#include "GlobalFitCommon.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TDirectory.h>
#include <TH3D.h>
#include <TObjArray.h>
#include <TVector3.h>

namespace {
using FibreKeyType = std::tuple<int,int,int>; // projection and transverse coordinates in 0.5 mm
struct Observation {
    double x=0,y=0,z=0,time=0,q=0;
    unsigned int geomId=0;
    int projection=0,u=0,v=0;
};
struct MapHit { TVector3 position; double fraction; };

FibreKeyType FibreKey(int projection, const TVector3& p) {
    if (projection == 0) return {0, std::lround(2*p.X()), std::lround(2*p.Z())};
    if (projection == 1) return {1, std::lround(2*p.Y()), std::lround(2*p.Z())};
    return {2, std::lround(2*p.X()), std::lround(2*p.Y())};
}

class LightMap {
public:
    explicit LightMap(const std::string& name,double minimumFraction) : file(name.c_str(), "READ") {
        if (file.IsZombie()) throw std::runtime_error("Cannot open light map " + name);
        hist = dynamic_cast<TH3D*>(file.Get("Starting Vertex Distribution"));
        centre = dynamic_cast<TVector3*>(file.Get("cubeCentre"));
        auto* pd = dynamic_cast<TDirectory*>(file.Get("hitPositions"));
        auto* fd = dynamic_cast<TDirectory*>(file.Get("lightFractions"));
        if (!hist || !centre || !pd || !fd) throw std::runtime_error("Incomplete light map");
        TIter next(pd->GetListOfKeys());
        while (auto* key = static_cast<::TKey*>(next())) {
            std::string n = key->GetName();
            const auto first = n.find('_') + 1, last = n.find('_', first);
            const int bin = std::stoi(n.substr(first, last-first));
            auto* positions = dynamic_cast<TObjArray*>(pd->Get(n.c_str()));
            auto* fractions = fd->Get<std::vector<double>>(
                ("Bin_"+std::to_string(bin)+"_light_fractions").c_str());
            if (!positions || !fractions) continue;
            auto& out = bins[bin];
            const int count = std::min<int>(positions->GetEntries(), fractions->size());
            for (int i=0;i<count;++i) if(fractions->at(i)>=minimumFraction)
                out.push_back({*static_cast<TVector3*>(positions->At(i)), fractions->at(i)});
        }
    }
    const std::vector<MapHit>& At(const TVector3& global, TVector3& offset) const {
        const TVector3 rel = global-*centre;
        TVector3 unit(Fold(rel.X()),Fold(rel.Y()),Fold(rel.Z()));
        offset.SetXYZ(global.X()-unit.X(),global.Y()-unit.Y(),global.Z()-unit.Z());
        const int bin=hist->FindBin(unit.X(),unit.Y(),unit.Z());
        auto it=bins.find(bin); return it==bins.end()?empty:it->second;
    }
private:
    static double Fold(double v) { return v>5?std::fmod(v+5,10)-5:v<-5?std::fmod(v-5,10)+5:v; }
    TFile file; TH3D* hist=nullptr; TVector3* centre=nullptr;
    std::map<int,std::vector<MapHit>> bins; std::vector<MapHit> empty;
};

std::vector<Observation> Read(TTree& t,int wanted,double cut) {
    int event,projection,u,v;unsigned int geomId;double x,y,z,time,q;
    t.SetBranchAddress("event",&event); t.SetBranchAddress("projection",&projection);
    t.SetBranchAddress("x",&x);t.SetBranchAddress("y",&y);t.SetBranchAddress("z",&z);
    t.SetBranchAddress("time",&time);t.SetBranchAddress("charge",&q);t.SetBranchAddress("geom_id",&geomId);
    t.SetBranchAddress("u",&u);t.SetBranchAddress("v",&v);std::vector<Observation> out;
    for(Long64_t i=0;i<t.GetEntries();++i){t.GetEntry(i);if(event<wanted)continue;if(event>wanted)break;
        if(q>=cut&&std::isfinite(q))out.push_back({x,y,z,time,q,geomId,projection,u,v});}
    return out;
}

std::pair<double,double> Transverse(const Observation& hit) {
    if(hit.projection==0) return {hit.x,hit.z}; // Y-directed fibre: XZ
    if(hit.projection==1) return {hit.y,hit.z}; // X-directed fibre: YZ
    return {hit.x,hit.y};                       // Z-directed fibre: XY
}

// DBSCAN is performed separately in each physical fibre projection. Keeping
// the highest-charge cluster avoids selecting a small but dense conversion.
std::vector<Observation> MainDbscanClusters(const std::vector<Observation>& input,
                                            double epsilon,int minPoints) {
    std::vector<Observation> output;const double epsilon2=epsilon*epsilon;
    for(int projection=0;projection<3;++projection) {
        std::vector<size_t> index;for(size_t i=0;i<input.size();++i)if(input[i].projection==projection)index.push_back(i);
        const int n=index.size();if(n==0)continue;std::vector<int> label(n,-1); // -1 unvisited, -2 noise
        auto neighbours=[&](int i){std::vector<int> found;auto a=Transverse(input[index[i]]);
            for(int j=0;j<n;++j){auto b=Transverse(input[index[j]]);double du=a.first-b.first,dv=a.second-b.second;
                if(du*du+dv*dv<=epsilon2+1e-6)found.push_back(j);}return found;};
        int clusters=0;
        for(int i=0;i<n;++i){if(label[i]!=-1)continue;auto nearby=neighbours(i);if(static_cast<int>(nearby.size())<minPoints){label[i]=-2;continue;}
            const int cluster=clusters++;label[i]=cluster;
            for(size_t cursor=0;cursor<nearby.size();++cursor){int j=nearby[cursor];if(label[j]==-2)label[j]=cluster;if(label[j]!=-1)continue;label[j]=cluster;
                auto expanded=neighbours(j);if(static_cast<int>(expanded.size())>=minPoints)for(int k:expanded)
                    if(std::find(nearby.begin(),nearby.end(),k)==nearby.end())nearby.push_back(k);}}
        if(clusters==0) continue;
        std::vector<double> charge(clusters,0.0);
        for(int i=0;i<n;++i) if(label[i]>=0) charge[label[i]]+=input[index[i]].q;
        const int keep=std::max_element(charge.begin(),charge.end())-charge.begin();for(int i=0;i<n;++i)if(label[i]==keep)output.push_back(input[index[i]]);
    }
    return output;
}

std::vector<Observation> LineCorridor(const std::vector<Observation>& input,
                                     const std::array<double,3>& direction,
                                     const TVector3& centre,int halfWidth) {
    std::vector<Observation> output;
    for(int projection=0;projection<3;++projection) {
        // Coordinate axes corresponding to Transverse(): XZ, YZ, and XY.
        const int axisA[3]={0,1,0},axisB[3]={2,2,1};
        double da=direction[axisA[projection]],db=direction[axisB[projection]];
        const double norm=std::hypot(da,db);if(norm<1e-9){for(auto& h:input)if(h.projection==projection)output.push_back(h);continue;}
        da/=norm;db/=norm;const double centreCoordinate[3]={centre.X(),centre.Y(),centre.Z()};
        // halfWidth=1 means the closest row and its immediate neighbours:
        // a maximum perpendicular centre-to-line distance of one pitch.
        const double maximumDistance=10.0*halfWidth+0.5;
        for(const auto& hit:input)if(hit.projection==projection){auto uv=Transverse(hit);double a=uv.first-centreCoordinate[axisA[projection]],b=uv.second-centreCoordinate[axisB[projection]];
            const double transverse=-a*db+b*da;if(std::abs(transverse)<=maximumDistance)output.push_back(hit);}
    }
    return output;
}

std::array<double,3> Seed(const std::vector<Observation>& o,const TVector3& centre) {
    std::array<double,3> f{{0,0,0}};
    const bool measured[3][3]={{true,false,true},{false,true,true},{true,true,false}};
    for(auto& h:o){double c[3]={h.x-centre.X(),h.y-centre.Y(),h.z-centre.Z()};for(int a=0;a<3;++a)
        if(measured[h.projection][a]&&std::abs(c[a])>std::abs(f[a]))f[a]=c[a];}
    double n=std::sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
    if(n==0) throw std::runtime_error("Zero direction seed");
    for(double& v:f) v/=n;
    return f;
}

class Likelihood {
public:
    Likelihood(const LightMap& m,const std::vector<Observation>& o,const std::array<double,3>& seed,
               const std::array<double,3>& lo,const std::array<double,3>& hi):map(m),low(lo),high(hi),origin(
                   .5*(lo[0]+hi[0]),.5*(lo[1]+hi[1]),.5*(lo[2]+hi[2])) {
        for(auto& h:o) observed[FibreKey(h.projection,TVector3(h.x,h.y,h.z))]+=h.q;
        total=0;for(auto& v:observed)total+=v.second;
        TVector3 d(seed[0],seed[1],seed[2]);
        TVector3 trial=std::abs(d.Z())<.9?TVector3(0,0,1):TVector3(0,1,0);
        e1=d.Cross(trial).Unit();e2=d.Cross(e1).Unit();
    }
    double operator()(const double* p) const {
        std::map<FibreKeyType,double> predicted;
        double sum=0;
        if(!Predict(p,predicted,sum)) return 1e100;
        double nll=0;for(auto& v:observed){double mu=1e-15;auto it=predicted.find(v.first);if(it!=predicted.end())mu=it->second/sum;
            nll-=v.second*std::log(mu);}return nll/total;
    }
    double Chi2(const double* p,int& ndof) const {
        std::map<FibreKeyType,double> predicted;double unusedSum=0;
        if(!Predict(p,predicted,unusedSum)){ndof=0;return 1e100;}
        // Condition the Pearson comparison on the same charge-selected fibre
        // set used by the likelihood. The overall normalization is not fitted.
        double selectedPrediction=0;
        for(const auto& v:observed){auto it=predicted.find(v.first);if(it!=predicted.end())selectedPrediction+=it->second;}
        if(selectedPrediction<=0){ndof=0;return 1e100;}
        double chi2=0;
        for(const auto& v:observed){auto it=predicted.find(v.first);double weight=it==predicted.end()?0.0:it->second;
            double expected=total*weight/selectedPrediction;
            double delta=v.second-expected;
            // Data-variance (Neyman) chi-square. The q>=cut selection makes
            // the Poisson variance estimate finite even for zero prediction.
            chi2+=delta*delta/v.second;}
        ndof=std::max(0,static_cast<int>(observed.size())-5);
        return chi2;
    }
    TVector3 Point(double a,double b)const{return origin+a*e1+b*e2;}
private:
    bool Predict(const double* p,std::map<FibreKeyType,double>& predicted,double& sum) const {
        TVector3 point=origin+p[0]*e1+p[1]*e2;
        TVector3 dir(std::sin(p[2])*std::cos(p[3]),std::sin(p[2])*std::sin(p[3]),std::cos(p[2]));
        double a=-1e30,b=1e30;
        for(int k=0;k<3;++k){double x=point[k],d=dir[k];if(std::abs(d)<1e-9){if(x<low[k]||x>high[k])return 1e100;continue;}
            double q1=(low[k]-x)/d,q2=(high[k]-x)/d;if(q1>q2)std::swap(q1,q2);a=std::max(a,q1);b=std::min(b,q2);}
        if(a>=b) return false;
        constexpr double step=2.0;
        for(double t=a+step/2;t<b;t+=step){TVector3 pos=point+t*dir,off;for(auto& h:map.At(pos,off)){
            TVector3 fp=h.position+off;if(fp.X()<low[0]||fp.X()>high[0]||fp.Y()<low[1]||fp.Y()>high[1]||fp.Z()<low[2]||fp.Z()>high[2])continue;
            int projection=std::abs(h.position.X())<.01?1:std::abs(h.position.Y())<.01?0:2;
            predicted[FibreKey(projection,fp)]+=h.fraction*step;}}
        sum=0;for(auto& v:predicted)sum+=v.second;
        return sum>0;
    }
    const LightMap& map; std::map<FibreKeyType,double> observed; double total=0;
    std::array<double,3> low,high; TVector3 origin,e1,e2;
};
}

int main(int argc,char** argv){
    auto usage=[](){std::cout<<R"(Usage:
  global_light_fit INPUT.root OUTPUT.root LIGHTMAP.root [KEY=value ...]

Options:
  EVENT=0|all                    Event number or all events (default: 0)
  TREE=homo_truth|homo_raw|fiber_hits
                                 Fibre input tree (default: homo_truth)
  MIN_CHARGE=10                  Minimum measured fibre charge
  DBSCAN=0|1                     Enable DBSCAN (default: 1)
  DBSCAN_EPSILON_MM=14.2         DBSCAN radius; includes diagonal fibres
  DBSCAN_MIN_POINTS=2            Minimum DBSCAN neighbourhood size
  CORRIDOR=0|1                   Enable straight seed corridor (default: 0)
  CORRIDOR_HALF_WIDTH_FIBRES=1   Corridor half-width; 1 means 10.5 mm
  MIN_MAP_FRACTION=0             Ignore smaller light-map fractions

Example:
  global_light_fit flat.root fitted.root lightmap.root EVENT=all \
    TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 \
    DBSCAN_MIN_POINTS=2 CORRIDOR=1 CORRIDOR_HALF_WIDTH_FIBRES=1 \
    MIN_MAP_FRACTION=1e-4
)";};
    if(argc==1||(argc==2&&std::string(argv[1])=="--help")){usage();return 0;}
    if(argc<4){usage();return 2;}
    try{
        std::map<std::string,std::string> options;const std::set<std::string> allowed={"EVENT","TREE","MIN_CHARGE","DBSCAN","DBSCAN_EPSILON_MM","DBSCAN_MIN_POINTS","CORRIDOR","CORRIDOR_HALF_WIDTH_FIBRES","MIN_MAP_FRACTION"};
        for(int i=4;i<argc;++i){std::string argument=argv[i];auto equal=argument.find('=');if(equal==std::string::npos)throw std::runtime_error("Optional arguments must use KEY=value: "+argument);
            std::string key=argument.substr(0,equal),value=argument.substr(equal+1);if(!allowed.count(key))throw std::runtime_error("Unknown option: "+key);if(value.empty())throw std::runtime_error("Empty value for "+key);options[key]=value;}
        auto get=[&](const std::string& key,const std::string& fallback){auto it=options.find(key);return it==options.end()?fallback:it->second;};
        auto boolean=[&](const std::string& key,bool fallback){std::string value=get(key,fallback?"1":"0");if(value=="1"||value=="true")return true;if(value=="0"||value=="false")return false;throw std::runtime_error(key+" must be 0, 1, true, or false");};
        std::string inName=argv[1],outName=argv[2],mapName=argv[3],treeName=get("TREE","homo_truth");
        std::string eventSpec=get("EVENT","0");double cut=std::stod(get("MIN_CHARGE","10"));
        const bool useDbscan=boolean("DBSCAN",true);
        const double dbscanEpsilon=std::stod(get("DBSCAN_EPSILON_MM","14.2"));
        const int dbscanMinPoints=std::stoi(get("DBSCAN_MIN_POINTS","2"));
        const bool useCorridor=boolean("CORRIDOR",false);
        const int corridorHalfWidth=std::stoi(get("CORRIDOR_HALF_WIDTH_FIBRES","1"));
        const double minimumMapFraction=std::stod(get("MIN_MAP_FRACTION","0"));
        if(dbscanEpsilon<=0||dbscanMinPoints<1||corridorHalfWidth<0||minimumMapFraction<0)throw std::runtime_error("Invalid selection parameters");
        if(inName==outName)throw std::runtime_error("Input and output files must differ");
        TFile input(inName.c_str(),"READ");if(input.IsZombie())throw std::runtime_error("Cannot open input");
        auto* tree=dynamic_cast<TTree*>(input.Get(treeName.c_str()));if(!tree)throw std::runtime_error("Missing input tree");
        std::vector<int> events;
        if(eventSpec=="all") { int value=0;tree->SetBranchAddress("event",&value);std::set<int> unique;
            for(Long64_t i=0;i<tree->GetEntries();++i){tree->GetEntry(i);unique.insert(value);}events.assign(unique.begin(),unique.end()); }
        else events.push_back(std::stoi(eventSpec));
        std::array<double,3> low,high; const char* axis[3]={"x","y","z"};
        for(int k=0;k<3;++k){low[k]=tree->GetMinimum(axis[k])-5;high[k]=tree->GetMaximum(axis[k])+5;}
        const TVector3 detectorCentre(.5*(low[0]+high[0]),.5*(low[1]+high[1]),.5*(low[2]+high[2]));
        LightMap map(mapName,minimumMapFraction);TFile output(outName.c_str(),"RECREATE");CloneFlatFile(input,output);output.cd();
        GlobalFitResult r;TTree result("global_fit","Global light-pattern straight-track fit");BranchGlobalFit(result,r);
        int selectedEvent=0,selectedProjection=0,selectedU=0,selectedV=0;unsigned int selectedGeomId=0;
        double selectedX=0,selectedY=0,selectedZ=0,selectedTime=0,selectedCharge=0;
        TTree selected("global_fit_fibres","Fibres selected for the global fit");selected.Branch("event",&selectedEvent);
        selected.Branch("x",&selectedX);selected.Branch("y",&selectedY);selected.Branch("z",&selectedZ);selected.Branch("time",&selectedTime);
        selected.Branch("charge",&selectedCharge);selected.Branch("geom_id",&selectedGeomId);selected.Branch("projection",&selectedProjection);
        selected.Branch("u",&selectedU);selected.Branch("v",&selectedV);
        for(int event:events){
            auto obs=Read(*tree,event,cut);const int observationsBefore=obs.size();if(useDbscan)obs=MainDbscanClusters(obs,dbscanEpsilon,dbscanMinPoints);
            const int observationsAfterDbscan=obs.size();if(obs.empty()){std::cerr<<"Skipping event "<<event<<": no observations survive DBSCAN\n";continue;}
            auto seed=Seed(obs,detectorCentre);if(useCorridor)obs=LineCorridor(obs,seed,detectorCentre,corridorHalfWidth);
            if(obs.empty()){std::cerr<<"Skipping event "<<event<<": no observations survive corridor\n";continue;}seed=Seed(obs,detectorCentre);
            for(const auto& hit:obs){selectedEvent=event;selectedX=hit.x;selectedY=hit.y;selectedZ=hit.z;selectedTime=hit.time;selectedCharge=hit.q;
                selectedGeomId=hit.geomId;selectedProjection=hit.projection;selectedU=hit.u;selectedV=hit.v;selected.Fill();}
            Likelihood likelihood(map,obs,seed,low,high);double theta=std::acos(std::clamp(seed[2],-1.0,1.0)),phi=std::atan2(seed[1],seed[0]);
            auto minimizer=std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex"));
            ROOT::Math::Functor functor(likelihood,4);minimizer->SetFunction(functor);minimizer->SetMaxFunctionCalls(2000);minimizer->SetTolerance(1e-3);
            minimizer->SetLimitedVariable(0,"impact_a",0,5,-500,500);minimizer->SetLimitedVariable(1,"impact_b",0,5,-500,500);
            minimizer->SetLimitedVariable(2,"theta",theta,.01,0.001,3.1405);minimizer->SetVariable(3,"phi",phi,.01);minimizer->Minimize();
            const double* p=minimizer->X();TVector3 point=likelihood.Point(p[0],p[1]);TVector3 dir(std::sin(p[2])*std::cos(p[3]),std::sin(p[2])*std::sin(p[3]),std::cos(p[2]));
            r=GlobalFitResult();r.event=event;r.status=minimizer->Status();r.observations=obs.size();r.observationsBeforeClustering=observationsBefore;r.observationsAfterDbscan=observationsAfterDbscan;
            r.dbscanEnabled=useDbscan;r.dbscanEpsilon=dbscanEpsilon;r.dbscanMinPoints=dbscanMinPoints;r.minimumCharge=cut;r.nll=minimizer->MinValue();r.chi2=likelihood.Chi2(p,r.ndof);r.inputTree=treeName;
            r.corridorEnabled=useCorridor;r.corridorHalfWidth=corridorHalfWidth;
            r.minimumMapFraction=minimumMapFraction;
            r.seedX=detectorCentre.X();r.seedY=detectorCentre.Y();r.seedZ=detectorCentre.Z();r.seedDx=seed[0];r.seedDy=seed[1];r.seedDz=seed[2];r.fitX=point.X();r.fitY=point.Y();r.fitZ=point.Z();r.fitDx=dir.X();r.fitDy=dir.Y();r.fitDz=dir.Z();
            r.errorA=minimizer->Errors()[0];r.errorB=minimizer->Errors()[1];r.errorTheta=minimizer->Errors()[2];r.errorPhi=minimizer->Errors()[3];result.Fill();
            std::cout<<"event "<<event<<" fibres "<<observationsBefore<<" -> "<<observationsAfterDbscan<<" -> "<<r.observations<<" status "<<r.status<<" nll "<<r.nll<<" chi2/ndof "<<r.chi2<<"/"<<r.ndof<<" = "<<(r.ndof>0?r.chi2/r.ndof:0)<<" direction "<<r.fitDx<<" "<<r.fitDy<<" "<<r.fitDz<<"\n";
        }
        output.cd();result.Write();selected.Write();output.Close();std::cout<<"output "<<outName<<"\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}return 0;
}
