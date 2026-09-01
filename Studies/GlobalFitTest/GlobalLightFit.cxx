#include "GlobalFitCommon.hxx"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
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
struct McSegment { int segment=0,primary=0;TVector3 start,stop; };
#ifdef GLOBAL_FIT_COLUMN_MODE
struct ColumnMipPriorConfig {
    bool enabled=false;double floorEnergy=0.5,mpvEnergy=1.5,strength=0;
    double chargePerEnergy=1,minPathFraction=0.8;int maximumIterations=500;
    double convergence=1e-6;
};
struct ColumnFitResult {
    int axis=0,index=0;double low=0,high=0,pathLength=0;
    double fittedCharge=0,chargePerMm=0,responseScale=0;
    double unregularizedFittedCharge=0,unregularizedResponseScale=0;
    int mipPriorEligible=0;double mipFloorCharge=0,mipMpvCharge=0;
};
struct ColumnFitSummary {
    std::vector<ColumnFitResult> columns;int converged=0,iterations=0,activeColumns=0,ndof=0;
    double maxRelativeChange=0,nll=0,chi2=0;
    int mipPriorEnabled=0,mipPriorConverged=0,mipPriorStalled=0,mipPriorIterations=0,mipPriorEligibleColumns=0;
    double unregularizedNll=0,mipPriorPenalty=0,mipPriorObjective=0,mipPriorMaxRelativeChange=0;
};
#endif

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

std::pair<std::set<int>,std::string> ParseFitViews(std::string value) {
    for(char& c:value)c=std::toupper(static_cast<unsigned char>(c));
    if(value=="ALL"||value=="3D")return {{0,1,2},"ALL"};
    if(value=="2D_X")return {{0,2},"XZ,XY"}; // exclude YZ / X-directed fibres
    if(value=="2D_Y")return {{1,2},"YZ,XY"}; // exclude XZ / Y-directed fibres
    if(value=="2D_Z")return {{0,1},"XZ,YZ"}; // exclude XY / Z-directed fibres
    std::set<int> views;std::string canonical;
    size_t begin=0;
    while(begin<value.size()){
        const size_t end=value.find(',',begin);const std::string name=value.substr(begin,end==std::string::npos?std::string::npos:end-begin);
        int projection=-1;if(name=="XZ"||name=="ZX")projection=0;else if(name=="YZ"||name=="ZY")projection=1;else if(name=="XY"||name=="YX")projection=2;
        else throw std::runtime_error("FIT_VIEWS entries must be XY, XZ, or YZ; got "+name);
        if(views.insert(projection).second){if(!canonical.empty())canonical+=",";canonical+=(projection==0?"XZ":projection==1?"YZ":"XY");}
        if(end==std::string::npos)break;
        begin=end+1;
    }
    if(views.size()<2)throw std::runtime_error("FIT_VIEWS needs at least two views to fit a 3D line");
    return {views,canonical};
}

std::vector<Observation> SelectFitViews(const std::vector<Observation>& input,const std::set<int>& views) {
    std::vector<Observation> output;for(const auto& hit:input)if(views.count(hit.projection))output.push_back(hit);return output;
}

double WeightedQuantile(std::vector<std::pair<double,double>> values,double fraction) {
    if(values.empty())throw std::runtime_error("Cannot determine fit range from an empty coordinate distribution");
    std::sort(values.begin(),values.end(),[](const auto& a,const auto& b){return a.first<b.first;});
    double totalWeight=0;for(const auto& value:values)totalWeight+=value.second;
    if(totalWeight<=0)throw std::runtime_error("Fit-range coordinate distribution has no positive charge");
    const double target=std::clamp(fraction,0.0,1.0)*totalWeight;double cumulative=0;
    for(const auto& value:values){cumulative+=value.second;if(cumulative>=target)return value.first;}
    return values.back().first;
}

std::pair<std::array<double,3>,std::array<double,3>> DataFitBounds(
        const std::vector<Observation>& observations,const std::array<double,3>& detectorLow,
        const std::array<double,3>& detectorHigh,double quantile,double padding) {
    std::array<std::vector<std::pair<double,double>>,3> coordinates;
    for(const auto& hit:observations){
        // A fibre constrains only coordinates transverse to its direction.
        if(hit.projection!=1)coordinates[0].push_back({hit.x,hit.q}); // XZ or XY
        if(hit.projection!=0)coordinates[1].push_back({hit.y,hit.q}); // YZ or XY
        if(hit.projection!=2)coordinates[2].push_back({hit.z,hit.q}); // XZ or YZ
    }
    std::array<double,3> low,high;
    for(int axis=0;axis<3;++axis){
        low[axis]=std::max(detectorLow[axis],WeightedQuantile(coordinates[axis],quantile)-padding);
        high[axis]=std::min(detectorHigh[axis],WeightedQuantile(coordinates[axis],1.0-quantile)+padding);
        if(high[axis]-low[axis]<10.0){const double centre=.5*(low[axis]+high[axis]);
            low[axis]=std::max(detectorLow[axis],centre-5.0);high[axis]=std::min(detectorHigh[axis],centre+5.0);}
    }
    return {low,high};
}

double SignedViewResidual(const TVector3& point,const TVector3& linePoint,
                          const TVector3& direction,int axisA,int axisB) {
    const double da=direction[axisA],db=direction[axisB],norm=std::hypot(da,db);
    if(norm<1e-12)return std::numeric_limits<double>::quiet_NaN();
    const double a=point[axisA]-linePoint[axisA],b=point[axisB]-linePoint[axisB];
    return (-db*a+da*b)/norm;
}

// Fit an orthogonal-regression (PCA) line to the muon segment start points.
// This benchmark is deliberately independent of fibres and of the light map.
void WriteMcSegmentBenchmark(TFile& output,TTree& segments,TTree& tracks,
                             const std::vector<int>& requestedEvents) {
    const std::set<int> wanted(requestedEvents.begin(),requestedEvents.end());
    int event=0,trackId=0,pdg=0;tracks.SetBranchAddress("event",&event);tracks.SetBranchAddress("track_id",&trackId);tracks.SetBranchAddress("pdg",&pdg);
    std::map<int,std::set<int>> muonIds;
    for(Long64_t i=0;i<tracks.GetEntries();++i){tracks.GetEntry(i);if(wanted.count(event)&&std::abs(pdg)==13)muonIds[event].insert(trackId);}
    tracks.ResetBranchAddresses();

    int detector=0,segment=0,primary=0;double sx=0,sy=0,sz=0,ex=0,ey=0,ez=0;
    segments.SetBranchAddress("event",&event);segments.SetBranchAddress("segment",&segment);segments.SetBranchAddress("detector",&detector);segments.SetBranchAddress("primary_id",&primary);
    segments.SetBranchAddress("start_x",&sx);segments.SetBranchAddress("start_y",&sy);segments.SetBranchAddress("start_z",&sz);
    segments.SetBranchAddress("stop_x",&ex);segments.SetBranchAddress("stop_y",&ey);segments.SetBranchAddress("stop_z",&ez);
    std::map<std::pair<int,int>,std::vector<McSegment>> byTrack;
    for(Long64_t i=0;i<segments.GetEntries();++i){segments.GetEntry(i);if(detector!=0||!wanted.count(event)||!muonIds[event].count(primary))continue;
        byTrack[{event,primary}].push_back({segment,primary,TVector3(sx,sy,sz),TVector3(ex,ey,ez)});}
    segments.ResetBranchAddresses();

    std::map<int,std::vector<McSegment>> selected;
    for(auto& item:byTrack){const int e=item.first.first;if(!selected.count(e)||item.second.size()>selected[e].size())selected[e]=std::move(item.second);}

    int fitEvent=0,fitPrimary=0,fitStarts=0,fitNdof=0;double fitX=0,fitY=0,fitZ=0,fitDx=0,fitDy=0,fitDz=0,fitSum2=0,fitRms=0;
    TTree fitTree("mc_segment_line_fit","Straight-line PCA fit to muon MC segment start positions");
    fitTree.Branch("event",&fitEvent);fitTree.Branch("primary_id",&fitPrimary);fitTree.Branch("n_start_points",&fitStarts);fitTree.Branch("ndof",&fitNdof);
    fitTree.Branch("fit_x",&fitX);fitTree.Branch("fit_y",&fitY);fitTree.Branch("fit_z",&fitZ);fitTree.Branch("fit_dx",&fitDx);fitTree.Branch("fit_dy",&fitDy);fitTree.Branch("fit_dz",&fitDz);
    fitTree.Branch("sum_squared_residual_mm2",&fitSum2);fitTree.Branch("rms_residual_mm",&fitRms);

    int residualEvent=0,residualSegment=0,residualPrimary=0;double startDistance=0,endDistance=0;
    double startRx=0,startRy=0,startRz=0,endRx=0,endRy=0,endRz=0;
    double startXY=0,startXZ=0,startYZ=0,endXY=0,endXZ=0,endYZ=0;
    TTree residualTree("mc_segment_line_residuals","MC segment start/end residuals to the MC-start straight-line fit");
    residualTree.Branch("event",&residualEvent);residualTree.Branch("segment",&residualSegment);residualTree.Branch("primary_id",&residualPrimary);
    residualTree.Branch("start_distance",&startDistance);residualTree.Branch("end_distance",&endDistance);
    residualTree.Branch("start_residual_x",&startRx);residualTree.Branch("start_residual_y",&startRy);residualTree.Branch("start_residual_z",&startRz);
    residualTree.Branch("end_residual_x",&endRx);residualTree.Branch("end_residual_y",&endRy);residualTree.Branch("end_residual_z",&endRz);
    residualTree.Branch("start_residual_xy",&startXY);residualTree.Branch("start_residual_xz",&startXZ);residualTree.Branch("start_residual_yz",&startYZ);
    residualTree.Branch("end_residual_xy",&endXY);residualTree.Branch("end_residual_xz",&endXZ);residualTree.Branch("end_residual_yz",&endYZ);
    TH1D startHistogram("mc_segment_start_distance",";MC segment start distance to MC-start line [mm];points",500,0,500);
    TH1D endHistogram("mc_segment_end_distance",";MC segment end distance to MC-start line [mm];points",500,0,500);

    for(auto& item:selected){auto& values=item.second;if(values.size()<3){std::cerr<<"Skipping MC line benchmark event "<<item.first<<": fewer than 3 muon segments\n";continue;}
        TVector3 centre;for(const auto& value:values)centre+=value.start;centre*=1.0/values.size();
        TMatrixDSym covariance(3);covariance.Zero();TVector3 travel;
        for(const auto& value:values){const TVector3 delta=value.start-centre;travel+=value.stop-value.start;for(int a=0;a<3;++a)for(int b=0;b<3;++b)covariance(a,b)+=delta[a]*delta[b];}
        TMatrixDSymEigen eigen(covariance);const TVectorD eigenvalues=eigen.GetEigenValues();const TMatrixD eigenvectors=eigen.GetEigenVectors();int largest=0;
        for(int k=1;k<3;++k)if(eigenvalues[k]>eigenvalues[largest])largest=k;
        TVector3 direction(eigenvectors(0,largest),eigenvectors(1,largest),eigenvectors(2,largest));direction=direction.Unit();
        if(travel.Dot(direction)<0)direction*=-1;
        fitEvent=item.first;fitPrimary=values.front().primary;fitStarts=values.size();fitNdof=std::max(0,2*fitStarts-4);fitX=centre.X();fitY=centre.Y();fitZ=centre.Z();fitDx=direction.X();fitDy=direction.Y();fitDz=direction.Z();fitSum2=0;
        for(const auto& value:values)fitSum2+=((value.start-centre).Cross(direction)).Mag2();
        fitRms=std::sqrt(fitSum2/fitStarts);fitTree.Fill();
        for(const auto& value:values){auto residual=[&](const TVector3& point){const TVector3 delta=point-centre;return delta-direction*delta.Dot(direction);};
            const TVector3 sr=residual(value.start),er=residual(value.stop);residualEvent=item.first;residualSegment=value.segment;residualPrimary=value.primary;
            startDistance=sr.Mag();endDistance=er.Mag();startRx=sr.X();startRy=sr.Y();startRz=sr.Z();endRx=er.X();endRy=er.Y();endRz=er.Z();
            startXY=SignedViewResidual(value.start,centre,direction,0,1);startXZ=SignedViewResidual(value.start,centre,direction,0,2);startYZ=SignedViewResidual(value.start,centre,direction,1,2);
            endXY=SignedViewResidual(value.stop,centre,direction,0,1);endXZ=SignedViewResidual(value.stop,centre,direction,0,2);endYZ=SignedViewResidual(value.stop,centre,direction,1,2);
            startHistogram.Fill(startDistance);endHistogram.Fill(endDistance);residualTree.Fill();}}
    output.cd();fitTree.Write();residualTree.Write();startHistogram.Write();endHistogram.Write();
    std::cout<<"MC segment benchmark: "<<fitTree.GetEntries()<<" line fits and "<<residualTree.GetEntries()<<" segment residual rows\n";
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

struct MedianDiameterSeedResult {
    std::array<double,3> direction{{0,0,1}};
    std::array<double,3> median{{0,0,0}}; // XZ, YZ, XY by projection number
    int observations=0;
};

MedianDiameterSeedResult MedianDiameterSeed(const std::vector<Observation>& input,double factor) {
    MedianDiameterSeedResult result;TMatrixDSym constraint(3);constraint.Zero();int usableViews=0;
    const int axisA[3]={0,1,0},axisB[3]={2,2,1};
    for(int projection=0;projection<3;++projection){
        std::vector<const Observation*> view;for(const auto& hit:input)if(hit.projection==projection)view.push_back(&hit);
        if(view.empty())continue;
        std::vector<double> charges;for(const auto* hit:view)charges.push_back(hit->q);
        const size_t middle=charges.size()/2;std::nth_element(charges.begin(),charges.begin()+middle,charges.end());double median=charges[middle];
        if(charges.size()%2==0){const auto lower=std::max_element(charges.begin(),charges.begin()+middle);median=.5*(median+*lower);}result.median[projection]=median;
        std::vector<const Observation*> retained;for(const auto* hit:view)if(hit->q>=factor*median)retained.push_back(hit);result.observations+=retained.size();
        if(retained.size()<2)continue;
        double maximum=-1;const Observation *first=nullptr,*second=nullptr;
        for(size_t i=0;i<retained.size();++i)for(size_t j=i+1;j<retained.size();++j){auto a=Transverse(*retained[i]),b=Transverse(*retained[j]);const double du=b.first-a.first,dv=b.second-a.second,d2=du*du+dv*dv;
            if(d2>maximum){maximum=d2;first=retained[i];second=retained[j];}}
        if(!first||!second||maximum<=0)continue;
        auto a=Transverse(*first),b=Transverse(*second);const double da=b.first-a.first,db=b.second-a.second,norm=std::hypot(da,db);
        double normal[3]={0,0,0};normal[axisA[projection]]=-db/norm;normal[axisB[projection]]=da/norm;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)constraint(i,j)+=normal[i]*normal[j];
        ++usableViews;
    }
    if(usableViews<2)throw std::runtime_error("Median-diameter seed needs at least two usable views");
    TMatrixDSymEigen eigen(constraint);const TVectorD values=eigen.GetEigenValues();const TMatrixD vectors=eigen.GetEigenVectors();int smallest=0;
    for(int i=1;i<3;++i)if(values[i]<values[smallest])smallest=i;
    TVector3 direction(vectors(0,smallest),vectors(1,smallest),vectors(2,smallest));direction=direction.Unit();
    int dominant=0;for(int i=1;i<3;++i)if(std::abs(direction[i])>std::abs(direction[dominant]))dominant=i;if(direction[dominant]<0)direction*=-1;
    result.direction={{direction.X(),direction.Y(),direction.Z()}};return result;
}

// Truth-assisted diagnostic seed. This deliberately uses MC information and
// must never be used as a reconstruction-performance result.
std::map<int,std::array<double,3>> McSegment10Seeds(
        TTree& segments,TTree& tracks,const std::vector<int>& requestedEvents) {
    const std::set<int> wanted(requestedEvents.begin(),requestedEvents.end());
    int event=0,trackId=0,pdg=0;
    tracks.SetBranchAddress("event",&event);tracks.SetBranchAddress("track_id",&trackId);tracks.SetBranchAddress("pdg",&pdg);
    std::map<int,std::set<int>> muonIds;
    for(Long64_t i=0;i<tracks.GetEntries();++i){tracks.GetEntry(i);if(wanted.count(event)&&std::abs(pdg)==13)muonIds[event].insert(trackId);}
    tracks.ResetBranchAddresses();

    int detector=0,segment=0,primary=0;double sx=0,sy=0,sz=0,ex=0,ey=0,ez=0;
    segments.SetBranchAddress("event",&event);segments.SetBranchAddress("segment",&segment);segments.SetBranchAddress("detector",&detector);segments.SetBranchAddress("primary_id",&primary);
    segments.SetBranchAddress("start_x",&sx);segments.SetBranchAddress("start_y",&sy);segments.SetBranchAddress("start_z",&sz);
    segments.SetBranchAddress("stop_x",&ex);segments.SetBranchAddress("stop_y",&ey);segments.SetBranchAddress("stop_z",&ez);
    std::map<std::pair<int,int>,std::vector<McSegment>> byTrack;
    for(Long64_t i=0;i<segments.GetEntries();++i){segments.GetEntry(i);if(detector!=0||!wanted.count(event)||!muonIds[event].count(primary))continue;
        byTrack[{event,primary}].push_back({segment,primary,TVector3(sx,sy,sz),TVector3(ex,ey,ez)});}
    segments.ResetBranchAddresses();

    std::map<int,std::vector<McSegment>> selected;
    for(auto& item:byTrack){auto& values=item.second;std::sort(values.begin(),values.end(),[](const auto& a,const auto& b){return a.segment<b.segment;});
        const int e=item.first.first;if(!selected.count(e)||values.size()>selected[e].size())selected[e]=values;}
    std::map<int,std::array<double,3>> result;
    for(const auto& item:selected){const auto& values=item.second;if(values.size()<2)continue;
        const McSegment* target=nullptr;for(const auto& value:values)if(value.segment==10){target=&value;break;}
        if(!target&&values.size()>=10)target=&values[9];
        if(!target)continue;
        TVector3 direction=target->start-values.front().start;if(direction.Mag2()<=0)continue;direction=direction.Unit();
        result[item.first]={{direction.X(),direction.Y(),direction.Z()}};}
    return result;
}

class Likelihood {
public:
    Likelihood(const LightMap& m,const std::vector<Observation>& o,const std::array<double,3>& seed,
               const std::array<double,3>& lo,const std::array<double,3>& hi,const std::set<int>& views,
               const std::array<double,3>& detectorLo
#ifdef GLOBAL_FIT_COLUMN_MODE
               ,const ColumnMipPriorConfig& columnPriorConfig
#endif
               ):map(m),enabledViews(views),low(lo),high(hi),columnGridLow(detectorLo),origin(
                   .5*(lo[0]+hi[0]),.5*(lo[1]+hi[1]),.5*(lo[2]+hi[2])) {
        for(auto& h:o) observed[FibreKey(h.projection,TVector3(h.x,h.y,h.z))]+=h.q;
        total=0;for(auto& v:observed)total+=v.second;
        TVector3 d(seed[0],seed[1],seed[2]);
        TVector3 trial=std::abs(d.Z())<.9?TVector3(0,0,1):TVector3(0,1,0);
        e1=d.Cross(trial).Unit();e2=d.Cross(e1).Unit();
#ifdef GLOBAL_FIT_COLUMN_MODE
        columnPrior=columnPriorConfig;
        columnAxis=0;if(std::abs(d.Y())>std::abs(d[columnAxis]))columnAxis=1;
        if(std::abs(d.Z())>std::abs(d[columnAxis]))columnAxis=2;
#endif
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
#ifdef GLOBAL_FIT_COLUMN_MODE
    ColumnFitSummary FitColumnCharges(const double* p) const {
        ColumnFitSummary summary;
        std::map<int,std::map<FibreKeyType,double>> components;std::map<int,double> lengths;
        if(!PredictColumns(p,components,lengths))return summary;
        std::vector<int> columns;for(const auto& item:components)columns.push_back(item.first);
        std::vector<double> amplitude(columns.size(),0.0),normalization(columns.size(),0.0);
        double initialResponse=0;
        for(size_t c=0;c<columns.size();++c){
            for(const auto& value:components.at(columns[c]))if(observed.count(value.first))normalization[c]+=value.second;
            initialResponse+=lengths.at(columns[c])*normalization[c];
        }
        if(initialResponse<=0)return summary;
        for(size_t c=0;c<columns.size();++c)amplitude[c]=total*lengths.at(columns[c])/initialResponse;
        // Multiplicative Poisson/EM updates give a non-negative solution for
        // fixed geometry. Only the charge-selected fibre set is conditioned
        // on, matching the geometry likelihood's observations.
        for(int iteration=0;iteration<2000;++iteration){
            std::map<FibreKeyType,double> expected;
            for(size_t c=0;c<columns.size();++c)for(const auto& value:components.at(columns[c]))
                if(observed.count(value.first))expected[value.first]+=amplitude[c]*value.second;
            double largestRelativeChange=0;std::vector<double> updated=amplitude;
            for(size_t c=0;c<columns.size();++c){
                if(normalization[c]<=0){updated[c]=0;continue;}double numerator=0;
                for(const auto& value:components.at(columns[c])){auto data=observed.find(value.first);if(data==observed.end())continue;
                    const double mu=expected[value.first];if(mu>0)numerator+=value.second*data->second/mu;}
                updated[c]=amplitude[c]*numerator/normalization[c];
                largestRelativeChange=std::max(largestRelativeChange,std::abs(updated[c]-amplitude[c])/std::max(amplitude[c],1e-12));
            }
            amplitude.swap(updated);summary.iterations=iteration+1;summary.maxRelativeChange=largestRelativeChange;
            if(largestRelativeChange<1e-6){summary.converged=1;break;}
        }
        const std::vector<double> unregularizedAmplitude=amplitude;
        auto expectedFrom=[&](const std::vector<double>& values){std::map<FibreKeyType,double> expected;
            for(size_t c=0;c<columns.size();++c)for(const auto& value:components.at(columns[c]))if(observed.count(value.first))expected[value.first]+=values[c]*value.second;
            return expected;};
        auto fullPoisson=[&](const std::vector<double>& values){const auto expected=expectedFrom(values);double value=0;
            for(const auto& data:observed){auto found=expected.find(data.first);const double mu=found==expected.end()?1e-15:std::max(found->second,1e-15);value+=mu-data.second*std::log(mu);}return value/std::max(total,1e-15);};
        auto normalizedShapeNll=[&](const std::vector<double>& values){const auto expected=expectedFrom(values);double predicted=0,nll=0;
            for(const auto& data:observed){auto found=expected.find(data.first);if(found!=expected.end())predicted+=found->second;}
            if(predicted<=0)return 1e100;
            for(const auto& data:observed){auto found=expected.find(data.first);const double mu=found==expected.end()?1e-15:std::max(found->second/predicted,1e-15);nll-=data.second*std::log(mu);}return nll/std::max(total,1e-15);};
        summary.unregularizedNll=normalizedShapeNll(unregularizedAmplitude);

        std::vector<int> eligible(columns.size(),0);std::vector<double> floorCharge(columns.size(),0),mpvCharge(columns.size(),0);
        if(columnPrior.enabled&&columns.size()>2){std::vector<double> internalLengths;for(size_t c=1;c+1<columns.size();++c)internalLengths.push_back(lengths.at(columns[c]));
            std::sort(internalLengths.begin(),internalLengths.end());const double medianLength=internalLengths[internalLengths.size()/2];
            for(size_t c=1;c+1<columns.size();++c){const double length=lengths.at(columns[c]);if(length+1e-12<columnPrior.minPathFraction*medianLength)continue;
                eligible[c]=1;++summary.mipPriorEligibleColumns;floorCharge[c]=columnPrior.floorEnergy*(length/10.0)*columnPrior.chargePerEnergy;
                mpvCharge[c]=columnPrior.mpvEnergy*(length/10.0)*columnPrior.chargePerEnergy;}}
        summary.mipPriorEnabled=columnPrior.enabled?1:0;
        auto penalty=[&](const std::vector<double>& values){double value=0;for(size_t c=0;c<columns.size();++c)if(eligible[c]){const double visible=values[c]*normalization[c];
                if(visible<mpvCharge[c]){const double width=std::max(mpvCharge[c]-floorCharge[c],1e-12);const double pull=(mpvCharge[c]-visible)/width;value+=pull*pull;}}
            return value/std::max(summary.mipPriorEligibleColumns,1);};
        if(columnPrior.enabled&&summary.mipPriorEligibleColumns>0&&columnPrior.strength>0){
            double objective=fullPoisson(amplitude)+columnPrior.strength*penalty(amplitude);
            auto scaledLargestChange=[&](const std::vector<double>& from,const std::vector<double>& to){double largest=0;
                for(size_t c=0;c<from.size();++c){double scale=std::max(std::abs(from[c]),1e-12);
                    // Near the non-negative boundary, dividing by a nearly
                    // zero amplitude makes an immaterial step look enormous.
                    // Use the prior's physical charge scale where available,
                    // and the unregularized solution elsewhere.
                    if(eligible[c]&&normalization[c]>0)scale=std::max(scale,mpvCharge[c]/normalization[c]);
                    else scale=std::max(scale,std::abs(unregularizedAmplitude[c]));
                    largest=std::max(largest,std::abs(to[c]-from[c])/std::max(scale,1e-12));}
                return largest;};
            for(int iteration=0;iteration<columnPrior.maximumIterations;++iteration){const auto expected=expectedFrom(amplitude);std::vector<double> gradient(columns.size(),0),hessian(columns.size(),0);
                for(size_t c=0;c<columns.size();++c){for(const auto& response:components.at(columns[c])){auto data=observed.find(response.first);if(data==observed.end())continue;auto muEntry=expected.find(response.first);const double mu=muEntry==expected.end()?1e-15:std::max(muEntry->second,1e-15);
                        gradient[c]+=response.second*(1-data->second/mu)/std::max(total,1e-15);hessian[c]+=response.second*response.second*data->second/(mu*mu*std::max(total,1e-15));}
                    if(eligible[c]){const double visible=amplitude[c]*normalization[c];if(visible<mpvCharge[c]){const double width=std::max(mpvCharge[c]-floorCharge[c],1e-12),scale=normalization[c]/width;
                            gradient[c]+=-columnPrior.strength*2*(mpvCharge[c]-visible)*normalization[c]/(width*width*summary.mipPriorEligibleColumns);
                            hessian[c]+=columnPrior.strength*2*scale*scale/summary.mipPriorEligibleColumns;}}}
                std::vector<double> proposal(amplitude.size());for(size_t c=0;c<amplitude.size();++c)proposal[c]=std::max(1e-15,amplitude[c]-gradient[c]/std::max(hessian[c],1e-12));
                double step=1;bool accepted=false,belowTolerance=false;std::vector<double> trial;double trialObjective=0,largest=0;
                for(int lineSearch=0;lineSearch<60;++lineSearch){trial.resize(amplitude.size());for(size_t c=0;c<amplitude.size();++c)trial[c]=std::max(1e-15,amplitude[c]+step*(proposal[c]-amplitude[c]));
                    largest=scaledLargestChange(amplitude,trial);
                    trialObjective=fullPoisson(trial)+columnPrior.strength*penalty(trial);if(trialObjective<objective-1e-12){accepted=true;break;}
                    if(largest<columnPrior.convergence){belowTolerance=true;break;}step*=.5;}
                if(!accepted){if(belowTolerance)summary.mipPriorConverged=1;else summary.mipPriorStalled=1;break;}
                amplitude.swap(trial);objective=trialObjective;summary.mipPriorIterations=iteration+1;summary.mipPriorMaxRelativeChange=largest;
                if(largest<columnPrior.convergence){summary.mipPriorConverged=1;break;}}
        }
        summary.mipPriorPenalty=penalty(amplitude);summary.mipPriorObjective=fullPoisson(amplitude)+columnPrior.strength*summary.mipPriorPenalty;
        std::map<FibreKeyType,double> finalExpected;double predictedTotal=0;
        for(size_t c=0;c<columns.size();++c)for(const auto& value:components.at(columns[c]))if(observed.count(value.first))
            finalExpected[value.first]+=amplitude[c]*value.second;
        for(const auto& value:finalExpected)predictedTotal+=value.second;
        if(predictedTotal>0)for(const auto& value:observed){const double mu=finalExpected[value.first];
            summary.nll-=value.second*std::log(std::max(mu/predictedTotal,1e-15))/total;
            const double delta=value.second-mu;summary.chi2+=delta*delta/value.second;}
        for(size_t c=0;c<columns.size();++c){const int index=columns[c];const double length=lengths.at(index);
            const double visible=amplitude[c]*normalization[c];
            if(visible>1e-6)++summary.activeColumns;
            const double physicalLow=columnGridLow[columnAxis]+10.0*index,physicalHigh=physicalLow+10.0;
            const int physicalIndex=static_cast<int>(std::lround(.5*(physicalLow+physicalHigh)/10.0));
            const double unregularizedVisible=unregularizedAmplitude[c]*normalization[c];
            summary.columns.push_back({columnAxis,physicalIndex,physicalLow,physicalHigh,length,
                              visible,length>0?visible/length:0.0,amplitude[c],unregularizedVisible,
                              unregularizedAmplitude[c],eligible[c],floorCharge[c],mpvCharge[c]});}
        summary.ndof=std::max(0,static_cast<int>(observed.size())-summary.activeColumns);
        return summary;
    }
#endif
private:
    bool Predict(const double* p,std::map<FibreKeyType,double>& predicted,double& sum) const {
#ifdef GLOBAL_FIT_COLUMN_MODE
        std::map<int,std::map<FibreKeyType,double>> columns;std::map<int,double> lengths;
        if(!PredictColumns(p,columns,lengths))return false;
        for(const auto& column:columns)for(const auto& value:column.second)predicted[value.first]+=value.second;
        sum=0;for(const auto& value:predicted)sum+=value.second;return sum>0;
#else
        TVector3 point=origin+p[0]*e1+p[1]*e2;
        TVector3 dir(std::sin(p[2])*std::cos(p[3]),std::sin(p[2])*std::sin(p[3]),std::cos(p[2]));
        double a=-1e30,b=1e30;
        for(int k=0;k<3;++k){double x=point[k],d=dir[k];if(std::abs(d)<1e-9){if(x<low[k]||x>high[k])return 1e100;continue;}
            double q1=(low[k]-x)/d,q2=(high[k]-x)/d;if(q1>q2)std::swap(q1,q2);a=std::max(a,q1);b=std::min(b,q2);}
        if(a>=b) return false;
        // Match the 1 mm light-map subvoxel pitch.  A 2 mm midpoint step can
        // systematically skip alternating response bins near a fibre.
        constexpr double step=1.0;
        for(double t=a+step/2;t<b;t+=step){TVector3 pos=point+t*dir,off;for(auto& h:map.At(pos,off)){
            TVector3 fp=h.position+off;if(fp.X()<low[0]||fp.X()>high[0]||fp.Y()<low[1]||fp.Y()>high[1]||fp.Z()<low[2]||fp.Z()>high[2])continue;
            int projection=std::abs(h.position.X())<.01?1:std::abs(h.position.Y())<.01?0:2;
            if(!enabledViews.count(projection))continue;
            predicted[FibreKey(projection,fp)]+=h.fraction*step;}}
        sum=0;for(auto& v:predicted)sum+=v.second;
        return sum>0;
#endif
    }
#ifdef GLOBAL_FIT_COLUMN_MODE
    bool PredictColumns(const double* p,std::map<int,std::map<FibreKeyType,double>>& columns,
                        std::map<int,double>& lengths) const {
        const TVector3 point=origin+p[0]*e1+p[1]*e2;
        const TVector3 dir(std::sin(p[2])*std::cos(p[3]),std::sin(p[2])*std::sin(p[3]),std::cos(p[2]));
        double trackBegin=-1e30,trackEnd=1e30;
        for(int k=0;k<3;++k){const double x=point[k],d=dir[k];if(std::abs(d)<1e-9){if(x<low[k]||x>high[k])return false;continue;}
            double q1=(low[k]-x)/d,q2=(high[k]-x)/d;if(q1>q2)std::swap(q1,q2);trackBegin=std::max(trackBegin,q1);trackEnd=std::min(trackEnd,q2);}
        if(trackBegin>=trackEnd)return false;
        const int firstColumn=static_cast<int>(std::floor((low[columnAxis]-columnGridLow[columnAxis])/10.0));
        const int lastColumn=static_cast<int>(std::ceil((high[columnAxis]-columnGridLow[columnAxis])/10.0))-1;
        for(int column=firstColumn;column<=lastColumn;++column){const double columnLow=std::max(low[columnAxis],columnGridLow[columnAxis]+10.0*column);
            const double columnHigh=std::min(high[columnAxis],columnGridLow[columnAxis]+10.0*(column+1));double begin=trackBegin,end=trackEnd;
            if(std::abs(dir[columnAxis])<1e-9){if(point[columnAxis]<columnLow||point[columnAxis]>=columnHigh)continue;}
            else {double q1=(columnLow-point[columnAxis])/dir[columnAxis],q2=(columnHigh-point[columnAxis])/dir[columnAxis];if(q1>q2)std::swap(q1,q2);begin=std::max(begin,q1);end=std::min(end,q2);}
            if(begin>=end)continue;
            lengths[column]=end-begin;
            constexpr double maximumStep=1.0;const int steps=std::max(1,static_cast<int>(std::ceil((end-begin)/maximumStep)));
            const double step=(end-begin)/steps;
            for(int sample=0;sample<steps;++sample){const TVector3 pos=point+(begin+(sample+.5)*step)*dir;TVector3 off;
                for(const auto& h:map.At(pos,off)){const TVector3 fp=h.position+off;
                    if(fp.X()<low[0]||fp.X()>high[0]||fp.Y()<low[1]||fp.Y()>high[1]||fp.Z()<low[2]||fp.Z()>high[2])continue;
                    const int projection=std::abs(h.position.X())<.01?1:std::abs(h.position.Y())<.01?0:2;
                    if(enabledViews.count(projection))columns[column][FibreKey(projection,fp)]+=h.fraction*step;}}
        }
        return !columns.empty();
    }
#endif
    const LightMap& map; const std::set<int> enabledViews; std::map<FibreKeyType,double> observed; double total=0;
    std::array<double,3> low,high,columnGridLow; TVector3 origin,e1,e2;
#ifdef GLOBAL_FIT_COLUMN_MODE
    int columnAxis=0;ColumnMipPriorConfig columnPrior;
#endif
};
}

#ifndef GLOBAL_LIGHT_FIT_NO_MAIN
int main(int argc,char** argv){
    auto usage=[](){std::cout<<"Usage:\n";
#ifdef GLOBAL_FIT_COLUMN_MODE
        std::cout<<"  global_light_fit_columns INPUT.root OUTPUT.root LIGHTMAP.root [KEY=value ...]\n";
#else
        std::cout<<"  global_light_fit INPUT.root OUTPUT.root LIGHTMAP.root [KEY=value ...]\n";
#endif
        std::cout<<R"(
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
  FIT_VIEWS=ALL|2D_X|2D_Y|2D_Z  Views entering selection and likelihood
                                 ALL uses XY+XZ+YZ (default)
                                 2D_X uses XZ+XY; excludes YZ/X fibres
                                 2D_Y uses YZ+XY; excludes XZ/Y fibres
                                 2D_Z uses XZ+YZ; excludes XY/Z fibres
                                 Explicit lists such as FIT_VIEWS=XY,XZ work too
  SEED_DIRECTION=FIBRE|VIEW_MEDIAN_DIAMETER|MC_SEGMENT_10
                                 Initial direction (default: FIBRE)
  SEED_MEDIAN_FACTOR=1           Per-view median multiplier used by
                                 VIEW_MEDIAN_DIAMETER after DBSCAN
  MAX_FUNCTION_CALLS=5000        Maximum likelihood evaluations per event
  TOLERANCE=1e-2                 Minuit2 minimization tolerance
  FIT_RANGE=0|1                  Use data-derived XYZ range (default: 1)
  FIT_RANGE_QUANTILE=0.01        Charge fraction trimmed at each range end
  FIT_RANGE_PADDING_MM=10        Padding around data-derived XYZ ranges
)";
#ifdef GLOBAL_FIT_COLUMN_MODE
        std::cout<<R"(
  COLUMN_MIP_PRIOR=0             Enable one-sided internal-column MIP prior
  COLUMN_MIP_FLOOR_ENERGY=0.5    Low-charge reference in calibrated energy units
  COLUMN_MIP_MPV_ENERGY=1.5      MIP most-probable charge; no penalty above it
  COLUMN_MIP_STRENGTH=0          Penalty weight (must be >0 when enabled)
  COLUMN_CHARGE_PER_ENERGY=1     Fixed independent charge/energy calibration
  COLUMN_MIP_MIN_PATH_FRACTION=0.8  Minimum path/median path; endpoints excluded
  COLUMN_MIP_MAX_ITERATIONS=500  Maximum regularized Newton updates
  COLUMN_MIP_CONVERGENCE=1e-6    Largest relative amplitude-change criterion
)";
#endif
        std::cout<<R"(

Example:
  global_light_fit flat.root fitted.root lightmap.root EVENT=all \
    TREE=fiber_hits MIN_CHARGE=10 DBSCAN=1 DBSCAN_EPSILON_MM=14.2 \
    DBSCAN_MIN_POINTS=2 CORRIDOR=1 CORRIDOR_HALF_WIDTH_FIBRES=1 \
    MIN_MAP_FRACTION=1e-4 FIT_VIEWS=2D_X

Truth-assisted seed diagnostic (not reconstruction performance):
  global_light_fit flat.root fitted_mc_seed.root lightmap.root EVENT=all \
    TREE=fiber_hits FIT_VIEWS=2D_X SEED_DIRECTION=MC_SEGMENT_10
)";};
    if(argc==1||(argc==2&&std::string(argv[1])=="--help")){usage();return 0;}
    if(argc<4){usage();return 2;}
    try{
        std::map<std::string,std::string> options;std::set<std::string> allowed={"EVENT","TREE","MIN_CHARGE","DBSCAN","DBSCAN_EPSILON_MM","DBSCAN_MIN_POINTS","CORRIDOR","CORRIDOR_HALF_WIDTH_FIBRES","MIN_MAP_FRACTION","FIT_VIEWS","SEED_DIRECTION","SEED_MEDIAN_FACTOR","MAX_FUNCTION_CALLS","TOLERANCE","FIT_RANGE","FIT_RANGE_QUANTILE","FIT_RANGE_PADDING_MM"};
#ifdef GLOBAL_FIT_COLUMN_MODE
        allowed.insert({"COLUMN_MIP_PRIOR","COLUMN_MIP_FLOOR_ENERGY","COLUMN_MIP_MPV_ENERGY","COLUMN_MIP_STRENGTH","COLUMN_CHARGE_PER_ENERGY","COLUMN_MIP_MIN_PATH_FRACTION","COLUMN_MIP_MAX_ITERATIONS","COLUMN_MIP_CONVERGENCE"});
#endif
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
        std::string seedMethod=get("SEED_DIRECTION","FIBRE");for(char& c:seedMethod)c=std::toupper(static_cast<unsigned char>(c));
        if(seedMethod!="FIBRE"&&seedMethod!="VIEW_MEDIAN_DIAMETER"&&seedMethod!="MC_SEGMENT_10")throw std::runtime_error("SEED_DIRECTION must be FIBRE, VIEW_MEDIAN_DIAMETER, or MC_SEGMENT_10");
        const double seedMedianFactor=std::stod(get("SEED_MEDIAN_FACTOR","1"));if(seedMedianFactor<=0)throw std::runtime_error("SEED_MEDIAN_FACTOR must be positive");
        const int maximumFunctionCalls=std::stoi(get("MAX_FUNCTION_CALLS","5000"));
        const double minimizerTolerance=std::stod(get("TOLERANCE","1e-2"));
        const bool useFitRange=boolean("FIT_RANGE",true);
        const double fitRangeQuantile=std::stod(get("FIT_RANGE_QUANTILE","0.01"));
        const double fitRangePadding=std::stod(get("FIT_RANGE_PADDING_MM","10"));
#ifdef GLOBAL_FIT_COLUMN_MODE
        ColumnMipPriorConfig columnPrior;columnPrior.enabled=boolean("COLUMN_MIP_PRIOR",false);
        columnPrior.floorEnergy=std::stod(get("COLUMN_MIP_FLOOR_ENERGY","0.5"));columnPrior.mpvEnergy=std::stod(get("COLUMN_MIP_MPV_ENERGY","1.5"));
        columnPrior.strength=std::stod(get("COLUMN_MIP_STRENGTH","0"));columnPrior.chargePerEnergy=std::stod(get("COLUMN_CHARGE_PER_ENERGY","1"));
        columnPrior.minPathFraction=std::stod(get("COLUMN_MIP_MIN_PATH_FRACTION","0.8"));columnPrior.maximumIterations=std::stoi(get("COLUMN_MIP_MAX_ITERATIONS","500"));
        columnPrior.convergence=std::stod(get("COLUMN_MIP_CONVERGENCE","1e-6"));
        if(columnPrior.floorEnergy<0||columnPrior.mpvEnergy<=columnPrior.floorEnergy||columnPrior.strength<0||columnPrior.chargePerEnergy<=0||columnPrior.minPathFraction<0||columnPrior.minPathFraction>1||columnPrior.maximumIterations<1||columnPrior.convergence<=0)throw std::runtime_error("Invalid COLUMN_MIP_* parameters");
        if(columnPrior.enabled&&columnPrior.strength<=0)throw std::runtime_error("COLUMN_MIP_STRENGTH must be positive when COLUMN_MIP_PRIOR=1");
#endif
        const auto fitViewSelection=ParseFitViews(get("FIT_VIEWS","ALL"));const auto& fitViews=fitViewSelection.first;const std::string& fitViewsName=fitViewSelection.second;
        if(dbscanEpsilon<=0||dbscanMinPoints<1||corridorHalfWidth<0||minimumMapFraction<0||maximumFunctionCalls<1||minimizerTolerance<=0||fitRangeQuantile<0||fitRangeQuantile>=0.5||fitRangePadding<0)throw std::runtime_error("Invalid selection, fit-range, or minimizer parameters");
        if(inName==outName)throw std::runtime_error("Input and output files must differ");
        TFile input(inName.c_str(),"READ");if(input.IsZombie())throw std::runtime_error("Cannot open input");
        auto* tree=dynamic_cast<TTree*>(input.Get(treeName.c_str()));if(!tree)throw std::runtime_error("Missing input tree");
        std::vector<int> events;
        if(eventSpec=="all") { int value=0;tree->SetBranchAddress("event",&value);std::set<int> unique;
            for(Long64_t i=0;i<tree->GetEntries();++i){tree->GetEntry(i);unique.insert(value);}events.assign(unique.begin(),unique.end()); }
        else events.push_back(std::stoi(eventSpec));
        std::array<double,3> low,high; const char* axis[3]={"x","y","z"};
        for(int k=0;k<3;++k){low[k]=tree->GetMinimum(axis[k])-5;high[k]=tree->GetMaximum(axis[k])+5;}
        LightMap map(mapName,minimumMapFraction);TFile output(outName.c_str(),"RECREATE");CloneFlatFile(input,output);output.cd();
        auto* mcSegments=dynamic_cast<TTree*>(input.Get("mc_virtual_segments"));auto* mcTracks=dynamic_cast<TTree*>(input.Get("mc_track_points"));
        if(!mcSegments||!mcTracks)throw std::runtime_error("Missing mc_virtual_segments or mc_track_points for MC line benchmark");
        const auto mcSeeds=seedMethod=="MC_SEGMENT_10"?McSegment10Seeds(*mcSegments,*mcTracks,events):std::map<int,std::array<double,3>>{};
        WriteMcSegmentBenchmark(output,*mcSegments,*mcTracks,events);
        GlobalFitResult r;TTree result("global_fit","Global light-pattern straight-track fit");BranchGlobalFit(result,r);
        int selectedEvent=0,selectedProjection=0,selectedU=0,selectedV=0;unsigned int selectedGeomId=0;
        double selectedX=0,selectedY=0,selectedZ=0,selectedTime=0,selectedCharge=0;
        TTree selected("global_fit_fibres","Fibres selected for the global fit");selected.Branch("event",&selectedEvent);
        selected.Branch("x",&selectedX);selected.Branch("y",&selectedY);selected.Branch("z",&selectedZ);selected.Branch("time",&selectedTime);
        selected.Branch("charge",&selectedCharge);selected.Branch("geom_id",&selectedGeomId);selected.Branch("projection",&selectedProjection);
        selected.Branch("u",&selectedU);selected.Branch("v",&selectedV);
        TTree dbscanSelected("global_fit_dbscan_fibres","Fibres retained after DBSCAN, before the optional corridor");
        dbscanSelected.Branch("event",&selectedEvent);dbscanSelected.Branch("x",&selectedX);dbscanSelected.Branch("y",&selectedY);
        dbscanSelected.Branch("z",&selectedZ);dbscanSelected.Branch("time",&selectedTime);dbscanSelected.Branch("charge",&selectedCharge);
        dbscanSelected.Branch("geom_id",&selectedGeomId);dbscanSelected.Branch("projection",&selectedProjection);
        dbscanSelected.Branch("u",&selectedU);dbscanSelected.Branch("v",&selectedV);
#ifdef GLOBAL_FIT_COLUMN_MODE
        int columnEvent=0,columnAxis=0,columnIndex=0;double columnLow=0,columnHigh=0,columnPathLength=0;
        double columnCharge=0,columnChargePerMm=0,columnResponseScale=0,columnUnregularizedCharge=0,columnUnregularizedScale=0,columnMipFloorCharge=0,columnMipMpvCharge=0;int columnMipEligible=0;
        TTree columnCharges("global_fit_column_charges","Non-negative column charges fitted after fixing the track geometry");
        columnCharges.Branch("event",&columnEvent);columnCharges.Branch("axis",&columnAxis);columnCharges.Branch("column_index",&columnIndex);
        columnCharges.Branch("column_low",&columnLow);columnCharges.Branch("column_high",&columnHigh);columnCharges.Branch("path_length",&columnPathLength);
        columnCharges.Branch("fitted_selected_fibre_charge",&columnCharge);columnCharges.Branch("charge_per_mm",&columnChargePerMm);
        columnCharges.Branch("response_scale",&columnResponseScale);
        columnCharges.Branch("unregularized_fitted_selected_fibre_charge",&columnUnregularizedCharge);
        columnCharges.Branch("unregularized_response_scale",&columnUnregularizedScale);
        columnCharges.Branch("mip_prior_eligible",&columnMipEligible);columnCharges.Branch("mip_prior_floor_charge",&columnMipFloorCharge);columnCharges.Branch("mip_prior_mpv_charge",&columnMipMpvCharge);
#endif
        for(int event:events){
            auto obs=SelectFitViews(Read(*tree,event,cut),fitViews);const int observationsBefore=obs.size();if(useDbscan)obs=MainDbscanClusters(obs,dbscanEpsilon,dbscanMinPoints);
            const int observationsAfterDbscan=obs.size();if(obs.empty()){std::cerr<<"Skipping event "<<event<<": no observations survive DBSCAN\n";continue;}
            for(const auto& hit:obs){selectedEvent=event;selectedX=hit.x;selectedY=hit.y;selectedZ=hit.z;selectedTime=hit.time;selectedCharge=hit.q;
                selectedGeomId=hit.geomId;selectedProjection=hit.projection;selectedU=hit.u;selectedV=hit.v;dbscanSelected.Fill();}
            const auto fitBounds=useFitRange?DataFitBounds(obs,low,high,fitRangeQuantile,fitRangePadding):std::make_pair(low,high);
            const auto& fitLow=fitBounds.first;const auto& fitHigh=fitBounds.second;
            const TVector3 fitCentre(.5*(fitLow[0]+fitHigh[0]),.5*(fitLow[1]+fitHigh[1]),.5*(fitLow[2]+fitHigh[2]));
            MedianDiameterSeedResult medianSeed;auto seed=Seed(obs,fitCentre);
            if(seedMethod=="VIEW_MEDIAN_DIAMETER"){medianSeed=MedianDiameterSeed(obs,seedMedianFactor);seed=medianSeed.direction;}
            else if(seedMethod=="MC_SEGMENT_10"){auto truth=mcSeeds.find(event);if(truth==mcSeeds.end()){std::cerr<<"Skipping event "<<event<<": no MC segment-10 seed\n";continue;}seed=truth->second;}
            if(useCorridor)obs=LineCorridor(obs,seed,fitCentre,corridorHalfWidth);
            if(obs.empty()){std::cerr<<"Skipping event "<<event<<": no observations survive corridor\n";continue;}if(seedMethod=="FIBRE")seed=Seed(obs,fitCentre);
            for(const auto& hit:obs){selectedEvent=event;selectedX=hit.x;selectedY=hit.y;selectedZ=hit.z;selectedTime=hit.time;selectedCharge=hit.q;
                selectedGeomId=hit.geomId;selectedProjection=hit.projection;selectedU=hit.u;selectedV=hit.v;selected.Fill();}
            Likelihood likelihood(map,obs,seed,fitLow,fitHigh,fitViews,low
#ifdef GLOBAL_FIT_COLUMN_MODE
                                  ,columnPrior
#endif
                                  );double theta=std::acos(std::clamp(seed[2],-1.0,1.0)),phi=std::atan2(seed[1],seed[0]);
            const double seedParameters[4]={0,0,theta,phi};const double seedNll=likelihood(seedParameters);int seedNdof=0;const double seedChi2=likelihood.Chi2(seedParameters,seedNdof);
            auto minimizer=std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex"));
            ROOT::Math::Functor functor(likelihood,4);minimizer->SetFunction(functor);minimizer->SetMaxFunctionCalls(maximumFunctionCalls);minimizer->SetTolerance(minimizerTolerance);
            minimizer->SetLimitedVariable(0,"impact_a",0,5,-500,500);minimizer->SetLimitedVariable(1,"impact_b",0,5,-500,500);
            minimizer->SetLimitedVariable(2,"theta",theta,.01,0.001,3.1405);minimizer->SetVariable(3,"phi",phi,.01);minimizer->Minimize();
            const double* p=minimizer->X();TVector3 point=likelihood.Point(p[0],p[1]);TVector3 dir(std::sin(p[2])*std::cos(p[3]),std::sin(p[2])*std::sin(p[3]),std::cos(p[2]));
#ifdef GLOBAL_FIT_COLUMN_MODE
            const auto columnFit=likelihood.FitColumnCharges(p);
            for(const auto& column:columnFit.columns){columnEvent=event;columnAxis=column.axis;columnIndex=column.index;
                columnLow=column.low;columnHigh=column.high;columnPathLength=column.pathLength;columnCharge=column.fittedCharge;
                columnChargePerMm=column.chargePerMm;columnResponseScale=column.responseScale;columnUnregularizedCharge=column.unregularizedFittedCharge;
                columnUnregularizedScale=column.unregularizedResponseScale;columnMipEligible=column.mipPriorEligible;columnMipFloorCharge=column.mipFloorCharge;columnMipMpvCharge=column.mipMpvCharge;columnCharges.Fill();}
#endif
            r=GlobalFitResult();r.event=event;r.status=minimizer->Status();r.observations=obs.size();r.observationsBeforeClustering=observationsBefore;r.observationsAfterDbscan=observationsAfterDbscan;
            r.dbscanEnabled=useDbscan;r.dbscanEpsilon=dbscanEpsilon;r.dbscanMinPoints=dbscanMinPoints;r.minimumCharge=cut;r.nll=minimizer->MinValue();
            r.fitRangeEnabled=useFitRange;r.fitRangeQuantile=fitRangeQuantile;r.fitRangePadding=fitRangePadding;
            r.fitLowX=fitLow[0];r.fitLowY=fitLow[1];r.fitLowZ=fitLow[2];r.fitHighX=fitHigh[0];r.fitHighY=fitHigh[1];r.fitHighZ=fitHigh[2];
            r.edm=minimizer->Edm();r.functionCalls=minimizer->NCalls();r.maximumFunctionCalls=maximumFunctionCalls;r.tolerance=minimizerTolerance;
            r.chi2=likelihood.Chi2(p,r.ndof);r.inputTree=treeName;
            r.corridorEnabled=useCorridor;r.corridorHalfWidth=corridorHalfWidth;
            r.minimumMapFraction=minimumMapFraction;
            r.fitViews=fitViewsName;
            r.seedMethod=seedMethod;
            r.seedObservations=seedMethod=="VIEW_MEDIAN_DIAMETER"?medianSeed.observations:obs.size();r.seedMedianFactor=seedMedianFactor;
            r.seedMedianXZ=medianSeed.median[0];r.seedMedianYZ=medianSeed.median[1];r.seedMedianXY=medianSeed.median[2];r.seedNll=seedNll;r.seedChi2=seedChi2;r.seedNdof=seedNdof;
            r.seedX=fitCentre.X();r.seedY=fitCentre.Y();r.seedZ=fitCentre.Z();r.seedDx=seed[0];r.seedDy=seed[1];r.seedDz=seed[2];r.fitX=point.X();r.fitY=point.Y();r.fitZ=point.Z();r.fitDx=dir.X();r.fitDy=dir.Y();r.fitDz=dir.Z();
#ifdef GLOBAL_FIT_COLUMN_MODE
            r.columnFitAvailable=1;r.columnFitConverged=columnFit.converged;r.columnFitIterations=columnFit.iterations;
            r.columnFitActiveColumns=columnFit.activeColumns;r.columnFitMaxRelativeChange=columnFit.maxRelativeChange;
            r.columnNll=columnFit.nll;r.columnChi2=columnFit.chi2;r.columnNdof=columnFit.ndof;
            r.columnMipPriorEnabled=columnFit.mipPriorEnabled;r.columnMipPriorConverged=columnFit.mipPriorConverged;r.columnMipPriorStalled=columnFit.mipPriorStalled;r.columnMipPriorIterations=columnFit.mipPriorIterations;
            r.columnMipPriorEligibleColumns=columnFit.mipPriorEligibleColumns;r.columnMipPriorFloorEnergy=columnPrior.floorEnergy;r.columnMipPriorMpvEnergy=columnPrior.mpvEnergy;
            r.columnMipPriorStrength=columnPrior.strength;r.columnMipPriorChargePerEnergy=columnPrior.chargePerEnergy;r.columnMipPriorMinPathFraction=columnPrior.minPathFraction;
            r.columnMipPriorPenalty=columnFit.mipPriorPenalty;r.columnMipPriorObjective=columnFit.mipPriorObjective;r.columnMipPriorMaxRelativeChange=columnFit.mipPriorMaxRelativeChange;
            r.columnUnregularizedNll=columnFit.unregularizedNll;
#endif
            r.errorA=minimizer->Errors()[0];r.errorB=minimizer->Errors()[1];r.errorTheta=minimizer->Errors()[2];r.errorPhi=minimizer->Errors()[3];result.Fill();
            std::cout<<"event "<<event<<" views "<<r.fitViews<<" fibres "<<observationsBefore<<" -> "<<observationsAfterDbscan<<" -> "<<r.observations<<" status "<<r.status<<" nll "<<r.nll<<" edm "<<r.edm<<" calls "<<r.functionCalls<<"/"<<r.maximumFunctionCalls<<" chi2/ndof "<<r.chi2<<"/"<<r.ndof<<" = "<<(r.ndof>0?r.chi2/r.ndof:0)<<" direction "<<r.fitDx<<" "<<r.fitDy<<" "<<r.fitDz<<"\n";
#ifdef GLOBAL_FIT_COLUMN_MODE
            std::cout<<"  column fit converged "<<r.columnFitConverged<<" iterations "<<r.columnFitIterations<<" max relative change "<<r.columnFitMaxRelativeChange
                     <<" active columns "<<r.columnFitActiveColumns<<" nll "<<r.columnNll<<" chi2/ndof "<<r.columnChi2<<"/"<<r.columnNdof
                     <<" = "<<(r.columnNdof>0?r.columnChi2/r.columnNdof:0)<<"\n";
            if(r.columnMipPriorEnabled)std::cout<<"  MIP prior converged "<<r.columnMipPriorConverged<<" stalled "<<r.columnMipPriorStalled<<" iterations "<<r.columnMipPriorIterations
                     <<" eligible columns "<<r.columnMipPriorEligibleColumns<<" penalty "<<r.columnMipPriorPenalty<<" objective "<<r.columnMipPriorObjective<<"\n";
#endif
        }
        output.cd();result.Write();selected.Write();dbscanSelected.Write();
#ifdef GLOBAL_FIT_COLUMN_MODE
        columnCharges.Write();
#endif
        output.Close();std::cout<<"output "<<outName<<"\n";
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}return 0;
}
#endif
