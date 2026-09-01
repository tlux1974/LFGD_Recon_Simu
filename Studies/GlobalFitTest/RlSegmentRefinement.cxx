// Offline third-stage refinement of GlobalLightFitColumns results.
// The light-map response is the Richardson-Lucy forward operator.  A
// sample-wide reco--MC residual density supplies the displacement prior and a
// discrete second-difference prior suppresses column-to-column zigzags.
#define GLOBAL_LIGHT_FIT_NO_MAIN
#include "GlobalLightFit.cxx"

#include <TParameter.h>
#include <TTree.h>

#include <iomanip>
#include <numeric>
#include <sstream>

namespace {

struct ColumnInput { int axis=0,index=0;double low=0,high=0,path=0,visible=0,scale=0; };
struct EventInput { GlobalFitResult fit;std::vector<ColumnInput> columns;std::vector<Observation> observations; };
struct Candidate {
    int column=0;double a=0,b=0,amplitude=0,initial=0,norm=0,prior=0;
    std::vector<std::pair<int,double>> response;
};
struct IterationRecord { int event=0,iteration=0;double nll=0,maxChange=0,curvatureRms=0; };

bool Boolean(std::string value) {
    for(char& c:value)c=std::tolower(static_cast<unsigned char>(c));
    if(value=="1"||value=="true"||value=="yes"||value=="on")return true;
    if(value=="0"||value=="false"||value=="no"||value=="off")return false;
    throw std::runtime_error("Expected boolean value, got "+value);
}

std::pair<TVector3,TVector3> Basis(const TVector3& direction) {
    const TVector3 d=direction.Unit();const TVector3 trial=std::abs(d.Z())<.9?TVector3(0,0,1):TVector3(0,1,0);
    const TVector3 e1=d.Cross(trial).Unit();return {e1,d.Cross(e1).Unit()};
}

std::map<int,EventInput> ReadEvents(TFile& input) {
    auto* fitTree=dynamic_cast<TTree*>(input.Get("global_fit"));
    auto* columnTree=dynamic_cast<TTree*>(input.Get("global_fit_column_charges"));
    auto* fibreTree=dynamic_cast<TTree*>(input.Get("global_fit_fibres"));
    if(!fitTree||!columnTree||!fibreTree)throw std::runtime_error("Input needs global_fit, global_fit_column_charges, and global_fit_fibres");
    std::map<int,EventInput> events;GlobalFitResult fit;SetGlobalFitAddresses(*fitTree,fit);
    for(Long64_t i=0;i<fitTree->GetEntries();++i){fitTree->GetEntry(i);events[fit.event].fit=fit;}
    fitTree->ResetBranchAddresses();
    int event=0,axis=0,index=0;double low=0,high=0,path=0,visible=0,scale=0;
    columnTree->SetBranchAddress("event",&event);columnTree->SetBranchAddress("axis",&axis);columnTree->SetBranchAddress("column_index",&index);
    columnTree->SetBranchAddress("column_low",&low);columnTree->SetBranchAddress("column_high",&high);columnTree->SetBranchAddress("path_length",&path);
    columnTree->SetBranchAddress("fitted_selected_fibre_charge",&visible);columnTree->SetBranchAddress("response_scale",&scale);
    for(Long64_t i=0;i<columnTree->GetEntries();++i){columnTree->GetEntry(i);if(events.count(event))events[event].columns.push_back({axis,index,low,high,path,visible,scale});}
    columnTree->ResetBranchAddresses();
    int projection=0,u=0,v=0;unsigned int geomId=0;double x=0,y=0,z=0,time=0,charge=0;
    fibreTree->SetBranchAddress("event",&event);fibreTree->SetBranchAddress("projection",&projection);fibreTree->SetBranchAddress("u",&u);fibreTree->SetBranchAddress("v",&v);
    fibreTree->SetBranchAddress("geom_id",&geomId);fibreTree->SetBranchAddress("x",&x);fibreTree->SetBranchAddress("y",&y);fibreTree->SetBranchAddress("z",&z);
    fibreTree->SetBranchAddress("time",&time);fibreTree->SetBranchAddress("charge",&charge);
    for(Long64_t i=0;i<fibreTree->GetEntries();++i){fibreTree->GetEntry(i);if(events.count(event))events[event].observations.push_back({x,y,z,time,charge,geomId,projection,u,v});}
    fibreTree->ResetBranchAddresses();
    for(auto& item:events)std::sort(item.second.columns.begin(),item.second.columns.end(),[](const auto& a,const auto& b){return a.index<b.index;});
    return events;
}

std::vector<std::array<double,2>> ResidualSample(TFile& input,const std::map<int,EventInput>& events,double maximumDistance) {
    auto* tree=dynamic_cast<TTree*>(input.Get("mc_virtual_segments"));if(!tree)throw std::runtime_error("Missing mc_virtual_segments");
    if(!tree->GetBranch("contributor_pdgs"))throw std::runtime_error("mc_virtual_segments lacks contributor_pdgs; regenerate the flat tree and column fit");
    int event=0,detector=0,cubeX=0,cubeY=0,cubeZ=0;double x=0,y=0,z=0;std::vector<int>* pdgs=nullptr;
    tree->SetBranchAddress("event",&event);tree->SetBranchAddress("detector",&detector);tree->SetBranchAddress("cube_x",&cubeX);tree->SetBranchAddress("cube_y",&cubeY);tree->SetBranchAddress("cube_z",&cubeZ);
    tree->SetBranchAddress("start_x",&x);tree->SetBranchAddress("start_y",&y);tree->SetBranchAddress("start_z",&z);tree->SetBranchAddress("contributor_pdgs",&pdgs);
    std::map<std::tuple<int,int,int,int>,std::pair<double,std::array<double,2>>> closest;
    for(Long64_t i=0;i<tree->GetEntries();++i){tree->GetEntry(i);auto found=events.find(event);if(detector!=0||found==events.end()||!pdgs)continue;
        bool muon=false;for(int pdg:*pdgs)if(std::abs(pdg)==13){muon=true;break;}if(!muon)continue;
        const auto& f=found->second.fit;TVector3 point(f.fitX,f.fitY,f.fitZ),direction(f.fitDx,f.fitDy,f.fitDz),truth(x,y,z);auto basis=Basis(direction);
        const TVector3 delta=truth-point-direction.Unit()*((truth-point).Dot(direction.Unit()));const double distance=delta.Mag();if(distance>maximumDistance)continue;
        const auto key=std::make_tuple(event,cubeX,cubeY,cubeZ);std::array<double,2> offset{{delta.Dot(basis.first),delta.Dot(basis.second)}};
        auto old=closest.find(key);if(old==closest.end()||distance<old->second.first)closest[key]={distance,offset};
    }
    tree->ResetBranchAddresses();std::vector<std::array<double,2>> result;for(const auto& item:closest)result.push_back(item.second.second);return result;
}

double PriorDensity(double a,double b,const std::vector<std::array<double,2>>& sample,double bandwidth) {
    double sum=0;const double inverse=1.0/(2*bandwidth*bandwidth);
    for(const auto& value:sample){const double da=a-value[0],db=b-value[1];sum+=std::exp(-(da*da+db*db)*inverse);}
    return std::max(sum/std::max<size_t>(sample.size(),1),1e-15);
}

bool SegmentResponse(const LightMap& map,const GlobalFitResult& fit,const ColumnInput& column,const TVector3& point,
                     const TVector3& direction,double step,const std::map<FibreKeyType,int>& observed,
                     std::vector<std::pair<int,double>>& response,double& norm) {
    std::array<double,3> low{{fit.fitLowX,fit.fitLowY,fit.fitLowZ}},high{{fit.fitHighX,fit.fitHighY,fit.fitHighZ}};
    low[column.axis]=std::max(low[column.axis],column.low);high[column.axis]=std::min(high[column.axis],column.high);
    double begin=-1e30,end=1e30;
    for(int k=0;k<3;++k){const double d=direction[k];if(std::abs(d)<1e-9){if(point[k]<low[k]||point[k]>high[k])return false;continue;}
        double t1=(low[k]-point[k])/d,t2=(high[k]-point[k])/d;if(t1>t2)std::swap(t1,t2);begin=std::max(begin,t1);end=std::min(end,t2);}
    if(begin>=end)return false;
    const int steps=std::max(1,static_cast<int>(std::ceil((end-begin)/step)));const double ds=(end-begin)/steps;
    std::map<int,double> accumulated;
    for(int s=0;s<steps;++s){const TVector3 position=point+(begin+(s+.5)*ds)*direction,zero;TVector3 offset;
        for(const auto& hit:map.At(position,offset)){const TVector3 fibre=hit.position+offset;
            if(fibre.X()<fit.fitLowX||fibre.X()>fit.fitHighX||fibre.Y()<fit.fitLowY||fibre.Y()>fit.fitHighY||fibre.Z()<fit.fitLowZ||fibre.Z()>fit.fitHighZ)continue;
            const int projection=std::abs(hit.position.X())<.01?1:std::abs(hit.position.Y())<.01?0:2;auto key=observed.find(FibreKey(projection,fibre));
            if(key!=observed.end())accumulated[key->second]+=hit.fraction*ds;}}
    norm=0;for(const auto& value:accumulated){response.push_back(value);norm+=value.second;}return norm>0;
}

double Nll(const std::vector<double>& observed,const std::vector<double>& expected) {
    double value=0;for(size_t i=0;i<observed.size();++i)value+=expected[i]-observed[i]*std::log(std::max(expected[i],1e-15));return value;
}

double CurvatureRms(const std::vector<std::array<double,2>>& centroid) {
    if(centroid.size()<3)return 0;
    double sum=0;int count=0;for(size_t k=1;k+1<centroid.size();++k)for(int d=0;d<2;++d){const double c=centroid[k+1][d]-2*centroid[k][d]+centroid[k-1][d];sum+=c*c;++count;}return std::sqrt(sum/std::max(count,1));
}

}

int main(int argc,char** argv) {
    auto usage=[](){std::cout<<R"(Usage:
  rl_segment_refinement INPUT_COLUMN_FIT.root OUTPUT.root LIGHTMAP.root [KEY=value ...]

Options:
  EVENT=all|N                   Events to refine (default: all)
  GRID_RADIUS_MM=3              Transverse candidate-grid radius
  GRID_STEP_MM=1               Transverse candidate-grid spacing
  ITERATIONS=50                 Maximum RL iterations
  CONVERGENCE=1e-4              Maximum relative amplitude-change criterion
  POSITION_CONVERGENCE_MM=1e-3  Maximum centroid displacement criterion
  OBJECTIVE_CONVERGENCE=1e-7    Relative penalized-objective criterion
  PRIOR_BANDWIDTH_MM=0.5        KDE bandwidth for aggregate reco-MC residual PDF
  PRIOR_MAX_RESIDUAL_MM=6       Largest residual admitted to the aggregate PDF
  PRIOR_STRENGTH=0.02           Per-iteration displacement-prior exponent
  CURVATURE_SIGMA_MM=0.75       Scale of the discrete-curvature prior
  CURVATURE_STRENGTH=0.20       Per-iteration curvature-prior exponent
  RESPONSE_STEP_MM=1            Light-map integration step (default: map pitch)
  MIN_MAP_FRACTION=0            Ignore smaller light-map fractions
  SAVE_CANDIDATES=1             Write final candidate amplitudes
)";};
    if(argc<4){usage();return argc==1?0:2;}
    try {
        std::map<std::string,std::string> options;for(int i=4;i<argc;++i){std::string value=argv[i];const auto split=value.find('=');if(split==std::string::npos)throw std::runtime_error("Option must use KEY=value: "+value);options[value.substr(0,split)]=value.substr(split+1);}
        auto get=[&](const std::string& key,const std::string& fallback){auto found=options.find(key);return found==options.end()?fallback:found->second;};
        const std::string eventSpec=get("EVENT","all");const double radius=std::stod(get("GRID_RADIUS_MM","3")),gridStep=std::stod(get("GRID_STEP_MM","1"));
        const int maximumIterations=std::stoi(get("ITERATIONS","50"));const double convergence=std::stod(get("CONVERGENCE","1e-4"));
        const double positionConvergence=std::stod(get("POSITION_CONVERGENCE_MM","1e-3")),objectiveConvergence=std::stod(get("OBJECTIVE_CONVERGENCE","1e-7"));
        const double priorBandwidth=std::stod(get("PRIOR_BANDWIDTH_MM","0.5")),priorMaximum=std::stod(get("PRIOR_MAX_RESIDUAL_MM","6"));
        const double priorStrength=std::stod(get("PRIOR_STRENGTH","0.02")),curvatureSigma=std::stod(get("CURVATURE_SIGMA_MM","0.75"));
        const double curvatureStrength=std::stod(get("CURVATURE_STRENGTH","0.20")),responseStep=std::stod(get("RESPONSE_STEP_MM","1"));
        const double minimumMapFraction=std::stod(get("MIN_MAP_FRACTION","0"));const bool saveCandidates=Boolean(get("SAVE_CANDIDATES","1"));
        if(radius<0||gridStep<=0||maximumIterations<1||convergence<=0||positionConvergence<=0||objectiveConvergence<=0||priorBandwidth<=0||priorMaximum<=0||curvatureSigma<=0||responseStep<=0||priorStrength<0||curvatureStrength<0)throw std::runtime_error("Invalid refinement option");
        TFile input(argv[1],"READ");if(input.IsZombie())throw std::runtime_error("Cannot open input");auto events=ReadEvents(input);
        const auto residualSample=ResidualSample(input,events,priorMaximum);if(residualSample.size()<10)throw std::runtime_error("Too few primary-muon residuals for the aggregate prior");
        LightMap map(argv[3],minimumMapFraction);TFile output(argv[2],"RECREATE");if(output.IsZombie())throw std::runtime_error("Cannot create output");
        int outEvent=0,outConverged=0,outStalled=0,outIterations=0,outBestNllIteration=0,outBestObjectiveIteration=0,outObservations=0,outColumns=0;double initialNll=0,finalNll=0,outBestNll=0,initialObjective=0,outMaxChange=0,outCurvature=0,outObjective=0,outPriorPenalty=0,outCurvaturePenalty=0;
        TTree fitOutput("rl_fit","Offline Richardson-Lucy segment refinement summary");fitOutput.Branch("event",&outEvent);fitOutput.Branch("converged",&outConverged);fitOutput.Branch("iterations",&outIterations);
        fitOutput.Branch("observations",&outObservations);fitOutput.Branch("columns",&outColumns);fitOutput.Branch("initial_nll",&initialNll);fitOutput.Branch("final_nll",&finalNll);
        fitOutput.Branch("best_nll",&outBestNll);fitOutput.Branch("best_nll_iteration",&outBestNllIteration);fitOutput.Branch("best_objective_iteration",&outBestObjectiveIteration);
        fitOutput.Branch("stalled",&outStalled);fitOutput.Branch("max_relative_change",&outMaxChange);fitOutput.Branch("curvature_rms_mm",&outCurvature);
        fitOutput.Branch("initial_penalized_objective",&initialObjective);fitOutput.Branch("penalized_objective",&outObjective);fitOutput.Branch("displacement_prior_penalty",&outPriorPenalty);fitOutput.Branch("curvature_penalty",&outCurvaturePenalty);
        int cEvent=0,cAxis=0,cIndex=0;double cLow=0,cHigh=0,cA=0,cB=0,cX=0,cY=0,cZ=0,cInitial=0,cScale=0,cVisible=0,cCurvA=0,cCurvB=0,cBestA=0,cBestB=0,cBestX=0,cBestY=0,cBestZ=0,cBestScale=0,cBestVisible=0;
        TTree columnOutput("rl_columns","Final RL-refined segment positions and charges");columnOutput.Branch("event",&cEvent);columnOutput.Branch("axis",&cAxis);columnOutput.Branch("column_index",&cIndex);
        columnOutput.Branch("column_low",&cLow);columnOutput.Branch("column_high",&cHigh);columnOutput.Branch("offset_a_mm",&cA);columnOutput.Branch("offset_b_mm",&cB);
        columnOutput.Branch("centre_x",&cX);columnOutput.Branch("centre_y",&cY);columnOutput.Branch("centre_z",&cZ);columnOutput.Branch("initial_response_scale",&cInitial);
        columnOutput.Branch("refined_response_scale",&cScale);columnOutput.Branch("predicted_selected_fibre_charge",&cVisible);columnOutput.Branch("curvature_a_mm",&cCurvA);columnOutput.Branch("curvature_b_mm",&cCurvB);
        columnOutput.Branch("best_nll_offset_a_mm",&cBestA);columnOutput.Branch("best_nll_offset_b_mm",&cBestB);columnOutput.Branch("best_nll_centre_x",&cBestX);columnOutput.Branch("best_nll_centre_y",&cBestY);columnOutput.Branch("best_nll_centre_z",&cBestZ);
        columnOutput.Branch("best_nll_response_scale",&cBestScale);columnOutput.Branch("best_nll_predicted_selected_fibre_charge",&cBestVisible);
        int iEvent=0,iIteration=0,iAccepted=0;double iNll=0,iNllBefore=0,iNllAfterRl=0,iChange=0,iPositionChange=0,iCurvature=0,iPriorPenalty=0,iCurvaturePenalty=0,iObjective=0,iStep=0;
        TTree iterationOutput("rl_iterations","RL convergence and penalized-objective history");iterationOutput.Branch("event",&iEvent);iterationOutput.Branch("iteration",&iIteration);iterationOutput.Branch("accepted",&iAccepted);
        iterationOutput.Branch("nll",&iNll);iterationOutput.Branch("nll_before",&iNllBefore);iterationOutput.Branch("nll_after_rl",&iNllAfterRl);iterationOutput.Branch("max_relative_change",&iChange);
        iterationOutput.Branch("max_position_change_mm",&iPositionChange);iterationOutput.Branch("curvature_rms_mm",&iCurvature);iterationOutput.Branch("displacement_prior_penalty",&iPriorPenalty);
        iterationOutput.Branch("curvature_penalty",&iCurvaturePenalty);iterationOutput.Branch("penalized_objective",&iObjective);iterationOutput.Branch("accepted_step_fraction",&iStep);
        int pEvent=0,pColumn=0;double pA=0,pB=0,pAmplitude=0,pBestNllAmplitude=0,pPrior=0,pNorm=0;TTree candidateOutput("rl_candidates","Final transverse-grid amplitudes");candidateOutput.Branch("event",&pEvent);candidateOutput.Branch("column_index",&pColumn);
        candidateOutput.Branch("offset_a_mm",&pA);candidateOutput.Branch("offset_b_mm",&pB);candidateOutput.Branch("response_scale",&pAmplitude);candidateOutput.Branch("best_nll_response_scale",&pBestNllAmplitude);candidateOutput.Branch("prior_density",&pPrior);candidateOutput.Branch("response_normalization",&pNorm);

        for(auto& item:events){const int event=item.first;if(eventSpec!="all"&&event!=std::stoi(eventSpec))continue;auto& data=item.second;if(data.columns.empty()||data.observations.empty())continue;
            const TVector3 direction(data.fit.fitDx,data.fit.fitDy,data.fit.fitDz),linePoint(data.fit.fitX,data.fit.fitY,data.fit.fitZ);auto basis=Basis(direction);
            std::map<FibreKeyType,int> observationIndex;std::vector<double> observed;for(const auto& hit:data.observations){const auto key=FibreKey(hit.projection,TVector3(hit.x,hit.y,hit.z));auto found=observationIndex.find(key);
                if(found==observationIndex.end()){const int index=observed.size();observationIndex[key]=index;observed.push_back(hit.q);}else observed[found->second]+=hit.q;}
            std::vector<Candidate> candidates;std::vector<std::vector<int>> byColumn(data.columns.size());
            for(size_t k=0;k<data.columns.size();++k){const auto& column=data.columns[k];std::vector<double> priors;std::vector<int> local;
                for(double a=-radius;a<=radius+1e-9;a+=gridStep)for(double b=-radius;b<=radius+1e-9;b+=gridStep){Candidate candidate;candidate.column=k;candidate.a=a;candidate.b=b;candidate.prior=PriorDensity(a,b,residualSample,priorBandwidth);
                    const TVector3 shifted=linePoint+a*basis.first+b*basis.second;if(!SegmentResponse(map,data.fit,column,shifted,direction.Unit(),responseStep,observationIndex,candidate.response,candidate.norm))continue;
                    local.push_back(candidates.size());priors.push_back(candidate.prior);candidates.push_back(std::move(candidate));}
                const double priorSum=std::accumulate(priors.begin(),priors.end(),0.0);for(size_t j=0;j<local.size();++j){auto& candidate=candidates[local[j]];candidate.amplitude=column.scale*priors[j]/std::max(priorSum,1e-15);candidate.initial=candidate.amplitude;byColumn[k].push_back(local[j]);}}
            const double observedTotal=std::accumulate(observed.begin(),observed.end(),0.0);
            auto predictValues=[&](const std::vector<double>& amplitude){std::vector<double> prediction(observed.size(),0);for(size_t j=0;j<candidates.size();++j)for(const auto& response:candidates[j].response)prediction[response.first]+=amplitude[j]*response.second;return prediction;};
            auto centroidValues=[&](const std::vector<double>& amplitude){std::vector<std::array<double,2>> value(byColumn.size(),{{0,0}});for(size_t k=0;k<byColumn.size();++k){double total=0;for(int j:byColumn[k]){total+=amplitude[j];value[k][0]+=amplitude[j]*candidates[j].a;value[k][1]+=amplitude[j]*candidates[j].b;}if(total>0){value[k][0]/=total;value[k][1]/=total;}}return value;};
            auto penalties=[&](const std::vector<double>& amplitude,const std::vector<std::array<double,2>>& centre){double displacement=0;for(const auto& indices:byColumn){double total=0,maximumPrior=0;for(int j:indices){total+=amplitude[j];maximumPrior=std::max(maximumPrior,candidates[j].prior);}for(int j:indices)if(total>0)displacement+=(amplitude[j]/total)*(-std::log(std::max(candidates[j].prior/std::max(maximumPrior,1e-15),1e-15)));}displacement/=std::max<size_t>(byColumn.size(),1);
                double curvature=0;int terms=0;for(size_t k=1;k+1<centre.size();++k)for(int d=0;d<2;++d){const double second=centre[k+1][d]-2*centre[k][d]+centre[k-1][d];curvature+=second*second/(2*curvatureSigma*curvatureSigma);++terms;}curvature/=std::max(terms,1);return std::make_pair(displacement,curvature);};
            auto objective=[&](double nll,const std::pair<double,double>& penalty){return nll/std::max(observedTotal,1e-15)+priorStrength*penalty.first+curvatureStrength*penalty.second;};
            std::vector<double> amplitudes;for(const auto& candidate:candidates)amplitudes.push_back(candidate.amplitude);
            auto expected=predictValues(amplitudes);auto centroids=centroidValues(amplitudes);auto currentPenalty=penalties(amplitudes,centroids);double currentNll=Nll(observed,expected),currentObjective=objective(currentNll,currentPenalty);
            std::vector<double> bestNllAmplitudes=amplitudes;double bestNll=currentNll;int bestNllIteration=0,bestObjectiveIteration=0;
            initialNll=currentNll;initialObjective=currentObjective;outConverged=0;outStalled=0;outIterations=0;outMaxChange=0;
            for(int iteration=1;iteration<=maximumIterations;++iteration){const std::vector<double> previous=amplitudes;const auto previousCentroids=centroids;const double previousNll=currentNll,previousObjective=currentObjective;
                std::vector<double> proposal(candidates.size());for(size_t j=0;j<candidates.size();++j){double numerator=0;for(const auto& response:candidates[j].response)numerator+=response.second*observed[response.first]/std::max(expected[response.first],1e-15);proposal[j]=amplitudes[j]*numerator/std::max(candidates[j].norm,1e-15);}
                const auto rlExpected=predictValues(proposal);const double rlNll=Nll(observed,rlExpected);auto rlCentroids=centroidValues(proposal);
                for(size_t k=0;k<byColumn.size();++k){double before=0,after=0;for(int j:byColumn[k])before+=proposal[j];for(int j:byColumn[k]){double factor=std::pow(candidates[j].prior,priorStrength);
                        if(k>0&&k+1<byColumn.size()){const double targetA=.5*(rlCentroids[k-1][0]+rlCentroids[k+1][0]),targetB=.5*(rlCentroids[k-1][1]+rlCentroids[k+1][1]);const double da=candidates[j].a-targetA,db=candidates[j].b-targetB;factor*=std::exp(-curvatureStrength*(da*da+db*db)/(2*curvatureSigma*curvatureSigma));}
                        proposal[j]*=factor;after+=proposal[j];}if(after>0)for(int j:byColumn[k])proposal[j]*=before/after;}
                double acceptedStep=1;bool accepted=false;std::vector<double> trial;std::vector<double> trialExpected;std::vector<std::array<double,2>> trialCentroids;std::pair<double,double> trialPenalty;double trialNll=0,trialObjective=0;
                for(int lineSearch=0;lineSearch<20;++lineSearch){trial.resize(amplitudes.size());for(size_t j=0;j<trial.size();++j)trial[j]=std::max(previous[j]+acceptedStep*(proposal[j]-previous[j]),1e-15);trialExpected=predictValues(trial);trialNll=Nll(observed,trialExpected);trialCentroids=centroidValues(trial);trialPenalty=penalties(trial,trialCentroids);trialObjective=objective(trialNll,trialPenalty);
                    if(trialObjective<=previousObjective+1e-12){accepted=true;break;}acceptedStep*=.5;}
                double maximumChange=0,maximumPositionChange=0;if(accepted){for(const auto& indices:byColumn){double change=0,normalization=0;for(int j:indices){change+=std::abs(trial[j]-previous[j]);normalization+=previous[j];}maximumChange=std::max(maximumChange,change/std::max(normalization,1e-12));}
                    for(size_t k=0;k<trialCentroids.size();++k)maximumPositionChange=std::max(maximumPositionChange,std::hypot(trialCentroids[k][0]-previousCentroids[k][0],trialCentroids[k][1]-previousCentroids[k][1]));
                    amplitudes.swap(trial);expected.swap(trialExpected);centroids.swap(trialCentroids);currentPenalty=trialPenalty;currentNll=trialNll;currentObjective=trialObjective;}
                if(accepted){bestObjectiveIteration=iteration;if(currentNll<bestNll){bestNll=currentNll;bestNllIteration=iteration;bestNllAmplitudes=amplitudes;}}
                iEvent=event;iIteration=iteration;iAccepted=accepted?1:0;iNll=currentNll;iNllBefore=previousNll;iNllAfterRl=rlNll;iChange=maximumChange;iPositionChange=maximumPositionChange;iCurvature=CurvatureRms(centroids);iPriorPenalty=currentPenalty.first;iCurvaturePenalty=currentPenalty.second;iObjective=currentObjective;iStep=accepted?acceptedStep:0;iterationOutput.Fill();outIterations=iteration;outMaxChange=maximumChange;
                if(!accepted){outStalled=1;break;}const double relativeObjectiveChange=std::abs(currentObjective-previousObjective)/std::max(std::abs(previousObjective),1e-12);
                if(maximumChange<convergence||(maximumPositionChange<positionConvergence&&relativeObjectiveChange<objectiveConvergence)){outConverged=1;break;}}
            for(size_t j=0;j<candidates.size();++j)candidates[j].amplitude=amplitudes[j];
            const auto bestNllCentroids=centroidValues(bestNllAmplitudes);finalNll=currentNll;outBestNll=bestNll;outBestNllIteration=bestNllIteration;outBestObjectiveIteration=bestObjectiveIteration;outCurvature=CurvatureRms(centroids);outObjective=currentObjective;outPriorPenalty=currentPenalty.first;outCurvaturePenalty=currentPenalty.second;outEvent=event;outObservations=observed.size();outColumns=data.columns.size();fitOutput.Fill();
            for(size_t k=0;k<data.columns.size();++k){const auto& column=data.columns[k];double totalScale=0,visiblePrediction=0;for(int j:byColumn[k]){totalScale+=candidates[j].amplitude;visiblePrediction+=candidates[j].amplitude*candidates[j].norm;}
                double bestScale=0,bestVisible=0;for(int j:byColumn[k]){bestScale+=bestNllAmplitudes[j];bestVisible+=bestNllAmplitudes[j]*candidates[j].norm;}
                const double coordinate=.5*(column.low+column.high);double t=(coordinate-linePoint[column.axis])/direction[column.axis];TVector3 centre=linePoint+t*direction+centroids[k][0]*basis.first+centroids[k][1]*basis.second;
                const TVector3 bestCentre=linePoint+t*direction+bestNllCentroids[k][0]*basis.first+bestNllCentroids[k][1]*basis.second;
                cEvent=event;cAxis=column.axis;cIndex=column.index;cLow=column.low;cHigh=column.high;cA=centroids[k][0];cB=centroids[k][1];cX=centre.X();cY=centre.Y();cZ=centre.Z();cInitial=column.scale;cScale=totalScale;cVisible=visiblePrediction;
                cBestA=bestNllCentroids[k][0];cBestB=bestNllCentroids[k][1];cBestX=bestCentre.X();cBestY=bestCentre.Y();cBestZ=bestCentre.Z();cBestScale=bestScale;cBestVisible=bestVisible;
                cCurvA=k>0&&k+1<centroids.size()?centroids[k+1][0]-2*centroids[k][0]+centroids[k-1][0]:0;cCurvB=k>0&&k+1<centroids.size()?centroids[k+1][1]-2*centroids[k][1]+centroids[k-1][1]:0;columnOutput.Fill();}
            if(saveCandidates)for(size_t j=0;j<candidates.size();++j){const auto& candidate=candidates[j];pEvent=event;pColumn=data.columns[candidate.column].index;pA=candidate.a;pB=candidate.b;pAmplitude=candidate.amplitude;pBestNllAmplitude=bestNllAmplitudes[j];pPrior=candidate.prior;pNorm=candidate.norm;candidateOutput.Fill();}
            std::cout<<"event "<<event<<" RL converged "<<outConverged<<" stalled "<<outStalled<<" iterations "<<outIterations<<" nll "<<initialNll<<" -> "<<finalNll
                     <<" objective "<<outObjective<<" prior penalty "<<outPriorPenalty<<" curvature penalty "<<outCurvaturePenalty<<" curvature RMS "<<outCurvature<<" mm\n";
        }
        output.cd();fitOutput.Write();columnOutput.Write();iterationOutput.Write();if(saveCandidates)candidateOutput.Write();
        TParameter<Long64_t>("rl_prior_residual_entries",residualSample.size()).Write();TParameter<double>("rl_grid_radius_mm",radius).Write();TParameter<double>("rl_grid_step_mm",gridStep).Write();
        TParameter<double>("rl_prior_bandwidth_mm",priorBandwidth).Write();TParameter<double>("rl_prior_strength",priorStrength).Write();TParameter<double>("rl_curvature_sigma_mm",curvatureSigma).Write();TParameter<double>("rl_curvature_strength",curvatureStrength).Write();
        TParameter<double>("rl_amplitude_convergence",convergence).Write();TParameter<double>("rl_position_convergence_mm",positionConvergence).Write();TParameter<double>("rl_objective_convergence",objectiveConvergence).Write();
        output.Close();std::cout<<"output "<<argv[2]<<"\n";
    } catch(const std::exception& error){std::cerr<<"ERROR: "<<error.what()<<"\n";return 1;}return 0;
}
