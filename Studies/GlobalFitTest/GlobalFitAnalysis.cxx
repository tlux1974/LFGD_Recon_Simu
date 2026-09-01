#include "GlobalFitCommon.hxx"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TTree.h>
#include <TVector3.h>

namespace {
std::map<int,TVector3> ReadTruthDirections(TTree& tree) {
    int event=0,detector=0,primary=0;double x1,y1,z1,x2,y2,z2,length=0;
    tree.SetBranchAddress("event",&event);tree.SetBranchAddress("detector",&detector);
    tree.SetBranchAddress("primary_id",&primary);tree.SetBranchAddress("track_length",&length);
    tree.SetBranchAddress("start_x",&x1);tree.SetBranchAddress("start_y",&y1);tree.SetBranchAddress("start_z",&z1);
    tree.SetBranchAddress("stop_x",&x2);tree.SetBranchAddress("stop_y",&y2);tree.SetBranchAddress("stop_z",&z2);
    std::map<std::pair<int,int>,double> lengths;std::map<std::pair<int,int>,TVector3> sums;
    for(Long64_t i=0;i<tree.GetEntries();++i){tree.GetEntry(i);if(detector!=0)continue;auto key=std::make_pair(event,primary);
        lengths[key]+=length;sums[key]+=TVector3(x2-x1,y2-y1,z2-z1);}
    std::map<int,std::pair<double,int>> dominant;
    for(auto& item:lengths)if(!dominant.count(item.first.first)||item.second>dominant[item.first.first].first)
        dominant[item.first.first]={item.second,item.first.second};
    std::map<int,TVector3> result;for(auto& item:dominant){auto v=sums[{item.first,item.second.second}];if(v.Mag()>0)result[item.first]=v.Unit();}
    tree.ResetBranchAddresses();
    return result;
}
double Angle(const TVector3& a,const TVector3& b){return std::acos(std::clamp(a.Dot(b),-1.0,1.0))*180.0/M_PI;}

TVector3 ResidualVector(const TVector3& point,const TVector3& linePoint,const TVector3& direction) {
    const TVector3 unit=direction.Unit(),delta=point-linePoint;return delta-unit*delta.Dot(unit);
}

double SignedProjectionResidual(const TVector3& point,const TVector3& linePoint,
                                const TVector3& direction,int axisA,int axisB) {
    const double da=direction[axisA],db=direction[axisB],norm=std::hypot(da,db);
    if(norm<1e-12)return std::numeric_limits<double>::quiet_NaN();
    const double a=point[axisA]-linePoint[axisA],b=point[axisB]-linePoint[axisB];
    return (-db*a+da*b)/norm;
}

void FitPeak(TH1D& histogram,TF1& gaussian) {
    if(histogram.GetEntries()<10)return;
    int peak=histogram.GetMaximumBin(),left=peak,right=peak;
    const double half=.5*histogram.GetBinContent(peak);
    while(left>1&&histogram.GetBinContent(left)>half)--left;
    while(right<histogram.GetNbinsX()&&histogram.GetBinContent(right)>half)++right;
    const double centre=histogram.GetBinCenter(peak);
    const double sigma=std::max(histogram.GetBinWidth(peak),
        (histogram.GetBinCenter(right)-histogram.GetBinCenter(left))/2.355);
    gaussian.SetRange(centre-2*sigma,centre+2*sigma);histogram.Fit(&gaussian,"RQ");
}

std::string GaussianLabel(const char* prefix,const TF1& function) {
    std::ostringstream text;text<<std::fixed<<std::setprecision(3)<<prefix
        <<": #mu="<<function.GetParameter(1)<<" mm, #sigma="<<std::abs(function.GetParameter(2))<<" mm";
    return text.str();
}

void DrawResidualOverlay(TH1D& fibreInput,TH1D& segmentInput,const char* canvasName,
                         const char* title) {
    TH1D fibre(fibreInput),segment(segmentInput);fibre.SetDirectory(nullptr);segment.SetDirectory(nullptr);
    const bool fitFibre=fibre.Integral()>=10,fitSegment=segment.Integral()>=10;
    if(fibre.Integral()>0)fibre.Scale(1.0/fibre.Integral());
    if(segment.Integral()>0)segment.Scale(1.0/segment.Integral());
    fibre.SetLineColor(kRed+1);fibre.SetLineWidth(2);segment.SetLineColor(kBlue+1);segment.SetLineWidth(2);
    fibre.SetTitle(title);fibre.GetYaxis()->SetTitle("normalized entries");
    TF1 fibreGaussian((std::string(canvasName)+"_fibre_gaussian").c_str(),"gaus",-10,10);
    TF1 segmentGaussian((std::string(canvasName)+"_segment_gaussian").c_str(),"gaus",-10,10);
    if(fitFibre)FitPeak(fibre,fibreGaussian);
    if(fitSegment)FitPeak(segment,segmentGaussian);
    fibreGaussian.SetLineColor(kRed+1);fibreGaussian.SetLineStyle(2);segmentGaussian.SetLineColor(kBlue+1);segmentGaussian.SetLineStyle(2);
    fibre.SetMaximum(1.2*std::max(fibre.GetMaximum(),segment.GetMaximum()));fibre.DrawCopy("HIST");segment.DrawCopy("HIST SAME");
    if(fitFibre)fibreGaussian.DrawClone("SAME");
    if(fitSegment)segmentGaussian.DrawClone("SAME");
    TLegend legend(.12,.70,.88,.89);legend.SetBorderSize(0);legend.SetFillStyle(0);
    const std::string fibreLabel=fitFibre?GaussianLabel("Fibre/light fit",fibreGaussian):"Fibre/light fit: insufficient entries in #pm10 mm";
    const std::string segmentLabel=fitSegment?GaussianLabel("MC-start line",segmentGaussian):"MC-start line: insufficient entries in #pm10 mm";
    legend.AddEntry(&fibre,fibreLabel.c_str(),"l");legend.AddEntry(&segment,segmentLabel.c_str(),"l");legend.DrawClone();
    fibreInput.Write();segmentInput.Write();fibreGaussian.Write();segmentGaussian.Write();
}

std::map<int,std::set<int>> ReadMuonTrackIds(TTree& tree) {
    int event=0,trackId=0,pdg=0;tree.SetBranchAddress("event",&event);
    tree.SetBranchAddress("track_id",&trackId);tree.SetBranchAddress("pdg",&pdg);
    std::map<int,std::set<int>> result;for(Long64_t i=0;i<tree.GetEntries();++i){tree.GetEntry(i);if(std::abs(pdg)==13)result[event].insert(trackId);}
    tree.ResetBranchAddresses();
    return result;
}
}

int main(int argc,char** argv){
    auto usage=[](){std::cout<<R"(Usage:
  global_fit_analysis FITTED.root [OUTPUT=global_fit_analysis.root]

Creates direction, chi2, signed/unsigned per-view muon residual, and Gaussian
summary plots. OUTPUT must use the shown KEY=value form.
)";};
    if(argc==1||(argc==2&&std::string(argv[1])=="--help")){usage();return 0;}
    if(argc<2||argc>3){usage();return 2;}
    std::string outputName="global_fit_analysis.root";
    if(argc==3){std::string argument=argv[2];const std::string prefix="OUTPUT=";if(argument.rfind(prefix,0)!=0){std::cerr<<"Expected OUTPUT=filename.root\n";return 2;}outputName=argument.substr(prefix.size());}
    TFile input(argv[1],"READ");
    auto* fits=dynamic_cast<TTree*>(input.Get("global_fit"));auto* segments=dynamic_cast<TTree*>(input.Get("mc_virtual_segments"));
    auto* tracks=dynamic_cast<TTree*>(input.Get("mc_track_points"));auto* mcLineFits=dynamic_cast<TTree*>(input.Get("mc_segment_line_fit"));
    auto* mcLineResiduals=dynamic_cast<TTree*>(input.Get("mc_segment_line_residuals"));
    if(!fits||!segments||!tracks||!mcLineFits||!mcLineResiduals){std::cerr<<"Missing global_fit, MC input trees, or MC segment benchmark trees; rerun global_light_fit with the current version\n";return 1;}auto truth=ReadTruthDirections(*segments);
    TFile output(outputName.c_str(),"RECREATE");const char* axis[3]={"x","y","z"};
    TH2D* seedComponent[3];TH2D* fitComponent[3];for(int k=0;k<3;++k){
        seedComponent[k]=new TH2D((std::string("mc_vs_seed_d")+axis[k]).c_str(),(std::string(";")+"MC d"+axis[k]+";seed d"+axis[k]).c_str(),100,-1,1,100,-1,1);
        fitComponent[k]=new TH2D((std::string("mc_vs_fit_d")+axis[k]).c_str(),(std::string(";")+"MC d"+axis[k]+";fit d"+axis[k]).c_str(),100,-1,1,100,-1,1);}
    TH1D seedAngle("mc_seed_angle",";angle(MC, seed) [deg];events",180,0,180),fitAngle("mc_fit_angle",";angle(MC, fit) [deg];events",180,0,180),seedFitAngle("seed_fit_angle",";angle(seed, fit) [deg];events",180,0,180);
    // ROOT 6.32 can crash while extending a file-owned histogram over many
    // orders of magnitude. Determine finite ranges before creating them.
    const double maximumChi2=std::max(1.0,1.05*fits->GetMaximum("chi2"));
    const double maximumNdof=std::max(1.0,1.05*fits->GetMaximum("ndof"));
    double maximumReducedChi2=1.0;GlobalFitResult rangeResult;SetGlobalFitAddresses(*fits,rangeResult);
    for(Long64_t i=0;i<fits->GetEntries();++i){fits->GetEntry(i);if(rangeResult.ndof>0&&std::isfinite(rangeResult.chi2))
        maximumReducedChi2=std::max(maximumReducedChi2,1.05*rangeResult.chi2/rangeResult.ndof);}
    TH1D chi2Plot("fit_chi2",";#chi^{2};events",200,0,maximumChi2);
    TH1D ndofPlot("fit_ndof",";degrees of freedom;events",250,0,maximumNdof);
    TH1D reducedChi2Plot("fit_chi2_ndof",";#chi^{2}/ndof;events",200,0,maximumReducedChi2);
    int event,status,ndof;double mcDx,mcDy,mcDz,seedDx,seedDy,seedDz,fitDx,fitDy,fitDz,seedResidualDx,seedResidualDy,seedResidualDz,fitResidualDx,fitResidualDy,fitResidualDz,angleSeed,angleFit,angleSeedFit,chi2,chi2Ndof,seedNll,fitNll,nllImprovement,seedChi2,chi2Improvement;TTree comparison("direction_comparison","MC, seed and fit direction comparison");
    comparison.Branch("event",&event);comparison.Branch("mc_dx",&mcDx);comparison.Branch("mc_dy",&mcDy);comparison.Branch("mc_dz",&mcDz);
    comparison.Branch("seed_dx",&seedDx);comparison.Branch("seed_dy",&seedDy);comparison.Branch("seed_dz",&seedDz);
    comparison.Branch("fit_dx",&fitDx);comparison.Branch("fit_dy",&fitDy);comparison.Branch("fit_dz",&fitDz);
    comparison.Branch("seed_residual_dx",&seedResidualDx);comparison.Branch("seed_residual_dy",&seedResidualDy);comparison.Branch("seed_residual_dz",&seedResidualDz);
    comparison.Branch("fit_residual_dx",&fitResidualDx);comparison.Branch("fit_residual_dy",&fitResidualDy);comparison.Branch("fit_residual_dz",&fitResidualDz);
    comparison.Branch("mc_seed_angle",&angleSeed);comparison.Branch("mc_fit_angle",&angleFit);comparison.Branch("seed_fit_angle",&angleSeedFit);comparison.Branch("status",&status);
    comparison.Branch("seed_nll",&seedNll);comparison.Branch("fit_nll",&fitNll);comparison.Branch("nll_improvement",&nllImprovement);
    comparison.Branch("seed_chi2",&seedChi2);comparison.Branch("chi2_improvement",&chi2Improvement);
    comparison.Branch("chi2",&chi2);comparison.Branch("ndof",&ndof);comparison.Branch("chi2_ndof",&chi2Ndof);
    GlobalFitResult r;SetGlobalFitAddresses(*fits,r);for(Long64_t i=0;i<fits->GetEntries();++i){fits->GetEntry(i);auto it=truth.find(r.event);if(it==truth.end())continue;
        event=r.event;TVector3 mc=it->second,seed(r.seedDx,r.seedDy,r.seedDz),fit(r.fitDx,r.fitDy,r.fitDz);if(seed.Mag()==0||fit.Mag()==0)continue;seed=seed.Unit();fit=fit.Unit();
        mcDx=mc.X();mcDy=mc.Y();mcDz=mc.Z();seedDx=seed.X();seedDy=seed.Y();seedDz=seed.Z();fitDx=fit.X();fitDy=fit.Y();fitDz=fit.Z();
        if(mc.Dot(seed)<0)seed*=-1;
        if(mc.Dot(fit)<0)fit*=-1;
        seedDx=seed.X();seedDy=seed.Y();seedDz=seed.Z();fitDx=fit.X();fitDy=fit.Y();fitDz=fit.Z();
        seedResidualDx=seedDx-mcDx;seedResidualDy=seedDy-mcDy;seedResidualDz=seedDz-mcDz;fitResidualDx=fitDx-mcDx;fitResidualDy=fitDy-mcDy;fitResidualDz=fitDz-mcDz;
        angleSeed=Angle(mc,seed);angleFit=Angle(mc,fit);angleSeedFit=Angle(seed,fit);status=r.status;seedNll=r.seedNll;fitNll=r.nll;nllImprovement=seedNll-fitNll;seedChi2=r.seedChi2;chi2Improvement=seedChi2-r.chi2;chi2=r.chi2;ndof=r.ndof;chi2Ndof=ndof>0?chi2/ndof:0;double m[3]={mcDx,mcDy,mcDz},s[3]={seedDx,seedDy,seedDz},f[3]={fitDx,fitDy,fitDz};
        for(int k=0;k<3;++k){seedComponent[k]->Fill(m[k],s[k]);fitComponent[k]->Fill(m[k],f[k]);}seedAngle.Fill(angleSeed);fitAngle.Fill(angleFit);
        seedFitAngle.Fill(angleSeedFit);chi2Plot.Fill(chi2);ndofPlot.Fill(ndof);if(ndof>0)reducedChi2Plot.Fill(chi2Ndof);comparison.Fill();}

    // Build fitted lines and associate segment primary IDs with muon PDGs via
    // the trajectory tree. Shared segment endpoints are quantized and filled once.
    std::map<int,std::pair<TVector3,TVector3>> fittedLines,seedLines;SetGlobalFitAddresses(*fits,r);
    for(Long64_t i=0;i<fits->GetEntries();++i){fits->GetEntry(i);TVector3 d(r.fitDx,r.fitDy,r.fitDz),s(r.seedDx,r.seedDy,r.seedDz);if(d.Mag()>0)fittedLines[r.event]={TVector3(r.fitX,r.fitY,r.fitZ),d.Unit()};if(s.Mag()>0)seedLines[r.event]={TVector3(r.seedX,r.seedY,r.seedZ),s.Unit()};}
    auto muonIds=ReadMuonTrackIds(*tracks);TH1D muonDistance("muon_fit_point_distance",";minimum distance from unique muon segment endpoint to fitted line [mm];endpoints",200,0,100);
    TH1D residualXY("muon_fit_residual_xy",";signed MC endpoint residual in XY [mm];endpoints",400,-100,100);
    TH1D residualXZ("muon_fit_residual_xz",";signed MC endpoint residual in XZ [mm];endpoints",400,-100,100);
    TH1D residualYZ("muon_fit_residual_yz",";signed MC endpoint residual in YZ [mm];endpoints",400,-100,100);
    TH1D distanceXY("muon_fit_distance_xy",";minimum MC endpoint distance to fit in XY [mm];endpoints",200,0,100);
    TH1D distanceXZ("muon_fit_distance_xz",";minimum MC endpoint distance to fit in XZ [mm];endpoints",200,0,100);
    TH1D distanceYZ("muon_fit_distance_yz",";minimum MC endpoint distance to fit in YZ [mm];endpoints",200,0,100);
    const char* viewName[3]={"xy","xz","yz"};
    std::array<std::unique_ptr<TH1D>,3> fibreStartResidual,fibreEndResidual,segmentStartResidual,segmentEndResidual;
    for(int view=0;view<3;++view){
        fibreStartResidual[view]=std::make_unique<TH1D>((std::string("fibre_fit_start_residual_")+viewName[view]).c_str(),(";signed start residual in "+std::string(viewName[view])+" [mm];points").c_str(),200,-10,10);
        fibreEndResidual[view]=std::make_unique<TH1D>((std::string("fibre_fit_end_residual_")+viewName[view]).c_str(),(";signed end residual in "+std::string(viewName[view])+" [mm];points").c_str(),200,-10,10);
        segmentStartResidual[view]=std::make_unique<TH1D>((std::string("mc_line_start_residual_")+viewName[view]).c_str(),(";signed start residual in "+std::string(viewName[view])+" [mm];points").c_str(),200,-10,10);
        segmentEndResidual[view]=std::make_unique<TH1D>((std::string("mc_line_end_residual_")+viewName[view]).c_str(),(";signed end residual in "+std::string(viewName[view])+" [mm];points").c_str(),200,-10,10);
    }
    int benchmarkEvent=0,benchmarkPrimary=0;std::map<int,int> benchmarkPrimaryByEvent;
    mcLineFits->SetBranchAddress("event",&benchmarkEvent);mcLineFits->SetBranchAddress("primary_id",&benchmarkPrimary);
    for(Long64_t i=0;i<mcLineFits->GetEntries();++i){mcLineFits->GetEntry(i);benchmarkPrimaryByEvent[benchmarkEvent]=benchmarkPrimary;}mcLineFits->ResetBranchAddresses();
    int distanceEvent=0,distanceTrack=0;double pointX=0,pointY=0,pointZ=0,distance=0,signedXY=0,signedXZ=0,signedYZ=0,viewDistanceXY=0,viewDistanceXZ=0,viewDistanceYZ=0;
    TTree distanceTree("muon_fit_distances","Unique muon segment endpoint distances to fitted trajectory");distanceTree.Branch("event",&distanceEvent);distanceTree.Branch("track_id",&distanceTrack);
    distanceTree.Branch("x",&pointX);distanceTree.Branch("y",&pointY);distanceTree.Branch("z",&pointZ);distanceTree.Branch("distance",&distance);
    distanceTree.Branch("residual_xy",&signedXY);distanceTree.Branch("residual_xz",&signedXZ);distanceTree.Branch("residual_yz",&signedYZ);
    distanceTree.Branch("distance_xy",&viewDistanceXY);distanceTree.Branch("distance_xz",&viewDistanceXZ);distanceTree.Branch("distance_yz",&viewDistanceYZ);
    int residualEvent=0,residualPrimary=0,residualSegment=0,residualEndpoint=0;double residualPointX=0,residualPointY=0,residualPointZ=0;
    double seedRx=0,seedRy=0,seedRz=0,seedDistance=0,fitRx=0,fitRy=0,fitRz=0,fitDistance=0;
    TTree xyzResidualTree("mc_seed_fit_residuals","MC segment endpoint residual vectors to seed and fitted lines");
    xyzResidualTree.Branch("event",&residualEvent);xyzResidualTree.Branch("primary_id",&residualPrimary);xyzResidualTree.Branch("segment",&residualSegment);xyzResidualTree.Branch("endpoint",&residualEndpoint);
    xyzResidualTree.Branch("x",&residualPointX);xyzResidualTree.Branch("y",&residualPointY);xyzResidualTree.Branch("z",&residualPointZ);
    xyzResidualTree.Branch("seed_residual_x",&seedRx);xyzResidualTree.Branch("seed_residual_y",&seedRy);xyzResidualTree.Branch("seed_residual_z",&seedRz);xyzResidualTree.Branch("seed_distance",&seedDistance);
    xyzResidualTree.Branch("fit_residual_x",&fitRx);xyzResidualTree.Branch("fit_residual_y",&fitRy);xyzResidualTree.Branch("fit_residual_z",&fitRz);xyzResidualTree.Branch("fit_distance",&fitDistance);
    std::array<std::unique_ptr<TH1D>,3> seedXyzResidual,fitXyzResidual;for(int k=0;k<3;++k){seedXyzResidual[k]=std::make_unique<TH1D>((std::string("mc_seed_residual_")+axis[k]).c_str(),(std::string(";MC residual ")+axis[k]+" to seed line [mm];segment endpoints").c_str(),400,-100,100);fitXyzResidual[k]=std::make_unique<TH1D>((std::string("mc_fit_residual_")+axis[k]).c_str(),(std::string(";MC residual ")+axis[k]+" to fitted line [mm];segment endpoints").c_str(),400,-100,100);}
    int segmentEvent=0,detector=0,primary=0,segmentNumber=0;double sx,sy,sz,ex,ey,ez;segments->SetBranchAddress("event",&segmentEvent);segments->SetBranchAddress("detector",&detector);segments->SetBranchAddress("segment",&segmentNumber);
    segments->SetBranchAddress("primary_id",&primary);segments->SetBranchAddress("start_x",&sx);segments->SetBranchAddress("start_y",&sy);segments->SetBranchAddress("start_z",&sz);
    segments->SetBranchAddress("stop_x",&ex);segments->SetBranchAddress("stop_y",&ey);segments->SetBranchAddress("stop_z",&ez);
    std::set<std::tuple<int,int,long long,long long,long long>> used;
    for(Long64_t i=0;i<segments->GetEntries();++i){segments->GetEntry(i);if(detector!=0||!fittedLines.count(segmentEvent)||!muonIds[segmentEvent].count(primary)||!benchmarkPrimaryByEvent.count(segmentEvent)||benchmarkPrimaryByEvent[segmentEvent]!=primary)continue;
        const auto& comparisonLine=fittedLines[segmentEvent];const TVector3 startPoint(sx,sy,sz),endPoint(ex,ey,ez);
        if(seedLines.count(segmentEvent))for(int endpoint=0;endpoint<2;++endpoint){const TVector3 point=endpoint?endPoint:startPoint;const TVector3 sr=ResidualVector(point,seedLines[segmentEvent].first,seedLines[segmentEvent].second),fr=ResidualVector(point,comparisonLine.first,comparisonLine.second);
            residualEvent=segmentEvent;residualPrimary=primary;residualSegment=segmentNumber;residualEndpoint=endpoint;residualPointX=point.X();residualPointY=point.Y();residualPointZ=point.Z();seedRx=sr.X();seedRy=sr.Y();seedRz=sr.Z();seedDistance=sr.Mag();fitRx=fr.X();fitRy=fr.Y();fitRz=fr.Z();fitDistance=fr.Mag();
            const double sv[3]={seedRx,seedRy,seedRz},fv[3]={fitRx,fitRy,fitRz};for(int k=0;k<3;++k){seedXyzResidual[k]->Fill(sv[k]);fitXyzResidual[k]->Fill(fv[k]);}xyzResidualTree.Fill();}
        const double fibreStart[3]={SignedProjectionResidual(startPoint,comparisonLine.first,comparisonLine.second,0,1),SignedProjectionResidual(startPoint,comparisonLine.first,comparisonLine.second,0,2),SignedProjectionResidual(startPoint,comparisonLine.first,comparisonLine.second,1,2)};
        const double fibreEnd[3]={SignedProjectionResidual(endPoint,comparisonLine.first,comparisonLine.second,0,1),SignedProjectionResidual(endPoint,comparisonLine.first,comparisonLine.second,0,2),SignedProjectionResidual(endPoint,comparisonLine.first,comparisonLine.second,1,2)};
        for(int view=0;view<3;++view){if(std::isfinite(fibreStart[view]))fibreStartResidual[view]->Fill(fibreStart[view]);if(std::isfinite(fibreEnd[view]))fibreEndResidual[view]->Fill(fibreEnd[view]);}
        for(const TVector3& point:{TVector3(sx,sy,sz),TVector3(ex,ey,ez)}){auto key=std::make_tuple(segmentEvent,primary,std::llround(point.X()*1000),std::llround(point.Y()*1000),std::llround(point.Z()*1000));
            if(!used.insert(key).second)continue;
            const auto& line=fittedLines[segmentEvent];distance=((point-line.first).Cross(line.second)).Mag();
            signedXY=SignedProjectionResidual(point,line.first,line.second,0,1);signedXZ=SignedProjectionResidual(point,line.first,line.second,0,2);
            signedYZ=SignedProjectionResidual(point,line.first,line.second,1,2);viewDistanceXY=std::abs(signedXY);viewDistanceXZ=std::abs(signedXZ);viewDistanceYZ=std::abs(signedYZ);
            distanceEvent=segmentEvent;distanceTrack=primary;pointX=point.X();pointY=point.Y();pointZ=point.Z();
            muonDistance.Fill(distance);if(std::isfinite(signedXY))residualXY.Fill(signedXY);if(std::isfinite(signedXZ))residualXZ.Fill(signedXZ);
            if(std::isfinite(signedYZ))residualYZ.Fill(signedYZ);
            if(std::isfinite(viewDistanceXY))distanceXY.Fill(viewDistanceXY);
            if(std::isfinite(viewDistanceXZ))distanceXZ.Fill(viewDistanceXZ);
            if(std::isfinite(viewDistanceYZ))distanceYZ.Fill(viewDistanceYZ);
            distanceTree.Fill();}}
    double benchmarkStartXY=0,benchmarkStartXZ=0,benchmarkStartYZ=0,benchmarkEndXY=0,benchmarkEndXZ=0,benchmarkEndYZ=0;
    mcLineResiduals->SetBranchAddress("start_residual_xy",&benchmarkStartXY);mcLineResiduals->SetBranchAddress("start_residual_xz",&benchmarkStartXZ);mcLineResiduals->SetBranchAddress("start_residual_yz",&benchmarkStartYZ);
    mcLineResiduals->SetBranchAddress("end_residual_xy",&benchmarkEndXY);mcLineResiduals->SetBranchAddress("end_residual_xz",&benchmarkEndXZ);mcLineResiduals->SetBranchAddress("end_residual_yz",&benchmarkEndYZ);
    for(Long64_t i=0;i<mcLineResiduals->GetEntries();++i){mcLineResiduals->GetEntry(i);const double starts[3]={benchmarkStartXY,benchmarkStartXZ,benchmarkStartYZ},ends[3]={benchmarkEndXY,benchmarkEndXZ,benchmarkEndYZ};
        for(int view=0;view<3;++view){if(std::isfinite(starts[view]))segmentStartResidual[view]->Fill(starts[view]);if(std::isfinite(ends[view]))segmentEndResidual[view]->Fill(ends[view]);}}
    mcLineResiduals->ResetBranchAddresses();
    output.cd();comparison.Write();distanceTree.Write();xyzResidualTree.Write();seedAngle.Write();fitAngle.Write();seedFitAngle.Write();chi2Plot.Write();ndofPlot.Write();reducedChi2Plot.Write();for(int k=0;k<3;++k){seedComponent[k]->Write();fitComponent[k]->Write();seedXyzResidual[k]->Write();fitXyzResidual[k]->Write();}
    TCanvas canvas("direction_summary","Direction comparison",1500,500);canvas.Divide(3,1);canvas.cd(1);seedAngle.Draw();canvas.cd(2);fitAngle.Draw();canvas.cd(3);seedFitAngle.Draw();canvas.Write();canvas.SaveAs("direction_comparison.png");
    TCanvas xyzCanvas("mc_seed_fit_xyz_residual_summary","MC segment endpoint residual components",1500,900);xyzCanvas.Divide(3,2);
    for(int k=0;k<3;++k){xyzCanvas.cd(k+1);seedXyzResidual[k]->SetLineColor(kBlue+1);seedXyzResidual[k]->Draw("HIST");xyzCanvas.cd(k+4);fitXyzResidual[k]->SetLineColor(kRed+1);fitXyzResidual[k]->Draw("HIST");}xyzCanvas.Write();xyzCanvas.SaveAs("mc_seed_fit_xyz_residuals.png");
    TCanvas quality("fit_quality_summary","Fit quality",1500,500);quality.Divide(3,1);quality.cd(1);chi2Plot.Draw();quality.cd(2);ndofPlot.Draw();quality.cd(3);reducedChi2Plot.Draw();quality.Write();quality.SaveAs("fit_quality.png");
    TCanvas distanceCanvas("muon_distance_summary","Muon endpoint distance to fit",900,700);TF1 gaussian("muon_distance_gaussian","gaus",0,100);FitPeak(muonDistance,gaussian);
    muonDistance.Draw();gaussian.SetLineColor(kRed);if(muonDistance.GetEntries()>=10)gaussian.Draw("SAME");muonDistance.Write();gaussian.Write();distanceCanvas.Write();distanceCanvas.SaveAs("muon_fit_distance.png");
    TF1 gaussianXY("muon_residual_xy_gaussian","gaus",-100,100),gaussianXZ("muon_residual_xz_gaussian","gaus",-100,100),gaussianYZ("muon_residual_yz_gaussian","gaus",-100,100);
    FitPeak(residualXY,gaussianXY);FitPeak(residualXZ,gaussianXZ);FitPeak(residualYZ,gaussianYZ);
    TCanvas residualCanvas("muon_signed_residual_summary","Signed muon endpoint residuals",1500,500);residualCanvas.Divide(3,1);
    residualCanvas.cd(1);residualXY.Draw();gaussianXY.SetLineColor(kRed);gaussianXY.Draw("SAME");residualCanvas.cd(2);residualXZ.Draw();gaussianXZ.SetLineColor(kRed);gaussianXZ.Draw("SAME");
    residualCanvas.cd(3);residualYZ.Draw();gaussianYZ.SetLineColor(kRed);gaussianYZ.Draw("SAME");residualXY.Write();residualXZ.Write();residualYZ.Write();gaussianXY.Write();gaussianXZ.Write();gaussianYZ.Write();
    residualCanvas.Write();residualCanvas.SaveAs("muon_fit_signed_residuals.png");
    TF1 distanceXYGaussian("muon_distance_xy_gaussian","gaus",0,100),distanceXZGaussian("muon_distance_xz_gaussian","gaus",0,100),distanceYZGaussian("muon_distance_yz_gaussian","gaus",0,100);
    FitPeak(distanceXY,distanceXYGaussian);FitPeak(distanceXZ,distanceXZGaussian);FitPeak(distanceYZ,distanceYZGaussian);
    TCanvas viewDistanceCanvas("muon_view_distance_summary","Muon endpoint distances by view",1500,500);viewDistanceCanvas.Divide(3,1);
    viewDistanceCanvas.cd(1);distanceXY.Draw();distanceXYGaussian.SetLineColor(kRed);distanceXYGaussian.Draw("SAME");
    viewDistanceCanvas.cd(2);distanceXZ.Draw();distanceXZGaussian.SetLineColor(kRed);distanceXZGaussian.Draw("SAME");
    viewDistanceCanvas.cd(3);distanceYZ.Draw();distanceYZGaussian.SetLineColor(kRed);distanceYZGaussian.Draw("SAME");
    distanceXY.Write();distanceXZ.Write();distanceYZ.Write();distanceXYGaussian.Write();distanceXZGaussian.Write();distanceYZGaussian.Write();
    viewDistanceCanvas.Write();viewDistanceCanvas.SaveAs("muon_fit_view_distances.png");
    TCanvas startOverlay("start_residual_overlay","Segment-start residual comparison",1800,600);startOverlay.Divide(3,1);
    for(int view=0;view<3;++view){startOverlay.cd(view+1);const std::string title=std::string("Segment starts, ")+viewName[view]+";signed residual [mm];normalized entries";
        DrawResidualOverlay(*fibreStartResidual[view],*segmentStartResidual[view],(std::string("start_")+viewName[view]).c_str(),title.c_str());}
    startOverlay.Write();startOverlay.SaveAs("fibre_vs_mc_line_start_residuals.png");
    TCanvas endOverlay("end_residual_overlay","Segment-end residual comparison",1800,600);endOverlay.Divide(3,1);
    for(int view=0;view<3;++view){endOverlay.cd(view+1);const std::string title=std::string("Segment ends, ")+viewName[view]+";signed residual [mm];normalized entries";
        DrawResidualOverlay(*fibreEndResidual[view],*segmentEndResidual[view],(std::string("end_")+viewName[view]).c_str(),title.c_str());}
    endOverlay.Write();endOverlay.SaveAs("fibre_vs_mc_line_end_residuals.png");
    std::cout<<"wrote "<<comparison.GetEntries()<<" comparisons to "<<outputName<<" and direction_comparison.png\n";return 0;
}
