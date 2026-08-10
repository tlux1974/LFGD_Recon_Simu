#include "GlobalFitCommon.hxx"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
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
    auto* tracks=dynamic_cast<TTree*>(input.Get("mc_track_points"));
    if(!fits||!segments||!tracks){std::cerr<<"Missing global_fit, mc_virtual_segments, or mc_track_points\n";return 1;}auto truth=ReadTruthDirections(*segments);
    TFile output(outputName.c_str(),"RECREATE");const char* axis[3]={"x","y","z"};
    TH2D* seedComponent[3];TH2D* fitComponent[3];for(int k=0;k<3;++k){
        seedComponent[k]=new TH2D((std::string("mc_vs_seed_d")+axis[k]).c_str(),(std::string(";")+"MC d"+axis[k]+";seed d"+axis[k]).c_str(),100,-1,1,100,-1,1);
        fitComponent[k]=new TH2D((std::string("mc_vs_fit_d")+axis[k]).c_str(),(std::string(";")+"MC d"+axis[k]+";fit d"+axis[k]).c_str(),100,-1,1,100,-1,1);}
    TH1D seedAngle("mc_seed_angle",";angle(MC, seed) [deg];events",180,0,180),fitAngle("mc_fit_angle",";angle(MC, fit) [deg];events",180,0,180);
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
    int event,ndof;double mcDx,mcDy,mcDz,seedDx,seedDy,seedDz,fitDx,fitDy,fitDz,angleSeed,angleFit,chi2,chi2Ndof;TTree comparison("direction_comparison","MC, seed and fit direction comparison");
    comparison.Branch("event",&event);comparison.Branch("mc_dx",&mcDx);comparison.Branch("mc_dy",&mcDy);comparison.Branch("mc_dz",&mcDz);
    comparison.Branch("seed_dx",&seedDx);comparison.Branch("seed_dy",&seedDy);comparison.Branch("seed_dz",&seedDz);
    comparison.Branch("fit_dx",&fitDx);comparison.Branch("fit_dy",&fitDy);comparison.Branch("fit_dz",&fitDz);
    comparison.Branch("mc_seed_angle",&angleSeed);comparison.Branch("mc_fit_angle",&angleFit);
    comparison.Branch("chi2",&chi2);comparison.Branch("ndof",&ndof);comparison.Branch("chi2_ndof",&chi2Ndof);
    GlobalFitResult r;SetGlobalFitAddresses(*fits,r);for(Long64_t i=0;i<fits->GetEntries();++i){fits->GetEntry(i);auto it=truth.find(r.event);if(it==truth.end())continue;
        event=r.event;TVector3 mc=it->second,seed(r.seedDx,r.seedDy,r.seedDz),fit(r.fitDx,r.fitDy,r.fitDz);if(seed.Mag()==0||fit.Mag()==0)continue;seed=seed.Unit();fit=fit.Unit();
        mcDx=mc.X();mcDy=mc.Y();mcDz=mc.Z();seedDx=seed.X();seedDy=seed.Y();seedDz=seed.Z();fitDx=fit.X();fitDy=fit.Y();fitDz=fit.Z();
        angleSeed=Angle(mc,seed);angleFit=Angle(mc,fit);chi2=r.chi2;ndof=r.ndof;chi2Ndof=ndof>0?chi2/ndof:0;double m[3]={mcDx,mcDy,mcDz},s[3]={seedDx,seedDy,seedDz},f[3]={fitDx,fitDy,fitDz};
        for(int k=0;k<3;++k){seedComponent[k]->Fill(m[k],s[k]);fitComponent[k]->Fill(m[k],f[k]);}seedAngle.Fill(angleSeed);fitAngle.Fill(angleFit);
        chi2Plot.Fill(chi2);ndofPlot.Fill(ndof);if(ndof>0)reducedChi2Plot.Fill(chi2Ndof);comparison.Fill();}

    // Build fitted lines and associate segment primary IDs with muon PDGs via
    // the trajectory tree. Shared segment endpoints are quantized and filled once.
    std::map<int,std::pair<TVector3,TVector3>> fittedLines;SetGlobalFitAddresses(*fits,r);
    for(Long64_t i=0;i<fits->GetEntries();++i){fits->GetEntry(i);TVector3 d(r.fitDx,r.fitDy,r.fitDz);if(d.Mag()>0)fittedLines[r.event]={TVector3(r.fitX,r.fitY,r.fitZ),d.Unit()};}
    auto muonIds=ReadMuonTrackIds(*tracks);TH1D muonDistance("muon_fit_point_distance",";minimum distance from unique muon segment endpoint to fitted line [mm];endpoints",200,0,100);
    TH1D residualXY("muon_fit_residual_xy",";signed MC endpoint residual in XY [mm];endpoints",400,-100,100);
    TH1D residualXZ("muon_fit_residual_xz",";signed MC endpoint residual in XZ [mm];endpoints",400,-100,100);
    TH1D residualYZ("muon_fit_residual_yz",";signed MC endpoint residual in YZ [mm];endpoints",400,-100,100);
    TH1D distanceXY("muon_fit_distance_xy",";minimum MC endpoint distance to fit in XY [mm];endpoints",200,0,100);
    TH1D distanceXZ("muon_fit_distance_xz",";minimum MC endpoint distance to fit in XZ [mm];endpoints",200,0,100);
    TH1D distanceYZ("muon_fit_distance_yz",";minimum MC endpoint distance to fit in YZ [mm];endpoints",200,0,100);
    int distanceEvent=0,distanceTrack=0;double pointX=0,pointY=0,pointZ=0,distance=0,signedXY=0,signedXZ=0,signedYZ=0,viewDistanceXY=0,viewDistanceXZ=0,viewDistanceYZ=0;
    TTree distanceTree("muon_fit_distances","Unique muon segment endpoint distances to fitted trajectory");distanceTree.Branch("event",&distanceEvent);distanceTree.Branch("track_id",&distanceTrack);
    distanceTree.Branch("x",&pointX);distanceTree.Branch("y",&pointY);distanceTree.Branch("z",&pointZ);distanceTree.Branch("distance",&distance);
    distanceTree.Branch("residual_xy",&signedXY);distanceTree.Branch("residual_xz",&signedXZ);distanceTree.Branch("residual_yz",&signedYZ);
    distanceTree.Branch("distance_xy",&viewDistanceXY);distanceTree.Branch("distance_xz",&viewDistanceXZ);distanceTree.Branch("distance_yz",&viewDistanceYZ);
    int segmentEvent=0,detector=0,primary=0;double sx,sy,sz,ex,ey,ez;segments->SetBranchAddress("event",&segmentEvent);segments->SetBranchAddress("detector",&detector);
    segments->SetBranchAddress("primary_id",&primary);segments->SetBranchAddress("start_x",&sx);segments->SetBranchAddress("start_y",&sy);segments->SetBranchAddress("start_z",&sz);
    segments->SetBranchAddress("stop_x",&ex);segments->SetBranchAddress("stop_y",&ey);segments->SetBranchAddress("stop_z",&ez);
    std::set<std::tuple<int,int,long long,long long,long long>> used;
    for(Long64_t i=0;i<segments->GetEntries();++i){segments->GetEntry(i);if(detector!=0||!fittedLines.count(segmentEvent)||!muonIds[segmentEvent].count(primary))continue;
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
    output.cd();comparison.Write();distanceTree.Write();seedAngle.Write();fitAngle.Write();chi2Plot.Write();ndofPlot.Write();reducedChi2Plot.Write();for(int k=0;k<3;++k){seedComponent[k]->Write();fitComponent[k]->Write();}
    TCanvas canvas("direction_summary","Direction comparison",1200,500);canvas.Divide(2,1);canvas.cd(1);seedAngle.Draw();canvas.cd(2);fitAngle.Draw();canvas.Write();canvas.SaveAs("direction_comparison.png");
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
    std::cout<<"wrote "<<comparison.GetEntries()<<" comparisons to "<<outputName<<" and direction_comparison.png\n";return 0;
}
