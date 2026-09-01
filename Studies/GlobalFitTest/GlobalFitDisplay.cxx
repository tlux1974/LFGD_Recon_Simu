#include "GlobalFitCommon.hxx"

#include <algorithm>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH2D.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

namespace {
struct Hit { double x,y,z,q; int projection; };
struct Segment { double x1,y1,z1,x2,y2,z2; };
struct RlPoint { int column;double x,y,z,bestX,bestY,bestZ; };

std::vector<RlPoint> ReadRlPoints(TTree& tree,int wanted) {
    int event=0,column=0;double x=0,y=0,z=0,bx=0,by=0,bz=0;
    tree.SetBranchAddress("event",&event);tree.SetBranchAddress("column_index",&column);
    tree.SetBranchAddress("centre_x",&x);tree.SetBranchAddress("centre_y",&y);tree.SetBranchAddress("centre_z",&z);
    tree.SetBranchAddress("best_nll_centre_x",&bx);tree.SetBranchAddress("best_nll_centre_y",&by);tree.SetBranchAddress("best_nll_centre_z",&bz);
    std::vector<RlPoint> out;
    for(Long64_t i=0;i<tree.GetEntries();++i){tree.GetEntry(i);if(event==wanted)out.push_back({column,x,y,z,bx,by,bz});}
    tree.ResetBranchAddresses();std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){return a.column<b.column;});return out;
}

std::vector<Segment> ReadSegments(TTree& tree,int wanted) {
    int event=0,detector=0,primary=0;double x1,y1,z1,x2,y2,z2;
    tree.SetBranchAddress("event",&event);tree.SetBranchAddress("detector",&detector);
    tree.SetBranchAddress("primary_id",&primary);
    tree.SetBranchAddress("start_x",&x1);tree.SetBranchAddress("start_y",&y1);tree.SetBranchAddress("start_z",&z1);
    tree.SetBranchAddress("stop_x",&x2);tree.SetBranchAddress("stop_y",&y2);tree.SetBranchAddress("stop_z",&z2);
    std::vector<Segment> out;
    for(Long64_t i=0;i<tree.GetEntries();++i){tree.GetEntry(i);if(event<wanted)continue;if(event>wanted)break;
        if(detector==0&&primary==1)out.push_back({x1,y1,z1,x2,y2,z2});}
    return out;
}

std::vector<Hit> ReadHits(TTree& tree,int wanted,double cut) {
    int event=0,projection=0;double x=0,y=0,z=0,q=0;
    tree.SetBranchAddress("event",&event);tree.SetBranchAddress("x",&x);
    tree.SetBranchAddress("y",&y);tree.SetBranchAddress("z",&z);
    tree.SetBranchAddress("charge",&q);tree.SetBranchAddress("projection",&projection);
    std::vector<Hit> out;
    for(Long64_t i=0;i<tree.GetEntries();++i) {
        tree.GetEntry(i);
        if(event<wanted) continue;
        if(event>wanted) break;
        if(q>=cut) out.push_back({x,y,z,q,projection});
    }
    return out;
}

bool ReadFit(TTree& tree,int wanted,GlobalFitResult& fit) {
    SetGlobalFitAddresses(tree,fit);std::string* inputTree=nullptr;tree.SetBranchAddress("input_tree",&inputTree);
    for(Long64_t i=0;i<tree.GetEntries();++i){tree.GetEntry(i);
        if(fit.event==wanted){if(inputTree)fit.inputTree=*inputTree;return true;}}return false;
}

void DrawProjection(int pad,TCanvas& canvas,const std::vector<Hit>& hits,
                    const std::vector<Hit>& selectedHits,bool showSelected,
                    const std::vector<Hit>& dbscanHits,bool showDbscan,
                    const std::vector<Segment>& segments,
                    const std::vector<RlPoint>& rlPoints,const std::string& rlSolution,
                    const GlobalFitResult& fit,int a,int b,int requiredProjection,
                    const char* title,const char* an,const char* bn,
                    const std::string& colourMode) {
    canvas.cd(pad);std::vector<double> av,bv,qv;
    for(auto& h:hits)if(h.projection==requiredProjection){double c[3]={h.x,h.y,h.z};av.push_back(c[a]);bv.push_back(c[b]);qv.push_back(h.q);}
    if(colourMode=="off" || colourMode=="all" || av.empty()) {
        auto* graph=new TGraph(av.size(),av.data(),bv.data());graph->SetTitle((std::string(title)+";"+an+" [mm];"+bn+" [mm]").c_str());
        graph->SetMarkerStyle(20);graph->SetMarkerSize(.55);graph->SetMarkerColor(kBlue+1);graph->Draw("AP");
    }
    else {
        const auto aminmax=std::minmax_element(av.begin(),av.end()),bminmax=std::minmax_element(bv.begin(),bv.end());
        const double alo=*aminmax.first-5,ahi=*aminmax.second+5,blo=*bminmax.first-5,bhi=*bminmax.second+5;
        const int na=std::max(1,static_cast<int>(std::lround((ahi-alo)/10.0)));
        const int nb=std::max(1,static_cast<int>(std::lround((bhi-blo)/10.0)));
        auto* charge=new TH2D(("charge_view_"+std::to_string(pad)).c_str(),
            (std::string(title)+";"+an+" [mm];"+bn+" [mm];fibre charge").c_str(),na,alo,ahi,nb,blo,bhi);
        charge->SetDirectory(nullptr);for(size_t i=0;i<av.size();++i)charge->Fill(av[i],bv[i],qv[i]);
        gPad->SetLogz(colourMode=="log");charge->SetMinimum(colourMode=="log"?std::max(1.0,fit.minimumCharge):0.0);charge->Draw("COLZ");
    }
    std::vector<double> selectedA,selectedB;if(showSelected)for(const auto& h:selectedHits)if(h.projection==requiredProjection){double c[3]={h.x,h.y,h.z};selectedA.push_back(c[a]);selectedB.push_back(c[b]);}
    if(showSelected&&!selectedA.empty()){auto* accepted=new TGraph(selectedA.size(),selectedA.data(),selectedB.data());accepted->SetMarkerStyle(24);
        accepted->SetMarkerSize(.65);accepted->SetMarkerColor(kBlack);accepted->Draw("P SAME");}
    std::vector<double> dbscanA,dbscanB;if(showDbscan)for(const auto& h:dbscanHits)if(h.projection==requiredProjection){double c[3]={h.x,h.y,h.z};dbscanA.push_back(c[a]);dbscanB.push_back(c[b]);}
    if(showDbscan&&!dbscanA.empty()){auto* accepted=new TGraph(dbscanA.size(),dbscanA.data(),dbscanB.data());accepted->SetMarkerStyle(24);
        accepted->SetMarkerSize(.8);accepted->SetMarkerColor(kMagenta+2);accepted->SetLineColor(kMagenta+2);accepted->Draw("P SAME");}
    double p[3]={fit.fitX,fit.fitY,fit.fitZ},d[3]={fit.fitDx,fit.fitDy,fit.fitDz};
    double t1=-2500,t2=2500,la[2]={p[a]+t1*d[a],p[a]+t2*d[a]},lb[2]={p[b]+t1*d[b],p[b]+t2*d[b]};
    auto* line=new TGraph(2,la,lb);line->SetLineColor(kRed);line->SetLineWidth(3);line->Draw("L SAME");
    for(const auto& s:segments){double c1[3]={s.x1,s.y1,s.z1},c2[3]={s.x2,s.y2,s.z2};double ma[2]={c1[a],c2[a]},mb[2]={c1[b],c2[b]};
        auto* mc=new TGraph(2,ma,mb);mc->SetLineColor(kGreen+2);mc->SetLineWidth(2);mc->Draw("L SAME");}
    auto drawRl=[&](bool best,int colour,int style){std::vector<double> ra,rb;for(const auto& p:rlPoints){double c[3]={best?p.bestX:p.x,best?p.bestY:p.y,best?p.bestZ:p.z};ra.push_back(c[a]);rb.push_back(c[b]);}
        if(!ra.empty()){auto* graph=new TGraph(ra.size(),ra.data(),rb.data());graph->SetLineColor(colour);graph->SetMarkerColor(colour);graph->SetLineWidth(3);graph->SetLineStyle(style);graph->SetMarkerStyle(20);graph->SetMarkerSize(.45);graph->Draw("LP SAME");}};
    if(rlSolution=="regularized"||rlSolution=="both")drawRl(false,kMagenta+1,1);
    if(rlSolution=="best-nll"||rlSolution=="both")drawRl(true,kOrange+7,rlSolution=="both"?2:1);
}

void DrawEvent(TCanvas& canvas,TTree& hitsTree,TTree& fitTree,TTree* segmentTree,TTree* selectedTree,TTree* dbscanTree,TTree* rlTree,int event,const std::string& colourMode,const std::string& rlSolution) {
    GlobalFitResult fit;if(!ReadFit(fitTree,event,fit)){std::cerr<<"No fit for event "<<event<<"\n";return;}
    const bool forceAll=colourMode=="all-log"||colourMode=="all-linear";
    const std::string drawMode=colourMode=="all-log"?"log":colourMode=="all-linear"?"linear":colourMode;
    auto hits=ReadHits(hitsTree,event,fit.minimumCharge);auto selectedHits=selectedTree?ReadHits(*selectedTree,event,0):std::vector<Hit>();
    auto dbscanHits=dbscanTree?ReadHits(*dbscanTree,event,0):std::vector<Hit>();
    const bool showSelected=!forceAll&&drawMode!="all"&&drawMode!="off"&&drawMode!="selected"&&drawMode!="dbscan"&&!selectedHits.empty()&&selectedHits.size()<hits.size();
    const bool showDbscan=drawMode=="dbscan"&&!dbscanHits.empty();
    if(drawMode=="selected") hits=selectedHits;
    auto segments=segmentTree?ReadSegments(*segmentTree,event):std::vector<Segment>();auto rlPoints=rlTree?ReadRlPoints(*rlTree,event):std::vector<RlPoint>();canvas.Clear();canvas.Divide(2,2);
    const char* title=showDbscan?"charge, DBSCAN=magenta circles, fit=red, MC=green":"charge, selected=black circles, fit=red, MC=green";
    DrawProjection(1,canvas,hits,selectedHits,showSelected,dbscanHits,showDbscan,segments,rlPoints,rlSolution,fit,0,1,2,(std::string("XY: ")+title).c_str(),"x","y",drawMode);
    DrawProjection(2,canvas,hits,selectedHits,showSelected,dbscanHits,showDbscan,segments,rlPoints,rlSolution,fit,0,2,0,(std::string("XZ: ")+title).c_str(),"x","z",drawMode);
    DrawProjection(3,canvas,hits,selectedHits,showSelected,dbscanHits,showDbscan,segments,rlPoints,rlSolution,fit,1,2,1,(std::string("YZ: ")+title).c_str(),"y","z",drawMode);
    canvas.cd(4);std::vector<double>x,y,z;for(auto& h:hits){x.push_back(h.x);y.push_back(h.y);z.push_back(h.z);}
    auto* points=new TGraph2D(x.size(),x.data(),y.data(),z.data());points->SetTitle("Fibre representatives and fit;x [mm];y [mm];z [mm]");
    points->SetMarkerStyle(20);points->SetMarkerSize(drawMode=="off"||drawMode=="all"?.35:0.0);points->SetMarkerColor(kBlue+1);points->Draw("P0");
    if(drawMode!="off"&&drawMode!="all"&&!hits.empty()) {
        double qmin=hits.front().q,qmax=qmin;for(const auto& h:hits){qmin=std::min(qmin,h.q);qmax=std::max(qmax,h.q);}
        constexpr int groups=20;std::vector<TPolyMarker3D*> markers(groups,nullptr);for(int g=0;g<groups;++g){markers[g]=new TPolyMarker3D();markers[g]->SetMarkerStyle(20);markers[g]->SetMarkerSize(.45);markers[g]->SetMarkerColor(gStyle->GetColorPalette(g*(gStyle->GetNumberOfColors()-1)/(groups-1)));}
        const double lo=drawMode=="log"?std::log(std::max(qmin,1.0)):qmin,hi=drawMode=="log"?std::log(std::max(qmax,1.0)):qmax;
        for(const auto& h:hits){double value=drawMode=="log"?std::log(std::max(h.q,1.0)):h.q;int group=hi>lo?std::clamp(static_cast<int>((value-lo)/(hi-lo)*groups),0,groups-1):0;
            markers[group]->SetNextPoint(h.x,h.y,h.z);}for(auto* marker:markers)if(marker->GetN()>0)marker->Draw();
    }
    if(showSelected){auto* accepted3=new TPolyMarker3D();accepted3->SetMarkerStyle(24);accepted3->SetMarkerSize(.6);accepted3->SetMarkerColor(kBlack);
        for(const auto& h:selectedHits) accepted3->SetNextPoint(h.x,h.y,h.z);
        accepted3->Draw();}
    if(showDbscan){auto* accepted3=new TPolyMarker3D();accepted3->SetMarkerStyle(24);accepted3->SetMarkerSize(.8);accepted3->SetMarkerColor(kMagenta+2);
        for(const auto& h:dbscanHits) accepted3->SetNextPoint(h.x,h.y,h.z);
        accepted3->Draw();}
    double px[2],py[2],pz[2];for(int i=0;i<2;++i){double t=i?2500:-2500;px[i]=fit.fitX+t*fit.fitDx;py[i]=fit.fitY+t*fit.fitDy;pz[i]=fit.fitZ+t*fit.fitDz;}
    auto* line3=new TPolyLine3D(2,px,py,pz);line3->SetLineColor(kRed);line3->SetLineWidth(4);line3->Draw();
    for(const auto& s:segments){double mx[2]={s.x1,s.x2},my[2]={s.y1,s.y2},mz[2]={s.z1,s.z2};auto* mc=new TPolyLine3D(2,mx,my,mz);
        mc->SetLineColor(kGreen+2);mc->SetLineWidth(2);mc->Draw();}
    auto drawRl3=[&](bool best,int colour,int style){auto* line=new TPolyLine3D();for(const auto& p:rlPoints)line->SetNextPoint(best?p.bestX:p.x,best?p.bestY:p.y,best?p.bestZ:p.z);line->SetLineColor(colour);line->SetLineWidth(4);line->SetLineStyle(style);line->Draw();};
    if(rlSolution=="regularized"||rlSolution=="both")drawRl3(false,kMagenta+1,1);
    if(rlSolution=="best-nll"||rlSolution=="both")drawRl3(true,kOrange+7,rlSolution=="both"?2:1);
    canvas.SetTitle(("Global fit event "+std::to_string(event)).c_str());canvas.Modified();canvas.Update();gSystem->ProcessEvents();
    std::cout<<"event "<<event<<", status "<<fit.status<<", NLL "<<fit.nll
             <<", EDM "<<fit.edm<<", calls "<<fit.functionCalls<<"/"<<fit.maximumFunctionCalls
             <<", chi2/ndof "<<fit.chi2<<"/"<<fit.ndof<<" = "
             <<(fit.ndof>0?fit.chi2/fit.ndof:0.0)
             <<", selected fibres "<<fit.observationsBeforeClustering<<" -> "<<fit.observationsAfterDbscan<<" -> "<<fit.observations
             <<(fit.dbscanEnabled?" (DBSCAN)":" (no DBSCAN)")
             <<(fit.corridorEnabled?" + corridor":"")<<"\n";
}
}

int main(int argc,char** argv){
    auto usage=[](){std::cout<<R"(Usage:
  global_fit_display FITTED.root [KEY=value ...]

Options:
  EVENT=first|NUMBER   Initial fitted event (default: first)
  MODE=log             Charge colours plus selected-fibre outline
  MODE=linear          Linear charge colours plus selection outline
  MODE=all             All input fibres, uniform colour, no selection
  MODE=all-log         All input fibres, log charge colours, no selection
  MODE=all-linear      All input fibres, linear charge colours, no selection
  MODE=selected        Only fibres used by the fit
  MODE=dbscan          All threshold-passing fibres, DBSCAN-retained outlined
  RL=FILE.root         Overlay rl_columns from an RL-refinement output
  RL_SOLUTION=regularized|best-nll|both  RL state to draw (default: both)

Overlay colours: MC truth=green, global straight fit=red,
regularized RL=magenta, best-raw-NLL RL=orange dashed.
)";};
    if(argc==1||(argc==2&&std::string(argv[1])=="--help")){usage();return 0;}
    if(argc<2){usage();return 2;}
    // TApplication consumes ROOT-style command-line arguments and may rewrite
    // argc/argv. Preserve this program's positional arguments first and give
    // ROOT an independent minimal argument list.
    const std::string inputName=argv[1];
    std::string eventOption="first",colourMode="log",rlName,rlSolution="both";
    for(int i=2;i<argc;++i){std::string argument=argv[i];auto equal=argument.find('=');if(equal==std::string::npos){std::cerr<<"Options must use KEY=value: "<<argument<<"\n";return 2;}
        std::string key=argument.substr(0,equal),value=argument.substr(equal+1);if(key=="EVENT")eventOption=value;else if(key=="MODE")colourMode=value;else if(key=="RL")rlName=value;else if(key=="RL_SOLUTION")rlSolution=value;else{std::cerr<<"Unknown option: "<<key<<"\n";return 2;}}
    const bool requestedEvent=eventOption!="first";
    const int initialEvent=requestedEvent?std::stoi(eventOption):0;
    if(colourMode!="log"&&colourMode!="linear"&&colourMode!="off"&&colourMode!="all"&&colourMode!="all-log"&&colourMode!="all-linear"&&colourMode!="selected"&&colourMode!="dbscan") {std::cerr<<"Invalid display MODE\n";return 2;}
    if(rlSolution!="regularized"&&rlSolution!="best-nll"&&rlSolution!="both"){std::cerr<<"Invalid RL_SOLUTION\n";return 2;}
    int rootArgc=1;
    char applicationName[]="global_fit_display";
    char* rootArgv[]={applicationName,nullptr};
    TApplication app("global_fit_display",&rootArgc,rootArgv);
    TFile file(inputName.c_str(),"READ");
    auto* fits=dynamic_cast<TTree*>(file.Get("global_fit"));if(!fits){std::cerr<<"Missing global_fit tree\n";return 1;}
    std::set<int> available;int ev=0;fits->SetBranchAddress("event",&ev);for(Long64_t i=0;i<fits->GetEntries();++i){fits->GetEntry(i);available.insert(ev);}
    if(available.empty()){std::cerr<<"global_fit is empty\n";return 1;}int current=requestedEvent?initialEvent:*available.begin();
    GlobalFitResult first;if(!ReadFit(*fits,current,first)){current=*available.begin();ReadFit(*fits,current,first);}
    auto* hits=dynamic_cast<TTree*>(file.Get(first.inputTree.c_str()));if(!hits){std::cerr<<"Missing "<<first.inputTree<<" tree\n";return 1;}
    auto* segments=dynamic_cast<TTree*>(file.Get("mc_virtual_segments"));
    auto* selected=dynamic_cast<TTree*>(file.Get("global_fit_fibres"));
    auto* dbscan=dynamic_cast<TTree*>(file.Get("global_fit_dbscan_fibres"));
    std::unique_ptr<TFile> rlFile;TTree* rlTree=nullptr;if(!rlName.empty()){rlFile=std::make_unique<TFile>(rlName.c_str(),"READ");if(rlFile->IsZombie()){std::cerr<<"Cannot open RL file "<<rlName<<"\n";return 1;}rlTree=dynamic_cast<TTree*>(rlFile->Get("rl_columns"));if(!rlTree){std::cerr<<"Missing rl_columns in "<<rlName<<"\n";return 1;}}
    if(colourMode=="dbscan"&&!dbscan){std::cerr<<"Missing global_fit_dbscan_fibres; rerun global_light_fit with the updated executable\n";return 1;}
    TCanvas canvas("global_fit_display","Global light fit display",1400,900);DrawEvent(canvas,*hits,*fits,segments,selected,dbscan,rlTree,current,colourMode,rlSolution);
    std::cout<<"Commands: n=next, p=previous, EVENT=jump, q=quit\n";std::string command;
    while(std::cout<<"> "&&std::cin>>command&&command!="q"){
        auto it=available.find(current);if(command=="n"){if(++it!=available.end())current=*it;}
        else if(command=="p"){if(it!=available.begin())current=*--it;}
        else try{int requested=std::stoi(command);if(available.count(requested))current=requested;else std::cerr<<"No fit for event "<<requested<<"\n";}catch(...){std::cerr<<"Unknown command\n";}
        DrawEvent(canvas,*hits,*fits,segments,selected,dbscan,rlTree,current,colourMode,rlSolution);
    }return 0;
}
