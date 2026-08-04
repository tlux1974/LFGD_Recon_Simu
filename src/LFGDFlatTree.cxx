#include <iostream>
#include <string>

#include <TAlgorithmResult.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <TND280Event.hxx>
#include <TND280Log.hxx>
#include <TParametersOptionManager.hxx>
#include <TReconNode.hxx>
#include <TReconTrack.hxx>
#include <TTrackState.hxx>
#include <TG4Trajectory.hxx>
#include <ND280GeomId.hxx>
#include <nd280EventLoop.hxx>

#include <TFile.h>
#include <TTree.h>

class LFGDFlatTree : public ND::TND280EventLoopFunction {
public:
    LFGDFlatTree() : fEvent(0), fOutput(nullptr), fFiberTree(nullptr),
                     fHitTree(nullptr), fTrackTree(nullptr), fMCTrackTree(nullptr) {}

    bool SetOption(std::string option, std::string value = "") override {
        if (option == "outfile") { fOutputName = value; return true; }
        if (fOptions.IsRelevantOption(option)) {
            fOptions.UseRelevantOption(value); return true;
        }
        return false;
    }

    void Usage() override {
        std::cout << "Create flat fibre_hits, hits3d, track_nodes, and mc_track_points trees.\n"
                  << "  -O outfile=flat.root\n";
    }

    void Initialize() override {
        if (fOutputName.empty()) throw std::runtime_error("Missing -O outfile=...");
        fOutput = TFile::Open(fOutputName.c_str(), "RECREATE");
        if (!fOutput || fOutput->IsZombie()) throw std::runtime_error("Cannot create output");

        fFiberTree = new TTree("fiber_hits", "Calibrated fibre hits");
        BranchCommon(fFiberTree);
        fFiberTree->Branch("projection", &fProjection);
        fFiberTree->Branch("u", &fU);
        fFiberTree->Branch("v", &fV);

        fHitTree = new TTree("hits3d", "Reconstructed 3D hits");
        BranchCommon(fHitTree);

        fTrackTree = new TTree("track_nodes", "Reconstructed HFG track nodes");
        fTrackTree->Branch("event", &fEvent);
        fTrackTree->Branch("track", &fTrack);
        fTrackTree->Branch("node", &fNode);
        fTrackTree->Branch("x", &fX);
        fTrackTree->Branch("y", &fY);
        fTrackTree->Branch("z", &fZ);
        fTrackTree->Branch("t", &fTime);

        fMCTrackTree = new TTree("mc_track_points", "Geant4 truth trajectory points");
        fMCTrackTree->Branch("event", &fEvent);
        fMCTrackTree->Branch("track_id", &fMCTrackId);
        fMCTrackTree->Branch("parent_id", &fMCParentId);
        fMCTrackTree->Branch("pdg", &fMCPdg);
        fMCTrackTree->Branch("particle", &fMCParticle);
        fMCTrackTree->Branch("point", &fMCPoint);
        fMCTrackTree->Branch("x", &fX);
        fMCTrackTree->Branch("y", &fY);
        fMCTrackTree->Branch("z", &fZ);
        fMCTrackTree->Branch("t", &fTime);
        fMCTrackTree->Branch("px", &fPx);
        fMCTrackTree->Branch("py", &fPy);
        fMCTrackTree->Branch("pz", &fPz);
    }

    bool operator()(ND::TND280Event& event) override {
        const char* inputName = event.GetHitSelection("homo") ? "homo" : "hfg";
        auto fibers = event.GetHitSelection(inputName);
        if (fibers) for (const auto& hit : *fibers) {
            fX = hit->GetPosition().X(); fY = hit->GetPosition().Y();
            fZ = hit->GetPosition().Z(); fTime = hit->GetTime();
            fCharge = hit->GetCharge(); fGeomId = hit->GetGeomId().AsInt();
            if (std::string(inputName) == "homo") {
                fProjection = ND::GeomId::Homo::GetFiberDirection(hit->GetGeomId());
                fU = ND::GeomId::Homo::GetFiberU(hit->GetGeomId());
                fV = ND::GeomId::Homo::GetFiberV(hit->GetGeomId());
            } else {
                fProjection = ND::GeomId::HFG::GetFiberProjection(hit->GetGeomId());
                fU = ND::GeomId::HFG::GetFiberU(hit->GetGeomId());
                fV = ND::GeomId::HFG::GetFiberV(hit->GetGeomId());
            }
            fFiberTree->Fill();
        }

        auto hits3d = event.GetHitSelection("hfg_3d");
        if (hits3d) for (const auto& hit : *hits3d) {
            bool is3d = ND::GeomId::Homo::IsCube(hit->GetGeomId())
                     || ND::GeomId::HFG::IsCube(hit->GetGeomId());
            if (!is3d) continue;
            fX = hit->GetPosition().X(); fY = hit->GetPosition().Y();
            fZ = hit->GetPosition().Z(); fTime = hit->GetTime();
            fCharge = hit->GetCharge(); fGeomId = hit->GetGeomId().AsInt();
            fHitTree->Fill();
        }

        auto result = event.GetFit("THFGRecon");
        if (!result) result = event.GetFit("hfgRecon");
        auto objects = result ? result->GetResultsContainer("final")
                              : ND::THandle<ND::TReconObjectContainer>();
        fTrack = 0;
        if (objects) for (const auto& object : *objects) {
            ND::THandle<ND::TReconTrack> track = object;
            // The final THFGRecon objects are normally TReconPID wrappers.
            // Recover the reconstructed track stored as their constituent.
            if (!track) {
                auto constituents = object->GetConstituents();
                if (constituents) for (const auto& constituent : *constituents) {
                    track = constituent;
                    if (track) break;
                }
            }
            if (!track) continue;
            fNode = 0;
            for (const auto& node : track->GetNodes()) {
                ND::THandle<ND::TTrackState> state = node->GetState();
                if (!state) continue;
                auto position = state->GetPosition();
                fX = position.X(); fY = position.Y(); fZ = position.Z();
                fTime = position.T();
                fTrackTree->Fill();
                ++fNode;
            }
            ++fTrack;
        }

        auto trajectories = event.Get<ND::TG4TrajectoryContainer>(
            "truth/G4Trajectories");
        if (trajectories) for (const auto& entry : *trajectories) {
            const auto& trajectory = entry.second;
            fMCTrackId = trajectory.GetTrackId();
            fMCParentId = trajectory.GetParentId();
            fMCPdg = trajectory.GetPDGEncoding();
            fMCParticle = trajectory.GetParticleName();
            fMCPoint = 0;
            for (const auto& point : trajectory.GetTrajectoryPoints()) {
                auto position = point.GetPosition();
                auto momentum = point.GetMomentum();
                fX = position.X(); fY = position.Y(); fZ = position.Z();
                fTime = position.T();
                fPx = momentum.X(); fPy = momentum.Y(); fPz = momentum.Z();
                fMCTrackTree->Fill();
                ++fMCPoint;
            }
        }
        ++fEvent;
        return true;
    }

    void Finalize(ND::TND280Output* const) override {
        fOutput->cd();
        fFiberTree->Write(); fHitTree->Write(); fTrackTree->Write();
        fMCTrackTree->Write();
        fOutput->Close();
    }

private:
    void BranchCommon(TTree* tree) {
        tree->Branch("event", &fEvent); tree->Branch("x", &fX);
        tree->Branch("y", &fY); tree->Branch("z", &fZ);
        tree->Branch("time", &fTime); tree->Branch("charge", &fCharge);
        tree->Branch("geom_id", &fGeomId);
    }

    ND::TParametersOptionManager fOptions;
    std::string fOutputName;
    int fEvent, fTrack, fNode, fProjection, fU, fV;
    int fMCTrackId, fMCParentId, fMCPdg, fMCPoint;
    std::string fMCParticle;
    unsigned int fGeomId;
    double fX, fY, fZ, fTime, fCharge, fPx, fPy, fPz;
    TFile* fOutput;
    TTree *fFiberTree, *fHitTree, *fTrackTree, *fMCTrackTree;
};

int main(int argc, char** argv) {
    LFGDFlatTree loop;
    ND::TND280Log::SetLogLevel(ND::TND280Log::LogLevel);
    nd280EventLoop(argc, argv, loop);
}
