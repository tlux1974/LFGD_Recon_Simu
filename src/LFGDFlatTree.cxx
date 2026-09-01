#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include <TAlgorithmResult.hxx>
#include <THandle.hxx>
#include <THitSelection.hxx>
#include <TMCHit.hxx>
#include <TND280Event.hxx>
#include <TND280Log.hxx>
#include <TParametersOptionManager.hxx>
#include <TReconNode.hxx>
#include <TReconCluster.hxx>
#include <TReconTrack.hxx>
#include <TTrackState.hxx>
#include <TG4Trajectory.hxx>
#include <TG4HitSegment.hxx>
#include <TG4VHit.hxx>
#include <ND280GeomId.hxx>
#include <TGeomInfo.hxx>
#include <THomoGeom.hxx>
#include <THFGGeom.hxx>
#include <nd280EventLoop.hxx>

#include <TFile.h>
#include <TTree.h>

class LFGDFlatTree : public ND::TND280EventLoopFunction {
public:
    LFGDFlatTree() : fEvent(0), fOutput(nullptr), fFiberTree(nullptr),
                     fHomoRawTree(nullptr), fHomoTruthTree(nullptr),
                     fHitTree(nullptr), fHitViewTree(nullptr),
                     fTrackTree(nullptr), fTrackHitTree(nullptr),
                     fMCTrackTree(nullptr),
                     fMCSegmentTree(nullptr) {}

    bool SetOption(std::string option, std::string value = "") override {
        if (option == "outfile") { fOutputName = value; return true; }
        if (fOptions.IsRelevantOption(option)) {
            fOptions.UseRelevantOption(value); return true;
        }
        return false;
    }

    void Usage() override {
        std::cout << "Create flat fibre, 3D-hit, track, trajectory, and per-cube MC segment trees.\n"
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

        fHomoRawTree = new TTree(
            "homo_raw", "HOMO fibre photons before MPPC and electronics");
        BranchFiber(fHomoRawTree);

        fHomoTruthTree = new TTree(
            "homo_truth", "HOMO light-map fibre truth hits");
        BranchFiber(fHomoTruthTree);

        fHitTree = new TTree("hits3d", "Reconstructed 3D hits");
        BranchCommon(fHitTree);

        fHitViewTree = new TTree(
            "hit3d_views", "2D-view and fibre composition of 3D hits");
        fHitViewTree->Branch("event", &fEvent);
        fHitViewTree->Branch("hit", &fHit3D);
        fHitViewTree->Branch("geom_id", &fGeomId);
        fHitViewTree->Branch("hit_charge", &fCharge);
        fHitViewTree->Branch("view", &fProjection);
        fHitViewTree->Branch("view_geom_id", &fViewGeomId);
        fHitViewTree->Branch("view_charge", &fViewCharge);
        fHitViewTree->Branch("fiber_count", &fFiberCount);
        fHitViewTree->Branch("fiber_charge_sum", &fFiberChargeSum);

        fTrackTree = new TTree("track_nodes", "Reconstructed HFG track nodes");
        fTrackTree->Branch("event", &fEvent);
        fTrackTree->Branch("track", &fTrack);
        fTrackTree->Branch("node", &fNode);
        fTrackTree->Branch("x", &fX);
        fTrackTree->Branch("y", &fY);
        fTrackTree->Branch("z", &fZ);
        fTrackTree->Branch("t", &fTime);

        // A track-node state is the fitted position and need not lie in the
        // same voxel as the measurement that created the node.  Preserve the
        // actual node-to-hit relation separately for reconstruction audits.
        fTrackHitTree = new TTree(
            "track_node_hits", "3D hits attached to reconstructed track nodes");
        BranchCommon(fTrackHitTree);
        fTrackHitTree->Branch("track", &fTrack);
        fTrackHitTree->Branch("node", &fNode);
        fTrackHitTree->Branch("node_hit", &fNodeHit);

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

        // Keep the historical tree name for compatibility.  It contains
        // virtual-cube segments for HOMO and physical-cube segments for HFG.
        fMCSegmentTree = new TTree(
            "mc_virtual_segments", "Per-cube Geant4 truth segments");
        fMCSegmentTree->Branch("event", &fEvent);
        fMCSegmentTree->Branch("segment", &fMCSegment);
        fMCSegmentTree->Branch("detector", &fMCDetector);
        fMCSegmentTree->Branch("primary_id", &fMCPrimaryId);
        fMCSegmentTree->Branch("primary_pdg", &fMCPrimaryPdg);
        fMCSegmentTree->Branch("contributors", &fMCContributors);
        fMCSegmentTree->Branch("contributor_track_ids", &fMCContributorTrackIds);
        fMCSegmentTree->Branch("contributor_pdgs", &fMCContributorPdgs);
        fMCSegmentTree->Branch("cube_x", &fMCCubeX);
        fMCSegmentTree->Branch("cube_y", &fMCCubeY);
        fMCSegmentTree->Branch("cube_z", &fMCCubeZ);
        fMCSegmentTree->Branch("geom_id", &fGeomId);
        fMCSegmentTree->Branch("start_x", &fStartX);
        fMCSegmentTree->Branch("start_y", &fStartY);
        fMCSegmentTree->Branch("start_z", &fStartZ);
        fMCSegmentTree->Branch("start_t", &fStartT);
        fMCSegmentTree->Branch("stop_x", &fStopX);
        fMCSegmentTree->Branch("stop_y", &fStopY);
        fMCSegmentTree->Branch("stop_z", &fStopZ);
        fMCSegmentTree->Branch("stop_t", &fStopT);
        fMCSegmentTree->Branch("energy_deposit", &fEnergyDeposit);
        fMCSegmentTree->Branch("track_length", &fMCTrackLength);
    }

    bool operator()(ND::TND280Event& event) override {
        const char* inputName = event.GetHitSelection("homo") ? "homo" : "hfg";
        auto fibers = event.GetHitSelection(inputName);
        if (fibers) for (const auto& hit : *fibers) {
            TVector3 position = hit->GetPosition();
            // HOMO fibre IDs are backed by representative virtual-cube nodes,
            // not physical TGeo fibre nodes.  Obtain their real staggered-grid
            // position from THomoGeom instead of serializing the representative
            // cube position (or the zero default from an older geometry map).
            if (std::string(inputName) == "homo") {
                position = ND::TGeomInfo::Get().HOMO()
                    .GetFiber(hit->GetGeomId()).GetPosition();
            }
            fX = position.X(); fY = position.Y();
            fZ = position.Z(); fTime = hit->GetTime();
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

        FillHomoFiberTree(event.GetHitSelection("homo_raw"), fHomoRawTree);
        FillHomoFiberTree(event.GetHitSelection("homo_truth"), fHomoTruthTree);

        auto hits3d = event.GetHitSelection("hfg_3d");
        fHit3D = 0;
        if (hits3d) for (const auto& hit : *hits3d) {
            bool is3d = ND::GeomId::Homo::IsCube(hit->GetGeomId())
                     || ND::GeomId::HFG::IsCube(hit->GetGeomId());
            if (!is3d) continue;
            fX = hit->GetPosition().X(); fY = hit->GetPosition().Y();
            fZ = hit->GetPosition().Z(); fTime = hit->GetTime();
            fCharge = hit->GetCharge(); fGeomId = hit->GetGeomId().AsInt();
            fHitTree->Fill();

            // Contributor zero is the fake cube seed.  The remaining
            // contributors are the reconstructed 2D view hits.  For HOMO a
            // view hit can itself contain a central fibre plus neighbouring
            // fibres assigned by the local 2D clustering method.
            for (int contributor = 1;
                 contributor < hit->GetContributorCount(); ++contributor) {
                ND::THandle<ND::THit> viewHit
                    = hit->GetContributor(contributor);
                if (!viewHit) continue;
                const auto viewId = viewHit->GetGeomId();
                if (ND::GeomId::Homo::IsFiber(viewId)) {
                    fProjection
                        = ND::GeomId::Homo::GetFiberDirection(viewId);
                }
                else if (ND::GeomId::HFG::IsFiber(viewId)) {
                    fProjection
                        = ND::GeomId::HFG::GetFiberProjection(viewId);
                }
                else continue;
                fViewGeomId = viewId.AsInt();
                fViewCharge = viewHit->GetCharge();
                fFiberCount = 0;
                fFiberChargeSum = 0.0;
                for (int fibre = 0;
                     fibre < viewHit->GetContributorCount(); ++fibre) {
                    ND::THandle<ND::THit> fiberHit
                        = viewHit->GetContributor(fibre);
                    if (!fiberHit) continue;
                    const auto fiberId = fiberHit->GetGeomId();
                    if (!ND::GeomId::Homo::IsFiber(fiberId)
                        && !ND::GeomId::HFG::IsFiber(fiberId)) continue;
                    ++fFiberCount;
                    fFiberChargeSum += fiberHit->GetCharge();
                }
                // A raw fibre contributor has no lower-level constituents.
                if (fFiberCount == 0) {
                    fFiberCount = 1;
                    fFiberChargeSum = viewHit->GetCharge();
                }
                fHitViewTree->Fill();
            }
            ++fHit3D;
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

                ND::THandle<ND::TReconCluster> cluster = node->GetObject();
                auto nodeHits = cluster ? cluster->GetHits()
                                        : ND::THandle<ND::THitSelection>();
                fNodeHit = 0;
                if (nodeHits) for (const auto& hit : *nodeHits) {
                    const auto geomId = hit->GetGeomId();
                    if (!ND::GeomId::Homo::IsCube(geomId)
                        && !ND::GeomId::HFG::IsCube(geomId)) continue;
                    const auto hitPosition = hit->GetPosition();
                    fX = hitPosition.X(); fY = hitPosition.Y();
                    fZ = hitPosition.Z(); fTime = hit->GetTime();
                    fCharge = hit->GetCharge(); fGeomId = geomId.AsInt();
                    fTrackHitTree->Fill();
                    ++fNodeHit;
                }
                ++fNode;
            }
            ++fTrack;
        }

        auto trajectories = event.Get<ND::TG4TrajectoryContainer>(
            "truth/G4Trajectories");
        fMCTrackPdgById.clear();
        if (trajectories) for (const auto& entry : *trajectories) {
            const auto& trajectory = entry.second;
            fMCTrackId = trajectory.GetTrackId();
            fMCParentId = trajectory.GetParentId();
            fMCPdg = trajectory.GetPDGEncoding();
            fMCTrackPdgById[fMCTrackId] = fMCPdg;
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

        fMCSegment = 0;
        if (std::string(inputName) == "homo") {
            auto virtualSegments = event.Get<ND::TG4HitContainer>(
                "truth/g4Hits/homoVirtualCube");
            if (virtualSegments) for (const auto& baseHit : *virtualSegments) {
                const auto* segment
                    = dynamic_cast<const ND::TG4HitSegment*>(baseHit);
                if (!segment) continue;
                const TVector3 midpoint(
                    0.5*(segment->GetStartX()+segment->GetStopX()),
                    0.5*(segment->GetStartY()+segment->GetStopY()),
                    0.5*(segment->GetStartZ()+segment->GetStopZ()));
                const TVector3 local = ND::TGeomInfo::Get().MasterToLocal(
                    ND::GeomId::Homo::Detector(),
                    midpoint.X(),midpoint.Y(),midpoint.Z());
                const ND::THomoGeom& homo = ND::TGeomInfo::Get().HOMO();
                const int cubeX = static_cast<int>(std::floor(
                    (local.X()+5.0*homo.GetXCubes())/10.0));
                const int cubeY = static_cast<int>(std::floor(
                    (local.Y()+5.0*homo.GetYCubes())/10.0));
                const int cubeZ = static_cast<int>(std::floor(
                    (local.Z()+5.0*homo.GetZCubes())/10.0));
                if (cubeX < 0 || cubeX >= homo.GetXCubes()
                    || cubeY < 0 || cubeY >= homo.GetYCubes()
                    || cubeZ < 0 || cubeZ >= homo.GetZCubes()) continue;
                FillMCSegment(*segment, ND::GeomId::Homo::Cube(
                    cubeX,cubeY,cubeZ), cubeX,cubeY,cubeZ, 0);
            }
        }
        else {
            // hfg_truth supplies the exact physical cube identifier for each
            // contributing Geant4 segment.  This avoids a fragile coordinate
            // rounding operation near cube boundaries.
            auto truthHits = event.GetHitSelection("hfg_truth");
            if (truthHits) for (const auto& hit : *truthHits) {
                const auto cubeId = hit->GetGeomId();
                if (!ND::GeomId::HFG::IsCube(cubeId)) continue;
                ND::THandle<ND::TMCHit> truth = hit;
                if (!truth) continue;
                std::set<const ND::TG4HitSegment*> uniqueSegments;
                for (const auto* baseHit : truth->GetContributors()) {
                    const auto* segment
                        = dynamic_cast<const ND::TG4HitSegment*>(baseHit);
                    if (!segment || !uniqueSegments.insert(segment).second)
                        continue;
                    const auto cube = ND::TGeomInfo::Get().HFG().GetCube(cubeId);
                    FillMCSegment(*segment, cubeId, cube.GetCubeX(),
                                  cube.GetCubeY(), cube.GetCubeZ(), 1);
                }
            }
        }
        ++fEvent;
        return true;
    }

    void Finalize(ND::TND280Output* const) override {
        fOutput->cd();
        fFiberTree->Write();
        fHomoRawTree->Write(); fHomoTruthTree->Write();
        fHitTree->Write(); fHitViewTree->Write();
        fTrackTree->Write();
        fTrackHitTree->Write();
        fMCTrackTree->Write();
        fMCSegmentTree->Write();
        fOutput->Close();
    }

private:
    void BranchFiber(TTree* tree) {
        BranchCommon(tree);
        tree->Branch("projection", &fProjection);
        tree->Branch("u", &fU);
        tree->Branch("v", &fV);
    }

    void FillHomoFiberTree(ND::THandle<ND::THitSelection> hits,
                           TTree* tree) {
        if (!hits) return;
        for (const auto& hit : *hits) {
            const auto geomId = hit->GetGeomId();
            if (!ND::GeomId::Homo::IsFiber(geomId)) continue;
            const TVector3 position
                = ND::TGeomInfo::Get().HOMO().GetFiber(geomId).GetPosition();
            fX = position.X(); fY = position.Y(); fZ = position.Z();
            fTime = hit->GetTime(); fCharge = hit->GetCharge();
            fGeomId = geomId.AsInt();
            fProjection = ND::GeomId::Homo::GetFiberDirection(geomId);
            fU = ND::GeomId::Homo::GetFiberU(geomId);
            fV = ND::GeomId::Homo::GetFiberV(geomId);
            tree->Fill();
        }
    }

    void FillMCSegment(const ND::TG4HitSegment& segment,
                       ND::TGeometryId cubeId,
                       int cubeX, int cubeY, int cubeZ, int detector) {
        fMCDetector = detector;
        fMCCubeX = cubeX; fMCCubeY = cubeY; fMCCubeZ = cubeZ;
        fMCPrimaryId = segment.GetPrimaryId();
        auto primaryPdg=fMCTrackPdgById.find(fMCPrimaryId);
        fMCPrimaryPdg=primaryPdg==fMCTrackPdgById.end()?0:primaryPdg->second;
        fMCContributorTrackIds.clear();fMCContributorPdgs.clear();
        for(const auto contributor:segment.GetContributors()){
            fMCContributorTrackIds.push_back(contributor);auto pdg=fMCTrackPdgById.find(contributor);
            fMCContributorPdgs.push_back(pdg==fMCTrackPdgById.end()?0:pdg->second);}
        fMCContributors = fMCContributorTrackIds.size();
        fGeomId = cubeId.AsInt();
        fStartX = segment.GetStartX(); fStartY = segment.GetStartY();
        fStartZ = segment.GetStartZ(); fStartT = segment.GetStartT();
        fStopX = segment.GetStopX(); fStopY = segment.GetStopY();
        fStopZ = segment.GetStopZ(); fStopT = segment.GetStopT();
        fEnergyDeposit = segment.GetEnergyDeposit();
        fMCTrackLength = segment.GetTrackLength();
        fMCSegmentTree->Fill();
        ++fMCSegment;
    }

    void BranchCommon(TTree* tree) {
        tree->Branch("event", &fEvent); tree->Branch("x", &fX);
        tree->Branch("y", &fY); tree->Branch("z", &fZ);
        tree->Branch("time", &fTime); tree->Branch("charge", &fCharge);
        tree->Branch("geom_id", &fGeomId);
    }

    ND::TParametersOptionManager fOptions;
    std::string fOutputName;
    int fEvent, fTrack, fNode, fNodeHit, fHit3D, fProjection, fU, fV;
    int fFiberCount;
    int fMCTrackId, fMCParentId, fMCPdg, fMCPoint;
    int fMCSegment, fMCDetector, fMCPrimaryId, fMCPrimaryPdg, fMCContributors;
    int fMCCubeX, fMCCubeY, fMCCubeZ;
    std::string fMCParticle;
    std::map<int,int> fMCTrackPdgById;
    std::vector<int> fMCContributorTrackIds,fMCContributorPdgs;
    unsigned int fGeomId, fViewGeomId;
    double fX, fY, fZ, fTime, fCharge, fPx, fPy, fPz;
    double fViewCharge, fFiberChargeSum;
    double fStartX, fStartY, fStartZ, fStartT;
    double fStopX, fStopY, fStopZ, fStopT;
    double fEnergyDeposit, fMCTrackLength;
    TFile* fOutput;
    TTree *fFiberTree, *fHomoRawTree, *fHomoTruthTree;
    TTree *fHitTree, *fHitViewTree, *fTrackTree;
    TTree *fTrackHitTree, *fMCTrackTree;
    TTree *fMCSegmentTree;
};

int main(int argc, char** argv) {
    LFGDFlatTree loop;
    ND::TND280Log::SetLogLevel(ND::TND280Log::LogLevel);
    nd280EventLoop(argc, argv, loop);
}
