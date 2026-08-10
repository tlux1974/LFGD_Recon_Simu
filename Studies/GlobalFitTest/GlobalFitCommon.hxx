#ifndef GLOBAL_FIT_COMMON_HXX
#define GLOBAL_FIT_COMMON_HXX

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TKey.h>
#include <TObject.h>
#include <TTree.h>

struct GlobalFitResult {
    int event = 0;
    int status = -1;
    int observations = 0;
    int observationsBeforeClustering = 0;
    int observationsAfterDbscan = 0;
    int dbscanEnabled = 0;
    double dbscanEpsilon = 14.2;
    int dbscanMinPoints = 2;
    int corridorEnabled = 0;
    int corridorHalfWidth = 1;
    double minimumMapFraction = 0.0;
    double minimumCharge = 10.0;
    double nll = 0.0;
    double chi2 = 0.0;
    int ndof = 0;
    double seedX = 0.0, seedY = 0.0, seedZ = 0.0;
    double seedDx = 0.0, seedDy = 0.0, seedDz = 1.0;
    double fitX = 0.0, fitY = 0.0, fitZ = 0.0;
    double fitDx = 0.0, fitDy = 0.0, fitDz = 1.0;
    double errorA = 0.0, errorB = 0.0, errorTheta = 0.0, errorPhi = 0.0;
    std::string inputTree;
};

inline void BranchGlobalFit(TTree& tree, GlobalFitResult& r) {
    tree.Branch("event", &r.event);
    tree.Branch("status", &r.status);
    tree.Branch("observations", &r.observations);
    tree.Branch("observations_before_clustering", &r.observationsBeforeClustering);
    tree.Branch("observations_after_dbscan", &r.observationsAfterDbscan);
    tree.Branch("dbscan_enabled", &r.dbscanEnabled);
    tree.Branch("dbscan_epsilon", &r.dbscanEpsilon);
    tree.Branch("dbscan_min_points", &r.dbscanMinPoints);
    tree.Branch("corridor_enabled", &r.corridorEnabled);
    tree.Branch("corridor_half_width", &r.corridorHalfWidth);
    tree.Branch("minimum_map_fraction", &r.minimumMapFraction);
    tree.Branch("minimum_charge", &r.minimumCharge);
    tree.Branch("nll", &r.nll);
    tree.Branch("chi2", &r.chi2);
    tree.Branch("ndof", &r.ndof);
    tree.Branch("seed_x", &r.seedX); tree.Branch("seed_y", &r.seedY);
    tree.Branch("seed_z", &r.seedZ);
    tree.Branch("seed_dx", &r.seedDx); tree.Branch("seed_dy", &r.seedDy);
    tree.Branch("seed_dz", &r.seedDz);
    tree.Branch("fit_x", &r.fitX); tree.Branch("fit_y", &r.fitY);
    tree.Branch("fit_z", &r.fitZ);
    tree.Branch("fit_dx", &r.fitDx); tree.Branch("fit_dy", &r.fitDy);
    tree.Branch("fit_dz", &r.fitDz);
    tree.Branch("error_a", &r.errorA); tree.Branch("error_b", &r.errorB);
    tree.Branch("error_theta", &r.errorTheta); tree.Branch("error_phi", &r.errorPhi);
    tree.Branch("input_tree", &r.inputTree);
}

inline void SetGlobalFitAddresses(TTree& tree, GlobalFitResult& r) {
#define GF_ADDRESS(name, member) do { if (tree.GetBranch(name)) tree.SetBranchAddress(name, &r.member); } while(false)
    GF_ADDRESS("event", event); GF_ADDRESS("status", status);
    GF_ADDRESS("observations", observations); GF_ADDRESS("minimum_charge", minimumCharge);
    GF_ADDRESS("observations_before_clustering", observationsBeforeClustering);
    GF_ADDRESS("observations_after_dbscan", observationsAfterDbscan);
    GF_ADDRESS("dbscan_enabled", dbscanEnabled); GF_ADDRESS("dbscan_epsilon", dbscanEpsilon);
    GF_ADDRESS("dbscan_min_points", dbscanMinPoints);
    GF_ADDRESS("corridor_enabled", corridorEnabled);
    GF_ADDRESS("corridor_half_width", corridorHalfWidth);
    GF_ADDRESS("minimum_map_fraction", minimumMapFraction);
    GF_ADDRESS("nll", nll);
    GF_ADDRESS("chi2", chi2); GF_ADDRESS("ndof", ndof);
    GF_ADDRESS("seed_x", seedX); GF_ADDRESS("seed_y", seedY); GF_ADDRESS("seed_z", seedZ);
    GF_ADDRESS("seed_dx", seedDx); GF_ADDRESS("seed_dy", seedDy); GF_ADDRESS("seed_dz", seedDz);
    GF_ADDRESS("fit_x", fitX); GF_ADDRESS("fit_y", fitY); GF_ADDRESS("fit_z", fitZ);
    GF_ADDRESS("fit_dx", fitDx); GF_ADDRESS("fit_dy", fitDy); GF_ADDRESS("fit_dz", fitDz);
    GF_ADDRESS("error_a", errorA); GF_ADDRESS("error_b", errorB);
    GF_ADDRESS("error_theta", errorTheta); GF_ADDRESS("error_phi", errorPhi);
    // ROOT object branches such as std::string require a pointer-to-pointer
    // address. Readers that need input_tree attach it locally.
#undef GF_ADDRESS
}

// Copy the newest cycle of every top-level object. Trees are fast-cloned so
// the result is self-contained without deserializing all entries.
inline void CloneFlatFile(TFile& input, TFile& output) {
    TIter next(input.GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        if (key->GetCycle() != input.GetKey(key->GetName())->GetCycle()) continue;
        if (std::string(key->GetName()) == "global_fit") continue;
        TObject* object = key->ReadObj();
        output.cd();
        if (auto* tree = dynamic_cast<TTree*>(object)) {
            TTree* clone = tree->CloneTree(-1, "fast");
            clone->Write(tree->GetName(), TObject::kOverwrite);
            delete clone;
        }
        else object->Write(object->GetName(), TObject::kOverwrite);
        delete object;
    }
}

#endif
