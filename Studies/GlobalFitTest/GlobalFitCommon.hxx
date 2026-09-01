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
    int fitRangeEnabled = 1;
    double fitRangeQuantile = 0.01;
    double fitRangePadding = 10.0;
    double fitLowX = 0.0, fitLowY = 0.0, fitLowZ = 0.0;
    double fitHighX = 0.0, fitHighY = 0.0, fitHighZ = 0.0;
    std::string fitViews = "ALL";
    double nll = 0.0;
    double edm = 0.0;
    int functionCalls = 0;
    int maximumFunctionCalls = 5000;
    double tolerance = 1e-2;
    double chi2 = 0.0;
    int ndof = 0;
    double seedX = 0.0, seedY = 0.0, seedZ = 0.0;
    double seedDx = 0.0, seedDy = 0.0, seedDz = 1.0;
    std::string seedMethod = "FIBRE";
    int seedObservations = 0;
    double seedMedianFactor = 1.0;
    double seedMedianXY = 0.0, seedMedianXZ = 0.0, seedMedianYZ = 0.0;
    double seedNll = 0.0, seedChi2 = 0.0;
    int seedNdof = 0;
    double fitX = 0.0, fitY = 0.0, fitZ = 0.0;
    double fitDx = 0.0, fitDy = 0.0, fitDz = 1.0;
    double errorA = 0.0, errorB = 0.0, errorTheta = 0.0, errorPhi = 0.0;
    int columnFitAvailable = 0, columnFitConverged = 0, columnFitIterations = 0, columnFitActiveColumns = 0;
    double columnFitMaxRelativeChange = 0.0, columnNll = 0.0, columnChi2 = 0.0;
    int columnNdof = 0;
    int columnMipPriorEnabled = 0, columnMipPriorConverged = 0, columnMipPriorStalled = 0, columnMipPriorIterations = 0, columnMipPriorEligibleColumns = 0;
    double columnMipPriorFloorEnergy = 0.0, columnMipPriorMpvEnergy = 0.0, columnMipPriorStrength = 0.0;
    double columnMipPriorChargePerEnergy = 0.0, columnMipPriorMinPathFraction = 0.0;
    double columnMipPriorPenalty = 0.0, columnMipPriorObjective = 0.0, columnMipPriorMaxRelativeChange = 0.0;
    double columnUnregularizedNll = 0.0;
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
    tree.Branch("fit_range_enabled", &r.fitRangeEnabled);
    tree.Branch("fit_range_quantile", &r.fitRangeQuantile);
    tree.Branch("fit_range_padding", &r.fitRangePadding);
    tree.Branch("fit_low_x", &r.fitLowX); tree.Branch("fit_low_y", &r.fitLowY); tree.Branch("fit_low_z", &r.fitLowZ);
    tree.Branch("fit_high_x", &r.fitHighX); tree.Branch("fit_high_y", &r.fitHighY); tree.Branch("fit_high_z", &r.fitHighZ);
    tree.Branch("fit_views", &r.fitViews);
    tree.Branch("nll", &r.nll);
    tree.Branch("edm", &r.edm);
    tree.Branch("function_calls", &r.functionCalls);
    tree.Branch("maximum_function_calls", &r.maximumFunctionCalls);
    tree.Branch("tolerance", &r.tolerance);
    tree.Branch("chi2", &r.chi2);
    tree.Branch("ndof", &r.ndof);
    tree.Branch("seed_x", &r.seedX); tree.Branch("seed_y", &r.seedY);
    tree.Branch("seed_z", &r.seedZ);
    tree.Branch("seed_dx", &r.seedDx); tree.Branch("seed_dy", &r.seedDy);
    tree.Branch("seed_dz", &r.seedDz);
    tree.Branch("seed_method", &r.seedMethod);
    tree.Branch("seed_observations", &r.seedObservations);
    tree.Branch("seed_median_factor", &r.seedMedianFactor);
    tree.Branch("seed_median_xy", &r.seedMedianXY);
    tree.Branch("seed_median_xz", &r.seedMedianXZ);
    tree.Branch("seed_median_yz", &r.seedMedianYZ);
    tree.Branch("seed_nll", &r.seedNll);
    tree.Branch("seed_chi2", &r.seedChi2);
    tree.Branch("seed_ndof", &r.seedNdof);
    tree.Branch("fit_x", &r.fitX); tree.Branch("fit_y", &r.fitY);
    tree.Branch("fit_z", &r.fitZ);
    tree.Branch("fit_dx", &r.fitDx); tree.Branch("fit_dy", &r.fitDy);
    tree.Branch("fit_dz", &r.fitDz);
    tree.Branch("error_a", &r.errorA); tree.Branch("error_b", &r.errorB);
    tree.Branch("error_theta", &r.errorTheta); tree.Branch("error_phi", &r.errorPhi);
    tree.Branch("column_fit_available", &r.columnFitAvailable);
    tree.Branch("column_fit_converged", &r.columnFitConverged);
    tree.Branch("column_fit_iterations", &r.columnFitIterations);
    tree.Branch("column_fit_active_columns", &r.columnFitActiveColumns);
    tree.Branch("column_fit_max_relative_change", &r.columnFitMaxRelativeChange);
    tree.Branch("column_nll", &r.columnNll); tree.Branch("column_chi2", &r.columnChi2); tree.Branch("column_ndof", &r.columnNdof);
    tree.Branch("column_mip_prior_enabled", &r.columnMipPriorEnabled);
    tree.Branch("column_mip_prior_converged", &r.columnMipPriorConverged);
    tree.Branch("column_mip_prior_stalled", &r.columnMipPriorStalled);
    tree.Branch("column_mip_prior_iterations", &r.columnMipPriorIterations);
    tree.Branch("column_mip_prior_eligible_columns", &r.columnMipPriorEligibleColumns);
    tree.Branch("column_mip_prior_floor_energy", &r.columnMipPriorFloorEnergy);
    tree.Branch("column_mip_prior_mpv_energy", &r.columnMipPriorMpvEnergy);
    tree.Branch("column_mip_prior_strength", &r.columnMipPriorStrength);
    tree.Branch("column_mip_prior_charge_per_energy", &r.columnMipPriorChargePerEnergy);
    tree.Branch("column_mip_prior_min_path_fraction", &r.columnMipPriorMinPathFraction);
    tree.Branch("column_mip_prior_penalty", &r.columnMipPriorPenalty);
    tree.Branch("column_mip_prior_objective", &r.columnMipPriorObjective);
    tree.Branch("column_mip_prior_max_relative_change", &r.columnMipPriorMaxRelativeChange);
    tree.Branch("column_unregularized_nll", &r.columnUnregularizedNll);
    tree.Branch("input_tree", &r.inputTree);
}

inline void SetGlobalFitAddresses(TTree& tree, GlobalFitResult& r) {
#define GF_ADDRESS(name, member) do { if (tree.GetBranch(name)) tree.SetBranchAddress(name, &r.member); } while(false)
    GF_ADDRESS("event", event); GF_ADDRESS("status", status);
    GF_ADDRESS("observations", observations); GF_ADDRESS("minimum_charge", minimumCharge);
    GF_ADDRESS("fit_range_quantile", fitRangeQuantile); GF_ADDRESS("fit_range_padding", fitRangePadding);
    GF_ADDRESS("fit_low_x", fitLowX); GF_ADDRESS("fit_low_y", fitLowY); GF_ADDRESS("fit_low_z", fitLowZ);
    GF_ADDRESS("fit_high_x", fitHighX); GF_ADDRESS("fit_high_y", fitHighY); GF_ADDRESS("fit_high_z", fitHighZ);
    GF_ADDRESS("observations_before_clustering", observationsBeforeClustering);
    GF_ADDRESS("observations_after_dbscan", observationsAfterDbscan);
    GF_ADDRESS("dbscan_enabled", dbscanEnabled); GF_ADDRESS("dbscan_epsilon", dbscanEpsilon);
    GF_ADDRESS("dbscan_min_points", dbscanMinPoints);
    GF_ADDRESS("corridor_enabled", corridorEnabled);
    GF_ADDRESS("corridor_half_width", corridorHalfWidth);
    GF_ADDRESS("minimum_map_fraction", minimumMapFraction);
    GF_ADDRESS("fit_range_enabled", fitRangeEnabled);
    GF_ADDRESS("nll", nll);
    GF_ADDRESS("edm", edm); GF_ADDRESS("function_calls", functionCalls);
    GF_ADDRESS("maximum_function_calls", maximumFunctionCalls); GF_ADDRESS("tolerance", tolerance);
    GF_ADDRESS("chi2", chi2); GF_ADDRESS("ndof", ndof);
    GF_ADDRESS("seed_x", seedX); GF_ADDRESS("seed_y", seedY); GF_ADDRESS("seed_z", seedZ);
    GF_ADDRESS("seed_dx", seedDx); GF_ADDRESS("seed_dy", seedDy); GF_ADDRESS("seed_dz", seedDz);
    GF_ADDRESS("seed_observations", seedObservations);
    GF_ADDRESS("seed_median_factor", seedMedianFactor);
    GF_ADDRESS("seed_median_xy", seedMedianXY); GF_ADDRESS("seed_median_xz", seedMedianXZ);
    GF_ADDRESS("seed_median_yz", seedMedianYZ);
    GF_ADDRESS("seed_nll", seedNll); GF_ADDRESS("seed_chi2", seedChi2); GF_ADDRESS("seed_ndof", seedNdof);
    GF_ADDRESS("fit_x", fitX); GF_ADDRESS("fit_y", fitY); GF_ADDRESS("fit_z", fitZ);
    GF_ADDRESS("fit_dx", fitDx); GF_ADDRESS("fit_dy", fitDy); GF_ADDRESS("fit_dz", fitDz);
    GF_ADDRESS("error_a", errorA); GF_ADDRESS("error_b", errorB);
    GF_ADDRESS("error_theta", errorTheta); GF_ADDRESS("error_phi", errorPhi);
    GF_ADDRESS("column_fit_available", columnFitAvailable); GF_ADDRESS("column_fit_converged", columnFitConverged);
    GF_ADDRESS("column_fit_iterations", columnFitIterations); GF_ADDRESS("column_fit_active_columns", columnFitActiveColumns);
    GF_ADDRESS("column_fit_max_relative_change", columnFitMaxRelativeChange);
    GF_ADDRESS("column_nll", columnNll); GF_ADDRESS("column_chi2", columnChi2); GF_ADDRESS("column_ndof", columnNdof);
    GF_ADDRESS("column_mip_prior_enabled", columnMipPriorEnabled); GF_ADDRESS("column_mip_prior_converged", columnMipPriorConverged); GF_ADDRESS("column_mip_prior_stalled", columnMipPriorStalled);
    GF_ADDRESS("column_mip_prior_iterations", columnMipPriorIterations); GF_ADDRESS("column_mip_prior_eligible_columns", columnMipPriorEligibleColumns);
    GF_ADDRESS("column_mip_prior_floor_energy", columnMipPriorFloorEnergy); GF_ADDRESS("column_mip_prior_mpv_energy", columnMipPriorMpvEnergy);
    GF_ADDRESS("column_mip_prior_strength", columnMipPriorStrength); GF_ADDRESS("column_mip_prior_charge_per_energy", columnMipPriorChargePerEnergy);
    GF_ADDRESS("column_mip_prior_min_path_fraction", columnMipPriorMinPathFraction); GF_ADDRESS("column_mip_prior_penalty", columnMipPriorPenalty);
    GF_ADDRESS("column_mip_prior_objective", columnMipPriorObjective); GF_ADDRESS("column_mip_prior_max_relative_change", columnMipPriorMaxRelativeChange);
    GF_ADDRESS("column_unregularized_nll", columnUnregularizedNll);
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
