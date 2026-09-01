#include <TDirectory.h>
#include <TFile.h>
#include <TH3.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

void AddLightMapEfficiencyTree(const char* input, int launchedPhotons = 100000) {
    if (launchedPhotons <= 0) {
        std::cerr << "launchedPhotons must be positive\n";
        gSystem->Exit(2);
        return;
    }
    TFile file(input, "UPDATE");
    if (file.IsZombie()) {
        std::cerr << "Cannot open light map for update: " << input << '\n';
        gSystem->Exit(3);
        return;
    }
    if (file.Get("position_efficiency")) {
        std::cerr << input << " already contains position_efficiency; unchanged\n";
        gSystem->Exit(4);
        return;
    }
    auto* starts = file.Get<TH3>("Starting Vertex Distribution");
    auto* fractions = file.Get<TDirectory>("lightFractions");
    auto* centre = file.Get<TVector3>("cubeCentre");
    if (!starts || !fractions || !centre) {
        std::cerr << "Missing legacy light-map objects\n";
        gSystem->Exit(5);
        return;
    }

    int bin = -1;
    int launched_photons = launchedPhotons;
    long long fibre_hits = 0;
    double x_mm = 0.0, y_mm = 0.0, z_mm = 0.0;
    double fibre_hit_efficiency = 0.0;
    TTree tree("position_efficiency", "Per-position fibre-hit efficiency");
    tree.Branch("bin", &bin);
    tree.Branch("x_mm", &x_mm);
    tree.Branch("y_mm", &y_mm);
    tree.Branch("z_mm", &z_mm);
    tree.Branch("launched_photons", &launched_photons);
    tree.Branch("fibre_hits", &fibre_hits);
    tree.Branch("fibre_hit_efficiency", &fibre_hit_efficiency);

    for (int ix = 1; ix <= starts->GetNbinsX(); ++ix)
        for (int iy = 1; iy <= starts->GetNbinsY(); ++iy)
            for (int iz = 1; iz <= starts->GetNbinsZ(); ++iz) {
                bin = starts->GetBin(ix, iy, iz);
                if (starts->GetBinContent(bin) <= 0) continue;
                auto* values = fractions->Get<std::vector<double>>(
                    ("Bin_" + std::to_string(bin) + "_light_fractions").c_str());
                if (!values) {
                    std::cerr << "Missing light fractions for populated bin "
                              << bin << '\n';
                    gSystem->Exit(6);
                    return;
                }
                fibre_hit_efficiency = 0.0;
                for (double value : *values) fibre_hit_efficiency += value;
                fibre_hits = std::llround(
                    fibre_hit_efficiency * launched_photons);
                x_mm = centre->X() + starts->GetXaxis()->GetBinCenter(ix);
                y_mm = centre->Y() + starts->GetYaxis()->GetBinCenter(iy);
                z_mm = centre->Z() + starts->GetZaxis()->GetBinCenter(iz);
                tree.Fill();
            }

    file.cd();
    tree.Write();
    file.Close();
    std::cout << "Added position_efficiency with " << tree.GetEntries()
              << " positions to " << input << '\n';
}
