#include <TFile.h>
#include <TH3.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>
#include <vector>

#include <fstream>
#include <iomanip>
#include <iostream>

void ExportLightMapEfficiency(const char* input, const char* output) {
    TFile file(input, "READ");
    if (file.IsZombie()) {
        std::cerr << "Cannot open light map: " << input << '\n';
        gSystem->Exit(2);
        return;
    }
    auto* tree = file.Get<TTree>("position_efficiency");

    std::ofstream csv(output);
    if (!csv) {
        std::cerr << "Cannot create CSV: " << output << '\n';
        gSystem->Exit(4);
        return;
    }
    csv << "entry,bin,x_mm,y_mm,z_mm,launched_photons,fibre_hits,fibre_hit_efficiency\n";
    csv << std::setprecision(10);

    if (!tree) {
        auto* starts = file.Get<TH3>("Starting Vertex Distribution");
        auto* fractions = file.Get<TDirectory>("lightFractions");
        auto* centre = file.Get<TVector3>("cubeCentre");
        if (!starts || !fractions || !centre) {
            std::cerr << "Missing legacy light-map objects\n";
            gSystem->Exit(3);
            return;
        }
        Long64_t entry = 0;
        for (int ix = 1; ix <= starts->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= starts->GetNbinsY(); ++iy)
                for (int iz = 1; iz <= starts->GetNbinsZ(); ++iz) {
                    const int bin = starts->GetBin(ix, iy, iz);
                    if (starts->GetBinContent(bin) <= 0) continue;
                    auto* values = fractions->Get<std::vector<double>>(
                        ("Bin_" + std::to_string(bin) + "_light_fractions").c_str());
                    if (!values) continue;
                    double efficiency = 0.0;
                    for (double value : *values) efficiency += value;
                    const long long hits = std::llround(efficiency * 100000.0);
                    csv << entry++ << ',' << bin << ','
                        << centre->X() + starts->GetXaxis()->GetBinCenter(ix) << ','
                        << centre->Y() + starts->GetYaxis()->GetBinCenter(iy) << ','
                        << centre->Z() + starts->GetZaxis()->GetBinCenter(iz) << ','
                        << 100000 << ',' << hits << ',' << efficiency << '\n';
                }
        std::cout << "Wrote " << entry << " legacy-map positions to " << output << '\n';
        return;
    }

    int bin = -1;
    int launched = 0;
    long long hits = 0;
    double x = 0.0, y = 0.0, z = 0.0, efficiency = 0.0;
    tree->SetBranchAddress("bin", &bin);
    tree->SetBranchAddress("x_mm", &x);
    tree->SetBranchAddress("y_mm", &y);
    tree->SetBranchAddress("z_mm", &z);
    tree->SetBranchAddress("launched_photons", &launched);
    tree->SetBranchAddress("fibre_hits", &hits);
    tree->SetBranchAddress("fibre_hit_efficiency", &efficiency);

    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);
        csv << entry << ',' << bin << ',' << x << ',' << y << ',' << z << ','
            << launched << ',' << hits << ',' << efficiency << '\n';
    }
    std::cout << "Wrote " << tree->GetEntries() << " positions to " << output << '\n';
}
