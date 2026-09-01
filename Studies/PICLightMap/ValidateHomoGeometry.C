#include <TFile.h>
#include <TGeoManager.h>
#include <TKey.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

void ValidateHomoGeometry(const char* filename) {
    TFile file(filename, "READ");
    if (file.IsZombie()) {
        std::cerr << "VALIDATION FAILED: cannot open " << filename << '\n';
        gSystem->Exit(10);
        return;
    }

    TGeoManager* geometry = nullptr;
    TIter next(file.GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        if (TString(key->GetClassName()) == "TGeoManager") {
            geometry = dynamic_cast<TGeoManager*>(key->ReadObj());
            break;
        }
    }
    if (!geometry) {
        std::cerr << "VALIDATION FAILED: output contains no TGeoManager\n";
        gSystem->Exit(11);
        return;
    }
    if (!geometry->FindVolumeFast("HOMO")) {
        std::cerr << "VALIDATION FAILED: HOMO volume is absent; the wrong detector was simulated\n";
        gSystem->Exit(12);
        return;
    }
    if (geometry->FindVolumeFast("HFG")) {
        std::cerr << "VALIDATION FAILED: HFG is present together with HOMO\n";
        gSystem->Exit(13);
        return;
    }
    // A streamer is written only when this class was actually serialized in
    // the file. This check works without loading the oaEvent dictionaries and
    // prevents geometry-only files with no optical fibre hits being accepted.
    if (!file.GetStreamerInfoList()->FindObject("ND::TG4HitContainer") ||
        !file.GetStreamerInfoList()->FindObject("ND::TG4HitSegment")) {
        std::cerr << "VALIDATION FAILED: no serialized optical hit collection was found\n";
        gSystem->Exit(14);
        return;
    }
    std::cout << "VALIDATION PASSED: HOMO geometry and optical fibre hits are present\n";
}
