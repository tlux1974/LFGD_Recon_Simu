#include <TND280Input.hxx>
#include <TND280Event.hxx>
#include <TDataVector.hxx>
#include <TG4VHit.hxx>
#include <TG4HitSegment.hxx>

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>

void CompareOpticalHitCounting(const char* filename, int eventNumber = 0) {
    ND::TND280Input input(filename);
    if (!input.IsOpen()) throw std::runtime_error("Cannot open input file");
    ND::TND280Event* event = input.ReadEvent(eventNumber);
    if (!event) throw std::runtime_error("Cannot read requested event");
    ND::THandle<ND::TDataVector> collections =
        event->Get<ND::TDataVector>("truth/g4Hits");
    if (!collections) throw std::runtime_error("No truth/g4Hits");

    long long segments = 0;
    long long contributorReferences = 0;
    double depositedEV = 0.0;
    std::unordered_set<int> contributorIds;
    for (auto item = collections->begin(); item != collections->end(); ++item) {
        ND::THandle<ND::TG4HitContainer> hits = (*item)->Get<ND::TG4HitContainer>(".");
        if (!hits || std::string(hits->GetName()).find("opticalhomo") != 0) continue;
        for (const auto* hit : *hits) {
            const auto* segment = dynamic_cast<const ND::TG4HitSegment*>(hit);
            if (!segment) continue;
            ++segments;
            // oaEvent persists energies in MeV; one eV is 1e-6 MeV.
            depositedEV += segment->GetEnergyDeposit() / 1.0e-6;
            contributorReferences += segment->GetContributors().size();
            contributorIds.insert(segment->GetContributors().begin(),
                                  segment->GetContributors().end());
        }
    }
    std::cout << "event " << eventNumber
              << " optical hit segments " << segments
              << " contributor references " << contributorReferences
              << " unique optical tracks " << contributorIds.size()
              << " deposited energy/eV " << depositedEV << '\n';
}
