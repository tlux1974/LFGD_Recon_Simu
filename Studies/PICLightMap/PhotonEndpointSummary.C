#include <TND280Input.hxx>
#include <TND280Event.hxx>
#include <TG4Trajectory.hxx>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using EndpointKey = std::tuple<int, int, std::string>;

const char* ProcessTypeName(int type) {
    switch (type) {
        case 0: return "undefined";
        case 1: return "transportation";
        case 2: return "electromagnetic";
        case 3: return "optical";
        case 4: return "hadronic";
        case 5: return "photo-lepton-hadron";
        case 6: return "decay";
        case 7: return "general";
        case 8: return "parameterization";
        case 9: return "user-defined";
        default: return "unknown";
    }
}
}

void PhotonEndpointSummary(const char* filename, int eventNumber = 0,
                           int maximumRows = 30) {
    ND::TND280Input input(filename);
    if (!input.IsOpen()) throw std::runtime_error("Cannot open input file");
    if (eventNumber < 0 || eventNumber >= input.GetEventsInFile())
        throw std::runtime_error("Event number outside input range");

    ND::TND280Event* event = input.ReadEvent(eventNumber);
    if (!event) throw std::runtime_error("Cannot read requested event");
    ND::THandle<ND::TG4TrajectoryContainer> trajectories =
        event->Get<ND::TG4TrajectoryContainer>("truth/G4Trajectories");
    if (!trajectories) throw std::runtime_error("No truth/G4Trajectories");

    long long optical = 0;
    long long zeroPoints = 0;
    long long onePoint = 0;
    long long multiplePoints = 0;
    std::map<EndpointKey, long long> endpoints;

    for (const auto& item : *trajectories) {
        const ND::TG4Trajectory& trajectory = item.second;
        if (trajectory.GetParticleName() != "opticalphoton") continue;
        ++optical;
        const auto& points = trajectory.GetTrajectoryPoints();
        if (points.empty()) {
            ++zeroPoints;
            continue;
        }
        if (points.size() == 1) ++onePoint;
        else ++multiplePoints;
        const auto& endpoint = points.back();
        endpoints[{endpoint.GetProcessType(), endpoint.GetProcessSubtype(),
                   endpoint.GetVolumeName()}]++;
    }

    std::vector<std::pair<EndpointKey, long long>> sorted(endpoints.begin(),
                                                          endpoints.end());
    std::sort(sorted.begin(), sorted.end(), [](const auto& left,
                                                const auto& right) {
        return left.second > right.second;
    });

    std::cout << "file: " << filename << '\n'
              << "event: " << eventNumber << '\n'
              << "all trajectories: " << trajectories->size() << '\n'
              << "optical trajectories: " << optical << '\n'
              << "zero/one/multiple trajectory points: " << zeroPoints << '/'
              << onePoint << '/' << multiplePoints << "\n\n"
              << "count fraction process_type process_subtype terminal_volume\n";
    int row = 0;
    for (const auto& item : sorted) {
        if (maximumRows >= 0 && row++ >= maximumRows) break;
        const auto& [key, count] = item;
        const auto& [type, subtype, volume] = key;
        std::cout << count << ' ' << std::fixed << std::setprecision(6)
                  << (optical ? double(count) / optical : 0.0) << ' '
                  << ProcessTypeName(type) << '(' << type << ") " << subtype
                  << ' ' << (volume.empty() ? "<empty>" : volume) << '\n';
    }
}
