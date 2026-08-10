#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

namespace {

struct FibreObservation {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double charge = 0.0;
    int projection = -1;
};

struct TrackSeed {
    std::array<double, 3> start{{0.0, 0.0, 0.0}};
    std::array<double, 3> direction{{0.0, 0.0, 1.0}};
};

// Projection 0 is Y-directed, 1 is X-directed and 2 is Z-directed.
// A fibre constrains only the two coordinates transverse to its direction.
bool MeasuresCoordinate(int projection, int coordinate) {
    static const bool measured[3][3] = {
        {true,  false, true },
        {false, true,  true },
        {true,  true,  false}
    };
    return projection >= 0 && projection < 3 && measured[projection][coordinate];
}

std::vector<FibreObservation> ReadEvent(TTree& tree, int requestedEvent,
                                        double minimumCharge) {
    int event = 0;
    int projection = 0;
    double x = 0.0, y = 0.0, z = 0.0, charge = 0.0;

    tree.SetBranchAddress("event", &event);
    tree.SetBranchAddress("x", &x);
    tree.SetBranchAddress("y", &y);
    tree.SetBranchAddress("z", &z);
    tree.SetBranchAddress("charge", &charge);
    tree.SetBranchAddress("projection", &projection);

    std::vector<FibreObservation> result;
    for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        tree.GetEntry(entry);
        if (event < requestedEvent) continue;
        if (event > requestedEvent) break; // Flat trees are event ordered.
        if (charge < minimumCharge || !std::isfinite(charge)) continue;
        result.push_back({x, y, z, charge, projection});
    }
    return result;
}

TrackSeed MakeSeed(const std::vector<FibreObservation>& observations,
                   const std::array<double, 3>& detectorCentre) {
    TrackSeed seed;
    seed.start = detectorCentre;
    std::array<double, 3> farthest{{0.0, 0.0, 0.0}};

    for (const auto& hit : observations) {
        const std::array<double, 3> coordinate{{hit.x, hit.y, hit.z}};
        for (int axis = 0; axis < 3; ++axis) {
            if (!MeasuresCoordinate(hit.projection, axis)) continue;
            const double displacement = coordinate[axis] - seed.start[axis];
            if (std::abs(displacement) > std::abs(farthest[axis])) {
                farthest[axis] = displacement;
            }
        }
    }

    const double norm = std::sqrt(farthest[0] * farthest[0]
                                + farthest[1] * farthest[1]
                                + farthest[2] * farthest[2]);
    if (!(norm > 0.0)) throw std::runtime_error("Cannot determine a direction seed");
    for (int axis = 0; axis < 3; ++axis) seed.direction[axis] = farthest[axis] / norm;
    return seed;
}

void PlotDistributions(const std::vector<FibreObservation>& observations,
                       int event, const std::string& prefix) {
    std::array<double, 3> low{{ std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::max() }};
    std::array<double, 3> high{{ -std::numeric_limits<double>::max(),
                                -std::numeric_limits<double>::max(),
                                -std::numeric_limits<double>::max() }};
    for (const auto& hit : observations) {
        const std::array<double, 3> coordinate{{hit.x, hit.y, hit.z}};
        for (int axis = 0; axis < 3; ++axis) {
            if (!MeasuresCoordinate(hit.projection, axis)) continue;
            low[axis] = std::min(low[axis], coordinate[axis]);
            high[axis] = std::max(high[axis], coordinate[axis]);
        }
    }

    const char* axisName[3] = {"x", "y", "z"};
    std::array<TH1D*, 3> histogram{};
    for (int axis = 0; axis < 3; ++axis) {
        // Fibre coordinates are on a 10 mm lattice. Pad by half a pitch.
        const int bins = std::max(1, static_cast<int>(std::lround(
            (high[axis] - low[axis]) / 10.0)) + 1);
        const std::string name = std::string("charge_") + axisName[axis];
        const std::string title = "Event " + std::to_string(event)
            + ";" + axisName[axis] + " [mm];summed fibre charge";
        histogram[axis] = new TH1D(name.c_str(), title.c_str(), bins,
                                   low[axis] - 5.0, high[axis] + 5.0);
    }

    for (const auto& hit : observations) {
        const std::array<double, 3> coordinate{{hit.x, hit.y, hit.z}};
        for (int axis = 0; axis < 3; ++axis) {
            if (MeasuresCoordinate(hit.projection, axis))
                histogram[axis]->Fill(coordinate[axis], hit.charge);
        }
    }

    TFile output((prefix + ".root").c_str(), "RECREATE");
    TCanvas canvas("coordinate_distributions", "Coordinate distributions", 1500, 500);
    canvas.Divide(3, 1);
    for (int axis = 0; axis < 3; ++axis) {
        canvas.cd(axis + 1);
        histogram[axis]->Draw("HIST");
        histogram[axis]->Write();
    }
    canvas.Write();
    canvas.SaveAs((prefix + ".png").c_str());
}

} // namespace

int main(int argc, char** argv) {
    auto usage=[](){std::cout<<R"(Usage:
  global_fit_seed FLAT.root [KEY=value ...]

Options:
  EVENT=0                 Event number
  TREE=homo_truth         Input fibre tree (also homo_raw or fiber_hits)
  OUTPUT_PREFIX=seed_event_N
                          ROOT/PNG output prefix
  MIN_CHARGE=10           Minimum fibre charge
)";};
    if(argc==1||(argc==2&&std::string(argv[1])=="--help")){usage();return 0;}
    if(argc<2){usage();return 2;}

    const std::string inputName = argv[1];
    std::map<std::string,std::string> options;
    for(int i=2;i<argc;++i){std::string argument=argv[i];auto equal=argument.find('=');if(equal==std::string::npos){std::cerr<<"Options must use KEY=value: "<<argument<<"\n";return 2;}
        std::string key=argument.substr(0,equal);if(key!="EVENT"&&key!="TREE"&&key!="OUTPUT_PREFIX"&&key!="MIN_CHARGE"){std::cerr<<"Unknown option: "<<key<<"\n";return 2;}options[key]=argument.substr(equal+1);}
    auto get=[&](const std::string& key,const std::string& fallback){auto it=options.find(key);return it==options.end()?fallback:it->second;};
    const int event = std::stoi(get("EVENT","0"));
    const std::string treeName = get("TREE","homo_truth");
    const std::string prefix = get("OUTPUT_PREFIX","seed_event_"+std::to_string(event));
    const double minimumCharge = std::stod(get("MIN_CHARGE","10"));
    if (minimumCharge < 0.0) {
        std::cerr << "ERROR: MIN_CHARGE must not be negative\n";
        return 2;
    }

    try {
        TFile input(inputName.c_str(), "READ");
        if (input.IsZombie()) throw std::runtime_error("Cannot open " + inputName);
        auto* tree = dynamic_cast<TTree*>(input.Get(treeName.c_str()));
        if (!tree) throw std::runtime_error("Cannot find tree " + treeName);
        for (const char* branch : {"event", "x", "y", "z", "charge", "projection"}) {
            if (!tree->GetBranch(branch))
                throw std::runtime_error(std::string("Missing branch ") + branch);
        }

        const auto observations = ReadEvent(*tree, event, minimumCharge);
        if (observations.empty())
            throw std::runtime_error("No observations pass the charge cut");
        const std::array<double, 3> detectorCentre{{
            0.5 * (tree->GetMinimum("x") + tree->GetMaximum("x")),
            0.5 * (tree->GetMinimum("y") + tree->GetMaximum("y")),
            0.5 * (tree->GetMinimum("z") + tree->GetMaximum("z"))}};
        const TrackSeed seed = MakeSeed(observations, detectorCentre);
        PlotDistributions(observations, event, prefix);

        std::cout << std::setprecision(9)
                  << "tree: " << treeName << "\n"
                  << "event: " << event << "\n"
                  << "minimum charge: " << minimumCharge << "\n"
                  << "fibre observations: " << observations.size() << "\n"
                  << "start seed [mm]: " << seed.start[0] << " "
                  << seed.start[1] << " " << seed.start[2] << " (free in fit)\n"
                  << "direction seed: " << seed.direction[0] << " "
                  << seed.direction[1] << " " << seed.direction[2] << "\n"
                  << "plots: " << prefix << ".png and " << prefix << ".root\n";
    }
    catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
