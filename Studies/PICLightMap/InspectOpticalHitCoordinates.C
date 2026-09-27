#include <TND280Input.hxx>
#include <TND280Event.hxx>
#include <TDataVector.hxx>
#include <TG4VHit.hxx>
#include <TG4HitSegment.hxx>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
double residue(double value, double origin) {
    double result = std::fmod(value-origin+5.0,10.0);
    if (result < 0.0) result += 10.0;
    return result-5.0;
}
}

void InspectOpticalHitCoordinates(const char* filename, int eventNumber=0,
                                  int maximumSegments=30) {
    ND::TND280Input input(filename);
    if (!input.IsOpen()) throw std::runtime_error("Cannot open input file");
    ND::TND280Event* event=input.ReadEvent(eventNumber);
    if (!event) throw std::runtime_error("Cannot read requested event");
    auto collections=event->Get<ND::TDataVector>("truth/g4Hits");
    if (!collections) throw std::runtime_error("No truth/g4Hits");

    int shown=0;
    std::cout << "Expected lattice origins relative to global coordinates:\n"
              << "  axis 0: x=+2.5 mod 10, z=907.5 mod 10\n"
              << "  axis 1: y=27.5 mod 10, z=912.5 mod 10\n"
              << "  axis 2: x=-2.5 mod 10, y=32.5 mod 10\n";
    std::cout << "start_xyz -> stop_xyz; contributor_count; minimum transverse distances [mm]\n";
    for (auto item=collections->begin(); item!=collections->end(); ++item) {
        auto hits=(*item)->Get<ND::TG4HitContainer>(".");
        if (!hits || std::string(hits->GetName()).find("opticalhomo")!=0) continue;
        for (const auto* hit:*hits) {
            const auto* s=dynamic_cast<const ND::TG4HitSegment*>(hit);
            if (!s) continue;
            const double x=s->GetStopX(), y=s->GetStopY(), z=s->GetStopZ();
            const double d0=std::hypot(residue(x,2.5),residue(z,907.5));
            const double d1=std::hypot(residue(y,27.5),residue(z,912.5));
            const double d2=std::hypot(residue(x,-2.5),residue(y,32.5));
            std::cout << std::fixed << std::setprecision(4)
                      << s->GetStartX() << ',' << s->GetStartY() << ',' << s->GetStartZ()
                      << " -> " << x << ',' << y << ',' << z << "; "
                      << s->GetContributors().size() << "; "
                      << d0 << ',' << d1 << ',' << d2 << '\n';
            if (++shown>=maximumSegments) return;
        }
    }
    std::cout << "Only " << shown << " optical HOMO segments found\n";
}
