#include "Combiners.h"
#include "core/Particle.h"
#include "core/PartonKDTree.h"
#include <cmath>
#include <limits>
#include <unordered_set>
#include <array>

std::vector<Hadron*> KDTreeDualGreedy::Combine(const std::vector<Parton*>& partons) {
    std::vector<Hadron*> hadrons;
    PartonKDTree tree(partons);

    for (auto* a : partons) {
        if (a->IsUsed()) continue;

        // Try meson (opposite baryon number)
        auto mesonNeighbor = tree.FindNearestOpposite(a);
        double mesonDist = std::numeric_limits<double>::max();
        if (mesonNeighbor && !mesonNeighbor->IsUsed()) {
            mesonDist = a->DistanceTo(*mesonNeighbor);
        }

        // Try baryon (same baryon number, form triplet)
        auto sameNeighbors = tree.FindNearestSame(a, 2);
        double baryonDist = std::numeric_limits<double>::max();
        if (sameNeighbors.size() == 2 && !sameNeighbors[0]->IsUsed() && !sameNeighbors[1]->IsUsed()) {
            auto d1 = a->DistanceTo(*sameNeighbors[0]);
            auto d2 = a->DistanceTo(*sameNeighbors[1]);
            baryonDist = (d1 + d2) * m_r;
        }

        if (mesonDist < baryonDist && mesonNeighbor && !mesonNeighbor->IsUsed()) {
            // Form meson
            double px = a->Px() + mesonNeighbor->Px();
            double py = a->Py() + mesonNeighbor->Py();
            double pz = a->Pz() + mesonNeighbor->Pz();
            double x = (a->X() + mesonNeighbor->X()) / 2;
            double y = (a->Y() + mesonNeighbor->Y()) / 2;
            double z = (a->Z() + mesonNeighbor->Z()) / 2;
            Hadron* h = new Hadron(
                x, y, z,
                px, py, pz,
                0, mesonDist
            );
            h->AddConstituentID(a->UniqueID());
            h->AddConstituentID(mesonNeighbor->UniqueID());
            hadrons.push_back(h);
            a->MarkUsed();
            mesonNeighbor->MarkUsed();
        } else if (baryonDist < std::numeric_limits<double>::max()) {
            auto* b = sameNeighbors[0];
            auto* c = sameNeighbors[1];
            double px = a->Px() + b->Px() + c->Px();
            double py = a->Py() + b->Py() + c->Py();
            double pz = a->Pz() + b->Pz() + c->Pz();
            double x = (a->X() + b->X() + c->X()) / 3;
            double y = (a->Y() + b->Y() + c->Y()) / 3;
            double z = (a->Z() + b->Z() + c->Z()) / 3;
            int baryonNumber = std::round(a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber());
            Hadron* h = new Hadron(
                x, y, z,
                px, py, pz,
                baryonNumber,
                baryonDist
            );
            h->AddConstituentID(a->UniqueID());
            h->AddConstituentID(b->UniqueID());
            h->AddConstituentID(c->UniqueID());
            hadrons.push_back(h);
            a->MarkUsed();
            b->MarkUsed();
            c->MarkUsed();
        }
    }

    return hadrons;
}
