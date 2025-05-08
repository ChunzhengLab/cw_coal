#include <algorithm>
#include <cmath>
#include <limits>

#include "Combiners.h"
#include "core/Particle.h"
#include "core/PartonKDTree.h"

std::vector<Hadron*> KDTreeDualGreedy::Combine(const std::vector<Parton*>& partons) {
  std::vector<Hadron*> hadrons;
  PartonKDTree tree(partons);

  for (auto* a : partons) {
    if (a->IsUsed()) continue;

    // Collect all neighbors sorted by distance
    auto neighbors = tree.kNearestSearch(a, partons.size());
    // Find best meson partner
    Parton* bestOpp = nullptr;
    double bestOppDist = std::numeric_limits<double>::max();
    // Find best baryon pair
    Parton* bestSame1 = nullptr;
    Parton* bestSame2 = nullptr;
    double bestTripletDist = std::numeric_limits<double>::max();

    for (auto& [b, d_ab] : neighbors) {
      if (b == a || b->IsUsed()) continue;
      int sum2 = std::round(a->GetBaryonNumber() + b->GetBaryonNumber());
      if (sum2 == 0) {
        if (d_ab < bestOppDist) {
          bestOppDist = d_ab;
          bestOpp = b;
        }
      } else if (b->GetBaryonNumber() == a->GetBaryonNumber()) {
        for (auto& [c, d_ac] : neighbors) {
          if (c == a || c == b || c->IsUsed() || c->GetBaryonNumber() != a->GetBaryonNumber()) continue;
          int bsum = std::round(a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber());
          if (std::abs(bsum) != 1) continue;
          double d_bc = b->DistanceTo(*c);
          double td = d_ab + d_ac + d_bc;
          if (td < bestTripletDist) {
            bestTripletDist = td;
            bestSame1 = b;
            bestSame2 = c;
          }
        }
      }
    }
    // Assign meson and baryon results
    Parton* mesonNeighbor = bestOpp;
    double mesonDist = bestOppDist;
    std::vector<Parton*> sameNeighbors;
    double baryonDist = std::numeric_limits<double>::max();
    if (bestSame1 && bestSame2) {
      sameNeighbors.push_back(bestSame1);
      sameNeighbors.push_back(bestSame2);
      baryonDist = bestTripletDist / m_r;
    }

    if (mesonDist < baryonDist && mesonNeighbor && !mesonNeighbor->IsUsed()) {
      // Form meson
      double px = a->Px() + mesonNeighbor->Px();
      double py = a->Py() + mesonNeighbor->Py();
      double pz = a->Pz() + mesonNeighbor->Pz();
      double x = (a->X() + mesonNeighbor->X()) / 2;
      double y = (a->Y() + mesonNeighbor->Y()) / 2;
      double z = (a->Z() + mesonNeighbor->Z()) / 2;
      // Kinematic calculation for meson
      double m1 = a->GetMassFromPDG();
      double m2 = mesonNeighbor->GetMassFromPDG();
      double E1 = std::sqrt(a->Px() * a->Px() + a->Py() * a->Py() + a->Pz() * a->Pz() + m1 * m1);
      double E2 = std::sqrt(mesonNeighbor->Px() * mesonNeighbor->Px() + mesonNeighbor->Py() * mesonNeighbor->Py() +
                            mesonNeighbor->Pz() * mesonNeighbor->Pz() + m2 * m2);
      double E_sum = E1 + E2;
      double invMass = std::sqrt(std::max(0.0, E_sum * E_sum - (px * px + py * py + pz * pz)));

      auto* h = new Hadron(x, y, z, px, py, pz, 0, mesonDist);
      h->SetMass(invMass);
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
      // Kinematic calculation for baryon
      double m1 = a->GetMassFromPDG();
      double m2 = b->GetMassFromPDG();
      double m3 = c->GetMassFromPDG();
      double E1 = std::sqrt(a->Px() * a->Px() + a->Py() * a->Py() + a->Pz() * a->Pz() + m1 * m1);
      double E2 = std::sqrt(b->Px() * b->Px() + b->Py() * b->Py() + b->Pz() * b->Pz() + m2 * m2);
      double E3 = std::sqrt(c->Px() * c->Px() + c->Py() * c->Py() + c->Pz() * c->Pz() + m3 * m3);
      double E_sum = E1 + E2 + E3;
      double invMass = std::sqrt(std::max(0.0, E_sum * E_sum - (px * px + py * py + pz * pz)));

      // Use raw sum of three edges as fill distance
      double rawBaryonDist = bestTripletDist;
      auto* h = new Hadron(x, y, z, px, py, pz, baryonNumber, rawBaryonDist);
      h->SetMass(invMass);
      h->AddConstituentID(a->UniqueID());
      h->AddConstituentID(b->UniqueID());
      h->AddConstituentID(c->UniqueID());
      hadrons.push_back(h);
      a->MarkUsed();
      b->MarkUsed();
      c->MarkUsed();
    }
  }

  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}
