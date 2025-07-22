#include <algorithm>
#include <cmath>

#include "Combiners.h"

auto ExhaustiveCompetitive::Combine(const std::vector<Parton*>& partons) -> std::vector<Hadron*> {
  std::vector<Hadron*> hadrons;

  for (auto* a : partons) {
    if (a->IsUsed()) continue;

    Parton* closestOpp = nullptr;
    double closestOppDist = std::numeric_limits<double>::max();
    Parton* closestSame1 = nullptr;
    Parton* closestSame2 = nullptr;
    double bestTripletDist = std::numeric_limits<double>::max();

    for (auto* b : partons) {
      if (b == a || b->IsUsed()) continue;
      int sum = std::round(a->GetBaryonNumber() + b->GetBaryonNumber());
      double d = a->DistanceTo(*b);

      if (sum == 0 && d < closestOppDist) {
        closestOpp = b;
        closestOppDist = d;
      } else if (sum != 0 && a->GetBaryonNumber() == b->GetBaryonNumber()) {
        for (auto* c : partons) {
          if (c == a || c == b || c->IsUsed()) continue;
          if (c->GetBaryonNumber() != a->GetBaryonNumber()) continue;

          int bsum = std::round(a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber());
          if (std::abs(bsum) != 1) continue;

          double td = a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c);
          if (td < bestTripletDist) {
            bestTripletDist = td;
            closestSame1 = b;
            closestSame2 = c;
          }
        }
      }
    }

    if (closestOpp && (!closestSame1 || closestOppDist < (bestTripletDist / m_r))) {
      double x = (a->X() + closestOpp->X()) / 2;
      double y = (a->Y() + closestOpp->Y()) / 2;
      double z = (a->Z() + closestOpp->Z()) / 2;
      double px = a->Px() + closestOpp->Px();
      double py = a->Py() + closestOpp->Py();
      double pz = a->Pz() + closestOpp->Pz();

      // Kinematic calculation for meson
      double m1 = a->GetMassFromPDG();
      double m2 = closestOpp->GetMassFromPDG();
      double E1 = std::sqrt(a->Px() * a->Px() + a->Py() * a->Py() + a->Pz() * a->Pz() + m1 * m1);
      double E2 = std::sqrt(closestOpp->Px() * closestOpp->Px() + closestOpp->Py() * closestOpp->Py() +
                            closestOpp->Pz() * closestOpp->Pz() + m2 * m2);
      double E_sum = E1 + E2;
      double invMass = std::sqrt(std::max(0.0, E_sum * E_sum - (px * px + py * py + pz * pz)));

      auto h = new Hadron(x, y, z, px, py, pz, 0, closestOppDist);
      h->SetMass(invMass);
      h->AddConstituentID(a->UniqueID());
      h->AddConstituentID(closestOpp->UniqueID());
      hadrons.push_back(h);
      a->MarkUsed();
      closestOpp->MarkUsed();
    } else if (closestSame1 && closestSame2) {
      double x = (a->X() + closestSame1->X() + closestSame2->X()) / 3;
      double y = (a->Y() + closestSame1->Y() + closestSame2->Y()) / 3;
      double z = (a->Z() + closestSame1->Z() + closestSame2->Z()) / 3;
      double px = a->Px() + closestSame1->Px() + closestSame2->Px();
      double py = a->Py() + closestSame1->Py() + closestSame2->Py();
      double pz = a->Pz() + closestSame1->Pz() + closestSame2->Pz();
      int baryonNumber =
          std::round(a->GetBaryonNumber() + closestSame1->GetBaryonNumber() + closestSame2->GetBaryonNumber());

      // formation distance for baryon
      double formationDist = bestTripletDist;

      // Kinematic calculation for baryon
      double m1 = a->GetMassFromPDG();
      double m2 = closestSame1->GetMassFromPDG();
      double m3 = closestSame2->GetMassFromPDG();
      double E1 = std::sqrt(a->Px() * a->Px() + a->Py() * a->Py() + a->Pz() * a->Pz() + m1 * m1);
      double E2 = std::sqrt(closestSame1->Px() * closestSame1->Px() + closestSame1->Py() * closestSame1->Py() +
                            closestSame1->Pz() * closestSame1->Pz() + m2 * m2);
      double E3 = std::sqrt(closestSame2->Px() * closestSame2->Px() + closestSame2->Py() * closestSame2->Py() +
                            closestSame2->Pz() * closestSame2->Pz() + m3 * m3);
      double E_sum = E1 + E2 + E3;
      double invMass = std::sqrt(std::max(0.0, E_sum * E_sum - (px * px + py * py + pz * pz)));

      auto h = new Hadron(x, y, z, px, py, pz, baryonNumber, formationDist);
      h->SetMass(invMass);
      h->AddConstituentID(a->UniqueID());
      h->AddConstituentID(closestSame1->UniqueID());
      h->AddConstituentID(closestSame2->UniqueID());
      hadrons.push_back(h);
      a->MarkUsed();
      closestSame1->MarkUsed();
      closestSame2->MarkUsed();
    }
  }

  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}
