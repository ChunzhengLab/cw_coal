#include "Combiners.h"
#include <cmath>
#include <algorithm>


std::vector<Hadron*> BruteForceDualGreedy::Combine(const std::vector<Parton*>& partons) {
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
      Hadron* h = new Hadron(x, y, z, px, py, pz, 0, closestOppDist);
      h->AddConstituentID(a->UniqueID());
      h->AddConstituentID(closestOpp->UniqueID());
      hadrons.push_back(h);
      a->MarkUsed(); closestOpp->MarkUsed();
    } else if (closestSame1 && closestSame2) {
      double x = (a->X() + closestSame1->X() + closestSame2->X()) / 3;
      double y = (a->Y() + closestSame1->Y() + closestSame2->Y()) / 3;
      double z = (a->Z() + closestSame1->Z() + closestSame2->Z()) / 3;
      double px = a->Px() + closestSame1->Px() + closestSame2->Px();
      double py = a->Py() + closestSame1->Py() + closestSame2->Py();
      double pz = a->Pz() + closestSame1->Pz() + closestSame2->Pz();
      int baryonNumber = std::round(a->GetBaryonNumber() + closestSame1->GetBaryonNumber() + closestSame2->GetBaryonNumber());
      Hadron* h = new Hadron(x, y, z, px, py, pz, baryonNumber, bestTripletDist);
      h->AddConstituentID(a->UniqueID());
      h->AddConstituentID(closestSame1->UniqueID());
      h->AddConstituentID(closestSame2->UniqueID());
      hadrons.push_back(h);
      a->MarkUsed(); closestSame1->MarkUsed(); closestSame2->MarkUsed();
    }
  }

  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}
