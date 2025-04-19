#include "Combiners.h"
#include "core/PartonKDTree.h"
#include <algorithm>
#include <cmath>

struct Candidate {
  double distance;
  bool isBaryon;
  Parton* a;
  Parton* b;
  Parton* c;

  Candidate(double d, Parton* x, Parton* y)
      : distance(d), isBaryon(false), a(x), b(y), c(nullptr) {}

  Candidate(double d, Parton* x, Parton* y, Parton* z)
      : distance(d), isBaryon(true), a(x), b(y), c(z) {}

  bool operator<(const Candidate& other) const { return distance < other.distance; }
};

std::vector<Hadron*> KDTreeGlobal::Combine(const std::vector<Parton*>& partons) {
  std::vector<Hadron*> hadrons;
  PartonKDTree tree(partons);
  std::vector<Candidate> candidates;

  // meson candidates
  for (auto* a : partons) {
    auto neighbors = tree.FindNeighbors(a, 20);
    for (auto& [b, dist] : neighbors) {
      if (a == b) continue;
      if (std::round(a->GetBaryonNumber() + b->GetBaryonNumber()) != 0) continue;
      candidates.emplace_back(dist, a, b);
    }
  }

  // baryon candidates
  for (auto* a : partons) {
    auto neighbors = tree.FindNeighbors(a, 10);
    for (size_t i = 0; i < neighbors.size(); ++i) {
      Parton* b = neighbors[i].first;
      if (a == b) continue;
      for (size_t j = i + 1; j < neighbors.size(); ++j) {
        Parton* c = neighbors[j].first;
        if (c == a || c == b ) continue;
        double baryonNumber = a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber();
        if (std::round(baryonNumber) != 1 && std::round(baryonNumber) != -1) continue;
        double dist = a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c);
        dist *= m_r;
        candidates.emplace_back(dist, a, b, c);
      }
    }
  }

  std::sort(candidates.begin(), candidates.end());

  for (const auto& candi : candidates) {
    if (candi.isBaryon) {
      if (candi.a->IsUsed() || candi.b->IsUsed() || candi.c->IsUsed()) continue;
      double px = candi.a->Px() + candi.b->Px() + candi.c->Px();
      double py = candi.a->Py() + candi.b->Py() + candi.c->Py();
      double pz = candi.a->Pz() + candi.b->Pz() + candi.c->Pz();
      double x = (candi.a->X() + candi.b->X() + candi.c->X()) / 3;
      double y = (candi.a->Y() + candi.b->Y() + candi.c->Y()) / 3;
      double z = (candi.a->Z() + candi.b->Z() + candi.c->Z()) / 3;
      double formation = candi.a->DistanceTo(*candi.b) + candi.a->DistanceTo(*candi.c) + candi.b->DistanceTo(*candi.c);

      auto* h = new Hadron(x,y,z,px,py,pz,
                          std::round(candi.a->GetBaryonNumber() +
                          candi.b->GetBaryonNumber() +
                          candi.c->GetBaryonNumber()),
                          formation);
      
      h->SetPosition(x, y, z);
      h->AddConstituentID(candi.a->UniqueID());
      h->AddConstituentID(candi.b->UniqueID());
      h->AddConstituentID(candi.c->UniqueID());
      hadrons.push_back(h);
      candi.a->MarkUsed(); candi.b->MarkUsed(); candi.c->MarkUsed();
    } else {
      if (candi.a->IsUsed() || candi.b->IsUsed()) continue;
      double px = candi.a->Px() + candi.b->Px();
      double py = candi.a->Py() + candi.b->Py();
      double pz = candi.a->Pz() + candi.b->Pz();
      double x = (candi.a->X() + candi.b->X()) / 2;
      double y = (candi.a->Y() + candi.b->Y()) / 2;
      double z = (candi.a->Z() + candi.b->Z()) / 2;
      double formation = candi.a->DistanceTo(*candi.b);

      auto* h = new Hadron(x,y,z,px,py,pz,
                          std::round(candi.a->GetBaryonNumber() +
                          candi.b->GetBaryonNumber()),
                          formation);
      h->SetPosition(x, y, z);
      h->AddConstituentID(candi.a->UniqueID());
      h->AddConstituentID(candi.b->UniqueID());
      hadrons.push_back(h);
      candi.a->MarkUsed(); candi.b->MarkUsed();
    }
  }

  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}
