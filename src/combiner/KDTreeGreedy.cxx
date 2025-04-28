#include "Particle.h"
#include <vector>
#include <nanoflann.hpp>
#include "Combiners.h"
#include "core/PartonKDTree.h"
#include <cmath>
#include <algorithm>
#include <random>

std::vector<Hadron*> KDTreeGreedy::Combine(const std::vector<Parton*>& partons) {
  // Baryon preference factor m_r: controls probability to reject meson
  std::mt19937 gen(std::random_device{}());
  std::bernoulli_distribution rejectMesonDist(m_r / (1.0 + m_r));

  PartonKDTree searcher(partons);
  std::vector<Hadron*> hadrons;

  std::unordered_set<const Parton*> used;

  // Try meson (2-parton) formation
  for (Parton* a : partons) {
    if (a->IsUsed()) continue;

    auto neighbors = searcher.kNearestSearch(a, 50);
    for (const auto& [b, dist] : neighbors) {
      if (b == a || b->IsUsed()) continue;
      if (std::round(a->GetBaryonNumber() + b->GetBaryonNumber()) != 0) continue;

      double px = a->Px() + b->Px();
      double py = a->Py() + b->Py();
      double pz = a->Pz() + b->Pz();
      double x = (a->X() + b->X()) / 2;
      double y = (a->Y() + b->Y()) / 2;
      double z = (a->Z() + b->Z()) / 2;
      // Kinematic calculation for meson
      double m1 = a->GetMassFromPDG();
      double m2 = b->GetMassFromPDG();
      double E1 = std::sqrt(a->Px()*a->Px() + a->Py()*a->Py() + a->Pz()*a->Pz() + m1*m1);
      double E2 = std::sqrt(b->Px()*b->Px() + b->Py()*b->Py() + b->Pz()*b->Pz() + m2*m2);
      double E_sum = E1 + E2;
      double invMass = std::sqrt(std::max(0.0, E_sum*E_sum - (px*px + py*py + pz*pz)));

      if (rejectMesonDist(gen)) {
        continue;  // probabilistically reject this meson
      }

      auto* h = new Hadron(x, y, z, px, py, pz, 0, dist);
      h->SetMass(invMass);
      h->AddConstituentID(a->UniqueID());
      h->AddConstituentID(b->UniqueID());
      hadrons.push_back(h);
      a->MarkUsed();
      b->MarkUsed();
      break;
    }
  }

  // Try baryon (3-parton) formation
  for (Parton* a : partons) {
    if (a->IsUsed()) continue;

    auto neighbors = searcher.kNearestSearch(a, 50);
    bool formed = false;

    for (size_t i = 0; i < neighbors.size(); ++i) {
      Parton* b = neighbors[i].first;
      if (b->IsUsed() || b == a) continue;
      for (size_t j = i + 1; j < neighbors.size(); ++j) {
        Parton* c = neighbors[j].first;
        if (c->IsUsed() || c == a || c == b) continue;

        double baryonSum = a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber();
        if (std::round(baryonSum) != 1 && std::round(baryonSum) != -1) continue;

        double dist = a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c);
        // No normalization: use raw sum of three edges
        // dist = dist / m_r;

        double px = a->Px() + b->Px() + c->Px();
        double py = a->Py() + b->Py() + c->Py();
        double pz = a->Pz() + b->Pz() + c->Pz();
        double x = (a->X() + b->X() + c->X()) / 3;
        double y = (a->Y() + b->Y() + c->Y()) / 3;
        double z = (a->Z() + b->Z() + c->Z()) / 3;

        // Kinematic calculation for baryon
        double m1 = a->GetMassFromPDG();
        double m2 = b->GetMassFromPDG();
        double m3 = c->GetMassFromPDG();
        double E1 = std::sqrt(a->Px()*a->Px() + a->Py()*a->Py() + a->Pz()*a->Pz() + m1*m1);
        double E2 = std::sqrt(b->Px()*b->Px() + b->Py()*b->Py() + b->Pz()*b->Pz() + m2*m2);
        double E3 = std::sqrt(c->Px()*c->Px() + c->Py()*c->Py() + c->Pz()*c->Pz() + m3*m3);
        double E_sum = E1 + E2 + E3;
        double invMass = std::sqrt(std::max(0.0, E_sum*E_sum - (px*px + py*py + pz*pz)));

        auto* h = new Hadron(x, y, z, px, py, pz, std::round(baryonSum), dist);
        h->SetMass(invMass);
        h->AddConstituentID(a->UniqueID());
        h->AddConstituentID(b->UniqueID());
        h->AddConstituentID(c->UniqueID());
        hadrons.push_back(h);
        a->MarkUsed();
        b->MarkUsed();
        c->MarkUsed();
        formed = true;
        break;
      }
      if (formed) break;
    }
  }

  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}