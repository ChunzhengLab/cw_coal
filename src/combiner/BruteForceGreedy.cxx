#include "Combiners.h"
#include <algorithm>
#include <cmath>
#include <random>

class PartonNeighborSearcherBrute {
public:
  PartonNeighborSearcherBrute(const std::vector<Parton*>& partons);
  std::vector<std::pair<Parton*, double>> nearestNeighbors(const Parton* query, size_t maxNeighbors = 50) const;
private:
  const std::vector<Parton*>& m_partons;
};

std::vector<Hadron*> BruteForceGreedy::Combine(const std::vector<Parton*>& partons) {
  std::vector<Hadron*> hadrons;
  PartonNeighborSearcherBrute searcher(partons);

  // Baryon preference factor m_r: controls probability to reject meson
  static std::mt19937 gen(std::random_device{}());
  std::bernoulli_distribution rejectMesonDist(m_r / (1.0 + m_r));

  for (auto* a : partons) {
    if (a->IsUsed()) continue;
    bool matched = false;

    auto neighbors = searcher.nearestNeighbors(a);
    for (auto& [b, dist] : neighbors) {
      if (b->IsUsed()) continue;

      int sum = std::round(a->GetBaryonNumber() + b->GetBaryonNumber());
      if (sum == 0) {
        // Probabilistically reject some meson combinations based on m_r
        if (rejectMesonDist(gen)) {
          continue;  // skip this meson pairing
        }
        // Kinematic calculation
        double x_mid = (a->X() + b->X()) / 2;
        double y_mid = (a->Y() + b->Y()) / 2;
        double z_mid = (a->Z() + b->Z()) / 2;
        double px = a->Px() + b->Px();
        double py = a->Py() + b->Py();
        double pz = a->Pz() + b->Pz();
        double m1 = a->GetMassFromPDG();
        double m2 = b->GetMassFromPDG();
        double E1 = std::sqrt(a->Px()*a->Px() + a->Py()*a->Py() + a->Pz()*a->Pz() + m1 * m1);
        double E2 = std::sqrt(b->Px()*b->Px() + b->Py()*b->Py() + b->Pz()*b->Pz() + m2 * m2);
        double E_sum = E1 + E2;
        double invMass = std::sqrt(std::max(0.0, E_sum * E_sum - (px*px + py*py + pz*pz)));
        auto* h = new Hadron(x_mid, y_mid, z_mid, px, py, pz, sum, dist);
        h->SetMass(invMass);
        h->AddConstituentID(a->UniqueID());
        h->AddConstituentID(b->UniqueID());
        hadrons.push_back(h);
        a->MarkUsed(); b->MarkUsed();
        matched = true;
        break;
      }

      for (auto* c : partons) {
        if (c == a || c == b || c->IsUsed()) continue;
        double baryonSum = a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber();
        if (std::round(baryonSum) != 1 && std::round(baryonSum) != -1) continue;

        // formation distance for baryon
        double formationDist = a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c);

        // Kinematic calculation for baryon
        double x_mid = (a->X() + b->X() + c->X()) / 3;
        double y_mid = (a->Y() + b->Y() + c->Y()) / 3;
        double z_mid = (a->Z() + b->Z() + c->Z()) / 3;
        double px = a->Px() + b->Px() + c->Px();
        double py = a->Py() + b->Py() + c->Py();
        double pz = a->Pz() + b->Pz() + c->Pz();
        double m1 = a->GetMassFromPDG();
        double m2 = b->GetMassFromPDG();
        double m3 = c->GetMassFromPDG();
        double E1 = std::sqrt(a->Px()*a->Px() + a->Py()*a->Py() + a->Pz()*a->Pz() + m1*m1);
        double E2 = std::sqrt(b->Px()*b->Px() + b->Py()*b->Py() + b->Pz()*b->Pz() + m2*m2);
        double E3 = std::sqrt(c->Px()*c->Px() + c->Py()*c->Py() + c->Pz()*c->Pz() + m3*m3);
        double E_sum = E1 + E2 + E3;
        double invMass = std::sqrt(std::max(0.0, E_sum*E_sum - (px*px + py*py + pz*pz)));
        auto* h = new Hadron(x_mid, y_mid, z_mid, px, py, pz, std::round(baryonSum), formationDist);
        h->SetMass(invMass);
        h->AddConstituentID(a->UniqueID());
        h->AddConstituentID(b->UniqueID());
        h->AddConstituentID(c->UniqueID());
        hadrons.push_back(h);
        a->MarkUsed(); b->MarkUsed(); c->MarkUsed();
        matched = true;
        break;  // break out of c-loop
      }
      if (matched) {
          break;  // break out of b-loop
      }
    }
  }

  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}

PartonNeighborSearcherBrute::PartonNeighborSearcherBrute(const std::vector<Parton*>& partons)
  : m_partons(partons) {}

std::vector<std::pair<Parton*, double>>
PartonNeighborSearcherBrute::nearestNeighbors(const Parton* query, size_t maxResults) const {
  std::vector<std::pair<Parton*, double>> results;

  for (auto* a : m_partons) {
    if (a == query || a->IsUsed()) continue;
    double dist = query->DistanceTo(*a);
    results.emplace_back(a, dist);
  }

  std::sort(results.begin(), results.end(),
            [](const auto& a, const auto& b) { return a.second < b.second; });

  if (results.size() > maxResults)
    results.resize(maxResults);
  return results;
}
