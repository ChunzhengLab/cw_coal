#include "Combiners.h"
#include <algorithm>
#include <cmath>

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

  for (auto* a : partons) {
    if (a->IsUsed()) continue;

    auto neighbors = searcher.nearestNeighbors(a);
    for (auto& [b, dist] : neighbors) {
      if (b->IsUsed()) continue;

      int sum = std::round(a->GetBaryonNumber() + b->GetBaryonNumber());
      if (sum == 0) {
        auto* h = new Hadron(
          (a->X() + b->X()) / 2, (a->Y() + b->Y()) / 2, (a->Z() + b->Z()) / 2,
          a->Px() + b->Px(), a->Py() + b->Py(), a->Pz() + b->Pz(),
          sum, a->DistanceTo(*b)
        );
        hadrons.push_back(h);
        a->MarkUsed(); b->MarkUsed();
        break;
      }

      for (auto* c : partons) {
        if (c == a || c == b || c->IsUsed()) continue;
        double baryonSum = a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber();
        if (std::round(baryonSum) != 1 && std::round(baryonSum) != -1) continue;

        auto* h = new Hadron(
          (a->X() + b->X() + c->X()) / 3, (a->Y() + b->Y() + c->Y()) / 3, (a->Z() + b->Z() + c->Z()) / 3,
          a->Px() + b->Px() + c->Px(), a->Py() + b->Py() + c->Py(), a->Pz() + b->Pz() + c->Pz(),
          std::round(baryonSum),
          a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c)
        );
        hadrons.push_back(h);
        a->MarkUsed(); b->MarkUsed(); c->MarkUsed();
        goto next_parton;
      }
    }
  next_parton:;
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
