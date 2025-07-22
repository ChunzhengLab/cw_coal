#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_set>

#include "Combiners.h"

struct Candidate {
  double distance;
  bool isBaryon;
  Parton* a;
  Parton* b;
  Parton* c;

  Candidate(double dist, Parton* x, Parton* y) : distance(dist), isBaryon(false), a(x), b(y), c(nullptr) {
  }

  Candidate(double dist, Parton* x, Parton* y, Parton* z) : distance(dist), isBaryon(true), a(x), b(y), c(z) {
  }

  bool operator<(const Candidate& other) const {
    return distance < other.distance;
  }
};

std::vector<Hadron*> ExhaustiveSorted::Combine(const std::vector<Parton*>& partons) {
  std::vector<Hadron*> hadrons;
  if (partons.empty()) return hadrons;
  
  // Build time frames using TimeFrameManager
  timeFrameManager_->BuildFrames(partons);
  std::set<Parton*> leftover; // Rolling buffer for unused partons
  
  // Helper lambda: combine partons within a single frame
  auto combineFrame = [&](const std::vector<Parton*>& frameParts) {
    std::vector<Candidate> candidates;

    // Generate meson candidates
    for (size_t i = 0; i < frameParts.size(); ++i) {
      if (frameParts[i]->IsUsed()) continue;
      for (size_t j = i + 1; j < frameParts.size(); ++j) {
        if (frameParts[j]->IsUsed()) continue;

        int bn = std::round(frameParts[i]->GetBaryonNumber() + frameParts[j]->GetBaryonNumber());
        if (bn != 0) continue;

        double dist = frameParts[i]->DistanceTo(*frameParts[j]);
        candidates.emplace_back(dist, frameParts[i], frameParts[j]);
      }
    }

    // Generate baryon candidates
    for (size_t i = 0; i < frameParts.size(); ++i) {
      if (frameParts[i]->IsUsed()) continue;
      for (size_t j = i + 1; j < frameParts.size(); ++j) {
        if (frameParts[j]->IsUsed()) continue;
        for (size_t k = j + 1; k < frameParts.size(); ++k) {
          if (frameParts[k]->IsUsed()) continue;

          double bn = frameParts[i]->GetBaryonNumber() + frameParts[j]->GetBaryonNumber() + frameParts[k]->GetBaryonNumber();
          if (std::round(bn) != 1 && std::round(bn) != -1) continue;

          double dist = frameParts[i]->DistanceTo(*frameParts[j]) + frameParts[i]->DistanceTo(*frameParts[k]) +
                        frameParts[j]->DistanceTo(*frameParts[k]);
          dist = dist / m_r;
          candidates.emplace_back(dist, frameParts[i], frameParts[j], frameParts[k]);
        }
      }
    }

    std::sort(candidates.begin(), candidates.end());

    // Process candidates
    for (const auto& candi : candidates) {
      if (candi.isBaryon) {
        if (candi.a->IsUsed() || candi.b->IsUsed() || candi.c->IsUsed()) continue;
        
        double bn = candi.a->GetBaryonNumber() + candi.b->GetBaryonNumber() + candi.c->GetBaryonNumber();
        double px = candi.a->Px() + candi.b->Px() + candi.c->Px();
        double py = candi.a->Py() + candi.b->Py() + candi.c->Py();
        double pz = candi.a->Pz() + candi.b->Pz() + candi.c->Pz();
        double x = (candi.a->X() + candi.b->X() + candi.c->X()) / 3;
        double y = (candi.a->Y() + candi.b->Y() + candi.c->Y()) / 3;
        double z = (candi.a->Z() + candi.b->Z() + candi.c->Z()) / 3;

        double m1 = candi.a->GetMassFromPDG();
        double m2 = candi.b->GetMassFromPDG();
        double m3 = candi.c->GetMassFromPDG();
        double E1 = std::hypot(candi.a->Px(), std::hypot(candi.a->Py(), candi.a->Pz()), m1);
        double E2 = std::hypot(candi.b->Px(), std::hypot(candi.b->Py(), candi.b->Pz()), m2);
        double E3 = std::hypot(candi.c->Px(), std::hypot(candi.c->Py(), candi.c->Pz()), m3);
        double E_sum = E1 + E2 + E3;
        double p2 = px * px + py * py + pz * pz;
        double invMass = (E_sum * E_sum > p2) ? std::sqrt(E_sum * E_sum - p2) : 0.0;

        double rawDist = candi.a->DistanceTo(*candi.b) + candi.a->DistanceTo(*candi.c) + candi.b->DistanceTo(*candi.c);
        auto* h = new Hadron(x, y, z, px, py, pz, std::round(bn), rawDist);
        h->SetMass(invMass);
        h->AddConstituentID(candi.a->UniqueID());
        h->AddConstituentID(candi.b->UniqueID());
        h->AddConstituentID(candi.c->UniqueID());
        hadrons.push_back(h);
        candi.a->MarkUsed();
        candi.b->MarkUsed();
        candi.c->MarkUsed();
      } else {
        if (candi.a->IsUsed() || candi.b->IsUsed()) continue;
        
        double px = candi.a->Px() + candi.b->Px();
        double py = candi.a->Py() + candi.b->Py();
        double pz = candi.a->Pz() + candi.b->Pz();
        double x = (candi.a->X() + candi.b->X()) / 2;
        double y = (candi.a->Y() + candi.b->Y()) / 2;
        double z = (candi.a->Z() + candi.b->Z()) / 2;

        double m1 = candi.a->GetMassFromPDG();
        double m2 = candi.b->GetMassFromPDG();
        double E1 = std::hypot(candi.a->Px(), std::hypot(candi.a->Py(), candi.a->Pz()), m1);
        double E2 = std::hypot(candi.b->Px(), std::hypot(candi.b->Py(), candi.b->Pz()), m2);
        double E_sum = E1 + E2;
        double p2 = px * px + py * py + pz * pz;
        double invMass = (E_sum * E_sum > p2) ? std::sqrt(E_sum * E_sum - p2) : 0.0;

        double rawDist = candi.a->DistanceTo(*candi.b);
        auto* h = new Hadron(x, y, z, px, py, pz, 0, rawDist);
        h->SetMass(invMass);
        h->AddConstituentID(candi.a->UniqueID());
        h->AddConstituentID(candi.b->UniqueID());
        hadrons.push_back(h);
        candi.a->MarkUsed();
        candi.b->MarkUsed();
      }
    }
  };

  // Loop over time frames with rolling leftover strategy
  for (size_t f = 0; f < timeFrameManager_->GetNumFrames(); ++f) {
    // Get partons in current frame
    std::vector<Parton*> frameParts = timeFrameManager_->GetPartonsInFrame(partons, f);
    
    // Add leftovers from previous frame
    for (auto* p : leftover) {
      if (!p->IsUsed()) frameParts.push_back(p);
    }
    leftover.clear();

    if (frameParts.empty()) continue;
    combineFrame(frameParts);

    // Collect new leftovers and move them to next timeframe
    for (auto* p : frameParts) {
      if (!p->IsUsed()) {
        leftover.insert(p);
      }
    }
    
    // Move leftover partons to next timeframe
    timeFrameManager_->MovePartonsToNextFrame(std::vector<Parton*>(leftover.begin(), leftover.end()), f);
  }

  // Final afterburner pass
  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}