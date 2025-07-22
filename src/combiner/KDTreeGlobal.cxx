#include <algorithm>
#include <cmath>
#include <set>

#include "Combiners.h"
#include "core/PartonKDTree.h"

struct Candidate {
  double distance;
  bool isBaryon;
  Parton* a;
  Parton* b;
  Parton* c;

  Candidate(double d, Parton* x, Parton* y) : distance(d), isBaryon(false), a(x), b(y), c(nullptr) {
  }

  Candidate(double d, Parton* x, Parton* y, Parton* z) : distance(d), isBaryon(true), a(x), b(y), c(z) {
  }

  bool operator<(const Candidate& other) const {
    return distance < other.distance;
  }
};

std::vector<Hadron*> KDTreeGlobal::Combine(const std::vector<Parton*>& partons)
{
  const double rInvScale = m_r;                // keep existing scale factor

  std::vector<Hadron*> hadrons;                // final output
  std::set<Parton*>    leftover;               // scheme‑A rolling buffer

  // ------------------------------------------------------------------
  // 1. build time frames using TimeFrameManager
  // ------------------------------------------------------------------
  if (partons.empty()) return hadrons;
  timeFrameManager_->BuildFrames(partons);
  const auto& timeFrames = timeFrameManager_->GetFrameBoundaries();

  // ------------------------------------------------------------------
  // 2. helper lambda: build & consume candidates inside a frame
  // ------------------------------------------------------------------
  auto combineFrame = [&](const std::vector<Parton*>& frameParts)
  {
    std::vector<Candidate> candidates;
    PartonKDTree tree(frameParts);

    // ----- meson candidates -----
    for (auto* a : frameParts) {
      if (a->IsUsed()) continue;
      auto neighbors = tree.FindNeighbors(a, 20);
      for (auto& [b, dist] : neighbors) {
        if (a == b || b->IsUsed()) continue;
        if (std::round(a->GetBaryonNumber() + b->GetBaryonNumber()) != 0) continue;
        candidates.emplace_back(dist, a, b);
      }
    }

    // ----- baryon candidates -----
    for (auto* a : frameParts) {
      if (a->IsUsed()) continue;
      auto neighbors = tree.FindNeighbors(a, 10);
      for (size_t i = 0; i < neighbors.size(); ++i) {
        Parton* b = neighbors[i].first;
        if (a == b || b->IsUsed()) continue;
        for (size_t j = i + 1; j < neighbors.size(); ++j) {
          Parton* c = neighbors[j].first;
          if (c == a || c == b || c->IsUsed()) continue;
          double baryonNumber = a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber();
          if (std::round(baryonNumber) != 1 && std::round(baryonNumber) != -1) continue;
          double dist = (a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c)) / rInvScale;
          candidates.emplace_back(dist, a, b, c);
        }
      }
    }

    std::sort(candidates.begin(), candidates.end());

    // ----- consume candidates -----
    for (const auto& candi : candidates) {
      if (candi.isBaryon) {
        if (candi.a->IsUsed() || candi.b->IsUsed() || candi.c->IsUsed()) continue;
        // build baryon
        double m1 = candi.a->GetMassFromPDG();
        double m2 = candi.b->GetMassFromPDG();
        double m3 = candi.c->GetMassFromPDG();
        double E1 = std::hypot(candi.a->Px(), std::hypot(candi.a->Py(), candi.a->Pz()), m1);
        double E2 = std::hypot(candi.b->Px(), std::hypot(candi.b->Py(), candi.b->Pz()), m2);
        double E3 = std::hypot(candi.c->Px(), std::hypot(candi.c->Py(), candi.c->Pz()), m3);
        double px = candi.a->Px() + candi.b->Px() + candi.c->Px();
        double py = candi.a->Py() + candi.b->Py() + candi.c->Py();
        double pz = candi.a->Pz() + candi.b->Pz() + candi.c->Pz();
        double Etot = E1 + E2 + E3;
        double p2   = px * px + py * py + pz * pz;
        double invM = (Etot * Etot > p2) ? std::sqrt(Etot * Etot - p2) : 0.0;
        double x = (candi.a->X() + candi.b->X() + candi.c->X()) / 3;
        double y = (candi.a->Y() + candi.b->Y() + candi.c->Y()) / 3;
        double z = (candi.a->Z() + candi.b->Z() + candi.c->Z()) / 3;
        double formation = candi.a->DistanceTo(*candi.b) + candi.a->DistanceTo(*candi.c) +
                           candi.b->DistanceTo(*candi.c);

        auto* h = new Hadron(x, y, z, px, py, pz,
                             std::round(candi.a->GetBaryonNumber() + candi.b->GetBaryonNumber() +
                                        candi.c->GetBaryonNumber()),
                             formation);
        h->SetMass(invM);
        h->AddConstituentID(candi.a->UniqueID());
        h->AddConstituentID(candi.b->UniqueID());
        h->AddConstituentID(candi.c->UniqueID());
        hadrons.push_back(h);
        candi.a->MarkUsed();
        candi.b->MarkUsed();
        candi.c->MarkUsed();
      } else {
        if (candi.a->IsUsed() || candi.b->IsUsed()) continue;
        // build meson
        double m1 = candi.a->GetMassFromPDG();
        double m2 = candi.b->GetMassFromPDG();
        double E1 = std::hypot(candi.a->Px(), std::hypot(candi.a->Py(), candi.a->Pz()), m1);
        double E2 = std::hypot(candi.b->Px(), std::hypot(candi.b->Py(), candi.b->Pz()), m2);
        double px = candi.a->Px() + candi.b->Px();
        double py = candi.a->Py() + candi.b->Py();
        double pz = candi.a->Pz() + candi.b->Pz();
        double Etot = E1 + E2;
        double p2   = px * px + py * py + pz * pz;
        double invM = (Etot * Etot > p2) ? std::sqrt(Etot * Etot - p2) : 0.0;
        double x = (candi.a->X() + candi.b->X()) / 2;
        double y = (candi.a->Y() + candi.b->Y()) / 2;
        double z = (candi.a->Z() + candi.b->Z()) / 2;
        double formation = candi.a->DistanceTo(*candi.b);

        auto* h = new Hadron(x, y, z, px, py, pz,
                             std::round(candi.a->GetBaryonNumber() + candi.b->GetBaryonNumber()),
                             formation);
        h->SetMass(invM);
        h->AddConstituentID(candi.a->UniqueID());
        h->AddConstituentID(candi.b->UniqueID());
        hadrons.push_back(h);
        candi.a->MarkUsed();
        candi.b->MarkUsed();
      }
    }
  };
  // ------------------------------------------------------------------
  // 3. loop over frames with rolling leftover strategy (scheme A)
  // ------------------------------------------------------------------
  for (size_t f = 0; f < timeFrameManager_->GetNumFrames(); ++f) {
    // Get partons in current frame (including leftovers)
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

  // => scheme says we ignore any final leftovers (remain deconfined) or you
  //    can push them to afterburner; here we simply call Afterburner as before.
  auto afterburned = Afterburner(partons);
  hadrons.insert(hadrons.end(), afterburned.begin(), afterburned.end());

  return hadrons;
}
