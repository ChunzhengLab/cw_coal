#include "Combiners.h"
#include <cmath>
#include <unordered_set>
#include <algorithm>

struct Candidate {
    double distance;
    bool isBaryon;
    size_t idxA, idxB, idxC;

    Candidate(double dist, bool isB, size_t a, size_t b)
        : distance(dist), isBaryon(isB), idxA(a), idxB(b), idxC(0) {}

    Candidate(double dist, bool isB, size_t a, size_t b, size_t c)
        : distance(dist), isBaryon(isB), idxA(a), idxB(b), idxC(c) {}

    bool operator<(const Candidate& other) const { return distance < other.distance; }
};

std::vector<Hadron*> BruteForceGlobal::Combine(const std::vector<Parton*>& partons) {
    std::vector<Candidate> candidates;
    std::unordered_set<size_t> used;
    std::vector<Hadron*> hadrons;

    for (size_t i = 0; i < partons.size(); ++i) {
        if (partons[i]->IsUsed()) continue;
        for (size_t j = i + 1; j < partons.size(); ++j) {
            if (partons[j]->IsUsed()) continue;

            int bn = std::round(partons[i]->GetBaryonNumber() + partons[j]->GetBaryonNumber());
            if (bn != 0) continue;

            double dist = partons[i]->DistanceTo(*partons[j]);
            candidates.emplace_back(dist, false, i, j);
        }
    }

    for (size_t i = 0; i < partons.size(); ++i) {
        if (partons[i]->IsUsed()) continue;
        for (size_t j = i + 1; j < partons.size(); ++j) {
            if (partons[j]->IsUsed()) continue;
            for (size_t k = j + 1; k < partons.size(); ++k) {
                if (partons[k]->IsUsed()) continue;

                double bn = partons[i]->GetBaryonNumber() + partons[j]->GetBaryonNumber() + partons[k]->GetBaryonNumber();
                if (std::round(bn) != 1 && std::round(bn) != -1) continue;

                double dist = partons[i]->DistanceTo(*partons[j]) +
                              partons[i]->DistanceTo(*partons[k]) +
                              partons[j]->DistanceTo(*partons[k]);
                dist = dist / m_r;
                candidates.emplace_back(dist, true, i, j, k);
            }
        }
    }

    std::sort(candidates.begin(), candidates.end());

    for (const auto& candi : candidates) {
        if (candi.isBaryon) {
            if (used.count(candi.idxA) || used.count(candi.idxB) || used.count(candi.idxC)) continue;
            Parton *a = partons[candi.idxA], *b = partons[candi.idxB], *c = partons[candi.idxC];
            double bn = a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber();
            double px = a->Px() + b->Px() + c->Px();
            double py = a->Py() + b->Py() + c->Py();
            double pz = a->Pz() + b->Pz() + c->Pz();
            double x = (a->X() + b->X() + c->X()) / 3;
            double y = (a->Y() + b->Y() + c->Y()) / 3;
            double z = (a->Z() + b->Z() + c->Z()) / 3;

            double m1 = partons[candi.idxA]->GetMassFromPDG();
            double m2 = partons[candi.idxB]->GetMassFromPDG();
            double m3 = partons[candi.idxC]->GetMassFromPDG();
            double E1 = std::sqrt(partons[candi.idxA]->Px()*partons[candi.idxA]->Px()
                        + partons[candi.idxA]->Py()*partons[candi.idxA]->Py()
                        + partons[candi.idxA]->Pz()*partons[candi.idxA]->Pz()
                        + m1*m1);
            double E2 = std::sqrt(partons[candi.idxB]->Px()*partons[candi.idxB]->Px()
                        + partons[candi.idxB]->Py()*partons[candi.idxB]->Py()
                        + partons[candi.idxB]->Pz()*partons[candi.idxB]->Pz()
                        + m2*m2);
            double E3 = std::sqrt(partons[candi.idxC]->Px()*partons[candi.idxC]->Px()
                        + partons[candi.idxC]->Py()*partons[candi.idxC]->Py()
                        + partons[candi.idxC]->Pz()*partons[candi.idxC]->Pz()
                        + m3*m3);
            double E_sum = E1 + E2 + E3;
            double invMass = std::sqrt(std::max(0.0, E_sum*E_sum - (px*px + py*py + pz*pz)));

            double rawDist = a->DistanceTo(*b) 
                           + a->DistanceTo(*c) 
                           + b->DistanceTo(*c);
            auto* h = new Hadron(x, y, z, px, py, pz, std::round(bn), rawDist);
            h->SetMass(invMass);
            h->AddConstituentID(a->UniqueID());
            h->AddConstituentID(b->UniqueID());
            h->AddConstituentID(c->UniqueID());
            hadrons.push_back(h);
            used.insert(candi.idxA); used.insert(candi.idxB); used.insert(candi.idxC);
            a->MarkUsed(); b->MarkUsed(); c->MarkUsed();
        } else {
            if (used.count(candi.idxA) || used.count(candi.idxB)) continue;
            Parton *a = partons[candi.idxA], *b = partons[candi.idxB];
            double px = a->Px() + b->Px();
            double py = a->Py() + b->Py();
            double pz = a->Pz() + b->Pz();
            double x = (a->X() + b->X()) / 2;
            double y = (a->Y() + b->Y()) / 2;
            double z = (a->Z() + b->Z()) / 2;

            double m1 = partons[candi.idxA]->GetMassFromPDG();
            double m2 = partons[candi.idxB]->GetMassFromPDG();
            double E1 = std::sqrt(partons[candi.idxA]->Px()*partons[candi.idxA]->Px()
                        + partons[candi.idxA]->Py()*partons[candi.idxA]->Py()
                        + partons[candi.idxA]->Pz()*partons[candi.idxA]->Pz()
                        + m1*m1);
            double E2 = std::sqrt(partons[candi.idxB]->Px()*partons[candi.idxB]->Px()
                        + partons[candi.idxB]->Py()*partons[candi.idxB]->Py()
                        + partons[candi.idxB]->Pz()*partons[candi.idxB]->Pz()
                        + m2*m2);
            double E_sum = E1 + E2;
            double invMass = std::sqrt(std::max(0.0, E_sum*E_sum - (px*px + py*py + pz*pz)));

            double rawDist = a->DistanceTo(*b);
            auto* h = new Hadron(x, y, z, px, py, pz, 0, rawDist);
            h->SetMass(invMass);
            h->AddConstituentID(a->UniqueID());
            h->AddConstituentID(b->UniqueID());
            hadrons.push_back(h);
            used.insert(candi.idxA); used.insert(candi.idxB);
            a->MarkUsed(); b->MarkUsed();
        }
    }

    return hadrons;
}
