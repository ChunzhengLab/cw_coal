#include "CombinerBase.h"
#include "Particle.h"
#include <algorithm>
#include <vector>


std::vector<Hadron*> CombinerBase::Afterburner(const std::vector<Parton*>& partons) {
    std::vector<Parton*> unused;
    for (auto* p : partons) {
        if (!p->IsUsed()) {
            unused.push_back(p);
        }
    }

    std::vector<Hadron*> result;

    std::vector<Parton*> plus, minus;
    for (auto* p : unused) {
        if (p->GetBaryonNumber() > 0) plus.push_back(p);
        else minus.push_back(p);
    }

    size_t mesonCount = std::min(plus.size(), minus.size());
    for (size_t i = 0; i < mesonCount; ++i) {
        Parton* a = plus[i];
        Parton* b = minus[i];

        double x = (a->X() + b->X()) / 2;
        double y = (a->Y() + b->Y()) / 2;
        double z = (a->Z() + b->Z()) / 2;
        double px = a->Px() + b->Px();
        double py = a->Py() + b->Py();
        double pz = a->Pz() + b->Pz();
        double formation = a->DistanceTo(*b);

        Hadron* h = new Hadron(x, y, z, px, py, pz, 0, formation);
        h->AddConstituentID(a->UniqueID());
        h->AddConstituentID(b->UniqueID());
        a->MarkUsed();
        b->MarkUsed();
        result.push_back(h);
    }

    std::vector<Parton*> remain;
    for (auto* p : unused) {
        if (!p->IsUsed()) remain.push_back(p);
    }

    for (size_t i = 0; i + 2 < remain.size(); i += 3) {
        Parton* a = remain[i];
        Parton* b = remain[i + 1];
        Parton* c = remain[i + 2];

        double x = (a->X() + b->X() + c->X()) / 3;
        double y = (a->Y() + b->Y() + c->Y()) / 3;
        double z = (a->Z() + b->Z() + c->Z()) / 3;
        double px = a->Px() + b->Px() + c->Px();
        double py = a->Py() + b->Py() + c->Py();
        double pz = a->Pz() + b->Pz() + c->Pz();
        double formation = a->DistanceTo(*b) + a->DistanceTo(*c) + b->DistanceTo(*c);
        int baryonNumber = std::round(a->GetBaryonNumber() + b->GetBaryonNumber() + c->GetBaryonNumber());

        Hadron* h = new Hadron(x, y, z, px, py, pz, baryonNumber, formation);
        h->AddConstituentID(a->UniqueID());
        h->AddConstituentID(b->UniqueID());
        h->AddConstituentID(c->UniqueID());
        a->MarkUsed();
        b->MarkUsed();
        c->MarkUsed();
        result.push_back(h);
    }

    return result;
}