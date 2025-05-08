#include "io/EventRandomGen.h"
#include "core/PhysicsConstants.h"
#include "core/Event.h"
#include <cstdlib>
#include <cmath>
#include "core/Particle.h"
#include <vector>

EventRandomGen::EventRandomGen(const std::string& histFilePath)
  : m_histFilePath(
      std::getenv("CW_COAL_PARTON_HIST")
        ? std::getenv("CW_COAL_PARTON_HIST")
        : histFilePath
    )
{
}

void EventRandomGen::GenerateEvent(Event& out, int nParts, int sumBaryonNumber,
                                   SamplingMode mode) const {
    out.Reset();

    // Determine number of partons
    int parts = (nParts < 0)
        ? static_cast<int>(PhysicsConstants::GetMultiplicityHistogram().GetRandom())
        : nParts;

    std::vector<Parton*> partons;
    partons.reserve(parts);

    // Sampler lambda: Toy or file-based
    auto sampler = [this, mode]() {
        return (mode == kToyMode)
            ? Parton::Random(nullptr)
            : Parton::RandomFromHists(m_histFilePath.c_str());
    };

    // Initial sample
    int baryonSum3 = 0;
    for (int i = 0; i < parts; ++i) {
        Parton* p = sampler();
        partons.push_back(p);
        baryonSum3 += static_cast<int>(std::round(p->GetBaryonNumber() * 3.0));
    }

    // Adjust to match desired baryon number
    int target3 = sumBaryonNumber * 3;
    while (baryonSum3 != target3) {
        Parton* p = sampler();
        int bn3 = static_cast<int>(std::round(p->GetBaryonNumber() * 3.0));
        if ((baryonSum3 < target3 && bn3 > 0) ||
            (baryonSum3 > target3 && bn3 < 0)) {
            baryonSum3 += bn3;
            partons.push_back(p);
        } else {
            delete p;
        }
    }

    // Add to event
    for (auto* p : partons) {
        out.AddParton(p);
    }
}
