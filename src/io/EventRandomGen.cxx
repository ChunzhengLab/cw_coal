#include "io/EventRandomGen.h"
#include "core/PhysicsConstants.h"
#include "core/Event.h"
#include <cstdlib>
#include <cmath>
#include "core/Particle.h"

EventRandomGen::EventRandomGen(const std::string& histFilePath)
  : m_histFilePath(
      std::getenv("CW_COAL_PARTON_HIST")
        ? std::getenv("CW_COAL_PARTON_HIST")
        : histFilePath
    )
{
}

Parton* EventRandomGen::RandomParton() const {
  return Parton::RandomFromHists(m_histFilePath.c_str());
}

Event EventRandomGen::GenerateEvent(int nParts, int sumBaryonNumber) const {
    int parts = nParts;
    if (parts < 0) {
        // 如果未指定，则根据多重性直方图确定初始 Parton 数量
        parts = static_cast<int>(
            PhysicsConstants::GetMultiplicityHistogram().GetRandom()
        );
    }

    Event event;
    std::vector<Parton*> partons;
    partons.reserve(parts);
    int baryonSum3 = 0; // 以 1/3 单位累加

    // 先抽取 parts 个 Parton
    for (int i = 0; i < parts; ++i) {
        Parton* p = RandomParton();
        int bn3 = static_cast<int>(std::round(p->GetBaryonNumber() * 3.0));
        baryonSum3 += bn3;
        partons.push_back(p);
    }

    // 目标以 1/3 单位累加
    int targetSum3 = sumBaryonNumber * 3;

    // 继续抽取或丢弃，直到达到目标 baryonSum3
    while (baryonSum3 != targetSum3) {
        Parton* p = RandomParton();
        int bn3 = static_cast<int>(std::round(p->GetBaryonNumber() * 3.0));

        if ((baryonSum3 < targetSum3 && bn3 > 0) ||
            (baryonSum3 > targetSum3 && bn3 < 0)) {
            baryonSum3 += bn3;
            partons.push_back(p);
        } else {
            delete p;
        }
    }

    // 将 Parton 添加到 Event 中
    for (auto* p : partons) {
        event.AddParton(p);
    }

    return event;
}