#ifndef PIDASSIGNER_H
#define PIDASSIGNER_H

#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <vector>

#include "core/Event.h"
#include "core/PIDInference.h"
#include "core/Particle.h"

class PIDAssigner {
 public:
  /**
   * @brief Assigns PDG codes to all hadrons in the event.
   *
   * This will iterate over each hadron in the event, infer its PDG code
   * using PIDInference, and update the hadron's PID.
   *
   * @param event The Event containing hadrons with quark composition and mass.
   */
  static void Assign(Event& event) {
    // 获取所有 hadron 指针
    const std::vector<Hadron*>& hadronList = event.GetHadrons();

    // Build UID -> Parton* map
    std::unordered_map<unsigned int, Parton*> uidToPartonMap;
    for (Parton* parton : event.GetPartons()) { uidToPartonMap[parton->UniqueID()] = parton; }

    // Phase 1: Assign non-diagonal mesons and all baryons
    for (Hadron* hadronPtr : hadronList) {
      auto constituentUIDs = hadronPtr->GetConstituentIDs();
      std::vector<int> quarkFlavorCodes;
      quarkFlavorCodes.reserve(constituentUIDs.size());
      for (auto uid : constituentUIDs) { quarkFlavorCodes.push_back(uidToPartonMap[uid]->GetPID()); }
      // diagonal light mesons are size==2 & q0==-q1 & |q|<=2
      if (!(quarkFlavorCodes.size() == 2 && quarkFlavorCodes[0] == -quarkFlavorCodes[1] &&
            std::abs(quarkFlavorCodes[0]) <= 2)) {
        int pdg = PIDInference::InferPID(quarkFlavorCodes, hadronPtr->GetMass());
        hadronPtr->SetPID(pdg);
      }
    }

    // Count charged pions (±211) and rhos (±213)
    int numChargedPions = 0;
    int numChargedRhos = 0;
    for (Hadron* hadronPtr : hadronList) {
      int pid = hadronPtr->GetPID();
      int absPid = std::abs(pid);
      if (absPid == 211) ++numChargedPions;
      if (absPid == 213) ++numChargedRhos;
    }

    // Phase 2: Batch assign diagonal light mesons using event counts
    std::vector<size_t> diagonalMesonIndices;
    std::vector<double> diagonalMesonMasses;
    for (size_t i = 0; i < hadronList.size(); ++i) {
      Hadron* hadronPtr = hadronList[i];
      auto constituentUIDs = hadronPtr->GetConstituentIDs();
      std::vector<int> quarkFlavorCodes;
      quarkFlavorCodes.reserve(constituentUIDs.size());
      for (auto uid : constituentUIDs) { quarkFlavorCodes.push_back(uidToPartonMap[uid]->GetPID()); }
      if (quarkFlavorCodes.size() == 2 && quarkFlavorCodes[0] == -quarkFlavorCodes[1] &&
          std::abs(quarkFlavorCodes[0]) <= 2) {
        diagonalMesonIndices.push_back(i);
        diagonalMesonMasses.push_back(hadronPtr->GetMass());
      }
    }
    if (!diagonalMesonIndices.empty()) {
      std::vector<int> diagonalMesonPDGs;
      PIDInference::BatchAssignDiagonalLightMesons(diagonalMesonMasses, diagonalMesonPDGs, numChargedPions,
                                                   numChargedRhos);
      for (size_t k = 0; k < diagonalMesonIndices.size(); ++k) {
        Hadron* hadronPtr = hadronList[diagonalMesonIndices[k]];
        hadronPtr->SetPID(diagonalMesonPDGs[k]);
      }
    }
  }
};

#endif  // PIDASSIGNER_H
