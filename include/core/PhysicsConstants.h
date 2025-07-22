#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

#include <optional>
#include <unordered_map>
#include <vector>

#include "TH1.h"
#include "TPDGCode.h"  // Correct path for ROOT's TPDGCode

namespace PhysicsConstants {
inline const std::unordered_map<int, double>& GetMassTable() {
  static const std::unordered_map<int, double> massTable = {
      {kDown, 0.325},
      {kUp, 0.325},
      {kStrange, 0.5},
      {kCharm, 1.60},
      {kBottom, 5.0}, // constituent mass from pythia8
      {kTop, 172.5},
      {kGluon, 0.0},
      {kElectron, 0.000510999},
      {kPositron, 0.000510999},
      {kNuE, 0.0},
      {kNuEBar, 0.0},
      {kMuonMinus, 0.1056584},
      {kMuonPlus, 0.1056584},
      {kNuMu, 0.0},
      {kNuMuBar, 0.0},
      {kTauMinus, 1.77686},
      {kTauPlus, 1.77686},
      {kNuTau, 0.0},
      {kNuTauBar, 0.0},
      {kGamma, 0.0},
      {kZ0, 91.1876},
      {kWPlus, 80.377},
      {kWMinus, 80.377},
      {kPi0, 0.1349768},
      {kPiPlus, 0.1395704},
      {kPiMinus, 0.1395704},
      {kK0, 0.497611},
      {kK0Bar, 0.497611},
      {kK0Short, 0.497611},
      {kK0Long, 0.497611},
      {kKPlus, 0.493677},
      {kKMinus, 0.493677},
      {kProton, 0.9382721},
      {kProtonBar, 0.9382721},
      {kNeutron, 0.9395654},
      {kNeutronBar, 0.9395654},
      {kLambda0, 1.115683},
      {kLambda0Bar, 1.115683},
      {kLambda1520, 1.519},
      {kSigmaMinus, 1.197449},
      {kSigmaPlus, 1.18937},
      {kSigma0, 1.192642},
      {kXiMinus, 1.32171},
      {kXiPlusBar, 1.32171},
      {kOmegaMinus, 1.67245},
      {kOmegaPlusBar, 1.67245},

      // Charmonium states (in GeV)
      {441, 2.9839},       // ηc(1S)
      {443, 3.09690},      // J/ψ(1S)
      {10441, 3.41475},    // χc0(1P)
      {10443, 3.51066},    // χc1(1P)
      {445, 3.55620},      // χc2(1P)
      {100441, 3.63990},   // ηc(2S)
      {100443, 3.68610},   // ψ(2S)
      {30443, 3.77313},    // ψ(3770)
      {100445, 4.15300},   // ψ(4160)
      {9000443, 4.03900},  // χc0(4500) [if used]
      {9010443, 4.19100},  // χc1(4500) [if used]
      {9020443, 4.42100},  // χc2(4500) [if used]
      // Bottomonium states (in GeV)
      {551, 9.3987},     // ηb(1S)
      {553, 9.46030},    // ϒ(1S)
      {10551, 10.5794},  // ϒ(4S) mapped to code 10551
  };
  return massTable;
}

inline std::optional<double> GetMass(int pdgCode) {
  const auto& table = GetMassTable();
  auto it = table.find(pdgCode);
  if (it != table.end()) { return it->second; }
  return std::nullopt;
}

// Controls the relative likelihood of forming baryons vs mesons
constexpr double kBaryonPreferenceFactor = 1.0;

//----------------------------------------------------------------------------//
// Parameters used by PIDInference for meson and baryon formation
// Vector-to-pseudoscalar meson ratio (V/P)
constexpr double kMesonVectorToPseudoscalarRatio = 1.0 / 3.0;
// ρ⁰/π⁰ production ratio (ρ⁰/π⁰) initialization
constexpr double kRhoToPionRatio = 0.36;
// ω/ρ⁰ production ratio (ω/ρ⁰) initialization
constexpr double kOmegaToRhoRatio = 1.90;
// K*/K production ratio (K*∶K) initialization
constexpr double kKStarToKRatio = 0.50;

// Returns a reference to the predefined multiplicity histogram
inline const TH1I& GetMultiplicityHistogram() {
  static bool initialized = false;
  static TH1I mult("mult", "Multiplicity;Multiplicity;Counts", 100, 0, 20000);
  if (!initialized) {
    initialized = true;
    // Fill bin contents (example entries; add the rest similarly)
    mult.SetBinContent(18, 1);
    mult.SetBinContent(22, 1);
    mult.SetBinContent(23, 7);
    mult.SetBinContent(24, 8);
    mult.SetBinContent(25, 11);
    mult.SetBinContent(26, 5);
    mult.SetBinContent(27, 17);
    mult.SetBinContent(28, 27);
    mult.SetBinContent(29, 30);
    mult.SetBinContent(30, 67);
    mult.SetBinContent(31, 84);
    mult.SetBinContent(32, 95);
    mult.SetBinContent(33, 110);
    mult.SetBinContent(34, 163);
    mult.SetBinContent(35, 169);
    mult.SetBinContent(36, 254);
    mult.SetBinContent(37, 323);
    mult.SetBinContent(38, 333);
    mult.SetBinContent(39, 407);
    mult.SetBinContent(40, 430);
    mult.SetBinContent(41, 513);
    mult.SetBinContent(42, 541);
    mult.SetBinContent(43, 634);
    mult.SetBinContent(44, 676);
    mult.SetBinContent(45, 670);
    mult.SetBinContent(46, 754);
    mult.SetBinContent(47, 740);
    mult.SetBinContent(48, 773);
    mult.SetBinContent(49, 831);
    mult.SetBinContent(50, 777);
    mult.SetBinContent(51, 868);
    mult.SetBinContent(52, 887);
    mult.SetBinContent(53, 834);
    mult.SetBinContent(54, 831);
    mult.SetBinContent(55, 833);
    mult.SetBinContent(56, 885);
    mult.SetBinContent(57, 831);
    mult.SetBinContent(58, 774);
    mult.SetBinContent(59, 791);
    mult.SetBinContent(60, 728);
    mult.SetBinContent(61, 687);
    mult.SetBinContent(62, 604);
    mult.SetBinContent(63, 611);
    mult.SetBinContent(64, 587);
    mult.SetBinContent(65, 501);
    mult.SetBinContent(66, 473);
    mult.SetBinContent(67, 460);
    mult.SetBinContent(68, 416);
    mult.SetBinContent(69, 343);
    mult.SetBinContent(70, 332);
    mult.SetBinContent(71, 276);
    mult.SetBinContent(72, 259);
    mult.SetBinContent(73, 233);
    mult.SetBinContent(74, 199);
    mult.SetBinContent(75, 134);
    mult.SetBinContent(76, 141);
    mult.SetBinContent(77, 124);
    mult.SetBinContent(78, 89);
    mult.SetBinContent(79, 71);
    mult.SetBinContent(80, 83);
    mult.SetBinContent(81, 56);
    mult.SetBinContent(82, 47);
    mult.SetBinContent(83, 41);
    mult.SetBinContent(84, 32);
    mult.SetBinContent(85, 26);
    mult.SetBinContent(86, 23);
    mult.SetBinContent(87, 14);
    mult.SetBinContent(88, 12);
    mult.SetBinContent(89, 4);
    mult.SetBinContent(90, 10);
    mult.SetBinContent(91, 2);
    mult.SetBinContent(92, 6);
    mult.SetBinContent(93, 6);
    mult.SetBinContent(94, 3);
    mult.SetBinContent(96, 1);
    mult.SetBinContent(97, 2);
    mult.SetBinContent(99, 3);
    mult.SetBinContent(101, 1);
    mult.SetEntries(24625);
  }
  return mult;
}

//----------------------------------------------------------------------------//
// PID weight distribution for Parton sampling: {pid, weight}
// Approximate relative parton PID weights:
// Based on scaled input numbers
inline const std::vector<std::pair<int, double>>& GetPartonPIDWeights() {
  static const std::vector<std::pair<int, double>> pid_weights = {
      {-3, 3},  /* sbar: strange antiquark */
      {-2, 10}, /* dbar: down antiquark */
      {-1, 10}, /* ubar: up antiquark */
      {1, 10},  /* u: up quark */
      {2, 10},  /* d: down quark */
      {3, 3},   /* s: strange quark */
  };
  return pid_weights;
}

}  // namespace PhysicsConstants

#endif  // PHYSICS_CONSTANTS_H
