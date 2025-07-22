#include "ana/AnalyzerCVE.h"

#include <cmath>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TProfile.h"
#include "TVector2.h"
#include "core/Event.h"
#include "core/Particle.h"

// Labels for the six charge combinations
static const char* comboLabels[6] = {"Baryon_Baryon",         "Baryon_AntiBaryon", "Baryon_Meson",
                                     "AntiBaryon_AntiBaryon", "AntiBaryon_Meson",  "Meson_Meson"};

static constexpr int idxMap[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};

// map baryon number to type index: 0=baryon,1=antibaryon,2=meson
static auto partonType = [](double b) { return b > 0 ? 0 : (b < 0 ? 1 : 2); };
// normalize delta-phi into [-0.5π,1.5π)
static auto rangePhi = [](double d) {
  d = TVector2::Phi_mpi_pi(d);
  if (d < -0.5 * TMath::Pi()) d += 2 * TMath::Pi();
  return d;
};

AnalyzerCVE::AnalyzerCVE() {
  for (int i = 0; i < 6; ++i) {
    hCdPhiP[i] = nullptr;
    hCdPhiM[i] = nullptr;
    deltaP[i] = nullptr;
    gammaP[i] = nullptr;
    deltaS[i] = nullptr;
    gammaS[i] = nullptr;
    deltaMNtrk[i] = nullptr;
    gammaMNtrk[i] = nullptr;
    deltaSNtrk[i] = nullptr;
    gammaSNtrk[i] = nullptr;
  }
}

AnalyzerCVE::~AnalyzerCVE() {
  // Histograms and profiles are owned by ROOT; do not delete here to avoid double-free.
}

void AnalyzerCVE::Init() {
  for (int i = 0; i < 6; ++i) {
    const std::string suffix = isProcessMixed ? "_MixEvt" : "";
    const std::string lab = comboLabels[i] + suffix;
    const std::string titleSuffix = isProcessMixed ? " Mix Event" : "";
    // Δφ histograms
    hCdPhiP[i] = new TH1D(("hCdPhiP_" + lab).c_str(),
                          ("#Delta#phi Position " + std::string(comboLabels[i]) + titleSuffix).c_str(), 64,
                          -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    hCdPhiP[i]->SetXTitle("#Delta#phi");
    hCdPhiP[i]->SetYTitle("Counts");

    hCdPhiM[i] = new TH1D(("hCdPhiM_" + lab).c_str(),
                          ("#Delta#phi Momentum " + std::string(comboLabels[i]) + titleSuffix).c_str(), 64,
                          -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    hCdPhiM[i]->SetXTitle("#Delta#phi");
    hCdPhiM[i]->SetYTitle("Counts");

    // sumφ
    hSdPhiP[i] =
        new TH1D(("hSdPhiP_" + lab).c_str(), ("Sum#phi Position " + std::string(comboLabels[i]) + titleSuffix).c_str(),
                 64, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    hSdPhiP[i]->SetXTitle("Sum#phi");
    hSdPhiP[i]->SetYTitle("Counts");

    hSdPhiM[i] =
        new TH1D(("hSdPhiM_" + lab).c_str(), ("Sum#phi Momentum " + std::string(comboLabels[i]) + titleSuffix).c_str(),
                 64, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
    hSdPhiM[i]->SetXTitle("Sum#phi");
    hSdPhiM[i]->SetYTitle("Counts");

    // delta and gamma
    deltaP[i] =
        new TProfile(("pDeltaP_" + lab).c_str(),
                     ("#LTcos(#Delta#phi)#GT Position " + std::string(comboLabels[i]) + titleSuffix).c_str(), 1, 0, 1);
    deltaP[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

    gammaP[i] = new TProfile(
        ("pGammaP_" + lab).c_str(),
        ("#LTcos(#phi_{1}+#phi_{2})#GT Position " + std::string(comboLabels[i]) + titleSuffix).c_str(), 1, 0, 1);
    gammaP[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");

    deltaS[i] =
        new TProfile(("pDeltaM_" + lab).c_str(),
                     ("#LTcos(#Delta#phi)#GT Momentum " + std::string(comboLabels[i]) + titleSuffix).c_str(), 1, 0, 1);
    deltaS[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

    gammaS[i] = new TProfile(
        ("pGammaM_" + lab).c_str(),
        ("#LTcos(#phi_{1}+#phi_{2})#GT Momentum " + std::string(comboLabels[i]) + titleSuffix).c_str(), 1, 0, 1);
    gammaS[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");

    // vs N_tracks profiles
    deltaMNtrk[i] = new TProfile(
        ("pDeltaPNtrk_" + lab).c_str(),
        ("#LTcos(#Delta#phi)#GT vs N_{tracks} Position " + std::string(comboLabels[i]) + titleSuffix).c_str(), 150, 0,
        15000);
    deltaMNtrk[i]->SetXTitle("N_{tracks}");
    deltaMNtrk[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

    gammaMNtrk[i] = new TProfile(
        ("pGammaPNtrk_" + lab).c_str(),
        ("#LTcos(#phi_{1}+#phi_{2})#GT vs N_{tracks} Position " + std::string(comboLabels[i]) + titleSuffix).c_str(),
        150, 0, 15000);
    gammaMNtrk[i]->SetXTitle("N_{tracks}");
    gammaMNtrk[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");

    deltaSNtrk[i] = new TProfile(
        ("pDeltaMNtrk_" + lab).c_str(),
        ("#LTcos(#Delta#phi)#GT vs N_{tracks} Momentum " + std::string(comboLabels[i]) + titleSuffix).c_str(), 150, 0,
        15000);
    deltaSNtrk[i]->SetXTitle("N_{tracks}");
    deltaSNtrk[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

    gammaSNtrk[i] = new TProfile(
        ("pGammaMNtrk_" + lab).c_str(),
        ("#LTcos(#phi_{1}+#phi_{2})#GT vs N_{tracks} Momentum " + std::string(comboLabels[i]) + titleSuffix).c_str(),
        150, 0, 15000);
    gammaSNtrk[i]->SetXTitle("N_{tracks}");
    gammaSNtrk[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");
  }
}

void AnalyzerCVE::AnalyzePair(const Hadron* h1, const Hadron* h2, int nTrk) {
  // determine combination index
  int idx = idxMap[partonType(h1->GetBaryonNumber())][partonType(h2->GetBaryonNumber())];
  if (idx < 0) return;

  // kinematic cuts: 0.2 < pT < 3 GeV/c, -0.8 < eta < 0.8
  double pt1 = std::hypot(h1->Px(), h1->Py());
  double pt2 = std::hypot(h2->Px(), h2->Py());
  double eta1 = std::asinh(h1->Pz() / pt1);
  double eta2 = std::asinh(h2->Pz() / pt2);
  if (pt1 < 0.2 || pt1 > 8.0 || pt2 < 0.2 || pt2 > 8.0 || eta1 < -0.8 || eta1 > 0.8 || eta2 < -0.8 || eta2 > 0.8)
    return;

  // position-space angles
  double phi_p1 = TMath::ATan2(h1->Y(), h1->X());
  double phi_p2 = TMath::ATan2(h2->Y(), h2->X());
  double dphi_p = phi_p1 - phi_p2;
  double sumphi_p = phi_p1 + phi_p2;
  double delta_p = std::cos(dphi_p);
  double gamma_p = std::cos(sumphi_p);

  // momentum-space angles
  double phi_m1 = TMath::ATan2(h1->Py(), h1->Px());
  double phi_m2 = TMath::ATan2(h2->Py(), h2->Px());
  double dphi_m = phi_m1 - phi_m2;
  double sumphi_m = phi_m1 + phi_m2;
  double delta_m = std::cos(dphi_m);
  double gamma_m = std::cos(sumphi_m);

  // fill histograms
  // position
  hCdPhiP[idx]->Fill(rangePhi(dphi_p));
  hSdPhiP[idx]->Fill(rangePhi(sumphi_p));
  deltaP[idx]->Fill(0.5, delta_p);
  gammaP[idx]->Fill(0.5, gamma_p);
  deltaMNtrk[idx]->Fill(nTrk, delta_p);
  gammaMNtrk[idx]->Fill(nTrk, gamma_p);
  // momentum
  hCdPhiM[idx]->Fill(rangePhi(dphi_m));
  hSdPhiM[idx]->Fill(rangePhi(sumphi_m));
  deltaS[idx]->Fill(0.5, delta_m);
  gammaS[idx]->Fill(0.5, gamma_m);
  deltaSNtrk[idx]->Fill(nTrk, delta_m);
  gammaSNtrk[idx]->Fill(nTrk, gamma_m);
}

void AnalyzerCVE::Process(const Event& evt) {
  const auto& hadrons = evt.GetHadrons();
  const int nTrk = evt.GetMultiplicity();
  for (size_t i = 0; i < hadrons.size(); ++i) {
    for (size_t j = i + 1; j < hadrons.size(); ++j) {
      Hadron* h1 = hadrons[i];
      Hadron* h2 = hadrons[j];
      if (h1->IsAfterBurned()) continue;
      if (h2->IsAfterBurned()) continue;
      AnalyzePair(h1, h2, nTrk);
    }
  }
}

void AnalyzerCVE::ProcessMixed(const Event& signalEvt, const std::vector<const Event*>& backgroundEvts) {
  const auto& hadrons1 = signalEvt.GetHadrons();
  int nTrk = signalEvt.GetMultiplicity();
  for (const Event* bgEvt : backgroundEvts) {
    const auto& hadrons2 = bgEvt->GetHadrons();
    int nTrk2 = bgEvt->GetMultiplicity();
    for (size_t i = 0; i < hadrons1.size(); ++i) {
      for (size_t j = 0; j < hadrons2.size(); ++j) {
        if (hadrons1[i]->IsAfterBurned()) continue;
        if (hadrons2[j]->IsAfterBurned()) continue;
        AnalyzePair(hadrons1[i], hadrons2[j], nTrk + nTrk2);
      }
    }
  }
}

void AnalyzerCVE::Finish(const std::string& outFileName) {
  TFile f(outFileName.c_str(), "RECREATE");
  for (int i = 0; i < 6; ++i) {
    if (hCdPhiP[i]) hCdPhiP[i]->Write();
    if (hCdPhiM[i]) hCdPhiM[i]->Write();
    if (hSdPhiP[i]) hSdPhiP[i]->Write();
    if (hSdPhiM[i]) hSdPhiM[i]->Write();
    if (deltaP[i]) deltaP[i]->Write();
    if (gammaP[i]) gammaP[i]->Write();
    if (deltaS[i]) deltaS[i]->Write();
    if (gammaS[i]) gammaS[i]->Write();
    if (deltaMNtrk[i]) deltaMNtrk[i]->Write();
    if (gammaMNtrk[i]) gammaMNtrk[i]->Write();
    if (deltaSNtrk[i]) deltaSNtrk[i]->Write();
    if (gammaSNtrk[i]) gammaSNtrk[i]->Write();
  }
  f.Close();
}
