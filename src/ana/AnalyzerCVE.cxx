#include "ana/AnalyzerCVE.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVector2.h"
#include "core/Event.h"
#include "core/Particle.h"
#include <cmath>

// Labels for the six charge combinations
static const char* comboLabels[6] = {
    "Baryon-Baryon",
    "Baryon-AntiBaryon",
    "Baryon-Meson",
    "AntiBaryon-AntiBaryon",
    "AntiBaryon-Meson",
    "Meson-Meson"
};

static constexpr int idxMap[3][3] = {
    {0, 1, 2},
    {1, 3, 4},
    {2, 4, 5}
};

AnalyzerCVE::AnalyzerCVE() {
    for (int i = 0; i < 6; ++i) {
        hCdPhiP[i]     = nullptr;
        hCdPhiM[i]     = nullptr;
        deltaP[i]      = nullptr;
        gammaP[i]      = nullptr;
        deltaS[i]      = nullptr;
        gammaS[i]      = nullptr;
        deltaMNtrk[i]  = nullptr;
        gammaMNtrk[i]  = nullptr;
        deltaSNtrk[i]  = nullptr;
        gammaSNtrk[i]  = nullptr;
    }
}

AnalyzerCVE::~AnalyzerCVE() {
    // Histograms and profiles are owned by ROOT; do not delete here to avoid double-free.
}

void AnalyzerCVE::Init() {
    for (int i = 0; i < 6; ++i) {
        const std::string lab = comboLabels[i];
        // Δφ histograms
        hCdPhiP[i] = new TH1D(("hCdPhiP_" + lab).c_str(),
                              ("#Delta#phi Position " + lab).c_str(),
                              64, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hCdPhiP[i]->SetXTitle("#Delta#phi");
        hCdPhiP[i]->SetYTitle("Counts");

        hCdPhiM[i] = new TH1D(("hCdPhiM_" + lab).c_str(),
                              ("#Delta#phi Momentum " + lab).c_str(),
                              64, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hCdPhiM[i]->SetXTitle("#Delta#phi");
        hCdPhiM[i]->SetYTitle("Counts");

        // sumφ
        hSdPhiP[i] = new TH1D(("hSdPhiP_" + lab).c_str(),
                              ("Sum#phi Position " + lab).c_str(),
                              64, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hSdPhiP[i]->SetXTitle("Sum#phi");
        hSdPhiP[i]->SetYTitle("Counts");

        hSdPhiM[i] = new TH1D(("hSdPhiM_" + lab).c_str(),
                              ("Sum#phi Momentum " + lab).c_str(),
                              64, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hSdPhiM[i]->SetXTitle("Sum#phi");
        hSdPhiM[i]->SetYTitle("Counts");


        // delta and gamma
        deltaP[i] = new TProfile(("pDeltaP_" + lab).c_str(),
                                 ("#LTcos(#Delta#phi)#GT Position " + lab).c_str(),
                                 1, 0, 1);
        deltaP[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

        gammaP[i] = new TProfile(("pGammaP_" + lab).c_str(),
                                 ("#LTcos(#phi_{1}+#phi_{2})#GT Position " + lab).c_str(),
                                 1, 0, 1);
        gammaP[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");

        deltaS[i] = new TProfile(("pDeltaM_" + lab).c_str(),
                                 ("#LTcos(#Delta#phi)#GT Momentum " + lab).c_str(),
                                 1, 0, 1);
        deltaS[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

        gammaS[i] = new TProfile(("pGammaM_" + lab).c_str(),
                                 ("#LTcos(#phi_{1}+#phi_{2})#GT Momentum " + lab).c_str(),
                                 1, 0, 1);
        gammaS[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");

        // vs N_tracks profiles
        deltaMNtrk[i] = new TProfile(("pDeltaPNtrk_" + lab).c_str(),
                                     ("#LTcos(#Delta#phi)#GT vs N_{tracks} Position " + lab).c_str(),
                                     150, 0, 15000);
        deltaMNtrk[i]->SetXTitle("N_{tracks}");
        deltaMNtrk[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

        gammaMNtrk[i] = new TProfile(("pGammaPNtrk_" + lab).c_str(),
                                     ("#LTcos(#phi_{1}+#phi_{2})#GT vs N_{tracks} Position " + lab).c_str(),
                                     150, 0, 15000);
        gammaMNtrk[i]->SetXTitle("N_{tracks}");
        gammaMNtrk[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");

        deltaSNtrk[i] = new TProfile(("pDeltaMNtrk_" + lab).c_str(),
                                     ("#LTcos(#Delta#phi)#GT vs N_{tracks} Momentum " + lab).c_str(),
                                     150, 0, 15000);
        deltaSNtrk[i]->SetXTitle("N_{tracks}");
        deltaSNtrk[i]->SetYTitle("#LTcos(#Delta#phi)#GT");

        gammaSNtrk[i] = new TProfile(("pGammaMNtrk_" + lab).c_str(),
                                     ("#LTcos(#phi_{1}+#phi_{2})#GT vs N_{tracks} Momentum " + lab).c_str(),
                                     150, 0, 15000);
        gammaSNtrk[i]->SetXTitle("N_{tracks}");
        gammaSNtrk[i]->SetYTitle("#LTcos(#phi_{1}+#phi_{2})#GT");
    }
}

void AnalyzerCVE::Process(const Event& evt) {
    const auto& hadrons = evt.GetHadrons();
    const int nTrk = evt.GetMultiplicity();
    auto partonType = [](int b){ return b>0 ? 0 : (b<0 ? 1 : 2); }; // 0: baryon, 1: antibaryon, 2: meson
    for (size_t i = 0; i < hadrons.size(); ++i) {
        for (size_t j = i + 1; j < hadrons.size(); ++j) {
            Hadron* h1 = hadrons[i];
            Hadron* h2 = hadrons[j];
            auto range = [](double d) {
                d = TVector2::Phi_mpi_pi(d);
                if (d < -0.5*TMath::Pi()) d += 2*TMath::Pi();
                return d;
            };
            int b1 = h1->GetBaryonNumber();
            int b2 = h2->GetBaryonNumber();
            int idx = idxMap[partonType(b1)][partonType(b2)];
            if (idx < 0) continue;
            // kinematic cuts: 0.2 < pT < 3 GeV/c, -0.8 < eta < 0.8
            double pt1 = std::hypot(h1->Px(), h1->Py());
            double pt2 = std::hypot(h2->Px(), h2->Py());
            double eta1 = std::asinh(h1->Pz()/pt1);
            double eta2 = std::asinh(h2->Pz()/pt2);
            if (pt1 < 0.2 || pt1 > 3.0 || pt2 < 0.2 || pt2 > 3.0 ||
                eta1 < -0.8 || eta1 > 0.8 || eta2 < -0.8 || eta2 > 0.8) continue;

            // compute azimuthal angles in position and momentum space
            double phi_p1   = TMath::ATan2(h1->Y(),  h1->X());
            double phi_p2   = TMath::ATan2(h2->Y(),  h2->X());
            double dphi_p   = phi_p1 - phi_p2;
            double sumphi_p = phi_p1 + phi_p2;
            double delta_p = TMath::Cos(dphi_p);
            double gamma_p  = TMath::Cos(sumphi_p);

            double phi_m1   = TMath::ATan2(h1->Py(), h1->Px());
            double phi_m2   = TMath::ATan2(h2->Py(), h2->Px());
            double dphi_m   = phi_m1 - phi_m2;
            double sumphi_m = phi_m1 + phi_m2;
            double delta_m = TMath::Cos(dphi_m);
            double gamma_m  = TMath::Cos(sumphi_m);

            // Position-space
            hCdPhiP[idx]->Fill(range(dphi_p));
            hSdPhiP[idx]->Fill(range(sumphi_p));
            deltaP[idx]->Fill(0.5, delta_p);
            gammaP[idx]->Fill(0.5, gamma_p);
            deltaMNtrk[idx]->Fill(nTrk, delta_p);
            gammaMNtrk[idx]->Fill(nTrk, gamma_p);

            // Momentum-space
            hCdPhiM[idx]->Fill(range(dphi_m));
            hSdPhiM[idx]->Fill(range(sumphi_m));
            deltaS[idx]->Fill(0.5, delta_m);
            gammaS[idx]->Fill(0.5, gamma_m);
            deltaSNtrk[idx]->Fill(nTrk, delta_m);
            gammaSNtrk[idx]->Fill(nTrk, gamma_m);
        }
    }
}

void AnalyzerCVE::Finish(const std::string& outFileName) {
    TFile f(outFileName.c_str(), "RECREATE");
    for (int i = 0; i < 6; ++i) {
        if (hCdPhiP[i])     hCdPhiP[i]->Write();
        if (hCdPhiM[i])     hCdPhiM[i]->Write();
        if (hSdPhiP[i])     hSdPhiP[i]->Write();
        if (hSdPhiM[i])     hSdPhiM[i]->Write();
        if (deltaP[i])      deltaP[i]->Write();
        if (gammaP[i])      gammaP[i]->Write();
        if (deltaS[i])      deltaS[i]->Write();
        if (gammaS[i])      gammaS[i]->Write();
        if (deltaMNtrk[i])  deltaMNtrk[i]->Write();
        if (gammaMNtrk[i])  gammaMNtrk[i]->Write();
        if (deltaSNtrk[i])  deltaSNtrk[i]->Write();
        if (gammaSNtrk[i])  gammaSNtrk[i]->Write();
    }
    f.Close();
}