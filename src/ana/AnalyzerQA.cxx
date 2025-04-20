#include "ana/AnalyzerQA.h"
#include "TFile.h"
#include "TH1D.h"
#include <cmath>

AnalyzerQA::AnalyzerQA()
  : hPt_b(nullptr), hPt_ab(nullptr), hPt_m(nullptr),
    hEta_b(nullptr), hEta_ab(nullptr), hEta_m(nullptr),
    hPhiM_b(nullptr), hPhiM_ab(nullptr), hPhiM_m(nullptr),
    hPhiP_b(nullptr), hPhiP_ab(nullptr), hPhiP_m(nullptr) {}

AnalyzerQA::~AnalyzerQA() {
    // Histograms and profiles are owned by ROOT; do not delete here to avoid double-free.
    // delete hPt_b;   delete hPt_ab;   delete hPt_m;
    // delete hEta_b;  delete hEta_ab;  delete hEta_m;
    // delete hPhiM_b; delete hPhiM_ab; delete hPhiM_m;
    // delete hPhiP_b; delete hPhiP_ab; delete hPhiP_m;
}

void AnalyzerQA::Init() {
    // pT histograms
    hPt_b  = new TH1D("hPt_b",  "p_{T} - Baryons; p_{T}; Counts",    100, 0, 10);
    hPt_ab = new TH1D("hPt_ab", "p_{T} - Anti-Baryons; p_{T}; Counts",100, 0, 10);
    hPt_m  = new TH1D("hPt_m",  "p_{T} - Mesons; p_{T}; Counts",       100, 0, 10);
    // eta histograms
    hEta_b  = new TH1D("hEta_b",  "#eta - Baryons; #eta; Counts",   100, -5, 5);
    hEta_ab = new TH1D("hEta_ab", "#eta - Anti-Baryons; #eta; Counts",100, -5, 5);
    hEta_m  = new TH1D("hEta_m",  "#eta - Mesons; #eta; Counts",    100, -5, 5);
    // momentum-space phi
    hPhiM_b  = new TH1D("hPhiM_b",  "#phi_{m} - Baryons; #phi_{m}; Counts",64, 0, 2*M_PI);
    hPhiM_ab = new TH1D("hPhiM_ab", "#phi_{m} - Anti-Baryons; #phi_{m}; Counts",64, 0, 2*M_PI);
    hPhiM_m  = new TH1D("hPhiM_m",  "#phi_{m} - Mesons; #phi_{m}; Counts",64, 0, 2*M_PI);
    // position-space phi
    hPhiP_b  = new TH1D("hPhiP_b",  "#phi_{p} - Baryons; #phi_{p}; Counts",64, 0, 2*M_PI);
    hPhiP_ab = new TH1D("hPhiP_ab", "#phi_{p} - Anti-Baryons; #phi_{p}; Counts",64, 0, 2*M_PI);
    hPhiP_m  = new TH1D("hPhiP_m",  "#phi_{p} - Mesons; #phi_{p}; Counts",64, 0, 2*M_PI);
}

void AnalyzerQA::Process(const Event& evt) {
    for (auto* h : evt.GetHadrons()) {
        double px = h->Px(), py = h->Py(), pz = h->Pz();
        double x  = h->X(),  y  = h->Y();
        double pt = std::hypot(px, py);
        double eta = std::asinh(pz/pt);
        double phi_m = std::atan2(py, px); if (phi_m < 0) phi_m += 2*M_PI;
        double phi_p = std::atan2(y, x);   if (phi_p < 0) phi_p += 2*M_PI;

        double bn = h->GetBaryonNumber();
        if      (bn > 0) {
            hPt_b->Fill(pt);  hEta_b->Fill(eta);
            hPhiM_b->Fill(phi_m); hPhiP_b->Fill(phi_p);
        } else if (bn < 0) {
            hPt_ab->Fill(pt); hEta_ab->Fill(eta);
            hPhiM_ab->Fill(phi_m); hPhiP_ab->Fill(phi_p);
        } else {
            hPt_m->Fill(pt);  hEta_m->Fill(eta);
            hPhiM_m->Fill(phi_m); hPhiP_m->Fill(phi_p);
        }
    }
}

void AnalyzerQA::Finish(const std::string& outFileName) {
    TFile outFile(outFileName.c_str(), "RECREATE");
    hPt_b->Write();   hPt_ab->Write();   hPt_m->Write();
    hEta_b->Write();  hEta_ab->Write();  hEta_m->Write();
    hPhiM_b->Write(); hPhiM_ab->Write(); hPhiM_m->Write();
    hPhiP_b->Write(); hPhiP_ab->Write(); hPhiP_m->Write();
    outFile.Close();
}
