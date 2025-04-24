#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include "TH1.h"  // for TH1::kCanRebin
#include "ana/AnalyzerQA.h"
#include "TFile.h"
#include "TH1D.h"
#include <cmath>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

AnalyzerQA::AnalyzerQA()
  : hPt_b(nullptr), hPt_ab(nullptr), hPt_m(nullptr),
    hEta_b(nullptr), hEta_ab(nullptr), hEta_m(nullptr),
    hPhiM_b(nullptr), hPhiM_ab(nullptr), hPhiM_m(nullptr),
    hPhiP_b(nullptr), hPhiP_ab(nullptr), hPhiP_m(nullptr)
  , hPIDUnsort(nullptr), hPID(nullptr), hPIDName(nullptr)
  , hRatio(nullptr) {}

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

    // PID histogram with dynamic string labels
    hPIDUnsort = new TH1D("hPIDUnsort", "Unsorted Hadron PID labels;PID;Counts", 1, 0, 0);

    // Ratio histogram
    hRatio = new TH1D("hRatio", "Hadron Ratios;Ratio;Value", 7, 0.5, 7.5);
    {
      auto* ax = hRatio->GetXaxis();
      ax->SetBinLabel(1, "(#bar{B}+B)/M");
      ax->SetBinLabel(2, "#bar{B}/B");
      ax->SetBinLabel(3, "p/#pi^{+}");
      ax->SetBinLabel(4, "#bar{p}/p");
      ax->SetBinLabel(5, "#Lambda/p");
      ax->SetBinLabel(6, "K^{+}/#pi^{+}");
      ax->SetBinLabel(7, "#rho^{+}/#pi^{+}");
    }
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

    // Fill PID histogram
    for (auto* hadron : evt.GetHadrons()) {
        int pid = hadron->GetPID();
        std::string pidLabel = std::to_string(pid);
        hPIDUnsort->Fill(pidLabel.c_str(), 1.0);
    }
}

void AnalyzerQA::Finish(const std::string& outFileName) {
    TFile outFile(outFileName.c_str(), "RECREATE");
    hPt_b->Write();   hPt_ab->Write();   hPt_m->Write();
    hEta_b->Write();  hEta_ab->Write();  hEta_m->Write();
    hPhiM_b->Write(); hPhiM_ab->Write(); hPhiM_m->Write();
    hPhiP_b->Write(); hPhiP_ab->Write(); hPhiP_m->Write();

    // Prepare PID counts for sorting
    int nbins = hPIDUnsort->GetNbinsX();
    std::vector<std::pair<int,double>> pidCounts;
    pidCounts.reserve(nbins);
    for (int ib = 1; ib <= nbins; ++ib) {
        const char* lbl = hPIDUnsort->GetXaxis()->GetBinLabel(ib);
        if (!lbl || lbl[0]=='\0') continue;
        int pid = std::stoi(lbl);
        double cnt = hPIDUnsort->GetBinContent(ib);
        if (cnt <= 0) continue;
        pidCounts.emplace_back(pid, cnt);
    }
    // Sort by count descending
    std::sort(pidCounts.begin(), pidCounts.end(),
              [](auto &a, auto &b){ return a.second > b.second; });

    // 1) Histogram sorted by count with PID labels
    hPID = new TH1D("hPID", "PID Sorted by Count;PID;Counts",
                    pidCounts.size(), 0.5, pidCounts.size()+0.5);
    for (size_t i = 0; i < pidCounts.size(); ++i) {
        int pid = pidCounts[i].first;
        double cnt = pidCounts[i].second;
        hPID->SetBinContent(i+1, cnt);
        hPID->GetXaxis()->SetBinLabel(i+1, Form("%d", pid));
    }
    hPID->Write();

    // 2) Histogram sorted by count with particle names
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    hPIDName = new TH1D("hPIDName", "PID Sorted by Count with Names;Name;Counts",
                        pidCounts.size(), 0.5, pidCounts.size()+0.5);
    for (size_t i = 0; i < pidCounts.size(); ++i) {
        int pid = pidCounts[i].first;
        double cnt = pidCounts[i].second;
        TParticlePDG* part = pdgDB->GetParticle(pid);
        const char* name = part ? part->GetName() : "Unknown";
        hPIDName->SetBinContent(i+1, cnt);
        hPIDName->GetXaxis()->SetBinLabel(i+1, name);
    }
    hPIDName->Write();

    // Compute counts from hPIDUnsort using PDG
    double nBaryon = 0, nAntiBaryon = 0, nMeson = 0;
    double nProton = 0, nAntiProton = 0, nLambda = 0, nKaonPlus = 0, nRhoPlus = 0, nPionPlus = 0;
    for (int ib = 1; ib <= nbins; ++ib) {
      const char* lbl = hPIDUnsort->GetXaxis()->GetBinLabel(ib);
      if (!lbl || !*lbl) continue;
      int pid = std::stoi(lbl);
      double cnt = hPIDUnsort->GetBinContent(ib);
      if (cnt <= 0) continue;
      // Classify by PID using PDG ID range
      int abspid = std::abs(pid);
      if (abspid >= 1000 && abspid < 6000) {
          // Baryons have |PID| in [1000,6000)
          if (pid > 0) nBaryon += cnt;
          else         nAntiBaryon += cnt;
      } else {
          // Mesons have |PID| < 1000
          nMeson += cnt;
      }
      if      (pid == 2212)   nProton     = cnt;
      else if (pid == -2212)  nAntiProton = cnt;
      else if (pid == 3122)   nLambda     = cnt;
      else if (pid == 321)    nKaonPlus   = cnt;
      else if (pid == 213)    nRhoPlus    = cnt;
      else if (pid == 211)    nPionPlus   = cnt;
    }
    // Fill ratio bins
    if (nMeson > 0)      hRatio->SetBinContent(1, (nAntiBaryon + nBaryon) / nMeson);
    if (nBaryon > 0)     hRatio->SetBinContent(2, nAntiBaryon / nBaryon);
    if (nPionPlus > 0) {
      hRatio->SetBinContent(3, nProton     / nPionPlus);
      hRatio->SetBinContent(6, nKaonPlus   / nPionPlus);
      hRatio->SetBinContent(7, nRhoPlus    / nPionPlus);
    }
    if (nProton > 0)     hRatio->SetBinContent(4, nAntiProton / nProton);
    if (nProton > 0)     hRatio->SetBinContent(5, nLambda     / nProton);
    hRatio->Write();

    outFile.Close();
}
