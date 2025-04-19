#include <vector>
#include <algorithm>
#include <cmath>
#include "TLatex.h"

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TH1D.h"
#include "TLegend.h"

#include "../include/core/Event.h"
#include "../include/core/Particle.h"

void visualize_event(const char* filename = "test_KDTreeGlobal.root", int eventIndex = 0) {

    gSystem->Load("libcw_coal.dylib");

    // Open file
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        ::fprintf(stderr, "Error: cannot open file %s\n", filename);
        return;
    }

    // Read tree via SetBranchAddress
    TTree* tree = nullptr;
    file->GetObject("events", tree);
    if (!tree) {
        ::fprintf(stderr, "Error: TTree 'events' not found in %s\n", filename);
        return;
    }
    Event* evt = nullptr;
    tree->SetBranchAddress("event", &evt);
    tree->GetEntry(eventIndex);

    // Extract collections
    const auto& partons = evt->GetPartons();
    const auto& hadrons = evt->GetHadrons();

    // Gather coordinates
    std::vector<double> xs, ys, zs;
    xs.reserve(partons.size()+hadrons.size());
    ys.reserve(partons.size()+hadrons.size());
    zs.reserve(partons.size()+hadrons.size());
    for (auto* p : partons) { xs.push_back(p->X()); ys.push_back(p->Y()); zs.push_back(p->Z()); }
    for (auto* h : hadrons) { xs.push_back(h->X()); ys.push_back(h->Y()); zs.push_back(h->Z()); }

    auto [xminIt, xmaxIt] = std::minmax_element(xs.begin(), xs.end());
    auto [yminIt, ymaxIt] = std::minmax_element(ys.begin(), ys.end());
    auto [zminIt, zmaxIt] = std::minmax_element(zs.begin(), zs.end());
    double xAbs = 1.5 * std::max(std::abs(*xminIt), std::abs(*xmaxIt));
    double yAbs = 1.5 * std::max(std::abs(*yminIt), std::abs(*ymaxIt));
    double zAbs = 1.5 * std::max(std::abs(*zminIt), std::abs(*zmaxIt));

    // Common color codes
    constexpr int colQ  = kRed, colAQ = kBlue;
    constexpr int colB  = kGreen+2, colAB = kMagenta, colM = kBlack;

    // --- Canvas 1: 2D scatter + composition ---
    TCanvas* c2d = new TCanvas("c2d", "2D Scatter", 1200, 600);
    c2d->Divide(2,1);

    auto drawProjection = [&](int padIndex, 
                              auto coordX, auto coordY) {
        c2d->cd(padIndex);
        // draw axes frame with symmetric limits
        double x1, x2, y1, y2;
        if (padIndex == 1) {
            x1 = -xAbs; x2 = xAbs; y1 = -yAbs; y2 = yAbs;
        } else {
            x1 = -zAbs; x2 = zAbs; y1 = -yAbs; y2 = yAbs;
        }
        if (padIndex == 1) {
            // XY projection: label axes X and Y
            gPad->DrawFrame(x1, y1, x2, y2, ";X;Y");
        } else {
            // ZY projection: label axes Z and Y
            gPad->DrawFrame(x1, y1, x2, y2, ";Z;Y");
        }

        // Parton graphs
        TGraph* gQ  = new TGraph();
        TGraph* gAQ = new TGraph();
        for (auto* p : partons) {
            double x = coordX(p), y = coordY(p);
            if (p->GetBaryonNumber() > 0) gQ->AddPoint(x, y);
            else                          gAQ->AddPoint(x, y);
        }
        gQ->SetMarkerStyle(20);  gQ->SetMarkerColor(colQ);
        gAQ->SetMarkerStyle(21); gAQ->SetMarkerColor(colAQ);
        
        // Hadron graphs
        TGraph* gB  = new TGraph();
        TGraph* gAB = new TGraph();
        TGraph* gM  = new TGraph();
        for (auto* h : hadrons) {
            double x = coordX(h), y = coordY(h);
            if      (h->GetBaryonNumber() > 0) gB->AddPoint(x, y);
            else if (h->GetBaryonNumber() < 0) gAB->AddPoint(x, y);
            else                               gM->AddPoint(x, y);
        }
        gB->SetMarkerStyle(22);  gB->SetMarkerColor(colB);
        gAB->SetMarkerStyle(23); gAB->SetMarkerColor(colAB);
        gM->SetMarkerStyle(33);  gM->SetMarkerColor(colM);

        // draw parton graphs
        gQ->Draw("P same");
        gAQ->Draw("P same");
        // draw hadron graphs
        gB->Draw("P same");
        gAB->Draw("P same");
        gM->Draw("P same");

        // Composition: connect partons for each hadron
        for (auto* h : hadrons) {
            std::vector<Parton*> comps;
            for (auto id : h->GetConstituentIDs()) {
                auto it = std::find_if(partons.begin(), partons.end(),
                    [&](Parton* p){ return p->UniqueID() == id; });
                if (it != partons.end()) comps.push_back(*it);
            }
            if (comps.size() == 3) {
                TPolyLine* poly = new TPolyLine(4);
                for (int j=0; j<3; ++j) poly->SetPoint(j, coordX(comps[j]), coordY(comps[j]));
                poly->SetPoint(3, coordX(comps[0]), coordY(comps[0]));
                poly->SetLineColor(colB);
                poly->Draw("same L");
            }
            else if (comps.size() == 2) {
                TLine* line = new TLine(
                    coordX(comps[0]), coordY(comps[0]),
                    coordX(comps[1]), coordY(comps[1]) );
                line->SetLineColor(colM);
                line->Draw("same");
            }
        }
    };

    // XY projection: X->X, Y->Y
    drawProjection(1, 
                   [](auto* p){ return p->X(); },
                   [](auto* p){ return p->Y(); });
    // ZY projection: X->Z, Y->Y
    drawProjection(2,
                   [](auto* p){ return p->Z(); },
                   [](auto* p){ return p->Y(); });

    // --- Canvas 2: 1D distributions ---
    constexpr int bins1d = 100;
    // 1D distributions separated by particle type
    TH1D* hX_q   = new TH1D("hX_q","X Distribution - Quarks;X;Counts", bins1d, -xAbs, xAbs);
    TH1D* hX_aq  = new TH1D("hX_aq",";X;Counts", bins1d, -xAbs, xAbs);
    TH1D* hX_b   = new TH1D("hX_b",";X;Counts", bins1d, -xAbs, xAbs);
    TH1D* hX_ab  = new TH1D("hX_ab",";X;Counts", bins1d, -xAbs, xAbs);
    TH1D* hX_m   = new TH1D("hX_m",";X;Counts", bins1d, -xAbs, xAbs);

    TH1D* hY_q   = new TH1D("hY_q","Y Distribution - Quarks;Y;Counts", bins1d, -yAbs, yAbs);
    TH1D* hY_aq  = new TH1D("hY_aq",";Y;Counts", bins1d, -yAbs, yAbs);
    TH1D* hY_b   = new TH1D("hY_b",";Y;Counts", bins1d, -yAbs, yAbs);
    TH1D* hY_ab  = new TH1D("hY_ab",";Y;Counts", bins1d, -yAbs, yAbs);
    TH1D* hY_m   = new TH1D("hY_m",";Y;Counts", bins1d, -yAbs, yAbs);

    TH1D* hZ_q   = new TH1D("hZ_q","Z Distribution - Quarks;Z;Counts", bins1d, -zAbs, zAbs);
    TH1D* hZ_aq  = new TH1D("hZ_aq",";Z;Counts", bins1d, -zAbs, zAbs);
    TH1D* hZ_b   = new TH1D("hZ_b",";Z;Counts", bins1d, -zAbs, zAbs);
    TH1D* hZ_ab  = new TH1D("hZ_ab",";Z;Counts", bins1d, -zAbs, zAbs);
    TH1D* hZ_m   = new TH1D("hZ_m",";Z;Counts", bins1d, -zAbs, zAbs);

    // Fill histograms
    for (auto* p : partons) {
        if (p->GetBaryonNumber() > 0) {
            hX_q->Fill(p->X()); hY_q->Fill(p->Y()); hZ_q->Fill(p->Z());
        } else {
            hX_aq->Fill(p->X()); hY_aq->Fill(p->Y()); hZ_aq->Fill(p->Z());
        }
    }
    for (auto* h : hadrons) {
        if      (h->GetBaryonNumber() > 0) {
            hX_b->Fill(h->X()); hY_b->Fill(h->Y()); hZ_b->Fill(h->Z());
        } else if (h->GetBaryonNumber() < 0) {
            hX_ab->Fill(h->X()); hY_ab->Fill(h->Y()); hZ_ab->Fill(h->Z());
        } else {
            hX_m->Fill(h->X()); hY_m->Fill(h->Y()); hZ_m->Fill(h->Z());
        }
    }

    // Draw on canvas
    TCanvas* c1d = new TCanvas("c1d","1D Distributions", 1200, 800);
    c1d->Divide(2,2);

    // X-axis distributions
    c1d->cd(1);
    hX_q->SetLineColor(colQ);    hX_q->Draw("same");
    hX_aq->SetLineColor(colAQ);  hX_aq->Draw("SAME");
    hX_b->SetLineColor(colB);    hX_b->Draw("SAME");
    hX_ab->SetLineColor(colAB);  hX_ab->Draw("SAME");
    hX_m->SetLineColor(colM);    hX_m->Draw("SAME");
    gPad->Update();

    // Y-axis distributions
    c1d->cd(2);
    hY_q->SetLineColor(colQ);    hY_q->Draw("same");
    hY_aq->SetLineColor(colAQ);  hY_aq->Draw("SAME");
    hY_b->SetLineColor(colB);    hY_b->Draw("SAME");
    hY_ab->SetLineColor(colAB);  hY_ab->Draw("SAME");
    hY_m->SetLineColor(colM);    hY_m->Draw("SAME");
    gPad->Update();

    // Z-axis distributions
    c1d->cd(3);
    hZ_q->SetLineColor(colQ);    hZ_q->Draw("same");
    hZ_aq->SetLineColor(colAQ);  hZ_aq->Draw("SAME");
    hZ_b->SetLineColor(colB);    hZ_b->Draw("SAME");
    hZ_ab->SetLineColor(colAB);  hZ_ab->Draw("SAME");
    hZ_m->SetLineColor(colM);    hZ_m->Draw("SAME");
    gPad->Update();
    // Pad 4: Legend
    c1d->cd(4);
    TLegend* legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->AddEntry(hX_q, "Parton (B>0)", "l");
    legend->AddEntry(hX_aq, "Parton (B<0)", "l");
    legend->AddEntry(hX_b,  "Hadron (B>0)", "l");
    legend->AddEntry(hX_ab, "Hadron (B<0)", "l");
    legend->AddEntry(hX_m,  "Meson (B=0)",   "l");
    legend->Draw();
    gPad->Update();

    // --- Canvas 4: Global Ratios ---
    TCanvas* cRatio = new TCanvas("cRatio", "Global Ratios", 600, 400);
    // Count hadron types
    int countB = 0, countAB = 0, countM = 0;
    for (auto* h : hadrons) {
        double bn = h->GetBaryonNumber();
        if      (bn > 0)  ++countB;
        else if (bn < 0)  ++countAB;
        else              ++countM;
    }
    // Calculate ratios
    double ratioBM = (countM > 0 ? double(countB) / countM : 0);
    double ratioBAB = (countAB > 0 ? double(countB) / countAB : 0);
    // Create a 2-bin histogram: bin1=B/M, bin2=B/#bar{B}
    TH1D* hRatio = new TH1D("hRatio","Global Ratios;Type;Value",2,0.5,2.5);
    hRatio->GetXaxis()->SetBinLabel(1,"B/M");
    hRatio->GetXaxis()->SetBinLabel(2,"B/#bar{B}");
    hRatio->SetBinContent(1, ratioBM);
    hRatio->SetBinContent(2, ratioBAB);
    hRatio->SetLineColor(kBlack);
    hRatio->SetFillColorAlpha(kGray, 0.5);
    hRatio->Draw("B");
    // Draw bar values
    for (int ib=1; ib<=2; ++ib) {
        double x = hRatio->GetBinCenter(ib);
        double y = hRatio->GetBinContent(ib);
        TLatex txt;
        txt.SetTextAlign(22);
        txt.DrawLatex(x, y + 0.03*y, Form("%.2f", y));
    }
    cRatio->Update();
}
