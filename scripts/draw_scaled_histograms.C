void draw_scaled_histograms() {
    TFile* f = TFile::Open("cve_KDTreeGlobal_r1.50.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file!" << std::endl;
        return;
    }

    const char* names[] = {
        "hCdPhiM_Baryon-Baryon",
        "hCdPhiM_Baryon-AntiBaryon",
        "hCdPhiM_AntiBaryon-AntiBaryon"
    };

    const int colors[] = {kRed+1, kBlue+1, kGreen+2};

    TCanvas* c = new TCanvas("c", "Scaled Histograms", 800, 600);
    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    bool first = true;

    for (int i = 0; i < 3; ++i) {
        auto h = dynamic_cast<TH1*>(f->Get(names[i]));
        if (!h) {
            std::cerr << "Cannot find histogram: " << names[i] << std::endl;
            continue;
        }

        h->SetDirectory(nullptr);
        int nbins = h->GetNbinsX();
        double scale = h->GetEntries() > 0 ? nbins / h->GetEntries() : 1.0;
        h->Scale(scale);
        h->SetLineColor(colors[i]);
        h->SetLineWidth(2);

        if (first) {
            h->Draw("hist");
            first = false;
        } else {
            h->Draw("hist same");
        }

        leg->AddEntry(h, names[i], "l");
    }

    leg->Draw();
    c->SaveAs("scaled_histos.pdf");
}
