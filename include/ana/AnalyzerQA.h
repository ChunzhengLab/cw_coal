#ifndef ANALYZERQA_H
#define ANALYZERQA_H

#include <string>
#include "core/Event.h"

class TH1D;

/// \brief Quality-assurance analyzer for hadron distributions.
/// Separately fills histograms for baryons, anti-baryons, and mesons.
class AnalyzerQA {
public:
    /// Constructor: initializes pointers to nullptr.
    AnalyzerQA();

    /// Destructor: deletes any allocated histograms.
    ~AnalyzerQA();

    /// Create all histograms. Must be called before processing events.
    void Init();

    /// Process a single Event: fill histograms based on hadron baryon number.
    void Process(const Event& evt);

    /// Write histograms into a ROOT file and close it.
    /// \param outFileName name of the output ROOT file.
    void Finish(const std::string& outFileName);

private:
    // Histograms for pT distributions
    TH1D* hPt_b;
    TH1D* hPt_ab;
    TH1D* hPt_m;

    // Histograms for eta distributions
    TH1D* hEta_b;
    TH1D* hEta_ab;
    TH1D* hEta_m;

    // Histograms for momentum-space phi
    TH1D* hPhiM_b;
    TH1D* hPhiM_ab;
    TH1D* hPhiM_m;

    // Histograms for position-space phi
    TH1D* hPhiP_b;
    TH1D* hPhiP_ab;
    TH1D* hPhiP_m;
};

#endif // ANALYZERQA_H
