#ifndef ANALYZERQA_H
#define ANALYZERQA_H

#include <string>

#include "core/Event.h"

class TH1D;
class TProfile;

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
  TH1D* hMult_a;
  TH1D* hMult_ab;
  TH1D* hMult_m;

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

  TProfile* v2_pt_b;
  TProfile* v2_pt_ab;
  TProfile* v2_pt_m;

  TH1D* hPIDUnsort;
  TH1D* hPID;
  TH1D* hPIDName;

  TH1D* hRatio;

  TProfile* hAfterBurnedFlagRatio;

  // Global hadron counters, updated in Process()
  double m_nBaryonCount;
  double m_nAntiBaryonCount;
  double m_nMesonCount;
  double m_nProtonCount;
  double m_nAntiProtonCount;
  double m_nLambdaCount;
  double m_nKaonPlusCount;
  double m_nRhoPlusCount;
  double m_nPionPlusCount;
};

#endif  // ANALYZERQA_H
