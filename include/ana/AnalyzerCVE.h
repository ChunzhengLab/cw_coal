#ifndef ANALYZERCVE_H
#define ANALYZERCVE_H

#include <string>
#include "core/Event.h"

class TH1D;
class TProfile;

/// \brief Quality-assurance analyzer for hadron distributions.
/// Separately fills histograms for baryons, anti-baryons, and mesons.
class AnalyzerCVE {
public:
    /// Constructor: initializes pointers to nullptr.
    AnalyzerCVE();

    /// Destructor: deletes any allocated histograms.
    ~AnalyzerCVE();

    /// Create all histograms. Must be called before processing events.
    void Init();

    /// Process a single Event: fill histograms based on hadron baryon number.
    void Process(const Event& evt);

    /// Write histograms into a ROOT file and close it.
    /// \param outFileName name of the output ROOT file.
    void Finish(const std::string& outFileName);

private:
    //             baryon, antibaryon, meson
    // baryon,       0,       1,        2
    // antibaryon,            3,        4
    // meson,                           5

    TH1D* hCdPhiP[6]; // -> position space
    TH1D* hCdPhiM[6]; // -> momentum space

    TH1D* hSdPhiP[6]; // -> position space
    TH1D* hSdPhiM[6]; // -> momentum space

    TProfile* deltaP[6]; // delta = #LTcos(#Delta#phi)#GT (a !=b)
    TProfile* gammaP[6]; // gamma = #LTcos(#phi_{1}+#phi_{2})#GT (a !=b)

    TProfile* deltaS[6]; // delta = #LTcos(#Delta#phi)#GT (a !=b)
    TProfile* gammaS[6]; // gamma = #LTcos(#phi_{1}+#phi_{2})#GT (a !=b)

    TProfile* deltaMNtrk[6]; // delta = #LTcos(#Delta#phi)#GT (a !=b)
    TProfile* gammaMNtrk[6]; // gamma = #LTcos(#phi_{1}+#phi_{2})#GT (a !=b)

    TProfile* deltaSNtrk[6]; // delta = #LTcos(#Delta#phi)#GT (a !=b)
    TProfile* gammaSNtrk[6]; // gamma = #LTcos(#phi_{1}+#phi_{2})#GT (a !=b)
};

#endif // ANALYZERCVE_H
