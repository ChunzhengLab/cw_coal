#include "Particle.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <cmath>

ClassImp(Particle)
ClassImp(Parton)
ClassImp(Hadron)

unsigned int Particle::sNextID = 1;

Parton* Parton::Random(TRandom3* rng) {
  static TRandom3 defaultRng(0);
  if (!rng) rng = &defaultRng;

  double x = rng->Gaus(0, 0.5);
  double y = rng->Gaus(0, 0.5);
  double z = rng->Gaus(0, 2.0);

  double pt = rng->Exp(0.4);
  double phi = rng->Uniform(0, 2 * TMath::Pi());
  double px = pt * std::cos(phi);
  double py = pt * std::sin(phi);
  double pz = rng->Gaus(0, 0.6);

  double baryonNumber = rng->Rndm() < 0.5 ? +1.0/3.0 : -1.0/3.0;
  
  Parton* p = new Parton(x, y, z, px, py, pz, baryonNumber);
  return p;
}

Parton* Parton::RandomFromHists(const char* filename, TRandom3* rng) {
  // Use provided RNG or a static default
  static TRandom3 defaultRng(0);
  if (!rng) rng = &defaultRng;

  // Open the histogram file once
  static TFile* fin = nullptr;
  if (!fin) fin = TFile::Open(filename);
  if (!fin || fin->IsZombie()) {
    std::fprintf(stderr, "Error: cannot open histogram file %s\n", filename);
    return nullptr;
  }
  // Retrieve histograms
  auto hZ   = dynamic_cast<TH1D*>(fin->Get("h_z"));
  auto hXY  = dynamic_cast<TH2D*>(fin->Get("h_x_y"));
  auto hPxPy = dynamic_cast<TH2D*>(fin->Get("h_px_py"));
  auto hPz = dynamic_cast<TH1D*>(fin->Get("h_pz"));
  if (!hZ || !hXY || !hPxPy || !hPz) {
    std::fprintf(stderr, "Error: missing histograms in file %s\n", filename);
    return nullptr;
  }
  // Sample position
  double z = hZ->GetRandom();
  double x, y;
  hXY->GetRandom2(x, y);
  // Sample momentum
  double pz = hPz->GetRandom();
  double px, py;
  hPxPy->GetRandom2(px, py);
  // Assign baryon number randomly
  double baryonNumber = rng->Rndm() < 0.5 ? +1.0/3.0 : -1.0/3.0;
  // Create and return new Parton
  return new Parton(x, y, z, px, py, pz, baryonNumber);
}
