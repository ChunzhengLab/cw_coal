#include "Particle.h"

#include <cmath>
#include <utility>

#include "PhysicsConstants.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"

ClassImp(Particle);
ClassImp(Parton);
ClassImp(Hadron);

    unsigned int Particle::sNextID = 1;

Parton* Parton::Random(TRandom3* rng) {
  static TRandom3 defaultRng(0);
  if (!rng) rng = &defaultRng;

  // Position: uniform in unit disk
  double u = rng->Rndm();  // uniform [0,1)
  double phi_pos = rng->Uniform(0, 2 * TMath::Pi());
  double rpos = std::sqrt(u);  // ensures uniform area
  double x = rpos * std::cos(phi_pos);
  double y = rpos * std::sin(phi_pos);
  double z = 0.0;

  // pT: Tsallis-like distribution
  const double T = 0.7;  // "temperature" parameter in GeV
  const double n = 4.0;  // power index (controls high-pT tail)
  double pT;
  do {
    pT = rng->Uniform(0, 5.0);              // Max pT up to 5 GeV
    double f = std::pow(1.0 + pT / T, -n);  // Tsallis shape
    double accept = rng->Uniform();
    if (accept < f) break;  // accept-reject method
  } while (true);

  // Momentum angle
  double phi = rng->Uniform(0, 2 * TMath::Pi());
  double px = pT * std::cos(phi);
  double py = pT * std::sin(phi);
  double pz = 0.0;  // no longitudinal momentum

  // Baryon number: ±1/3 randomly
  double baryonNumber = (rng->Rndm() < 0.5) ? +1.0 / 3.0 : -1.0 / 3.0;

  // Create the parton
  Parton* p = new Parton(x, y, z, px, py, pz, baryonNumber);
  // Sample PID based on predefined weighted distribution
  const auto& pid_weights = PhysicsConstants::GetPartonPIDWeights();
  double total_weight = 0;
  for (const auto& pw : pid_weights) total_weight += pw.second;
  double r1 = rng->Rndm() * total_weight;
  for (const auto& pw : pid_weights) {
    if (r1 < pw.second) {
      p->SetPID(pw.first);
      break;
    }
    r1 -= pw.second;
  }
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
  auto hZ = dynamic_cast<TH1D*>(fin->Get("h_z"));
  auto hXY = dynamic_cast<TH2D*>(fin->Get("h_x_y"));
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
  double baryonNumber = rng->Rndm() < 0.5 ? +1.0 / 3.0 : -1.0 / 3.0;
  // Create the parton
  Parton* p = new Parton(x, y, z, px, py, pz, baryonNumber);
  // Sample PID based on predefined weighted distribution
  const auto& pid_weights = PhysicsConstants::GetPartonPIDWeights();
  double total_weight = 0;
  for (const auto& pw : pid_weights) total_weight += pw.second;
  double r1 = rng->Rndm() * total_weight;
  for (const auto& pw : pid_weights) {
    if (r1 < pw.second) {
      p->SetPID(pw.first);
      break;
    }
    r1 -= pw.second;
  }
  return p;
}
