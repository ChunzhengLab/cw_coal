#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <vector>
#include "PhysicsConstants.h"
#include "TPDGCode.h"
#include <TRandom3.h> // Added include for TRandom3

class Particle : public TObject {
public:
  Particle() : fUniqueID(sNextID++) {}
  Particle(double x, double y, double z,
           double px, double py, double pz, double baryonNumber)
    : fUniqueID(sNextID++), fX(x), fY(y), fZ(z),
      fPx(px), fPy(py), fPz(pz), fBaryonNumber(baryonNumber) {}

  virtual ~Particle() {}

  void SetMomentum(double px, double py, double pz) {
    fPx = px; fPy = py; fPz = pz;
  }

  void SetPosition(double x, double y, double z) {
    fX = x; fY = y; fZ = z;
  }

  double Px() const { return fPx; }
  double Py() const { return fPy; }
  double Pz() const { return fPz; }
  double X() const { return fX; }
  double Y() const { return fY; }
  double Z() const { return fZ; }

  void GetPosition(double* xyz) const {
    xyz[0] = fX; xyz[1] = fY; xyz[2] = fZ;
  }

  std::array<double, 3> GetPosition() const {
    return {fX, fY, fZ};
  }

  double DistanceTo(const Particle& other) const {
    double dx = fX - other.fX;
    double dy = fY - other.fY;
    double dz = fZ - other.fZ;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  double GetMassFromPDG() {
      auto m = PhysicsConstants::GetMass(fPID);
      return m.value_or(0.0);
  }

  void SetMass(double mass) {fMass = mass;}
  double GetMass() const {return fMass;}

  void SetPID(int pid) {fPID = pid;}
  int GetPID() const {return fPID;} 

  double GetBaryonNumber() const { return fBaryonNumber; }
  void SetBaryonNumber(double b) { fBaryonNumber = b; }
  
  unsigned int UniqueID() const { return fUniqueID; }
  void SetUID(unsigned int uid) { fUniqueID = uid; }

private:
  int fPID{0};
  double fPx{0}, fPy{0}, fPz{0};
  double fX{0}, fY{0}, fZ{0};
  double fMass{0};
  double fBaryonNumber{0.0};
  unsigned int fUniqueID{0}; // Changed to unsigned int
  static unsigned int sNextID; // Added static member

  ClassDef(Particle, 3)
};

// ================= Parton ==================

class Parton : public Particle {
public:
  /// Default constructor needed for ROOT I/O
  Parton() : Particle() {}
  Parton(double x, double y, double z,
         double px, double py, double pz, double baryonNumber)
    : Particle(x, y, z, px, py, pz, baryonNumber) {}

  void MarkUsed() { fUsed = true; }
  bool IsUsed() const { return fUsed; }

  // 可选：从模型或hist生成随机Parton
  static Parton* Random(TRandom3* rng = nullptr); // Modified declaration
  static Parton* RandomFromHists(const char* filename,
                                  TRandom3* rng = nullptr); // Modified declaration

private:
  bool fUsed{false};

  ClassDef(Parton, 2)
};

// ================= Hadron ==================

class Hadron : public Particle {
public:
  /// Default constructor needed for ROOT I/O
  Hadron() : Particle(), fFormationDistance(0.0) {}
  Hadron(double x, double y, double z,
         double px, double py, double pz, double baryonNumber, double formationDistance)
    : Particle(x, y, z, px, py, pz, baryonNumber), fFormationDistance(formationDistance) {}

  void SetFormationDistance(double d) { fFormationDistance = d; }
  double GetFormationDistance() const { return fFormationDistance; }

  void AddConstituentID(unsigned int id) { fConstituentIDs.push_back(id); }
  const std::vector<unsigned int>& GetConstituentIDs() const { return fConstituentIDs; }

private:
  double fFormationDistance{0.0};
  std::vector<unsigned int> fConstituentIDs; // Unique IDs of partons forming this hadron

  ClassDef(Hadron, 3)
};

#endif // PARTICLE_H