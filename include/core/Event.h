#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "Particle.h"
#include <vector>

class Event : public TObject {
public:
  Event() : fEvtID(sNextID++) {}
  virtual ~Event() {}

  static void ResetIDCounter(unsigned int next = 1) { sNextID = next; }

  void AddParton(Parton* p) { fPartons.push_back(p); }
  void AddHadron(Hadron* h) { fHadrons.push_back(h); }

  void SetReactionPlane(double psi) { fReactionPlane = psi; }
  void SetUniqueID(unsigned int id) { fEvtID = id; }

  const std::vector<Parton*>& GetPartons() const { return fPartons; }
  const std::vector<Hadron*>& GetHadrons() const { return fHadrons; }
  double GetReactionPlane() const { return fReactionPlane; }
  unsigned int GetUniqueID() const { return fEvtID; }

  int GetMultiplicity() const { return static_cast<int>(fHadrons.size()); }

  virtual void Clear(Option_t* opt = "") override;

private:
  static unsigned int sNextID;
  unsigned int fEvtID{0};
  std::vector<Parton*> fPartons;
  std::vector<Hadron*> fHadrons;
  double fReactionPlane{0.0};

  ClassDef(Event, 4)
};

#endif // EVENT_H
