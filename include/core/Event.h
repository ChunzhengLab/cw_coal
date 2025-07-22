#ifndef EVENT_H
#define EVENT_H

#include <vector>

#include "Particle.h"
#include "TObject.h"

class Event : public TObject {
 public:
  Event() : fEvtID(sNextID++) {
  }
  virtual ~Event() {
      Reset();
  }

  static void ResetIDCounter(unsigned int next = 1) {
    sNextID = next;
  }

  enum class ShuffleLevel : std::uint8_t {
    kLevel1,
    kLevel2,
    kLevel3,
    kLevel4
  };

  void ShufflePartons(ShuffleLevel level);
  /// Shuffle a fraction of parton positions; fraction must be between 0 and 1.
  void ShufflePartons(double fraction);

  void AddParton(Parton* p) {
    fPartons.push_back(p);
  }
  void AddHadron(Hadron* h) {
    fHadrons.push_back(h);
  }

  void SetReactionPlane(double psi) {
    fReactionPlane = psi;
  }
  void SetUID(unsigned int id) {
    fEvtID = id;
  }

  const std::vector<Parton*>& GetPartons() const {
    return fPartons;
  }
  const std::vector<Hadron*>& GetHadrons() const {
    return fHadrons;
  }
  double GetReactionPlane() const {
    return fReactionPlane;
  }
  unsigned int GetUID() const {
    return fEvtID;
  }

  int GetMultiplicity() const {
    return static_cast<int>(fHadrons.size());
  }

  void Reset(Option_t* opt = "");

 private:
  static unsigned int sNextID;
  unsigned int fEvtID{0};
  std::vector<Parton*> fPartons;
  std::vector<Hadron*> fHadrons;
  double fReactionPlane{0.0};

  ClassDef(Event, 4)
};

#endif  // EVENT_H
