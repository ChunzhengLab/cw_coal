#ifndef COMBINER_BASE_H
#define COMBINER_BASE_H

#include <vector>

#include "Particle.h"

class CombinerBase {
 public:
  virtual ~CombinerBase() = default;

  /// Combine partons into hadrons using a specific strategy
  virtual std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) = 0;

  /// Final sweep to combine all unused partons into mesons or baryons
  virtual std::vector<Hadron*> Afterburner(const std::vector<Parton*>& partons);

 protected:
  CombinerBase() = default;
};

#endif  // COMBINER_BASE_H
