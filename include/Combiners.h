#ifndef COMBINERS_H
#define COMBINERS_H

#include <vector>

#include "core/CombinerBase.h"
#include "core/Particle.h"

class ExhaustiveSorted : public CombinerBase {
 public:
  explicit ExhaustiveSorted(double baryonPreference = 1.0) : m_r(3 * baryonPreference) {
  }
  std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

 private:
  double m_r;
};

class ExhaustiveSequential : public CombinerBase {
 public:
  explicit ExhaustiveSequential(double baryonPreference = 1.0) : m_r(baryonPreference) {
  }
  std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

 private:
  double m_r;
};

class ExhaustiveCompetitive : public CombinerBase {
 public:
  explicit ExhaustiveCompetitive(double baryonPreference = 1.0) : m_r(3 * baryonPreference) {
  }
  std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

 private:
  double m_r;
};

class KDTreeSorted : public CombinerBase {
 public:
  explicit KDTreeSorted(double baryonPreference = 1.0) : m_r(3 * baryonPreference) {
  }
  std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

 private:
  double m_r;
};

class KDTreeSequential : public CombinerBase {
 public:
  explicit KDTreeSequential(double baryonPreference = 1.0) : m_r(baryonPreference) {
  }
  std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

 private:
  double m_r;
};

class KDTreeCompetitive : public CombinerBase {
 public:
  explicit KDTreeCompetitive(double baryonPreference = 1.0) : m_r(3 * baryonPreference) {
  }
  std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

 private:
  double m_r;
};

#endif  // COMBINERS_H
