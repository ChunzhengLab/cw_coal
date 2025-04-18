#ifndef COMBINERS_H
#define COMBINERS_H

#include "core/CombinerBase.h"
#include "core/Particle.h"
#include <vector>

class BruteForceGlobal : public CombinerBase {
public:
    explicit BruteForceGlobal(double baryonPreference = 1.0) : m_r(baryonPreference) {}
    std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

private:
    double m_r;
};

class BruteForceGreedy : public CombinerBase {
public:
    explicit BruteForceGreedy(double baryonPreference = 1.0) : m_r(baryonPreference) {}
    std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

private:
    double m_r;
};

class BruteForceDualGreedy : public CombinerBase {
public:
    explicit BruteForceDualGreedy(double baryonPreference = 1.0) : m_r(baryonPreference) {}
    std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

private:
    double m_r;
};

class KDTreeGlobal : public CombinerBase {
public:
    explicit KDTreeGlobal(double baryonPreference = 1.0) : m_r(baryonPreference) {}
    std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

private:
    double m_r;
};

class KDTreeGreedy : public CombinerBase {
public:
    explicit KDTreeGreedy(double baryonPreference = 1.0) : m_r(baryonPreference) {}
    std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

private:
    double m_r;
};

class KDTreeDualGreedy : public CombinerBase {
public:
    explicit KDTreeDualGreedy(double baryonPreference = 1.0) : m_r(baryonPreference) {}
    std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) override;

private:
    double m_r;
};

#endif // COMBINERS_H
