#include "core/PIDInference.h"
#include <algorithm>
#include <iostream>
#include "core/PhysicsConstants.h"
#include <random>
#include <vector>


// 定义内部 RNG
std::mt19937 PIDInference::rng_{std::random_device{}()};
std::uniform_real_distribution<double> PIDInference::dist_{0.0, 1.0};

// 根据两个夸克/反夸克的 flavor code 和系统质量，推断介子 PDG 码
// 对角态返回 π0 (111)、φ (333) 或调用 InferQuarkoniumPDG
// 非对角态按 100*kmax+10*kmin+(2S+1) 构造并加符号与相位
int PIDInference::InferMesonPDG(int q1, int q2, double mass)
{
    // 先处理所有对角态 quark–antiquark
    if (q1 == -q2) {
        int af = std::abs(q1);
        if (af <= 2) {
            return PIDInference::kPidPi0;   // π0 for uū/dd̄
        } else if (af == 3) {
            return PIDInference::kPidPhi;   // φ for s s̄
        } else {
            return InferQuarkoniumPDG(af, mass);  // heavy quarkonium
        }
    }
    // 非对角态按 100*kmax+10*kmin+(2S+1) 构造并加符号与相位
    int qmax = std::max(std::abs(q1), std::abs(q2));
    int qmin = std::min(std::abs(q1), std::abs(q2));
    int pdg = 100 * qmax + 10 * qmin + 1;
    int sign = ((q1 + q2) > 0) ? 1 : -1;
    int phase = ((qmax % 2 == 0) ? 1 : -1);
    return pdg * sign * phase;
}

// 根据三个夸克的 flavor code、系统质量 和自旋多重度，推断重子/反重子 PDG 码
// 降序排序 quarks，octet(2) 与 decuplet(4) 自旋多重度，
// octet 情况下按质量接近度择优两种排列，负总荷判定反重子
int PIDInference::InferBaryonPDG(int q1, int q2, int q3, double mass, int spinMult) {
    // Sort quark flavors in descending order for PDG construction
    std::vector<int> quarks = { std::abs(q1), std::abs(q2), std::abs(q3) };
    std::sort(quarks.begin(), quarks.end(), std::greater<int>());
    int k1 = quarks[0], k2 = quarks[1], k3 = quarks[2];

    // uds: choose Λ0 vs Σ0 by mass closeness
    if (k1 == 3 && k2 == 2 && k3 == 1) {
        double mLambda = PhysicsConstants::GetMass(3122).value_or(0.0);
        double mSigma0 = PhysicsConstants::GetMass(3212).value_or(0.0);
        bool wantSigma = std::abs(mass - mSigma0) < std::abs(mass - mLambda);
        int pdg = wantSigma ? 3212 : 3122;
        return (q1 + q2 + q3 > 0 ? pdg : -pdg);
    }

    // Determine spin multiplicity (2S+1): decuplet (4) vs octet (2)
    int mult = spinMult;
    if (spinMult < 0) {
        // auto-detect: all identical → decuplet (spin 3/2 → 4), else octet (spin 1/2 → 2)
        mult = (k1 == k2 && k2 == k3) ? 4 : 2;
    }

    int pdg = 0;
    if (mult == 4) {
        // decuplet: only one ordering
        pdg = 1000 * k1 + 100 * k2 + 10 * k3 + mult;
    } else {
        // octet: two possible orderings, pick by mass closeness
        int pdg1 = 1000 * k1 + 100 * k2 + 10 * k3 + mult;
        int pdg2 = 1000 * k1 + 100 * k3 + 10 * k2 + mult;
        double m1 = PhysicsConstants::GetMass(pdg1).value_or(0.0);
        double m2 = PhysicsConstants::GetMass(pdg2).value_or(0.0);
        // 当与两种排列的差值相等或更小时，优先选择 pdg1
        pdg = (std::abs(mass - m1) <= std::abs(mass - m2)) ? pdg1 : pdg2;
    }

    // If it's an antibaryon (sum of quark charges negative), flip sign
    if ((q1 + q2 + q3) < 0) {
        pdg = -pdg;
    }
    return pdg;
}

// 根据重夸克 flavor (4=c,5=b,6=t) 与质量，在预定义偶素列表中选最接近的 PDG 码
int PIDInference::InferQuarkoniumPDG(int flavor, double mass) {
    if (flavor == 4) { // charm
        // Local charmonium PDG codes
        static const std::vector<int> pdgList = {
            441, 443, 10441, 20443, 10443, 445,
            100441, 100443, 30443, 100445, 9000443, 9010443, 9020443
        };
        int bestPdg = pdgList.front();
        double bestDiff = std::abs(mass - PhysicsConstants::GetMass(bestPdg).value_or(0.0));
        for (int pdg : pdgList) {
            double m = PhysicsConstants::GetMass(pdg).value_or(0.0);
            double diff = std::abs(mass - m);
            if (diff < bestDiff) {
                bestDiff = diff;
                bestPdg = pdg;
            }
        }
        return bestPdg;
    } else if (flavor == 5) { // bottom
        int pdg = 100 * flavor + 10 * flavor + 1;
        if (std::abs(mass - PhysicsConstants::GetMass(pdg).value_or(0.0)) > std::abs(mass - PhysicsConstants::GetMass(pdg + 2).value_or(0.0))) {
            pdg += 2;
            if (std::abs(mass - PhysicsConstants::GetMass(553).value_or(0.0)) > std::abs(mass - PhysicsConstants::GetMass(10551).value_or(0.0)))
                pdg = 10551;
        }
        return pdg;
    } else if (flavor == 6) { // top
        return (std::abs(mass - PhysicsConstants::GetMass(661).value_or(0.0)) > std::abs(mass - PhysicsConstants::GetMass(663).value_or(0.0))) ? 663 : 661;
    }

    return 0;
}

// 基于矢量/伪标量产额比 (vpratio) 与随机数 rnd，判断介子自旋：0 (伪标量) 或 1 (向量)
int PIDInference::InferMesonSpin(double vpratio, double rnd) {
    double pPseudo = 1.0 / (1.0 + vpratio);
    return (rnd < pPseudo) ? 0 : 1;
}

// 针对轻味对角态 (uū、dd̄)，根据 rrhopi、romrho0 计算 π⁰/η/ρ⁰/ω 概率，并用 rnd 随机分配
int PIDInference::ResolveDiagonalLightMeson(double rnd, double rrhopi, double romrho0) {
    // Based on AMPT logic: probabilities for π0, η, ρ0, ω 
    // p_pi0 = 1/(2*(1+rrhopi))
    // p_rho0 = rrhopi/(2*(1+rrhopi))
    // p_omega = rrhopi*romrho0/(2*(1+rrhopi))
    // p_eta = (1 + rrhopi - rrhopi*romrho0)/(2*(1+rrhopi))
    double denom = 2.0 * (1.0 + rrhopi);
    double p_pi0   = 1.0 / denom;
    double p_rho0  = rrhopi / denom;
    double p_omega = rrhopi * romrho0 / denom;
    double p_eta   = (1.0 + rrhopi - rrhopi * romrho0) / denom;

    if (rnd < p_pi0) {
        return PIDInference::kPidPi0;   // π0
    } else if (rnd < p_pi0 + p_eta) {
        return PIDInference::kPidEta;   // η
    } else if (rnd < p_pi0 + p_eta + p_rho0) {
        return PIDInference::kPidRho0;   // ρ0
    } else {
        return PIDInference::kPidOmega;   // ω
    }
}

void PIDInference::BatchAssignDiagonalLightMesons(const std::vector<double>& masses,
                                                  std::vector<int>& outPDG,
                                                  int numChargedPions,
                                                  int numChargedRhos) {
    double rrhopi  = PhysicsConstants::kRhoToPionRatio;
    double romrho0 = PhysicsConstants::kOmegaToRhoRatio;

    size_t N = masses.size();
    outPDG.resize(N);

    // 1. Compute expected counts
    double xnpi0   = (numChargedPions + numChargedRhos) / (1.0 + rrhopi) / 2.0;
    double xnrho0  = xnpi0 * rrhopi;
    double xnomega = xnrho0 * romrho0;
    double xneta   = static_cast<double>(N) - xnpi0 - xnrho0 - xnomega;
    if (xneta < 0) xneta = 0;

    // 2. Convert to probabilities
    double p_pi0   = xnpi0   / N;
    double p_eta   = xneta   / N;
    double p_rho0  = xnrho0  / N;
    double p_omega = xnomega / N;

    // 3. Randomly assign each candidate
    for (size_t i = 0; i < N; ++i) {
        double r = dist_(rng_);
        if (r < p_pi0) {
            outPDG[i] = kPidPi0;
        } else if (r < p_pi0 + p_eta) {
            outPDG[i] = kPidEta;
        } else if (r < p_pi0 + p_eta + p_rho0) {
            outPDG[i] = kPidRho0;
        } else {
            outPDG[i] = kPidOmega;
        }
    }
}

// 通用 PID 推断入口
// quarks: 夸克列表，长度为2走介子逻辑，长度为3走重子逻辑
// mass: 系统质量，rnd/vpratio/rrhopi/romrho0: 相关随机数与比率参数
int PIDInference::InferPIDWithRNG(const std::vector<int>& quarks, double mass, double rnd, double vpratio, double rrhopi, double romrho0) {
    if (quarks.size() == 2) {
        int q1 = quarks[0], q2 = quarks[1];
        // 处理所有 quark–antiquark 对角态
        if (q1 == -q2) {
            int af = std::abs(q1);
            if (af <= 2) {
                return ResolveDiagonalLightMeson(rnd, rrhopi, romrho0);
            } else if (af == 3) {
                return PIDInference::kPidPhi;
            } else {
                return InferQuarkoniumPDG(af, mass);
            }
        }
        // 非对角态走自旋 + 基式构造
        int spin = InferMesonSpin(vpratio, rnd);
        int pdg = InferMesonPDG(q1, q2, mass);
        if (spin == 1 && std::abs(pdg % 10) == 1) {
            pdg += 2 * ((pdg > 0) ? 1 : -1);
        }
        return pdg;
    } else if (quarks.size() == 3) {
        return InferBaryonPDG(quarks[0], quarks[1], quarks[2], mass);
    } else {
        return 0;
    }
}

// 生成随机数版本的入口，内部维护 RNG
int PIDInference::InferPID(const std::vector<int>& quarks, double mass) {
    double rnd = dist_(rng_);
    return InferPIDWithRNG(quarks, mass, rnd,
                           kDefaultVPRatio,
                           kDefaultRhoPiRatio,
                           kDefaultOmegaRhoRatio);
}