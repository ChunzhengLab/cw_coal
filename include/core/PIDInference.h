#ifndef PIDINFERENCE_H
#define PIDINFERENCE_H

#include <random>
#include <vector>

#include "core/PhysicsConstants.h"

class PIDInference {
 public:
  // Common PDG codes for light mesons
  static constexpr int kPidPi0 = 111;
  static constexpr int kPidEta = 221;
  static constexpr int kPidRho0 = 113;
  static constexpr int kPidOmega = 223;
  static constexpr int kPidPhi = 333;

  // Default ratio parameters (from PhysicsConstants)
  static constexpr double kDefaultVPRatio = PhysicsConstants::kMesonVectorToPseudoscalarRatio;
  static constexpr double kDefaultRhoPiRatio = PhysicsConstants::kRhoToPionRatio;
  static constexpr double kDefaultOmegaRhoRatio = PhysicsConstants::kOmegaToRhoRatio;

  // 根据两个夸克/反夸克的 flavor code 和系统质量，推断生成的介子 PDG 码
  // q1, q2: ±1(d/d̄), ±2(u/ū), ±3(s/ s̄), ...
  // mass: 两体系统的不变量质量
  static int InferMesonPDG(int q1, int q2, double mass);

  // 根据三个夸克的 flavor code、系统质量 和自旋多重度，推断生成的重子或反重子 PDG 码
  // q1, q2, q3: 夸克或反夸克代码；mass: 三体不变量质量；spinMult: 自旋多重度(2S+1)，-1 表示自动判断
  static int InferBaryonPDG(int q1, int q2, int q3, double mass, int spinMult = -1);

  // 根据重夸克 flavo(u)r 和质量，从预定义偶素列表中选最接近的 quarkonia PDG 码
  // flavor: 4 (c), 5 (b), 6 (t)；mass: 对角态质量
  static int InferQuarkoniumPDG(int flavor, double mass);

  // 基于矢量/伪标量介子产额比 (V/P) 和随机数，推断介子的自旋态
  // vpratio: 矢量介子/标量介子产额比；rnd: [0,1) 随机数；返回 0 (伪标量) 或 1 (向量)
  static int InferMesonSpin(double vpratio, double rnd);

  // 针对轻味对角态介子 (uū、d d̄)，根据 ρ⁰/π⁰、ω/ρ⁰ 比率及随机数，分配 π⁰、η、ρ⁰ 或 ω
  // rnd: 随机数；rrhopi: ρ⁰/π⁰ 比；romrho0: ω/ρ⁰ 比
  static int ResolveDiagonalLightMeson(double rnd, double rrhopi, double romrho0);

  // 通用 PID 推断入口，根据夸克数目自动调用介子或重子推断逻辑
  // quarks: 输入夸克列表；mass: 系统质量
  // 生成随机数版本的 InferPID，内部维护 RNG，不需要传入比率参数
  static int InferPID(const std::vector<int>& quarks, double mass);

  // Batch-assign light diagonal mesons using event-level charged π and ρ counts
  static void BatchAssignDiagonalLightMesons(const std::vector<double>& masses, std::vector<int>& outPDG,
                                             int numChargedPions, int numChargedRhos);

 private:
  // Internal PID inference with explicit random and ratio parameters
  static int InferPIDWithRNG(const std::vector<int>& quarks, double mass, double rnd, double vpratio, double rrhopi,
                             double romrho0);
  // 内部随机数生成器和分布
  static std::mt19937 rng_;
  static std::uniform_real_distribution<double> dist_;
};

#endif  // PIDINFERENCE_H
