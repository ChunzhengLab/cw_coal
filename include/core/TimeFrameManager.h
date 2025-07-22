#ifndef TIMEFRAME_MANAGER_H
#define TIMEFRAME_MANAGER_H

#include <vector>
#include "Particle.h"

class TimeFrameManager {
 public:
  enum class Strategy {
    kFixedTimeStep,  // 模式1: 用户直接设置每一帧的固定时间步长
    kEqualTime,      // 模式2: 给定帧数，自动提取最大时间，平均分配时间帧
    kAdaptive        // 模式3: 给定帧数，自动提取最大时间，动态调整时间帧
  };

  TimeFrameManager(int numFrames = 10, Strategy strategy = Strategy::kEqualTime, double fixedTimeStep = 1.0);

  // 建立时间帧
  void BuildFrames(const std::vector<Parton*>& partons);
  
  // 获取时间帧边界
  const std::vector<float>& GetFrameBoundaries() const { return frameBoundaries_; }
  
  // 获取指定帧的时间范围
  std::pair<float, float> GetFrameRange(size_t frameIndex) const;
  
  // 计算时间步长
  float GetTimeStep(size_t frameIndex) const;
  
  // 获取指定时间帧内的粒子
  std::vector<Parton*> GetPartonsInFrame(const std::vector<Parton*>& partons, 
                                        size_t frameIndex) const;
  
  // 将粒子移动到下一时间帧
  void MovePartonsToNextFrame(const std::vector<Parton*>& partons, 
                             size_t frameIndex) const;
  
  // 获取时间帧数量
  size_t GetNumFrames() const { return numFrames_; }
  
  // 设置时间帧数量
  void SetNumFrames(int numFrames) { numFrames_ = numFrames; }
  
  // 设置分帧策略
  void SetStrategy(Strategy strategy) { strategy_ = strategy; }
  
  // 设置固定时间步长 (仅用于 kAbsoluteTime 模式)
  void SetFixedTimeStep(double timeStep) { fixedTimeStep_ = timeStep; }
  double GetFixedTimeStep() const { return fixedTimeStep_; }

 private:
  int numFrames_;
  Strategy strategy_;
  double fixedTimeStep_;  // 固定时间步长 (fm/c)
  std::vector<float> frameBoundaries_;
  
  // 不同的分帧策略实现
  void BuildFixedTimeStepFrames(const std::vector<Parton*>& partons);
  void BuildEqualTimeFrames(const std::vector<Parton*>& partons);
  void BuildAdaptiveFrames(const std::vector<Parton*>& partons);
};

#endif // TIMEFRAME_MANAGER_H