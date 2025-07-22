#include "core/TimeFrameManager.h"
#include <algorithm>
#include <cmath>

TimeFrameManager::TimeFrameManager(int numFrames, Strategy strategy, double fixedTimeStep) 
    : numFrames_(numFrames), strategy_(strategy), fixedTimeStep_(fixedTimeStep) {
}

void TimeFrameManager::BuildFrames(const std::vector<Parton*>& partons) {
  frameBoundaries_.clear();
  if (partons.empty() || numFrames_ <= 0) return;
  
  switch (strategy_) {
    case Strategy::kFixedTimeStep:
      BuildFixedTimeStepFrames(partons);
      break;
    case Strategy::kEqualTime:
      BuildEqualTimeFrames(partons);
      break;
    case Strategy::kAdaptive:
      BuildAdaptiveFrames(partons);
      break;
  }
}

std::pair<float, float> TimeFrameManager::GetFrameRange(size_t frameIndex) const {
  if (frameIndex >= frameBoundaries_.size() - 1) {
    return {0.0f, 0.0f};
  }
  return {frameBoundaries_[frameIndex], frameBoundaries_[frameIndex + 1]};
}

float TimeFrameManager::GetTimeStep(size_t frameIndex) const {
  if (frameIndex >= frameBoundaries_.size() - 1) {
    return 0.0f;
  }
  return frameBoundaries_[frameIndex + 1] - frameBoundaries_[frameIndex];
}

std::vector<Parton*> TimeFrameManager::GetPartonsInFrame(
    const std::vector<Parton*>& partons, size_t frameIndex) const {
  auto [tLow, tHigh] = GetFrameRange(frameIndex);
  
  std::vector<Parton*> framePartons;
  for (auto* p : partons) {
    if (!p->IsUsed() && p->T() >= tLow && p->T() < tHigh) {
      framePartons.push_back(p);
    }
  }
  return framePartons;
}

void TimeFrameManager::MovePartonsToNextFrame(
    const std::vector<Parton*>& partons, size_t frameIndex) const {
  float deltaTime = GetTimeStep(frameIndex);
  for (auto* p : partons) {
    if (!p->IsUsed()) {
      p->MoveOn(deltaTime);
    }
  }
}

void TimeFrameManager::BuildFixedTimeStepFrames(const std::vector<Parton*>& partons) {
  // 模式1: 用户直接设置每一帧的固定时间步长
  if (partons.empty()) return;
  
  // 找到最小时间作为起点
  float minTime = static_cast<float>(partons[0]->T());
  for (auto* p : partons) {
    minTime = std::min(minTime, static_cast<float>(p->T()));
  }
  
  // 按固定时间步长建立时间帧
  // 每帧的时间间隔都是 fixedTimeStep_
  for (int i = 0; i <= numFrames_; ++i) {
    frameBoundaries_.push_back(minTime + i * fixedTimeStep_);
  }
}

void TimeFrameManager::BuildEqualTimeFrames(const std::vector<Parton*>& partons) {
  // 模式2: 给定帧数，自动提取最大时间，平均分配时间帧
  if (partons.empty()) return;
  
  float minTime = static_cast<float>(partons[0]->T());
  float maxTime = static_cast<float>(partons[0]->T());
  
  for (auto* p : partons) {
    minTime = std::min(minTime, static_cast<float>(p->T()));
    maxTime = std::max(maxTime, static_cast<float>(p->T()));
  }
  
  // 按时间均匀分配
  float timeRange = maxTime - minTime;
  float timeStep = timeRange / numFrames_;
  
  for (int i = 0; i <= numFrames_; ++i) {
    frameBoundaries_.push_back(minTime + i * timeStep);
  }
}

void TimeFrameManager::BuildAdaptiveFrames(const std::vector<Parton*>& partons) {
  // 模式3: 给定帧数，自动提取最大时间，动态调整时间帧
  if (partons.empty()) return;
  
  // 先按时间排序
  std::vector<Parton*> sorted(partons.begin(), partons.end());
  std::sort(sorted.begin(), sorted.end(),
            [](Parton* a, Parton* b) { return a->T() < b->T(); });
  
  float minTime = sorted[0]->T();
  float maxTime = sorted.back()->T();
  float timeRange = maxTime - minTime;
  
  // 动态调整策略：根据粒子时间分布密度调整时间帧大小
  // 使用变分时间步长：粒子密集的区域用更小的时间步长
  
  frameBoundaries_.push_back(minTime);
  
  const size_t totalPartons = sorted.size();
  size_t processedPartons = 0;
  
  for (int frame = 1; frame <= numFrames_; ++frame) {
    if (frame == numFrames_) {
      // 最后一帧包含所有剩余时间
      frameBoundaries_.push_back(maxTime + 1.e-6f);
      break;
    }
    
    // 计算目标粒子数 (理想情况下每帧相等)
    size_t targetPartons = (totalPartons * frame) / numFrames_;
    
    // 找到对应的时间点
    if (targetPartons < totalPartons) {
      float targetTime = sorted[targetPartons]->T();
      
      // 动态调整：如果时间间隔太小，向后推移
      float minTimeStep = timeRange / (numFrames_ * 10); // 最小时间步长
      if (targetTime - frameBoundaries_.back() < minTimeStep) {
        targetTime = frameBoundaries_.back() + minTimeStep;
      }
      
      frameBoundaries_.push_back(targetTime);
    } else {
      // 均匀分配剩余时间
      float remainingTime = maxTime - frameBoundaries_.back();
      int remainingFrames = numFrames_ - frame + 1;
      frameBoundaries_.push_back(frameBoundaries_.back() + remainingTime / remainingFrames);
    }
  }
}