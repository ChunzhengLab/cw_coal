#ifndef COMBINER_BASE_H
#define COMBINER_BASE_H

#include <memory>
#include <vector>

#include "Particle.h"
#include "core/TimeFrameManager.h"

class CombinerBase {
 public:
  virtual ~CombinerBase() = default;

  /// Combine partons into hadrons using a specific strategy
  virtual std::vector<Hadron*> Combine(const std::vector<Parton*>& partons) = 0;

  /// Final sweep to combine all unused partons into mesons or baryons
  virtual std::vector<Hadron*> Afterburner(const std::vector<Parton*>& partons);

  // Timeframe management
  void SetTimeFrameCount(int numFrames) { timeFrameManager_->SetNumFrames(numFrames); }
  void SetTimeFrameStrategy(TimeFrameManager::Strategy strategy) { timeFrameManager_->SetStrategy(strategy); }
  void SetFixedTimeStep(double timeStep) { timeFrameManager_->SetFixedTimeStep(timeStep); }
  TimeFrameManager* GetTimeFrameManager() { return timeFrameManager_.get(); }

 protected:
  std::unique_ptr<TimeFrameManager> timeFrameManager_;
  CombinerBase() : timeFrameManager_(std::make_unique<TimeFrameManager>()) {}
};

#endif  // COMBINER_BASE_H
