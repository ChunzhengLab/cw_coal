#ifndef EVENT_READER_H
#define EVENT_READER_H

#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "core/Event.h"

class EventReader {
 public:
  enum CopyMode { kShallowCopy, kDeepCopy };

  explicit EventReader(const std::string& filename, CopyMode copyMode = kShallowCopy);
  ~EventReader();

  /// Returns a pointer to the next Event in the file.
  /// For shallow copy mode, returns internal pointer (valid until next call).
  /// For deep copy mode, returns a new Event pointer (caller must delete).
  Event* NextEvent();

  void Close();

  /// Get total number of events available in the file or chain
  Long64_t GetTotalEvents() const {
    return nEntries_;
  }

 private:
  CopyMode copyMode_;
  TFile* file_;
  TTree* tree_;
  Event* eventBuffer_;  // buffer for shallow copy
  Long64_t nEntries_;
  Long64_t currentEntry_;
};

#endif  // EVENT_READER_H
