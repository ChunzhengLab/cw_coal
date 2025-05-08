#include "io/EventReader.h"

#include <fstream>

#include "TChain.h"
#include "core/Particle.h"

EventReader::EventReader(const std::string& filename, EventReader::CopyMode copyMode)
    : file_(nullptr), tree_(nullptr), eventBuffer_(nullptr), copyMode_(copyMode), nEntries_(0), currentEntry_(0) {
  // Use TChain to support reading a .list of files or a single ROOT file
  TChain* chain = new TChain("events");
  if (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".list") {
    std::ifstream in(filename);
    std::string fn;
    while (std::getline(in, fn)) {
      if (fn.empty()) continue;
      chain->Add(fn.c_str());
    }
    // no TFile opened
  } else {
    chain->Add(filename.c_str());
  }
  tree_ = chain;

  eventBuffer_ = nullptr;
  tree_->SetBranchAddress("event", &eventBuffer_);
  nEntries_ = tree_->GetEntries();
}

EventReader::~EventReader() {
  Close();
  // If file_ is not used, tree_ is a TChain allocated here, so delete it
  if (!file_ && tree_) {
    delete tree_;
    tree_ = nullptr;
  }
}

Event* EventReader::NextEvent() {
  if (currentEntry_ >= nEntries_) { return nullptr; }

  tree_->GetEntry(currentEntry_);
  ++currentEntry_;

  if (copyMode_ == EventReader::kShallowCopy) {
    return eventBuffer_;
  } else {
    // Deep copy: create new Event and clone particles
    Event* copy = new Event();
    copy->SetUID(eventBuffer_->GetUID());
    copy->SetReactionPlane(eventBuffer_->GetReactionPlane());

    for (auto* p : eventBuffer_->GetPartons()) { copy->AddParton(new Parton(*p)); }
    for (auto* h : eventBuffer_->GetHadrons()) { copy->AddHadron(new Hadron(*h)); }
    return copy;
  }
}

void EventReader::Close() {
  if (file_) {
    file_->Close();
    delete file_;
    file_ = nullptr;
  }
  eventBuffer_ = nullptr;
  tree_ = nullptr;
  currentEntry_ = 0;
  nEntries_ = 0;
}
