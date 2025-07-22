#ifndef EVENTREADERAMPT_H
#define EVENTREADERAMPT_H

#include <array>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "core/Event.h"

class EventReaderAMPT {
 public:
  // 打开文件并绑定所有分支
  EventReaderAMPT(const std::string& filename);
  ~EventReaderAMPT();

  // 读取下一条事件，返回 false 表示读完了
  bool NextEvent(Event& out);
  /// Return the total number of events available in the input
  Long64_t GetTotalEvents() const {
    return fNentries;
  }

 private:
  void SetupBranches();

  TFile* fFile = nullptr;
  TTree* fTree = nullptr;
  Long64_t fNentries = 0;
  Long64_t fCurrent = 0;

  // 临时存放从 ROOT 读到的分支数据
  Int_t fEvent = 0;
  Int_t fNparton = 0;
  std::array<Int_t, 2> fEventBuf{};
  std::vector<Int_t> fID;
  std::vector<Float_t> fPx, fPy, fPz;
  std::vector<Float_t> fX, fY, fZ;
  std::vector<Float_t> fT;
};

#endif  // EVENTREADERAMPT_H
