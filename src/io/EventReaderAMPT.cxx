#include <fstream>
#include <iostream>
#include "TChain.h"
#include "io/EventReaderAMPT.h"
#include "core/Particle.h"

void EventReaderAMPT::SetupBranches() {
  // Disable all branches by default
  fTree->SetBranchStatus("*", 0);

  // Enable only the branches we need
  fTree->SetBranchStatus("Event", 1);
  fTree->SetBranchStatus("ID",    1);
  fTree->SetBranchStatus("Px",    1);
  fTree->SetBranchStatus("Py",    1);
  fTree->SetBranchStatus("Pz",    1);
  fTree->SetBranchStatus("X",     1);
  fTree->SetBranchStatus("Y",     1);
  fTree->SetBranchStatus("Z",     1);

  // Bind branches to our buffers
  fTree->SetBranchAddress("Event", fEventBuf.data());
  fTree->SetBranchAddress("ID",    fID.data());
  fTree->SetBranchAddress("Px",    fPx.data());
  fTree->SetBranchAddress("Py",    fPy.data());
  fTree->SetBranchAddress("Pz",    fPz.data());
  fTree->SetBranchAddress("X",     fX.data());
  fTree->SetBranchAddress("Y",     fY.data());
  fTree->SetBranchAddress("Z",     fZ.data());
}

EventReaderAMPT::EventReaderAMPT(const std::string& fn) {
  // 支持单文件或.list文件列表
  if (fn.size() >= 5 && fn.substr(fn.size() - 5) == ".list") {
    // 读取文件列表
    TChain* chain = nullptr;
    std::ifstream inputStream(fn);
    if (!inputStream) {
      std::cerr << "Cannot open list file " << fn << std::endl;
      std::exit(1);
    }
    std::string filePath;
    std::string treeName;
    while (std::getline(inputStream, filePath)) {
      if (filePath.empty()) continue;

      // 检测树的名称
      TFile testFile(filePath.c_str(), "READ");
      if (!testFile.IsOpen()) {
        std::cerr << "Cannot open file to check tree: " << filePath << std::endl;
        std::exit(1);
      }

      if (testFile.Get("AMPT")) {
        treeName = "AMPT";
      } else if (testFile.Get("AMPT_I")) {
        treeName = "AMPT_I";
      } else {
        std::cerr << "Neither 'AMPT' nor 'AMPT_I' tree found in " << filePath << std::endl;
        std::exit(1);
      }
      testFile.Close();

      if (!chain) {
        // 第一个文件确定tree名字后，再初始化chain
        chain = new TChain(treeName.c_str());
      }
      std::cout << "Adding file (" << treeName << "): " << filePath << std::endl;
      chain->Add(filePath.c_str());
    }
    fTree = chain;
  } else {
    // 单个 ROOT 文件
    fFile = TFile::Open(fn.c_str());
    if (!fFile || fFile->IsZombie()) {
      std::cerr << "Cannot open ROOT file " << fn << std::endl;
      std::exit(1);
    }
    if (fFile->Get("AMPT")) {
      fTree = static_cast<TTree*>(fFile->Get("AMPT"));
    } else if (fFile->Get("AMPT_I")) {
      fTree = static_cast<TTree*>(fFile->Get("AMPT_I"));
    } else {
      std::cerr << "Neither 'AMPT' nor 'AMPT_I' tree found in " << fn << std::endl;
      std::exit(1);
    }
  }
  
  fNentries = fTree->GetEntries();
  // 调整 vector 大小
  const std::size_t maxN = 100000;
  fID.resize(maxN);
  fPx.resize(maxN);
  fPy.resize(maxN);
  fPz.resize(maxN);
  fX.resize(maxN);
  fY.resize(maxN);
  fZ.resize(maxN);
  SetupBranches();
}

bool EventReaderAMPT::NextEvent(Event& out) {
  if (fCurrent >= fNentries) return false;
  fTree->GetEntry(fCurrent);
  fNparton = fEventBuf[1];
  out.SetUID(static_cast<unsigned int>(fCurrent++)); // 用序号当 ID
  // 清除上一次的数据
  out.Clear();
  // 填充所有 parton
  for (int i = 0; i < fNparton; ++i) {
    // PDG < 0 表示反夸克，B = -1/3；否则表示正夸克，B = +1/3
    int pdg = fID[i];
    double bn = (pdg < 0 ? -1.0/3.0 : 1.0/3.0);
    auto p = new Parton(
      fX[i], fY[i], fZ[i],
      fPx[i], fPy[i], fPz[i],
      bn
    );
    p->SetPID(pdg);
    out.AddParton(p);
  }
  return true;
}

EventReaderAMPT::~EventReaderAMPT() {
  if (fFile) { fFile->Close(); delete fFile; }
  else if (TChain* chain = dynamic_cast<TChain*>(fTree)) {
    delete chain;
  }
}