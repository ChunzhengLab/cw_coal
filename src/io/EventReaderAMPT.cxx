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
  fTree = nullptr;
  fFile = nullptr;
  fNentries = 0;

  if (fn.size() >= 5 && fn.substr(fn.size() - 5) == ".list") {
    TChain* chain = nullptr;
    std::ifstream inputStream(fn);
    if (!inputStream) {
      std::cerr << "[ERROR] Cannot open list file: " << fn << std::endl;
      return;
    }

    std::string filePath;
    std::string treeName;
    while (std::getline(inputStream, filePath)) {
      if (filePath.empty()) continue;

      TFile testFile(filePath.c_str(), "READ");
      if (!testFile.IsOpen() || testFile.IsZombie()) {
        std::cerr << "[WARN] Skipping unreadable or zombie file: " << filePath << std::endl;
        continue;
      }

      if (testFile.GetNkeys() == 0) {
        std::cerr << "[WARN] Skipping empty file (no keys): " << filePath << std::endl;
        testFile.Close();
        continue;
      }

      TTree* testTree = nullptr;
      if ((testTree = static_cast<TTree*>(testFile.Get("AMPT")))) {
        treeName = "AMPT";
      } else if ((testTree = static_cast<TTree*>(testFile.Get("AMPT_I")))) {
        treeName = "AMPT_I";
      } else {
        std::cerr << "[WARN] No AMPT/AMPT_I tree found in: " << filePath << std::endl;
        testFile.Close();
        continue;
      }

      testFile.Close();

      if (!chain) {
        chain = new TChain(treeName.c_str());
      }
      std::cout << "[INFO] Adding file (" << treeName << "): " << filePath << std::endl;
      chain->Add(filePath.c_str());
    }

    if (!chain || chain->GetNtrees() == 0) {
      std::cerr << "[ERROR] No valid ROOT files found in list." << std::endl;
      return;
    }

    fTree = chain;
  } else {
    // 单个文件模式
    fFile = TFile::Open(fn.c_str(), "READ");
    if (!fFile || fFile->IsZombie()) {
      std::cerr << "[ERROR] Cannot open ROOT file: " << fn << std::endl;
      return;
    }

    if (fFile->GetNkeys() == 0) {
      std::cerr << "[ERROR] ROOT file is empty (no keys): " << fn << std::endl;
      return;
    }

    if ((fTree = static_cast<TTree*>(fFile->Get("AMPT")))) {
      // OK
    } else if ((fTree = static_cast<TTree*>(fFile->Get("AMPT_I")))) {
      // OK
    } else {
      std::cerr << "[ERROR] No AMPT/AMPT_I tree found in: " << fn << std::endl;
      fTree = nullptr;
      return;
    }
  }

  if (!fTree) {
    std::cerr << "[ERROR] No valid tree loaded." << std::endl;
    return;
  }

  fNentries = fTree->GetEntries();

  const std::size_t maxN = 100000;
  fID.resize(maxN);
  fPx.resize(maxN);
  fPy.resize(maxN);
  fPz.resize(maxN);
  fX.resize(maxN);
  fY.resize(maxN);
  fZ.resize(maxN);

  SetupBranches();
  std::cout << "[INFO] EventReaderAMPT initialized with " << fNentries << " entries." << std::endl;
}

bool EventReaderAMPT::NextEvent(Event& out) {
  if (fCurrent >= fNentries) return false;
  fTree->GetEntry(fCurrent);
  fNparton = fEventBuf[1];
  out.SetUID(static_cast<unsigned int>(fCurrent++)); // 用序号当 ID
  // 清除上一次的数据
  out.Reset(nullptr);
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
