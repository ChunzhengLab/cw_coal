#include <getopt.h>

#include <deque>
#include <iostream>
#include <string>
#include <vector>

#include "ana/AnalyzerCVE.h"
#include "ana/AnalyzerQA.h"
#include "core/Event.h"
#include "io/EventReader.h"

static void PrintUsage() {
  std::cout << "Usage: analysis [options]\n"
            << "Options:\n"
            << "  -h, --help               Show this help message and exit\n"
            << "  -i, --data-input <file>  Input data file or list\n"
            << "  -s, --savedir <dir>      Output directory for all analysis files (default current directory)\n"
            << "  -m, --is-mix             Enable event mixing\n"
            << "  -p, --mixpool-size <N>   Size of the mixing pool (default 5)\n";
}

int main(int argc, char** argv) {
  std::string dataInput;
  bool doMix = false;
  std::string saveDir = ".";
  size_t mixPoolSize = 2;
  size_t totalReadEvents = 0;

  static struct option longOpts[] = {
      {"help", no_argument, nullptr, 'h'},          {"data-input", required_argument, nullptr, 'i'},
      {"is-mix", no_argument, nullptr, 'm'},        {"mixpool-size", required_argument, nullptr, 'p'},
      {"savedir", required_argument, nullptr, 's'}, {nullptr, 0, nullptr, 0}};

  int opt;
  int longIndex = 0;
  while ((opt = getopt_long(argc, argv, "hi:mp:s:", longOpts, &longIndex)) != -1) {
    switch (opt) {
      case 'h':
        PrintUsage();
        return 0;
      case 'i':
        dataInput = optarg;
        break;
      case 'm':
        doMix = true;
        break;
      case 'p':
        mixPoolSize = std::stoul(optarg);
        break;
      case 's':
        saveDir = optarg;
        break;
      default:
        PrintUsage();
        return 1;
    }
  }

  if (dataInput.empty()) {
    PrintUsage();
    return 1;
  }

  // Reader for cwcoal-generated data
  // Use deep‑copy mode so we can store events in mixPool without them being overwritten
  EventReader reader(dataInput, doMix ? EventReader::kDeepCopy : EventReader::kShallowCopy);

  // Startup banner
  std::cout << "==================================================================================\n";
  std::cout << "              Chunzheng Wang's Quark Coalescence Model Analyzer \n";
  std::cout << "             Author: Chunzheng Wang (chunzheng.wang@icloud.com) \n";
  std::cout << "==================================================================================\n";

  // Determine total events for progress bar
  size_t nEvents = reader.GetTotalEvents();  // replace with actual method if different

  std::cout << ">>>Input file: " << dataInput << std::endl;
  std::cout << ">>>Number of events to process: " << nEvents << std::endl;
  if (doMix) {
    std::cout << ">>> Mix event enabled: Yes"
              << ", mix pool size: " << mixPoolSize << std::endl;
  } else {
    std::cout << ">>> Mix event enabled: No" << std::endl;
  }
  std::cout << ">>>Save directory: " << saveDir << std::endl;

  // Progress bar
  auto printProgress = [&](int current) {
    int width = 80;
    int pos = static_cast<int>(width * current / nEvents);
    std::cout << "\r[";
    for (int i = 0; i < width; ++i) { std::cout << (i < pos ? '=' : ' '); }
    std::cout << "] " << static_cast<int>(100.0 * current / nEvents) << "% (" << current << "/" << nEvents << ")"
              << std::flush;
    if (current == nEvents) std::cout << std::endl;
  };

  // Offline analyzers
  AnalyzerCVE cveSame;
  cveSame.Init();
  AnalyzerQA qa;
  qa.Init();

  AnalyzerCVE cveMix;
  std::deque<Event*> mixPool;

  if (doMix) {
    cveMix.SetProcessMixed(true);
    cveMix.Init();
  }

  // Main event loop: handle same-event, mixed-event, and QA
  Event* evt = nullptr;
  while ((evt = reader.NextEvent())) {
    ++totalReadEvents;
    printProgress(totalReadEvents);
    // same-event analysis
    cveSame.Process(*evt);
    // QA histogram filling
    qa.Process(*evt);
    if (doMix) {
      // mixed-event analysis against events in pool
      std::vector<const Event*> bgEvts(mixPool.begin(), mixPool.end());
      cveMix.ProcessMixed(*evt, bgEvts);
      // add current event to pool (deep-copy mode ensures independence)
      mixPool.push_back(evt);
      if (mixPool.size() > mixPoolSize) {
        delete mixPool.front();
        mixPool.pop_front();
      }
    }
  }

  // clean up remaining events in pool
  if (doMix) {
    for (auto* e : mixPool) { delete e; }
  }

  // Write analysis outputs
  cveSame.Finish(saveDir + "/cve_single_offline.root");
  if (doMix) cveMix.Finish(saveDir + "/cve_mix_offline.root");
  qa.Finish(saveDir + "/qa_offline.root");

  if (doMix && totalReadEvents < mixPoolSize + 1) {
    std::cerr << "Warning: only " << totalReadEvents << " events read, less than mixpool-size+1 (" << (mixPoolSize + 1)
              << "), no mixed-event analysis was performed.\n";
  }

  return 0;
}