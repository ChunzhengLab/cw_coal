#include <TROOT.h>
#include <getopt.h>

#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

#include "Combiners.h"
#include "ana/AnalyzerCVE.h"
#include "ana/AnalyzerQA.h"
#include "core/Event.h"
#include "core/PIDAssigner.h"
#include "core/TimeFrameManager.h"
#include "io/EventRandomGen.h"
#include "io/EventReaderAMPT.h"
#include "io/EventWriter.h"

using SamplingMode = EventRandomGen::SamplingMode;

static void PrintUsage() {
  std::cout << "Usage: cwcoal [options]\n"
            << "  -h, --help               Show this help message\n"
            << "  -i, --data-input <file>  AMPT input ROOT file or list (if omitted, random generation mode)\n"
            << "  -o, --data-output <file> Output ROOT file for hadrons (if omitted, no output will be written)\n"
            << "  -a, --algorithm <name>   Combiner algorithm: KDTreeGlobal, KDTreeGreedy, BruteForceGlobal, "
               "BruteForceGreedy (default: KDTreeGlobal)\n"
            << "  -n, --events <N>         Number of events to process/generate (if omitted, process all events)\n"
            << "  -b, --bn <B>             Target total baryon number per event (default: 0)\n"
            << "  -p, --partons <P>        Number of partons per event (-1 to sample from histogram)\n"
            << "  -s, --savedir <dir>      Output directory for all files (default: current working directory)\n"
            << "  -r, --baryon-preference <R>  Baryon preference factor (default: 1.0)\n"
            << "  -F, --shuffle-fraction <F>   Shuffle fraction of parton positions (0.0–1.0, default: 0.0)\n"
            << "  -T, --toymode            Use Toy event generation mode (uniform disk + simple pT)\n"
            << "  --timeframes <N>         Number of time frames for evolution (default: 10)\n"
            << "  --timeframe-strategy <S> Time frame strategy: FixedTimeStep, EqualTime, Adaptive (default: EqualTime)\n"
            << "  --fixed-timestep <dt>    Fixed time step in fm/c (only for FixedTimeStep strategy, default: 1.0)\n";
}

int main(int argc, char** argv) {
  // Variables for input and output file management
  std::string dataInput;
  std::string dataOutput;
  bool outputEnabled = false;  // Only write events if --data-output is specified
  std::string algorithm = "KDTreeGlobal";
  std::string saveDir = ".";  // Output directory for all files
  int nEvents = 10;
  int sumBn = 0;
  int partonCount = -1;  // number of partons per event (-1 to sample from histogram)
  double baryonPreference = 1.0;
  bool eventsLimited = false;
  bool toyMode = false;

  // Shuffle feature
  double shuffleFraction = -1;
  
  // Timeframe features
  int timeFrames = 10;
  std::string timeFrameStrategy = "EqualTime";
  double fixedTimeStep = 1.0;

  const struct option longOpts[] = {{"help", no_argument, nullptr, 'h'},
                                    {"data-input", required_argument, nullptr, 'i'},
                                    {"data-output", required_argument, nullptr, 'o'},
                                    {"algorithm", required_argument, nullptr, 'a'},
                                    {"events", required_argument, nullptr, 'n'},
                                    {"bn", required_argument, nullptr, 'b'},
                                    {"partons", required_argument, nullptr, 'p'},
                                    {"savedir", required_argument, nullptr, 's'},
                                    {"baryon-preference", required_argument, nullptr, 'r'},
                                    {"shuffle-fraction", required_argument, nullptr, 'F'},
                                    {"toymode", no_argument, nullptr, 'T'},
                                    {"timeframes", required_argument, nullptr, 1000},
                                    {"timeframe-strategy", required_argument, nullptr, 1001},
                                    {"fixed-timestep", required_argument, nullptr, 1002},
                                    {nullptr, 0, nullptr, 0}};

  int opt;
  while ((opt = getopt_long(argc, argv, "hi:o:a:n:b:p:s:r:F:T", longOpts, nullptr)) != -1) {
    switch (opt) {
      case 'h':
        PrintUsage();
        return 0;
      case 'i':
        dataInput = optarg;
        break;
      case 'o':
        dataOutput = optarg;
        outputEnabled = true;
        break;
      case 'a':
        algorithm = optarg;
        break;
      case 'n':
        nEvents = std::stoi(optarg);
        eventsLimited = true;
        break;
      case 'b':
        sumBn = std::stoi(optarg);
        break;
      case 'p':
        partonCount = std::stoi(optarg);
        break;
      case 's':
        saveDir = optarg;
        break;
      case 'r':
        baryonPreference = std::stod(optarg);
        break;
      case 'F':
        shuffleFraction = std::stod(optarg);
        break;
      case 'T':
        toyMode = true;
        break;
      case 1000:
        timeFrames = std::stoi(optarg);
        break;
      case 1001:
        timeFrameStrategy = optarg;
        break;
      case 1002:
        fixedTimeStep = std::stod(optarg);
        break;
      default:
        PrintUsage();
        return 1;
    }
  }

  // Select combiner
  std::unique_ptr<CombinerBase> combiner;
  if (algorithm == "KDTreeGlobal") {
    combiner = std::make_unique<KDTreeGlobal>(baryonPreference);
  } else if (algorithm == "KDTreeGreedy") {
    combiner = std::make_unique<KDTreeGreedy>(baryonPreference);
  } else if (algorithm == "BruteForceGlobal") {
    combiner = std::make_unique<BruteForceGlobal>(baryonPreference);
  } else if (algorithm == "BruteForceGreedy") {
    combiner = std::make_unique<BruteForceGreedy>(baryonPreference);
  } else if (algorithm == "BruteForceDualGreedy") {
    combiner = std::make_unique<BruteForceDualGreedy>(baryonPreference);
  } else if (algorithm == "KDTreeDualGreedy") {
    combiner = std::make_unique<KDTreeDualGreedy>(baryonPreference);
  } 
  else {
    std::cerr << "Unknown algorithm: " << algorithm << std::endl;
    return 1;
  }

  // Configure timeframe manager
  combiner->SetTimeFrameCount(timeFrames);
  combiner->SetFixedTimeStep(fixedTimeStep);
  
  // Set timeframe strategy
  TimeFrameManager::Strategy strategy;
  if (timeFrameStrategy == "FixedTimeStep") {
    strategy = TimeFrameManager::Strategy::kFixedTimeStep;
  } else if (timeFrameStrategy == "EqualTime") {
    strategy = TimeFrameManager::Strategy::kEqualTime;
  } else if (timeFrameStrategy == "Adaptive") {
    strategy = TimeFrameManager::Strategy::kAdaptive;
  } else {
    std::cerr << "Unknown timeframe strategy: " << timeFrameStrategy << std::endl;
    return 1;
  }
  combiner->SetTimeFrameStrategy(strategy);

  // Initialize writer and analyzers
  std::unique_ptr<EventWriter> writer;
  if (outputEnabled) { writer = std::make_unique<EventWriter>(saveDir + "/" + dataOutput); }
  AnalyzerQA qa;
  AnalyzerCVE cve;
  qa.Init();
  cve.Init();

  // Startup banner
  std::cout << "==================================================================================\n";
  std::cout << "                  Chunzheng Wang's Quark Coalescence Model \n";
  std::cout << "             Author: Chunzheng Wang (chunzheng.wang@icloud.com) \n";
  std::cout << "==================================================================================\n";

  // Mode and event count summary
  std::cout << ">>>Mode: "
            << (dataInput.empty() ? "Random generation mode" : std::string("AMPT input mode (file: ") + dataInput + ")")
            << std::endl;

  Long64_t totalEvents = -1;
  std::function<bool(Event&)> fetch;
  if (!dataInput.empty()) {
    auto reader = std::make_shared<EventReaderAMPT>(dataInput);
    totalEvents = reader->GetTotalEvents();
    if (eventsLimited) {
      std::cout << ">>>Number of events to process: " << nEvents << std::endl;
    } else {
      std::cout << ">>>Number of events to process: all available events (" << totalEvents << ")" << std::endl;
      nEvents = static_cast<int>(totalEvents);
    }
    fetch = [reader](Event& e) { return reader->NextEvent(e); };
  } else {
    auto gen = std::make_shared<EventRandomGen>();
    // Determine sampling mode
    SamplingMode mode = toyMode ? SamplingMode::kToyMode : SamplingMode::kSampleFromFile;
    // Single lambda for both modes
    fetch = [gen, sumBn, mode, partonCount](Event& e) {
      gen->GenerateEvent(e, partonCount, sumBn, mode);
      return true;
    };
  }

  std::cout << ">>>Baryon preference factor: " << baryonPreference << std::endl;
  std::cout << ">>>Save directory: " << saveDir << std::endl;
  std::cout << ">>>Algorithm: " << algorithm << std::endl;
  std::cout << ">>>Time frames: " << timeFrames << " (" << timeFrameStrategy << ")" << std::endl;
  if (timeFrameStrategy == "FixedTimeStep") {
    std::cout << ">>>Fixed time step: " << fixedTimeStep << " fm/c" << std::endl;
  }
  if (writer) { std::cout << ">>>Hadrons output file: " << saveDir << "/" << dataOutput << std::endl; }
  if (shuffleFraction >= 0.0) {
    std::cout << ">>>Shuffle fraction of parton positions: " << shuffleFraction << std::endl;
  }

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

  // Main loop
  for (int ie = 0; ie < nEvents; ++ie) {
    Event evt;
    if (!fetch(evt)) break;
    if (shuffleFraction > 0.0) { evt.ShufflePartons(shuffleFraction); }
    auto partons = evt.GetPartons();
    auto hadrons = combiner->Combine(partons);
    for (auto* h : hadrons) evt.AddHadron(h);
    // Assign PDG codes to hadrons
    PIDAssigner::Assign(evt);
    if (writer) writer->WriteEvent(&evt);
    qa.Process(evt);
    cve.Process(evt);
    printProgress(ie + 1);
  }

  // Finalize
  std::ostringstream ossQA, ossCVE;
  ossQA << saveDir << "/qa_" << algorithm << "_r" << std::fixed << std::setprecision(2) << baryonPreference;
  ossCVE << saveDir << "/cve_" << algorithm << "_r" << std::fixed << std::setprecision(2) << baryonPreference;
  if (toyMode) {
    ossQA << "_n" << nEvents << "_p" << partonCount << "_bn" << sumBn;
    ossCVE << "_n" << nEvents << "_p" << partonCount << "_bn" << sumBn;
    if (shuffleFraction >= 0.0) {
      ossQA << "_sf" << std::fixed << std::setprecision(2) << shuffleFraction;
      ossCVE << "_sf" << std::fixed << std::setprecision(2) << shuffleFraction;
    }
  }
  ossQA << ".root";
  ossCVE << ".root";
  qa.Finish(ossQA.str());
  cve.Finish(ossCVE.str());
  if (writer) writer->Close();

  return 0;
}
