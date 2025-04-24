#include "core/Particle.h"
#include "core/Event.h"
#include "io/EventReaderAMPT.h"
#include "io/EventRandomGen.h"
#include "io/EventWriter.h"
#include "ana/AnalyzerQA.h"
#include "ana/AnalyzerCVE.h"
#include "Combiners.h"
#include <TROOT.h>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <functional>

static void PrintUsage() {
    std::cout << "Usage: cwcoal [options]\n"
              << "  -h, --help               Show this help message\n"
              << "  -i, --data-input <file>  AMPT input ROOT file or list (if omitted, random generation mode)\n"
              << "  -o, --data-output <file> Output ROOT file for hadrons (default: output.root)\n"
              << "  -a, --algorithm <name>   Combiner algorithm: KDTreeGlobal, KDTreeGreedy, BruteForceGlobal, BruteForceGreedy (default: KDTreeGlobal)\n"
              << "  -n, --events <N>         Number of events to process/generate (if omitted, process all events)\n"
              << "  -b, --bn <B>             Target total baryon number per event (default: 0)\n"
              << "  -s, --savedir <dir>      Output directory for all files (default: current working directory)\n"
              << "  -r, --baryon-preference <R>  Baryon preference factor (default: 1.0)\n";
}

int main(int argc, char** argv) {
    // Variables for input and output file management
    std::string dataInput;
    std::string dataOutput = "output.root";
    std::string algorithm = "KDTreeGlobal";
    std::string saveDir = "."; // Output directory for all files
    int nEvents = 10;
    int sumBn = 0;
    double baryonPreference = 1.0;
    bool eventsLimited = false;

    const struct option longOpts[] = {
        {"help",      no_argument,       nullptr, 'h'},
        {"data-input", required_argument, nullptr, 'i'},
        {"data-output", required_argument, nullptr, 'o'},
        {"algorithm", required_argument, nullptr, 'a'},
        {"events",    required_argument, nullptr, 'n'},
        {"bn",        required_argument, nullptr, 'b'},
        {"savedir",   required_argument, nullptr, 's'},
        {"baryon-preference", required_argument, nullptr, 'r'},
        {nullptr,     0,                 nullptr,  0 }
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "h:i:o:a:n:b:s:r:", longOpts, nullptr)) != -1) {
        switch (opt) {
            case 'h': PrintUsage(); return 0;
            case 'i': dataInput = optarg; break;
            case 'o': dataOutput = optarg; break;
            case 'a': algorithm = optarg; break;
            case 'n': nEvents = std::stoi(optarg); eventsLimited = true; break;
            case 'b': sumBn = std::stoi(optarg); break;
            case 's': saveDir = optarg; break;
            case 'r': baryonPreference = std::stod(optarg); break;
            default:  PrintUsage(); return 1;
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
    } else {
        std::cerr << "Unknown algorithm: " << algorithm << std::endl;
        return 1;
    }

    // Initialize writer and analyzers
    std::unique_ptr<EventWriter> writer = std::make_unique<EventWriter>(saveDir + "/" + dataOutput);
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
              << (dataInput.empty()
                    ? "Random generation mode"
                    : std::string("AMPT input mode (file: ") + dataInput + ")")
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
        std::cout << ">>>Number of events to process: " << nEvents << std::endl;
        auto gen = std::make_shared<EventRandomGen>();
        fetch = [gen](Event& e) { e = gen->GenerateEvent(-1, 0); return true; };
    }

    std::cout << ">>>Baryon preference factor: " << baryonPreference << std::endl;
    std::cout << ">>>Save directory: " << saveDir << std::endl;
    std::cout << ">>>Algorithm: " << algorithm << std::endl;
    if (writer) {
        std::cout << ">>>Hadrons output file: " << saveDir << "/" << dataOutput << std::endl;
    }

    // Progress bar
    auto printProgress = [&](int current) {
        int width = 80;
        int pos = static_cast<int>(width * current / nEvents);
        std::cout << "\r[";
        for (int i = 0; i < width; ++i) {
            std::cout << (i < pos ? '=' : ' ');
        }
        std::cout << "] " << static_cast<int>(100.0 * current / nEvents)
                  << "% (" << current << "/" << nEvents << ")" << std::flush;
        if (current == nEvents) std::cout << std::endl;
    };

    // Main loop
    for (int ie = 0; ie < nEvents; ++ie) {
        Event evt;
        if (!fetch(evt)) break;
        auto partons = evt.GetPartons();
        auto hadrons = combiner->Combine(partons);
        for (auto* h : hadrons) evt.AddHadron(h);
        if (writer) writer->WriteEvent(&evt);
        qa.Process(evt);
        cve.Process(evt);
        printProgress(ie + 1);
    }

    // Finalize
    qa.Finish(saveDir + "/qa_" + algorithm + ".root");
    cve.Finish(saveDir + "/cve_" + algorithm + ".root");
    if (writer) writer->Close();

    return 0;
}