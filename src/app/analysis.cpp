#include <iostream>
#include <string>
#include <deque>
#include <vector>
#include <getopt.h>
#include "io/EventReader.h"
#include "ana/AnalyzerCVE.h"
#include "ana/AnalyzerQA.h"
#include "core/Event.h"

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
    size_t mixPoolSize = 5;
    size_t totalReadEvents = 0;

    static struct option longOpts[] = {
        {"help",        no_argument,       nullptr, 'h'},
        {"data-input",  required_argument, nullptr, 'i'},
        {"is-mix",      no_argument,       nullptr, 'm'},
        {"mixpool-size",required_argument, nullptr, 'p'},
        {"savedir",     required_argument, nullptr, 's'},
        {nullptr,       0,                 nullptr, 0}
    };

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
    EventReader reader(dataInput,
        doMix ? EventReader::kDeepCopy : EventReader::kShallowCopy);

    // Offline analyzers
    AnalyzerCVE cveSame;
    cveSame.Init();
    AnalyzerQA qa;
    qa.Init();

    AnalyzerCVE cveMix;
    std::deque<Event*> mixPool;

    if (doMix) {
        cveMix.Init();
    }

    // Main event loop: handle same-event, mixed-event, and QA
    Event* evt = nullptr;
    while ((evt = reader.NextEvent())) {
        ++totalReadEvents;
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
        for (auto* e : mixPool) {
            delete e;
        }
    }

    // Write analysis outputs
    cveSame.Finish(saveDir + "/cve_single_offline.root");
    if (doMix) cveMix.Finish(saveDir + "/cve_mix_offline.root");
    qa.Finish(saveDir + "/qa_offline.root");

    if (doMix && totalReadEvents < mixPoolSize + 1) {
        std::cerr << "Warning: only " << totalReadEvents
                  << " events read, less than mixpool-size+1 ("
                  << (mixPoolSize + 1)
                  << "), no mixed-event analysis was performed.\n";
    }

    return 0;
}