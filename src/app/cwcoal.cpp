#include "core/Particle.h"
#include "core/Event.h"
#include "io/EventReaderAMPT.h"
#include "io/EventRandomGen.h"
#include "io/EventWriter.h"
#include "ana/AnalyzerQA.h"
#include "ana/AnalyzerCVE.h"
#include "Combiners.h"
#include <TROOT.h>
#include <TApplication.h>
#include <getopt.h>
#include <iostream>
#include <memory>

// Variables for input and output file management
std::string dataInput;
std::string dataOutput;
std::string algorithm = "KDTreeGlobal";
std::string workDir = "."; // Output directory for all files
int nEvents = 10;
int sumBn = 0;
double baryonPreference = 1.0;

static void PrintUsage() {
    std::cout << "Usage: cwcoal [options]\n"
              << "  -h, --help               Show this help message\n"
              << "  -i, --data-input <file>  AMPT input ROOT file or list (if omitted, random generation mode)\n"
              << "  -o, --data-output <file> Output ROOT file for hadrons (default: output.root)\n"
              << "  -a, --algorithm <name>   Combiner algorithm: KDTreeGlobal, KDTreeGreedy, BruteForceGlobal, BruteForceGreedy (default: KDTreeGlobal)\n"
              << "  -n, --events <N>         Number of events to process/generate (default: 10)\n"
              << "  -b, --bn <B>             Target total baryon number per event (default: 0)\n"
              << "  -w, --workdir <dir>      Output directory for all files (default: current working directory)\n"
              << "  -r, --baryon-preference <R>  Baryon preference factor (default: 1.0)\n";
}

int main(int argc, char** argv) {
    const struct option longOpts[] = {
        {"help",      no_argument,       nullptr, 'h'},
        {"data-input", required_argument, nullptr, 'i'},
        {"data-output", required_argument, nullptr, 'o'},
        {"algorithm", required_argument, nullptr, 'a'},
        {"events",    required_argument, nullptr, 'n'},
        {"bn",        required_argument, nullptr, 'b'},
        {"workdir",   required_argument, nullptr, 'w'},
        {"baryon-preference", required_argument, nullptr, 'r'},
        {nullptr,     0,                 nullptr,  0 }
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "hi:o:a:n:b:w:r:", longOpts, nullptr)) != -1) {
        switch (opt) {
            case 'h': PrintUsage(); return 0;
            case 'i': dataInput = optarg; break;
            case 'o': dataOutput = optarg; break;
            case 'a': algorithm = optarg; break;
            case 'n': nEvents = std::stoi(optarg); break;
            case 'b': sumBn = std::stoi(optarg); break;
            case 'w': workDir = optarg; break;
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

    std::unique_ptr<EventWriter> writer;
    if (!dataOutput.empty()) {
        writer = std::make_unique<EventWriter>(workDir + "/" + dataOutput);
    }
    AnalyzerQA qa;
    AnalyzerCVE cve;
    qa.Init();
    cve.Init();

    if (!dataInput.empty()) {
        EventReaderAMPT reader(dataInput);
        // Read and process events from AMPT
        for (int ie = 0; ie < nEvents; ++ie) {
            Event evt;
            if (!reader.NextEvent(evt)) {
                std::cerr << "No more events available in input file\n";
                break;
            }
            // Extract partons, combine into hadrons
            auto partons = evt.GetPartons();
            auto hadrons = combiner->Combine(partons);
            for (auto* h : hadrons) evt.AddHadron(h);
            if (writer) writer->WriteEvent(&evt);
            qa.Process(evt);
            cve.Process(evt);
        }
    } else {
        // Random generation mode
        EventRandomGen gen;
        for (int ie = 0; ie < nEvents; ++ie) {
            Event evt = gen.GenerateEvent(-1, 0); // -1 uses default multiplicity sampling
            // Combine
            auto partons = evt.GetPartons();
            auto hadrons = combiner->Combine(partons);
            for (auto* h : hadrons) evt.AddHadron(h);
            if (writer) writer->WriteEvent(&evt);
            qa.Process(evt);
            cve.Process(evt);
        }
    }

    qa.Finish(workDir + "/qa_" + algorithm + ".root");
    cve.Finish(workDir + "/cve_" + algorithm + ".root");
    if (writer) writer->Close();

    return 0;
}
