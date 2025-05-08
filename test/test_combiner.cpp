#include "core/Event.h"
#include "io/EventWriter.h"
#include "Combiners.h"
#include "core/PhysicsConstants.h"
#include "ana/AnalyzerQA.h"
#include "io/EventRandomGen.h"
#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <string>

// Path to parton histogram file: CLI override > env var > installed DATA_INSTALL_DIR
static const std::string defaultHistPath = std::string(DATA_INSTALL_DIR) + "/dist_parton_afART.root";
static const char* kPartonHistFile =
    std::getenv("CW_COAL_PARTON_HIST") ?
    std::getenv("CW_COAL_PARTON_HIST") :
    defaultHistPath.c_str();

// Lambda to deep-copy an Event (clone all Partons)
auto cloneEvent = [](const Event& src) {
    Event clone;
    for (auto* p : src.GetPartons()) {
        clone.AddParton(new Parton(*p));
    }
    return clone;
};

// Run a single test: receives a reference to base Event, deep-copies partons, applies combiner, and prints stats
void RunTest(const std::string& label, CombinerBase& combiner,
             Event& event) {
    auto localPartons = event.GetPartons();

    // Time the combine operation
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Hadron*> hadrons = combiner.Combine(localPartons);
    for (auto* h : hadrons) {
        event.AddHadron(h);
    }
    auto end   = std::chrono::high_resolution_clock::now();
    double ms  = std::chrono::duration<double, std::milli>(end - start).count();

    // Print statistics
    std::cout << label << ":\n";
    std::cout << localPartons.size() << " partons formed into " << hadrons.size() << " hadrons in " << ms << " ms\n";
    size_t unused = std::count_if(localPartons.begin(), localPartons.end(),
                                  [](Parton* p){ return !p->IsUsed(); });
    std::cout << "  Unused partons: " << unused << "\n";
    size_t nBaryons = 0, nMesons = 0;
    for (auto* h : hadrons) {
        if (h->GetBaryonNumber() == 0) ++nMesons;
        else ++nBaryons;
    }
    std::cout << "  Baryons: " << nBaryons << ", Mesons: " << nMesons << "\n";
    size_t totalUsed = 3 * nBaryons + 2 * nMesons + unused;
    if (totalUsed != localPartons.size()) {
        std::cerr << " Consistency check failed: "
                  << "3 * baryons + 2 * mesons + unused = " << totalUsed
                  << " != " << localPartons.size() << "\n";
    } else {
        std::cout << "Consistency check passed.\n";
    }
}

int main() {
    // Prepare combiner instances
    // BruteForceGlobal bfGlobal;
    // BruteForceGreedy bfGreedy;
    KDTreeGlobal    kdGlobal;
    KDTreeGreedy   kdGreedy;
    // BruteForceDualGreedy bfDualGreedy;
    // KDTreeDualGreedy    kdDualGreedy;

    // Create an EventRandomGen to produce parton lists
    EventRandomGen eventGen(kPartonHistFile);

    // Define test cases: label and combiner reference
    struct TestCase { const char* label; CombinerBase& combiner; };
    TestCase tests[] = {
      // { "BruteForceGlobal", bfGlobal },
      // { "BruteForceGreedy", bfGreedy },
      { "KDTreeGlobal",    kdGlobal },
      { "KDTreeGreedy",   kdGreedy }
      // , { "BruteForceDualGreedy", bfDualGreedy }
      // , { "KDTreeDualGreedy",    kdDualGreedy }
    };

    // Number of events to simulate per algorithm
    const int nEvents = 5;
    // Pre-generate base events so all tests use identical Parton sets
    std::vector<Event> baseEvents;
    baseEvents.reserve(nEvents);
    for (int ie = 0; ie < nEvents; ++ie) {
        int nParts = static_cast<int>(PhysicsConstants::GetMultiplicityHistogram().GetRandom());
        std::cout << "Generating base event " << (ie+1)
                  << " of " << nEvents << ": " << nParts << " partons" << std::endl;
        Event evt;
        eventGen.GenerateEvent(evt, nParts, /*sumBaryonNumber=*/0);
        baseEvents.push_back(evt);
    }

    // Run each test: each writer is locally created per test
    for (auto& tc : tests) {
        // Initialize QA analyzer for this combiner
        AnalyzerQA qa;
        qa.Init();
        EventWriter writer(std::string("test_") + tc.label + ".root");
        for (int ie = 0; ie < nEvents; ++ie) {
            std::cout << ">>>>>================================="<<std::endl;
            auto& evt = baseEvents[ie];
            Event event = cloneEvent(evt);
            std::cout << "Using base event " << (ie+1)
                      << " of " << nEvents << " with " << evt.GetPartons().size()
                      << " partons for test \"" << tc.label << "\"" << std::endl;
            RunTest(tc.label, tc.combiner, event);
            writer.WriteEvent(&event);
            qa.Process(event);
            event.Reset();
            std::cout << "=================================<<<<<"<<std::endl << std::endl;
        }
        // After all events, write out QA histograms
        qa.Finish(std::string("qa_") + tc.label + ".root");
    }

    for (auto &evt : baseEvents) {
      evt.Reset();
    }

    return 0;
}
