#include "core/Event.h"
#include "io/EventWriter.h"
#include "Combiners.h"
#include "core/PhysicsConstants.h"
#include "ana/AnalyzerQA.h"
#include <iostream>
#include <string>
#include <memory>
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstdlib> 
#include <string>

// Path to parton histogram file: CLI override > env var > installed DATA_INSTALL_DIR
static const std::string defaultHistPath = std::string(DATA_INSTALL_PREFIX) + "/" + std::string(DATA_INSTALL_DATADIR) + "/cw_coal/refdata" + "/dist_parton_afART.root";
static const char* kPartonHistFile =
    std::getenv("CW_COAL_PARTON_HIST") ?
    std::getenv("CW_COAL_PARTON_HIST") :
    defaultHistPath.c_str();

std::vector<Parton*> GenerateRandomPartons(size_t initialCount, int sumBaryonNumber) {
  std::vector<Parton*> partons;
  int baryonSum3 = 0;  // 以三分之一单位整数累加

  for (size_t i = 0; i < initialCount; ++i) {
    Parton* p = Parton::RandomFromHists(kPartonHistFile);
    int bn3 = std::round(p->GetBaryonNumber() * 3); // ±1
    baryonSum3 += bn3;
    partons.push_back(p);
  }

  int targetSum3 = sumBaryonNumber * 3;

  // 自动补充直到达到目标总和
  while (baryonSum3 != targetSum3) {
    Parton* p = Parton::RandomFromHists(kPartonHistFile);
    int bn3 = std::round(p->GetBaryonNumber() * 3); // ±1

    // 只接受能让当前和靠近目标的粒子
    if ((baryonSum3 < targetSum3 && bn3 > 0) ||
        (baryonSum3 > targetSum3 && bn3 < 0)) {
      baryonSum3 += bn3;
      partons.push_back(p);
    } else {
      delete p; // ❌ 不需要的粒子直接丢掉
    }
  }

  return partons;
}

// Helper to clone partons (deep copy)
std::vector<Parton*> ClonePartons(const std::vector<Parton*>& original) {
    std::vector<Parton*> copy;
    for (const auto* p : original) {
        copy.push_back(new Parton(*p));
    }
    return copy;
}

// Run a single test
void RunTest(const std::string& label, CombinerBase& combiner, const std::vector<Parton*>& partons, EventWriter& writer) {
    std::vector<Parton*> localPartons = ClonePartons(partons);

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Hadron*> hadrons = combiner.Combine(localPartons);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> duration = end - start;

    Event event;
    for (auto* p : localPartons) event.AddParton(p);
    for (auto* h : hadrons) event.AddHadron(h);

    writer.WriteEvent(&event);
    std::cout << label << ":\n";
    std::cout << "  Formed " << hadrons.size() << " hadrons in " << duration.count() << " ms\n";
    std::cout << "  Original partons: " << partons.size() << "\n";
    size_t unused = std::count_if(localPartons.begin(), localPartons.end(), [](Parton* p) { return !p->IsUsed(); });
    std::cout << "  Unused partons: " << unused << "\n";

    size_t nBaryons = 0, nMesons = 0;
    for (const auto* h : hadrons) {
        int b = h->GetBaryonNumber();
        if (b == 0) ++nMesons;
        else ++nBaryons;
    }
    std::cout << "  Baryons: " << nBaryons << ", Mesons: " << nMesons << "\n";
    size_t totalUsed = 3 * nBaryons + 2 * nMesons + unused;
    if (totalUsed != partons.size()) {
        std::cerr << " Consistency check failed: "
                  << "3*baryons + 2*mesons + unused = " << totalUsed
                  << " != " << partons.size() << "\n";
    } else {
        std::cout << "  ✅ Consistency check passed.\n";
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

    // Run each test: each writer is locally created per test
    for (auto& tc : tests) {
        // Initialize QA analyzer for this combiner
        AnalyzerQA qa;
        qa.Init();
        EventWriter writer(std::string("test_") + tc.label + ".root");
        for (int ie = 0; ie < nEvents; ++ie) {
            // Generate fresh partons for each event
            // Sample multiplicity from predefined histogram
            int nParts = static_cast<int>(PhysicsConstants::GetMultiplicityHistogram().GetRandom());
            std::cout << "Generating " << nParts 
                      << " partons for test \"" << tc.label 
                      << "\", event " << (ie+1) 
                      << " of " << nEvents << std::endl;
            auto basePartons = GenerateRandomPartons(nParts, 0);
            RunTest(tc.label, tc.combiner, basePartons, writer);
            // Feed the produced event into QA histograms
            // Note: RunTest creates a local Event 'event'; replicate its filling here
            {
                // Clone local event structure for QA (after RunTest)
                Event event;
                for (auto* p : basePartons) event.AddParton(p);
                auto hadrons = tc.combiner.Combine(ClonePartons(basePartons));
                for (auto* h : hadrons) event.AddHadron(h);
                qa.Process(event);
                // clean up hadrons from temporary combine
                for (auto* h : hadrons) delete h;
            }
            // Clean up partons
            for (auto* p : basePartons) delete p;
        }
        // After all events, write out QA histograms
        qa.Finish(std::string("qa_") + tc.label + ".root");
    }

    return 0;
}
