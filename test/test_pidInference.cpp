#include "core/PIDInference.h"
#include "core/PhysicsConstants.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <map>

// Helper to run a single test case
bool runTest(const std::vector<int>& quarks, double mass,
             int expected, const std::string& desc) {
    int result = PIDInference::InferPID(quarks, mass);
    bool pass = (result == expected);
    std::cout << (pass ? "[PASS]" : "[FAIL]") << " " << desc
              << " => Got " << result << ", expected " << expected << "\n";
    return pass;
}

int main() {
    int failures = 0;

    // PDG code to particle name mapping for statistics output
    auto pdgName = [](int pdg) {
        switch (pdg) {
            case 211:   return "pi+";
            case 213:   return "rho+";
            case 111:   return "pi0";
            case 113:   return "rho0";
            case 221:   return "eta";
            case 223:   return "omega";
            case 333:   return "phi";
            case 443:   return "J/psi";
            default:    return "unknown";
        }
    };

    // Meson tests
    // s + s̄ -> φ (333)
    failures += !runTest({3, -3}, 1.019, 333, "s s̄ -> φ");

    // c + c̄ -> J/ψ (443)
    failures += !runTest({4, -4}, 3.097, 443, "c c̄ -> J/ψ");

    // Baryon tests
    // u + u + d -> proton (2212)
    failures += !runTest({2, 2, 1}, 0.938, 2212, "u u d -> p");

    // ū + ū + d̄ -> anti-proton (-2212)
    failures += !runTest({-2, -2, -1}, 0.938, -2212, "ū ū d̄ -> anti-p");

    // Spin test for mesons
    // rnd < 1/(1+vpratio) => pseudoscalar (0)
    {
        double vpr = PhysicsConstants::kMesonVectorToPseudoscalarRatio;
        double thresh = 1.0 / (1.0 + vpr);
        bool pass = (PIDInference::InferMesonSpin(vpr, thresh / 2.0) == 0);
        std::cout << (pass ? "[PASS]" : "[FAIL]") << " InferMesonSpin pseudoscalar\n";
        failures += !pass;
    }
    // rnd > 1/(1+vpratio) => vector (1)
    {
        double vpr = PhysicsConstants::kMesonVectorToPseudoscalarRatio;
        double thresh = 1.0 / (1.0 + vpr);
        bool pass = (PIDInference::InferMesonSpin(vpr, thresh + (1.0 - thresh)/2.0) == 1);
        std::cout << (pass ? "[PASS]" : "[FAIL]") << " InferMesonSpin vector\n";
        failures += !pass;
    }

    // Statistical test for random outcomes
    {
        const int N = 10000;
        std::map<int,int> counts;
        // Test u + d̄ random distribution
        for(int i = 0; i < N; ++i) {
            int pdg = PIDInference::InferPID({2, -1}, 0.140);
            counts[pdg]++;
        }
        std::cout << "Statistics for u + d̄ over " << N << " trials:\n";
        for (auto& kv : counts) {
            std::cout << "  PDG " << kv.first << " (" << pdgName(kv.first) << "): "
                      << (100.0 * kv.second / N) << "%\n";
        }
    }

    {
        const int N = 10000;
        std::map<int,int> counts;
        // Test uū diagonal distribution
        for(int i = 0; i < N; ++i) {
            int pdg = PIDInference::InferPID({2, -2}, 0.135);
            counts[pdg]++;
        }
        std::cout << "Statistics for uū diagonal over " << N << " trials:\n";
        for (auto& kv : counts) {
            std::cout << "  PDG " << kv.first << " (" << pdgName(kv.first) << "): "
                      << (100.0 * kv.second / N) << "%\n";
        }
    }

    if (failures == 0) {
        std::cout << "All PIDInference tests passed!\n";
        return 0;
    } else {
        std::cout << failures << " PIDInference tests failed.\n";
        return 1;
    }
}