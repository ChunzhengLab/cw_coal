/// check_root.C
/// ------------
///  用 ROOT 直接测试  test_output.root  结构与反序列化有效性
///
///  root -l -b -q 'check_root.C("test_output.root", 3)'
///                     └──────── ② 可选：最多打印多少个事件的详细内容

#include "TFile.h"
#include "TTree.h"
#include <iostream>

// 如果头文件已通过字典编译，下面三行可省略；
// 若 ROOT 找不到类，请把 include 路径改成你的实际位置。
#include "../include/core/Event.h"
#include "../include/core/Particle.h"

void check_root(const char* fname = "test_output.root", int maxPrint = 5)
{
    //── 打开文件
    TFile f(fname, "READ");
    if (f.IsZombie()) {
        std::cerr << "❌  无法打开文件  " << fname << std::endl;
        return;
    }

    //── 取出 events 树
    TTree* t = dynamic_cast<TTree*>(f.Get("events"));
    if (!t) {
        std::cerr << "❌  文件中没有名为 \"events\" 的 TTree\n";
        return;
    }

    //── 绑定分支
    Event* evt = nullptr;
    t->SetBranchAddress("event", &evt);

    Long64_t nEvt = t->GetEntries();
    std::cout << "📊  事件条目数: " << nEvt << "\n";

    Long64_t totalPartons  = 0;
    Long64_t totalHadrons  = 0;
    int       badPartons   = 0;
    int       badHadrons   = 0;

    for (Long64_t i = 0; i < nEvt; ++i) {
        t->GetEntry(i);

        const auto& partons = evt->GetPartons();
        const auto& hadrons = evt->GetHadrons();

        totalPartons += partons.size();
        totalHadrons += hadrons.size();

        // 打印前 maxPrint 个事件的概要
        if (i < maxPrint) {
            std::cout << "─────────────────────────────────────────────\n";
            std::cout << "Event " << i
                      << "   Partons: " << partons.size()
                      << "   Hadrons: " << hadrons.size()
                      << "   ReactionPlane: " << evt->GetReactionPlane()
                      << "\n";

            if (!partons.empty() && partons.at(0)) {
                const Parton* p0 = partons.at(0);
                std::cout << "   • First Parton  XYZ=("
                          << p0->X() << "," << p0->Y() << "," << p0->Z() << ")  "
                          << "B=" << p0->GetBaryonNumber() << "\n";
            }
            if (!hadrons.empty() && hadrons.at(0)) {
                const Hadron* h0 = hadrons.at(0);
                std::cout << "   • First Hadron  XYZ=("
                          << h0->X() << "," << h0->Y() << "," << h0->Z() << ")  "
                          << "B=" << h0->GetBaryonNumber()
                          << "   Formation=" << h0->GetFormationDistance() << "\n";
            }
        }

        // 简单合法性检查
        for (auto* p : partons)
            if (!p) ++badPartons;
        for (auto* h : hadrons)
            if (!h) ++badHadrons;
    }

    std::cout << "─────────────────────────────────────────────\n";
    std::cout << "✅  总 Parton 数: " << totalPartons
              << "    (空指针: " << badPartons << ")\n";
    std::cout << "✅  总 Hadron 数: " << totalHadrons
              << "    (空指针: " << badHadrons << ")\n";
    std::cout << "🍀  完成检查\n";
}