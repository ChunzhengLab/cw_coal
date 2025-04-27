#!/usr/bin/env python3
import argparse
import subprocess
import os
import ROOT
import matplotlib.pyplot as plt
import re
plt.rcParams.update({"font.family": "serif"})

r_values = [0.1, 0.5, 1, 1.5, 2, 3, 999.99]
r_labels = [r'$(\overline{B}+B)/M$', r'$\overline{B}/B$', r'$p/\pi^+$',
            r'$\overline{p}/p$', r'$\Lambda/p$', r'$K^+/\pi^+$', r'$\rho^+/\pi^+$']
ampt_reference = [0.220177, 1./1.74087, 0.272094, 0.431481, 0.23647, 0.348683, 0.359315]
qa_prefix = "./qa_KDTreeGlobal_r{:.2f}.root"

def get_r_from_existing_files():
    files = os.listdir(".")
    r_list = []
    pattern = re.compile(r"qa_KDTreeGlobal_r([\d.]+)\.root")
    for fname in files:
        match = pattern.match(fname)
        if match:
            try:
                r = float(match.group(1))
                r_list.append(r)
            except ValueError:
                continue
    return sorted(r_list)

def run_all(exe_path, r_list, events):
    for r in r_list:
        print(f"Running {exe_path} with r={r:.2f}")
        cmd = [
            exe_path,
            "-i", "/Users/wangchunzheng/works/Models/CWCoalProject/zpc-1.root",
            "-n", str(events),
            "-r", f"{r:.2f}"
        ]
        subprocess.run(cmd, check=True)

def analyze_all(r_list):
    ratios = [[] for _ in range(7)]  # 7 bins
    for r in r_list:
        fname = qa_prefix.format(r)
        if not os.path.exists(fname):
            print(f"Warning: {fname} not found, skipping.")
            continue
        f = ROOT.TFile(fname, "READ")
        h = f.Get("hRatio")
        if not h:
            print(f"Warning: hRatio not found in {fname}, skipping.")
            continue
        for i in range(7):
            ratios[i].append(h.GetBinContent(i + 1))
        f.Close()
    return ratios

def plot_ratios(ratios, r_list):
    colors = plt.colormaps.get_cmap('tab10').resampled(7)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True)

    left_indices = [0, 1, 2, 3]  # 左图：B相关和p相关
    right_indices = [4, 5, 6]   # 右图：Lambda和轻介子

    for i in left_indices:
        ax1.plot(r_list, ratios[i], marker='o', label=r_labels[i], color=colors(i))
        ax1.hlines(ampt_reference[i], min(r_list), max(r_list),
                   linestyles='dashed', colors=colors(i), alpha=0.8)
        ax1.text(max(r_list), ampt_reference[i], f"{ampt_reference[i]:.3f}", 
                 va='center', ha='left', fontsize=8, color=colors(i))

    for i in right_indices:
        ax2.plot(r_list, ratios[i], marker='o', label=r_labels[i], color=colors(i))
        ax2.hlines(ampt_reference[i], min(r_list), max(r_list),
                   linestyles='dashed', colors=colors(i), alpha=0.8)
        ax2.text(max(r_list), ampt_reference[i], f"{ampt_reference[i]:.3f}", 
                 va='center', ha='left', fontsize=8, color=colors(i))

    ax1.set_title("Ratios: Baryons and Protons")
    ax2.set_title("Ratios: Lambda and Mesons")

    for ax in (ax1, ax2):
        ax.set_xlabel(r"Baryon Preference $r$")
        ax.set_ylabel("Hadron Yield Ratios")
        ax.grid(True)
        ax.legend()

    plt.tight_layout()
    plt.savefig("ratios_vs_r_split.pdf")
    print("Saved figure to ratios_vs_r_split.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run and analyze CWCoal hadron ratios.")
    parser.add_argument("-p", "--process", action="store_true", help="Run cwcoal with specified r values.")
    parser.add_argument("-a", "--ana", action="store_true", help="Analyze generated ROOT files.")
    parser.add_argument("-e", "--exe", type=str, default="~/works/Models/CWCoalProject/bin/cwcoal",
                        help="Path to cwcoal executable")
    parser.add_argument("-d", "--dir", type=str, default=".",
                        help="Working directory for running and analyzing")
    parser.add_argument("-r", "--baryon-preference", nargs="+", type=float,
                        default=r_values,
                        help="List of baryon preference r values (e.g. -r 0.1 0.5 1)")
    parser.add_argument("-n", "--events", type=int, default=100,
                        help="Number of events to process")
    args = parser.parse_args()

    # Change to working directory if specified
    if args.dir:
        os.chdir(args.dir)
    exe_path = args.exe
    r_list = args.baryon_preference
    events = args.events

    if args.process:
        run_all(exe_path, r_list, events)
    if args.ana:
        r_existing = get_r_from_existing_files()
        print("Detected r values from files:", r_existing)
        data = analyze_all(r_existing)
        plot_ratios(data, r_existing)
    if not args.process and not args.ana:
        parser.print_help()