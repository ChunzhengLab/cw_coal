#!/usr/bin/env python3
import os
import argparse
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import re

# Labels for combinations
comboLabels = [
    "Baryon-Baryon", "Baryon-AntiBaryon", "Baryon-Meson",
    "AntiBaryon-AntiBaryon", "AntiBaryon-Meson", "Meson-Meson"
]

# Histogram prefixes and titles
hist_names = {
    "hCdPhiP": "Δφ Position",
    "hCdPhiM": "Δφ Momentum",
    "hSdPhiP": "Σφ Position",
    "hSdPhiM": "Σφ Momentum"
}

# Color map: one color per combo
colors = plt.colormaps.get_cmap('tab10').resampled(len(comboLabels))

def extract_hist_array(h):
    """Return bin centers and normalized contents (scaled by Nbins)."""
    nbins = h.GetNbinsX()
    counts = np.array([h.GetBinContent(i+1) for i in range(nbins)], dtype=float)
    centers = np.array([h.GetBinCenter(i+1) for i in range(nbins)], dtype=float)
    total = counts.sum()
    scale = nbins / total if total > 0 else 1.0
    counts *= scale
    return centers, counts

def plot_single_or_mix(rootfile, suffix, outname):
    """Plot single or mix histograms from one ROOT file."""
    f = ROOT.TFile(rootfile, "READ")
    keys = [k.GetName() for k in f.GetListOfKeys()]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (prefix, title) in zip(axs.flatten(), hist_names.items()):
        # find all matching hist names for this prefix and suffix
        pattern = re.compile(rf"^{prefix}.*{re.escape(suffix)}$")
        names = [name for name in keys if pattern.match(name)]
        if not names:
            print(f"No histograms found for prefix '{prefix}' with suffix '{suffix}'")
        for name in sorted(names):
            h = f.Get(name)
            x, y = extract_hist_array(h)
            ax.plot(x, y, marker='o', label=name, color=None)
        ax.set_title(title)
        ax.set_xlabel("φ")
        ax.set_ylabel("Normalized Counts × Nbins")
        ax.grid(True)
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(outname)
    print(f"Saved: {outname}")
    f.Close()

def plot_ratio(single_file, mix_file):
    """Plot ratio mix/single from two ROOT files."""
    f1 = ROOT.TFile(single_file, "READ")
    f2 = ROOT.TFile(mix_file, "READ")
    keys1 = [k.GetName() for k in f1.GetListOfKeys()]
    mix_suffix = "_MixEvt"
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (prefix, title) in zip(axs.flatten(), hist_names.items()):
        # find single-event hist names for this prefix
        names1 = sorted(n for n in keys1 if n.startswith(prefix))
        if not names1:
            print(f"No single-event histograms for prefix '{prefix}'")
        for name1 in names1:
            # determine corresponding mix-event histogram name
            name2 = name1 + mix_suffix
            if not f2.Get(name2):
                name2 = name1
            h1 = f1.Get(name1)
            h2 = f2.Get(name2)
            if not h1 or not h2:
                print(f"Missing histogram pair: '{name1}' , '{name2}'")
                continue
            x1, y1 = extract_hist_array(h1)
            _, y2 = extract_hist_array(h2)
            # Compute single/mix ratio instead of mix/single
            ratio = np.divide(y1, y2, out=np.zeros_like(y1), where=y2!=0)
            ax.plot(x1, ratio, marker='o', label=name1)
        ax.set_title(f"{title} Ratio")
        ax.set_xlabel("φ")
        ax.set_ylabel("Single/Mix")
        ax.grid(True)
        ax.legend(fontsize=8)
    plt.tight_layout()
    outname = "CVE_phi_ratio.pdf"
    plt.savefig(outname)
    print(f"Saved: {outname}")
    f1.Close()
    f2.Close()

def main():
    parser = argparse.ArgumentParser(description="Plot CVE φ distributions.")
    parser.add_argument("-s", "--single", type=str, metavar="FILE",
                        help="Single-event ROOT file")
    parser.add_argument("-m", "--mix", type=str, metavar="FILE",
                        help="Mix-event ROOT file")
    args = parser.parse_args()

    if args.single and not args.mix:
        plot_single_or_mix(args.single, suffix="", outname="CVE_phi_single.pdf")
    elif args.mix and not args.single:
        plot_single_or_mix(args.mix, suffix="_MixEvt", outname="CVE_phi_mix.pdf")
    elif args.single and args.mix:
        plot_ratio(args.single, args.mix)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()