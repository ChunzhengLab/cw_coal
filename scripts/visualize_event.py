#!/usr/bin/env python3
"""
Visualize a single event from a ROOT file using PyROOT and matplotlib,
with same canvas layouts and improved styling.

Usage:
    python visualize_event.py -i INPUT_FILE.root -e EVENT_INDEX -l LIB_PATH
"""
import argparse
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def load_event(input_file, event_index, lib_path):
    ROOT.gSystem.Load(lib_path)
    f = ROOT.TFile.Open(input_file, 'READ')
    tree = f.Get('events')
    if not tree:
        raise RuntimeError(f"Error: 'events' tree not found in {input_file}")
    evt = ROOT.Event()
    tree.SetBranchAddress('event', evt)
    tree.GetEntry(event_index)
    v_partons = evt.GetPartons()
    v_hadrons = evt.GetHadrons()
    return [v_partons[i] for i in range(v_partons.size())], [v_hadrons[i] for i in range(v_hadrons.size())]

def determine_limits(coords):
    if coords.size == 0:
        return 1.0
    return 1.5 * max(abs(coords.min()), abs(coords.max()))

def visualize_event(input_file, event_index, lib_path, output_dir):
    # Load data
    partons, hadrons = load_event(input_file, event_index, lib_path)

    # Count and print composition
    num_q = sum(p.GetBaryonNumber()>0 for p in partons)
    num_aq = sum(p.GetBaryonNumber()<0 for p in partons)
    num_baryons = sum(h.GetBaryonNumber()>0 for h in hadrons)
    num_antibaryons = sum(h.GetBaryonNumber()<0 for h in hadrons)
    num_mesons = sum(h.GetBaryonNumber()==0 for h in hadrons)
    net_pre = (num_q - num_aq) / 3.0
    net_post = num_baryons - num_antibaryons
    print(f"Quarks: {num_q}, Anti-quarks: {num_aq} (net before: {net_pre:.2f})")
    print(f"Mesons: {num_mesons}, Baryons: {num_baryons}, Anti-baryons: {num_antibaryons} (net after: {net_post})")

    # Prepare arrays
    xs = np.array([p.X() for p in partons] + [h.X() for h in hadrons])
    ys = np.array([p.Y() for p in partons] + [h.Y() for h in hadrons])
    zs = np.array([p.Z() for p in partons] + [h.Z() for h in hadrons])
    types = np.array([
        'q' if p.GetBaryonNumber()>0 else 'aq' for p in partons
    ] + [
        'b' if h.GetBaryonNumber()>0 else 'ab' if h.GetBaryonNumber()<0 else 'm'
        for h in hadrons
    ])

    # Axis limits and zoom
    xL, yL, zL = determine_limits(xs), determine_limits(ys), determine_limits(zs)
    zoom_x, zoom_y, zoom_z = 1.0, 1.0, 2.0

    # Styling: marker, color, size
    style = {
        'q':   {'marker':'.','color':'#e67c73','size': 3 * 5}, # 红
        'aq':  {'marker':'.','color':'#7baaf7','size': 3 * 5}, # 蓝
        'm':   {'marker':'o','color':'#f7cb4d','size': 5 * 5}, # 黄
        'b':   {'marker':'^','color':'#41b375','size': 7 * 5}, # 绿
        'ab':  {'marker':'v','color':'#ba67c8','size': 7 * 5}  # 紫
    }
    # Full names for legend
    names = {
        'q': 'Quark',
        'aq': 'Anti-quark',
        'm': 'Meson',
        'b': 'Baryon',
        'ab': 'Anti-baryon'
    }

    # Draw constituent connections with lighter color based on hadron type
    def draw_connections(ax, coordX, coordY):
        pid_map = {p.UniqueID(): p for p in partons}
        for h in hadrons:
            # Determine hadron type key
            if h.GetBaryonNumber()>0:
                key = 'b'
            elif h.GetBaryonNumber()<0:
                key = 'ab'
            else:
                key = 'm'
            base_color = style[key]['color']
            # Use alpha to lighten
            line_kwargs = {'color': base_color, 'alpha':0.3, 'linewidth':0.8}
            comps = [pid_map[i] for i in h.GetConstituentIDs() if i in pid_map]
            if len(comps)==3:
                pts = [(coordX(p),coordY(p)) for p in comps] + [(coordX(comps[0]),coordY(comps[0]))]
                xs_p, ys_p = zip(*pts)
                ax.plot(xs_p, ys_p, **line_kwargs)
            elif len(comps)==2:
                x0,y0 = coordX(comps[0]), coordY(comps[0])
                x1,y1 = coordX(comps[1]), coordY(comps[1])
                ax.plot([x0,x1],[y0,y1], **line_kwargs)

    # --- Figure 1: full + zoomed projections ---
    fig1, axes1 = plt.subplots(2,2,figsize=(12,12))
    ax1, ax2 = axes1[0]
    ax3, ax4 = axes1[1]

    # Full XY
    for t,p in style.items():
        mask = types==t
        ax1.scatter(xs[mask], ys[mask], s=p['size'], c=p['color'], marker=p['marker'], label=names[t])
    draw_connections(ax1, lambda p:p.X(), lambda p:p.Y())
    ax1.set(xlim=(-xL,xL), ylim=(-yL,yL), xlabel='X', ylabel='Y', title='XY Projection')
    ax1.legend(fontsize='small', loc='best')

    # Full ZY
    for t,p in style.items():
        mask = types==t
        ax2.scatter(zs[mask], ys[mask], s=p['size'], c=p['color'], marker=p['marker'], label=names[t])
    draw_connections(ax2, lambda p:p.Z(), lambda p:p.Y())
    ax2.set(xlim=(-zL,zL), ylim=(-yL,yL), xlabel='Z', ylabel='Y', title='ZY Projection')
    ax2.legend(fontsize='small', loc='best')

    # Zoomed XY
    for t,p in style.items():
        ax3.scatter(xs[types==t], ys[types==t], s=p['size'], c=p['color'], marker=p['marker'], label=names[t])
    draw_connections(ax3, lambda p:p.X(), lambda p:p.Y())
    ax3.set(xlim=(-zoom_x,zoom_x), ylim=(-zoom_y,zoom_y), xlabel='X', ylabel='Y', title='Zoomed XY')
    ax3.legend(fontsize='small', loc='best')

    # Zoomed ZY
    for t,p in style.items():
        ax4.scatter(zs[types==t], ys[types==t], s=p['size'], c=p['color'], marker=p['marker'], label=names[t])
    draw_connections(ax4, lambda p:p.Z(), lambda p:p.Y())
    ax4.set(xlim=(-zoom_z,zoom_z), ylim=(-zoom_y,zoom_y), xlabel='Z', ylabel='Y', title='Zoomed ZY')
    ax4.legend(fontsize='small', loc='best')

    fig1.tight_layout()

    filename_base = f"event_{event_index:05d}"
    fig1.savefig(os.path.join(output_dir, f"{filename_base}_projections.pdf"))

    # --- Figure 2: 1D distributions ---
    fig2, axes2 = plt.subplots(2,2,figsize=(12,8))
    axx, ayy = axes2[0]
    azz, alegend = axes2[1]
    for t,p in style.items():
        axx.hist(xs[types==t], bins=100, range=(-xL,xL), histtype='step', color=p['color'], linewidth=1, label=names[t])
    axx.set(xlabel='X', ylabel='Counts', title='X Distribution')
    axx.legend(fontsize='small', loc='best')
    for t,p in style.items():
        ayy.hist(ys[types==t], bins=100, range=(-yL,yL), histtype='step', color=p['color'], linewidth=1, label=names[t])
    ayy.set(xlabel='Y', ylabel='Counts', title='Y Distribution')
    ayy.legend(fontsize='small', loc='best')
    for t,p in style.items():
        azz.hist(zs[types==t], bins=100, range=(-zL,zL), histtype='step', color=p['color'], linewidth=1, label=names[t])
    azz.set(xlabel='Z', ylabel='Counts', title='Z Distribution')
    azz.legend(fontsize='small', loc='best')
    alegend.axis('off')
    fig2.tight_layout()

    fig2.savefig(os.path.join(output_dir, f"{filename_base}_distributions.pdf"))

    # --- Figure 3: Global Ratios ---
    fig3, axr = plt.subplots(figsize=(6,4))
    mask_h = np.isin(types, ['b','ab','m'])
    th = types[mask_h]
    cB, cAB, cM = np.sum(th=='b'), np.sum(th=='ab'), np.sum(th=='m')
    rBM = (cB + cAB)/cM if cM else 0
    rBAB = cB/cAB if cAB else 0
    axr.bar(['(Anti)Baryon/Meson', 'Baryon/Anti-Baryon'], [rBM, rBAB], color='gray', alpha=0.5)
    for i,v in enumerate([rBM, rBAB]): axr.text(i, v*1.03, f"{v:.2f}", ha='center')
    axr.set(ylabel='Value', title='Global Ratios')
    fig3.tight_layout()
    fig3.savefig(os.path.join(output_dir, f"{filename_base}_ratios.pdf"))
#    plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-data', default='test_KDTreeGlobal.root', help='ROOT file')
    parser.add_argument('-e', '--event-index', type=int, default=0, help='Event index')
    parser.add_argument('-l', '--lib', dest='lib_path', default='libcw_coal.dylib', help='Shared lib path')
    parser.add_argument('-o', '--outputdir', default='.', help='Output directory for saved plots')
    args = parser.parse_args()
    visualize_event(args.input_data, args.event_index, args.lib_path, args.outputdir)
