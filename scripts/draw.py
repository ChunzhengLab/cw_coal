import ROOT
import matplotlib.pyplot as plt

ROOT.gSystem.Load("libcw_coal.dylib")
file = ROOT.TFile("test_output.root")
event_tree = file.Get("events")

print("📦 ROOT file opened: ", file.GetName())

parton_dict = {}
all_partons = []
all_hadrons = []

failed_entries = 0
for i, entry in enumerate(event_tree):
    try:
        evt = entry.event
        partons = evt.GetPartons()
        if partons is not None:
            for i in range(partons.size()):
                p = partons.at(i)
                if p:
                    parton_dict[p.UniqueID()] = (p.X(), p.Y(), p.Z(), p.GetBaryonNumber())
                    all_partons.append(p)
        hadrons = evt.GetHadrons()
        if hadrons is not None:
            for i in range(hadrons.size()):
                h = hadrons.at(i)
                if h:
                    all_hadrons.append(h)
    except Exception as e:
        print(f"⚠️ Failed to read event {i}: {e}")
        failed_entries += 1

if failed_entries > 0:
    print(f"❌ Total failed events: {failed_entries}")

print(f"🔬 Partons: {len(all_partons)}")
print(f"🔬 Hadrons: {len(all_hadrons)}")

n_meson = 0
n_baryon = 0
n_antibaryon = 0
n_missing_constituents = 0
n_total_constituents = 0

for h in all_hadrons:
    b = h.GetBaryonNumber()
    if b == 0:
        n_meson += 1
    elif b > 0:
        n_baryon += 1
    else:
        n_antibaryon += 1

    ids = h.GetConstituentIDs()
    n_total_constituents += ids.size()
    for i in range(ids.size()):
        if ids[i] not in parton_dict:
            n_missing_constituents += 1

print(f"  ✅ Mesons: {n_meson}")
print(f"  ✅ Baryons: {n_baryon}")
print(f"  ✅ Anti-Baryons: {n_antibaryon}")
print(f"  🔗 Total Constituent IDs: {n_total_constituents}")
print(f"  ⚠️ Missing Constituent Partons: {n_missing_constituents}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

for uid, (x, y, z, b) in parton_dict.items():
    color = 'red' if b > 0 else 'blue'
    ax1.plot(x, y, 'o', color=color, markersize=4, alpha=0.6)
    ax2.plot(z, y, 'o', color=color, markersize=4, alpha=0.6)

for h in all_hadrons:
    x, y, z = h.X(), h.Y(), h.Z()
    b = h.GetBaryonNumber()
    ids = h.GetConstituentIDs()

    h_color = 'green' if b == 0 else 'orange' if b > 0 else 'purple'

    ax1.plot(x, y, '*', color=h_color, markersize=10)
    ax2.plot(z, y, '*', color=h_color, markersize=10)

    points = []
    for i in range(ids.size()):
        pid = ids[i]
        if pid in parton_dict:
            px, py, pz, _ = parton_dict[pid]
            points.append((px, py, pz))

    if len(points) == 2:
        x0, y0, z0 = points[0]
        x1, y1, z1 = points[1]
        ax1.plot([x0, x1], [y0, y1], '-', color='gray', linewidth=1)
        ax2.plot([z0, z1], [y0, y1], '-', color='gray', linewidth=1)
    elif len(points) == 3:
        for i in range(3):
            j = (i + 1) % 3
            x0, y0, z0 = points[i]
            x1, y1, z1 = points[j]
            ax1.plot([x0, x1], [y0, y1], '-', color='gray', linewidth=1)
            ax2.plot([z0, z1], [y0, y1], '-', color='gray', linewidth=1)

ax1.set_title("XY Projection")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.grid(True)

ax2.set_title("ZY Projection")
ax2.set_xlabel("Z")
ax2.set_ylabel("Y")
ax2.grid(True)

plt.tight_layout()
plt.show()