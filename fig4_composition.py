#!/usr/bin/env python3
"""Figure 4. What a run does with the leaflet composition it was built with.

Data: CG_M3/leaflet_analysis/dipc_8us.json, written by leaflet_series.py.
4001 frames at 2 ns over 8 us, 2905 lipid molecules. The build gave the upper
leaflet 501 cholesterol molecules of 969 (system.top), and that number, not the
first frame of the production run, is the reference every trace is drawn against.

Panel A  the whole 8 us, all four species
Panel B  the first 500 ns, cholesterol alone, with the windows of Section 3.3
"""
import json, os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as style

# Where the figure is written, and where its data is read from. Both are
# directories, both default to this file's own directory, and both can be set
# from the environment so that the scripts run anywhere the deposit is unpacked.
HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.environ.get("MEMBLE_FIG_OUT", HERE)
DATA = os.environ.get("MEMBLE_DATA", HERE)

style.use()
D = json.load(open(os.path.join(DATA, "dipc_8us.json")))
t = np.asarray(D["time_ps"], float) / 1000.0          # ns
counts = D["counts"]
BUILT = {"CHOL": 501, "DLPC": 500, "PSM": 500, "DOPS": 0}
SETTLED_FROM = 1000.0


def running(y, win_ns):
    """Centred running mean that shortens the window at the two ends, so the
    first and last frames are not pulled toward zero by padding."""
    n = max(1, int(round(win_ns / (t[1] - t[0]))))
    h = n // 2
    c = np.concatenate(([0.0], np.cumsum(y)))
    out = np.empty_like(y)
    for i in range(len(y)):
        a = max(0, i - h); b = min(len(y), i + h + 1)
        out[i] = (c[b] - c[a]) / (b - a)
    return out


fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 2.7))

# ---- A: every species over 8 us
for sp, col, lab in (("DLPC", style.GREY, "DLPC"),
                     ("PSM", style.GREY, "PSM"),
                     ("DOPS", style.GREY, "DOPS")):
    y = np.asarray(counts[sp]["upper"], float) - BUILT[sp]
    axA.plot(t / 1000.0, y, color=col, lw=0.9, alpha=0.9, zorder=2)
y = np.asarray(counts["CHOL"]["upper"], float) - BUILT["CHOL"]
axA.plot(t / 1000.0, y, color=style.BLUE, lw=0.4, alpha=0.22, zorder=3)
axA.plot(t / 1000.0, running(y, 50.0), color=style.BLUE, lw=1.6, zorder=4)

settled = y[t >= SETTLED_FROM].mean()
axA.axhline(settled, color=style.INK, lw=0.8, ls=(0, (5, 2)), zorder=5)
axA.axhline(0, color=style.MUTED, lw=0.6, zorder=1)
axA.text(8.12, settled, "settled\n+%.1f" % settled, ha="left", va="center",
         fontsize=7.5, color=style.INK, clip_on=False)
axA.text(7.9, 2.0, "DLPC, PSM, DOPS", ha="right", va="bottom",
         fontsize=7.5, color=style.GREY)
axA.text(4.2, 11.0, "cholesterol", ha="left", va="bottom",
         fontsize=7.5, color=style.BLUE)
axA.set_xlabel("time (μs)")
axA.set_ylabel("molecules gained by the\nupper leaflet since the build")
axA.set_xlim(0, 8); axA.set_ylim(-6, 58)
style.panel_label(axA, "A")

# ---- B: the first 500 ns, cholesterol alone
m = t <= 500
axB.plot(t[m], y[m], color=style.BLUE, lw=0.4, alpha=0.22)
axB.plot(t[m], running(y, 18.0)[m], color=style.BLUE, lw=1.6)
axB.axhline(settled, color=style.INK, lw=0.8, ls=(0, (5, 2)))
axB.axhline(0, color=style.MUTED, lw=0.6)
axB.axvspan(100, 250, color=style.BLUE, alpha=0.07, lw=0)
axB.text(175, 3.0, "100–250 ns", ha="center", va="bottom", fontsize=7,
         color=style.MUTED)
axB.axvline(9.5, color=style.MUTED, lw=0.8, ls=(0, (1.5, 1.5)))
axB.text(14, 54, "9.5 ns, the whole\nstandard equilibration", ha="left", va="top",
         fontsize=7, color=style.MUTED)
axB.set_xlabel("time (ns)")
axB.set_ylabel("cholesterol molecules gained\nby the upper leaflet")
axB.set_xlim(0, 500); axB.set_ylim(-6, 58)
style.panel_label(axB, "B")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "Figure4_composition.%s" % ext))
print("settled %.1f molecules, %.2f%% of 969" % (settled, 100 * (501 + settled) / 969))
print("windows:", [(a, b, round(float(y[(t >= a) & (t < b)].mean()), 1))
                   for a, b in ((0, 50), (50, 100), (100, 250))])
