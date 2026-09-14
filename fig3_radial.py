#!/usr/bin/env python3
"""Figure 3. Cholesterol against the distance from the protein.

Data: CG_M3/leaflet_analysis/chol_radial.json, 1751 frames, 206 transmembrane
beads as the reference. Each shell holds the lipids whose head bead lies within
that distance of the nearest transmembrane bead, and the value is the share of
those lipids that are cholesterol.

The upper leaflet separates into two phases over the run. The protein sits in
the cholesterol-poor one. The lower leaflet shows no such gradient.
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
D = json.load(open(os.path.join(DATA, "chol_radial.json")))
e = np.asarray(D["edges"], float)
mid = np.array([(e[i] + min(e[i + 1], 8.0)) / 2 for i in range(len(e) - 1)])
mid[-1] = 7.6

fig, ax = plt.subplots(figsize=(3.6, 2.9))
for side, col, mk, lab in (("upper", style.BLUE, "o", "upper leaflet\n(cholesterol, DLPC, sphingomyelin)"),
                           ("lower", style.ORANGE, "s", "lower leaflet\n(cholesterol, DLPC, DOPS)")):
    y = 100.0 * np.asarray(D["chol"][side], float) / np.asarray(D["count"][side], float)
    ax.plot(mid, y, marker=mk, color=col, mfc="white", mew=1.3, ms=4.5, lw=1.4,
            label=lab)
    print(side, " ".join("%.1f" % v for v in y))

ax.annotate("22.0%", xy=(mid[0], 22.0), xytext=(mid[0] + 0.35, 19.6),
            fontsize=7.5, color=style.BLUE, ha="left",
            arrowprops=dict(arrowstyle="-", lw=0.6, color=style.BLUE))
ax.annotate("39.0%", xy=(mid[-1], 39.0), xytext=(mid[-1] - 0.35, 41.2),
            fontsize=7.5, color=style.BLUE, ha="right",
            arrowprops=dict(arrowstyle="-", lw=0.6, color=style.BLUE))
ax.set_xlabel("distance from the nearest\ntransmembrane bead (nm)")
ax.set_ylabel("cholesterol of the lipids\nin that shell (mol%)")
ax.set_xlim(0, 8.2); ax.set_ylim(17, 45)
ax.set_xticks([0, 2, 4, 6, 8])
ax.set_xticklabels(["0", "2", "4", "6", "> 7"])
ax.legend(loc="lower right", handletextpad=0.5, labelspacing=0.6)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "Figure3_radial.%s" % ext))
