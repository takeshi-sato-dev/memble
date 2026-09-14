#!/usr/bin/env python3
"""Figure S1. The composition of each leaflet against the distance from the protein.

Data: radial_species.json, 1751 frames from 1 to 8 us of the 8 us system of
Section 3.3, 206 transmembrane beads as the reference. Each lipid is given the
distance in the membrane plane to the nearest transmembrane bead, and the value
of a shell is the share of the lipids of that leaflet in that shell that are the
named species. Figure 3 of the article draws the cholesterol line of both
leaflets; this figure draws every species of each leaflet.

Panel A  the upper leaflet, which holds the sphingomyelin
Panel B  the lower leaflet, which holds the DOPS
"""
import json
import os
import sys
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
D = json.load(open(os.path.join(DATA, "radial_species.json")))
e = np.asarray(D["edges"], float)
mid = np.array([(e[i] + min(e[i + 1], 8.0)) / 2 for i in range(len(e) - 1)])
mid[-1] = 7.6
total = {k: np.asarray(v, float) for k, v in D["count"].items()}

COLOUR = {"CHOL": style.BLUE, "DLPC": style.ORANGE,
          "PSM": style.GREEN, "DOPS": style.GREEN}
MARKER = {"CHOL": "o", "DLPC": "s", "PSM": "^", "DOPS": "^"}
LABEL  = {"CHOL": "cholesterol", "DLPC": "DLPC",
          "PSM": "sphingomyelin", "DOPS": "DOPS"}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, side, title in ((axA, "upper", "upper leaflet, which holds the sphingomyelin"),
                        (axB, "lower", "lower leaflet, which holds the DOPS")):
    for sp in ("CHOL", "DLPC", "PSM", "DOPS"):
        counts = np.asarray(D["species"][sp][side], float)
        if counts.sum() == 0:
            continue
        y = 100.0 * counts / total[side]
        ax.plot(mid, y, marker=MARKER[sp], color=COLOUR[sp], mfc="white",
                mew=1.3, ms=4.5, lw=1.4, label=LABEL[sp])
        print("%-6s %-6s %s" % (side, sp, " ".join("%.1f" % v for v in y)))
    ax.set_xlabel("distance from the nearest\ntransmembrane bead (nm)")
    ax.set_xlim(0, 8.2)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_xticklabels(["0", "2", "4", "6", "> 7"])
    ax.set_title(title, loc="left", fontsize=8.5, color=style.MUTED)
    ax.legend(loc="upper right", handletextpad=0.5, labelspacing=0.5)
axA.set_ylabel("share of the lipids of that\nleaflet in that shell (mol%)")
axA.set_ylim(0, 75)
style.panel_label(axA, "A")
style.panel_label(axB, "B")
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "FigureS1_species.%s" % ext))
print("written to", OUT)
