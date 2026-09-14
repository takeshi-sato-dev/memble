#!/usr/bin/env python3
"""Figure 1. What MEMBLE does, and what it refuses.

A schematic. The left column is the order of the steps; the right column names
the quantity each step leaves behind. The eight checks are the ones
verify_system.py applies to the packed system, and run.sh is written only when
all eight hold.
"""
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as style

# Where the figure is written, and where its data is read from. Both are
# directories, both default to this file's own directory, and both can be set
# from the environment so that the scripts run anywhere the deposit is unpacked.
HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.environ.get("MEMBLE_FIG_OUT", HERE)
DATA = os.environ.get("MEMBLE_DATA", HERE)

style.use()

STEPS = [
    ("all-atom structure, composition,\nbox, the transmembrane range", "input"),
    ("orient the transmembrane range\nalong z", "tilt that remains"),
    ("coarse-grain with martinize2", "the secondary structure string"),
    ("hand COBY the number of molecules\nfor each leaflet", "coby.log keeps the seed"),
    ("solvate and neutralise", "water depth on each side"),
    ("measure eight properties of the\npacked system", "memble_report.json"),
    ("write run.sh only when all eight hold", "the run files"),
]
CHECKS = ["secondary structure", "composition", "protein placement", "charge",
          "periodic image", "water layer", "overlap", "leaflet area"]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.2, 3.9),
                               gridspec_kw={"width_ratios": [1.55, 1.0]})
for ax in (axL, axR):
    ax.set_axis_off()

n = len(STEPS)
h, gap = 0.92, 0.34
for i, (txt, note) in enumerate(STEPS):
    y = n - i * (h + gap) / 1.0
    face = "#eef3fa" if i in (5, 6) else "#f6f7f9"
    edge = style.BLUE if i in (5, 6) else style.MUTED
    axL.add_patch(FancyBboxPatch((0.03, y - h / 2), 0.66, h,
                                 boxstyle="round,pad=0.012,rounding_size=0.06",
                                 lw=1.0, ec=edge, fc=face, mutation_aspect=0.6))
    axL.text(0.36, y, txt, ha="center", va="center", fontsize=7.6,
             color=style.INK)
    axL.text(0.73, y, note, ha="left", va="center", fontsize=7,
             color=style.MUTED, style="italic")
    if i < n - 1:
        axL.add_patch(FancyArrowPatch((0.36, y - h / 2), (0.36, y - h / 2 - gap + 0.03),
                                      arrowstyle="-|>", mutation_scale=8,
                                      lw=0.9, color=style.MUTED))
axL.set_xlim(0, 1.30); axL.set_ylim(n - (n - 1) * (h + gap) - h, n + h)
style.panel_label(axL, "A", dx=0.0, dy=0.99)

axR.text(0.0, 0.965, "the eight properties measured on\nthe packed system",
         fontsize=8, color=style.INK, va="top", ha="left")
for k, c in enumerate(CHECKS):
    yy = 0.845 - k * 0.072
    axR.plot([0.035], [yy + 0.012], marker="o", ms=4, mfc="white",
             mec=style.BLUE, mew=1.3)
    axR.text(0.085, yy, "%d. %s" % (k + 1, c), fontsize=7.6, color=style.INK,
             va="center", ha="left")
axR.add_patch(FancyBboxPatch((0.0, 0.055), 0.98, 0.185,
                             boxstyle="round,pad=0.012,rounding_size=0.03",
                             lw=1.0, ec=style.ORANGE, fc="#fdf3ec"))
axR.text(0.49, 0.148,
         "one failure and run.sh is not written.\nThe report names the measured value,\n"
         "the expected value, and what to change.",
         fontsize=7.4, color=style.INK, va="center", ha="center")
axR.set_xlim(0, 1); axR.set_ylim(0, 1)
style.panel_label(axR, "B", dx=0.0, dy=0.99)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "Figure1_pipeline.%s" % ext))
print("Figure 1 written")
