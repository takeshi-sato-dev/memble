#!/usr/bin/env python3
"""Figure 5. Where cholesterol settles, against the phospholipid imbalance.

Every number here is the one curve_stats.py returns from the six points. The
points are the mean over the last 40% of a 250 ns run of the number of
cholesterol molecules in the upper leaflet, over the 134 cholesterol molecules
of the system. The imbalance f of each point is computed from the phospholipid
numbers in that system's system.top, never from delta.

Panel A  the six points, the quadratic, and the imbalance that holds the build
Panel B  the three builds of Table 4, which were not used in the fit
"""
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

PL_UPPER = np.array([120, 128, 134, 140, 146, 152])
PL_LOWER = np.array([153, 145, 139, 132, 126, 120])
F = 100.0 * (PL_UPPER - PL_LOWER) / (PL_UPPER + PL_LOWER)
C = np.array([64.78, 63.20, 61.17, 56.86, 50.26, 46.36])
SD = np.array([2.69, np.nan, np.nan, 1.71, np.nan, np.nan])
N = np.array([3, 1, 1, 3, 1, 1])
BUILT = 100.0 * 71.0 / 134.0
Q = np.polyfit(F, C, 2)
FIX, FIXSD = 6.25, 0.94

ROUTE = [("EqN", -0.37, 57.59, 2.26), ("areas measured on the packed system", 1.10, 57.43, 0.55),
         ("SA, from a table of areas", 2.94, 54.70, 1.84)]
SHORT = ["EqN", "measured areas", "SA (table)"]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.0))
x = np.linspace(-13.5, 13.5, 400)

# ---------------- A
axA.plot(x, np.polyval(Q, x), color=style.BLUE, lw=1.6, zorder=3)
err = np.where(N == 3, SD / np.sqrt(3.0), np.nan)
axA.errorbar(F, C, yerr=err, fmt="o", ms=5, color=style.BLUE, mfc="white",
             mew=1.4, ecolor=style.BLUE, elinewidth=1.0, capsize=2.5, zorder=4)
axA.axhline(BUILT, color=style.MUTED, lw=0.9, ls=(0, (5, 2)), zorder=2)
axA.plot([FIX], [BUILT], "o", ms=7, mfc="white", mec=style.INK, mew=1.5, zorder=5)
axA.errorbar([FIX], [BUILT], xerr=[FIXSD], fmt="none", ecolor=style.INK,
             elinewidth=1.0, capsize=2.5, zorder=5)
axA.annotate("f = +6.25 ± 0.94%\nholds the build" , xy=(FIX, BUILT),
             xytext=(8.2, 60.5), fontsize=7.5, color=style.INK, ha="left",
             arrowprops=dict(arrowstyle="-", lw=0.6, color=style.INK,
                             shrinkA=0, shrinkB=4))
axA.text(-13.2, BUILT + 0.5, "the build gave 52.99%", fontsize=7.5,
         color=style.MUTED, va="bottom")
axA.set_xlabel("phospholipid imbalance  f  (%)")
axA.set_ylabel("cholesterol of the upper leaflet\nafter the run (%)")
axA.set_xlim(-13.5, 13.5); axA.set_ylim(43, 69)
style.panel_label(axA, "A")

# ---------------- B, in the units the curve carries
axB.plot(x, np.polyval(Q, x), color=style.BLUE, lw=1.6, zorder=3)
axB.axhline(BUILT, color=style.MUTED, lw=0.9, ls=(0, (5, 2)), zorder=2)
for i, (lab, f, mv, se) in enumerate(ROUTE):
    axB.errorbar([f], [mv], yerr=[se], fmt=style.MARKERS[i], ms=6,
                 color=style.ORANGE, mfc="white", mew=1.4, ecolor=style.ORANGE,
                 elinewidth=1.0, capsize=2.5, zorder=5, label=SHORT[i])
axB.plot([FIX], [BUILT], "o", ms=7, mfc="white", mec=style.INK, mew=1.5, zorder=6)
axB.text(FIX + 0.4, BUILT + 0.8, "f = +6.25%", fontsize=7.5, color=style.INK)
axB.set_xlabel("phospholipid imbalance  f  (%)")
axB.set_ylabel("cholesterol of the upper leaflet\nafter the run (%)")
axB.set_xlim(-4.5, 9); axB.set_ylim(51.5, 62.0)
axB.legend(loc="lower left", handletextpad=0.4, labelspacing=0.4,
           borderaxespad=0.2)
style.panel_label(axB, "B")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "Figure5_curve.%s" % ext))

print("C = %.3f %+.4f f %+.5f f^2" % (Q[2], Q[1], Q[0]))
print("residual rms %.2f pt, max %.2f pt"
      % (np.sqrt(np.mean((C - np.polyval(Q, F)) ** 2)),
         np.max(np.abs(C - np.polyval(Q, F)))))
print("value at f = 0: %.2f%%   build: %.2f%%" % (np.polyval(Q, 0.0), BUILT))
for lab, f, mv, se in ROUTE:
    print("  %-34s curve %5.2f%%  measured %5.2f%% +- %.2f   diff %+.2f"
          % (lab.replace("\n", " "), np.polyval(Q, f), mv, se, mv - np.polyval(Q, f)))
