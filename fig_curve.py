"""Figure 4. Where cholesterol settles against the phospholipid imbalance.

A: the six 12 by 12 nm systems, the curve through them, the composition the
   build assigned, and the one 32 by 32 nm system.
B: the three routes of Table 3 placed on the same curve.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Written beside this script unless a directory is given on the command line.
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

SURF = "#fcfcfb"; INK = "#0b0b0b"; INK2 = "#52514e"
BLUE = "#2a78d6"; GREY = "#9a9892"; ORANGE = "#d1701f"; RED = "#b5322e"

# ---- the six points of the curve, 12 by 12 by 16 nm, 250 ns each ----------
X = np.array([-11.76, -6.23, -1.83, 2.94, 7.35, 11.76])
Y = np.array([64.78, 63.20, 61.17, 56.86, 50.26, 46.36])
SD = np.array([2.69, np.nan, np.nan, 1.71, np.nan, np.nan])
N = np.array([3, 1, 1, 3, 1, 1])
SE = SD / np.sqrt(N)

BUILT = 100 * 71 / 134.0            # the share the build gave every system
NCHOL = 134.0

C = np.polyfit(X, Y, 2)             # 59.193 - 0.8201 f - 0.02805 f^2
FIX = 6.239                         # where the curve crosses BUILT

# ---- the one 32 by 32 by 16 nm system, 300 ns -----------------------------
XB, YB, EB = -11.93, 74.47, 0.53

# ---- the three routes of Table 3, 500 ns, three seeds each ---------------
RX = np.array([-0.37, 1.10, 2.94])
RMOV = np.array([8.8, 7.1, 2.1])
RERR = np.array([2.8, 0.3, 2.4])
RNAME = ["equal\nnumbers", "areas from the\nmeasurement", "areas from\na table"]

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2,
                     "axes.linewidth": 0.8})

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.1), facecolor=SURF,
                             gridspec_kw={"width_ratios": [1.5, 1]})

# =========================== A ============================================
ax.set_facecolor(SURF)
g = np.linspace(-13.5, 13.0, 400)
ax.plot(g, np.polyval(C, g), color=BLUE, lw=1.8, zorder=3)
ax.axhline(BUILT, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2)
ax.plot([FIX, FIX], [40, BUILT], color=INK2, lw=0.8, ls=(0, (1, 2)), zorder=2)
ax.plot([FIX], [BUILT], marker="o", ms=5, mfc=SURF, mec=INK2, mew=1.2, zorder=6)

ax.errorbar(X, Y, yerr=SE, fmt="o", ms=5, color=BLUE, mfc=BLUE, mec=BLUE,
            ecolor=BLUE, elinewidth=1.2, capsize=2.5, zorder=5)
ax.errorbar([XB], [YB], yerr=[EB], fmt="s", ms=6.5, color=RED, mfc=SURF,
            mec=RED, mew=1.6, ecolor=RED, elinewidth=1.2, capsize=2.5, zorder=7)

ax.annotate("", xy=(XB, YB - 0.7), xytext=(XB, np.polyval(C, XB) + 0.7),
            arrowprops=dict(arrowstyle="<->", color=RED, lw=1.0,
                            shrinkA=0, shrinkB=0), zorder=6)
ax.text(XB + 0.7, 0.5 * (YB + np.polyval(C, XB)), "9.5 pt", color=RED,
        fontsize=8, fontweight="bold", va="center")
ax.text(-13.2, 75.6, "32 × 32 nm, one system", color=RED, fontsize=8,
        fontweight="bold")
ax.text(-13.2, 43.4, "12 × 12 nm, six systems", color=BLUE, fontsize=8,
        fontweight="bold")
ax.text(12.6, BUILT + 0.9, "as built, 52.99%", color=INK2, fontsize=7,
        ha="right")
ax.text(FIX + 0.5, 41.6, "+6.2%\nnothing moves", color=INK2, fontsize=7)

ax.set_xlim(-13.5, 13.0); ax.set_ylim(41, 78)
ax.set_xlabel("phospholipid imbalance, upper minus lower (% of the total)",
              color=INK)
ax.set_ylabel("cholesterol in the upper leaflet (%)", color=INK)
ax.spines[["top", "right"]].set_visible(False)
ax.text(-0.16, 1.04, "A", transform=ax.transAxes, fontsize=11,
        fontweight="bold", color=INK)

# =========================== B ============================================
bx.set_facecolor(SURF)
g2 = np.linspace(-4.0, 6.0, 200)
bx.plot(g2, (np.polyval(C, g2) - BUILT) * NCHOL / 100.0, color=BLUE, lw=1.8,
        zorder=3)
bx.axhline(0, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2)
bx.errorbar(RX, RMOV, yerr=RERR, fmt="D", ms=5.5, color=ORANGE, mfc=SURF,
            mec=ORANGE, mew=1.6, ecolor=ORANGE, elinewidth=1.2, capsize=2.5,
            zorder=5)
for x, y, nm in zip(RX, RMOV, RNAME):
    bx.annotate(nm, (x, y), textcoords="offset points", xytext=(7, 4),
                fontsize=6.5, color=INK2)
bx.plot([FIX], [0], marker="o", ms=5, mfc=SURF, mec=INK2, mew=1.2, zorder=6)
bx.set_xlim(-4.0, 6.6); bx.set_ylim(-1.5, 13.5)
bx.set_xlabel("phospholipid imbalance (%)", color=INK)
bx.set_ylabel("cholesterol molecules that changed leaflet", color=INK)
bx.spines[["top", "right"]].set_visible(False)
bx.text(-0.22, 1.04, "B", transform=bx.transAxes, fontsize=11,
        fontweight="bold", color=INK)
bx.text(-3.7, 2.4, "the curve of panel A", color=BLUE, fontsize=7.5,
        fontweight="bold")
bx.text(-3.7, 1.3, "the routes of Table 3", color=ORANGE, fontsize=7.5,
        fontweight="bold")

fig.tight_layout(pad=0.8)
fig.savefig(os.path.join(OUT, "Figure4_curve.svg"), facecolor=SURF)
fig.savefig(os.path.join(OUT, "Figure4_curve.png"), dpi=300, facecolor=SURF)
print("curve  C = %.3f %+.4f f %+.5f f^2" % (C[2], C[1], C[0]))
print("predicted at the three routes:",
      np.round((np.polyval(C, RX) - BUILT) * NCHOL / 100.0, 1))
print("written to %s" % OUT)
