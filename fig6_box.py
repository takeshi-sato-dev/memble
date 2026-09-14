#!/usr/bin/env python3
"""Figure 6. The curve depends on the size of the box.

Five systems of the base composition in a 32 by 32 by 16 nm box holding seven
copies of the construct, 250 ns each, against the curve measured at 12 by 12 nm.
The point at f = -11.93% is the first of three runs that differ only in their
starting velocities (74.04, 74.13, 73.36; mean 73.84 +- 0.42). Every one of the
five was built with the overlap check in force.

The quadratic was fitted over f = -12.09 to +11.76%, and all five points of the
larger box lie inside that range.

Five points carry no fit, so the response of the larger box is drawn as the
chord between adjacent measured points, and the chord of the smaller box is
taken over the same four intervals. box_slopes.py returns the numbers.
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
Q = [-0.02829, -0.8182, 59.202]        # curve_stats.py, the six 12 x 12 nm points
BUILT = 100.0 * 71.0 / 134.0
F = np.array([11.62, 2.74, -1.60, -6.04, -11.93])   # from the topology of each 32 x 32 system
BIG = np.array([46.66, 56.42, 61.62, 67.84, 74.04])
BIG_SD = np.array([np.nan, np.nan, np.nan, np.nan, 0.42])
SMALL = np.polyval(Q, F)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.0),
                               gridspec_kw={"width_ratios": [1.35, 1.0]})
FIT_LO, FIT_HI = -12.09, 11.76      # the range the quadratic was fitted over
xin = np.linspace(FIT_LO, FIT_HI, 400)
xlo = np.linspace(-13.5, FIT_LO, 60)
xhi = np.linspace(FIT_HI, 13.5, 60)

axA.plot(xin, np.polyval(Q, xin), color=style.BLUE, lw=1.6, zorder=3,
         label="12 × 12 nm, six systems and the fit")
for xo in (xlo, xhi):
    axA.plot(xo, np.polyval(Q, xo), color=style.BLUE, lw=1.2, zorder=3,
             ls=(0, (2, 2)), alpha=0.55)
axA.axvline(FIT_LO, color=style.MUTED, lw=0.7, ls=(0, (1, 3)), zorder=1)
axA.axvline(FIT_HI, color=style.MUTED, lw=0.7, ls=(0, (1, 3)), zorder=1)
axA.annotate("", xy=(FIT_LO, 77.2), xytext=(FIT_HI, 77.2),
             arrowprops=dict(arrowstyle="<->", color=style.MUTED, lw=0.8))
axA.text((FIT_LO + FIT_HI) / 2, 77.7, "the range the fit covers", fontsize=7, color=style.MUTED,
         ha="center", va="bottom")
axA.errorbar(F, BIG, yerr=BIG_SD, fmt="s", ms=6, color=style.GREEN, mfc="white",
             mew=1.4, ecolor=style.GREEN, elinewidth=1.0, capsize=2.5, zorder=5,
             label="32 × 32 nm, five systems")
for f, b, s_ in zip(F, BIG, SMALL):
    axA.plot([f, f], [s_, b], color=style.GREEN, lw=0.8, ls=(0, (1.5, 1.5)), zorder=4)
    if b - s_ > 1:
        axA.text(f + 0.55, (b + s_) / 2, "%+.2f" % (b - s_), fontsize=7,
                 color=style.GREEN, ha="left", va="center")
axA.axhline(BUILT, color=style.MUTED, lw=0.9, ls=(0, (5, 2)), zorder=2)
axA.text(13.2, BUILT + 0.5, "the build gave 52.99%", fontsize=7.5,
         color=style.MUTED, va="bottom", ha="right")
axA.set_xlabel("phospholipid imbalance  f  (%)")
axA.set_ylabel("cholesterol of the upper leaflet\nafter the run (%)")
axA.set_xlim(-13.5, 13.5); axA.set_ylim(43, 80.5)
axA.legend(loc="lower left", handletextpad=0.5, labelspacing=0.4)
style.panel_label(axA, "A", dx=-0.13)

# ---- B: the chords, taken the same way in both boxes
lab, sb, ss = [], [], []
for i in range(len(F) - 1):
    df = F[i + 1] - F[i]
    sb.append((BIG[i + 1] - BIG[i]) / df)
    ss.append((SMALL[i + 1] - SMALL[i]) / df)
    lab.append("%+.2f%% to\n%+.2f%%" % (F[i], F[i + 1]))
idx = np.arange(len(lab)); w = 0.34
axB.bar(idx - w / 2 - 0.01, np.abs(sb), w, color=style.GREEN, label="32 × 32 nm")
axB.bar(idx + w / 2 + 0.01, np.abs(ss), w, color=style.BLUE, label="12 × 12 nm")
for k in idx:
    axB.text(k - w / 2 - 0.01, abs(sb[k]) + 0.03, "%.2f" % abs(sb[k]),
             ha="center", va="bottom", fontsize=7, color=style.INK)
    axB.text(k + w / 2 + 0.01, abs(ss[k]) + 0.03, "%.2f" % abs(ss[k]),
             ha="center", va="bottom", fontsize=7, color=style.INK)
    axB.text(k, -0.115, "ratio %.1f" % (abs(sb[k]) / abs(ss[k])), ha="center",
             va="top", fontsize=7.5, color=style.MUTED, transform=axB.get_xaxis_transform())
axB.set_xticks(idx); axB.set_xticklabels(lab)
axB.tick_params(axis="x", length=0, pad=4)
axB.set_ylabel("points of cholesterol lost per\npoint of f, over that interval")
axB.set_ylim(0, 1.75)
axB.legend(loc="upper right", handletextpad=0.5, labelspacing=0.4)
style.panel_label(axB, "B", dx=-0.2)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "Figure6_box.%s" % ext))
for k in idx:
    print("%s  32x32 %.2f  12x12 %.2f  ratio %.1f"
          % (lab[k].replace("\n", " "), sb[k], ss[k], abs(sb[k] / ss[k])))
print("differences:", [round(float(b - s_), 2) for b, s_ in zip(BIG, SMALL)])
