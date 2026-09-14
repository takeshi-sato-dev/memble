#!/usr/bin/env python3
"""Figure 2. Where the packing step leaves the membrane-spanning helix.

Two builds of one system, a bilayer of cholesterol, DLPC and palmitoyl
sphingomyelin at 1:1:1 in a 10 by 10 by 16 nm box, holding residues 1 to 117 of
the receptor. The membrane-spanning range is residues 66 to 88.

A  the packing step was not told which residues cross the membrane, so it
   centred the protein on the centroid of the whole molecule
B  MEMBLE named the range 65 to 88 to the packing step as the residues to
   centre on

Coordinates are the packed systems as the build report saw them, memble_bench/figA and figB. The
midplane of each is taken as the mean z of the two phosphate planes, and every z
is drawn against it.
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
TM_FIRST, TM_LAST = 65, 88          # the range MEMBLE was given, as in the report


def read_gro(path):
    lines = open(path).read().splitlines()
    n = int(lines[1])
    resid, resn, name, z = [], [], [], []
    for L in lines[2:2 + n]:
        resid.append(int(L[0:5])); resn.append(L[5:10].strip())
        name.append(L[10:15].strip()); z.append(float(L[36:44]))
    box = [float(v) for v in lines[2 + n].split()[:3]]
    return (np.array(resid), np.array(resn), np.array(name), np.array(z), box)


def panel(ax, path, title, tag):
    resid, resn, name, z, box = read_gro(path)
    lipid = np.isin(resn, ["CHOL", "DLPC", "PSM"])
    prot = ~lipid & ~np.isin(resn, ["W", "NA", "CL", "ION"])
    po4 = lipid & (name == "PO4")
    up = po4 & (z > np.median(z[po4]))
    lo = po4 & ~up
    # MEMBLE takes the midplane as the mean z of every lipid bead, because a
    # sterol has no head bead and a mean over the two head planes is a height at
    # which nothing lies. verify_system.py line 326. The figure uses the same
    # definition, so the offset drawn here is the offset the report prints.
    mid = float(np.mean(z[lipid]))
    zz = z - mid
    tm = prot & (resid >= TM_FIRST) & (resid <= TM_LAST)
    bb = tm   # verify_system.py averages every bead of the range, line 420

    rng = np.random.default_rng(0)
    ax.scatter(rng.uniform(0, 1, lipid.sum()), zz[lipid], s=1.2,
               color=style.GREY, alpha=0.30, lw=0, zorder=1)
    ax.scatter(rng.uniform(0, 1, po4.sum()), zz[po4], s=5, color=style.MUTED,
               alpha=0.75, lw=0, zorder=2)
    ax.scatter(rng.uniform(0.28, 0.72, (prot & ~tm).sum()), zz[prot & ~tm], s=5,
               color=style.ORANGE, alpha=0.55, lw=0, zorder=3)
    ax.scatter(rng.uniform(0.32, 0.68, tm.sum()), zz[tm], s=9, color=style.BLUE,
               alpha=0.95, lw=0, zorder=4)

    for v in (zz[up].mean(), zz[lo].mean()):
        ax.axhline(v, color=style.MUTED, lw=0.8, ls=(0, (5, 2)), zorder=5)
    ax.axhline(0, color=style.INK, lw=0.7, zorder=5)
    off = zz[bb].mean()
    ax.plot([0.06], [off], marker="_", ms=26, mew=2.0, color=style.BLUE, zorder=6)
    ax.annotate("%.2f nm" % off, xy=(0.06, off), xytext=(0.03, off - 1.05),
                fontsize=8, color=style.BLUE, ha="left",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=style.BLUE))
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-4.6, 5.4)
    ax.set_xticks([])
    ax.set_title(title, fontsize=8, pad=6)
    ax.spines["bottom"].set_visible(False)
    style.panel_label(ax, tag, dx=-0.1)
    return (zz[up].mean(), zz[lo].mean(), off,
            int((tm & (zz > zz[up].mean())).sum()), int(tm.sum()), int(bb.sum()))


fig, (axA, axB) = plt.subplots(1, 2, figsize=(5.6, 3.4), sharey=True)
a = panel(axA, os.path.join(DATA, "figA_system.gro"),
          "centred on the whole molecule", "A")
b = panel(axB, os.path.join(DATA, "figB_system.gro"),
          "centred on residues 65–88", "B")
axA.set_ylabel("z against the bilayer midplane (nm)")
axB.text(1.02, a[0], "phosphate\nplanes", transform=axB.get_yaxis_transform(),
         fontsize=7.5, color=style.MUTED, va="center", ha="left")
h = [plt.Line2D([], [], ls="", marker="o", ms=4, color=c, label=l)
     for c, l in ((style.BLUE, "residues 65–88"), (style.ORANGE, "the rest of the protein"),
                  (style.MUTED, "PO4"), (style.GREY, "other lipid beads"))]
fig.legend(handles=h, loc="lower center", ncol=4, frameon=False,
           handletextpad=0.4, columnspacing=1.4, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "Figure2_placement.%s" % ext))
for tag, r in (("A", a), ("B", b)):
    print("%s upper PO4 %+.2f  lower PO4 %+.2f  centroid of the range %+.3f nm  "
          "beads of the range above the upper plane %d of %d  (%d beads averaged)"
          % ((tag,) + r))
