#!/usr/bin/env python3
"""The phospholipid imbalance a production run carries, taken over the run.

The leaflet rule of Section 2.6 reads a molecule from one frame, and in a
strained bilayer some molecules tilt far enough that one frame does not decide
which leaflet they belong to. This script decides each molecule by the leaflet
it occupies in most frames of the settled part of the production run, which is
the same window Section 2.6 averages the cholesterol over, and then counts the
two leaflets.

For every molecule it also reports how firmly it is assigned: a molecule on one
side in every frame is certain, and a molecule near half is not. The imbalance
is given three ways, over all molecules, over the molecules that are on one side
in at least 90% of frames, and over those above 70%, so that a reader can see
whether the answer depends on the molecules that are not decided.

Usage:
    python3 f_from_run.py prod.gro prod.xtc
    python3 f_from_run.py prod.gro prod.xtc --last-fraction 0.6 --stride 5
    python3 f_from_run.py --work ~/counts_run/curve/d-20/d-20_work
"""
import argparse
import glob
import os
import re

import numpy as np

TAIL = re.compile(r"^[CD]\d[AB]$")
HEADS = "PO4 P1 P2 P3 CNO"
PL = ("DLPC", "DLiPC", "DIPC", "PSM", "DPSM", "DOPS", "POP2", "POP2_45",
      "PIP2", "POPI")


def pick(work):
    for g, x in (("step7_production.gro", "step7_production.xtc"),
                 ("step7.gro", "step7.xtc"),
                 ("prod.gro", "prod.xtc"),
                 ("md.gro", "md.xtc")):
        a = os.path.join(work, g)
        b = os.path.join(work, x)
        if os.path.exists(a) and os.path.exists(b):
            return a, b
    gs = sorted(glob.glob(os.path.join(work, "*.gro")))
    xs = sorted(glob.glob(os.path.join(work, "*.xtc")))
    if gs and xs:
        return gs[0], xs[0]
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure", nargs="?")
    ap.add_argument("trajectory", nargs="?")
    ap.add_argument("--work", default=None)
    ap.add_argument("--last-fraction", type=float, default=0.6)
    ap.add_argument("--stride", type=int, default=5)
    a = ap.parse_args()

    if a.work:
        s, t = pick(os.path.expanduser(a.work))
        if not s:
            raise SystemExit("f_from_run: no structure and trajectory under %s"
                             % a.work)
    else:
        s, t = a.structure, a.trajectory
    print("f_from_run: %s" % os.path.basename(s))
    print("f_from_run: %s" % os.path.basename(t))

    import MDAnalysis as mda
    u = mda.Universe(s, t)
    n = len(u.trajectory)
    start = int((1.0 - a.last_fraction) * n)
    print("  %d frames, the last %d averaged, stride %d"
          % (n, n - start, a.stride))

    mols = []
    for name in PL:
        sel = u.select_atoms("resname %s" % name)
        if not len(sel):
            continue
        for r in sel.residues:
            h = r.atoms.select_atoms("name %s" % HEADS)
            tl = [at.index for at in r.atoms if TAIL.match(at.name)]
            if len(h) and tl:
                mols.append((name, h.indices[0], tl))
    if not mols:
        raise SystemExit("f_from_run: no phospholipid matched")
    hi = np.array([m[1] for m in mols])
    print("  %d phospholipid molecules" % len(mols))

    up = []
    for ts in u.trajectory[start::a.stride]:
        p = u.atoms.positions
        hz = p[hi, 2]
        tz = np.array([p[m[2], 2].mean() for m in mols])
        up.append(hz > tz)
    U = np.array(up)
    share = U.mean(axis=0)
    print("  %d frames read" % U.shape[0])

    names = np.array([m[0] for m in mols])
    print("")
    print("%-8s %-14s %-14s %s" % ("species", "upper/lower", "firm 0.9",
                                   "undecided 0.3-0.7"))
    for name in sorted(set(names)):
        k = names == name
        sh = share[k]
        u1 = int((sh > 0.5).sum())
        firm = int(((sh > 0.9) | (sh < 0.1)).sum())
        mid = int(((sh > 0.3) & (sh < 0.7)).sum())
        print("%-8s %6d /%6d  %6d of %-5d %d"
              % (name, u1, len(sh) - u1, firm, len(sh), mid))

    def imbalance(mask):
        sh = share[mask]
        if not len(sh):
            return float("nan"), 0, 0
        u1 = int((sh > 0.5).sum())
        l1 = len(sh) - u1
        return 100.0 * (u1 - l1) / (u1 + l1), u1, l1

    print("")
    for label, mask in (("all molecules", np.ones(len(share), bool)),
                        ("firm, 0.9 or above", (share > 0.9) | (share < 0.1)),
                        ("firm, 0.7 or above", (share > 0.7) | (share < 0.3))):
        f, u1, l1 = imbalance(mask)
        print("  %-22s upper %3d  lower %3d   f = %+.2f%%" % (label, u1, l1, f))
    print("")
    print("  molecules whose leaflet is not decided by the run: %d of %d"
          % (int(((share > 0.3) & (share < 0.7)).sum()), len(share)))


if __name__ == "__main__":
    main()
