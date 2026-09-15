#!/usr/bin/env python3
"""Whether a phospholipid that reads as crossed has crossed, or is misread.

Section 3.3 reports that no phospholipid changed leaflet. The 12 by 12 nm
systems report a few molecules on the side the build did not put them on, and a
count alone cannot say which of three things happened:

  the rule misreads a molecule the protein has pushed out of the plane, and the
  same molecule reads wrong in every frame, from the first one

  the rule flickers on a molecule that sits near the midplane, and that molecule
  changes side and returns within one or two frames

  the molecule crossed, once, and stayed

The home leaflet of a molecule is the leaflet it occupies in the first frame,
taken molecule by molecule. The first version of this script took the home
leaflet of a species from the majority of its molecules, which is wrong for a
species built into one leaflet only: sphingomyelin sits in the upper leaflet and
DOPS in the lower, so the majority test flagged an entire species as odd.

For every phospholipid the script records the leaflet in each frame and reports

  share away    the share of frames the molecule spent away from its own first
                frame. A molecule at 1/N changed side for one frame and returned
  changes       how many times the assignment changed
  ends away     whether the last frame has it away from home
  distance      the lateral distance from the nearest protein bead

A misreading gives share 1.000 with 0 changes, because the molecule reads wrong
from the first frame on. A flicker gives a share of a few frames with an even
number of changes. A crossing gives one change, ends away, and a share set by
when it happened.

Usage:
    python3 -u leaflet_identity2.py prod.gro prod.xtc
    python3 -u leaflet_identity2.py prod.gro prod.xtc --stride 5
    python3 -u leaflet_identity2.py prod.gro prod.xtc --build memble_build.json
"""
import argparse
import json
import os
import re

import numpy as np

PL = ("DLPC", "DLiPC", "DIPC", "PSM", "DPSM", "DOPS", "POP2", "PIP2", "POPI")
HEADS = ("PO4", "P1", "P2", "P3", "CNO")
TAIL = re.compile(r"^[CD]\d[AB]$")


def as_counts(v):
    if isinstance(v, dict):
        return {str(k): int(x) for k, x in v.items()}
    out = {}
    for tok in str(v).replace(",", " ").split():
        if ":" not in tok:
            continue
        name, _, num = tok.partition(":")
        num = num.split(":")[0]
        try:
            out[name] = int(float(num))
        except ValueError:
            continue
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure")
    ap.add_argument("trajectory")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--near-nm", type=float, default=2.0)
    ap.add_argument("--build", default=None,
                    help="memble_build.json, to compare the first frame with "
                         "what the build assigned")
    a = ap.parse_args()

    import MDAnalysis as mda
    u = mda.Universe(a.structure, a.trajectory)
    prot = u.select_atoms("name BB")
    print("leaflet_identity: %d frames, %d protein backbone beads"
          % (len(u.trajectory), len(prot)))

    species = {}
    for name in PL:
        sel = u.select_atoms("resname %s" % name)
        if not len(sel):
            continue
        mols = []
        for r in sel.residues:
            h = r.atoms.select_atoms("name %s" % " ".join(HEADS))
            t = [at.index for at in r.atoms if TAIL.match(at.name)]
            if len(h) >= 1 and t:
                mols.append((h.indices[0], t))
        if mols:
            species[name] = mols
            print("  %-6s %d molecules" % (name, len(mols)))
    if not species:
        raise SystemExit("leaflet_identity: no phospholipid matched")

    built = None
    if a.build and os.path.exists(a.build):
        d = json.load(open(a.build))
        for key in ("composition", "leaflets", "lipids", None):
            cand = d.get(key) if key else d
            if isinstance(cand, dict) and "upper" in cand and "lower" in cand:
                built = {"upper": as_counts(cand["upper"]),
                         "lower": as_counts(cand["lower"])}
                break

    frames = {k: [] for k in species}
    dist = {k: [] for k in species}
    for ts in u.trajectory[::a.stride]:
        pxy = prot.positions[:, :2]
        for name, mols in species.items():
            hz = np.array([u.atoms.positions[h, 2] for h, _ in mols])
            tz = np.array([u.atoms.positions[t, 2].mean() for _, t in mols])
            frames[name].append(hz > tz)
            hxy = np.array([u.atoms.positions[h, :2] for h, _ in mols])
            d = np.sqrt(((hxy[:, None, :] - pxy[None, :, :]) ** 2).sum(-1))
            dist[name].append(d.min(axis=1) / 10.0)

    print("")
    for name, mols in species.items():
        U = np.array(frames[name])
        D = np.array(dist[name]).mean(axis=0)
        n = U.shape[1]
        home = U[0]
        away = U != home[None, :]
        share = away.mean(axis=0)
        rev = (U[1:] != U[:-1]).sum(axis=0)
        ended = away[-1]

        nu = int(home.sum())
        line = "%-6s %d molecules | first frame %d upper %d lower" \
               % (name, n, nu, n - nu)
        if built is not None:
            bu = built["upper"].get(name, 0)
            bl = built["lower"].get(name, 0)
            line += " | build %d upper %d lower" % (bu, bl)
        print(line)
        print("       %d end on the other side | %d change side at least once | "
              "mean distance from the protein %.2f nm"
              % (int(ended.sum()), int((rev > 0).sum()), D.mean()))

        watch = np.where((rev > 0) | ended)[0]
        if not len(watch):
            print("       every molecule stays in the leaflet the first frame "
                  "put it in")
            print("")
            continue
        order = watch[np.argsort(-share[watch])]
        print("       %-6s %-13s %-9s %-10s %s"
              % ("mol", "share away", "changes", "ends away", "distance nm"))
        for i in order[:20]:
            print("       %-6d %-13.3f %-9d %-10s %.2f"
                  % (i, share[i], int(rev[i]),
                     "yes" if ended[i] else "no", D[i]))
        if len(order) > 20:
            print("       ... %d more" % (len(order) - 20))

        crossed = np.where((rev == 1) & ended)[0]
        flick = np.where((rev >= 2) & ~ended)[0]
        print("       %d cross once and stay, %d change side and return, "
              "%d sit within %.1f nm of the protein"
              % (len(crossed), len(flick), int((D[watch] < a.near_nm).sum()),
                 a.near_nm))
        if len(crossed):
            print("       ** at least one molecule changed side and stayed. "
                  "Report it as a crossing.")
        else:
            print("       ** no molecule changed side and stayed: the far-side "
                  "molecules are misread from the first frame, and the rest "
                  "are the midplane flicker")
        print("")


if __name__ == "__main__":
    main()
