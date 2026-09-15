#!/usr/bin/env python3
"""Whether a phospholipid changed leaflet during a run, molecule by molecule.

The rule is the rule of Section 2.6: a molecule belongs to the leaflet its own
head bead sits toward, taken against the tail beads of that same molecule. No
plane is drawn through the bilayer.

The home leaflet of a molecule is the leaflet it occupies in the first frame,
taken molecule by molecule. An earlier version took the home leaflet of a
species from the majority of its molecules, which is wrong for a species that
was built into one leaflet only: sphingomyelin sits in the upper leaflet and
DOPS in the lower, so the majority test flagged whole species as odd.

For each species the script reports how many molecules end the run on the side
the first frame did not put them on, how many change side at least once, and
for each such molecule the share of frames spent away from home, the number of
changes, and the lateral distance from the nearest protein bead.

Usage:
    python3 -u ident_mdtraj2.py traj.xtc system.gro [stride]
"""
import sys

import numpy as np
import mdtraj as md

PL = ("DLPC", "DLiPC", "DIPC", "PSM", "DPSM", "DOPS", "POP2", "PIP2", "POPI")
HEADS = ("PO4", "P1", "P2", "P3", "CNO")
TAILS = set()
for c in "CD":
    for i in range(1, 7):
        for ab in "AB":
            TAILS.add("%s%d%s" % (c, i, ab))


def main():
    traj_file = sys.argv[1]
    top_file = sys.argv[2]
    stride = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    top = md.load(top_file).topology
    print("ident: %s  %d atoms" % (top_file, top.n_atoms))

    species = {}
    for res in top.residues:
        if res.name not in PL:
            continue
        h = [a.index for a in res.atoms if a.name in HEADS]
        t = [a.index for a in res.atoms if a.name in TAILS]
        if not h or not t:
            continue
        species.setdefault(res.name, []).append((h[0], t))
    if not species:
        raise SystemExit("ident: no phospholipid matched")
    for k in sorted(species):
        print("  %-6s %d molecules" % (k, len(species[k])))

    prot = [a.index for a in top.atoms if a.name == "BB"]
    prot = np.array(prot, dtype=int)
    print("  protein backbone beads %d" % len(prot))
    print("  stride %d" % stride)

    up = {k: [] for k in species}
    dis = {k: [] for k in species}
    nframe = 0
    for chunk in md.iterload(traj_file, top=top_file, chunk=20, stride=stride):
        xyz = chunk.xyz
        nframe += len(xyz)
        pxy = xyz[:, prot, :2]
        for k, mols in species.items():
            hi = np.array([m[0] for m in mols], dtype=int)
            hz = xyz[:, hi, 2]
            tz = np.array([xyz[:, m[1], 2].mean(axis=1) for m in mols]).T
            up[k].append(hz > tz)
            hxy = xyz[:, hi, :2]
            d = np.sqrt(((hxy[:, :, None, :] - pxy[:, None, :, :]) ** 2).sum(-1))
            dis[k].append(d.min(axis=2))
    print("  frames read %d" % nframe)
    print("")

    for k in sorted(species):
        U = np.concatenate(up[k], axis=0)          # frames x molecules, bool
        D = np.concatenate(dis[k], axis=0).mean(axis=0)
        n = U.shape[1]
        home = U[0]                                # each molecule's own start
        away = U != home[None, :]
        share = away.mean(axis=0)
        rev = (U[1:] != U[:-1]).sum(axis=0)
        ended_away = away[-1]
        print("%-6s %d molecules | %d end on the other side | %d change side at "
              "least once | mean distance from protein %.2f nm"
              % (k, n, int(ended_away.sum()), int((rev > 0).sum()), D.mean()))
        watch = np.where((rev > 0) | ended_away)[0]
        if not len(watch):
            print("        every molecule stays in the leaflet the first frame "
                  "put it in")
            print("")
            continue
        order = watch[np.argsort(-share[watch])]
        print("        %-6s %-14s %-9s %-9s %s"
              % ("mol", "share away", "changes", "ends away", "distance nm"))
        for i in order[:20]:
            print("        %-6d %-14.3f %-9d %-9s %.2f"
                  % (i, share[i], int(rev[i]), "yes" if ended_away[i] else "no",
                     D[i]))
        if len(order) > 20:
            print("        ... %d more" % (len(order) - 20))
        real = np.where((share > 0.2) & ended_away)[0]
        print("        of these %d, %d spend more than a fifth of the run away "
              "from home and end there" % (len(watch), len(real)))
        print("")


if __name__ == "__main__":
    main()
