#!/usr/bin/env python3
"""
check_min_distance.py

Report the smallest inter-molecular bead-bead distance (minimum image / PBC) in
a built system, so the build can confirm there is no overlap that would give an
infinite force at minimisation. Coincident or near-coincident beads from
different molecules are what make GROMACS report Fmax = inf; this catches them
before the user ever runs gmx.

Prints the minimum distance and the worst offending pair. Exit code is non-zero
if anything is below --min (default 0.15 nm), so the build can warn loudly.

Usage:
  check_min_distance.py --gro system.gro --lipids "CHOL DLPC PSM DOPS" [--min 0.15]
"""

import argparse
import sys
import numpy as np
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--lipids", default="")
    ap.add_argument("--min", type=float, default=0.15)
    ap.add_argument("--exclude-beads", default="ROH R3")
    args = ap.parse_args()

    L = open(args.gro).read().splitlines()
    n = int(L[1])
    body = L[2:2 + n]
    box = np.array([float(v) for v in L[2 + n].split()[:3]])

    resid = np.array([int(b[0:5]) for b in body])
    resn = [b[5:10].strip() for b in body]
    aname = [b[10:15].strip() for b in body]
    X = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                  for b in body])

    known = set(args.lipids.split()) | {"W", "WF", "NA", "CL", "ION", "NA+", "CL-"}
    # molecule id: lipids/solvent/ions by residue, protein (rest) all -1
    mol = np.array([resid[i] if resn[i] in known else -1 for i in range(n)])
    skip = set(args.exclude_beads.split())
    use = np.array([aname[i] not in skip for i in range(n)])

    # also flag the worst INTRA-molecule close contact (same molecule, different
    # bead): a template that places two beads on top of each other shows up here
    # even though it is not an inter-molecular clash.
    intra_best = 1e9
    intra_pair = None
    bymol = defaultdict(list)
    for i in range(n):
        if use[i]:
            bymol[mol[i]].append(i)
    for m, ids in bymol.items():
        if m == -1 or len(ids) < 2:
            continue
        # only need one representative molecule per residue id; residue ids are
        # unique per lipid, so each group is a single molecule
        a = np.array(ids)
        P = X[a]
        d = P[:, None, :] - P[None, :, :]
        r = np.sqrt(np.einsum("abc,abc->ab", d, d))
        np.fill_diagonal(r, 1e9)
        if r.min() < intra_best:
            ii, jj = np.unravel_index(np.argmin(r), r.shape)
            intra_best = r.min()
            intra_pair = (a[ii], a[jj])
        if len(bymol) > 4000 and intra_best < 0.05:
            break

    # cell list at the query distance
    cell = max(args.min * 3, 0.3)
    ncells = np.maximum(np.floor(box / cell).astype(int), 1)
    keys = np.floor(X / cell).astype(np.int64)
    keys[:, 0] %= ncells[0]; keys[:, 1] %= ncells[1]; keys[:, 2] %= ncells[2]
    grid = defaultdict(list)
    for i in range(n):
        if use[i]:
            grid[(keys[i, 0], keys[i, 1], keys[i, 2])].append(i)

    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)]
    best = 1e9
    bestpair = None
    for (cx, cy, cz), ai in grid.items():
        cand = []
        for dx, dy, dz in offs:
            cand += grid.get(((cx + dx) % ncells[0], (cy + dy) % ncells[1],
                              (cz + dz) % ncells[2]), ())
        ai = np.array(ai); cand = np.array(cand)
        d = X[ai][:, None, :] - X[cand][None, :, :]
        d -= np.round(d / box) * box
        r = np.sqrt(np.einsum("abc,abc->ab", d, d))
        mi = mol[ai][:, None]; mc = mol[cand][None, :]
        r[(mi == mc)] = 1e9
        r[r == 0] = 1e9
        if r.size and r.min() < best:
            a, b = np.unravel_index(np.argmin(r), r.shape)
            best = r.min()
            bestpair = (ai[a], cand[b])

    if bestpair is None:
        print("check_min_distance: no inter-molecular pairs found")
        return
    i, j = bestpair
    msg = ("check_min_distance: smallest inter-molecular distance = %.3f nm "
           "between %s %d:%s and %s %d:%s"
           % (best, resn[i], resid[i], aname[i], resn[j], resid[j], aname[j]))
    print(msg)
    if intra_pair is not None:
        ii, jj = intra_pair
        print("check_min_distance: smallest intra-molecular contact  = %.3f nm "
              "between %s %d %s and %s"
              % (intra_best, resn[ii], resid[ii], aname[ii], aname[jj]))
    worst = min(best, intra_best)
    if worst < args.min:
        print("check_min_distance: WARNING below %.2f nm -- minimisation may see "
              "a very large force here" % args.min)
        sys.exit(1)
    else:
        print("check_min_distance: OK, no overlaps below %.2f nm "
              "(system should minimise without infinite forces)" % args.min)


if __name__ == "__main__":
    main()
