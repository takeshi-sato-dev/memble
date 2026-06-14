#!/usr/bin/env python3
"""
declash_gro.py

Remove inter-molecular atom overlaps from an assembled coarse-grained system
(GROMACS .gro) by gently pushing apart only pairs of beads that belong to
DIFFERENT molecules and are closer than a target distance. Intramolecular
geometry is never touched, so the carefully built single-molecule structures
(correct bond/constraint lengths) are preserved exactly; the only changes are
small rigid displacements of whole-bead positions that energy minimisation then
relaxes.

Dense membrane packing (e.g. tail ends from neighbouring lipids meeting near the
bilayer midplane) can leave bead pairs <0.05 nm apart, which produces effectively
infinite Lennard-Jones forces and aborts minimisation. Declashing to ~0.3 nm
gives finite forces so steepest-descent can proceed.

Molecule identity is taken from the gro residue number column (each lipid/ion/
water is its own residue; a multi-residue protein is treated as one molecule via
--protein-resids so its internal contacts are not disturbed).

Usage:
  declash_gro.py --gro system.gro [--target 0.30] [--iters 100]
                 [--protein-natoms N]
"""

import argparse
import math
from collections import defaultdict

import numpy as np


def read_gro(path):
    lines = open(path).read().splitlines()
    title = lines[0]
    n = int(lines[1])
    body = lines[2:2 + n]
    box = lines[2 + n]
    resid = np.empty(n, dtype=np.int64)
    coords = np.empty((n, 3), dtype=np.float64)
    meta = []
    for i, ln in enumerate(body):
        resid[i] = int(ln[0:5])
        meta.append((ln[0:5], ln[5:10], ln[10:15], ln[15:20]))
        coords[i] = (float(ln[20:28]), float(ln[28:36]), float(ln[36:44]))
    return title, n, meta, coords, box


def write_gro(path, title, meta, coords, box):
    n = len(coords)
    with open(path, "w") as fh:
        fh.write(title + "\n")
        fh.write("%5d\n" % n)
        for (r, rn, a, idx), xyz in zip(meta, coords):
            fh.write("%5s%5s%5s%5s%8.3f%8.3f%8.3f\n"
                     % (r, rn, a, idx, xyz[0], xyz[1], xyz[2]))
        fh.write(box + "\n")


def molecule_id(meta, resid, lipids, solvions):
    """Molecule id per atom. Lipids/solvent/ions are each their own molecule
    (id = residue number). Everything else (the protein, possibly multi-residue)
    shares one id so it is declashed as a single rigid body and never distorted
    internally."""
    n = len(meta)
    mol = np.empty(n, dtype=np.int64)
    known = set(lipids) | set(solvions)
    for i in range(n):
        rn = meta[i][1].strip()
        mol[i] = resid[i] if rn in known else -1
    return mol


def declash(coords, mol, target, iters, active=None, box=None, frozen=None):
    """Rigid-molecule declash: when two beads from different molecules are closer
    than `target`, translate BOTH whole molecules apart along the contact vector.
    `active` is an optional boolean mask; beads set False (e.g. reconstructed
    virtual sites) are ignored when detecting contacts but still move with their
    molecule. `frozen` is an optional set of molecule ids that must not move
    (e.g. the protein, mol id -1): clashes against them are resolved by moving
    only the other molecule, so the protein is never distorted or re-split.
    If `box` (bx,by,bz) is given, distances use the minimum image and the grid is
    toroidal, so contacts that cross the periodic boundary (a bead pushed just
    outside the box clashing with the opposite face) are caught too.
    Vectorised with a numpy cell list so it is fast on ~10^5 beads."""
    n = len(coords)
    if active is None:
        active = np.ones(n, dtype=bool)
    if frozen is None:
        frozen = set()
    mol_atoms = defaultdict(list)
    for i in range(n):
        mol_atoms[mol[i]].append(i)
    mol_atoms = {m: np.array(idx) for m, idx in mol_atoms.items()}

    use_pbc = box is not None
    if use_pbc:
        boxv = np.array(box[:3], dtype=float)
        # wrap each molecule rigidly so its centroid sits inside [0,box)
        for m, idx in mol_atoms.items():
            c = coords[idx].mean(axis=0)
            shift = np.floor(c / boxv) * boxv
            if np.any(shift):
                coords[idx] -= shift[None, :]
        ncells = np.maximum(np.floor(boxv / target).astype(int), 1)

    def mindiff(a, b):
        d = a - b
        if use_pbc:
            d -= np.round(d / boxv) * boxv
        return d

    cell = target
    t2 = target * target
    max_step = 0.5 * target
    remaining = 0
    import sys as _sys
    rng = np.random.default_rng(0)
    for it in range(iters):
        keys = np.floor(coords / cell).astype(np.int64)
        if use_pbc:
            keys[:, 0] %= ncells[0]; keys[:, 1] %= ncells[1]; keys[:, 2] %= ncells[2]
        grid = defaultdict(list)
        for i in range(n):
            if active[i]:
                grid[(keys[i, 0], keys[i, 1], keys[i, 2])].append(i)
        trans = defaultdict(lambda: np.zeros(3))
        cnt = defaultdict(int)
        nclash = 0
        offsets = [(dx, dy, dz) for dx in (-1, 0, 1)
                   for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
        for (cx, cy, cz), ai in grid.items():
            cand = []
            for dx, dy, dz in offsets:
                if use_pbc:
                    key = ((cx + dx) % ncells[0], (cy + dy) % ncells[1],
                           (cz + dz) % ncells[2])
                else:
                    key = (cx + dx, cy + dy, cz + dz)
                cand += grid.get(key, ())
            if not cand:
                continue
            cand = np.array(cand)
            ai = np.array(ai)
            Pi = coords[ai]
            Pc = coords[cand]
            d = Pi[:, None, :] - Pc[None, :, :]
            if use_pbc:
                d -= np.round(d / boxv) * boxv
            r2 = np.einsum("abc,abc->ab", d, d)
            mi = mol[ai][:, None]
            mc = mol[cand][None, :]
            mask = (r2 < t2) & (r2 > 0) & (mi != mc)
            if not mask.any():
                continue
            aidx, bidx = np.where(mask)
            for k in range(len(aidx)):
                i = ai[aidx[k]]
                j = cand[bidx[k]]
                if j <= i:
                    continue
                vec = mindiff(coords[j], coords[i])
                r = math.sqrt(float(vec @ vec))
                if r < 1e-9:
                    vec = rng.normal(size=3); vec /= np.linalg.norm(vec); r = 1e-3
                push = 0.5 * (target - r) / r
                trans[mol[i]] -= push * vec
                trans[mol[j]] += push * vec
                cnt[mol[i]] += 1
                cnt[mol[j]] += 1
                nclash += 1
        if nclash == 0:
            remaining = 0
            print("declash: resolved after %d iterations" % (it + 1))
            break
        remaining = nclash
        for m, t in trans.items():
            if cnt[m] == 0 or m in frozen:
                continue
            step = t
            nrm = float(np.linalg.norm(step))
            if nrm > max_step:
                step = step * (max_step / nrm)
            elif nrm < 0.05 * target:
                jit = rng.normal(size=3)
                step = step + jit / np.linalg.norm(jit) * (0.3 * target)
            coords[mol_atoms[m]] += step[None, :]
        if it % 10 == 0:
            print("declash: iter %d, %d contacts remaining" % (it, nclash))
            _sys.stdout.flush()
    # final wrap into the box
    if use_pbc:
        for m, idx in mol_atoms.items():
            if m in frozen:
                continue
            c = coords[idx].mean(axis=0)
            shift = np.floor(c / boxv) * boxv
            if np.any(shift):
                coords[idx] -= shift[None, :]
    return coords, remaining


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--target", type=float, default=0.30,
                    help="minimum inter-molecular distance to enforce (nm)")
    ap.add_argument("--iters", type=int, default=100)
    ap.add_argument("--lipids", default="",
                    help="space-separated lipid resnames present (each treated "
                         "as its own rigid molecule)")
    ap.add_argument("--exclude-beads", default="",
                    help="space-separated bead names to ignore in clash detection "
                         "(reconstructed virtual sites such as CHOL ROH/R3 are "
                         "rebuilt by mdrun each step and only pollute contacts)")
    ap.add_argument("--freeze-protein", action="store_true",
                    help="never move the protein (mol id -1); resolve clashes "
                         "against it by moving only the other molecule, so the "
                         "protein is not distorted or re-split across PBC")
    args = ap.parse_args()

    title, n, meta, coords, box = read_gro(args.gro)
    resid = np.array([int(m[0]) for m in meta], dtype=np.int64)
    solvions = ["W", "WF", "NA", "CL", "ION", "NA+", "CL-"]
    mol = molecule_id(meta, resid, args.lipids.split(), solvions)

    skip = set(args.exclude_beads.split())
    active = np.array([m[1].strip() not in skip for m in meta]) \
        if skip else np.ones(n, dtype=bool)

    frozen = {-1} if args.freeze_protein else set()
    boxvals = [float(v) for v in box.split()[:3]]
    coords, remaining = declash(coords, mol, args.target, args.iters, active,
                                boxvals, frozen)
    write_gro(args.gro, title, meta, coords, box)
    if remaining:
        print("declash: %d inter-molecular contacts still < %.2f nm after %d "
              "iters (minimisation will finish the job)"
              % (remaining, args.target, args.iters))
    else:
        print("declash: no inter-molecular contacts < %.2f nm remain"
              % args.target)


if __name__ == "__main__":
    main()
