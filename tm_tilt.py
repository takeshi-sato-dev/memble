#!/usr/bin/env python3
"""
tm_tilt.py

Measure the transmembrane helix tilt angle of each protein chain in a Martini
system, for comparison with a tilt angle measured by polarized ATR-IR.

The tilt is the angle between the TM-core helix long axis and the membrane
normal (the box z axis). The helix long axis is the largest-variance direction
(first principal axis) of the TM-core backbone beads of one chain. With a
trajectory, the per-frame tilt is averaged so the value matches an
order-parameter measurement, which is itself a time and ensemble average.

This reports the same quantity for the built system that an experiment reports,
so a build can be checked against an independent measurement rather than only
against whether it runs.

Usage:
  tm_tilt.py --gro system.gro --core 12-35 [--traj prod.xtc --top system.top]
  tm_tilt.py --gro system.gro --core 12-35           # single frame
The --core range is the TM-core residue range in the itp-local numbering that
the build prints (original residue minus the kept-range start plus one).
"""

import argparse
import sys
import numpy as np


AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def parse_core(spec):
    lo, hi = spec.split("-")
    return int(lo), int(hi)


def read_gro_backbone(path, core_lo, core_hi):
    """Return a list of arrays, one per chain, of TM-core BB coordinates (nm).

    Chains are split on a resid reset (the per-chain numbering restarts), which
    is how martinize names assembled chains.
    """
    lines = open(path).read().splitlines()
    n = int(lines[1])
    body = lines[2:2 + n]
    chains = []
    cur = []
    prev_resid = None
    for b in body:
        resname = b[5:10].strip().upper()
        atom = b[10:15].strip()
        if resname not in AA or atom != "BB":
            continue
        resid = int(b[0:5].strip())
        if prev_resid is not None and resid < prev_resid:
            if cur:
                chains.append(cur)
            cur = []
        if core_lo <= resid <= core_hi:
            x = float(b[20:28])
            y = float(b[28:36])
            z = float(b[36:44])
            cur.append((x, y, z))
        prev_resid = resid
    if cur:
        chains.append(cur)
    return [np.array(c) for c in chains if len(c) >= 2]


def tilt_of_chain(coords):
    """Angle in degrees between the helix long axis and the z axis."""
    c = coords - coords.mean(axis=0)
    _u, _s, vt = np.linalg.svd(c, full_matrices=False)
    axis = vt[0]
    axis = axis / np.linalg.norm(axis)
    cosang = abs(float(axis[2]))            # axis sign is arbitrary, fold to 0-90
    cosang = min(1.0, max(0.0, cosang))
    return np.degrees(np.arccos(cosang))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--core", required=True, help="TM-core itp-local range, e.g. 12-35")
    ap.add_argument("--traj", default="", help="optional xtc/trr for a time average")
    ap.add_argument("--top", default="", help="topology, needed with --traj")
    args = ap.parse_args()
    core_lo, core_hi = parse_core(args.core)

    if args.traj:
        try:
            import mdtraj as md
        except ImportError:
            sys.exit("mdtraj is needed for --traj; install it or omit --traj")
        t = md.load(args.traj, top=args.top if args.top else args.gro)
        # identify TM-core BB atoms per chain from the topology
        top = t.topology
        chains_idx = []
        for chain in top.chains:
            idx = []
            resids = sorted({r.resSeq for r in chain.residues})
            for res in chain.residues:
                if res.name.upper() in AA and core_lo <= res.resSeq <= core_hi:
                    for a in res.atoms:
                        if a.name == "BB":
                            idx.append(a.index)
            if len(idx) >= 2:
                chains_idx.append(idx)
        if not chains_idx:
            sys.exit("no TM-core BB atoms found; check --core range and naming")
        per_chain = [[] for _ in chains_idx]
        for f in range(t.n_frames):
            for ci, idx in enumerate(chains_idx):
                coords = t.xyz[f, idx, :]
                per_chain[ci].append(tilt_of_chain(coords))
        print("TM tilt angle (deg), time-averaged over %d frames:" % t.n_frames)
        allm = []
        for ci, vals in enumerate(per_chain):
            vals = np.array(vals)
            allm.append(vals.mean())
            print("  chain %d: %.1f +/- %.1f" % (ci, vals.mean(), vals.std()))
        allm = np.array(allm)
        print("  all chains: %.1f +/- %.1f" % (allm.mean(), allm.std()))
    else:
        chains = read_gro_backbone(args.gro, core_lo, core_hi)
        if not chains:
            sys.exit("no TM-core BB beads found; check --core range")
        print("TM tilt angle (deg), single frame:")
        tilts = []
        for ci, coords in enumerate(chains):
            tl = tilt_of_chain(coords)
            tilts.append(tl)
            print("  chain %d: %.1f" % (ci, tl))
        tilts = np.array(tilts)
        print("  all chains: %.1f +/- %.1f" % (tilts.mean(), tilts.std()))
        print("compare with the ATR-IR tilt angle for the same composition")


if __name__ == "__main__":
    main()
