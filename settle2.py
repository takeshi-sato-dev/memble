#!/usr/bin/env python3
"""Where the cholesterol of one finished run settled.

The rule is the rule of Section 2.6: a molecule belongs to the leaflet its own
head bead sits toward, taken against the rest of the beads of that same
molecule. No plane is drawn through the bilayer, because the protein bends the
membrane and a plane then misassigns the lipids beside it.

  cholesterol   head = ROH
  phospholipid  head = PO4, or the first phosphate the molecule carries

The numbers the build assigned are read from memble_build.json and are never
counted from the packed system, for the reason Section 3.3 gives: counting
misreads a few molecules wherever the protein bends the membrane.

Usage:
    python3 settle.py ~/counts_run/hold/hold_work
    python3 settle.py ~/counts_run/bracket_work --out relax_bracket.json

Prints what the build assigned, what the run carried over the last 40%, and
both in the mol% of each leaflet, which is the unit a composition is stated in.
"""
import argparse
import json
import os
import re

import numpy as np

STEROL = "CHOL"
STEROL_HEAD = ("ROH",)
PL_HEAD = ("PO4", "P1", "P2", "P3", "CNO", "PO1", "PO2")
SOLVENT = ("W", "WF", "ION", "NA", "CL", "NA+", "CL-")


def as_counts(v):
    """The build record gives a leaflet either as a dict or as its spec string.

    memble writes UPPER and LOWER back as the strings it was given, for example
    "CHOL:76 DLPC:70 PSM:70". Both forms are accepted so that a build record
    from any version is read the same way.
    """
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


def lookup(counts, name):
    """Match a residue name against a build key, allowing a truncated name."""
    if name in counts:
        return counts[name]
    for k, v in counts.items():
        if k.startswith(name) or name.startswith(k):
            return v
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work", help="the *_work directory of the build")
    ap.add_argument("--structure", default=None)
    ap.add_argument("--trajectory", default=None)
    ap.add_argument("--last-fraction", type=float, default=0.4)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    import MDAnalysis as mda

    def pick(names):
        for n in names:
            p = os.path.join(a.work, n)
            if os.path.exists(p):
                return p
        return None

    struct = a.structure or pick(("step7_production.gro", "prod.gro", "md.gro"))
    traj = a.trajectory or pick(("step7_production.xtc", "prod.xtc", "md.xtc"))
    if not struct or not traj:
        raise SystemExit("settle: structure or trajectory not found under %s" % a.work)
    print("settle: %s" % os.path.basename(struct))
    print("settle: %s" % os.path.basename(traj))

    u = mda.Universe(struct, traj)
    lipid = u.select_atoms("not protein and not resname %s" % " ".join(SOLVENT))

    mols = {}
    for r in lipid.residues:
        heads = STEROL_HEAD if r.resname == STEROL else PL_HEAD
        h = r.atoms.select_atoms("name %s" % " ".join(heads))
        body = [at.index for at in r.atoms if at.name not in heads]
        if not len(h) or not body:
            continue
        mols.setdefault(r.resname, []).append((h.indices[0], body))
    if STEROL not in mols:
        raise SystemExit("settle: no %s found" % STEROL)
    for k in sorted(mols):
        print("  %-8s %d molecules" % (k, len(mols[k])))

    n = len(u.trajectory)
    start = int((1.0 - a.last_fraction) * n)
    counts = {k: {"upper": [], "lower": []} for k in mols}
    first = {}
    for i, ts in enumerate(u.trajectory):
        for k, ms in mols.items():
            hz = np.array([u.atoms.positions[h, 2] for h, _ in ms])
            bz = np.array([u.atoms.positions[b, 2].mean() for _, b in ms])
            up = int((hz > bz).sum())
            counts[k]["upper"].append(up)
            counts[k]["lower"].append(len(ms) - up)
            if i == 0:
                first[k] = (up, len(ms) - up)

    built = None
    mb = os.path.join(a.work, "memble_build.json")
    if os.path.exists(mb):
        d = json.load(open(mb))
        for key in ("composition", "leaflets", "lipids", None):
            cand = d.get(key) if key else d
            if isinstance(cand, dict) and "upper" in cand and "lower" in cand:
                built = {"upper": as_counts(cand["upper"]),
                         "lower": as_counts(cand["lower"])}
                break
    if built is None:
        raise SystemExit("settle: memble_build.json not found or has no "
                         "upper/lower block under %s" % a.work)

    print("")
    print("%-8s %-22s %-22s %-10s" % ("", "built upper/lower", "run, last 40%", "moved"))
    for k in sorted(mols):
        bu = lookup(built["upper"], k)
        bl = lookup(built["lower"], k)
        su = float(np.mean(counts[k]["upper"][start:]))
        sl = float(np.mean(counts[k]["lower"][start:]))
        print("%-8s %-22s %-22s %+8.1f" % (k, "%d / %d" % (bu, bl),
                                           "%.1f / %.1f" % (su, sl), su - bu))

    pu = sum(v for kk, v in built["upper"].items() if kk != STEROL)
    pl = sum(v for kk, v in built["lower"].items() if kk != STEROL)
    cu = lookup(built["upper"], STEROL)
    cl = lookup(built["lower"], STEROL)
    su = float(np.mean(counts[STEROL]["upper"][start:]))
    sl = float(np.mean(counts[STEROL]["lower"][start:]))
    tot = cu + cl

    print("")
    print("  phospholipids            %d upper, %d lower" % (pu, pl))
    print("  imbalance f              %+.2f%%" % (100.0 * (pu - pl) / (pu + pl)))
    print("")
    print("  cholesterol built        %d / %d  = %.2f%% above" % (cu, cl, 100.0 * cu / tot))
    print("  cholesterol settled      %.1f / %.1f  = %.2f%% above" % (su, sl, 100.0 * su / (su + sl)))
    print("  moved                    %+.2f percentage points" % (100.0 * su / (su + sl) - 100.0 * cu / tot))
    print("")
    print("  upper leaflet cholesterol  built %.2f mol%%   run %.2f mol%%"
          % (100.0 * cu / (pu + cu), 100.0 * su / (pu + su)))
    print("  lower leaflet cholesterol  built %.2f mol%%   run %.2f mol%%"
          % (100.0 * cl / (pl + cl), 100.0 * sl / (pl + sl)))

    out = a.out or os.path.join(a.work, "relax.json")
    json.dump({"counts": counts, "built_counts": built,
               "last_fraction": a.last_fraction, "n_frames": n,
               "structure": struct, "trajectory": traj},
              open(out, "w"))
    print("")
    print("written to %s" % out)


if __name__ == "__main__":
    main()
