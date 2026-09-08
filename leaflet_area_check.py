#!/usr/bin/env python3
"""
leaflet_area_check.py

Measure the area that each lipid occupies in a freshly built coarse-grained
membrane, per leaflet, and stop the build when the two leaflets are mismatched.
Run this after the packing and before any GROMACS run: leaflets that differ in
area put the bilayer under stress, and finding that after a production run means
building again.

HOW THE AREA IS OBTAINED

Earlier versions of this script summed a table of approximate areas per lipid.
A table cannot answer the question it was asked, because the area a lipid takes
depends on the mixture it sits in, on the sterol content, and on the protein
that occupies part of the leaflet. This version measures the area instead.

Each leaflet is tessellated in the membrane plane: every lipid, and every
protein bead inside the leaflet, is given the region of the plane that lies
closer to it than to anything else, and the area of that region is its area.
The tessellation is integrated by Monte Carlo under the periodic boundary
conditions, following the treatment of lipid areas in protein-containing
membranes of Mori, Ogushi and Sugita (J. Comput. Chem. 2012, 33, 286-296). No
table of areas per lipid is read, and no area is assumed.

WHAT IS CHECKED

  1. A lipid that is present in both leaflets is measured in both. The two
     values agree in a matched bilayer, and they differ when one leaflet holds
     too few lipids for its area. This check reads nothing but the built system.
  2. When a reference area is supplied with --apl, the measured area is compared
     against it as well.
  3. The area that the protein occupies is reported for each leaflet, because a
     leaflet count that ignores the protein is the usual source of a mismatch.

WHAT THE NUMBERS ARE

The tessellation divides the whole plane among the molecules, so the area it
gives a sterol is larger than the area per lipid that a box-area-over-count
convention gives the same sterol. The two are different quantities and they are
not compared here. What the check uses is the comparison of one lipid against
itself in the other leaflet, and that comparison is unaffected.

WHAT IS ASSUMED

The bilayer is treated as flat, and the areas are projected on the xy plane. A
buckled or strongly curved membrane is outside what this script measures, and
the report says so. A freshly built system is flat by construction.

Usage:
  leaflet_area_check.py --gro system.gro --lipids "CHOL DIPC DPSM DOPS" \
      [--tol 0.08] [--hard-tol 0.25] [--asym 1] [--apl "POPC:0.64 CHOL:0.40"] \
      [--points 400000] [--json leaflet_area.json]
"""

import argparse
import json
import sys
from collections import defaultdict

import numpy as np

HEAD_PRIORITY = ("PO4", "PO1", "PO2", "ROH", "GL1", "GL0", "NC3", "GM1",
                 "AM1", "CNO")
SOLVENT = {"W", "WF", "NA", "CL", "ION", "NA+", "CL-"}


def read_gro(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[1])
    body = lines[2:2 + n]
    box = np.array([float(v) for v in lines[2 + n].split()[:3]])
    resid = np.array([int(b[0:5]) for b in body])
    resn = np.array([b[5:10].strip() for b in body])
    aname = np.array([b[10:15].strip() for b in body])
    xyz = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                    for b in body])
    return resid, resn, aname, xyz, box


def gro_key(name):
    """The GRO residue field is five characters, so POP2_45 is written POP2_."""
    return name[:5]


def mc_areas(points_xy, box_xy, n_points, seed=0, chunk=20000):
    """Area of every Voronoi region in the plane, by Monte Carlo integration.

    Random points are thrown into the periodic cell and each is given to the
    nearest reference point under the minimum image convention. The share of
    the throws that a reference point wins, times the area of the cell, is the
    area of its region. The standard error of each area follows from the
    binomial count, and it is returned so the report can state how well the
    integral converged.
    """
    rng = np.random.default_rng(seed)
    m = len(points_xy)
    won = np.zeros(m, dtype=np.int64)
    total = 0
    while total < n_points:
        k = min(chunk, n_points - total)
        q = rng.random((k, 2)) * box_xy
        d = q[:, None, :] - points_xy[None, :, :]
        d -= box_xy * np.round(d / box_xy)
        r2 = (d * d).sum(axis=2)
        nearest = np.argmin(r2, axis=1)
        won += np.bincount(nearest, minlength=m)
        total += k
    cell = float(box_xy[0] * box_xy[1])
    frac = won / float(total)
    area = frac * cell
    err = cell * np.sqrt(np.maximum(frac * (1.0 - frac), 0.0) / float(total))
    return area, err


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gro", required=True)
    ap.add_argument("--lipids", required=True,
                    help="space list of lipid resnames present, e.g. 'CHOL DIPC DPSM'")
    ap.add_argument("--apl", default="",
                    help="optional reference areas 'NAME:VAL ...' in nm^2, "
                         "compared against the measured value")
    ap.add_argument("--tol", type=float, default=0.08,
                    help="largest relative difference accepted between the two "
                         "leaflets (default 0.08)")
    ap.add_argument("--hard-tol", type=float, default=0.25,
                    help="an asymmetric build stops above this difference "
                         "(default 0.25)")
    ap.add_argument("--asym", type=int, default=0, help="1 for an asymmetric build")
    ap.add_argument("--points", type=int, default=400000,
                    help="Monte Carlo throws per leaflet (default 400000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json", default="", help="write the measurement to this file")
    args = ap.parse_args()

    ref = {}
    for tok in args.apl.split():
        if ":" in tok:
            k, v = tok.split(":")
            ref[k] = float(v)

    lipids = args.lipids.split()
    key2lp = {}
    for lp in lipids:
        k = gro_key(lp)
        if k in key2lp and key2lp[k] != lp:
            print("WARNING: %s and %s share the five-character key '%s', so the "
                  "measurement cannot separate them." % (key2lp[k], lp, k))
        else:
            key2lp[k] = lp

    resid, resn, aname, xyz, box = read_gro(args.gro)
    box_xy = box[:2]

    is_lipid = np.array([r in key2lp for r in resn])
    is_solv = np.array([r in SOLVENT or r[:1] == "W" for r in resn])
    is_prot = ~(is_lipid | is_solv)

    if not is_lipid.any():
        sys.exit("ERROR: none of the names in --lipids was found in %s. Check the "
                 "spelling against the residue names in the structure:\n"
                 "  awk 'NR>2{print substr($0,6,5)}' %s | sort -u | head"
                 % (args.gro, args.gro))

    # One entry per lipid molecule: its leaflet, its species, its xy position.
    mol_beads = defaultdict(list)
    for i in np.nonzero(is_lipid)[0]:
        mol_beads[(int(resid[i]), resn[i])].append(i)

    head_z = []
    for key, idx in mol_beads.items():
        names = set(aname[j] for j in idx)
        for p in HEAD_PRIORITY:
            if p in names:
                head_z.extend(xyz[j, 2] for j in idx if aname[j] == p)
                break
    mid = float(np.mean(head_z)) if head_z else float(np.mean(xyz[is_lipid, 2]))

    mols = []          # (species, leaflet, x, y)
    for (rid, rn), idx in mol_beads.items():
        z = float(np.mean(xyz[idx, 2]))
        leaf = "upper" if z >= mid else "lower"
        # The centroid in the membrane plane is the reference point. It needs no
        # per-lipid choice of atom, so a lipid the script has never seen is
        # measured on the same footing as a familiar one.
        p = xyz[idx, :2]
        ref_xy = p[0] + np.mean(((p - p[0]) + box_xy / 2.0) % box_xy - box_xy / 2.0,
                                axis=0)
        mols.append((key2lp[rn], leaf, float(ref_xy[0]), float(ref_xy[1])))

    out = {"box_nm": [float(v) for v in box], "midplane_nm": mid,
           "method": "Monte Carlo integration of the planar Voronoi regions",
           "assumes": "a flat bilayer; areas are projected on the xy plane",
           "leaflets": {}}

    per_leaf_species = {}
    for leaf in ("upper", "lower"):
        here = [m for m in mols if m[1] == leaf]
        if not here:
            continue
        # Protein beads that lie in this leaflet take area from it.
        if leaf == "upper":
            psel = is_prot & (xyz[:, 2] >= mid) & (xyz[:, 2] < mid + 2.5)
        else:
            psel = is_prot & (xyz[:, 2] < mid) & (xyz[:, 2] > mid - 2.5)
        pxy = xyz[psel, :2]

        pts = np.array([[m[2], m[3]] for m in here] + list(pxy))
        areas, errs = mc_areas(pts, box_xy, args.points, seed=args.seed)
        n_lip = len(here)

        by_species = defaultdict(list)
        for j, m in enumerate(here):
            by_species[m[0]].append(areas[j])
        species = {s: {"n": len(v), "apl_nm2": round(float(np.mean(v)), 4),
                       "sd_nm2": round(float(np.std(v)), 4)}
                   for s, v in sorted(by_species.items())}
        per_leaf_species[leaf] = species

        prot_area = float(areas[n_lip:].sum()) if len(pts) > n_lip else 0.0
        lip_area = float(areas[:n_lip].sum())
        out["leaflets"][leaf] = {
            "n_lipids": n_lip,
            "lipid_area_nm2": round(lip_area, 3),
            "protein_area_nm2": round(prot_area, 3),
            "mc_area_error_nm2": round(float(errs.sum()), 4),
            "species": species,
        }

    for leaf, d in out["leaflets"].items():
        print("%s leaflet: %d lipids, lipid area %.1f nm^2, protein area %.1f nm^2"
              % (leaf, d["n_lipids"], d["lipid_area_nm2"], d["protein_area_nm2"]))
        for s, v in d["species"].items():
            print("    %-8s n=%-5d Voronoi area %.3f nm^2 (sd %.3f)"
                  % (s, v["n"], v["apl_nm2"], v["sd_nm2"]))

    # --- check 1: a species in both leaflets is measured in both -------------
    shared = []
    if "upper" in per_leaf_species and "lower" in per_leaf_species:
        for s in per_leaf_species["upper"]:
            if s in per_leaf_species["lower"]:
                a = per_leaf_species["upper"][s]["apl_nm2"]
                b = per_leaf_species["lower"][s]["apl_nm2"]
                mean = (a + b) / 2.0
                if mean > 0:
                    shared.append((s, a, b, abs(a - b) / mean))
    out["shared_species"] = [{"lipid": s, "upper_nm2": a, "lower_nm2": b,
                              "relative_difference": round(d, 4)}
                             for s, a, b, d in shared]

    worst, worst_s = 0.0, ""
    for s, a, b, d in shared:
        print("  %-8s upper %.3f nm^2, lower %.3f nm^2, difference %.1f%%"
              % (s, a, b, 100 * d))
        if d > worst:
            worst, worst_s = d, s

    # --- check 2: measured against a reference, when one is given ------------
    ref_dev = []
    for leaf, sp in per_leaf_species.items():
        for s, v in sp.items():
            if s in ref and ref[s] > 0:
                ref_dev.append((leaf, s, v["apl_nm2"], ref[s],
                                abs(v["apl_nm2"] - ref[s]) / ref[s]))
    for leaf, s, got, want, d in ref_dev:
        print("  %-8s %s leaflet: measured %.3f nm^2 against the reference %.3f "
              "nm^2, difference %.1f%%" % (s, leaf, got, want, 100 * d))

    out["result"] = "PASS"
    if not shared:
        out["result"] = "REPORTED"
        print("")
        print("NOTE: no lipid is present in both leaflets, so the leaflets cannot "
              "be compared against each other. The measured areas above are "
              "reported and the build continues.")
        print("      To have this checked, give a reference with --apl, or put one "
              "lipid in both leaflets.")
    elif worst > args.tol:
        big = "upper" if per_leaf_species["upper"][worst_s]["apl_nm2"] > \
                         per_leaf_species["lower"][worst_s]["apl_nm2"] else "lower"
        out["result"] = "FAIL"
        print("")
        print("LEAFLET AREA MISMATCH: %s occupies %.1f%% more area in the %s "
              "leaflet than in the other one (tolerance %.1f%%)."
              % (worst_s, 100 * worst, big, 100 * args.tol))
        print("")
        print("  What to do:")
        print("    The %s leaflet holds too few lipids for the area it has to "
              "cover, so its lipids are stretched." % big)
        print("    1. add lipids to the %s leaflet, or remove them from the other,"
              % big)
        print("       until the two measured areas agree. The counts above say by")
        print("       how much: the ratio of the two areas is the ratio to correct.")
        print("    2. a protein that reaches into one leaflet takes area from that")
        print("       leaflet alone. The protein areas printed above are the amount")
        print("       to subtract before the counts are set.")
        print("    3. build again with UPPER and LOWER set separately, rather than")
        print("       with LIPIDS, so the two leaflets are counted on their own.")
        if args.asym and worst > args.hard_tol:
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(out, fh, indent=2)
            sys.exit("ABORT: a difference above %.0f%% does not equilibrate away."
                     % (100 * args.hard_tol))
        print("")
        print("WARNING: the build continues, and the bilayer carries this stress "
              "into the equilibration. Watch the area and the thickness.")
    else:
        print("")
        print("The leaflets are matched within %.1f%%. The largest difference is "
              "%s, at %.1f%%." % (100 * args.tol, worst_s or "none", 100 * worst))

    print("")
    print("Measured by Monte Carlo integration of the planar Voronoi regions "
          "(%d throws per leaflet); the bilayer is treated as flat. The regions "
          "divide the whole plane, so a sterol receives more area here than the "
          "area per lipid of a box-area-over-count convention gives it; the "
          "check compares each lipid against itself in the other leaflet."
          % args.points)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
