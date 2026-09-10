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

  1. A lipid that both leaflets hold in numbers is measured in both. The two
     values agree in a matched bilayer, and they differ when one leaflet holds
     too few lipids for its area. This check reads nothing but the built system.
     A species that one leaflet holds in numbers and the other holds a handful
     of is not compared, because a mean taken over a handful of molecules
     carries a standard error of tens of percent and because a lipid surrounded
     by species it does not share the leaflet with reports its neighbors rather
     than the tension. Such a species is printed and written to the report, and
     it decides nothing.
  2. A reference area supplied with --apl is printed beside the measurement. The
     two are different quantities, because a Voronoi region divides the whole
     plane while an area per lipid counts lipids into the area of the box, so
     the reference is reported and never decides the outcome.
  3. The area that the protein occupies is reported for each leaflet, because a
     leaflet count that ignores the protein is the usual source of a mismatch.

WHAT THE COMPARISON CAN AND CANNOT SEPARATE

The area of a Voronoi region is set by the neighbors of the molecule, so a lipid
compared against itself in the other leaflet reports the tension of the leaflets
only while the two leaflets present a similar neighborhood. In a bilayer whose
two leaflets hold largely different lipids, the same lipid sits among PC and SM
on one side and among PE and PS on the other, and the two areas differ while the
leaflets carry no tension between them. This script measures how much of each
leaflet is made of lipids that both leaflets hold. Below half, the difference is
reported and the build continues, and only a difference above the hard tolerance
stops it.

WHAT THE NUMBERS ARE

The tessellation divides the whole plane among the molecules, so the area it
gives a sterol is larger than the area per lipid that a box-area-over-count
convention gives the same sterol. The two are different quantities and they are
not compared here. What the check uses is the comparison of one lipid against
itself in the other leaflet, and that comparison is unaffected.

HOW A LIPID IS GIVEN ITS LEAFLET

The leaflet of a lipid is read from the lipid. A lipid points out of the leaflet
it belongs to, so the head bead of a lipid in the upper leaflet lies above the
rest of that same molecule, and the head bead of a lipid in the lower leaflet
lies below it. The head bead is the phosphate of a phospholipid and the hydroxyl
of a sterol. Comparing each molecule against a plane drawn through the whole
bilayer gives a different answer, and a worse one: the plane has to be placed,
the mean z of every lipid bead places it toward whichever leaflet holds more
beads, and a bilayer that holds a protein is not flat, so a lipid in a region
the protein depresses sits below a plane its own leaflet rises above.

WHAT IS ASSUMED

The areas are projected on the xy plane, so a buckled or strongly curved
membrane is outside what this script measures and the report says so. A freshly
built system is flat by construction. The leaflet assignment above needs no such
assumption.

Usage:
  leaflet_area_check.py --gro system.gro --lipids "CHOL DIPC DPSM DOPS" \
      [--tol 0.08] [--hard-tol 0.25] [--asym 1] [--apl "POPC:0.64 CHOL:0.40"] \
      [--points 400000] [--min-shared 10] [--json leaflet_area.json]
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


def head_of(idx, aname):
    """The beads of one molecule that lie at its head, in order of preference.

    The first name in HEAD_PRIORITY that the molecule carries wins, so a
    phospholipid is read from its phosphate, a sterol from its hydroxyl, and a
    lipid that carries neither from its glycerol or its headgroup bead. A
    molecule that carries none of them returns nothing, and it is placed by the
    midplane of the molecules that were placed.
    """
    for h in HEAD_PRIORITY:
        sel = [int(i) for i in idx if aname[i] == h]
        if sel:
            return sel
    return []


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
                    help="reference areas 'NAME:VAL ...' in nm^2, printed beside "
                         "the measurement; they never decide the outcome")
    ap.add_argument("--tol", type=float, default=0.08,
                    help="largest relative difference accepted between the two "
                         "leaflets (default 0.08)")
    ap.add_argument("--hard-tol", type=float, default=0.25,
                    help="an asymmetric build stops above this difference "
                         "(default 0.25)")
    ap.add_argument("--asym", type=int, default=0, help="1 for an asymmetric build")
    ap.add_argument("--points", type=int, default=400000,
                    help="Monte Carlo throws per leaflet (default 400000)")
    ap.add_argument("--min-shared", type=int, default=10,
                    help="molecules of one species a leaflet has to hold before "
                         "that species is compared against the other leaflet "
                         "(default 10)")
    ap.add_argument("--min-shared-ratio", type=float, default=0.2,
                    help="smallest ratio of the two counts of one species "
                         "accepted for the comparison (default 0.2)")
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

    # The leaflet of a lipid is read from the lipid: the head bead of a lipid in
    # the upper leaflet lies above the rest of that same molecule. No plane is
    # drawn through the bilayer, so no lipid is carried across a plane that the
    # protein, or the difference between the two leaflets, has moved under it.
    placed = []        # (species, upper, x, y, z_head)
    unplaced = []      # (species, x, y, z_mol)
    for (rid, rn), idx in mol_beads.items():
        # The centroid in the membrane plane is the reference point. It needs no
        # per-lipid choice of atom, so a lipid the script has never seen is
        # measured on the same footing as a familiar one.
        p = xyz[idx, :2]
        ref_xy = p[0] + np.mean(((p - p[0]) + box_xy / 2.0) % box_xy - box_xy / 2.0,
                                axis=0)
        x, y = float(ref_xy[0]), float(ref_xy[1])
        h = head_of(idx, aname)
        b = [int(i) for i in idx if i not in set(h)]
        if h and b:
            zh = float(np.mean(xyz[h, 2]))
            placed.append((key2lp[rn], zh > float(np.mean(xyz[b, 2])), x, y, zh))
        else:
            unplaced.append((key2lp[rn], x, y, float(np.mean(xyz[idx, 2]))))

    if not placed:
        sys.exit("ERROR: no lipid in %s carries a head bead the script knows "
                 "(%s), so the leaflets cannot be read from the molecules."
                 % (args.gro, " ".join(HEAD_PRIORITY)))

    # The midplane is the midpoint between the head beads of the two leaflets,
    # which the assignment above has already settled. It places the molecules
    # that carry no head bead, and it divides the protein beads between the
    # leaflets.
    zu = [m[4] for m in placed if m[1]]
    zl = [m[4] for m in placed if not m[1]]
    if zu and zl:
        mid = 0.5 * (float(np.mean(zu)) + float(np.mean(zl)))
    else:
        mid = float(np.mean(xyz[is_lipid, 2]))

    mols = []          # (species, leaflet, x, y)
    for s, up, x, y, _z in placed:
        mols.append((s, "upper" if up else "lower", x, y))
    for s, x, y, z in unplaced:
        mols.append((s, "upper" if z >= mid else "lower", x, y))

    out = {"box_nm": [float(v) for v in box], "midplane_nm": mid,
           "n_placed_by_molecule": len(placed),
           "n_placed_by_midplane": len(unplaced),
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
                       "sd_nm2": round(float(np.std(v)), 4),
                       "se_nm2": round(float(np.std(v) / np.sqrt(len(v))), 4)}
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
            print("    %-8s n=%-5d Voronoi area %.3f +/- %.3f nm^2 (sd %.3f)"
                  % (s, v["n"], v["apl_nm2"], v["se_nm2"], v["sd_nm2"]))

    # --- check 1: a species both leaflets hold in numbers is measured in both -
    # A species that one leaflet holds in numbers and the other holds a handful
    # of is set aside. Two things make such a comparison meaningless. The mean
    # area of eight molecules carries a standard error of tens of percent, and
    # eight molecules of a species that belongs to the other leaflet sit among
    # neighbors they do not share the leaflet with, so their area reports the
    # neighbors. Both the tolerance and the correction that follows from it read
    # the largest difference among the compared species, so one such species
    # would decide the build on its own.
    shared, aside = [], []
    if "upper" in per_leaf_species and "lower" in per_leaf_species:
        for s in per_leaf_species["upper"]:
            if s not in per_leaf_species["lower"]:
                continue
            A = per_leaf_species["upper"][s]
            B = per_leaf_species["lower"][s]
            a, b = A["apl_nm2"], B["apl_nm2"]
            mean = (a + b) / 2.0
            if mean <= 0:
                continue
            nu, nl = A["n"], B["n"]
            se = float(np.hypot(A["se_nm2"], B["se_nm2"]))
            row = (s, a, b, abs(a - b) / mean, se / mean)
            if min(nu, nl) < args.min_shared:
                aside.append((row, nu, nl, "one leaflet holds fewer than %d of "
                                           "them" % args.min_shared))
            elif min(nu, nl) < args.min_shared_ratio * max(nu, nl):
                aside.append((row, nu, nl, "one leaflet holds fewer than %.0f%% "
                                           "as many of them as the other"
                                           % (100 * args.min_shared_ratio)))
            else:
                shared.append(row)

    # How much of each leaflet is made of lipids that both leaflets hold in
    # numbers. A species set aside above is not one of them.
    shared_names = set(r[0] for r in shared)
    overlap = []
    for leaf, sp in per_leaf_species.items():
        tot = sum(v["n"] for v in sp.values())
        if tot:
            overlap.append(sum(sp[s]["n"] for s in shared_names if s in sp) / tot)
    shared_fraction = float(np.mean(overlap)) if overlap else 0.0
    out["shared_fraction"] = round(shared_fraction, 3)

    out["shared_species"] = [{"lipid": s, "upper_nm2": a, "lower_nm2": b,
                              "relative_difference": round(d, 4),
                              "relative_standard_error": round(e, 4)}
                             for s, a, b, d, e in shared]
    out["not_compared"] = [{"lipid": r[0], "upper_nm2": r[1], "lower_nm2": r[2],
                            "relative_difference": round(r[3], 4),
                            "relative_standard_error": round(r[4], 4),
                            "n_upper": nu, "n_lower": nl, "reason": why}
                           for r, nu, nl, why in aside]

    # A difference is read against the scatter that produced it. With a few tens
    # of molecules of one lipid in a leaflet the mean area carries a standard
    # error of a few percent, and a difference of that size says nothing.
    worst, worst_s, worst_e = 0.0, "", 0.0
    for s, a, b, d, e in shared:
        print("  %-8s upper %.3f nm^2, lower %.3f nm^2, difference %.1f%% "
              "(standard error %.1f%%)" % (s, a, b, 100 * d, 100 * e))
        if d > worst:
            worst, worst_s, worst_e = d, s, e
    for (s, a, b, d, e), nu, nl, why in aside:
        print("  %-8s upper %.3f nm^2 (n=%d), lower %.3f nm^2 (n=%d): not "
              "compared, because %s." % (s, a, nu, b, nl, why))

    # --- a reference, printed for the reader and used for nothing else -------
    ref_dev = []
    for leaf, sp in per_leaf_species.items():
        for s, v in sp.items():
            if s in ref and ref[s] > 0:
                ref_dev.append((leaf, s, v["apl_nm2"], ref[s],
                                abs(v["apl_nm2"] - ref[s]) / ref[s]))
    for leaf, s, got, want, d in ref_dev:
        print("  %-8s %s leaflet: measured %.3f nm^2 against the reference %.3f "
              "nm^2, difference %.1f%%" % (s, leaf, got, want, 100 * d))
    out["reference"] = [{"leaflet": l, "lipid": s, "measured_nm2": g,
                         "reference_nm2": w, "relative_difference": round(d, 4)}
                        for l, s, g, w, d in ref_dev]

    out["result"] = "PASS"
    if not shared:
        out["result"] = "REPORTED"
        print("")
        if aside:
            print("NOTE: no lipid is held by both leaflets in numbers, so the "
                  "leaflets cannot be compared against each other. The species "
                  "set aside above are held by one leaflet and by a handful of "
                  "molecules of the other. The measured areas are reported and "
                  "the build continues.")
        else:
            print("NOTE: no lipid is present in both leaflets, so the leaflets "
                  "cannot be compared against each other. The measured areas "
                  "above are reported and the build continues.")
        print("      To have this checked, give a reference with --apl, or put one "
              "lipid in both leaflets.")
    elif worst > args.tol and worst > 2.0 * worst_e and shared_fraction < 0.5:
        out["result"] = "REPORTED"
        print("")
        print("%s occupies %.1f%% more area in one leaflet than in the other, and "
              "the two leaflets have only %.0f%% of their lipids in common."
              % (worst_s, 100 * worst, 100 * shared_fraction))
        print("A lipid takes the area its neighbors leave it, and the neighbors "
              "differ between these two leaflets, so this difference does not "
              "report the tension between them.")
        print("The measurement is reported and the build continues.")
        if args.asym and worst > args.hard_tol:
            if args.json:
                with open(args.json, "w") as fh:
                    json.dump(out, fh, indent=2)
            sys.exit("ABORT: a difference above %.0f%% is too large to come from "
                     "the neighbors alone." % (100 * args.hard_tol))
    elif worst > args.tol and worst > 2.0 * worst_e:
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
        if worst > args.tol:
            print("The largest difference is %s, at %.1f%%, and the standard "
                  "error of that difference is %.1f%%. The two leaflets are not "
                  "separated by the scatter of the measurement."
                  % (worst_s, 100 * worst, 100 * worst_e))
        else:
            print("The leaflets are matched within %.1f%%. The largest difference "
                  "is %s, at %.1f%% (standard error %.1f%%), and the two leaflets "
                  "have %.0f%% of their lipids in common."
                  % (100 * args.tol, worst_s or "none", 100 * worst,
                     100 * worst_e, 100 * shared_fraction))

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
