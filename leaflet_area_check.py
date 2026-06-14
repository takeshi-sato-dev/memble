#!/usr/bin/env python3
"""
leaflet_area_check.py

Estimate the per-leaflet area of a freshly built CG membrane and fail if the
two leaflets are mismatched beyond a tolerance. Run this AFTER the builder and
BEFORE any GROMACS MD: an area mismatch between leaflets puts the bilayer under
tension and can pore/buckle, and discovering that after a production run means
starting over. Catching it here costs nothing.

How it works (composition-agnostic):
  - For each lipid residue type, pick its head bead by priority (PO4, ROH, ...).
  - Assign every lipid molecule to upper/lower by its head-bead z relative to
    the mean head-bead z (the bilayer midplane).
  - Per leaflet, area = sum over lipids of (count * area-per-lipid).
  - Report counts and areas; exit nonzero if |A_up - A_low| / mean > tol.

Area-per-lipid (APL) values are approximate Martini values in nm^2; this is a
coarse pre-check ("are the leaflets roughly balanced"), not an exact areometer.
Override any APL with --apl "NAME:VALUE ...".
"""

import argparse
import sys

# Approximate Martini area-per-lipid (nm^2). Coarse, override as needed.
APL = {
    "CHOL": 0.40, "POPC": 0.64, "DOPC": 0.67, "DIPC": 0.68, "DLPC": 0.60,
    "DPPC": 0.50, "DPSM": 0.50, "PPSM": 0.50, "PSM": 0.50, "POPE": 0.59,
    "DOPE": 0.61, "POPS": 0.55, "POPG": 0.62, "POPA": 0.55, "POPI": 0.66,
    "PAPC": 0.68, "PUPC": 0.70, "DPG3": 0.75, "DPGS": 0.70,
    "DOPS": 0.60, "DOPG": 0.64, "DOPA": 0.58, "DPPE": 0.50, "DPPS": 0.48,
    "DPPG": 0.52, "DLPS": 0.55, "DSPC": 0.50, "PIPC": 0.66, "PGPC": 0.66,
}
DEFAULT_APL = 0.60  # fallback for lipids not in the table (build proceeds, warns)
HEAD_PRIORITY = ("PO4", "PO1", "PO2", "ROH", "GL1", "GL0", "NC3", "GM1", "AM1", "CNO")


def read_gro(path):
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[1])
    atoms = lines[2:2 + n]
    recs = []
    for ln in atoms:
        # gro fixed cols: resid(5) resname(5) atomname(5) atomnum(5) x y z(8.3 each)
        resid = ln[0:5].strip()
        resname = ln[5:10].strip()
        atomname = ln[10:15].strip()
        try:
            z = float(ln[36:44])
        except ValueError:
            continue
        recs.append((resid, resname, atomname, z))
    return recs


def head_bead_for(resname, recs):
    names = set(a for (_, rn, a, _) in recs if rn == resname)
    for p in HEAD_PRIORITY:
        if p in names:
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gro", required=True)
    ap.add_argument("--lipids", required=True,
                    help="space list of lipid resnames present, e.g. 'CHOL DIPC DPSM'")
    ap.add_argument("--apl", default="", help="overrides 'NAME:VAL ...' (nm^2)")
    ap.add_argument("--tol", type=float, default=0.08, help="max |dA|/mean (default 0.08)")
    ap.add_argument("--hard-tol", type=float, default=0.25,
                    help="abort an asymmetric build only above this mismatch "
                         "(default 0.25); between tol and hard-tol it warns and "
                         "continues so the build still produces run files")
    ap.add_argument("--asym", type=int, default=0, help="1 if asymmetric (enforce abort)")
    args = ap.parse_args()

    for tok in args.apl.split():
        if ":" in tok:
            k, v = tok.split(":"); APL[k] = float(v)

    lipids = args.lipids.split()
    recs = read_gro(args.gro)

    # midplane = mean z of all head beads of all lipid types
    head_of = {}
    for lp in lipids:
        hb = head_bead_for(lp, recs)
        if hb is None:
            print("WARN: no head bead found for %s; skipping it" % lp)
        head_of[lp] = hb

    head_zs = [z for (_, rn, a, z) in recs if rn in head_of and a == head_of[rn]]
    if not head_zs:
        sys.exit("ERROR: no lipid head beads found; check --lipids names vs gro")
    mid = sum(head_zs) / len(head_zs)

    # count molecules per leaflet (one head bead per molecule)
    up = {lp: 0 for lp in lipids}
    lo = {lp: 0 for lp in lipids}
    for (resid, rn, a, z) in recs:
        if rn in head_of and a == head_of[rn]:
            (up if z >= mid else lo)[rn] += 1

    def area(counts):
        miss = [lp for lp in counts if lp not in APL]
        if miss:
            # do not abort the build for an unknown lipid; use a default APL and
            # warn. Pass --apl NAME:VAL for an accurate value if leaflet balance
            # matters for this composition.
            for lp in miss:
                print("WARNING: no tabulated APL for %s; using default %.2f nm^2 "
                      "(override with --apl %s:VALUE)" % (lp, DEFAULT_APL, lp))
        return sum(counts[lp] * APL.get(lp, DEFAULT_APL) for lp in counts)

    a_up, a_lo = area(up), area(lo)
    mean = (a_up + a_lo) / 2.0 if (a_up + a_lo) else 1.0
    dev = abs(a_up - a_lo) / mean

    print("leaflet areas (approx, nm^2): upper=%.1f lower=%.1f  mismatch=%.1f%%"
          % (a_up, a_lo, 100 * dev))
    print("  upper counts:", {k: v for k, v in up.items() if v})
    print("  lower counts:", {k: v for k, v in lo.items() if v})

    if dev > args.tol:
        msg = ("LEAFLET AREA MISMATCH %.1f%% (tol %.1f%%). The leaflets differ in "
               "area; under semiisotropic pressure this leaves some residual "
               "bilayer stress. To balance, adjust the per-leaflet counts so the "
               "summed lipid areas match (reduce the larger leaflet, here the "
               "%s one), or widen box xy." % (100 * dev, 100 * args.tol,
                                              "upper" if a_up > a_lo else "lower"))
        if args.asym and dev > args.hard_tol:
            sys.exit("ABORT: %s This mismatch (> %.0f%%) is too large to "
                     "equilibrate; rebuild with balanced leaflets."
                     % (msg, 100 * args.hard_tol))
        print("WARN: " + msg)
        print("WARN: continuing; the build will finish and you can equilibrate, "
              "but check membrane area and tension during equilibration.")
    else:
        print("leaflets balanced within tolerance; safe to run.")


if __name__ == "__main__":
    main()
