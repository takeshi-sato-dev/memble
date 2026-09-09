#!/usr/bin/env python3
"""
verify_system.py

Final gate of a memble build. Reads the finished system and confirms every
post-condition that a membrane protein system has to satisfy before it is worth
running. Writes memble_report.json and memble_report.txt, and exits non-zero if
any check fails, so a build that did not pass is never reported as a success.

The checks are deliberately the ones that a broken build still survives: a
system that fails any of them starts, minimizes and runs, and returns numbers
that look ordinary.

  1. secondary structure   the string handed to martinize2 covers every residue
  2. composition           the lipids placed match the composition requested
  3. protein placement     the transmembrane segment sits at the bilayer midplane
  4. charge                the system is neutral
  5. periodic image        the protein and its image in z are further apart than
                           the non-bonded cutoff
  6. water layer           the water above and below reaches the requested depth
  7. overlap               no two beads of different molecules are closer than
                           the minimum separation
  8. leaflet area          one lipid occupies the same area in both leaflets

Usage:
  verify_system.py --gro system.gro --top system.top --itp-dir . \
      --lipids "CHOL DLPC DPSM DOPS POPC" \
      [--ss-string HHHH...] [--ss-mode dssp] [--tm-resids 619-641] \
      [--expect-upper "CHOL:1 DLPC:1"] [--expect-lower "POPC:3 POPE:1"] \
      [--water-nm 2.5] [--cutoff 1.1] [--min-dist 0.12] \
      [--out-json memble_report.json] [--out-txt memble_report.txt]

Every threshold is an option, and --allow NAME skips one named check and records
the override in the report.
"""

import argparse
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

SOLVENT = {"W", "WF", "NA", "CL", "ION", "NA+", "CL-", "TIP3", "SOL"}
HEAD_PRIORITY = ("PO4", "PO1", "PO2", "ROH", "GL1", "GL0", "NC3", "GM1",
                 "AM1", "CNO")
CHECKS = ("secondary_structure", "composition", "protein_placement", "charge",
          "periodic_image", "water_layer", "overlap", "leaflet_area")

# What to do when a check fails. A gate that only says "FAIL" costs the user the
# afternoon that the gate was meant to save.
REMEDY = {
    "secondary_structure":
        "Build again with SS_MODE=dssp (the default), which runs DSSP in memble\n"
        "and records the string it hands to martinize2.\n"
        "If the assignment is given by hand, check that SS_OVERRIDE has one\n"
        "character per residue of the protein as martinize2 coarse-grained it,\n"
        "counting every chain and every copy:\n"
        "  awk '/^\\[ *atoms/{a=1;next} /^\\[/{a=0} a&&NF{print $3}' molecule_0.itp | sort -un | wc -l",
    "composition":
        "The lipids that were placed differ from the ratio that was asked for.\n"
        "  1. a lipid the packing could not fit usually means the box is too small:\n"
        "     raise BOX_X and BOX_Y by 1 nm and build again\n"
        "  2. check the spelling of every name in LIPIDS, UPPER and LOWER against\n"
        "     the lipidome:  ls $M3_DIR/*.itp\n"
        "  3. a large protein takes area out of one leaflet, so an asymmetric\n"
        "     composition needs UPPER and LOWER set separately, not LIPIDS\n"
        "  4. to accept the placed composition:  export MEMBLE_ALLOW=composition",
    "protein_placement":
        "The protein is not centered in the bilayer, so part of the transmembrane\n"
        "segment sits in water.\n"
        "  1. check that TM_RANGE names the residues that actually cross the\n"
        "     membrane, in the numbering of the input PDB\n"
        "  2. move the protein through the bilayer with Z_SHIFT, in nm:\n"
        "       export Z_SHIFT=0.4     (toward the upper leaflet)\n"
        "       export Z_SHIFT=-0.4    (toward the lower leaflet)\n"
        "  3. a protein oriented upside down is fixed with NTERM_SIDE=up or down",
    "charge":
        "The system is not neutral, so GROMACS applies a uniform background charge\n"
        "and every electrostatic number is shifted.\n"
        "  1. check that SALT_M is set, for example SALT_M=0.15\n"
        "  2. check that the ion itp is in M3_DIR:  ls $M3_DIR | grep -i ion\n"
        "  3. a charge that is not a whole number means an itp was read wrongly;\n"
        "     the moleculetype named in the message above is the one to look at",
    "periodic_image":
        "The protein is close enough to its own image in z to interact with it.\n"
        "  1. raise the water layer, which raises the box:  export WATER_NM=3.0\n"
        "  2. or set the box directly:  export BOX_Z=14.0\n"
        "  3. if the span looks far larger than the protein, a chain is still split\n"
        "     across the boundary; check that make_protein_whole ran without a\n"
        "     warning above",
    "water_layer":
        "The water above or below the bilayer is thinner than WATER_NM asked for.\n"
        "  1. raise the box in z:  export BOX_Z=<the value memble reported + 2>\n"
        "  2. or lower the request:  export WATER_NM=2.0\n"
        "  3. a tall protein in a small xy box leaves no room; raise BOX_X and BOX_Y",
    "leaflet_area":
        "One lipid occupies a different area in the two leaflets, so the leaflet\n"
        "with the larger area holds too few lipids and its lipids are stretched.\n"
        "  1. change the ratio of that lipid between UPPER and LOWER by the ratio\n"
        "     of the two measured areas, and build again\n"
        "  2. a protein that reaches into one leaflet takes area from that leaflet\n"
        "     alone; leaflet_area.json reports how much\n"
        "  3. to accept the leaflets as built:  export MEMBLE_ALLOW=leaflet_area",
    "overlap":
        "Two beads of different molecules are close enough to give an infinite\n"
        "force at minimization. The pair is named above.\n"
        "  1. give the packing room:    raise BOX_X and BOX_Y by 1 nm\n"
        "  2. lower the lipid density:  raise COBY_APL\n"
        "  3. for several copies:       raise SPACING_NM, or lower N_COPY\n"
        "  4. to keep this system and look at it:  export MEMBLE_ALLOW=overlap",
}


# ----------------------------------------------------------------- readers

def read_gro(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[1])
    body = lines[2:2 + n]
    box = [float(v) for v in lines[2 + n].split()[:3]]
    resid = np.array([int(b[0:5]) for b in body])
    resn = [b[5:10].strip() for b in body]
    aname = [b[10:15].strip() for b in body]
    xyz = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                    for b in body])
    return resid, resn, aname, xyz, np.array(box)


def top_include_dirs(path):
    """Directories named by the #include lines of a topology.

    The ion and water topologies usually sit in the lipidome directory while the
    protein topologies sit next to the system, so the charge cannot be summed
    from one directory alone.
    """
    out, base = [], os.path.dirname(os.path.abspath(path))
    try:
        fh = open(path)
    except OSError:
        return out
    with fh:
        for raw in fh:
            line = raw.strip()
            if not line.startswith("#include"):
                continue
            parts = line.split('"')
            if len(parts) < 2:
                continue
            inc = parts[1]
            d = os.path.dirname(inc if os.path.isabs(inc)
                                else os.path.join(base, inc))
            if d and d not in out:
                out.append(d)
    return out


def read_top_molecules(path):
    """Return [(moleculetype, count), ...] from the [ molecules ] section."""
    out, in_mol = [], False
    with open(path) as fh:
        for raw in fh:
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("["):
                in_mol = line.replace(" ", "").lower() == "[molecules]"
                continue
            if in_mol:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        out.append((parts[0], int(parts[1])))
                    except ValueError:
                        pass
    return out


def itp_charges(itp_dirs):
    """Map moleculetype name -> (total charge, bead count) over every itp found."""
    charge, beads = {}, {}
    seen = set()
    for d in itp_dirs:
        if not d or not os.path.isdir(d):
            continue
        for path in sorted(glob.glob(os.path.join(d, "*.itp"))):
            real = os.path.realpath(path)
            if real in seen:
                continue
            seen.add(real)
            name, section, q, n = None, None, 0.0, 0
            try:
                fh = open(path)
            except OSError:
                continue
            with fh:
                for raw in fh:
                    line = raw.split(";")[0].strip()
                    if not line:
                        continue
                    if line.startswith("["):
                        if section == "atoms" and name is not None:
                            charge.setdefault(name, q)
                            beads.setdefault(name, n)
                        section = line.strip("[] ").lower()
                        if section == "moleculetype":
                            name, q, n = None, 0.0, 0
                        continue
                    if section == "moleculetype" and name is None:
                        name = line.split()[0]
                    elif section == "atoms":
                        parts = line.split()
                        if len(parts) >= 7:
                            try:
                                q += float(parts[6])
                                n += 1
                            except ValueError:
                                pass
            if section == "atoms" and name is not None:
                charge.setdefault(name, q)
                beads.setdefault(name, n)
    return charge, beads


def parse_ratio(spec):
    """'CHOL:1 DLPC:2' -> {'CHOL': 1.0, 'DLPC': 2.0}"""
    out = {}
    for token in (spec or "").split():
        parts = token.split(":")
        if len(parts) >= 2:
            try:
                out[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return out


def parse_resids(spec):
    """'619-641' or '619:641' or '10,12,14' -> set of ints"""
    out = set()
    for chunk in (spec or "").replace(":", "-").split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            a, b = chunk.split("-")[:2]
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(chunk))
    return out


# ----------------------------------------------------------------- checks

class Report:
    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail, measured=None, expected=None,
            skipped=False, remedy=""):
        self.rows.append({"check": name, "pass": bool(ok), "skipped": skipped,
                          "measured": measured, "expected": expected,
                          "detail": detail,
                          "remedy": "" if (ok or skipped) else
                                    (remedy or REMEDY.get(name, ""))})

    def failed(self):
        return [r for r in self.rows if not r["pass"] and not r["skipped"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--itp-dir", default=None, action="append",
                    help="directory holding itp files; repeat for several. The "
                         "directories named by #include lines in the topology "
                         "are read as well.")
    ap.add_argument("--lipids", default="")
    ap.add_argument("--ss-string", default="")
    ap.add_argument("--ss-mode", default="")
    ap.add_argument("--tm-resids", default="")
    ap.add_argument("--expect-upper", default="")
    ap.add_argument("--expect-lower", default="")
    ap.add_argument("--water-nm", type=float, default=0.0)
    ap.add_argument("--cutoff", type=float, default=1.1)
    ap.add_argument("--min-dist", type=float, default=0.12)
    ap.add_argument("--placement-tol", type=float, default=0.5)
    ap.add_argument("--composition-tol", type=float, default=0.10)
    ap.add_argument("--charge-tol", type=float, default=1e-3)
    ap.add_argument("--allow", action="append", default=[],
                    help="skip one named check (repeatable)")
    ap.add_argument("--meta", default="", help="JSON file of build provenance")
    ap.add_argument("--leaflet-json", default="",
                    help="leaflet_area.json written by leaflet_area_check.py")
    ap.add_argument("--out-json", default="memble_report.json")
    ap.add_argument("--out-txt", default="memble_report.txt")
    args = ap.parse_args()

    for name in args.allow:
        if name not in CHECKS:
            sys.stderr.write("verify_system: unknown check %r; known checks are %s\n"
                             % (name, ", ".join(CHECKS)))
            return 2

    resid, resn, aname, xyz, box = read_gro(args.gro)
    lipid_keys = set(s[:5] for s in args.lipids.split())
    solvent_keys = set(s[:5] for s in SOLVENT)

    is_lipid = np.array([r[:5] in lipid_keys for r in resn])
    is_solv = np.array([r[:5] in solvent_keys for r in resn])
    is_prot = ~(is_lipid | is_solv)

    rep = Report()
    skip = set(args.allow)

    # --- bilayer midplane ------------------------------------------------
    # The midplane is the mean z of every lipid bead. Taking it from head beads
    # instead needs a head bead per species, and a sterol has none: its hydroxyl
    # sits well below the phosphate plane, so a mean over the two is a height at
    # which nothing lies. The tails dominate the bead count and they are centered
    # on the midplane, so the mean over all beads needs no table.
    midplane = (float(np.mean(xyz[is_lipid, 2])) if is_lipid.any()
                else float("nan"))

    # --- 1. secondary structure ------------------------------------------
    prot_resids = sorted(set(resid[i] for i in range(len(resn)) if is_prot[i]))
    n_prot_res = len(prot_resids)
    if "secondary_structure" in skip:
        rep.add("secondary_structure", True, "skipped by --allow", skipped=True)
    elif not args.ss_string:
        rep.add("secondary_structure", False,
                "no secondary structure string was recorded for this build. "
                "memble passes the string it gave martinize2; a build without "
                "one cannot be reproduced.",
                measured=None, expected="a string of %d characters" % n_prot_res)
    else:
        ok = len(args.ss_string) == n_prot_res
        rep.add("secondary_structure", ok,
                "the string covers every protein residue" if ok else
                "the string and the protein disagree on the number of residues; "
                "martinize2 assigned bonded parameters to a different set of "
                "residues from the ones in the system.",
                measured=len(args.ss_string), expected=n_prot_res)

    # --- 2. composition ---------------------------------------------------
    if "composition" in skip:
        rep.add("composition", True, "skipped by --allow", skipped=True)
    else:
        want_up = parse_ratio(args.expect_upper)
        want_lo = parse_ratio(args.expect_lower)
        counts = {"upper": Counter(), "lower": Counter()}
        # One entry per molecule, and the leaflet from the mean z of that
        # molecule. A lipid is counted whether or not it carries a bead the
        # script knows the name of.
        by_mol = defaultdict(list)
        for i in range(len(resn)):
            if is_lipid[i]:
                by_mol[(resid[i], resn[i])].append(i)
        for (rid_, rn_), ids in by_mol.items():
            side = ("upper" if float(np.mean(xyz[ids, 2])) >= midplane
                    else "lower")
            counts[side][rn_] += 1
        measured = {s: dict(counts[s]) for s in counts}
        if not want_up and not want_lo:
            total = sum(counts["upper"].values()) + sum(counts["lower"].values())
            rep.add("composition", total > 0,
                    "no requested composition was recorded, so the placed "
                    "lipids are reported without a comparison."
                    if total else "no lipid was placed.",
                    measured=measured, expected=None)
        else:
            bad = []
            for side, want in (("upper", want_up), ("lower", want_lo)):
                if not want:
                    continue
                tot_w = sum(want.values())
                tot_c = sum(counts[side].get(k, 0) for k in want)
                if tot_c == 0:
                    bad.append("%s leaflet holds none of the requested lipids" % side)
                    continue
                for lipid, w in want.items():
                    f_want = w / tot_w
                    f_got = counts[side].get(lipid, 0) / tot_c
                    if abs(f_got - f_want) > args.composition_tol:
                        bad.append("%s %s: asked %.3f, placed %.3f"
                                   % (side, lipid, f_want, f_got))
            rep.add("composition", not bad,
                    "the placed lipids match the requested ratio"
                    if not bad else "; ".join(bad),
                    measured=measured,
                    expected={"upper": want_up, "lower": want_lo})

    # --- 3. protein placement --------------------------------------------
    if "protein_placement" in skip:
        rep.add("protein_placement", True, "skipped by --allow", skipped=True)
    elif not is_prot.any():
        rep.add("protein_placement", True, "the system holds no protein",
                skipped=True)
    else:
        want_res = parse_resids(args.tm_resids)
        if want_res:
            sel = np.array([is_prot[i] and resid[i] in want_res
                            for i in range(len(resn))])
            label = "transmembrane range %s" % args.tm_resids
        else:
            lo, hi = midplane - 1.5, midplane + 1.5
            sel = is_prot & (xyz[:, 2] > lo) & (xyz[:, 2] < hi)
            label = "protein beads within 1.5 nm of the midplane"
        if not sel.any():
            rep.add("protein_placement", False,
                    "no protein bead lies in the bilayer. The protein is "
                    "outside the membrane, and a run of this system measures "
                    "a protein in water.",
                    measured=None, expected="a transmembrane segment")
        else:
            zc = float(np.mean(xyz[sel, 2]))
            d = abs(zc - midplane)
            rep.add("protein_placement", d <= args.placement_tol,
                    "the %s is centered on the bilayer midplane" % label
                    if d <= args.placement_tol else
                    "the %s sits %.2f nm off the bilayer midplane" % (label, d),
                    measured=round(d, 3), expected="<= %.2f nm" % args.placement_tol)

    # --- 4. charge ---------------------------------------------------------
    if "charge" in skip:
        rep.add("charge", True, "skipped by --allow", skipped=True)
    else:
        dirs = list(args.itp_dir) if args.itp_dir else ["."]
        dirs.append(os.path.dirname(os.path.abspath(args.top)))
        dirs.extend(top_include_dirs(args.top))
        seen_d, uniq = set(), []
        for d in dirs:
            r = os.path.realpath(d)
            if r not in seen_d:
                seen_d.add(r)
                uniq.append(d)
        dirs = uniq
        qmap, _ = itp_charges(dirs)
        mols = read_top_molecules(args.top)
        total, unknown = 0.0, []
        for name, count in mols:
            if name in qmap:
                total += qmap[name] * count
            else:
                unknown.append(name)
        if unknown:
            rep.add("charge", False,
                    "no topology was found for %s, so the net charge could not "
                    "be summed. The directories read were: %s"
                    % (", ".join(sorted(set(unknown))), ", ".join(dirs)),
                    measured=None, expected=0,
                    remedy="The topology of every moleculetype has to be readable.\n"
                           "  1. point at the lipidome as well as the build directory:\n"
                           "       verify_system.py --itp-dir . --itp-dir $M3_DIR ...\n"
                           "  2. or check that the #include lines of system.top resolve:\n"
                           "       grep '#include' system.top")
        else:
            ok = abs(total) <= args.charge_tol
            rep.add("charge", ok,
                    "the system is neutral" if ok else
                    "the system carries a net charge of %+.3f e" % total,
                    measured=round(total, 4), expected=0)

    # --- 5. periodic image in z -------------------------------------------
    if "periodic_image" in skip:
        rep.add("periodic_image", True, "skipped by --allow", skipped=True)
    elif not is_prot.any():
        rep.add("periodic_image", True, "the system holds no protein", skipped=True)
    else:
        span = float(xyz[is_prot, 2].max() - xyz[is_prot, 2].min())
        clearance = float(box[2]) - span
        need = 2.0 * args.cutoff
        rep.add("periodic_image", clearance >= need,
                "the protein and the image of the protein in z are %.2f nm apart"
                % clearance if clearance >= need else
                "the protein spans %.2f nm in a box of %.2f nm, leaving %.2f nm "
                "to its periodic image. A protein that sees its own image gives "
                "an ordinary-looking trajectory of the wrong system."
                % (span, box[2], clearance),
                measured=round(clearance, 3), expected=">= %.2f nm" % need)

    # --- 6. water layer ----------------------------------------------------
    if "water_layer" in skip or args.water_nm <= 0:
        rep.add("water_layer", True,
                "skipped by --allow" if "water_layer" in skip
                else "no water thickness was requested", skipped=True)
    else:
        wsel = np.array([resn[i][:5] in ("W", "WF") for i in range(len(resn))])
        if not wsel.any():
            rep.add("water_layer", False, "the system holds no water bead.",
                    measured=0, expected=">= %.2f nm per side" % args.water_nm)
        else:
            memb_hi = float(xyz[is_lipid, 2].max()) if is_lipid.any() else midplane
            memb_lo = float(xyz[is_lipid, 2].min()) if is_lipid.any() else midplane
            up = float(xyz[wsel, 2].max()) - memb_hi
            lo = memb_lo - float(xyz[wsel, 2].min())
            worst = min(up, lo)
            ok = worst >= 0.8 * args.water_nm
            rep.add("water_layer", ok,
                    "water reaches %.2f nm above and %.2f nm below the bilayer"
                    % (up, lo),
                    measured={"above": round(up, 3), "below": round(lo, 3)},
                    expected=">= %.2f nm per side" % args.water_nm)

    # --- 7. overlap --------------------------------------------------------
    if "overlap" in skip:
        rep.add("overlap", True, "skipped by --allow", skipped=True)
    else:
        worst, pair = grid_min_distance(resid, resn, xyz, box, args.min_dist)
        ok = worst >= args.min_dist
        rep.add("overlap", ok,
                "the closest pair of beads from different molecules is %.3f nm "
                "apart" % worst if ok else
                "beads of %s and %s are %.3f nm apart. Minimization of this "
                "system reports an infinite force, or moves the two molecules "
                "far enough to distort them." % (pair[0], pair[1], worst),
                measured=round(worst, 4), expected=">= %.3f nm" % args.min_dist)

    # --- 8. leaflet areas --------------------------------------------------
    # The areas are measured by leaflet_area_check.py during the build. Reading
    # the result back here keeps one report: a build cannot say PASS while a
    # leaflet was reported as mismatched somewhere earlier in the same log.
    if "leaflet_area" in skip:
        rep.add("leaflet_area", True, "skipped by --allow", skipped=True)
    elif not args.leaflet_json or not os.path.isfile(args.leaflet_json):
        rep.add("leaflet_area", True,
                "no leaflet measurement was recorded for this build", skipped=True)
    else:
        try:
            with open(args.leaflet_json) as fh:
                la = json.load(fh)
        except (OSError, ValueError):
            la = {}
        res = la.get("result", "")
        shared = la.get("shared_species", [])
        worst = max(shared, key=lambda x: x.get("relative_difference", 0.0),
                    default=None)
        if res == "FAIL" and worst:
            rep.add("leaflet_area", False,
                    "%s occupies %.1f%% more area in one leaflet than in the "
                    "other, against a standard error of %.1f%% on that difference"
                    % (worst["lipid"], 100 * worst["relative_difference"],
                       100 * worst.get("relative_standard_error", 0.0)),
                    measured=round(100 * worst["relative_difference"], 2),
                    expected="within the tolerance, or within twice the "
                             "standard error")
        elif res == "REPORTED":
            f = la.get("shared_fraction")
            if worst and f is not None:
                d = "%s differs by %.1f%% between the leaflets, which hold %.0f%% " \
                    "of their lipids in common, so the difference does not report " \
                    "the tension between them" % (worst["lipid"],
                        100 * worst["relative_difference"], 100 * f)
            else:
                d = "no lipid is present in both leaflets, so the two leaflets " \
                    "were not compared"
            rep.add("leaflet_area", True, d, skipped=True)
        else:
            d = (100 * worst["relative_difference"]) if worst else 0.0
            rep.add("leaflet_area", True,
                    "the leaflets are matched; the largest difference is %.1f%%"
                    % d, measured=round(d, 2), expected="within the tolerance")

    # ----------------------------------------------------------------- out
    meta = {}
    if args.meta and os.path.isfile(args.meta):
        try:
            with open(args.meta) as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            meta = {}
    meta.setdefault("ss_mode", args.ss_mode)
    meta.setdefault("ss_string", args.ss_string)
    meta.setdefault("tm_resids", args.tm_resids)
    meta["box_nm"] = [round(float(v), 4) for v in box]
    meta["n_beads"] = int(len(resn))
    meta["n_protein_residues"] = int(n_prot_res)

    failures = rep.failed()
    result = {"result": "PASS" if not failures else "FAIL",
              "build": meta, "checks": rep.rows}

    with open(args.out_json, "w") as fh:
        json.dump(result, fh, indent=2)

    lines = ["memble verification report", "=" * 60, ""]
    for key in ("ss_mode", "tm_resids", "box_nm", "n_beads",
                "n_protein_residues"):
        if key in meta:
            lines.append("%-22s %s" % (key, meta[key]))
    if meta.get("ss_string"):
        lines.append("%-22s %s" % ("ss_string", meta["ss_string"]))
    lines += ["", "%-22s %-6s %s" % ("check", "result", "detail"), "-" * 60]
    for r in rep.rows:
        mark = "SKIP" if r["skipped"] else ("PASS" if r["pass"] else "FAIL")
        lines.append("%-22s %-6s %s" % (r["check"], mark, r["detail"]))
        if r["measured"] is not None:
            lines.append("%-22s %-6s measured %s, expected %s"
                         % ("", "", r["measured"], r["expected"]))
    if failures:
        lines += ["", "=" * 60, "WHAT TO DO", "=" * 60]
        for r in failures:
            lines.append("")
            lines.append("%s:" % r["check"])
            for ln in (r["remedy"] or "no remedy recorded").split("\n"):
                lines.append("  " + ln)
    lines += ["", "-" * 60, "RESULT: %s" % result["result"], ""]
    text = "\n".join(lines)
    with open(args.out_txt, "w") as fh:
        fh.write(text)
    sys.stdout.write(text)

    return 0 if not failures else 1


def grid_min_distance(resid, resn, xyz, box, floor):
    """Smallest inter-molecular bead distance under the minimum image convention.

    Cells of the side of the floor distance keep this linear in the number of
    beads, so a system of a million beads is checked in seconds.
    """
    mol = np.array(["%d:%s" % (resid[i], resn[i]) for i in range(len(resn))])
    cell = max(floor, 0.05) * 2.0
    dims = np.maximum((box / cell).astype(int), 1)
    idx = np.floor(xyz / cell).astype(int) % dims
    buckets = defaultdict(list)
    for i in range(len(xyz)):
        buckets[(idx[i, 0], idx[i, 1], idx[i, 2])].append(i)

    best, pair = float("inf"), ("", "")
    offsets = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    for key, here in buckets.items():
        near = []
        for off in offsets:
            k = ((key[0] + off[0]) % dims[0], (key[1] + off[1]) % dims[1],
                 (key[2] + off[2]) % dims[2])
            near.extend(buckets.get(k, ()))
        if not near:
            continue
        near = np.array(near)
        for i in here:
            d = xyz[near] - xyz[i]
            d -= box * np.round(d / box)
            r = np.sqrt((d * d).sum(axis=1))
            same = mol[near] == mol[i]
            r[same] = np.inf
            j = int(np.argmin(r))
            if r[j] < best:
                best, pair = float(r[j]), (mol[i], mol[near[j]])
    if best == float("inf"):
        return float("inf"), ("", "")
    return best, pair


if __name__ == "__main__":
    sys.exit(main())
