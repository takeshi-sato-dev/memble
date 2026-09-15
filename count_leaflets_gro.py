#!/usr/bin/env python3
"""How many phospholipid molecules each leaflet holds, read from a .gro file.

The rule is the rule of Section 2.6: a molecule belongs to the leaflet its own
head bead sits toward, taken against the tail beads of that same molecule. No
plane is drawn through the bilayer, because the protein bends the membrane and a
plane then misassigns the lipids beside it.

The script reads the .gro format directly and imports nothing outside the
standard library, so that it runs on a machine that carries neither MDAnalysis
nor mdtraj. A .gro line is fixed width: five columns of residue number, five of
residue name, five of atom name, five of atom number, then three fields of eight
for the coordinates in nanometres. A new residue is taken to start wherever the
residue number or the residue name changes, which is correct because the
molecules of a .gro file are contiguous.

What it is for. Section 3.5 takes the phospholipid imbalance *f* from the
structure that begins the production run and not from the build, because the
equilibration removes part of a large imbalance. Run on step6.0 and on step6.6 of
the same system, this script gives both numbers.

Usage:
    python3 count_leaflets_gro.py step6.0.gro step6.6.gro
    python3 count_leaflets_gro.py --walk '~/ext_*/curve/*/*_work'
    python3 count_leaflets_gro.py --walk '~/counts_run/*/*/*_work'

With --walk the argument is a shell pattern of directories; each directory that
holds both a step6.0 and a step6.6 file is reported as one row.
"""
import argparse
import glob
import os
import re
import sys

TAIL = re.compile(r"^[CD]\d[AB]$")
HEADS = ("PO4", "P1", "P2", "P3", "CNO")
PL = ("DLPC", "DLiPC", "DIPC", "PSM", "DPSM", "DOPS", "POP2", "PIP2", "POPI")
STEROL = "CHOL"
STEROL_HEAD = "ROH"


def residues(path):
    """Yield (resname, head z, list of tail z) for every residue of the file."""
    with open(path) as fh:
        lines = fh.readlines()
    cur = None
    out = []
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        key = (line[0:5], line[5:10].strip())
        name = line[10:15].strip()
        try:
            z = float(line[36:44])
        except ValueError:
            continue
        if cur is None or cur[0] != key:
            cur = [key, [], []]
            out.append(cur)
        if name in HEADS or name == STEROL_HEAD:
            if not cur[1]:
                cur[1].append(z)
        else:
            cur[2].append(z)
    return [(k[1], h[0], t) for k, h, t in out if h and t]


def split(path, species):
    up = lo = 0
    per = {}
    for name, hz, tz in residues(path):
        if name not in species:
            continue
        side = 0 if hz - sum(tz) / len(tz) > 0 else 1
        per.setdefault(name, [0, 0])[side] += 1
        if side == 0:
            up += 1
        else:
            lo += 1
    return up, lo, per


def report(a, b, label, species, detail):
    u0, l0, p0 = split(a, species)
    u6, l6, p6 = split(b, species)
    n0, n6 = u0 + l0, u6 + l6
    if not n0 or not n6:
        print("%-40s no molecule of %s matched" % (label, "/".join(species)))
        return
    print("%-40s %4d /%4d  %4d /%4d  %+8.2f %+8.2f  %+d"
          % (label, u0, l0, u6, l6,
             100.0 * (u0 - l0) / n0, 100.0 * (u6 - l6) / n6, u6 - u0))
    if detail:
        for k in sorted(set(p0) | set(p6)):
            x = p0.get(k, [0, 0])
            y = p6.get(k, [0, 0])
            print("    %-8s %3d /%3d   %3d /%3d" % (k, x[0], x[1], y[0], y[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="step6.0.gro and step6.6.gro")
    ap.add_argument("--walk", default=None,
                    help="shell pattern of work directories to scan")
    ap.add_argument("--sterol", action="store_true",
                    help="count cholesterol instead of the phospholipids")
    ap.add_argument("--detail", action="store_true",
                    help="also print the split of each species")
    a = ap.parse_args()

    species = (STEROL,) if a.sterol else PL
    print("%-40s %-12s %-12s %-9s %-9s %s"
          % ("system", "step6.0", "step6.6", "f build", "f run", "moved"))

    if a.walk:
        for w in sorted(glob.glob(os.path.expanduser(a.walk))):
            first = sorted(glob.glob(os.path.join(w, "step6.0*.gro")))
            last = sorted(glob.glob(os.path.join(w, "step6.6*.gro")))
            if not first or not last:
                continue
            report(first[0], last[0],
                   os.path.relpath(w, os.path.expanduser("~")), species,
                   a.detail)
        return

    if len(a.files) != 2:
        print("give two .gro files, or use --walk", file=sys.stderr)
        raise SystemExit(2)
    report(a.files[0], a.files[1], os.path.basename(os.path.dirname(a.files[0])),
           species, a.detail)


if __name__ == "__main__":
    main()
