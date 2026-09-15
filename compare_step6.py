#!/usr/bin/env python3
"""Whether a molecule that changes leaflet during the equilibration has changed
leaflet, or is misread by the leaflet rule once the bilayer is strained.

Section 3.5 takes the phospholipid imbalance from the structure that begins the
production run, because the equilibration removes part of a large imbalance.
That step rests on the leaflet rule of Section 2.6 holding in a strained
bilayer, and the rule was checked on a bilayer that is not strained. This
script makes the check on any system, by comparing the packed structure, which
carries the numbers the build assigned, with the structure that ends the
equilibration.

The two files hold the same molecules in the same order, so every molecule can
be followed from one to the other without a trajectory. For each molecule that
changes leaflet the script reports

  dz          the head bead less the mean of the tail beads, along the normal,
              in the second structure, against the median of that species. A
              molecule lying in the plane of the bilayer has a dz near zero and
              its leaflet is not defined
  distance    the lateral distance from the nearest protein bead, against the
              mean of that species

A real transfer gives a dz of the usual size and no preference for the
neighbourhood of the protein. A failure of the rule gives a small dz, or puts
the molecules that changed beside the protein, where the bilayer is bent.

Usage:
    python3 compare_step6.py step6.0.gro step6.6.gro
    python3 compare_step6.py step6.0.gro step6.6.gro --near-nm 2.0
"""
import argparse
import math
import re

TAIL = re.compile(r"^[CD]\d[AB]$")
HEADS = ("PO4", "P1", "P2", "P3", "CNO")
PL = ("DLPC", "DLiPC", "DIPC", "PSM", "DPSM", "DOPS", "POP2", "POP2_",
      "PIP2", "POPI")
PROT = "BB"


def is_pl(name):
    return any(name == x or name.startswith(x) for x in PL)


def read(path):
    """Return the molecules of a .gro file and the protein bead positions."""
    with open(path) as fh:
        lines = fh.readlines()
    mols = []
    prot = []
    cur = None
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        key = (line[0:5], line[5:10].strip())
        name = line[10:15].strip()
        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue
        if name == PROT:
            prot.append((x, y))
            continue
        if not is_pl(key[1]):
            continue
        if cur is None or cur[0] != key:
            cur = [key, None, [], None]
            mols.append(cur)
        if name in HEADS and cur[1] is None:
            cur[1] = (x, y, z)
        elif TAIL.match(name):
            cur[2].append(z)
    out = []
    for key, head, tails, _ in mols:
        if head is None or not tails:
            continue
        out.append((key[1], head[0], head[1], head[2] - sum(tails) / len(tails)))
    return out, prot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("packed", help="step6.0 .gro, before any dynamics")
    ap.add_argument("equilibrated", help="step6.6 .gro, the start of production")
    ap.add_argument("--near-nm", type=float, default=2.0)
    a = ap.parse_args()

    first, _ = read(a.packed)
    second, prot = read(a.equilibrated)
    if len(first) != len(second):
        raise SystemExit("the two files hold %d and %d molecules"
                         % (len(first), len(second)))
    print("%d phospholipid molecules, %d protein beads" % (len(first), len(prot)))

    def dist(x, y):
        if not prot:
            return float("nan")
        return math.sqrt(min((x - px) ** 2 + (y - py) ** 2 for px, py in prot))

    species = {}
    for i, (name, x, y, dz) in enumerate(second):
        species.setdefault(name, []).append(i)

    for name in sorted(species):
        idx = species[name]
        med = sorted(abs(second[i][3]) for i in idx)[len(idx) // 2]
        dmean = sum(dist(second[i][1], second[i][2]) for i in idx) / len(idx)
        moved = [i for i in idx
                 if (first[i][3] > 0) != (second[i][3] > 0)]
        up0 = sum(1 for i in idx if first[i][3] > 0)
        up6 = sum(1 for i in idx if second[i][3] > 0)
        print("")
        print("%-6s %d molecules | packed %d/%d | equilibrated %d/%d | %d changed side"
              % (name, len(idx), up0, len(idx) - up0, up6, len(idx) - up6,
                 len(moved)))
        print("       all molecules: median |dz| %.2f nm, mean distance to the "
              "protein %.2f nm" % (med, dmean))
        if not moved:
            print("       every molecule stays where the build put it")
            continue
        print("       %-6s %-9s %-9s %s" % ("mol", "dz nm", "|dz|/median", "distance nm"))
        flat = near = 0
        for i in moved[:25]:
            dz = second[i][3]
            d = dist(second[i][1], second[i][2])
            if abs(dz) < 0.6:
                flat += 1
            if d < a.near_nm:
                near += 1
            print("       %-6d %+9.2f %-9.2f %.2f" % (i, dz, abs(dz) / med, d))
        if len(moved) > 25:
            print("       ... %d more" % (len(moved) - 25))
            for i in moved[25:]:
                if abs(second[i][3]) < 0.6:
                    flat += 1
                if dist(second[i][1], second[i][2]) < a.near_nm:
                    near += 1
        print("       of the %d that changed, %d lie flatter than 0.6 nm and %d "
              "sit within %.1f nm of the protein"
              % (len(moved), flat, near, a.near_nm))
        if flat == 0 and near <= len(moved) / 2:
            print("       ** these molecules stand as the others do and are not "
                  "gathered at the protein: the transfer is real")
        else:
            print("       ** these molecules lie flat or sit at the protein: "
                  "the leaflet rule is failing here, not the lipid moving")


if __name__ == "__main__":
    main()
