#!/usr/bin/env python3
"""
gen_ss.py

Generate a martinize2 secondary-structure string (-ss) from transmembrane
ranges. Residues inside a TM range are marked H (helix); everything else is
marked C (coil). The string is built by reading the exact PDB that martinize2
will consume (the oriented all-atom PDB), so it always matches the residue count
and chain order martinize2 sees.

This is for single-pass TM and TM-JM constructs, where the membrane-spanning
segment is helical and the flanking and juxtamembrane regions should stay
flexible. It is NOT appropriate for multi-helix proteins such as GPCRs; for
those, let DSSP assign the secondary structure (SS_MODE=dssp).

TM specification (--tm):
  per-chain : "A:65-88;B:65-88"      apply 65-88 as helix in chains A and B
  all-chain : "ALL:65-88"            apply 65-88 as helix in every chain
  multiple  : "A:34-60,72-96"        more than one helix range in a chain

The residue numbers are the numbers present in the PDB being read.
"""

import argparse
import sys
from collections import OrderedDict


def parse_tm(spec):
    """Return dict chain -> list of (start,end). Chain 'ALL' applies to every chain."""
    out = {}
    for part in spec.split(";"):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise SystemExit(f"bad --tm token (need CHAIN:start-end): {part}")
        chain, ranges = part.split(":", 1)
        spans = []
        for r in ranges.split(","):
            r = r.strip()
            if "-" not in r:
                raise SystemExit(f"bad range (need start-end): {r}")
            a, b = r.split("-", 1)
            spans.append((int(a), int(b)))
        out.setdefault(chain.strip(), []).extend(spans)
    return out


def read_chain_residues(pdb):
    """Return OrderedDict chain -> ordered list of residue numbers (unique, in
    file order). Uses CA atoms to define residues."""
    chains = OrderedDict()
    seen = set()
    for line in open(pdb):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        chain = line[21:22].strip() or "A"
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        key = (chain, resnum)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(chain, []).append(resnum)
    return chains


def in_any(resnum, spans):
    return any(a <= resnum <= b for a, b in spans)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="oriented all-atom PDB martinize2 will read")
    ap.add_argument("--tm", required=True, help='TM ranges, e.g. "A:65-88;B:65-88" or "ALL:65-88"')
    args = ap.parse_args()

    tm = parse_tm(args.tm)
    chains = read_chain_residues(args.pdb)
    if not chains:
        sys.exit("gen_ss: no CA atoms found in " + args.pdb)

    ss = []
    for chain, residues in chains.items():
        spans = tm.get(chain) or tm.get("ALL") or []
        if not spans:
            sys.stderr.write(f"gen_ss: warning: no TM range for chain {chain}; marking all coil\n")
        for r in residues:
            ss.append("H" if in_any(r, spans) else "C")

    sys.stdout.write("".join(ss))


if __name__ == "__main__":
    main()
