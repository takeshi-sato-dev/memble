#!/usr/bin/env python3
"""
ss_from_dssp.py

Emit a one-character-per-residue secondary-structure string (H for helix, C for
everything else) for an all-atom protein PDB, in CA order. Used by the multi-TM
orientation path so a 7-TM bundle (GPCR) can be oriented on its helix axes
rather than on a single principal component.

Helix detection, in order of preference:
  1. --dssp /path/to/mkdssp   run DSSP 2.2.1 / 3.x, parse its SS column
  2. otherwise                fall back to mdtraj.compute_dssp

Output: the SS string on stdout, nothing else. Letters G (3-10) and I (pi) are
folded into H, since for orientation any membrane-spanning helix counts.
"""

import argparse
import subprocess
import sys
import tempfile
import os


def ca_order_resids(pdb):
    """List of (chain, resSeq) in CA order, to align SS back to residues."""
    out = []
    seen = set()
    with open(pdb) as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
                key = (ln[21], ln[22:27])
                if key not in seen:
                    seen.add(key)
                    out.append(key)
    return out


def ss_via_dssp(pdb, dssp_bin):
    """Run mkdssp, return {(chain,resseq_field): 'H'/'C'}."""
    with tempfile.NamedTemporaryFile(suffix=".dssp", delete=False) as tf:
        outp = tf.name
    try:
        # DSSP 3.x: mkdssp -i in.pdb -o out.dssp ; some builds accept positional
        for args in ([dssp_bin, "-i", pdb, "-o", outp],
                     [dssp_bin, pdb, outp]):
            try:
                r = subprocess.run(args, capture_output=True, text=True)
                if r.returncode == 0 and os.path.getsize(outp) > 0:
                    break
            except Exception:
                continue
        else:
            return None
        mapping = {}
        started = False
        with open(outp) as fh:
            for ln in fh:
                if ln.startswith("  #  RESIDUE"):
                    started = True
                    continue
                if not started:
                    continue
                if len(ln) < 17 or ln[13] == "!":
                    continue
                ch = ln[11]
                resnum = ln[5:10].strip()
                sscode = ln[16]
                mapping[(ch, resnum)] = "H" if sscode in "HGI" else "C"
        return mapping or None
    finally:
        try:
            os.unlink(outp)
        except OSError:
            pass


def ss_via_mdtraj(pdb):
    try:
        import mdtraj as md
    except Exception:
        return None
    try:
        t = md.load(pdb)
        ss = md.compute_dssp(t, simplified=False)[0]
        # one entry per residue in topology order
        out = []
        for s in ss:
            out.append("H" if s in ("H", "G", "I") else "C")
        return out
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--dssp", default=None, help="path to mkdssp (3.x); omit to use mdtraj")
    args = ap.parse_args()

    order = ca_order_resids(args.pdb)
    if not order:
        sys.exit("ERROR: no CA atoms in %s" % args.pdb)

    ss = None
    if args.dssp and args.dssp != "mdtraj":
        mp = ss_via_dssp(args.pdb, args.dssp)
        if mp:
            # align to CA order; DSSP resnum field is resSeq without icode space
            ss = "".join(mp.get((ch, rs.strip().lstrip()), mp.get((ch, rs.strip()), "C"))
                         for (ch, rs) in order)
            # the key formats differ slightly; retry on plain integer match
            if ss.count("H") == 0:
                # rebuild mapping keyed by (chain, int)
                mp2 = {}
                for (ch, rs), v in mp.items():
                    try:
                        mp2[(ch, int(rs))] = v
                    except ValueError:
                        pass
                ss = "".join(mp2.get((ch, int(rs)), "C") for (ch, rs) in order)

    if ss is None or ss.count("H") == 0:
        seq = ss_via_mdtraj(args.pdb)
        if seq is not None:
            if len(seq) == len(order):
                ss = "".join(seq)
            else:
                ss = "".join(seq)[:len(order)].ljust(len(order), "C")

    if not ss or ss.count("H") == 0:
        sys.exit("ERROR: could not derive secondary structure (no helices found)")

    sys.stdout.write(ss)


if __name__ == "__main__":
    main()
