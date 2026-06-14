#!/usr/bin/env python3
"""
fix_protein_resid.py

Restore the original residue numbers (and chain identity) of the protein in a
built system. martinize2 renumbers each chain from 1, so after assembly all
protein chains share resids 1..N and overlap; CHARMM-GUI instead keeps the input
numbering (e.g. 54..103 per chain). This reads the TRUE per-chain residue numbers
from the oriented all-atom PDB that was fed to martinize2 (which preserves the
input resSeq) and writes them back onto the protein beads of system.gro, so the
PSF/PDB and any analysis see the correct numbering. Nothing is hardcoded: the
numbers come from the input PDB, the chain order from the topology.

Usage:
  fix_protein_resid.py --gro system.gro --oriented oriented_aa.pdb \
      --top system.top [--itp-dir DIR]
"""

import argparse
import os
import re

AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def chains_resseq_from_pdb(path):
    """Ordered list of (chain, [resSeq in order of appearance, unique])."""
    order = []
    seen = {}
    for ln in open(path):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        resn = ln[17:21].strip().upper()
        if resn not in AA:
            continue
        chain = ln[21]
        try:
            rs = int(ln[22:26])
        except ValueError:
            continue
        if chain not in seen:
            seen[chain] = []
            order.append(chain)
        if not seen[chain] or seen[chain][-1] != rs:
            # keep order; collapse consecutive duplicates (same residue, many atoms)
            if rs not in seen[chain]:
                seen[chain].append(rs)
    return [(c, seen[c]) for c in order]


def parse_itp_atoms(itp_path):
    """{molname: (natoms, [resnr per bead], [resname per bead])}."""
    out = {}
    mol = None
    resnr = []
    resn = []
    section = None

    def flush():
        if mol is not None:
            out[mol] = (len(resnr), list(resnr), list(resn))

    for raw in open(itp_path):
        s = raw.split(";")[0].strip()
        if not s:
            continue
        m = re.match(r"\[\s*([a-zA-Z_0-9]+)\s*\]", s)
        if m:
            sec = m.group(1).lower()
            if sec == "moleculetype":
                flush()
                mol, resnr, resn = None, [], []
            section = sec
            continue
        f = s.split()
        if section == "moleculetype":
            mol = f[0]
        elif section == "atoms" and len(f) >= 5:
            try:
                resnr.append(int(f[2]))
            except ValueError:
                resnr.append(len(resnr) + 1)
            resn.append(f[3])
    flush()
    return out


def collect_itps(top_path, itp_dir):
    paths = []
    base = os.path.dirname(os.path.abspath(top_path))
    for ln in open(top_path):
        m = re.match(r'\s*#include\s+"([^"]+)"', ln)
        if m:
            p = m.group(1)
            cand = p if os.path.isabs(p) else os.path.join(base, p)
            if os.path.exists(cand):
                paths.append(cand)
    if itp_dir and os.path.isdir(itp_dir):
        for fn in os.listdir(itp_dir):
            if fn.endswith(".itp"):
                paths.append(os.path.join(itp_dir, fn))
    return paths


def mol_order(top_path):
    order = []
    in_mol = False
    for ln in open(top_path):
        s = ln.split(";")[0].strip()
        if s.lower().startswith("[ molecules"):
            in_mol = True
            continue
        if in_mol and s:
            if s.startswith("["):
                break
            p = s.split()
            order.append((p[0], int(p[1])))
    return order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--oriented", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--itp-dir", default="")
    args = ap.parse_args()

    chains = chains_resseq_from_pdb(args.oriented)   # [(chain,[resseq...]),...]
    atoms_by_mol = {}
    for itp in collect_itps(args.top, args.itp_dir):
        try:
            for mol, v in parse_itp_atoms(itp).items():
                atoms_by_mol.setdefault(mol, v)
        except Exception:
            continue

    order = mol_order(args.top)

    # protein molecule instances, in order, each as (natoms, resnr-per-bead)
    prot_instances = []
    for molname, count in order:
        info = atoms_by_mol.get(molname)
        if not info:
            continue
        natoms, resnr, resn = info
        is_pro = sum(1 for r in resn if r.upper() in AA) >= max(2, 0.5 * len(resn))
        if is_pro:
            for _ in range(count):
                prot_instances.append((natoms, resnr))

    if not prot_instances:
        print("fix_protein_resid: no protein molecules found; nothing to do")
        return
    if len(prot_instances) != len(chains):
        print("fix_protein_resid: WARNING %d protein chains in topology but %d in "
              "%s; mapping the overlap in order"
              % (len(prot_instances), len(chains), os.path.basename(args.oriented)))

    # read gro
    L = open(args.gro).read().splitlines()
    title, n = L[0], int(L[1])
    body = L[2:2 + n]
    box = L[2 + n]

    # walk the leading protein beads, assign original resids per chain
    idx = 0
    nfixed = 0
    for ci, (natoms, resnr) in enumerate(prot_instances):
        if ci >= len(chains):
            idx += natoms
            continue
        orig = chains[ci][1]                  # original resSeq list for this chain
        # group beads of this molecule into residues by runs of equal resnr
        groups = []
        last = None
        for r in resnr:
            if r != last:
                groups.append(0)
                last = r
            groups[-1] += 1
        # assign: residue g -> orig[g]
        if len(groups) != len(orig):
            print("fix_protein_resid: WARNING chain #%d has %d CG residues but %d "
                  "original residues; assigning the overlap in order"
                  % (ci + 1, len(groups), len(orig)))
        b = 0
        for g, cnt in enumerate(groups):
            newid = orig[g] if g < len(orig) else (orig[-1] + 1 + g - len(orig))
            for _ in range(cnt):
                if idx + b < n:
                    ln = body[idx + b]
                    body[idx + b] = "%5d%s" % (newid % 100000, ln[5:])
                    nfixed += 1
                b += 1
        idx += natoms

    with open(args.gro, "w") as fh:
        fh.write(title + "\n%d\n" % n)
        for ln in body:
            fh.write(ln + "\n")
        fh.write(box + "\n")
    print("fix_protein_resid: restored original residue numbers on %d protein "
          "beads across %d chains" % (nfixed, len(prot_instances)))


if __name__ == "__main__":
    main()
