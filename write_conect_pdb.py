#!/usr/bin/env python3
"""
write_conect_pdb.py

Write a PDB with CONECT records (and per-molecule bonds) from a built CG system
so it displays WITH BONDS by default in VMD/PyMOL/ChimeraX, with no extra steps.

Coarse-grained .gro files carry no connectivity, and viewers guess bonds from
distance; CG bead spacing (~0.47 nm) is too large for the default guess, so the
system looks like loose dots. This reads the bonds/constraints from each
moleculetype in the topology and emits explicit CONECT records, the universal
way to ship connectivity inside a coordinate file.

Usage:
  write_conect_pdb.py --gro system.gro --top system.top --out system_view.pdb
                      [--itp-dir DIR]
"""

import argparse
import os
import re
import sys


def parse_mol_bonds(itp_path):
    """{molname: (natoms, [(i,j), ...])} using [bonds] and [constraints]
    (0-based atom indices within the molecule)."""
    out = {}
    mol = None
    natoms = 0
    bonds = []
    section = None

    def flush():
        if mol is not None:
            out[mol] = (natoms, list(bonds))

    for raw in open(itp_path):
        line = raw.split(";")[0]
        s = line.strip()
        if not s:
            continue
        m = re.match(r"\[\s*([a-zA-Z_0-9]+)\s*\]", s)
        if m:
            section = m.group(1).lower()
            if section == "moleculetype":
                flush()
                mol = None
                natoms = 0
                bonds = []
            continue
        f = s.split()
        if section == "moleculetype":
            mol = f[0]
        elif section == "atoms":
            natoms += 1
        elif section in ("bonds", "constraints"):
            try:
                i, j = int(f[0]) - 1, int(f[1]) - 1
                bonds.append((i, j))
            except (ValueError, IndexError):
                pass
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


def read_gro(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    body = L[2:2 + n]
    box = [float(v) for v in L[2 + n].split()[:3]]
    return n, body, box


AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def is_protein_moltype(info):
    """info = (natoms, bonds, anames, resns) -> protein if residues are AAs."""
    if not info or len(info) < 4 or not info[3]:
        return False
    return sum(1 for r in info[3] if r.upper() in AA) >= max(2, 0.5 * len(info[3]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--itp-dir", default="")
    args = ap.parse_args()

    molbonds = {}
    for itp in collect_itps(args.top, args.itp_dir):
        try:
            for mol, v in parse_mol_bonds(itp).items():
                molbonds.setdefault(mol, v)
        except Exception:
            continue

    n, body, box = read_gro(args.gro)

    # parse gro atoms
    resids, resns, anames, xyz = [], [], [], []
    for b in body:
        resids.append(b[0:5].strip())
        resns.append(b[5:10].strip())
        anames.append(b[10:15].strip())
        xyz.append((float(b[20:28]) * 10.0, float(b[28:36]) * 10.0,
                    float(b[36:44]) * 10.0))   # nm -> angstrom

    # build global bond list by walking molecules in topology order
    order = mol_order(args.top)
    bonds = []
    idx = 0
    for molname, count in order:
        info = molbonds.get(molname)
        if info is None:
            # advance by residue-run if unknown (single-residue molecules)
            a = guess_napm(body, idx)
            idx += a * count
            continue
        napm, mb = info
        for _ in range(count):
            for (i, j) in mb:
                if i < napm and j < napm:
                    bonds.append((idx + i, idx + j))
            idx += napm

    # write PDB
    with open(args.out, "w") as fh:
        fh.write("REMARK  CG system with explicit CONECT bonds (for viewers)\n")
        fh.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n"
                 % (box[0] * 10, box[1] * 10, box[2] * 10))
        ions = {"NA", "CL", "ION", "NA+", "CL-", "K", "MG", "CA"}
        waters = {"W", "WF", "WAT", "SOL"}
        for k in range(n):
            serial = (k + 1) % 100000
            name = anames[k][:4]
            resn = resns[k][:4]
            ru = resn.upper()
            if ru in AA:
                chain = "P"
            elif ru in waters:
                chain = "W"
            elif ru in ions:
                chain = "I"
            else:
                chain = "M"
            try:
                resid = int(resids[k]) % 10000
            except ValueError:
                resid = k % 10000
            x, y, z = xyz[k]
            # explicit PDB columns: name 13-16, resName 18-21, chainID 22,
            # resSeq 23-26, coords 31-54
            line = ("ATOM  "
                    + "%5d" % serial
                    + " "
                    + "%-4s" % name
                    + " "
                    + "%-4s" % resn
                    + "%1s" % chain
                    + "%4d" % resid
                    + "    "
                    + "%8.3f%8.3f%8.3f" % (x, y, z)
                    + "  1.00  0.00")
            fh.write(line + "\n")
        # CONECT (PDB serial is 1-based); group up to 4 partners per line
        from collections import defaultdict
        adj = defaultdict(list)
        for i, j in bonds:
            adj[i].append(j)
            adj[j].append(i)
        for i in sorted(adj):
            partners = adj[i]
            for s in range(0, len(partners), 4):
                chunk = partners[s:s + 4]
                fh.write("CONECT%5d" % ((i + 1) % 100000)
                         + "".join("%5d" % ((p + 1) % 100000) for p in chunk)
                         + "\n")
        fh.write("END\n")
    print("write_conect_pdb: %s (%d atoms, %d bonds)" % (args.out, n, len(bonds)))


def guess_napm(body, start):
    if start >= len(body):
        return 1
    rid0 = body[start][0:5]
    k = start
    while k < len(body) and body[k][0:5] == rid0:
        k += 1
    return max(k - start, 1)


if __name__ == "__main__":
    main()
