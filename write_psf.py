#!/usr/bin/env python3
"""
write_psf.py

Write an X-PLOR/CHARMM PSF with correct bonds for a built CG system, so the
system displays WITH BONDS in VMD/PyMOL exactly like CHARMM-GUI: load system.psf
as the topology, then read system.gro or a trajectory (step*.xtc) on top of it
and every frame shows connectivity.

Connectivity comes straight from each moleculetype's [bonds] and [constraints]
in the topology (ParmEd drops these for Martini because of nrexcl/funct quirks,
so we write the PSF directly). Atom names, residue names/ids and segids match
the .gro/.xtc atom order so the topology and coordinates line up.

Usage:
  write_psf.py --gro system.gro --top system.top --out system.psf [--itp-dir DIR]
"""

import argparse
import os
import re


def parse_moltypes(itp_path):
    """{molname: (natoms, [(i,j)...], [atomname...], [resname...])} from one itp,
    using [atoms] for names and [bonds]+[constraints] for connectivity (0-based)."""
    out = {}
    mol = None
    anames = []
    resns = []
    bonds = []
    section = None

    def flush():
        if mol is not None:
            out[mol] = (len(anames), list(bonds), list(anames), list(resns))

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
                anames, resns, bonds = [], [], []
            continue
        f = s.split()
        if section == "moleculetype":
            mol = f[0]
        elif section == "atoms":
            # id type resnr resname atomname cgnr [charge mass]
            anames.append(f[4] if len(f) > 4 else ("B%d" % len(anames)))
            resns.append(f[3] if len(f) > 3 else (mol or "MOL"))
        elif section in ("bonds", "constraints"):
            try:
                bonds.append((int(f[0]) - 1, int(f[1]) - 1))
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
    return n, body


AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def is_protein_moltype(info):
    """info = (natoms, bonds, anames, resns); protein if its residues are AAs."""
    if not info or len(info) < 4:
        return False
    resns = info[3]
    if not resns:
        return False
    aa_hits = sum(1 for r in resns if r.upper() in AA)
    return aa_hits >= max(2, 0.5 * len(resns))


# segid per molecule type group (protein chains get PROA.., else resname-based)
def seg_for(molname, resname):
    rn = (resname or molname)[:4].upper()
    if rn in ("W", "WF", "WAT", "SOL"):
        return "SOLV"
    if rn in ("NA", "CL", "ION", "NA+", "CL-", "K", "MG", "CA"):
        return "IONS"
    return "MEMB"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--itp-dir", default="")
    args = ap.parse_args()

    moltypes = {}
    for itp in collect_itps(args.top, args.itp_dir):
        try:
            for mol, v in parse_moltypes(itp).items():
                moltypes.setdefault(mol, v)
        except Exception:
            continue

    n, body = read_gro(args.gro)

    # gro atom fields (authoritative for names/resids/order)
    g_resid = [b[0:5].strip() for b in body]
    g_resn = [b[5:10].strip() for b in body]
    g_aname = [b[10:15].strip() for b in body]

    order = mol_order(args.top)
    bonds = []
    seg_of = [None] * n
    idx = 0
    pro_chain = 0
    for molname, count in order:
        info = moltypes.get(molname)
        napm = info[0] if info else guess_napm(body, idx)
        mb = info[1] if info else []
        is_pro = is_protein_moltype(info)
        for _ in range(count):
            if napm <= 0:
                break
            if is_pro:
                seg = "PRO%s" % chr(ord("A") + min(pro_chain, 25))
                pro_chain += 1
            else:
                seg = seg_for(molname, g_resn[idx] if idx < n else molname)
            for k in range(napm):
                if idx + k < n:
                    seg_of[idx + k] = seg
            for (i, j) in mb:
                if i < napm and j < napm:
                    bonds.append((idx + i, idx + j))
            idx += napm
    for k in range(n):
        if seg_of[k] is None:
            seg_of[k] = seg_for("", g_resn[k])

    # write PSF (X-PLOR style: named atom types)
    with open(args.out, "w") as fh:
        fh.write("PSF EXT XPLOR\n\n")
        fh.write("%10d !NTITLE\n" % 1)
        fh.write("* CG system PSF with bonds (written by memble write_psf)\n\n")
        fh.write("%10d !NATOM\n" % n)
        for k in range(n):
            try:
                rid = int(g_resid[k])
            except ValueError:
                rid = k + 1
            seg = seg_of[k]
            name = g_aname[k]
            resn = g_resn[k]
            # PSF columns: id seg resid resname name type charge mass
            fh.write("%10d %-8s %-8d %-8s %-8s %-6s %10.6f %13.4f %11d\n"
                     % (k + 1, seg[:8], rid, resn[:8], name[:8], name[:6],
                        0.0, 72.0, 0))
        fh.write("\n%10d !NBOND: bonds\n" % len(bonds))
        col = 0
        for (i, j) in bonds:
            fh.write("%10d%10d" % (i + 1, j + 1))
            col += 1
            if col % 4 == 0:
                fh.write("\n")
        if col % 4 != 0:
            fh.write("\n")
        # empty angle/dihedral/etc sections so VMD parses cleanly
        for tag, nm in [("NTHETA", "angles"), ("NPHI", "dihedrals"),
                        ("NIMPHI", "impropers"), ("NDON", "donors"),
                        ("NACC", "acceptors"), ("NNB", "")]:
            fh.write("\n%10d !%s: %s\n" % (0, tag, nm))
        fh.write("\n")
    print("write_psf: %s (%d atoms, %d bonds)" % (args.out, n, len(bonds)))


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
