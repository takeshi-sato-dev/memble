#!/usr/bin/env python3
"""
make_protein_whole.py

Make the protein contiguous across periodic boundaries and center the system on
it, so no chain is split across the z edge. A membrane-spanning TM-JM protein is
tall in z; if the assembled box places it across the z boundary, per-atom PBC
wrapping (done upstream by COBY/solvation) maps part of a chain to the opposite
face. Consecutive backbone beads then sit a full box apart, the bond constraint
(~0.31 nm) cannot be satisfied, LINCS blows up and minimisation reports an
infinite force.

This walks each protein chain in atom order and removes box jumps between bonded
neighbours (minimum-image unwrap), recenters the whole system so the protein sits
in the middle of the box, then rigidly wraps only the non-protein molecules
(lipids, water, ions) back into the box. The protein is left whole.

Nothing is hardcoded: protein chains and their sizes come from the topology.

Usage:
  make_protein_whole.py --gro system.gro --top system.top [--itp-dir DIR]
"""

import argparse
import os
import re
import numpy as np

AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def parse_itp_atoms(itp_path):
    out = {}
    mol = None
    resn = []
    section = None

    def flush():
        if mol is not None:
            out[mol] = (len(resn), list(resn))

    for raw in open(itp_path):
        s = raw.split(";")[0].strip()
        if not s:
            continue
        m = re.match(r"\[\s*([a-zA-Z_0-9]+)\s*\]", s)
        if m:
            sec = m.group(1).lower()
            if sec == "moleculetype":
                flush()
                mol, resn = None, []
            section = sec
            continue
        f = s.split()
        if section == "moleculetype":
            mol = f[0]
        elif section == "atoms" and len(f) >= 5:
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
    ap.add_argument("--top", required=True)
    ap.add_argument("--itp-dir", default="")
    args = ap.parse_args()

    atoms_by_mol = {}
    for itp in collect_itps(args.top, args.itp_dir):
        try:
            for mol, v in parse_itp_atoms(itp).items():
                atoms_by_mol.setdefault(mol, v)
        except Exception:
            continue
    order = mol_order(args.top)

    L = open(args.gro).read().splitlines()
    title, n = L[0], int(L[1])
    body = L[2:2 + n]
    box = np.array([float(v) for v in L[2 + n].split()[:3]])

    # coordinates
    X = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                  for b in body])

    # protein chains as (start, natoms); non-protein molecules as (start,natoms)
    prot_spans = []
    nonprot_spans = []
    idx = 0
    for molname, count in order:
        info = atoms_by_mol.get(molname)
        if info is None:
            # unknown: treat as single-residue molecules of guessed size
            a = guess_napm(body, idx)
            for _ in range(count):
                nonprot_spans.append((idx, a))
                idx += a
            continue
        natoms, resn = info
        is_pro = sum(1 for r in resn if r.upper() in AA) >= max(2, 0.5 * len(resn))
        for _ in range(count):
            if is_pro:
                prot_spans.append((idx, natoms))
            else:
                nonprot_spans.append((idx, natoms))
            idx += natoms

    if not prot_spans:
        print("make_protein_whole: no protein found; nothing to do")
        return

    # 1) unwrap the protein as ONE assembly, preserving the inter-chain geometry
    # of the input PDB. First make each chain internally contiguous (remove box
    # jumps between consecutive beads), then bring every chain into the same
    # periodic image as the first chain (minimum image of the chain centroid, in
    # all axes). Unwrapping chains independently would leave some chains a whole
    # box away from the others (the assembly was wrapped per-atom upstream), which
    # splits a 4-chain bundle into separate clusters that drift apart and leave
    # the membrane two at a time.
    for (st, na) in prot_spans:
        for k in range(st + 1, st + na):
            d = X[k] - X[k - 1]
            X[k] -= box * np.round(d / box)
    if prot_spans:
        st0, na0 = prot_spans[0]
        ref = X[st0:st0 + na0].mean(axis=0)
        for (st, na) in prot_spans[1:]:
            c = X[st:st + na].mean(axis=0)
            X[st:st + na] -= box * np.round((c - ref) / box)

    # 2) recenter in z. The protein is the tallest object, so center on its
    # z-midpoint: it then has equal water on both sides, which matters for an
    # asymmetric transmembrane-juxtamembrane protein whose long tail is on one
    # side only. The membrane sits inside the protein span, so it stays well
    # away from the box edge and is not displayed split. With no protein, fall
    # back to the membrane circular mean. x and y are left alone.
    lipid_idx = []
    for (st, na) in nonprot_spans:
        rn = body[st][5:10].strip().upper()
        if rn in ("W", "WF", "WAT", "SOL", "NA", "CL", "ION", "NA+", "CL-",
                  "K", "MG", "CA"):
            continue
        lipid_idx.extend(range(st, st + na))
    pro_idx = np.concatenate([np.arange(st, st + na) for (st, na) in prot_spans])
    if len(pro_idx):
        center_z = 0.5 * (X[pro_idx, 2].min() + X[pro_idx, 2].max())
    elif lipid_idx:
        ang = 2.0 * np.pi * X[lipid_idx, 2] / box[2]
        mean_ang = np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())
        center_z = (mean_ang % (2.0 * np.pi)) / (2.0 * np.pi) * box[2]
    else:
        center_z = X[:, 2].mean()
    X[:, 2] += box[2] / 2.0 - center_z

    # 3) wrap only non-protein molecules rigidly back into [0, box)
    for (st, na) in nonprot_spans:
        c = X[st:st + na].mean(axis=0)
        s = np.floor(c / box) * box
        if np.any(s != 0.0):
            X[st:st + na] -= s
    # keep the protein whole and inside the box: shift each chain rigidly by a
    # whole box vector if its centroid drifted out (does not split the chain)
    for (st, na) in prot_spans:
        c = X[st:st + na].mean(axis=0)
        s = np.floor(c / box) * box
        if np.any(s != 0.0):
            X[st:st + na] -= s

    # write back (preserve the gro fixed-width coordinate columns)
    with open(args.gro, "w") as fh:
        fh.write(title + "\n%d\n" % n)
        for i, b in enumerate(body):
            fh.write("%s%8.3f%8.3f%8.3f%s\n"
                     % (b[:20], X[i, 0], X[i, 1], X[i, 2], b[44:].rstrip()))
        fh.write("%10.5f%10.5f%10.5f\n" % (box[0], box[1], box[2]))

    span = X[pro_idx, 2].max() - X[pro_idx, 2].min()
    cushion = (box[2] - span) / 2.0
    print("make_protein_whole: unwrapped %d chains, protein z-span %.2f nm, "
          "centered in box_z %.2f nm" % (len(prot_spans), span, box[2]))
    if cushion < 1.5:
        print("make_protein_whole: WARNING water cushion beyond the protein is "
              "only %.2f nm per side; semiisotropic pressure will shrink box_z "
              "and the protein may contact its periodic image. Increase WATER_NM "
              "or set BOX_Z so box_z >= protein_span + 2*~2.5 nm." % cushion)
    else:
        print("make_protein_whole: water cushion %.2f nm per side beyond the "
              "protein" % cushion)


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
