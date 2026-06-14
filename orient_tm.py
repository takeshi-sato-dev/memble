#!/usr/bin/env python3
"""
orient_tm.py

Given an all-atom protein PDB, orient it so the transmembrane (TM) helix lies
along z and is centered at the origin, then write the oriented PDB. This removes
the manual "pre-orient the peptide" prerequisite: feed any AA structure and the
membrane builder receives a correctly oriented protein.

TM span detection (in order of preference):
  1. --tm-range START:END   explicit residue range (resSeq), if you know it
  2. --ss STRING            secondary-structure string (longest 'H' run = TM)
  3. automatic              Kyte-Doolittle hydrophobicity, most hydrophobic
                            window of --tm-window residues among CA atoms

The TM axis is the first principal component of the TM CA coordinates; the whole
structure is rotated to bring that axis onto z (Rodrigues rotation) and
translated so the TM centroid sits at the origin.

No external deps beyond numpy. Operates on ATOM/HETATM records only.
"""

import argparse
import sys
import numpy as np

# Kyte-Doolittle hydropathy (3-letter), higher = more hydrophobic
KD = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9, "TYR": -1.3,
    "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
    "LYS": -3.9, "ARG": -4.5,
}


def read_atoms(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append(line.rstrip("\n"))
    if not atoms:
        sys.exit("ERROR: no ATOM/HETATM records in %s" % path)
    return atoms


def xyz(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def set_xyz(line, v):
    return "%s%8.3f%8.3f%8.3f%s" % (line[:30], v[0], v[1], v[2], line[54:])


def resseq(line):
    return int(line[22:26])


def resname(line):
    return line[17:20].strip()


def atomname(line):
    return line[12:16].strip()


def ca_table(atoms):
    """Ordered unique residues by their CA atom: list of (resSeq, resName, coord)."""
    out = []
    seen = set()
    for ln in atoms:
        if atomname(ln) == "CA":
            rs = resseq(ln)
            if rs not in seen:
                seen.add(rs)
                out.append((rs, resname(ln), xyz(ln)))
    if not out:
        sys.exit("ERROR: no CA atoms found; is this an all-atom protein PDB?")
    return out


def tm_indices_from_ss(ss, cas):
    if len(ss) != len(cas):
        sys.exit("ERROR: --ss length %d != number of residues %d"
                 % (len(ss), len(cas)))
    best_i = best_len = cur_i = cur_len = 0
    for i, c in enumerate(ss):
        if c in "Hh":
            if cur_len == 0:
                cur_i = i
            cur_len += 1
            if cur_len > best_len:
                best_len, best_i = cur_len, cur_i
        else:
            cur_len = 0
    if best_len == 0:
        sys.exit("ERROR: no helix ('H') found in --ss string")
    return list(range(best_i, best_i + best_len))


def tm_indices_auto(cas, window):
    if len(cas) < window:
        window = len(cas)
    scores = [KD.get(rn, 0.0) for _, rn, _ in cas]
    best_i, best_sum = 0, -1e9
    for i in range(0, len(cas) - window + 1):
        s = sum(scores[i:i + window])
        if s > best_sum:
            best_sum, best_i = s, i
    return list(range(best_i, best_i + window))


def tm_indices_from_range(rng, cas):
    a, b = rng.split(":")
    a, b = int(a), int(b)
    idx = [i for i, (rs, _, _) in enumerate(cas) if a <= rs <= b]
    if not idx:
        sys.exit("ERROR: no residues in --tm-range %s" % rng)
    return idx


def rotation_to_z(axis):
    """Rotation matrix bringing unit vector `axis` onto +z (Rodrigues)."""
    a = axis / np.linalg.norm(axis)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(a, z)
    s = np.linalg.norm(v)
    c = float(np.dot(a, z))
    if s < 1e-8:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="in_pdb", required=True)
    p.add_argument("--out", dest="out_pdb", required=True)
    p.add_argument("--tm-range", default=None, help="resSeq START:END override")
    p.add_argument("--ss", default=None, help="secondary-structure string")
    p.add_argument("--tm-window", type=int, default=21,
                   help="residues in the auto hydrophobic window (default 21)")
    args = p.parse_args()

    atoms = read_atoms(args.in_pdb)
    cas = ca_table(atoms)

    if args.tm_range:
        tm = tm_indices_from_range(args.tm_range, cas)
        how = "range %s" % args.tm_range
    elif args.ss:
        tm = tm_indices_from_ss(args.ss, cas)
        how = "ss longest-helix"
    else:
        tm = tm_indices_auto(cas, args.tm_window)
        how = "auto KD window=%d" % args.tm_window

    tm_coords = np.array([cas[i][2] for i in tm])
    tm_centroid = tm_coords.mean(axis=0)

    # principal axis of the TM CA cloud
    centered = tm_coords - tm_centroid
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    axis = vt[0]

    R = rotation_to_z(axis)

    # apply: rotate all atoms about TM centroid, then put TM centroid at origin
    out = []
    for ln in atoms:
        v = R @ (xyz(ln) - tm_centroid)
        out.append(set_xyz(ln, v))

    with open(args.out_pdb, "w") as fh:
        fh.write("REMARK  oriented by orient_tm.py; TM span %d res (%s)\n"
                 % (len(tm), how))
        for ln in out:
            fh.write(ln + "\n")
        fh.write("END\n")

    # report resulting TM axis tilt (should be ~0 deg from z)
    new_tm = np.array([R @ (cas[i][2] - tm_centroid) for i in tm])
    _, _, vt2 = np.linalg.svd(new_tm - new_tm.mean(axis=0), full_matrices=False)
    tilt = np.degrees(np.arccos(abs(vt2[0][2])))
    print("oriented %s -> %s | TM detect: %s | %d residues | residual tilt %.2f deg"
          % (args.in_pdb, args.out_pdb, how, len(tm), tilt))


if __name__ == "__main__":
    main()
