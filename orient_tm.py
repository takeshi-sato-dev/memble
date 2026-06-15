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


def tm_runs_from_ss(ss, cas, min_len=12):
    """All helical runs (list of index-lists), for multi-pass TM proteins.

    A single-helix TM peptide gives one run; a 7-TM bundle gives ~7. Runs
    shorter than min_len (loops, short 3-10 turns, the amphipathic H8 when it
    is short) are dropped so only membrane-spanning helices set the axis.
    """
    if len(ss) != len(cas):
        sys.exit("ERROR: --ss length %d != number of residues %d"
                 % (len(ss), len(cas)))
    runs = []
    cur = []
    for i, c in enumerate(ss):
        if c in "Hh":
            cur.append(i)
        else:
            if len(cur) >= min_len:
                runs.append(cur)
            cur = []
    if len(cur) >= min_len:
        runs.append(cur)
    return runs


def helix_axis(coords):
    """Unit long axis of one helix = first principal component of its CA cloud."""
    c = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    a = vt[0]
    return a / np.linalg.norm(a)


def bundle_axis(cas, runs):
    """Membrane normal for a multi-helix bundle.

    Each TM helix has a long axis but an arbitrary sign. Anti-parallel packing
    (the rule in 7-TM receptors) means a plain average cancels out, and the
    first principal component of the pooled cloud points along the widest
    in-plane spread, not the normal. Instead take each helix axis, flip all of
    them into a common hemisphere, and average: that mean points along the
    shared spanning direction, i.e. the membrane normal.
    """
    axes = []
    for run in runs:
        coords = np.array([cas[i][2] for i in run])
        axes.append(helix_axis(coords))
    axes = np.array(axes)
    ref = axes[0]
    for k in range(len(axes)):
        if np.dot(axes[k], ref) < 0:
            axes[k] = -axes[k]
    mean = axes.mean(axis=0)
    n = np.linalg.norm(mean)
    if n < 1e-6:
        # degenerate (helices cancel); fall back to pooled principal axis
        pooled = np.array([cas[i][2] for run in runs for i in run])
        return helix_axis(pooled)
    return mean / n


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
    p.add_argument("--multi-tm", action="store_true",
                   help="treat the protein as a multi-pass TM bundle (e.g. a "
                        "7-TM GPCR): detect every helix from --ss and set the "
                        "membrane normal to the mean helix axis, not the first "
                        "principal component of the pooled cloud")
    p.add_argument("--multi-tm-minlen", type=int, default=12,
                   help="minimum helix length (residues) counted as a TM helix "
                        "in --multi-tm mode (default 12)")
    p.add_argument("--nterm-side", choices=["up", "down"], default=None,
                   help="after orienting, flip so the N-terminus faces +z (up) "
                        "or -z (down). The bundle axis sign is otherwise "
                        "arbitrary, so a 7-TM receptor can end up inverted. For "
                        "a GPCR the N-terminus is extracellular: use up if the "
                        "extracellular side should face the upper leaflet.")
    args = p.parse_args()

    atoms = read_atoms(args.in_pdb)
    cas = ca_table(atoms)

    # Decide the orientation axis.
    # Multi-TM (bundle) path: explicit --multi-tm, or --ss that contains several
    # membrane-length helices. Single-TM path: unchanged from the original tool.
    runs = None
    if args.ss and not args.tm_range:
        candidate = tm_runs_from_ss(args.ss, cas, args.multi_tm_minlen)
        if args.multi_tm or len(candidate) >= 2:
            runs = candidate

    if runs:
        tm = sorted(i for run in runs for i in run)
        tm_coords = np.array([cas[i][2] for i in tm])
        tm_centroid = tm_coords.mean(axis=0)
        axis = bundle_axis(cas, runs)
        how = "ss multi-TM bundle (%d helices)" % len(runs)
    else:
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

    # Optionally fix the up/down sense so a chosen terminus faces +z. The bundle
    # axis (and a single-helix principal axis) has an arbitrary sign, so without
    # this a multi-pass receptor can come out inverted (N-terminus on the wrong
    # leaflet). cas is in file order, so cas[0] is the N-terminus.
    if args.nterm_side:
        nz = float((R @ (cas[0][2] - tm_centroid))[2])
        cz = float((R @ (cas[-1][2] - tm_centroid))[2])
        need_flip = (args.nterm_side == "up" and nz < cz) or \
                    (args.nterm_side == "down" and nz > cz)
        if need_flip:
            flip = np.diag([1.0, -1.0, -1.0])   # 180 deg about x: swaps +z/-z, keeps axis vertical
            R = flip @ R

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
    if runs:
        new_axis = R @ axis
        tilt = np.degrees(np.arccos(min(1.0, abs(new_axis[2]))))
    else:
        new_tm = np.array([R @ (cas[i][2] - tm_centroid) for i in tm])
        _, _, vt2 = np.linalg.svd(new_tm - new_tm.mean(axis=0), full_matrices=False)
        tilt = np.degrees(np.arccos(abs(vt2[0][2])))
    print("oriented %s -> %s | TM detect: %s | %d residues | residual tilt %.2f deg"
          % (args.in_pdb, args.out_pdb, how, len(tm), tilt))


if __name__ == "__main__":
    main()
