#!/usr/bin/env python3
"""
prebuild_orient.py

Prepare a pre-assembled multi-chain all-atom PDB (for example four EGFR TM-JM
peptides already arranged by a prior equilibration) for membrane building,
WITHOUT disturbing their relative arrangement.

Two per-chain selections, both interpreted on the input PDB resSeq (numbering is
never changed):

  --keep   "A:54-103;B:54-103;C:50-99;D:54-103"
           residues to keep in the system (everything else is trimmed)

  --core   "A:60-82;B:60-82;C:56-78;D:60-82"
           the TM-core residues whose combined centroid is moved to z = 0 and
           whose principal axis is aligned to z. If omitted, the kept residues
           are used.

The whole assembly is moved as ONE rigid body, so the inter-chain geometry is
preserved exactly. resSeq is left untouched (remap later, at analysis time).
"""

import argparse
import sys

import numpy as np


def parse_ranges(spec):
    """'A:54-103;B:50-99' -> {'A': (54,103), 'B': (50,99)}."""
    out = {}
    if not spec:
        return out
    for entry in spec.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        chain, rng = entry.split(":")
        lo, hi = rng.split("-")
        out[chain.strip()] = (int(lo), int(hi))
    return out


def read_atoms(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                rows.append(ln.rstrip("\n"))
    if not rows:
        sys.exit("ERROR: no ATOM/HETATM records in %s" % path)
    return rows


def chain_of(ln):
    return ln[21]


def resseq_of(ln):
    return int(ln[22:26])


def xyz_of(ln):
    return (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))


def set_xyz(ln, x, y, z):
    return "%s%8.3f%8.3f%8.3f%s" % (ln[:30], x, y, z, ln[54:])


def in_range(ln, ranges):
    c = chain_of(ln)
    if c not in ranges:
        return False
    lo, hi = ranges[c]
    return lo <= resseq_of(ln) <= hi


def principal_axis(coords):
    c = coords - coords.mean(axis=0)
    cov = c.T @ c
    w, v = np.linalg.eigh(cov)
    return v[:, np.argmax(w)]


def bundle_normal(core_atoms):
    """Membrane normal for the assembly.

    For three or more chains the membrane normal is the normal of the plane that
    best fits the per-chain TM-core centroids: those centroids lie in the
    membrane midplane, so their plane normal IS the membrane normal. The TM
    helices themselves are usually tilted relative to the normal (EGFR TM tilt is
    real), so aligning the average HELIX axis to z instead rotates the whole
    laterally-extended bundle and lifts one side, spreading the cores in z. Using
    the centroid plane keeps the cores coplanar and preserves the helix tilt.

    With only one or two chains a plane is ill-defined, so fall back to the
    sign-aligned average of each chain helix axis.
    """
    from collections import defaultdict
    by_chain = defaultdict(list)
    for a in core_atoms:
        by_chain[chain_of(a)].append(xyz_of(a))
    centroids = [np.mean(np.array(pts), axis=0) for pts in by_chain.values()
                 if len(pts) >= 2]
    if len(centroids) >= 3:
        P = np.array(centroids)
        P = P - P.mean(axis=0)
        # plane normal = direction of least variance of the centroid cloud
        _u, _s, vt = np.linalg.svd(P, full_matrices=False)
        return vt[-1] / np.linalg.norm(vt[-1])
    axes = []
    for ch, pts in by_chain.items():
        if len(pts) >= 2:
            axes.append(principal_axis(np.array(pts)))
    if not axes:
        return principal_axis(np.array([xyz_of(a) for a in core_atoms]))
    ref = axes[0]
    aligned = [ax if float(np.dot(ax, ref)) >= 0 else -ax for ax in axes]
    m = np.mean(aligned, axis=0)
    nrm = np.linalg.norm(m)
    if nrm < 1e-8:
        return ref
    return m / nrm


def rotation_to_z(axis):
    """Rotation matrix sending unit `axis` onto +z."""
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
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--keep", required=True, help='e.g. "A:54-103;B:54-103"')
    p.add_argument("--core", default="", help='TM core, e.g. "A:60-82;B:60-82"')
    p.add_argument("--no-rotate", action="store_true",
                   help="translate the core to z=0 but do not align its axis")
    args = p.parse_args()

    keep = parse_ranges(args.keep)
    core = parse_ranges(args.core) if args.core else keep
    if not keep:
        sys.exit("ERROR: --keep is empty")

    atoms = read_atoms(args.inp)
    kept = [a for a in atoms if in_range(a, keep)]
    if not kept:
        sys.exit("ERROR: nothing kept; check --keep chains/ranges against the PDB")

    core_atoms = [a for a in kept if in_range(a, core)]
    if not core_atoms:
        sys.exit("ERROR: no atoms in --core range")

    kept_xyz = np.array([xyz_of(a) for a in kept])
    core_xyz = np.array([xyz_of(a) for a in core_atoms])

    # 1) rigid rotation: average per-chain TM-helix axis -> z (stands the
    #    helices up; pooled principal axis would lay a side-by-side bundle flat)
    if args.no_rotate:
        R = np.eye(3)
    else:
        R = rotation_to_z(bundle_normal(core_atoms))

    pivot = kept_xyz.mean(axis=0)
    kept_rot = (kept_xyz - pivot) @ R.T + pivot
    core_rot = (core_xyz - pivot) @ R.T + pivot

    # 2) rigid translation: core centroid -> z = 0 (x,y left as is)
    dz = -core_rot[:, 2].mean()
    kept_rot[:, 2] += dz

    out = [set_xyz(a, kept_rot[i][0], kept_rot[i][1], kept_rot[i][2])
           for i, a in enumerate(kept)]

    with open(args.out, "w") as fh:
        for ln in out:
            fh.write(ln + "\n")
        fh.write("END\n")

    chains = sorted(set(chain_of(a) for a in kept))
    core_z = (core_rot[:, 2] + dz)
    print("prebuild: kept %d atoms across chains %s (resSeq preserved); "
          "TM-core centroid z = %.3f nm after centering"
          % (len(kept), ",".join(chains), core_z.mean() / 10.0))


if __name__ == "__main__":
    main()
