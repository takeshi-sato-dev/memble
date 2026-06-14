#!/usr/bin/env python3
"""
membrane_props.py

Measure area per lipid and bilayer thickness of a Martini membrane, for checking
a built bilayer against the values reported for its lipid parameters.

Area per lipid is the in-plane box area divided by the number of lipids in one
leaflet. Thickness is the phosphate to phosphate distance, the mean z of the
phosphate beads in the upper leaflet minus that in the lower leaflet. Lipids are
assigned to a leaflet by the z of their head bead relative to the bilayer center.
With a trajectory the values are averaged over frames.

Usage:
  membrane_props.py --gro membrane.gro [--traj eq.xtc --top system.top]
  membrane_props.py --gro membrane.gro --lipids "POPC"          # restrict APL count
"""

import argparse
import sys
import numpy as np

NONLIPID = {"W", "WF", "WAT", "SOL", "NA", "CL", "ION", "NA+", "CL-", "K", "MG",
            "CA", "TIP3", "POT", "SOD", "CLA"}
AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}
PHOSPHATE = ("PO4", "PO1", "PO2", "PO3")
HEAD_PREF = ("PO4", "PO1", "PO2", "NC3", "CNO", "ROH", "GL0", "AM1", "GM1")


def read_gro_frame(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    body = L[2:2 + n]
    box = np.array([float(v) for v in L[2 + n].split()[:3]])
    recs = []
    for b in body:
        resid = b[0:5].strip()
        resn = b[5:10].strip().upper()
        atom = b[10:15].strip()
        x = float(b[20:28]); y = float(b[28:36]); z = float(b[36:44])
        recs.append((resid, resn, atom, x, y, z))
    return recs, box


def lipid_residues(recs, restrict):
    """Group lipid beads into residues keyed by (resid, resname), in order."""
    res = {}
    order = []
    for (resid, resn, atom, x, y, z) in recs:
        if resn in NONLIPID or resn in AA:
            continue
        if restrict and resn not in restrict:
            continue
        key = (resid, resn)
        if key not in res:
            res[key] = []
            order.append(key)
        res[key].append((atom, np.array([x, y, z])))
    return res, order


def head_bead(beads):
    names = {a: xyz for a, xyz in beads}
    for h in HEAD_PREF:
        if h in names:
            return names[h]
    # fall back to the bead with the highest |z| spread end is not robust; use mean
    return np.mean([xyz for _, xyz in beads], axis=0)


def frame_props(recs, box, restrict):
    res, order = lipid_residues(recs, restrict)
    if not res:
        return None
    heads = {k: head_bead(v) for k, v in res.items()}
    zmid = np.mean([h[2] for h in heads.values()])
    upper = [k for k, h in heads.items() if h[2] >= zmid]
    lower = [k for k, h in heads.items() if h[2] < zmid]
    n_leaflet = 0.5 * (len(upper) + len(lower))
    area = box[0] * box[1]
    apl = area / n_leaflet if n_leaflet else float("nan")
    # thickness from phosphate beads only
    def phos_z(keys):
        zs = []
        for k in keys:
            for a, xyz in res[k]:
                if a in PHOSPHATE:
                    zs.append(xyz[2])
        return np.mean(zs) if zs else None
    zu = phos_z(upper); zl = phos_z(lower)
    thick = (zu - zl) if (zu is not None and zl is not None) else float("nan")
    return apl, thick, len(upper), len(lower)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--traj", default="")
    ap.add_argument("--top", default="")
    ap.add_argument("--lipids", default="", help="restrict APL count to these resnames")
    ap.add_argument("--skip", type=int, default=0,
                    help="discard this many initial frames before averaging "
                         "(use to drop the unequilibrated start)")
    args = ap.parse_args()
    restrict = set(s.upper()[:5] for s in args.lipids.split()) if args.lipids else None

    if args.traj:
        try:
            import mdtraj as md
        except ImportError:
            sys.exit("mdtraj is needed for --traj")
        # mdtraj reads a structure file as topology, not a GROMACS .top. If
        # --top is a structure mdtraj understands, use it; otherwise fall back
        # to the gro, so passing system.top by habit still works.
        md_top_ext = (".gro", ".pdb", ".pdb.gz", ".h5", ".lh5", ".prmtop",
                      ".parm7", ".prm7", ".psf", ".mol2", ".hoomdxml",
                      ".hdf5", ".gsd", ".arc")
        if args.top and args.top.lower().endswith(md_top_ext):
            md_top = args.top
        else:
            md_top = args.gro
        t = md.load(args.traj, top=md_top)
        apls = []; thicks = []
        top = t.topology
        # precompute per-residue atom indices for lipids
        start = max(0, args.skip)
        if start >= t.n_frames:
            sys.exit("--skip %d discards all %d frames" % (args.skip, t.n_frames))
        for f in range(start, t.n_frames):
            recs = []
            box = t.unitcell_lengths[f]
            xyz = t.xyz[f]
            for atom in top.atoms:
                r = atom.residue
                recs.append((str(r.resSeq), r.name.upper(), atom.name,
                             xyz[atom.index][0], xyz[atom.index][1], xyz[atom.index][2]))
            pr = frame_props(recs, box, restrict)
            if pr:
                apls.append(pr[0]); thicks.append(pr[1])
        apls = np.array(apls); thicks = np.array(thicks)
        print("frames averaged: %d (skipped first %d)" % (len(apls), start))
        print("area per lipid: %.4f +/- %.4f nm^2" % (np.nanmean(apls), np.nanstd(apls)))
        print("thickness (P-P): %.4f +/- %.4f nm" % (np.nanmean(thicks), np.nanstd(thicks)))
    else:
        recs, box = read_gro_frame(args.gro)
        pr = frame_props(recs, box, restrict)
        if not pr:
            sys.exit("no lipids found; check --lipids or the file")
        apl, thick, nu, nl = pr
        print("single frame")
        print("lipids per leaflet: upper %d, lower %d" % (nu, nl))
        print("area per lipid: %.4f nm^2" % apl)
        print("thickness (P-P): %.4f nm" % thick)
        print("compare with the values reported for these Martini 3 parameters")


if __name__ == "__main__":
    main()
