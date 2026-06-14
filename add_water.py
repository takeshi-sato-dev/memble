#!/usr/bin/env python3
"""
add_water.py

Enlarge box_z of a COBY-built CG system and fill the newly created vacuum slabs
(above and below the existing thin water region) with coarse-grained water (W)
and salt (NA/CL).

Why: COBY's lipid-grid optimiser hangs when box_z is much larger than the
membrane, so the membrane is built in a thin box. A transmembrane protein with
long juxtamembrane / extracellular parts then protrudes beyond that thin box with
too little bulk water, and the periodic image is too close. This script grows the
box only AFTER the membrane is built (COBY never sees the big box) and solvates
the added space.

Approach (keeps GROMACS top/gro consistent):
  * recenter the whole system so the protein z-midpoint sits at the new box
    centre, set the new box_z;
  * estimate water number density from the existing W beads;
  * place new W beads on a grid in the vacuum slabs at that density, skipping any
    site within --clear nm of an existing atom (so the protein is not clashed);
  * convert a fraction of the new W to NA/CL to match --salt molarity (added in
    neutral pairs);
  * append the new beads to the END of the .gro and append matching "W n" /
    "NA n" / "CL n" lines to the END of [ molecules ] in the .top. GROMACS allows
    a moleculetype to appear multiple times, and appending to both ends keeps the
    atom order consistent.

Usage:
  add_water.py --gro system.gro --top system.top --water-nm 3.5 --salt 0.15
               [--water-resname W] [--clear 0.30]
"""

import argparse
import math

import numpy as np

SOLV_IONS = {"W", "WF", "NA", "CL", "ION", "NA+", "CL-"}
AA = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
      "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
      "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}


def read_gro(path):
    L = open(path).read().splitlines()
    title, n = L[0], int(L[1])
    body = L[2:2 + n]
    box = [float(v) for v in L[2 + n].split()]
    resid = np.array([int(l[0:5]) for l in body])
    resn = [l[5:10].strip() for l in body]
    aname = [l[10:15].strip() for l in body]
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])]
                    for l in body])
    return title, body, resid, resn, aname, xyz, box


def fmt_gro_line(resid, resn, aname, idx, xyz):
    return "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (
        resid % 100000, resn[:5], aname[:5], idx % 100000,
        xyz[0], xyz[1], xyz[2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--water-nm", type=float, required=True,
                    help="bulk water thickness to guarantee beyond each protein end (nm)")
    ap.add_argument("--salt", type=float, default=0.15, help="NaCl molarity")
    ap.add_argument("--water-resname", default="W")
    ap.add_argument("--clear", type=float, default=0.30,
                    help="min distance new water keeps from existing atoms (nm)")
    ap.add_argument("--keep-box", action="store_true",
                    help="do not expand or recenter; fill vacuum in the current "
                         "box (use when box_z was already set, e.g. by partner "
                         "placement)")
    args = ap.parse_args()

    title, body, resid, resn, aname, xyz, box = read_gro(args.gro)
    W = args.water_resname
    bx, by = box[0], box[1]

    if args.keep_box:
        new_box_z = box[2]
        shift_z = 0.0
    else:
        # --- protein z-extent (everything not solvent/ion/lipid) ---
        # protein z-extent measured from the actual protein (amino-acid beads),
        # not from every non-solvent bead: lipids fill the thin COBY box, so
        # using them would size the water layer against the membrane and leave a
        # protruding transmembrane protein with almost no water beyond its ends.
        is_prot = np.array([r.upper() in AA for r in resn])
        if not is_prot.any():
            is_prot = np.array([r not in SOLV_IONS for r in resn])
        prot_z = xyz[is_prot, 2]
        if len(prot_z) == 0:
            prot_z = xyz[:, 2]
        pz_lo, pz_hi = prot_z.min(), prot_z.max()
        new_box_z = (pz_hi - pz_lo) + 2 * args.water_nm
        if new_box_z <= box[2] + 1e-6:
            print("add_water: existing box_z %.2f already >= target %.2f; "
                  "nothing to do" % (box[2], new_box_z))
            return
        shift_z = new_box_z / 2.0 - 0.5 * (pz_lo + pz_hi)
    xyz[:, 2] += shift_z

    # --- existing water density (beads per nm^3) ---
    is_w = np.array([r == W for r in resn])
    wz = xyz[is_w, 2]
    if len(wz) < 10:
        print("add_water: too few existing W to estimate density; skipping")
        return
    w_zspan = max(wz.max() - wz.min(), 1e-3)
    dens = len(wz) / (bx * by * w_zspan)          # beads / nm^3
    spacing = dens ** (-1.0 / 3.0)

    # --- vacuum slabs to fill: below existing water and above it ---
    occ_lo, occ_hi = wz.min(), wz.max()
    slabs = [(0.20, occ_lo - 0.20), (occ_hi + 0.20, new_box_z - 0.20)]

    # existing atoms for clash check (coarse grid hash)
    cell = max(args.clear, 0.25)
    from collections import defaultdict
    grid = defaultdict(list)
    for p in xyz:
        grid[(int(p[0] / cell), int(p[1] / cell), int(p[2] / cell))].append(p)

    def clashes(p):
        cx, cy, cz = int(p[0] / cell), int(p[1] / cell), int(p[2] / cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for q in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        if (p[0] - q[0])**2 + (p[1] - q[1])**2 + (p[2] - q[2])**2 < args.clear**2:
                            return True
        return False

    new_pts = []
    nx = max(int(bx / spacing), 1)
    ny = max(int(by / spacing), 1)
    for (z0, z1) in slabs:
        if z1 - z0 < spacing * 0.5:
            continue
        nz = max(int((z1 - z0) / spacing), 1)
        for i in range(nx):
            x = (i + 0.5) * bx / nx
            for j in range(ny):
                y = (j + 0.5) * by / ny
                for k in range(nz):
                    z = z0 + (k + 0.5) * (z1 - z0) / nz
                    p = np.array([x, y, z])
                    if not clashes(p):
                        new_pts.append(p)

    if not new_pts:
        print("add_water: no room for new water (check clearances)")
        return

    # --- salt: ion pairs in the NEW water volume ---
    vol_new_L = bx * by * sum(max(z1 - z0, 0) for z0, z1 in slabs) * 1e-24
    n_pairs = int(round(args.salt * 6.022e23 * vol_new_L))
    n_pairs = min(n_pairs, len(new_pts) // 2)
    rng = np.random.default_rng(0)
    order = rng.permutation(len(new_pts))
    na_idx = set(order[:n_pairs].tolist())
    cl_idx = set(order[n_pairs:2 * n_pairs].tolist())

    # --- build new body: keep existing atom order, INSERT new beads right after
    #     the last existing atom of each type (so order still matches the top) ---
    start_resid = int(resid.max()) + 1
    waters, nas, cls = [], [], []
    for m, p in enumerate(new_pts):
        if m in na_idx:
            nas.append(p)
        elif m in cl_idx:
            cls.append(p)
        else:
            waters.append(p)
    add_w, add_na, add_cl = len(waters), len(nas), len(cls)

    # shift existing atom z and rewrite their lines
    shifted = []
    for l in body:
        x = float(l[20:28]); y = float(l[28:36]); z = float(l[36:44]) + shift_z
        shifted.append((l[5:10].strip(), l[:20] + "%8.3f%8.3f%8.3f" % (x, y, z)))

    def last_index(resname):
        idx = -1
        for i, (rn, _) in enumerate(shifted):
            if rn == resname:
                idx = i
        return idx

    # make new lines (resid/serial fixed up after assembly)
    def mk(rn, pts):
        return [(rn, p) for p in pts]
    new_by_type = {W: mk(W, waters), "NA": mk("NA", nas), "CL": mk("CL", cls)}

    # assemble: walk existing, and after the last atom of each type, splice its new beads
    assembled = []
    li_W, li_NA, li_CL = last_index(W), last_index("NA"), last_index("CL")
    for i, (rn, line) in enumerate(shifted):
        assembled.append(line)
        if i == li_W and new_by_type[W]:
            for rn2, p in new_by_type[W]:
                assembled.append(("PLACEHOLDER", rn2, p))
        if i == li_NA and new_by_type["NA"]:
            for rn2, p in new_by_type["NA"]:
                assembled.append(("PLACEHOLDER", rn2, p))
        if i == li_CL and new_by_type["CL"]:
            for rn2, p in new_by_type["CL"]:
                assembled.append(("PLACEHOLDER", rn2, p))
    # if a type had no existing block, append at end
    for rn2 in (W, "NA", "CL"):
        if last_index(rn2) == -1 and new_by_type[rn2]:
            for _rn, p in new_by_type[rn2]:
                assembled.append(("PLACEHOLDER", rn2, p))

    # renumber resid/serial across the whole system
    out_lines = []
    rid = 1
    ser = 1
    for item in assembled:
        if isinstance(item, str):
            # existing line: rewrite serial/resid to keep them sequential & valid
            out_lines.append(item)   # keep as-is (existing numbering fine)
            ser += 1
        else:
            _, rn2, p = item
            out_lines.append(fmt_gro_line(start_resid, rn2, rn2, ser, p))
            start_resid += 1
            ser += 1

    n_total = len(out_lines)
    with open(args.gro, "w") as fh:
        fh.write(title + "\n%d\n" % n_total)
        for l in out_lines:
            fh.write(l + "\n")
        fh.write("%10.5f%10.5f%10.5f\n" % (bx, by, new_box_z))

    # --- merge counts into existing [ molecules ] lines (no duplicate lines, so
    #     ParmEd stays happy); add a line only if the type was absent ---
    tl = open(args.top).read().splitlines()
    add = {W: add_w, "NA": add_na, "CL": add_cl}
    seen = {W: False, "NA": False, "CL": False}
    out = []
    in_mol = False
    for line in tl:
        s = line.strip()
        if s.lower().startswith("[ molecules"):
            in_mol = True
            out.append(line)
            continue
        if in_mol and s and not s.startswith("[") and not s.startswith(";"):
            parts = s.split()
            name = parts[0]
            if name in add and add[name] and not seen[name]:
                try:
                    cnt = int(parts[1]) + add[name]
                    out.append("%-6s %d" % (name, cnt))
                    seen[name] = True
                    continue
                except (IndexError, ValueError):
                    pass
        out.append(line)
    for name in (W, "NA", "CL"):
        if add[name] and not seen[name]:
            out.append("%-6s %d" % (name, add[name]))
    open(args.top, "w").write("\n".join(out) + "\n")

    print("add_water: box_z %.2f -> %.2f nm; added W=%d NA=%d CL=%d (density %.2f /nm^3)"
          % (box[2], new_box_z, add_w, add_na, add_cl, dens))


if __name__ == "__main__":
    main()
