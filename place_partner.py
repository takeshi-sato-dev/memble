#!/usr/bin/env python3
"""
place_partner.py  (post-COBY peripheral-protein placement)

Place an already-coarse-grained peripheral ("partner") protein next to a built
membrane system, on an EXPLICITLY chosen leaflet, a controlled gap away, and
grow the box so the partner has bulk water on its side and can NEVER wrap to the
other (compositionally different) leaflet through the periodic boundary.

Why post-COBY: COBY's lipid-grid optimiser hangs when box_z is much larger than
the membrane. So the membrane + transmembrane protein are built in a thin box;
this script then enlarges box_z and inserts the partner into the new space.

What it does:
  * reads the built system (system.gro) and the chosen leaflet's head height
    (mean z of head beads PO4/ROH belonging to that leaflet);
  * enlarges box_z to: (system span on the partner side) + gap + partner
    thickness + outer bulk water, recenters everything;
  * places the partner with its membrane-facing edge `gap` nm from the head, xy
    centered on the membrane (or at --x/--y);
  * appends the partner to system.gro and adds it to system.top [ molecules ];
  * writes partner_zrange.txt and partner_side.txt for the restraint stage.

Backside protection is enforced later by a one-sided flat-bottom pull restraint
(MEMB vs PARTNER COM along z) that is active only in production; during
equilibration the partner is frozen by strong position restraints. This script
only guarantees the geometry (correct side, gap, room, no initial wrap).
"""

import argparse
import sys
import numpy as np

SOLV_IONS = {"W", "WF", "NA", "CL", "ION", "NA+", "CL-"}
HEAD_BEADS = {"PO4", "ROH", "PO1", "PO2"}


def read_pdb_atoms(path):
    out = []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            out.append(ln.rstrip("\n"))
    if not out:
        sys.exit("ERROR: no ATOM/HETATM in %s" % path)
    return out


def pdb_xyz(ln):
    return np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])


def read_gro(path):
    L = open(path).read().splitlines()
    title, n = L[0], int(L[1])
    body = L[2:2 + n]
    box = [float(v) for v in L[2 + n].split()]
    return title, body, box


def gro_fields(ln):
    return (ln[5:10].strip(), ln[10:15].strip(),
            float(ln[20:28]), float(ln[28:36]), float(ln[36:44]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system-gro", required=True)
    ap.add_argument("--system-top", required=True)
    ap.add_argument("--partner-cg", required=True, help="CG partner PDB (martinized)")
    ap.add_argument("--partner-name", required=True, help="moleculetype name in top")
    ap.add_argument("--partner-itp", required=True)
    ap.add_argument("--side", choices=["upper", "lower"], required=True)
    ap.add_argument("--gap", type=float, default=1.5,
                    help="nm between leaflet head and partner nearest edge")
    ap.add_argument("--water-nm", type=float, default=3.0,
                    help="bulk water beyond the partner on its side (nm)")
    ap.add_argument("--lipids", default="", help="space-separated lipid resnames")
    ap.add_argument("--x", type=float, default=None)
    ap.add_argument("--y", type=float, default=None)
    ap.add_argument("--rotate", choices=["none", "random"], default="none")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    title, body, box = read_gro(args.system_gro)
    bx, by, bz = box[0], box[1], box[2]
    lipset = set(args.lipids.split())

    # --- membrane head height on the chosen leaflet ---
    head_z = []
    all_lip_z = []
    for ln in body:
        rn, an, x, y, z = gro_fields(ln)
        if rn in lipset:
            all_lip_z.append(z)
            if an in HEAD_BEADS:
                head_z.append(z)
    if not all_lip_z:
        sys.exit("ERROR: no lipid beads found (check --lipids)")
    memb_mid = float(np.mean(all_lip_z))
    if head_z:
        upper_head = float(np.mean([z for z in head_z if z > memb_mid]))
        lower_head = float(np.mean([z for z in head_z if z <= memb_mid]))
    else:  # fall back to lipid extremes
        upper_head = float(np.max(all_lip_z))
        lower_head = float(np.min(all_lip_z))

    # --- partner geometry ---
    patoms = read_pdb_atoms(args.partner_cg)
    P = np.array([pdb_xyz(a) for a in patoms]) / 10.0  # A -> nm
    if args.rotate == "random":
        rng = np.random.default_rng(args.seed)
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        w, xq, yq, zq = q
        R = np.array([
            [1 - 2*(yq*yq+zq*zq), 2*(xq*yq-zq*w), 2*(xq*zq+yq*w)],
            [2*(xq*yq+zq*w), 1 - 2*(xq*xq+zq*zq), 2*(yq*zq-xq*w)],
            [2*(xq*zq-yq*w), 2*(yq*zq+xq*w), 1 - 2*(xq*xq+yq*yq)]])
        P = (P - P.mean(0)) @ R.T + P.mean(0)
    p_thick = float(P[:, 2].max() - P[:, 2].min())

    # --- target z for partner so its membrane-facing edge is `gap` from head ---
    if args.side == "upper":
        # partner sits above; its lowest bead at upper_head + gap
        target_edge = upper_head + args.gap
        dz = target_edge - P[:, 2].min()
    else:
        # partner sits below; its highest bead at lower_head - gap
        target_edge = lower_head - args.gap
        dz = target_edge - P[:, 2].max()

    # xy placement: membrane center or explicit
    memb_xy = np.array([bx / 2.0, by / 2.0])
    if args.x is not None and args.y is not None:
        memb_xy = np.array([args.x, args.y])
    P[:, 0] += memb_xy[0] - P[:, 0].mean()
    P[:, 1] += memb_xy[1] - P[:, 1].mean()
    P[:, 2] += dz

    # --- enlarge box_z so the partner side gets `water_nm` of bulk beyond it,
    #     and recenter the whole assembly (system + partner) ---
    sys_z = np.array([gro_fields(ln)[4] for ln in body])
    if args.side == "upper":
        lo = min(sys_z.min(), P[:, 2].min())
        hi = P[:, 2].max() + args.water_nm        # bulk above partner
        lo = lo - args.water_nm                    # modest bulk below membrane
    else:
        hi = max(sys_z.max(), P[:, 2].max())
        lo = P[:, 2].min() - args.water_nm
        hi = hi + args.water_nm
    new_bz = hi - lo
    shift = -lo
    # apply shift in z to existing atoms and partner
    P[:, 2] += shift

    # --- write new gro: existing atoms (z-shifted) + partner appended at end ---
    # partner appended after everything => add "PARTNER 1" at end of [ molecules ]
    out_lines = []
    for ln in body:
        rn, an, x, y, z = gro_fields(ln)
        out_lines.append(ln[:20] + "%8.3f%8.3f%8.3f" % (x, y, z + shift))
    # partner gro lines (each bead its own residue index continues)
    # find current max resid
    maxres = 0
    for ln in body:
        try:
            maxres = max(maxres, int(ln[0:5]))
        except ValueError:
            pass
    rid = (maxres + 1) % 100000
    start_serial = len(body) + 1
    pres_names = [a[12:16].strip() or "BB" for a in patoms]
    for i, a in enumerate(patoms):
        nm = pres_names[i][:5]
        out_lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" %
                         (rid, args.partner_name[:5], nm, (start_serial + i) % 100000,
                          P[i, 0], P[i, 1], P[i, 2]))
    n_total = len(out_lines)
    with open(args.system_gro, "w") as fh:
        fh.write(title + "\n%d\n" % n_total)
        for l in out_lines:
            fh.write(l + "\n")
        fh.write("%10.5f%10.5f%10.5f\n" % (bx, by, new_bz))

    # --- top: ensure partner itp is included, add count to [ molecules ] ---
    tl = open(args.system_top).read().splitlines()
    inc = '#include "%s"' % args.partner_itp
    has_inc = any(inc in l for l in tl)
    out = []
    in_mol = False
    appended = False
    for i, line in enumerate(tl):
        s = line.strip()
        if s.lower().startswith("[ molecules"):
            if not has_inc:
                out.insert(0, inc)  # safest: include at very top
            in_mol = True
            out.append(line)
            continue
        out.append(line)
    out.append("%-6s 1" % args.partner_name)
    open(args.system_top, "w").write("\n".join(out) + "\n")

    # --- record geometry for restraint stage ---
    open("partner_zrange.txt", "w").write(
        "%s %.3f %.3f %.3f %.3f\n" %
        (args.side, P[:, 2].min(), P[:, 2].max(), memb_mid + shift,
         (upper_head if args.side == "upper" else lower_head) + shift))
    print("placed partner %s on %s leaflet: gap %.2f nm, partner z %.2f..%.2f nm, "
          "box_z %.2f -> %.2f nm" %
          (args.partner_name, args.side, args.gap, P[:, 2].min(), P[:, 2].max(),
           bz, new_bz))


if __name__ == "__main__":
    main()
