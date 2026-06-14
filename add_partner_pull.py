#!/usr/bin/env python3
"""
add_partner_pull.py

Add a one-sided flat-bottom pull restraint (MEMB vs PARTNER centre-of-mass along
z) to the production mdp, and a PARTNER index group, so that during production
the peripheral protein:
  * can freely approach and associate with its leaflet (separation may DECREASE),
  * is gently pushed back if it drifts too far / tries to wrap through the
    periodic boundary to the other leaflet (separation cannot exceed init),
giving full lateral and rotational freedom while making backside crossing
impossible. During equilibration the partner is instead frozen by strong
position restraints (handled separately), so it never touches the relaxing
membrane.

Usage:
  add_partner_pull.py --gro system.gro --ndx index.ndx --mdp step7_production.mdp
                      --partner-name PART0 [--margin 1.0] [--k 1000]
"""

import argparse
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--ndx", required=True)
    ap.add_argument("--mdp", required=True)
    ap.add_argument("--partner-name", required=True)
    ap.add_argument("--lipids", default="")
    ap.add_argument("--margin", type=float, default=1.0,
                    help="extra nm the partner COM may move outward before the "
                         "flat-bottom wall engages")
    ap.add_argument("--k", type=float, default=1000.0,
                    help="flat-bottom force constant kJ/mol/nm^2")
    args = ap.parse_args()

    L = open(args.gro).read().splitlines()
    n = int(L[1]); body = L[2:2 + n]
    lipset = set(args.lipids.split())

    part_idx = []
    part_z = []
    memb_z = []
    pname5 = args.partner_name[:5]
    for i, ln in enumerate(body, 1):
        rn = ln[5:10].strip()
        z = float(ln[36:44])
        if rn == pname5:
            part_idx.append(i); part_z.append(z)
        elif rn in lipset:
            memb_z.append(z)
    if not part_idx:
        print("add_partner_pull: no partner beads found; skipping")
        return

    memb_com_z = float(np.mean(memb_z))
    part_com_z = float(np.mean(part_z))
    init = abs(part_com_z - memb_com_z) + args.margin   # max allowed |z-sep|

    # --- append PARTNER group to index.ndx ---
    with open(args.ndx, "a") as fh:
        fh.write("[ PARTNER ]\n")
        for k in range(0, len(part_idx), 15):
            fh.write(" ".join("%d" % x for x in part_idx[k:k + 15]) + "\n")
        fh.write("\n")

    # --- append pull (flat-bottom) to production mdp ---
    pull = [
        "",
        "; --- peripheral-protein backside guard (one-sided flat-bottom) ---",
        "pull                 = yes",
        "pull-ngroups         = 2",
        "pull-group1-name     = MEMB",
        "pull-group2-name     = PARTNER",
        "pull-ncoords         = 1",
        "pull-coord1-type     = flat-bottom",   # force only when coord > init
        "pull-coord1-geometry = distance",
        "pull-coord1-dim      = N N Y",          # z-separation only
        "pull-coord1-groups   = 1 2",
        "pull-coord1-init     = %.3f" % init,
        "pull-coord1-k        = %.1f" % args.k,
        "pull-coord1-start    = no",
        "pull-pbc-ref-prev    = no",
        "",
    ]
    with open(args.mdp, "a") as fh:
        fh.write("\n".join(pull) + "\n")

    print("add_partner_pull: PARTNER group (%d beads), flat-bottom |z-sep| <= "
          "%.2f nm (k=%.0f) added to %s" %
          (len(part_idx), init, args.k, args.mdp))


if __name__ == "__main__":
    main()
