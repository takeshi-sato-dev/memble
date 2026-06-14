#!/usr/bin/env python3
"""
balance_apl.py

Compute a per-leaflet area-per-lipid for COBY so an asymmetric bilayer comes out
with matched leaflet areas. The user gives only composition ratios; different
lipids occupy different areas, so passing one apl to both leaflets (as COBY does
by default) leaves the leaflet with larger lipids over-area and the bilayer under
stress. This sizes each leaflet apl from its composition, automating the manual
trial-and-error of balancing leaflet counts.

The leaflet mean area is the ratio-weighted sum of per-lipid areas. The two
leaflet apls are scaled to keep their average equal to the requested packing apl
(which carries headroom against overlaps), while their ratio reflects the real
composition areas, so the leaflets end up with equal total area and the packing
density is preserved.

Prints: "<apl_upper> <apl_lower>"

Usage:
  balance_apl.py --upper "CHOL:1 DLPC:1 PSM:1" --lower "CHOL:1 DLPC:1 DOPS:1" \
                 --base-apl 0.70 [--apl "PSM:0.47 ..."]
"""

import argparse

# Approximate Martini 3 area-per-lipid (nm^2). Coarse reference; override with
# --apl. Kept in sync with leaflet_area_check.py.
APL = {
    "CHOL": 0.40, "POPC": 0.64, "DOPC": 0.67, "DIPC": 0.68, "DLPC": 0.60,
    "DPPC": 0.50, "DPSM": 0.50, "PPSM": 0.50, "PSM": 0.50, "POPE": 0.59,
    "DOPE": 0.61, "POPS": 0.55, "POPG": 0.62, "POPA": 0.55, "POPI": 0.66,
    "PAPC": 0.68, "PUPC": 0.70, "DPG3": 0.75, "DPGS": 0.70,
    "DOPS": 0.60, "DOPG": 0.64, "DOPA": 0.58, "DPPE": 0.50, "DPPS": 0.48,
    "DPPG": 0.52, "DLPS": 0.55, "DSPC": 0.50, "PIPC": 0.66, "PGPC": 0.66,
}
DEFAULT_APL = 0.60


def mean_area(comp):
    tot_r = 0.0
    tot_a = 0.0
    for tok in comp.split():
        parts = tok.split(":")
        name = parts[0]
        ratio = float(parts[1]) if len(parts) > 1 else 1.0
        tot_r += ratio
        tot_a += ratio * APL.get(name.upper(), DEFAULT_APL)
    return tot_a / tot_r if tot_r else DEFAULT_APL


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--upper", required=True)
    ap.add_argument("--lower", required=True)
    ap.add_argument("--base-apl", type=float, default=0.70)
    ap.add_argument("--apl", default="")
    args = ap.parse_args()

    for tok in args.apl.split():
        if ":" in tok:
            k, v = tok.split(":")
            APL[k.upper()] = float(v)

    a_up = mean_area(args.upper)
    a_lo = mean_area(args.lower)
    a_mean = 0.5 * (a_up + a_lo)
    # scale so the average apl stays at the requested packing value, only the
    # leaflet ratio changes
    apl_up = args.base_apl * a_up / a_mean
    apl_lo = args.base_apl * a_lo / a_mean
    print("%.4f %.4f" % (apl_up, apl_lo))


if __name__ == "__main__":
    main()
