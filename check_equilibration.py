#!/usr/bin/env python3
"""
check_equilibration.py

Say whether a membrane has stopped changing, instead of assuming that a run of a
stated length was long enough.

The area of a bilayer relaxes on its own time scale, and that time scale depends
on the composition, on the sterol content, and on the protein in the membrane. A
run of five nanoseconds equilibrates one system and leaves another still moving,
so the length of the run answers nothing on its own. This script reads the area
and the surface tension from the energy file, fits the second half of the run,
and reports whether the area is still drifting and whether the bilayer is under
tension.

  area          Lx * Ly from the energy file, in nm^2
  drift         the slope of the second half, as a percentage of the area per ns
  tension       the mean surface tension of the second half, in bar nm

A membrane that has equilibrated shows a drift that the scatter of the run
cannot separate from zero, and a surface tension that covers zero. A membrane
that fails either one is reported with the length of run its own drift implies.

Usage:
  check_equilibration.py --edr step6.6.edr [--gmx gmx]
  check_equilibration.py --xvg box.xvg --tension-xvg surften.xvg
  check_equilibration.py --edr step7.edr --lipids-per-leaflet 138 --last 0.5

Exit code is non-zero when the membrane has not settled, so a build script can
hold the production run back.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def read_xvg(path):
    t, cols = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#@":
                continue
            parts = line.split()
            try:
                vals = [float(v) for v in parts]
            except ValueError:
                continue
            t.append(vals[0])
            cols.append(vals[1:])
    if not t:
        return np.array([]), np.zeros((0, 0))
    return np.array(t), np.array(cols)


def gmx_energy(gmx, edr, terms, outdir):
    """Extract named energy terms into an xvg. Returns (time, {term: values})."""
    out = os.path.join(outdir, "energy.xvg")
    sel = "\n".join(terms) + "\n\n"
    proc = subprocess.run([gmx, "energy", "-f", edr, "-o", out],
                          input=sel, text=True, capture_output=True)
    if not os.path.isfile(out):
        sys.stderr.write(proc.stderr[-2000:] + "\n")
        return None, None
    t, cols = read_xvg(out)
    if cols.size == 0:
        return None, None
    # gmx writes the selected terms in the order it holds them, and it names
    # them in the legend. Read the legend rather than trusting the order.
    legend = []
    with open(out) as fh:
        for line in fh:
            if line.startswith("@ s") and "legend" in line:
                legend.append(line.split('"')[1])
            elif not line.startswith(("#", "@")):
                break
    got = {}
    for i, name in enumerate(legend):
        if i < cols.shape[1]:
            got[name] = cols[:, i]
    return t, got


def fit_drift(t, y):
    """Slope of a straight line, and the standard error of that slope."""
    n = len(t)
    if n < 4:
        return 0.0, float("inf")
    a, b = np.polyfit(t, y, 1)
    resid = y - (a * t + b)
    s2 = float((resid ** 2).sum()) / (n - 2)
    sxx = float(((t - t.mean()) ** 2).sum())
    err = float(np.sqrt(s2 / sxx)) if sxx > 0 else float("inf")
    return float(a), err


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edr", default="", help="GROMACS energy file")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--xvg", default="", help="xvg holding Box-X and Box-Y")
    ap.add_argument("--tension-xvg", default="", help="xvg holding the surface tension")
    ap.add_argument("--last", type=float, default=0.5,
                    help="fraction of the run to fit (default the second half)")
    ap.add_argument("--drift-tol", type=float, default=0.5,
                    help="largest area drift accepted, in percent per ns (default 0.5)")
    ap.add_argument("--tension-tol", type=float, default=50.0,
                    help="largest surface tension accepted, in bar nm (default 50)")
    ap.add_argument("--lipids-per-leaflet", type=int, default=0,
                    help="report the area per lipid as well as the area")
    ap.add_argument("--json", default="")
    args = ap.parse_args()

    t = None
    area = None
    tension = None

    if args.edr:
        if not shutil.which(args.gmx):
            sys.exit("ERROR: %s was not found.\n\n  What to do:\n"
                     "    Set --gmx to the GROMACS binary, for example\n"
                     "      --gmx /usr/local/gromacs/bin/gmx\n"
                     "    or extract the terms yourself and pass them with --xvg:\n"
                     "      gmx energy -f %s -o box.xvg   (select Box-X and Box-Y)"
                     % (args.gmx, args.edr))
        with tempfile.TemporaryDirectory() as d:
            t, got = gmx_energy(args.gmx, args.edr,
                                ["Box-X", "Box-Y", "#Surf*SurfTen"], d)
            if t is None:
                sys.exit("ERROR: no energy term could be read from %s.\n\n"
                         "  What to do:\n"
                         "    List what the file holds and pass the terms yourself:\n"
                         "      gmx energy -f %s" % (args.edr, args.edr))
            if "Box-X" in got and "Box-Y" in got:
                area = got["Box-X"] * got["Box-Y"]
            for k in got:
                if "Surf" in k:
                    tension = got[k]
    elif args.xvg:
        t, cols = read_xvg(args.xvg)
        if cols.shape[1] >= 2:
            area = cols[:, 0] * cols[:, 1]
        elif cols.shape[1] == 1:
            area = cols[:, 0]
        if args.tension_xvg:
            _, tc = read_xvg(args.tension_xvg)
            if tc.size:
                tension = tc[:, 0]
    else:
        sys.exit("ERROR: give --edr or --xvg.")

    if area is None or len(area) < 4:
        sys.exit("ERROR: the area could not be read.\n\n  What to do:\n"
                 "    The energy file has to hold Box-X and Box-Y. A run written\n"
                 "    with nstenergy = 0 holds nothing; set nstenergy in the mdp\n"
                 "    and run the stage again.")

    n0 = int(len(t) * (1.0 - args.last))
    tt, aa = t[n0:], area[n0:]
    span_ns = (tt[-1] - tt[0]) / 1000.0 if tt[-1] > 100 else (tt[-1] - tt[0])
    unit = "ns" if tt[-1] > 100 else "ns"
    # gmx writes time in ps.
    tt_ns = tt / 1000.0

    slope, slope_err = fit_drift(tt_ns, aa)
    mean_area = float(np.mean(aa))
    drift_pct = 100.0 * slope / mean_area if mean_area else 0.0
    drift_err = 100.0 * slope_err / mean_area if mean_area else float("inf")

    print("area of the membrane over the last %.0f%% of the run" % (100 * args.last))
    print("  interval fitted     %.2f to %.2f ns" % (tt_ns[0], tt_ns[-1]))
    print("  mean area           %.2f nm^2" % mean_area)
    if args.lipids_per_leaflet:
        print("  area per lipid      %.4f nm^2"
              % (mean_area / args.lipids_per_leaflet))
    print("  drift               %+.3f %% per ns (standard error %.3f)"
          % (drift_pct, drift_err))

    settled = abs(drift_pct) <= max(args.drift_tol, 2.0 * drift_err)
    result = {"mean_area_nm2": round(mean_area, 3),
              "drift_percent_per_ns": round(drift_pct, 4),
              "drift_standard_error": round(drift_err, 4),
              "fitted_from_ns": round(float(tt_ns[0]), 3),
              "fitted_to_ns": round(float(tt_ns[-1]), 3)}

    tension_ok = True
    if tension is not None and len(tension) > n0:
        ten = tension[n0:]
        mt = float(np.mean(ten))
        se = float(np.std(ten) / np.sqrt(max(len(ten), 1)))
        print("  surface tension     %+.1f bar nm (standard error %.1f)" % (mt, se))
        result["surface_tension_bar_nm"] = round(mt, 2)
        result["surface_tension_standard_error"] = round(se, 2)
        tension_ok = abs(mt) <= max(args.tension_tol, 2.0 * se)

    print("")
    if settled and tension_ok:
        result["result"] = "PASS"
        print("The membrane has settled. The area drift is within the scatter of "
              "the run, and the bilayer carries no tension that the run can "
              "separate from zero.")
    else:
        result["result"] = "FAIL"
        if not settled:
            # How much longer the area needs, if it keeps the drift it has.
            need = abs(drift_pct) / max(args.drift_tol, 1e-6)
            print("The area is still changing at %+.2f %% per ns, which the run "
                  "separates from zero." % drift_pct)
            print("")
            print("  What to do:")
            print("    1. run this stage for about %.0f times longer, and read it"
                  % need)
            print("       again. The drift falls as the membrane relaxes.")
            print("    2. a drift that does not fall is not a matter of length. It")
            print("       means the two leaflets hold areas that do not match, so")
            print("       run leaflet_area_check.py on the starting structure and")
            print("       correct the counts before running again.")
            print("    3. a composition whose lipids demix relaxes slowly by")
            print("       nature; report the drift rather than waiting it out.")
        if not tension_ok:
            print("")
            print("  The bilayer is under a surface tension that the run separates")
            print("  from zero. Under semi-isotropic coupling a tensionless bilayer")
            print("  is the intended state.")
            print("  What to do:")
            print("    1. check that the mdp couples x and y separately from z")
            print("       (semiisotropic), and that ref-p holds one value for each")
            print("    2. check the leaflet areas; a mismatch shows up here first")

    if args.json:
        import json
        with open(args.json, "w") as fh:
            json.dump(result, fh, indent=2)

    return 0 if (settled and tension_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
