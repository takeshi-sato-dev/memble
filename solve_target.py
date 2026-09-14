#!/usr/bin/env python3
"""solve_target.py : the leaflet numbers that hold a required sterol asymmetry.

A build sets how many cholesterol molecules each leaflet receives, and the run
does not keep those numbers. What the run keeps is set by how many phospholipid
molecules each leaflet holds. run_curve.sh measures that relation over a range
of phospholipid numbers, and this script reads the relation backwards: given the
share of the sterol the upper leaflet is required to hold, it returns the two
phospholipid numbers and the two sterol numbers a build has to be given.

The build is given four numbers and not two. The two phospholipid numbers are
the quantity the run cannot change, and they are what carries the sterol to the
target. The two sterol numbers are the target written into the build, and they
are what makes the system start where the run would otherwise take it. A build
given the phospholipid numbers alone reaches the target only after the sterol
has moved, and over that part of the run it carries a composition that a Methods
section does not state.

Usage:
    python3 solve_target.py <curve dir> --target 70
    python3 solve_target.py <curve dir> --target 70 --total-pl 1937
    python3 solve_target.py <curve dir> --target 70 --emit-env

<curve dir> holds the relax_d<delta>.json files run_curve.sh writes. Repeats of
a point (relax_d0_s2_250ns.json) are averaged into that point.

The script refuses a target the measured curve does not reach, and says what
range it does reach. Extrapolating a fitted quadratic past its data returns a
number with no measurement behind it, and a leaflet stretched far enough stops
behaving like a bilayer.
"""
import argparse
import glob
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

STEROL = "CHOL"


def read_points(d, prefix="relax_d", last_fraction=0.4):
    """Return [(imbalance %, settled sterol share %, n runs, sd)] per delta."""
    pat = re.compile(re.escape(prefix) + r"(-?\d+)(?:_s\d+)?(?:_\d+ns)?\.json$")
    by_delta = defaultdict(list)
    meta = {}
    for f in sorted(glob.glob(os.path.join(d, prefix + "*.json"))):
        m = pat.search(os.path.basename(f))
        if not m:
            continue
        rec = json.load(open(f))
        counts = rec["counts"][STEROL]
        n = len(counts["upper"])
        cut = int((1.0 - last_fraction) * n)
        up = float(np.mean(counts["upper"][cut:]))
        lo = float(np.mean(counts["lower"][cut:]))
        built = rec.get("built_counts")
        if not built:
            continue
        pl = {s: sum(v for k, v in built[s].items() if k != STEROL)
              for s in ("upper", "lower")}
        delta = int(m.group(1))
        by_delta[delta].append(100.0 * up / (up + lo))
        meta[delta] = (pl["upper"], pl["lower"],
                       built["upper"][STEROL], built["lower"][STEROL])
    out = []
    for delta in sorted(by_delta):
        pu, pl_, cu, cl = meta[delta]
        share = np.array(by_delta[delta], float)
        out.append({"delta": delta, "pl_upper": pu, "pl_lower": pl_,
                    "imbalance": 100.0 * (pu - pl_) / (pu + pl_),
                    "share": float(share.mean()), "n": len(share),
                    "sd": float(share.std(ddof=1)) if len(share) > 1 else None,
                    "chol_built": (cu, cl)})
    return out


def fit(points, degree=2):
    x = np.array([p["imbalance"] for p in points])
    y = np.array([p["share"] for p in points])
    if len(x) < degree + 1:
        degree = max(1, len(x) - 1)
    c = np.polyfit(x, y, degree)
    resid = y - np.polyval(c, x)
    return c, float(np.sqrt((resid ** 2).mean())), float(np.abs(resid).max())


def invert(c, target, lo, hi):
    """The imbalance inside [lo, hi] at which the curve gives `target`."""
    roots = np.roots(np.concatenate([c[:-1], [c[-1] - target]]))
    real = [float(r.real) for r in roots if abs(r.imag) < 1e-9]
    inside = [r for r in real if lo - 1e-9 <= r <= hi + 1e-9]
    if not inside:
        return None, real
    # the curve is monotone over the measured range, so one root is inside
    return min(inside, key=lambda r: abs(r)), real


def buildable(f, total):
    """The phospholipid pair nearest the imbalance `f`.

    The two leaflets share `total` molecules, so their difference has the parity
    of `total` and the imbalance moves in steps of 200 / total percentage points.
    """
    d_real = f / 100.0 * total
    best = None
    for d in range(int(np.floor(d_real)) - 2, int(np.ceil(d_real)) + 3):
        if (d - total) % 2:
            continue
        up, down = (total + d) // 2, (total - d) // 2
        if up < 0 or down < 0:
            continue
        got = 100.0 * (up - down) / total
        if best is None or abs(got - f) < abs(best[2] - f):
            best = (up, down, got)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", help="the directory run_curve.sh wrote its points to")
    ap.add_argument("--target", type=float, required=True,
                    help="the percentage of the sterol the upper leaflet is "
                         "required to hold after the run has settled")
    ap.add_argument("--total-pl", type=int, default=0,
                    help="the number of phospholipid molecules the two leaflets "
                         "will hold together in the system to be built "
                         "(default: the number the curve was measured on)")
    ap.add_argument("--total-chol", type=int, default=0,
                    help="the number of sterol molecules the system to be built "
                         "will hold (default: the number the curve was measured "
                         "on)")
    ap.add_argument("--prefix", default="relax_d")
    ap.add_argument("--last-fraction", type=float, default=0.4)
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--emit-env", action="store_true",
                    help="print the four numbers as shell variables")
    a = ap.parse_args()

    pts = read_points(a.dir, a.prefix, a.last_fraction)
    if len(pts) < 2:
        raise SystemExit("solve_target: fewer than two points under %s; "
                         "run run_curve.sh first" % a.dir)

    c, rms, worst = fit(pts, a.degree)
    xs = [p["imbalance"] for p in pts]
    lo, hi = min(xs), max(xs)
    ys = [p["share"] for p in pts]
    scatter = [p["sd"] for p in pts if p["sd"] is not None]

    print("solve_target: %d points, imbalance %.2f%% to %.2f%%, "
          "sterol share %.2f%% to %.2f%%" % (len(pts), lo, hi, min(ys), max(ys)))
    print("solve_target: fit residual %.2f points rms, %.2f worst" % (rms, worst))
    if scatter:
        print("solve_target: scatter between runs of one build, %s points"
              % ", ".join("%.2f" % s for s in scatter))

    f, real = invert(c, a.target, lo, hi)
    if f is None:
        print("solve_target: the curve does not reach %.2f%%. It covers %.2f%% "
              "to %.2f%%, and a target outside that range needs points measured "
              "there. Extrapolating the fit returns a number with no "
              "measurement behind it, and a leaflet stretched far enough stops "
              "behaving like a bilayer." % (a.target, min(ys), max(ys)))
        print("solve_target: a target the curve does not reach leaves this "
              "procedure. run_settled.sh builds a system, runs it, and builds "
              "it again at the sterol numbers that run settled to, which gives "
              "a composition the run keeps but not the one that was asked for.")
        sys.exit(2)

    # The points need not all hold the same number of phospholipid molecules:
    # the packing rounds a species present in a handful of copies up at one
    # imbalance and down at another. The default is the total of the measured
    # point nearest the solved imbalance, and the totals are named when they
    # differ so that a reader can pass --total-pl instead.
    totals = sorted({p["pl_upper"] + p["pl_lower"] for p in pts})
    near = min(pts, key=lambda p: abs(p["imbalance"] - f))
    total = a.total_pl or (near["pl_upper"] + near["pl_lower"])
    if not a.total_pl and len(totals) > 1:
        print("solve_target: the measured points hold %s phospholipid "
              "molecules; %d is used, from the point at %+.2f%%. Pass "
              "--total-pl for the system to be built."
              % (" and ".join(str(t) for t in totals), total,
                 near["imbalance"]))
    up, down, f_built = buildable(f, total)

    cu, cl = pts[0]["chol_built"]
    total_chol = a.total_chol or (cu + cl)
    # The phospholipid numbers set where the run takes the sterol. The sterol
    # numbers set where the build starts. Both are handed over, so that the
    # system starts at the distribution its phospholipid numbers demand.
    share_built = float(np.polyval(c, f_built))
    su = int(round(share_built / 100.0 * total_chol))
    sl = total_chol - su

    print("")
    print("  target sterol share      %.2f%% in the upper leaflet" % a.target)
    print("  imbalance that holds it  %+.2f%% of the phospholipid molecules" % f)
    print("  phospholipids            %d upper, %d lower (of %d), imbalance "
          "%+.2f%%" % (up, down, total, f_built))
    print("  sterol                   %d upper, %d lower (of %d), which is "
          "%.2f%%" % (su, sl, total_chol, 100.0 * su / total_chol))
    print("  the curve at that imbalance gives %.2f%%, and the build holds "
          "%.2f%%" % (share_built, 100.0 * su / total_chol))
    print("  the two leaflets share %d phospholipid molecules, so the imbalance "
          "moves in steps of %.2f points" % (total, 200.0 / total))
    if scatter:
        # how far the target moves if a point is off by the seed scatter
        slope = np.polyval(np.polyder(c), f)
        if abs(slope) > 1e-6:
            df = float(np.mean(scatter)) / abs(slope)
            print("  the scatter between runs of one build moves that "
                  "imbalance by %.2f%%, which is %.0f phospholipid molecules"
                  % (df, round(df / 100.0 * total)))
    print("  the curve was measured with the sterol at %d upper and %d lower"
          % (cu, cl))
    print("")
    print("  Build at all four numbers. The two phospholipid numbers are what "
          "the run keeps,")
    print("  and they are what carries the sterol to the target. The two sterol "
          "numbers are")
    print("  the target written into the build, and they are what makes the "
          "system start")
    print("  where the run would otherwise take it. A build given the "
          "phospholipid numbers")
    print("  alone reaches the target only after the sterol has moved, and over "
          "that part of")
    print("  the run it carries a composition that a Methods section does not "
          "state.")

    if a.emit_env:
        print("")
        print("PL_UPPER=%d" % up)
        print("PL_LOWER=%d" % down)
        print("CHOL_UPPER=%d" % su)
        print("CHOL_LOWER=%d" % sl)
        print("IMBALANCE_PCT=%.4f" % f_built)


if __name__ == "__main__":
    main()
