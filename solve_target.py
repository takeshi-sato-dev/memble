#!/usr/bin/env python3
"""solve_target.py : the leaflet numbers that hold a required cholesterol content.

A build sets how many cholesterol molecules each leaflet receives, and the run
does not keep those numbers. What the run keeps is set by how many phospholipid
molecules each leaflet holds. run_curve.sh measures that relation over a range
of phospholipid numbers, and this script reads the relation backwards: given the
cholesterol content a leaflet is required to hold once the run has settled, it
returns the two phospholipid numbers and the two cholesterol numbers a build has
to be given.

The target is stated the way a build states a composition, as the cholesterol
mol% of a leaflet:

    cholesterol mol% of a leaflet = 100 * CHOL / (CHOL + phospholipid)

counted over the molecules of that leaflet alone. That is the number a user of
any builder types, and it is the number a Methods section prints. The fit itself
runs on a different quantity, the share of the cholesterol of the whole system
that one leaflet holds, because that share is what the six measured systems have
in common. The script converts between the two and reports both, so that nobody
has to.

The build is given four numbers and not two. The two phospholipid numbers are
the quantity the run cannot change, and they are what carries the cholesterol to
the target. The two cholesterol numbers are the target written into the build,
and they are what makes the system start where the run would otherwise take it.
A build given the phospholipid numbers alone reaches the target only after the
cholesterol has moved, and over that part of the run it carries a composition
that a Methods section does not state.

Usage:
    python3 solve_target.py <curve dir> --chol-upper 40 --total-pl 272 --total-chol 134
    python3 solve_target.py <curve dir> --chol-lower 25 --total-pl 272 --total-chol 134
    python3 solve_target.py <curve dir> --chol-upper 40 --total-pl 272 --total-chol 134 --emit-env
    python3 solve_target.py <curve dir> --band --total-pl 272 --total-chol 134

    python3 solve_target.py <curve dir> --share 56.55        # the internal variable

<curve dir> holds the relax_d<delta>.json files run_curve.sh writes. Repeats of
a point (relax_d0_s2_250ns.json) are averaged into that point.

The script refuses a target the measured curve does not reach, and says in mol%
what it does reach. Extrapolating a fitted quadratic past its data returns a
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


# ---------------------------------------------------------------- mol%
#
# The two leaflets of a system holding `total_pl` phospholipid molecules and
# `total_chol` cholesterol molecules, at a phospholipid imbalance f, hold
#
#     phospholipid   upper  total_pl * (1 + f/100) / 2
#                    lower  total_pl * (1 - f/100) / 2
#     cholesterol    upper  total_chol * C(f) / 100
#                    lower  total_chol - that
#
# with C the fitted curve. The cholesterol mol% of a leaflet is the cholesterol
# of that leaflet over every molecule of that leaflet. Both are monotone in f
# over the range the curve was fitted on, the upper falling as f rises and the
# lower rising, so each is inverted by bisection.


def mol_percent(c, f, total_pl, total_chol):
    """(upper, lower) cholesterol mol% at imbalance f."""
    share = float(np.polyval(c, f))
    pu = total_pl * (1.0 + f / 100.0) / 2.0
    pl_ = total_pl - pu
    su = total_chol * share / 100.0
    sl = total_chol - su
    return 100.0 * su / (pu + su), 100.0 * sl / (pl_ + sl)


def invert_mol(c, target, leaflet, lo, hi, total_pl, total_chol):
    """The imbalance inside [lo, hi] at which `leaflet` holds `target` mol%.

    Returns (f, None) on success and (None, (lo_value, hi_value)) when the
    target lies outside what the measured points reach, so that the caller can
    name the range in the same unit the user asked in.
    """
    idx = 0 if leaflet == "upper" else 1
    a = mol_percent(c, lo, total_pl, total_chol)[idx]
    b = mol_percent(c, hi, total_pl, total_chol)[idx]
    if not (min(a, b) - 1e-9 <= target <= max(a, b) + 1e-9):
        return None, (min(a, b), max(a, b))
    x0, x1 = lo, hi
    for _ in range(200):
        mid = 0.5 * (x0 + x1)
        v = mol_percent(c, mid, total_pl, total_chol)[idx]
        if (v - target) * (mol_percent(c, x0, total_pl, total_chol)[idx]
                           - target) <= 0.0:
            x1 = mid
        else:
            x0 = mid
    return 0.5 * (x0 + x1), None


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
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--chol-upper", type=float,
                   help="the cholesterol mol%% the upper leaflet is required to "
                        "hold once the run has settled, counted over the "
                        "molecules of that leaflet alone")
    g.add_argument("--chol-lower", type=float,
                   help="the same for the lower leaflet")
    g.add_argument("--share", "--target", type=float, dest="share",
                   help="the target as the percentage of the cholesterol of the "
                        "whole system that the upper leaflet holds, which is the "
                        "quantity the curve is fitted on")
    g.add_argument("--band", action="store_true",
                   help="print the range of targets the measured points reach "
                        "and exit")
    ap.add_argument("--total-pl", type=int, default=0,
                    help="the number of phospholipid molecules the two leaflets "
                         "will hold together in the system to be built "
                         "(default: the number the curve was measured on)")
    ap.add_argument("--total-chol", type=int, default=0,
                    help="the number of cholesterol molecules the system to be "
                         "built will hold (default: the number the curve was "
                         "measured on)")
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

    # The points need not all hold the same number of phospholipid molecules:
    # the packing rounds a species present in a handful of copies up at one
    # imbalance and down at another. The default is the total of the measured
    # point nearest the middle of the range, and the totals are named when they
    # differ so that a reader can pass --total-pl instead.
    totals = sorted({p["pl_upper"] + p["pl_lower"] for p in pts})
    mid = min(pts, key=lambda p: abs(p["imbalance"] - 0.5 * (lo + hi)))
    total = a.total_pl or (mid["pl_upper"] + mid["pl_lower"])
    cu0, cl0 = pts[0]["chol_built"]
    total_chol = a.total_chol or (cu0 + cl0)

    up_lo, low_lo = mol_percent(c, lo, total, total_chol)
    up_hi, low_hi = mol_percent(c, hi, total, total_chol)

    print("solve_target: %d points, imbalance %.2f%% to %.2f%%" % (len(pts), lo, hi))
    print("solve_target: with %d phospholipid and %d cholesterol molecules, the "
          "measured points reach" % (total, total_chol))
    print("              cholesterol %.1f to %.1f mol%% in the upper leaflet, "
          "%.1f to %.1f mol%% in the lower"
          % (min(up_lo, up_hi), max(up_lo, up_hi),
             min(low_lo, low_hi), max(low_lo, low_hi)))
    print("solve_target: fit residual %.2f points rms, %.2f worst" % (rms, worst))
    if scatter:
        print("solve_target: scatter between runs of one build, %s points of "
              "the share" % ", ".join("%.2f" % s for s in scatter))
    if not a.total_pl and len(totals) > 1:
        print("solve_target: the measured points hold %s phospholipid "
              "molecules; %d is used. Pass --total-pl for the system to be built."
              % (" and ".join(str(t) for t in totals), total))

    if a.band:
        return

    # ------------------------------------------------------------ the target
    if a.share is not None:
        f, real = invert(c, a.share, lo, hi)
        asked = "%.2f%% of the cholesterol of the system, in the upper leaflet" % a.share
        reach = "%.2f%% to %.2f%%" % (min(ys), max(ys))
    else:
        leaflet = "upper" if a.chol_upper is not None else "lower"
        want = a.chol_upper if leaflet == "upper" else a.chol_lower
        f, bad = invert_mol(c, want, leaflet, lo, hi, total, total_chol)
        asked = "cholesterol %.2f mol%% in the %s leaflet" % (want, leaflet)
        reach = ("%.1f mol%% to %.1f mol%%" % bad) if bad else ""

    if f is None:
        print("")
        print("solve_target: the curve does not reach %s. It covers %s, and a "
              "target outside that range needs points measured there."
              % (asked, reach))
        print("solve_target: raising the cholesterol content of BOTH leaflets is "
              "a different request, and no curve is needed for it: build the "
              "system with more cholesterol. The curve divides the cholesterol "
              "a system holds between its two leaflets, and it has to be "
              "measured again at the new content.")
        print("solve_target: a target the curve does not reach leaves this "
              "procedure. run_settled.sh builds a system, runs it, and builds "
              "it again at the cholesterol numbers that run settled to, which "
              "gives a composition the run keeps but not the one that was "
              "asked for.")
        sys.exit(2)

    up, down, f_built = buildable(f, total)
    share_built = float(np.polyval(c, f_built))
    su = int(round(share_built / 100.0 * total_chol))
    sl = total_chol - su
    got_up = 100.0 * su / (up + su)
    got_lo = 100.0 * sl / (down + sl)

    print("")
    print("  asked for                %s" % asked)
    print("  phospholipids            %d upper, %d lower (of %d)" % (up, down, total))
    print("  cholesterol              %d upper, %d lower (of %d)" % (su, sl, total_chol))
    print("")
    print("  the build then holds     %.2f mol%% cholesterol in the upper leaflet"
          % got_up)
    print("                           %.2f mol%% cholesterol in the lower leaflet"
          % got_lo)
    print("  phospholipid imbalance   %+.2f%% (the variable of the fit)" % f_built)
    print("  share of the cholesterol %.2f%% in the upper leaflet" % share_built)
    print("")
    print("  the two leaflets share %d phospholipid molecules, so the imbalance "
          "moves in steps of %.2f points" % (total, 200.0 / total))
    if scatter:
        # how far the target moves if a point is off by the seed scatter
        slope = np.polyval(np.polyder(c), f)
        if abs(slope) > 1e-6:
            df = float(np.mean(scatter)) / abs(slope)
            m_a = mol_percent(c, f, total, total_chol)[0]
            m_b = mol_percent(c, f + df, total, total_chol)[0]
            print("  the scatter between runs of one build is worth %.2f mol%% "
                  "of cholesterol, which is the precision of a target"
                  % abs(m_b - m_a))
    print("  the curve was measured with the cholesterol at %d upper and %d lower"
          % (cu0, cl0))
    print("")
    print("  Build at all four numbers. The two phospholipid numbers are what "
          "the run keeps,")
    print("  and they are what carries the cholesterol to the target. The two "
          "cholesterol numbers")
    print("  are the target written into the build, and they are what makes the "
          "system start")
    print("  where the run would otherwise take it. A build given the "
          "phospholipid numbers")
    print("  alone reaches the target only after the cholesterol has moved, and "
          "over that part")
    print("  of the run it carries a composition that a Methods section does not "
          "state.")

    if a.emit_env:
        print("")
        print("PL_UPPER=%d" % up)
        print("PL_LOWER=%d" % down)
        print("CHOL_UPPER=%d" % su)
        print("CHOL_LOWER=%d" % sl)
        print("CHOL_MOLPCT_UPPER=%.2f" % got_up)
        print("CHOL_MOLPCT_LOWER=%.2f" % got_lo)
        print("IMBALANCE_PCT=%.4f" % f_built)


if __name__ == "__main__":
    main()
