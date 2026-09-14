#!/usr/bin/env python3
"""The settled cholesterol share of the systems built in a large box, against the curve.

Section 3.5 measures a curve on 12 by 12 by 16 nm bilayers and then measures the
same quantity on 32 by 32 by 16 nm bilayers at five imbalances. This
script reads the second set and prints each value beside the value the curve of
the first set returns at the same imbalance.

  python3 bigread.py --root ~/                      # the five points and the repeats
  python3 bigread.py --root ~/ --prefix ext_d       # the points of the extended range

The value of a point is the mean over the last 40% of the run of the number of
cholesterol molecules in the upper leaflet, divided by the number of cholesterol
molecules in the system, which is what leaflet_relax.py writes and what
curve_report.py reads. The curve is the one curve_stats.py fits.
"""
import argparse
import glob
import json
import os
import statistics as st

# The quadratic of Section 3.5, returned by curve_stats.py from the six points
# of the 12 by 12 by 16 nm series.
A2, A1, A0 = -0.028290, -0.818244, 59.202100


def curve(f):
    return A0 + A1 * f + A2 * f * f


# The phospholipid numbers each leaflet of each point received, read from the
# system.top of that build. The numbers are not a function of delta: the packing
# rounds the ratio of the minor species, and the two leaflets of the points at
# delta = -20, -12 and -6 hold 273 phospholipid molecules together where the
# other three hold 272. The one molecule is the phosphoinositide of the lower
# leaflet, which the packing rounds to eight there and to seven in the other
# three. The imbalance is therefore taken from these counts and never computed
# from delta.
COUNTS = {-20: (120, 153), -12: (128, 145), -6: (134, 139),
          0: (140, 132), 6: (146, 126), 12: (152, 120)}
STEROL = "CHOL"


def imbalance_from_counts(upper, lower):
    """The difference between the two phospholipid numbers over their sum."""
    return 100.0 * (upper - lower) / float(upper + lower)


def imbalance_of_point(path, delta):
    """The imbalance of a built point, from its own topology where possible.

    leaflet_relax.py writes the numbers the topology gave each leaflet into
    built_counts. Where that block is present the imbalance is taken from it, so
    that the value belongs to the system that was run. Where it is absent the
    counts above are used, and those were read from the same topologies.
    """
    d = json.load(open(path))
    b = d.get("built_counts")
    if b and "upper" in b and "lower" in b:
        up = sum(v for k, v in b["upper"].items() if k != STEROL)
        lo = sum(v for k, v in b["lower"].items() if k != STEROL)
        if up and lo:
            return imbalance_from_counts(up, lo), "built_counts of that run"
    up, lo = COUNTS[delta]
    return imbalance_from_counts(up, lo), "the counts table in this script"


def settled(path, last=0.4):
    d = json.load(open(path))
    up = d["counts"]["CHOL"]["upper"]
    lo = d["counts"]["CHOL"]["lower"]
    total = up[0] + lo[0]
    cut = int((1.0 - last) * len(up))
    return 100.0 * st.mean(up[cut:]) / total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.expanduser("~"))
    ap.add_argument("--prefix", default="bigbox_d")
    ap.add_argument("--last-fraction", type=float, default=0.4)
    a = ap.parse_args()

    print("%-24s %8s %9s %9s %8s   %s"
          % ("point", "f (%)", "measured", "curve", "diff", "where f comes from"))
    for d in (0, -6, -12, -20, 12):
        p = os.path.join(a.root, "%s%d" % (a.prefix, d), "curve",
                         "relax_d%d.json" % d)
        if not os.path.exists(p):
            continue
        f, where = imbalance_of_point(p, d)
        v = settled(p, a.last_fraction)
        print("%-24s %8.2f %9.2f %9.2f %+8.2f   %s"
              % ("delta %d" % d, f, v, curve(f), v - curve(f), where))

    reps = sorted(glob.glob(os.path.join(a.root, "%s-20" % a.prefix, "curve",
                                         "relax_d-20*.json")))
    if len(reps) > 1:
        print("")
        print("the runs of delta -20, one build and one velocity seed each")
        vals = []
        for p in reps:
            v = settled(p, a.last_fraction)
            vals.append(v)
            print("  %-38s %6.2f" % (os.path.basename(p), v))
        print("  mean %.2f, sample sd %.2f over %d runs"
              % (st.mean(vals), st.stdev(vals), len(vals)))
        f20 = imbalance_from_counts(*COUNTS[-20])
        print("  the curve at f = %.2f%% returns %.2f" % (f20, curve(f20)))


if __name__ == "__main__":
    main()
