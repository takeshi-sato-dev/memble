#!/usr/bin/env python3
"""The line of Section 3.5, fitted on the imbalance the production run starts from.

Every system is entered here as the pair it was measured as: the phospholipid
imbalance counted on the structure that begins its production run, from
`count_leaflets_gro.py`, and the share of the cholesterol of the system that the
upper leaflet holds over the settled part of that run, from `settle2.py` or from
the tables of Section 3.5 and Section S20.

Three groups of systems are entered and the script reports each against the
others, because they are three independent tests of the same relation:

  curve     the six systems of Section 3.5, which differ only in how the DLPC
            molecules are divided, 250 ns each
  ext       the nine systems of Section S20, built at a larger imbalance than
            the six and run 500 ns each. The three systems of one build settle
            at different imbalances, so each is entered separately
  second    the five systems of the second series, in which the sphingomyelin
            of the upper leaflet is moved with the DLPC, 250 ns each

Usage:
    python3 fit_curve_f.py
    python3 fit_curve_f.py --sd 2.83 --draws 20000
"""
import argparse

import numpy as np

CURVE = [  # f at the start of the run, settled share, label
    (-3.30, 64.78, "d-20"), (-2.56, 63.20, "d-12"), (-1.83, 61.17, "d-6"),
    (2.94, 56.86, "d0"), (7.35, 50.26, "d6"), (8.09, 46.36, "d12"),
]
EXT = [
    (-9.89, 74.42, "ext_d-26"), (-9.16, 75.13, "ext_d-26_b"),
    (-3.30, 65.91, "ext_d-26_c"), (-16.48, 83.07, "ext_d-32"),
    (-15.75, 81.20, "ext_d-32_b"), (-11.36, 77.71, "ext_d-32_c"),
    (-13.55, 79.55, "ext_d-40"), (-13.55, 80.08, "ext_d-40_b"),
    (-5.49, 71.82, "ext_d-40_c"),
]
SECOND = [
    (2.56, 56.77, "a-32"), (1.10, 60.52, "a-16"), (-1.10, 62.45, "a-4"),
    (5.15, 54.56, "a20"), (5.15, 53.41, "a32"),
]
BUILT_SHARE = 52.99


def fit(points, deg=1):
    f = np.array([p[0] for p in points])
    c = np.array([p[1] for p in points])
    k = np.polyfit(f, c, deg)
    r = c - np.polyval(k, f)
    return k, r, f, c


def line(name, points, sd, draws, rng):
    k, r, f, c = fit(points, 1)
    kq, rq, _, _ = fit(points, 2)
    F = ((r ** 2).sum() - (rq ** 2).sum()) / ((rq ** 2).sum() / (len(f) - 3))
    print("%s, %d systems, f from %+.2f to %+.2f" % (name, len(f), f.min(), f.max()))
    print("  C = %.2f %+.3f f     rms %.2f   largest %.2f"
          % (k[1], k[0], np.sqrt((r ** 2).mean()), np.abs(r).max()))
    print("  a quadratic term gives rms %.2f and F(1,%d) = %.2f"
          % (np.sqrt((rq ** 2).mean()), len(f) - 3, F))
    cr = []
    for _ in range(draws):
        kb = np.polyfit(f, c + rng.normal(0, sd, len(c)), 1)
        cr.append((kb[1] - BUILT_SHARE) / -kb[0])
    print("  value at f = 0   %.2f%%" % k[1])
    print("  crossing of %.2f%%   %+.2f +- %.2f%%"
          % (BUILT_SHARE, (k[1] - BUILT_SHARE) / -k[0], np.array(cr).std()))
    print("  r2 = %.3f" % np.corrcoef(f, c)[0, 1] ** 2)
    print("")
    return k


def against(name, points, k):
    f = np.array([p[0] for p in points])
    c = np.array([p[1] for p in points])
    r = c - np.polyval(k, f)
    print("%s against that line" % name)
    for i, p in enumerate(points):
        print("  %-12s f %+7.2f   measured %6.2f   line %6.2f   %+6.2f"
              % (p[2], p[0], p[1], np.polyval(k, p[0]), r[i]))
    print("  mean %+.2f   rms %.2f   largest %.2f"
          % (r.mean(), np.sqrt((r ** 2).mean()), np.abs(r).max()))
    print("")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sd", type=float, default=2.83,
                    help="the scatter between runs of one build, in points")
    ap.add_argument("--draws", type=int, default=20000)
    a = ap.parse_args()
    rng = np.random.default_rng(5)

    k6 = line("the six of Section 3.5", CURVE, a.sd, a.draws, rng)
    against("the nine of Section S20", EXT, k6)
    against("the second series", SECOND, k6)
    k15 = line("the six and the nine", CURVE + EXT, a.sd, a.draws, rng)
    against("the second series", SECOND, k15)
    line("all twenty", CURVE + EXT + SECOND, a.sd, a.draws, rng)


if __name__ == "__main__":
    main()
