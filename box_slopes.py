#!/usr/bin/env python3
"""How fast the settled cholesterol share changes with the imbalance, in each box.

Section 3.7 of the article compares the response of a 32 by 32 by 16 nm bilayer
with the response of the 12 by 12 by 16 nm curve of Section 3.5. The 32 by 32 nm
set holds five imbalances, which is too few to fit a quadratic, so no coefficient
is quoted for that box. What can be quoted is the slope of the chord between two
measured points, and the fair comparison is the chord of the 12 by 12 nm curve
over the same two imbalances rather than the derivative of the quadratic at some
other imbalance.

  python3 box_slopes.py

The five 32 by 32 nm values are the mean over the last 40% of a 250 ns run, read
from the relax file of each build; the value at f = -11.93% is the first of three
runs that differ only in their starting velocities (74.04, 74.13, 73.36). The
imbalance of each is computed from the topology of that 32 by 32 nm system and
not from delta. The 12 by 12 nm values are the quadratic of curve_stats.py
evaluated at the same five imbalances.
"""
import numpy as np

F = np.array([-11.93, -6.04, -1.60, 2.74, 11.62])     # imbalance, per cent
BIG = np.array([74.04, 67.84, 61.62, 56.42, 46.66])   # 32 by 32 nm, upper share
QUAD = [-0.02829, -0.8182, 59.202]                    # curve_stats.py, 12 by 12 nm
SMALL = np.polyval(QUAD, F)


def main():
    print("%9s %10s %10s %11s" % ("f (%)", "32 x 32", "12 x 12", "difference"))
    for f, b, s in zip(F, BIG, SMALL):
        print("%+8.2f%% %9.2f%% %9.2f%% %+10.2f" % (f, b, s, b - s))
    print("")
    print("slope of the chord between adjacent points, percentage points per unit of f")
    print("%-26s %10s %10s %7s" % ("interval in f", "32 x 32", "12 x 12", "ratio"))
    for i in range(len(F) - 1):
        df = F[i + 1] - F[i]
        sb = (BIG[i + 1] - BIG[i]) / df
        ss = (SMALL[i + 1] - SMALL[i]) / df
        print("%+7.2f%% to %+7.2f%%      %9.2f %10.2f %7.1f"
              % (F[i], F[i + 1], sb, ss, sb / ss))
    print("")
    print("derivative of the 12 by 12 nm quadratic, for reference")
    for f in (0.0, -4.0, -8.0):
        print("  at f = %+5.1f%%   %.2f points per unit of f"
              % (f, QUAD[1] + 2.0 * QUAD[0] * f))
    print("")
    print("Five points carry no fit, so the 32 by 32 nm response is reported as")
    print("these four chords and not as a coefficient.")


if __name__ == "__main__":
    main()
