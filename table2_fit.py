#!/usr/bin/env python3
"""What the leaflet area measurement returns against what the build imposed.

Table 2 of the article gives four builds of one composition that differ in the
pair of mean areas per lipid handed to the packing and in nothing else. For each
build the table gives:

  the imposed difference   the difference between the two mean areas per lipid
                           handed to the packing, over their mean, in per cent.
                           This is an input to the build.
  the measured differences the difference between the two leaflets in the mean
                           area per lipid, and the difference between the two
                           leaflets in the area that one lipid species occupies,
                           each over the mean of the two, measured by Monte
                           Carlo integration of the planar Voronoi regions of
                           the packed system by leaflet_area_check.py. These are
                           outputs of the build.

The packing places the molecules from a random start, so each of the four builds
was repeated, four or five times over. This script regresses each measured
difference on the imposed difference and reports the slope, the offset and the
standard error of each, together with the spread over the repeats of one build.

A slope of one says that the measurement returns its own input. An offset that
is not separated from zero by its own standard error says the measurement adds
no constant to that input. A spread over the repeats that is large beside the
imposed difference says the measurement adds noise instead.

  python3 table2_fit.py

The numbers below are the signed difference between the two leaflets, in per
cent, read from the leaflets block of leaflet_area.json of every build on the
workstation:

  counts_run/t2/r{1,2,3,4}{a,b,c}   the three repeats built for this table
  counts_run/eqn, measured, table   the original builds of rows 1 to 3, which
                                    are the three arms of compare_counts.sh
  counts_run/table_rebuild          a fourth repeat of row 3
  counts_run/chol7, chol7_b         two repeats of row 4

The original build of row 4 was bench/fixpt, which was not retained. Its values
(cholesterol 0.586 upper against 0.664 lower, DLPC 0.680 against 0.816, leaflet
means 141.1 nm^2 over 218 against 139.7 over 188) fall inside the range that the
five repeats cover.
"""
import numpy as np

NAME = ["equal areas per lipid (EqN)",
        "areas measured on the packed system (none of the four routes)",
        "areas from a table of reported values (SA)",
        "seven more cholesterol molecules in the upper leaflet"]

APL_UPPER = np.array([0.700, 0.6891, 0.6760, 0.65429])
APL_LOWER = np.array([0.700, 0.7110, 0.7240, 0.75096])

# the difference between the two leaflet mean areas per lipid, per build
LEAF = [[-0.06, 0.18, 0.03, -0.12],
        [3.27, 2.57, 2.80, 2.97],
        [6.73, 7.07, 6.91, 6.64, 6.69],
        [13.81, 13.68, 13.78, 13.69, 13.87]]
# the difference between the two leaflets for cholesterol, per build
CHOL = [[2.56, -1.48, -2.59, -4.73],
        [0.37, -1.95, 3.84, 2.19],
        [9.03, 6.34, 7.98, 2.88, 3.77],
        [15.00, 10.39, 7.91, 6.14, 12.63]]
# the difference between the two leaflets for DLPC, per build
DLPC = [[0.34, 0.04, -1.28, 2.74],
        [3.27, 2.60, -1.35, 6.01],
        [5.12, 4.87, 0.84, 10.29, 5.40],
        [10.68, 9.72, 15.59, 10.84, 13.55]]

# the tolerance in force when these builds were made
TOL = 8.0

# The mean area per lipid handed to the packing, held at this value in all four
# builds, so that the four differ in the difference of the pair and in nothing
# else. COBY_APL of memble.sh.
COBY_APL = 0.70

IMPOSED = 100.0 * (APL_LOWER - APL_UPPER) / COBY_APL


def fit(y):
    """Regress the four means on the four imposed differences."""
    y = np.asarray(y, float)
    b, a = np.polyfit(IMPOSED, y, 1)
    res = y - (a + b * IMPOSED)
    s = float(np.sqrt(np.sum(res ** 2) / (len(y) - 2)))
    sxx = float(np.sum((IMPOSED - IMPOSED.mean()) ** 2))
    return b, s / np.sqrt(sxx), a, s * np.sqrt(1.0 / len(y) + IMPOSED.mean() ** 2 / sxx), s, res


def main():
    print("%-62s %7s %7s %9s" % ("build", "upper", "lower", "imposed"))
    for i in range(4):
        print("%-62s %7.4f %7.4f %8.2f%%"
              % (NAME[i], APL_UPPER[i], APL_LOWER[i], IMPOSED[i]))
    print("  mean area per lipid handed to the packing in all four builds  %.2f nm^2"
          % COBY_APL)
    print("")

    print("the difference between the two leaflets, per build, in per cent")
    print("%-10s %-12s %s" % ("build", "comparison", "repeats"))
    for i in range(4):
        for label, D in (("leaflet means", LEAF), ("cholesterol", CHOL), ("DLPC", DLPC)):
            v = np.array(D[i], float)
            print("row %d     %-14s %s   mean %+6.2f   sd %5.2f   n %d"
                  % (i + 1, label, " ".join("%+7.2f" % x for x in v),
                     v.mean(), v.std(ddof=1), len(v)))
        print("")

    print("each measured difference regressed on the imposed difference")
    for label, D in (("leaflet means", LEAF), ("cholesterol", CHOL), ("DLPC", DLPC)):
        y = [float(np.mean(v)) for v in D]
        b, sb, a, sa, s, res = fit(y)
        print("  %-14s slope %6.3f +- %.3f   offset %+6.2f +- %.2f   residual sd %.2f"
              % (label, b, sb, a, sa, s))
        print("  %-14s residuals %s" % ("", "  ".join("%+.2f" % r for r in res)))
    print("")

    print("what a single build reports, against what the repeats show")
    for i in range(4):
        sl = np.std(LEAF[i], ddof=1)
        ss = max(np.std(CHOL[i], ddof=1), np.std(DLPC[i], ddof=1))
        print("  row %d   leaflet means +- %.2f   species +- %.2f points over %d builds"
              % (i + 1, sl, ss, len(LEAF[i])))
    print("")

    print("the verdict at the %.0f%% tolerance, build by build" % TOL)
    for i in range(4):
        larger = [c if abs(c) > abs(d) else d for c, d in zip(CHOL[i], DLPC[i])]
        who = ["cholesterol" if abs(c) > abs(d) else "DLPC"
               for c, d in zip(CHOL[i], DLPC[i])]
        over = [abs(x) > TOL for x in larger]
        print("  row %d   %s" % (i + 1, "  ".join(
            "%+.2f %s %s" % (x, w[:4], "refused" if o else "passed")
            for x, w, o in zip(larger, who, over))))
        print("          passed in %d of %d, refused in %d"
              % (len(over) - sum(over), len(over), sum(over)))


if __name__ == "__main__":
    main()
