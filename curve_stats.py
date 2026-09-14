#!/usr/bin/env python3
"""Every number the article quotes for the curve of Section 3.4, from the six points.

The six points are the settled cholesterol share of the six 12 by 12 by 16 nm
systems, against the phospholipid imbalance of each. Two of the six were run
three times from three sets of starting velocities, and the standard deviation
over those three runs is what this script uses as the uncertainty of a point.

What this script prints, and where the article uses each number:

  the quadratic fit and its residuals          Section 3.4, Figure 4
  the imbalance at which the curve crosses     Section 3.4, Section 3.5,
    the share the build assigned                 Conclusions, Figure 4 legend
  the uncertainty on that imbalance            the same four places
  the same crossing from a straight line       Section 3.4
  the F test that keeps the quadratic term     Section 3.4

Run with no arguments. The six points are written in the source below, and the
source of each of the six is named beside it.

  python3 curve_stats.py
"""
import numpy as np

# ---------------------------------------------------------------------------
# The six points. Each value is the mean over the last 40% of a 250 ns run of
# the number of cholesterol molecules in the upper leaflet, divided by the 134
# cholesterol molecules of the system. The file each value was read from is
# relax_d<delta>.json, written by leaflet_relax.py at the end of run_curve.sh.
# n is the number of runs of that build, and sd is the standard deviation over
# those runs where n is three.
# ---------------------------------------------------------------------------
DELTA = np.array([-20, -12, -6, 0, 6, 12])
# The phospholipid numbers each leaflet of each point received, read from the
# system.top of that build. The imbalance below is computed from these two
# numbers and is not a function of delta: the packing rounds the ratio of the
# minor species, and the points at delta = -20, -12 and -6 hold 273
# phospholipid molecules over the two leaflets where the other three hold 272.
# The three that hold 273 carry eight phosphoinositide molecules in the lower
# leaflet against seven in the other three; every other count is the same.
PL_UPPER = np.array([120, 128, 134, 140, 146, 152])
PL_LOWER = np.array([153, 145, 139, 132, 126, 120])
F = 100.0 * (PL_UPPER - PL_LOWER) / (PL_UPPER + PL_LOWER)   # imbalance, %
C = np.array([64.78, 63.20, 61.17, 56.86, 50.26, 46.36])  # upper share, %
N = np.array([3, 1, 1, 3, 1, 1])
SD = np.array([2.69, np.nan, np.nan, 1.71, np.nan, np.nan])

BUILT = 100.0 * 71.0 / 134.0     # the share the build gave every one of the six

# The uncertainty of a point. Where a build was run three times, the standard
# error of the mean is the standard deviation over the three runs divided by the
# square root of three. Where a build was run once, the uncertainty is taken as
# the pooled standard deviation of the two builds that were run three times,
# because one run carries no scatter of its own.
POOLED_SD = float(np.sqrt(np.nanmean(SD[np.isfinite(SD)] ** 2)))
SIGMA = np.where(N == 3, SD / np.sqrt(3.0), POOLED_SD)


def crossing(coef, level):
    """The imbalance at which the fitted curve equals `level`.

    A quadratic that opens downward meets a horizontal line twice, and the two
    meetings sit on either side of the maximum. The root wanted here is the one
    the measured points bracket, which is the real root of smallest absolute
    value, because every one of the six points lies within 12.09% of zero. A
    draw whose fit never reaches `level` returns nan and is counted separately.
    """
    c = np.array(coef, dtype=float)
    c[-1] -= level
    r = np.roots(c)
    r = r[np.abs(r.imag) < 1e-9].real
    if r.size == 0:
        return float("nan")
    return float(r[np.argmin(np.abs(r))])


def main():
    print("the six points")
    print("  %6s %6s %6s %8s %8s %4s %8s"
          % ("delta", "PL up", "PL lo", "f (%)", "C (%)", "n", "sigma"))
    for d, pu, pl, f, c, n, sg in zip(DELTA, PL_UPPER, PL_LOWER, F, C, N, SIGMA):
        print("  %6d %6d %6d %8.2f %8.2f %4d %8.2f" % (d, pu, pl, f, c, n, sg))
    print("  pooled standard deviation of the two repeated builds  %.2f" % POOLED_SD)
    print("  share the build assigned to the upper leaflet         %.2f%%" % BUILT)
    print("")

    q = np.polyfit(F, C, 2)
    l = np.polyfit(F, C, 1)
    rq = C - np.polyval(q, F)
    rl = C - np.polyval(l, F)
    print("quadratic fit   C = %.3f %+.4f f %+.5f f^2" % (q[2], q[1], q[0]))
    print("  rms residual %.2f pt, largest residual %.2f pt"
          % (np.sqrt(np.mean(rq ** 2)), np.max(np.abs(rq))))
    print("straight line   C = %.3f %+.4f f" % (l[1], l[0]))
    print("  rms residual %.2f pt, largest residual %.2f pt"
          % (np.sqrt(np.mean(rl ** 2)), np.max(np.abs(rl))))
    print("")

    print("value at zero imbalance   quadratic %.2f%%   line %.2f%%"
          % (np.polyval(q, 0.0), np.polyval(l, 0.0)))
    print("")

    fq = crossing(q, BUILT)
    fl = crossing(l, BUILT)
    print("imbalance at which the curve returns the share the build assigned")
    print("  from the quadratic  %+.2f%%" % fq)
    print("  from the line       %+.2f%%" % fl)

    # The uncertainty. Each point is drawn from a normal distribution centred on
    # the measured value with the standard deviation of the column sigma above,
    # the fit is repeated, and the crossing is found again.
    rng = np.random.default_rng(0)
    NDRAW = 20000
    sq, sl = [], []
    for _ in range(NDRAW):
        y = C + rng.normal(0.0, SIGMA)
        sq.append(crossing(np.polyfit(F, y, 2), BUILT))
        sl.append(crossing(np.polyfit(F, y, 1), BUILT))
    sq = np.array(sq); sl = np.array(sl)
    lost_q = int(np.sum(~np.isfinite(sq))); lost_l = int(np.sum(~np.isfinite(sl)))
    sq = sq[np.isfinite(sq)]; sl = sl[np.isfinite(sl)]
    print("  quadratic  %+.2f +- %.2f over %d draws, %d draws had no real crossing"
          % (sq.mean(), sq.std(), len(sq), lost_q))
    print("  line       %+.2f +- %.2f over %d draws, %d draws had no real crossing"
          % (sl.mean(), sl.std(), len(sl), lost_l))
    print("  draws whose crossing falls inside the measured range of %+.2f to %+.2f%%: %d of %d"
          % (F.min(), F.max(), int(np.sum((sq >= F.min()) & (sq <= F.max()))), len(sq)))
    inside = sq[(sq >= F.min()) & (sq <= F.max())]
    print("  over those draws alone, quadratic %+.2f +- %.2f" % (inside.mean(), inside.std()))
    print("")

    # Does the quadratic term earn its place? The residual sum of squares of the
    # line and of the quadratic are compared by an F test with one and three
    # degrees of freedom, six points carrying three parameters in the quadratic.
    rss_l = float(np.sum(rl ** 2))
    rss_q = float(np.sum(rq ** 2))
    Fstat = ((rss_l - rss_q) / 1.0) / (rss_q / (len(F) - 3))
    print("F test on the quadratic term")
    print("  residual sum of squares, line      %.3f" % rss_l)
    print("  residual sum of squares, quadratic %.3f" % rss_q)
    print("  F(1,3) = %.1f, against 10.13 at the 5%% level" % Fstat)
    print("")

    print("value the curve returns at the three imbalances of Table 4")
    for f in (-0.37, 1.10, 2.94):
        c = np.polyval(q, f)
        print("  f = %+5.2f%%   share %.2f%%   molecules moved %+.1f"
              % (f, c, (c - BUILT) * 134.0 / 100.0))


if __name__ == "__main__":
    main()
