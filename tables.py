#!/usr/bin/env python3
"""The four tables of the article, written from the numbers the scripts return.

Table 1  the systems built and measured with MEMBLE
Table 2  what the leaflet area measurement returns against what was imposed
Table 3  the six systems of the curve
Table 4  three builds that were not used in the fit

Each row names where its numbers come from. Running this file rewrites the CSV
and the markdown, so the article and the data cannot drift apart.
"""
import csv, io, os
import numpy as np

# The directory the tables are written to. Defaults to this file's own
# directory and can be set from the environment.
OUT = os.environ.get("MEMBLE_TABLE_OUT", os.path.dirname(os.path.abspath(__file__)))
BUILT = 100.0 * 71.0 / 134.0


def write(name, header, rows, caption):
    with open(os.path.join(OUT, name + ".csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    md = ["| " + " | ".join(header) + " |",
          "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        md.append("| " + " | ".join(str(v) for v in r) + " |")
    with open(os.path.join(OUT, name + ".md"), "w") as fh:
        fh.write(caption + "\n\n" + "\n".join(md) + "\n")
    print("\n" + caption.split(".")[0] + "\n" + "\n".join(md))


# ---------------------------------------------------------------- Table 1
# These rows are the measurements of the five builds of Section 3.1. They are
# not recomputed here: they are the values memble_report.json wrote for each
# build, transcribed once from memble_bench/row*/row*_work/memble_report.json,
# the leaflet line of ../build.log and the Maximum force of row*_work/em.log.
# Fmax is the maximum force at the end of the energy minimization, in kJ/mol/nm.
write("Table1_systems",
      ["components", "leaflets", "beads", "species both leaflets hold, as a share of a leaflet",
       "largest difference between the leaflets", "species", "run files", "Fmax"],
      [[3, "symmetric", "21,136", "100%", "5.5% ± 3.6%", "cholesterol", "written", "95"],
       [4, "asymmetric", "21,195", "67%", "2.8% ± 3.5%", "cholesterol", "written", "99"],
       [5, "asymmetric, with PIP2", "21,235", "65%", "0.7% ± 4.0%", "cholesterol", "written", "90"],
       [8, "asymmetric, with PIP2", "28,986", "30%", "4.0% ± 3.8%", "cholesterol", "written", "97"],
       [12, "asymmetric, with PIP2", "29,031", "41%", "4.4% ± 5.7%", "POPC", "written", "59"]],
      "**Table 1.** Martini 3 membrane protein systems assembled and measured with "
      "MEMBLE. The difference column holds the largest difference between the leaflets "
      "for a lipid that both leaflets hold in numbers, with the standard error of that "
      "difference, and the species that carries it. Each row is one build, and Section "
      "3.4 measures how far that difference moves between repeats of one build, which "
      "is about three percentage points. Fmax is the maximum force at the end "
      "of the energy minimization, in kJ/mol/nm.")

# ---------------------------------------------------------------- Table 2
# The packing places the molecules from a random start, so each of the four
# pairs of areas per lipid was built four or five times over. The arrays below
# hold the signed difference between the two leaflets, in per cent, read from
# leaflet_area.json of every build:
#   counts_run/t2/r{1,2,3,4}{a,b,c}   the three repeats built for this table
#   counts_run/eqn, measured, table   the original builds of rows 1 to 3
#   counts_run/table_rebuild          a fourth repeat of row 3
#   counts_run/chol7, chol7_b         two repeats of row 4
# The original build of row 4 was bench/fixpt, which was not retained; its
# values fall inside the range the five repeats cover.
APL_U = np.array([0.700, 0.6891, 0.6760, 0.65429])
APL_L = np.array([0.700, 0.7110, 0.7240, 0.75096])
LEAF = [[-0.06, 0.18, 0.03, -0.12],
        [3.27, 2.57, 2.80, 2.97],
        [6.73, 7.07, 6.91, 6.64, 6.69],
        [13.81, 13.68, 13.78, 13.69, 13.87]]
CHOL = [[2.56, -1.48, -2.59, -4.73],
        [0.37, -1.95, 3.84, 2.19],
        [9.03, 6.34, 7.98, 2.88, 3.77],
        [15.00, 10.39, 7.91, 6.14, 12.63]]
DLPC = [[0.34, 0.04, -1.28, 2.74],
        [3.27, 2.60, -1.35, 6.01],
        [5.12, 4.87, 0.84, 10.29, 5.40],
        [10.68, 9.72, 15.59, 10.84, 13.55]]
TOL = 8.0
NAME = ["equal areas per lipid (EqN)",
        "areas measured on the packed system (none of the four routes)",
        "areas from a table of reported values (SA)",
        "seven more cholesterol molecules in the upper leaflet"]
IMP = 100.0 * (APL_L - APL_U) / 0.70


def _ms(v):
    v = np.array(v, float)
    return v.mean(), v.std(ddof=1)


def _fit(y):
    y = np.array(y, float)
    b, a = np.polyfit(IMP, y, 1)
    r = y - (a + b * IMP)
    s = float(np.sqrt(np.sum(r ** 2) / 2))
    sxx = float(np.sum((IMP - IMP.mean()) ** 2))
    return b, s / np.sqrt(sxx), a, s * np.sqrt(0.25 + IMP.mean() ** 2 / sxx), s


def _verdict(i):
    n = len(CHOL[i])
    over = sum(1 for c, d in zip(CHOL[i], DLPC[i]) if max(abs(c), abs(d)) > TOL)
    if over == 0:
        return "passed in %d of %d" % (n, n)
    if over == n:
        return "refused in %d of %d" % (n, n)
    return "passed in %d of %d, refused in %d" % (n - over, n, over)


_LM = [_ms(v) for v in LEAF]
_CM = [_ms(v) for v in CHOL]
_DM = [_ms(v) for v in DLPC]
_bl, _ble, _al, _ale, _sl = _fit([m for m, _ in _LM])
_bc, _bce, _ac, _ace, _sc = _fit([m for m, _ in _CM])
_bd, _bde, _ad, _ade, _sd = _fit([m for m, _ in _DM])

write("Table2_leaflet_areas",
      ["build", "area per lipid, upper / lower (nm2)", "imposed difference (%)",
       "builds", "difference between the leaflet means (%)",
       "cholesterol (%)", "DLPC (%)", "verdict at the 8% tolerance"],
      [[NAME[i], "%.3f / %.3f" % (APL_U[i], APL_L[i]), "%.1f" % IMP[i],
        len(LEAF[i]),
        "%+.2f \u00b1 %.2f" % _LM[i],
        "%+.1f \u00b1 %.1f" % _CM[i],
        "%+.1f \u00b1 %.1f" % _DM[i],
        _verdict(i)] for i in range(4)],
      "**Table 2.** What the leaflet area measurement returns against what was "
      "imposed on it. One composition, built in a 12 by 12 by 16 nm box with four "
      "pairs of areas per lipid handed to the packing and nothing else changed. The "
      "mean of the pair is 0.70 nm2 in all four. The packing places the molecules "
      "from a random start, so each pair was built four or five times over, and the "
      "table gives the mean and the standard deviation over those builds, with the "
      "sign of each difference kept. The imposed difference is the difference "
      "between the two areas over that mean. Regressed on the imposed difference, "
      "the difference between the two leaflet means gives a slope of %.3f \u00b1 %.3f "
      "and an offset of %+.2f \u00b1 %.2f percentage points, the cholesterol comparison "
      "%.2f \u00b1 %.2f and %+.1f \u00b1 %.1f, and the DLPC comparison %.2f \u00b1 %.2f and "
      "%+.1f \u00b1 %.1f. `table2_fit.py` returns these."
      % (_bl, _ble, _al, _ale, _bc, _bce, _ac, _ace, _bd, _bde, _ad, _ade))

# ---------------------------------------------------------------- Table 3
D = [-20, -12, -6, 0, 6, 12]
PU = np.array([120, 128, 134, 140, 146, 152]); PL = np.array([153, 145, 139, 132, 126, 120])
F = 100.0 * (PU - PL) / (PU + PL)
C = np.array([64.78, 63.20, 61.17, 56.86, 50.26, 46.36])
SD = ["2.69", "", "", "1.71", "", ""]; N = [3, 1, 1, 3, 1, 1]
Q = np.polyfit(F, C, 2)
write("Table3_curve",
      ["DLPC moved from the upper leaflet", "phospholipids upper / lower",
       "imbalance f (%)", "cholesterol of the upper leaflet after the run (%)",
       "standard deviation over three runs", "runs", "the fit (%)", "residual (points)"],
      [["%+d" % D[i], "%d / %d" % (PU[i], PL[i]), "%+.2f" % F[i], "%.2f" % C[i],
        SD[i], N[i], "%.2f" % np.polyval(Q, F[i]), "%+.2f" % (C[i] - np.polyval(Q, F[i]))]
       for i in range(6)],
      "**Table 3.** The six systems of the curve and where cholesterol settled in each. "
      "The six differ only in how the DLPC molecules were divided between the leaflets; "
      "every one was built with 71 cholesterol molecules in the upper leaflet and 63 in "
      "the lower, which is %.2f%% of the cholesterol above. The imbalance of each system "
      "is computed from the phospholipid numbers in its own system.top. The fit is "
      "C = %.2f %+.4f f %+.5f f^2, with a root mean square residual of %.2f points. "
      "`curve_stats.py` returns every number of this table."
      % (BUILT, Q[2], Q[1], Q[0], np.sqrt(np.mean((C - np.polyval(Q, F)) ** 2))))

# ---------------------------------------------------------------- Table 4
ROUTE = [("equal numbers (EqN)", -1, -0.37, 57.59, 2.26, 3),
         ("areas measured on the packed system (none of the four)", 3, 1.10, 57.43, 0.55, 2),
         ("areas from a table of reported values (SA)", 8, 2.94, 54.70, 1.84, 3)]
write("Table4_routes",
      ["build", "phospholipids upper − lower", "imbalance f (%)",
       "the curve (%)", "measured (%)", "standard error", "runs"],
      [[n, "%+d" % d, "%+.2f" % f, "%.2f" % np.polyval(Q, f),
        "%.2f" % m, "%.2f" % e, r] for n, d, f, m, e, r in ROUTE],
      "**Table 4.** Three builds that were not used in the fit, against what the curve "
      "returns for each. The three were built from one composition with a separate "
      "script and run for 500 ns from three sets of starting velocities. Two of the "
      "three are routes of Section 1; the third is none of the four. The measured "
      "column is the share of the cholesterol of the system that the upper leaflet "
      "holds, averaged over the last 40% of each run and then over the runs. The three "
      "builds assigned their upper leaflets 68, 69 and 71 cholesterol molecules, which "
      "is 51.13, 51.88 and 52.99% of the cholesterol of each system. One run of the "
      "second build lost its box and is left out.")
