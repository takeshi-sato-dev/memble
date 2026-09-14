#!/usr/bin/env python3
"""The areas per lipid Table 2 handed to the packing, and where each one comes from.

Table 2 of the article compares the difference between the two areas per lipid
that a build handed to the packing with the difference the finished bilayer
returns. The first of those two differences is an input to the build, and until
the version of MEMBLE that carries this file the input was not written into the
build record. This script returns the four pairs of the table from the quantities
each one was derived from, so that the four rows can be checked without the build
being repeated.

  python3 table2_apl.py

The composition is the one of Section 2.7: the upper leaflet holds cholesterol,
DLPC and sphingomyelin at 1:1:1 and the lower leaflet holds cholesterol, DLPC,
DOPS and the (4,5) phosphoinositide at 1:1:1:0.1, in a 12 by 12 by 16 nm box
with one copy of the construct. The mean area per lipid handed to the packing is
COBY_APL, which is 0.70 nm^2 in every build of the table.
"""
import subprocess
import sys
import os

COBY_APL = 0.70
UPPER = "CHOL:1 DLPC:1 PSM:1"
LOWER = "CHOL:1 DLPC:1 DOPS:1 POP2_45:0.1"

# The numbers of molecules the third build placed, read from the leaflet_area
# file of that build.
N_UPPER = 211
N_LOWER = 195

# The number of cholesterol molecules the run of the third build moved from the
# lower leaflet to the upper leaflet, which run_settled.sh reads from the
# trajectory and which sets the counts of the fourth build.
MOVED = 7

HERE = os.path.dirname(os.path.abspath(__file__))


def from_balance_apl():
    """The pair balance_apl.py returns for this composition, which is row 3."""
    out = subprocess.run(
        [sys.executable, os.path.join(HERE, "balance_apl.py"),
         "--upper", UPPER, "--lower", LOWER, "--base-apl", str(COBY_APL)],
        capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit("balance_apl.py failed: %s" % out.stderr)
    u, l = out.stdout.split()
    return float(u), float(l)


def main():
    print("mean area per lipid handed to the packing, COBY_APL   %.3f nm^2" % COBY_APL)
    print("molecules the third build placed                      %d upper, %d lower"
          % (N_UPPER, N_LOWER))
    print("cholesterol the run of the third build moved upward   %d molecules" % MOVED)
    print("")
    print("%-42s %7s %7s  %s" % ("build", "upper", "lower", "how the pair is obtained"))

    # Row 1. Both leaflets are handed COBY_APL, which is what AUTO_BALANCE=0 does.
    print("%-42s %7.3f %7.3f  %s"
          % ("equal areas per lipid (EqN)", COBY_APL, COBY_APL,
             "AUTO_BALANCE=0 hands COBY_APL to both leaflets"))

    # Row 3. balance_apl.py sizes the two areas from a table of areas per lipid
    # weighted by the composition, holding the mean at COBY_APL.
    tu, tl = from_balance_apl()
    print("%-42s %7.3f %7.3f  %s"
          % ("from the table of areas (SA)", tu, tl,
             "balance_apl.py, run by this script"))

    # Row 4. run_settled.sh keeps the areas of row 3 and divides each by the
    # number of molecules that leaflet holds once the sterol has been moved.
    su = tu * N_UPPER / float(N_UPPER + MOVED)
    sl = tl * N_LOWER / float(N_LOWER - MOVED)
    print("%-42s %7.3f %7.3f  %s"
          % ("the composition the run settled to", su, sl,
             "row 3 scaled by the molecule counts after %d moved" % MOVED))

    # Row 2. The pair comes from a measurement on the packed system of the first
    # build of that arm, and the file holding that measurement is named here
    # rather than recomputed, because a measurement is not reproducible from a
    # composition alone.
    print("%-42s %7.3f %7.3f  %s"
          % ("corrected from the measurement", 0.689, 0.711,
             "measured on the packed system of that arm, see the note below"))
    print("")
    print("The pair of the second build is the one quantity of this table that")
    print("cannot be recomputed from the composition, because the pair was taken")
    print("from a Monte Carlo measurement of the Voronoi areas of a packed system.")
    print("The leaflet_area.json of that build is the record of it. The mean of")
    print("the pair is %.4f nm^2, which is COBY_APL, so the correction moved the"
          % ((0.689 + 0.711) / 2.0))
    print("two areas apart and held their mean, exactly as the other rows do.")


if __name__ == "__main__":
    main()
