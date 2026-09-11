#!/usr/bin/env python3
"""The points of a response curve, read against the builds they came from.

run_curve.sh and run_curve2.sh each build one system per point, run it, and
write the leaflet counts of every frame to relax_<tag>.json beside the build
directory. This script collects those files into the curve.

Both axes are taken from the build and not from the first frame of the run.
The number of phospholipids each leaflet was given is the quantity the curve
is drawn against, and the first frame does not report that number: a lipid
whose head bead sits near the midplane is assigned to a leaflet frame by
frame, and in a bilayer that holds a protein a handful of molecules are
assigned to the leaflet they did not come from. The build carries the numbers
that were imposed, and system.top holds them.

  python3 curve_report.py ~/counts_run/curve
  python3 curve_report.py ~/counts_run/curve2 --prefix relax_a
"""
import argparse
import glob
import importlib.util
import json
import os
import re

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "leaflet_relax", os.path.join(_here, "leaflet_relax.py"))
_lr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lr)

STEROL = "CHOL"


def work_dir(d, tag):
    """run_curve.sh puts the build of point <tag> in <d>/<tag>/<tag>_work."""
    return os.path.join(d, tag, "%s_work" % tag)


def built_of(d, tag, rec):
    """The numbers the build gave each leaflet.

    A file written by the present leaflet_relax.py carries them. A file
    written before that carries only the frames, and the build is read from
    system.top and memble_build.json of the point.
    """
    if rec.get("built_counts"):
        return rec["built_counts"], "the json"
    w = work_dir(d, tag)
    got = _lr.built_counts(os.path.join(w, "system.top"),
                           os.path.join(w, "memble_build.json"))
    return got, ("system.top" if got else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--prefix", default="relax_d",
                    help="relax_d for run_curve.sh, relax_a for run_curve2.sh")
    ap.add_argument("--last-fraction", type=float, default=0.4,
                    help="the fraction of each run averaged for the settled "
                         "value (default the last 40%%, as the runs were made)")
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    # One point of the curve is "<prefix><integer>.json" and nothing else. The
    # same directory also holds the repeats of a point, which carry a seed and
    # a run length in the name (relax_d0_s2_250ns.json), and those are not
    # points of the curve. Take the names that are a point and leave the rest,
    # rather than reading a seed as a delta or stopping on the first name that
    # does not parse.
    pat = re.compile(re.escape(a.prefix) + r"(-?\d+)\.json$")
    keyed = []
    skipped = []
    for f in glob.glob(os.path.join(a.dir, a.prefix + "*.json")):
        m = pat.search(os.path.basename(f))
        if m:
            keyed.append((int(m.group(1)), f))
        else:
            skipped.append(os.path.basename(f))
    files = [f for _, f in sorted(keyed)]
    if not files:
        raise SystemExit("no %s<integer>.json under %s" % (a.prefix, a.dir))
    if skipped:
        print("  (%d file%s in this directory are not points of the curve and "
              "were left out: %s)"
              % (len(skipped), "" if len(skipped) == 1 else "s",
                 ", ".join(sorted(skipped))))

    print("  PL up-lo    CHOL built      CHOL settled          upper share"
          "     moved   first frame")
    rows = []
    for f in files:
        tag = os.path.basename(f)[len("relax_"):-len(".json")]
        rec = json.load(open(f))
        counts = rec["counts"]
        n = len(rec["time_ps"])
        cut = int((1.0 - a.last_fraction) * n)
        u = counts[STEROL]["upper"]
        tot = u[0] + counts[STEROL]["lower"][0]
        settled = sum(u[cut:]) / len(u[cut:])

        built, src = built_of(a.dir, tag, rec)
        if built is None:
            print("  %-8s the build of this point was not found under %s"
                  % (tag, work_dir(a.dir, tag)))
            continue
        ch_u = built["upper"].get(STEROL, 0)
        ch_l = built["lower"].get(STEROL, 0)
        pl_u = sum(k for s, k in built["upper"].items() if s != STEROL)
        pl_l = sum(k for s, k in built["lower"].items() if s != STEROL)

        print("  %+8d  %4d/%-4d (%.1f%%)  %6.1f/%-6.1f (%.2f%%)  %+6.2f pt "
              "%+6.1f   %4d (%s)"
              % (pl_u - pl_l, ch_u, ch_l, 100.0 * ch_u / tot,
                 settled, tot - settled, 100.0 * settled / tot,
                 100.0 * settled / tot - 100.0 * ch_u / tot,
                 settled - ch_u, u[0], src))
        rows.append({"point": tag, "pl_upper": pl_u, "pl_lower": pl_l,
                     "pl_difference": pl_u - pl_l,
                     "chol_upper_built": ch_u, "chol_lower_built": ch_l,
                     "chol_total": tot,
                     "chol_upper_settled": round(settled, 2),
                     "upper_share_built": round(100.0 * ch_u / tot, 2),
                     "upper_share_settled": round(100.0 * settled / tot, 2),
                     "moved": round(settled - ch_u, 2),
                     "chol_upper_first_frame": u[0],
                     "baseline_from": src})

    print("")
    print("  PL up-lo is the difference between the numbers of phospholipid")
    print("  molecules the build gave the two leaflets. CHOL built is what the")
    print("  build gave each leaflet, and CHOL settled is the mean over the")
    print("  last %d%% of the run. moved is the second less the first. The last"
          % round(100 * a.last_fraction))
    print("  column is the first frame of the run, which the equilibration has")
    print("  already moved, and which is not the baseline used here.")

    if a.csv:
        import csv
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("written to %s" % a.csv)


if __name__ == "__main__":
    main()
