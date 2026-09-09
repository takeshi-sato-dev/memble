#!/usr/bin/env python3
"""Collect the runs of compare_counts.sh into one table.

  python3 summarise_counts.py /path/to/output/dir
"""
import glob
import json
import os
import sys

import numpy as np

ARM_LABEL = {"eqn": "equal numbers",
             "table": "numbers from a table",
             "measured": "numbers from the measurement"}


def load(d):
    runs = {}
    for f in sorted(glob.glob(os.path.join(d, "relax_*_s*.json"))):
        base = os.path.basename(f)[len("relax_"):-len(".json")]
        arm, seed = base.rsplit("_s", 1)
        runs.setdefault(arm, []).append((seed, json.load(open(f))))
    return runs


def mean_sem(v):
    v = np.asarray(v, dtype=float)
    if v.size < 2:
        return float(v.mean()), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(v.size))


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    runs = load(d)
    if not runs:
        raise SystemExit("no relax_*.json under %s" % d)

    sterols = set()
    for rs in runs.values():
        for _, r in rs:
            sterols |= {s for s in r["between_leaflets"] if s.startswith("CHOL")}
    sterol = sorted(sterols)[0] if sterols else None

    print("%-28s %8s %10s %10s %10s" % ("", "n", "sterol", "thickness", "order"))
    print("%-28s %8s %10s %10s %10s"
          % ("", "runs", "moved", "difference", "difference"))
    print("-" * 70)
    rows = {}
    for arm in ("eqn", "table", "measured"):
        if arm not in runs:
            continue
        moved, dth, dop = [], [], []
        for _, r in runs[arm]:
            if sterol:
                moved.append(abs(r["between_leaflets"][sterol]["moved"]))
            t = r["leaflet_thickness_nm"]
            dth.append(abs(t["upper"]["settled"] - t["lower"]["settled"]))
            o = r["tail_order"]
            dop.append(abs(o["upper"]["settled"] - o["lower"]["settled"]))
        rows[arm] = dict(moved=mean_sem(moved) if moved else (float("nan"),) * 2,
                         thickness=mean_sem(dth), order=mean_sem(dop),
                         n=len(runs[arm]))
        r = rows[arm]
        print("%-28s %8d %5.1f+-%3.1f %5.3f+-%5.3f %5.3f+-%5.3f"
              % (ARM_LABEL.get(arm, arm), r["n"],
                 r["moved"][0], r["moved"][1],
                 r["thickness"][0], r["thickness"][1],
                 r["order"][0], r["order"][1]))
    print("-" * 70)
    print("sterol moved: molecules that changed leaflet between the start and")
    print("  the second half of the run. thickness: the difference between the")
    print("  two leaflets, in nm. order: the difference in the tail order")
    print("  parameter between the two leaflets. Each is an absolute value,")
    print("  averaged over the runs of the arm, with the standard error.")

    with open(os.path.join(d, "counts_summary.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    print("written to %s" % os.path.join(d, "counts_summary.json"))


if __name__ == "__main__":
    main()
