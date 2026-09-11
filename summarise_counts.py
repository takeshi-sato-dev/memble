#!/usr/bin/env python3
"""Collect the runs of compare_counts.sh into one table.

  python3 summarise_counts.py /path/to/output/dir
"""
import glob
import importlib.util
import json
import os
import sys

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "leaflet_relax", os.path.join(_here, "leaflet_relax.py"))
_lr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lr)

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


def built_of(d, arm, rec):
    """The numbers the build of this arm gave each leaflet.

    A file written by the present leaflet_relax.py carries them. A file
    written before that carries only the frames, and the build is read from
    the system.top and memble_build.json compare_counts.sh left in
    <d>/<arm>/<arm>_work.

    The first frame of a production run is not the build: the equilibration
    has already moved cholesterol by the time that frame is written, and in a
    bilayer that holds a protein a handful of molecules are assigned to the
    leaflet they did not come from. Measuring the move against that frame
    reports the wrong number, and reports a different wrong number for each
    seed of one build.
    """
    if rec.get("built_counts"):
        return rec["built_counts"], "the json"
    w = os.path.join(d, arm, "%s_work" % arm)
    got = _lr.built_counts(os.path.join(w, "system.top"),
                           os.path.join(w, "memble_build.json"))
    return got, ("system.top" if got else None)


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
        sources = set()
        for _, r in runs[arm]:
            # Every quantity keeps its sign. Which leaflet gave up sterol, and
            # which leaflet ended thinner, is what separates one way of setting
            # the numbers from another, and an absolute value would discard it.
            # An absolute value would also report a run that did nothing as
            # having done something, because it turns the scatter of a
            # measurement centered on zero into a positive mean.
            if sterol:
                built, src = built_of(d, arm, r)
                sources.add(src)
                base = built["upper"].get(sterol) if built else None
                if base is None:
                    moved.append(r["between_leaflets"][sterol]["moved"])
                else:
                    moved.append(r["between_leaflets"][sterol]["settled"] - base)
            t = r["leaflet_thickness_nm"]
            dth.append(t["upper"]["settled"] - t["lower"]["settled"])
            o = r["tail_order"]
            dop.append(o["upper"]["settled"] - o["lower"]["settled"])
        rows[arm] = dict(moved=mean_sem(moved) if moved else (float("nan"),) * 2,
                         thickness=mean_sem(dth), order=mean_sem(dop),
                         n=len(runs[arm]),
                         baseline_from=sorted(x for x in sources if x))
        if None in sources:
            print("  NOTE: the build of the %s arm was not found under %s, and"
                  % (arm, os.path.join(d, arm, "%s_work" % arm)))
            print("  the sterol figure of that arm is taken against the first"
                  " frame of each run.")
        r = rows[arm]
        print("%-28s %8d %+5.1f+-%3.1f %+6.3f+-%5.3f %+6.3f+-%5.3f"
              % (ARM_LABEL.get(arm, arm), r["n"],
                 r["moved"][0], r["moved"][1],
                 r["thickness"][0], r["thickness"][1],
                 r["order"][0], r["order"][1]))
    print("-" * 70)
    print("sterol moved: the change in the number of sterol molecules the upper")
    print("  leaflet holds, between the number the build gave that leaflet and")
    print("  the second half of the production run. The first frame of that run")
    print("  is not the build, because the equilibration has already moved")
    print("  cholesterol by the time the frame is written.")
    print("  thickness: upper minus lower, in nm. order: upper")
    print("  minus lower, in the tail order parameter. Each is averaged over")
    print("  the runs of the arm, with the standard error, and each keeps its")
    print("  sign: a positive sterol figure means the upper leaflet took")
    print("  molecules from the lower one.")

    with open(os.path.join(d, "counts_summary.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    print("written to %s" % os.path.join(d, "counts_summary.json"))


if __name__ == "__main__":
    main()
