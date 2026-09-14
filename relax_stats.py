#!/usr/bin/env python3
"""Every number Section 3.3 quotes for the 8 microsecond run, from the counts file.

The counts file is the one `leaflet_relax.py` writes: for every frame it holds
the number of molecules of each species in each leaflet. Section 3.3 quotes a
first value, a settled value, a move, the course of that move over four windows,
a slope over the second half of the run, and the scatter between blocks of that
half. This script returns all of them from that one file.

Two conventions are fixed here and are the ones the article uses.

  A time window written "a to b ns" is half open: a frame is inside the window
  when its time is at least a and less than b. A frame at exactly b belongs to
  the next window and is counted once.

  A scatter over blocks is the sample standard deviation over those blocks,
  with one degree of freedom removed.

  python3 relax_stats.py dipc_8us.json --built 501 --species CHOL
"""
import argparse
import json

import numpy as np


def window(t, a, b):
    """Half open: a <= time < b."""
    return (t >= a) & (t < b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("counts_json")
    ap.add_argument("--species", default="CHOL")
    ap.add_argument("--built", type=float, required=True,
                    help="molecules of that species the build gave the upper "
                         "leaflet, read from system.top and not from frame one")
    ap.add_argument("--settled-from-ns", type=float, default=1000.0)
    ap.add_argument("--block-ns", type=float, default=500.0)
    ap.add_argument("--blocks-from-ns", type=float, default=2000.0)
    a = ap.parse_args()

    d = json.load(open(a.counts_json))
    t = np.asarray(d["time_ps"], dtype=float) / 1000.0
    up = np.asarray(d["counts"][a.species]["upper"], dtype=float)
    lo = np.asarray(d["counts"][a.species]["lower"], dtype=float)
    total = up[0] + lo[0]

    print("file        %s" % a.counts_json)
    print("species     %s, %d molecules in the system" % (a.species, total))
    print("frames      %d, spanning %.0f ns, spaced %.1f ns"
          % (len(t), t[-1] - t[0], t[1] - t[0]))
    print("")
    print("the build gave the upper leaflet          %d molecules, %.1f%%"
          % (a.built, 100.0 * a.built / total))
    print("the first frame of the run holds          %d molecules, %.2f%%"
          % (up[0], 100.0 * up[0] / total))

    settled = up[t >= a.settled_from_ns].mean()
    moved = settled - a.built
    print("mean from %.0f ns to the end               %.1f molecules, %.2f%%"
          % (a.settled_from_ns, settled, 100.0 * settled / total))
    print("moved against the build                   %.1f molecules, %.2f%% of the species"
          % (moved, 100.0 * moved / total))
    print("")

    print("the course of the move, each window half open")
    for lo_ns, hi_ns in ((0.0, 50.0), (50.0, 100.0), (100.0, 250.0)):
        w = window(t, lo_ns, hi_ns)
        print("  %6.0f to %6.0f ns   gained against the build %6.1f molecules"
              % (lo_ns, hi_ns, up[w].mean() - a.built))
    per_us = []
    k = int(np.ceil(t[-1] / 1000.0))
    for i in range(1, k):
        w = window(t, i * 1000.0, (i + 1) * 1000.0)
        if w.any():
            per_us.append(up[w].mean() - a.built)
    print("  each microsecond from 1 us on  %s"
          % " ".join("%.1f" % v for v in per_us))
    print("  the range over those microseconds   %.1f to %.1f molecules"
          % (min(per_us), max(per_us)))
    print("")

    w = t >= a.blocks_from_ns
    slope = np.polyfit(t[w] / 1000.0, up[w], 1)[0]
    print("a straight line from %.0f ns to the end    %+.2f molecules per microsecond"
          % (a.blocks_from_ns, slope))

    edges = np.arange(a.blocks_from_ns, t[-1], a.block_ns)
    blocks = [up[window(t, e, e + a.block_ns)].mean() for e in edges]
    blocks = [b for b in blocks if np.isfinite(b)]
    print("%.0f ns blocks after %.0f ns              %d blocks, sample sd %.1f molecules"
          % (a.block_ns, a.blocks_from_ns, len(blocks), np.std(blocks, ddof=1)))

    print("")
    print("every other species, first frame and settled mean")
    for sp in sorted(d["counts"]):
        u2 = np.asarray(d["counts"][sp]["upper"], dtype=float)
        print("  %-8s first %5d   settled %7.1f   change %+6.1f"
              % (sp, u2[0], u2[t >= a.settled_from_ns].mean(),
                 u2[t >= a.settled_from_ns].mean() - u2[0]))


if __name__ == "__main__":
    main()
