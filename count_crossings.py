#!/usr/bin/env python3
"""How often a lipid changes leaflet, counted one molecule at a time.

leaflet_relax.py reports how many molecules of each species each leaflet holds
in every frame, which is a net number: a molecule that goes up while another
goes down leaves it unchanged. This script counts the crossings themselves.

The two numbers answer two different questions. The net says where the
composition settled. The gross says whether the molecules are still moving
after it settled, and a species whose net stops changing while its molecules
keep crossing is a species at equilibrium, not a species that was frozen in
place by a barrier it cannot pass.

The leaflet of a molecule is read exactly as leaflet_relax.py reads it, from
the molecule itself: the head bead of a lipid in the upper leaflet lies above
the rest of that same molecule. No plane is drawn through the bilayer.

A crossing is counted only when a molecule stays on the new side for
--dwell frames. A molecule whose head and body sit almost level flickers
between the two answers from one frame to the next, and those flickers are not
crossings. The default dwell of 2 frames removes them; --dwell 1 counts every
change and is what a check of this filter compares against.

  python3 count_crossings.py --gro system.gro --xtc prod.xtc --json crossings.json
"""
import argparse
import json

import numpy as np

try:
    import mdtraj as md
except ImportError:                                    # pragma: no cover
    md = None

WATER = {"W", "WF", "ION", "NA", "CL", "NA+", "CL-"}
STEROL_HEAD = {"ROH"}
PHOSPHATE = {"PO4"}


def lipid_residues(top):
    """Every residue that is neither protein, water nor ion."""
    out = []
    for r in top.residues:
        nm = r.name.strip()
        if nm in WATER or r.is_protein:
            continue
        names = {a.name for a in r.atoms}
        if not (names & (PHOSPHATE | STEROL_HEAD)) and len(names) < 3:
            continue                                   # a free ion
        out.append(r)
    return out


def leaflet_of_every_lipid(t, res):
    """upper[frame, molecule]: True when that molecule sits in the upper leaflet."""
    z = t.xyz[:, :, 2]
    by_res, heads = [], []
    for r in res:
        by_res.append(np.array([a.index for a in r.atoms]))
        heads.append(np.array([a.index for a in r.atoms
                               if a.name in PHOSPHATE | STEROL_HEAD]))
    head_idx = [heads[i] if len(heads[i]) else by_res[i] for i in range(len(res))]
    body_idx = []
    for i in range(len(res)):
        h = set(int(x) for x in head_idx[i])
        b = np.array([int(x) for x in by_res[i] if int(x) not in h])
        body_idx.append(b if b.size else head_idx[i])
    zhead = np.stack([z[:, idx].mean(axis=1) for idx in head_idx], axis=1)
    zbody = np.stack([z[:, idx].mean(axis=1) for idx in body_idx], axis=1)
    return zhead > zbody


def settle(upper, dwell):
    """Remove the flicker: a side is taken only once it is held for `dwell` frames."""
    if dwell <= 1:
        return upper.copy()
    n_f, n_m = upper.shape
    out = np.empty_like(upper)
    state = upper[0].copy()
    run = np.ones(n_m, dtype=int)
    out[0] = state
    for f in range(1, n_f):
        same = upper[f] == upper[f - 1]
        run = np.where(same, run + 1, 1)
        flip = (upper[f] != state) & (run >= dwell)
        state = np.where(flip, upper[f], state)
        out[f] = state
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--xtc", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dwell", type=int, default=2,
                    help="frames a molecule must hold a side for the change to "
                         "count as a crossing (default 2; 1 counts every change)")
    ap.add_argument("--settled-from-ns", type=float, default=250.0,
                    help="the time after which the net is taken to have settled; "
                         "the rate is reported separately before and after it")
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    if md is None:
        raise SystemExit("count_crossings.py needs mdtraj:  pip install mdtraj")

    t = md.load(a.xtc, top=a.gro, stride=a.stride)
    res = lipid_residues(t.topology)
    if not res:
        raise SystemExit("no lipid residues found in %s" % a.gro)
    names = [r.name.strip() for r in res]
    time_ns = np.asarray(t.time) / 1000.0

    raw = leaflet_of_every_lipid(t, res)
    upper = settle(raw, a.dwell)

    changes = upper[1:] != upper[:-1]            # [frame-1, molecule]
    up_going = changes & upper[1:]               # lower -> upper
    down_going = changes & ~upper[1:]            # upper -> lower

    # the frame index at which the settled window begins
    mid_t = 0.5 * (time_ns[1:] + time_ns[:-1])
    late = mid_t >= a.settled_from_ns
    early = ~late
    span_early = float(mid_t[early].max() - time_ns[0]) if early.any() else 0.0
    span_late = float(time_ns[-1] - a.settled_from_ns) if late.any() else 0.0

    out = {"n_frames": int(t.n_frames), "n_lipids": len(res),
           "dwell_frames": a.dwell, "stride": a.stride,
           "frame_spacing_ns": float(np.median(np.diff(time_ns))) if t.n_frames > 1 else 0.0,
           "length_ns": float(time_ns[-1] - time_ns[0]),
           "settled_from_ns": a.settled_from_ns,
           "species": {}}

    print("%-8s %7s %9s %9s %9s   %s"
          % ("species", "n", "crossings", "per mol", "net", "crossings per us"))
    for s in sorted(set(names)):
        m = np.array([n == s for n in names])
        tot_up = int(up_going[:, m].sum())
        tot_dn = int(down_going[:, m].sum())
        n0 = int(upper[0, m].sum())
        n1 = int(upper[-1, m].sum())
        e_up = int(up_going[early][:, m].sum()); e_dn = int(down_going[early][:, m].sum())
        l_up = int(up_going[late][:, m].sum()); l_dn = int(down_going[late][:, m].sum())
        nmol = int(m.sum())
        rate_late = ((l_up + l_dn) / span_late * 1000.0) if span_late > 0 else float("nan")
        rate_early = ((e_up + e_dn) / span_early * 1000.0) if span_early > 0 else float("nan")
        out["species"][s] = {
            "n_molecules": nmol,
            "crossings_total": tot_up + tot_dn,
            "crossings_up": tot_up, "crossings_down": tot_dn,
            "crossings_per_molecule": (tot_up + tot_dn) / float(nmol),
            "upper_first_frame": n0, "upper_last_frame": n1,
            "net_change": n1 - n0,
            "before_settled": {"up": e_up, "down": e_dn, "net": e_up - e_dn,
                               "span_ns": span_early,
                               "crossings_per_us": rate_early},
            "after_settled": {"up": l_up, "down": l_dn, "net": l_up - l_dn,
                              "span_ns": span_late,
                              "crossings_per_us": rate_late},
        }
        print("%-8s %7d %9d %9.2f %9d   %8.1f before / %8.1f after"
              % (s, nmol, tot_up + tot_dn, (tot_up + tot_dn) / float(nmol),
                 n1 - n0, rate_early, rate_late))

    print("")
    print("net is where the composition settled. crossings is how often a molecule")
    print("changed side. A species whose net stops changing while its crossings")
    print("continue is at equilibrium; one with no crossings never had the chance")
    print("to reach one, and carries the number the build gave it.")

    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print("written %s" % a.json)


if __name__ == "__main__":
    main()
