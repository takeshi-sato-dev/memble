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

The trajectory is read in chunks and only the beads of the lipids are read, so
a run of several microseconds is counted on a machine that could not hold it.
Nothing but the state of each molecule is carried from one chunk to the next.

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


def head_and_body(res):
    """For every molecule, the beads of its head and the beads of the rest.

    A lipid points out of the leaflet it belongs to, so the leaflet is read
    from the molecule itself: the head bead of a lipid in the upper leaflet
    lies above the rest of that same molecule, and the head bead of a lipid in
    the lower leaflet lies below it. The head bead is the phosphate of a
    phospholipid and the hydroxyl of a sterol. Comparing each lipid against a
    plane drawn through the bilayer gives a different answer and a worse one,
    because the plane has to be placed and a bilayer holding a protein is not
    flat. Reading the molecule needs no plane at all. This is the rule
    leaflet_relax.py uses, and the two must agree.
    """
    head, body = [], []
    for r in res:
        allx = [a.index for a in r.atoms]
        h = [a.index for a in r.atoms if a.name in PHOSPHATE | STEROL_HEAD]
        if not h:
            h = allx
        hs = set(h)
        b = [x for x in allx if x not in hs] or h
        head.append(np.array(h))
        body.append(np.array(b))
    return head, body


def upper_of_chunk(z, head_local, body_local):
    """upper[frame, molecule] for one chunk of frames."""
    zhead = np.stack([z[:, idx].mean(axis=1) for idx in head_local], axis=1)
    zbody = np.stack([z[:, idx].mean(axis=1) for idx in body_local], axis=1)
    return zhead > zbody


class Crossings(object):
    """The state each molecule carries from one frame to the next.

    Only three arrays cross a chunk boundary: the side the molecule is held to
    be on, the side the last frame read gave, and how long that last answer has
    been repeated. The counts themselves are sums.
    """

    def __init__(self, n_mol, dwell):
        self.dwell = max(1, dwell)
        self.n = n_mol
        self.state = None
        self.prev_raw = None
        self.run = np.ones(n_mol, dtype=np.int64)
        self.up_early = np.zeros(n_mol, dtype=np.int64)
        self.dn_early = np.zeros(n_mol, dtype=np.int64)
        self.up_late = np.zeros(n_mol, dtype=np.int64)
        self.dn_late = np.zeros(n_mol, dtype=np.int64)
        self.first = None
        self.last = None
        self.prev_time = None
        self.n_frames = 0

    def feed(self, upper, times, settled_from_ns):
        for f in range(upper.shape[0]):
            raw = upper[f]
            t = times[f]
            if self.state is None:
                self.state = raw.copy()
                self.prev_raw = raw.copy()
                self.first = raw.copy()
                self.prev_time = t
                self.last = raw.copy()
                self.n_frames = 1
                continue
            same = raw == self.prev_raw
            self.run = np.where(same, self.run + 1, 1)
            flip = (raw != self.state) & (self.run >= self.dwell)
            if flip.any():
                mid = 0.5 * (t + self.prev_time)
                up = flip & raw
                dn = flip & ~raw
                if mid >= settled_from_ns:
                    self.up_late += up
                    self.dn_late += dn
                else:
                    self.up_early += up
                    self.dn_early += dn
                self.state = np.where(flip, raw, self.state)
            self.prev_raw = raw
            self.prev_time = t
            self.last = self.state.copy()
            self.n_frames += 1


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
    ap.add_argument("--chunk", type=int, default=200,
                    help="frames read at a time (default 200); lower it on a "
                         "machine with little memory")
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    if md is None:
        raise SystemExit("count_crossings.py needs mdtraj:  pip install mdtraj")

    top = md.load(a.gro).topology
    res = lipid_residues(top)
    if not res:
        raise SystemExit("no lipid residues found in %s" % a.gro)
    names = [r.name.strip() for r in res]
    head, body = head_and_body(res)

    # Read only the beads of the lipids. A trajectory that also holds the water
    # is several times larger than the part this measurement needs.
    sel = np.unique(np.concatenate([np.concatenate(head), np.concatenate(body)]))
    where = {int(g): i for i, g in enumerate(sel)}
    head_l = [np.array([where[int(x)] for x in h]) for h in head]
    body_l = [np.array([where[int(x)] for x in b]) for b in body]

    cr = Crossings(len(res), a.dwell)
    t_first = t_last = None
    for chunk in md.iterload(a.xtc, top=top, stride=a.stride,
                             chunk=a.chunk, atom_indices=sel):
        times = np.asarray(chunk.time) / 1000.0
        if t_first is None:
            t_first = float(times[0])
        t_last = float(times[-1])
        cr.feed(upper_of_chunk(chunk.xyz[:, :, 2], head_l, body_l),
                times, a.settled_from_ns)
    if cr.state is None:
        raise SystemExit("no frames were read from %s" % a.xtc)

    span_early = max(0.0, min(a.settled_from_ns, t_last) - t_first)
    span_late = max(0.0, t_last - max(a.settled_from_ns, t_first))

    out = {"n_frames": int(cr.n_frames), "n_lipids": len(res),
           "dwell_frames": a.dwell, "stride": a.stride,
           "length_ns": float(t_last - t_first),
           "settled_from_ns": a.settled_from_ns,
           "species": {}}

    print("%-8s %7s %9s %9s %9s   %s"
          % ("species", "n", "crossings", "per mol", "net", "crossings per us"))
    for s_name in sorted(set(names)):
        m = np.array([n == s_name for n in names])
        tot_up = int(cr.up_early[m].sum() + cr.up_late[m].sum())
        tot_dn = int(cr.dn_early[m].sum() + cr.dn_late[m].sum())
        n0 = int(cr.first[m].sum())
        n1 = int(cr.last[m].sum())
        e_up = int(cr.up_early[m].sum()); e_dn = int(cr.dn_early[m].sum())
        l_up = int(cr.up_late[m].sum()); l_dn = int(cr.dn_late[m].sum())
        nmol = int(m.sum())
        rate_late = ((l_up + l_dn) / span_late * 1000.0) if span_late > 0 else float("nan")
        rate_early = ((e_up + e_dn) / span_early * 1000.0) if span_early > 0 else float("nan")
        crossed = (cr.up_early + cr.dn_early + cr.up_late + cr.dn_late) > 0
        out["species"][s_name] = {
            "n_molecules": nmol,
            "n_molecules_that_crossed": int(crossed[m].sum()),
            "most_by_one_molecule": int(
                (cr.up_early + cr.dn_early + cr.up_late + cr.dn_late)[m].max()),
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
        print("%-8s %7d %9d %9.2f %9d   %8.1f before / %8.1f after   "
              "%d of %d molecules crossed, most %d times"
              % (s_name, nmol, tot_up + tot_dn, (tot_up + tot_dn) / float(nmol),
                 n1 - n0, rate_early, rate_late,
                 int(crossed[m].sum()), nmol,
                 int((cr.up_early + cr.dn_early + cr.up_late + cr.dn_late)[m].max())))

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
