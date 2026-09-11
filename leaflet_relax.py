#!/usr/bin/env python3
"""What a bilayer does with the lipid numbers it was built with.

The number of lipids each leaflet receives is fixed at the build. If the two
leaflets were given numbers that do not match their areas, the membrane relieves
the mismatch during the run: the sterol moves to the leaflet that has room, the
stretched leaflet thins, and its tails lose order. This script measures those
three things per leaflet over a trajectory, so two builds of the same
composition that differ only in the numbers can be compared.

No lateral pressure profile is computed, and no claim about differential stress
is made here. What is reported is what the membrane did.

  python3 leaflet_relax.py --gro system.gro --xtc step7.xtc --json relax.json
"""
import argparse
import json
import os

import numpy as np

try:
    import mdtraj as md
except ImportError:                                    # pragma: no cover
    # built_counts below reads two text files and needs no trajectory library.
    # summarise_counts.py and curve_report.py import this module for it alone,
    # and both must run where mdtraj is not installed. Reading a trajectory
    # still needs mdtraj, and main() says so.
    md = None

WATER = {"W", "WF", "ION", "NA", "CL", "NA+", "CL-"}
STEROL_HEAD = {"ROH"}
PHOSPHATE = {"PO4"}


def built_counts(top_path, build_path):
    """The number of molecules of each species the build gave each leaflet.

    memble writes the two leaflets one after the other in the [ molecules ]
    section of system.top, the upper leaflet first, and memble_build.json names
    the species of each leaflet. The number of species in each leaflet is what
    divides the lipid lines of system.top between them.

    This is the baseline the moved numbers below are measured against. The
    first frame of the trajectory is not that baseline, because the
    equilibration has already moved cholesterol by the time that frame is
    written.
    """
    def n_species(spec):
        """memble records a leaflet as "CHOL:1.0 DLPC:0.71 PSM:1"."""
        if isinstance(spec, dict):
            return len(spec)
        return len([x for x in str(spec).split() if ":" in x])

    try:
        with open(build_path) as fh:
            b = json.load(fh)
        n_up, n_lo = n_species(b["upper"]), n_species(b["lower"])
    except Exception:
        return None

    lines, inside = [], False
    try:
        with open(top_path) as fh:
            for ln in fh:
                s = ln.split(";")[0].strip()
                if not s:
                    continue
                if s.startswith("["):
                    inside = s.replace(" ", "").lower() == "[molecules]"
                    continue
                if not inside:
                    continue
                f = s.split()
                if len(f) == 2 and f[1].isdigit():
                    lines.append((f[0], int(f[1])))
    except Exception:
        return None

    lip = [(n, k) for n, k in lines
           if n.upper() not in WATER and not n.startswith("molecule_")]

    if not n_up or not n_lo:
        # A build given one composition for the whole bilayer records no
        # per-leaflet species list. Such a build writes the same species to
        # both leaflets, so the lipid lines divide in the middle.
        half = len(lip) // 2
        if not half or len(lip) % 2:
            return None
        if [n for n, _ in lip[:half]] != [n for n, _ in lip[half:]]:
            return None
        n_up = n_lo = half

    if len(lip) < n_up + n_lo:
        return None

    out = {"upper": {}, "lower": {}}
    for side, chunk in (("upper", lip[:n_up]), ("lower", lip[n_up:n_up + n_lo])):
        for n, k in chunk:
            out[side][n] = out[side].get(n, 0) + k
    return out


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


def tail_bonds(res):
    """Consecutive pairs inside each Martini tail, by the C/D naming."""
    tails = {}
    for a in res.atoms:
        n = a.name
        if len(n) >= 3 and n[0] in "CD" and n[1].isdigit():
            tails.setdefault(n[-1], []).append((int(n[1]), a.index))
    pairs = []
    for _, beads in tails.items():
        beads.sort()
        pairs += [(beads[i][1], beads[i + 1][1]) for i in range(len(beads) - 1)]
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--xtc", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--last-fraction", type=float, default=0.5,
                    help="the fraction of the trajectory averaged for the "
                         "settled values (default the second half)")
    ap.add_argument("--top", default="",
                    help="system.top of the build, which carries the number of "
                         "molecules each leaflet was given (default: system.top "
                         "beside --gro)")
    ap.add_argument("--build", default="",
                    help="memble_build.json of the build (default: beside --gro)")
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    if md is None:
        raise SystemExit("leaflet_relax.py needs mdtraj:  pip install mdtraj")

    work = os.path.dirname(os.path.abspath(a.gro))
    built = built_counts(a.top or os.path.join(work, "system.top"),
                         a.build or os.path.join(work, "memble_build.json"))

    t = md.load(a.xtc, top=a.gro, stride=a.stride)
    top = t.topology
    res = lipid_residues(top)
    if not res:
        raise SystemExit("no lipid was found in %s" % a.gro)

    lipid_atoms = np.array([x.index for r in res for x in r.atoms])
    by_res = [np.array([x.index for x in r.atoms]) for r in res]
    names = [r.name.strip() for r in res]
    heads = [np.array([x.index for x in r.atoms
                       if x.name in PHOSPHATE | STEROL_HEAD]) for r in res]
    bonds = [tail_bonds(r) for r in res]

    z = t.xyz[:, :, 2]                                  # (frame, atom)

    # --- leaflet of every lipid in every frame ---------------------------
    # A lipid points out of the leaflet it belongs to, so the leaflet is read
    # from the molecule itself: the head bead of a lipid in the upper leaflet
    # lies above the rest of that same molecule, and the head bead of a lipid
    # in the lower leaflet lies below it. The head bead is the phosphate of a
    # phospholipid and the hydroxyl of a sterol, and the rest of the molecule
    # is every other bead it carries.
    #
    # Comparing each lipid against a plane drawn through the whole bilayer
    # gives a different answer, and a worse one, for two reasons. The plane has
    # to be placed, and the mean z of every lipid bead does not place it: the
    # two leaflets of an asymmetric bilayer hold different numbers of beads, so
    # that mean sits toward the leaflet that holds more of them and moves as the
    # bilayer relaxes, which carries lipids across a plane that moved under
    # them. And a bilayer that holds a protein is not flat, so a lipid in a
    # region that the protein depresses sits below a plane its own leaflet
    # rises above. Reading the molecule needs no plane at all.
    head_idx = [heads[i] if len(heads[i]) else by_res[i] for i in range(len(res))]
    body_idx = []
    for i in range(len(res)):
        h = set(int(x) for x in head_idx[i])
        b = np.array([int(x) for x in by_res[i] if int(x) not in h])
        body_idx.append(b if b.size else head_idx[i])

    zhead = np.stack([z[:, idx].mean(axis=1) for idx in head_idx], axis=1)
    zbody = np.stack([z[:, idx].mean(axis=1) for idx in body_idx], axis=1)
    upper = zhead > zbody

    # The thickness below is measured from a plane, and that plane is the
    # midpoint between the head beads of the two leaflets, which the assignment
    # above has already settled.
    mid = np.empty(t.n_frames)
    for f in range(t.n_frames):
        u = upper[f]
        mid[f] = (0.5 * (zhead[f][u].mean() + zhead[f][~u].mean())
                  if u.any() and (~u).any() else float(z[f, lipid_atoms].mean()))

    out = {"n_frames": int(t.n_frames), "n_lipids": len(res),
           "time_ps": [float(x) for x in t.time]}

    # --- 1. the sterol, and every other species, per leaflet over time ---
    species = sorted(set(names))
    counts = {}
    for s in species:
        m = np.array([n == s for n in names])
        counts[s] = {"upper": upper[:, m].sum(axis=1).astype(int).tolist(),
                     "lower": (~upper[:, m]).sum(axis=1).astype(int).tolist()}
    out["counts"] = counts

    n0 = a.last_fraction
    cut = int((1.0 - n0) * t.n_frames)

    moved = {}
    for s in species:
        u = np.array(counts[s]["upper"], dtype=float)
        base = built["upper"].get(s) if built else None
        moved[s] = {"built": base,
                    "first_frame": int(u[0]),
                    "settled": float(u[cut:].mean()),
                    "moved": float(u[cut:].mean()
                                   - (base if base is not None else u[0])),
                    "baseline": "system.top" if base is not None else "first frame"}
    out["between_leaflets"] = moved
    out["built_counts"] = built

    # --- 2. the thickness of each leaflet --------------------------------
    po4 = [np.array([x.index for x in r.atoms if x.name in PHOSPHATE])
           for r in res]
    th = {"upper": [], "lower": []}
    for f in range(t.n_frames):
        for side, sel in (("upper", upper[f]), ("lower", ~upper[f])):
            d = [abs(z[f, po4[i]].mean() - mid[f])
                 for i in range(len(res)) if sel[i] and len(po4[i])]
            th[side].append(float(np.mean(d)) if d else float("nan"))
    out["leaflet_thickness_nm"] = {
        k: {"series": v, "settled": float(np.nanmean(v[cut:]))}
        for k, v in th.items()}

    # --- 3. the order of the tails in each leaflet -----------------------
    op = {"upper": [], "lower": []}
    for f in range(t.n_frames):
        for side, sel in (("upper", upper[f]), ("lower", ~upper[f])):
            vals = []
            for i in range(len(res)):
                if not sel[i] or not bonds[i]:
                    continue
                p = np.array(bonds[i])
                v = t.xyz[f, p[:, 1]] - t.xyz[f, p[:, 0]]
                v /= np.linalg.norm(v, axis=1)[:, None]
                vals.append(np.mean(0.5 * (3.0 * v[:, 2] ** 2 - 1.0)))
            op[side].append(float(np.mean(vals)) if vals else float("nan"))
    out["tail_order"] = {k: {"series": v, "settled": float(np.nanmean(v[cut:]))}
                         for k, v in op.items()}

    # --- 4. the area of the membrane -------------------------------------
    area = (t.unitcell_lengths[:, 0] * t.unitcell_lengths[:, 1]).astype(float)
    out["area_nm2"] = {"series": area.tolist(),
                       "settled": float(area[cut:].mean()),
                       "start": float(area[0])}

    print("frames %d, lipids %d, second half from %.1f ns"
          % (t.n_frames, len(res), t.time[cut] / 1000.0))
    if built is None:
        print("  NOTE: system.top and memble_build.json were not both read, so"
              " the moved numbers below are taken against the first frame of"
              " this trajectory and not against the build.")
    for s in species:
        m = moved[s]
        base = m["built"] if m["built"] is not None else m["first_frame"]
        print("  %-8s upper %3d at the build, %3d in the first frame, %6.1f "
              "settled, moved %+5.1f"
              % (s, base, m["first_frame"], m["settled"], m["moved"]))
    print("  thickness   upper %.3f nm, lower %.3f nm, difference %.3f nm"
          % (out["leaflet_thickness_nm"]["upper"]["settled"],
             out["leaflet_thickness_nm"]["lower"]["settled"],
             out["leaflet_thickness_nm"]["upper"]["settled"]
             - out["leaflet_thickness_nm"]["lower"]["settled"]))
    print("  tail order  upper %.3f, lower %.3f, difference %.3f"
          % (out["tail_order"]["upper"]["settled"],
             out["tail_order"]["lower"]["settled"],
             out["tail_order"]["upper"]["settled"]
             - out["tail_order"]["lower"]["settled"]))
    print("  area        %.1f nm^2 at the start, %.1f nm^2 settled"
          % (out["area_nm2"]["start"], out["area_nm2"]["settled"]))

    if a.json:
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=1)
        print("written to %s" % a.json)


if __name__ == "__main__":
    main()
