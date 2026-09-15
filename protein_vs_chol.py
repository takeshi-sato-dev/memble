#!/usr/bin/env python3
"""What the protein and each leaflet do across the systems of the curve.

The six systems of Section 3.5 differ only in how the phospholipid molecules
were divided between the two leaflets, and the run left them holding 29.0 to
42.0 mol% cholesterol in the leaflet that holds the sphingomyelin. Each holds
one copy of the construct. That makes the six a series in the cholesterol
content of a leaflet, built for another purpose and already run.

Version 2. Three things are different from version 1, and each was a fault of
version 1 rather than a result:

  Every selection reports how many beads it matched. A contact count of zero
  means no contact only when the selection found the beads to count, and
  version 1 could not tell the two apart.

  The quantities are taken per leaflet. Cholesterol leaves one leaflet and
  arrives in the other, so a quantity summed over the bilayer has the two
  changes cancelling inside it. The thickness of the whole bilayer was flat
  across the six systems for that reason and not because nothing happened.

  Replicate runs of one build are read when they are present, because the
  spread between them is the only yardstick for the spread across the six.

Measured per leaflet:
  thickness   the phosphate plane of that leaflet, above the midplane, in nm
  order       the second rank order parameter of the tail bonds of that
              leaflet against the bilayer normal
Measured on the protein:
  tilt        the angle between the membrane-spanning helix and the normal
  contacts    anionic head group beads within 0.7 nm of a basic side chain of
              the juxtamembrane region

Usage:
    python3 -u protein_vs_chol.py ~/counts_run/curve
    python3 -u protein_vs_chol.py ~/counts_run/curve --replicates
    python3 -u protein_vs_chol.py ~/counts_run/curve --tm 65:88 --jm 89:117
"""
import argparse
import glob
import json
import os
import re

import numpy as np

STEROL = "CHOL"
ANIONIC = ("DOPS", "POP2", "POPI", "PIP2", "POP1", "POP3")
BASIC = ("ARG", "LYS")
HEADNAMES = ("PO4", "P1", "P2", "P3", "CNO", "PO1", "PO2")
CUTOFF_NM = 0.7
TAIL = re.compile(r"^[CD]\d[AB]$")


def parse_range(s):
    a, b = s.split(":")
    return int(a), int(b)


def find_systems(root, prefix="relax_d", replicates=False):
    """Match every relax json to the run it belongs to.

    The layout is  curve/d<delta>/d<delta>_work/  and inside it

        prod.xtc              the run the curve was fitted on
        s<n>_prod.xtc         replicate n
        s<n>_prod250.xtc      replicate n cut to 250 ns, where both lengths
                              were kept

    so relax_d0_s2_250ns.json belongs to s2_prod250.xtc and relax_d0_s2_500ns
    to s2_prod.xtc. A json whose run cannot be found is reported rather than
    silently paired with the base run, which is what made every replicate
    return the numbers of its own base system.
    """
    pat = re.compile(re.escape(prefix) + r"(-?\d+)(?:_s(\d+))?(?:_(\d+)ns)?\.json$")
    out = []
    for js in sorted(glob.glob(os.path.join(root, prefix + "*.json"))):
        m = pat.search(os.path.basename(js))
        if not m:
            continue
        delta, rep, ns = int(m.group(1)), m.group(2), m.group(3)
        if rep and not replicates:
            continue
        work = os.path.join(root, "d%d" % delta, "d%d_work" % delta)
        if rep and ns:
            tnames = ["s%s_prod%s.xtc" % (rep, ns), "s%s_prod.xtc" % rep]
            if ns == "500":
                tnames = ["s%s_prod.xtc" % rep, "s%s_prod500.xtc" % rep]
        elif rep:
            tnames = ["s%s_prod.xtc" % rep, "s%s_prod250.xtc" % rep]
        else:
            tnames = ["prod.xtc", "prod250.xtc", "md.xtc", "traj_comp.xtc"]
        snames = ([] if not rep else ["s%s_prod.gro" % rep]) + \
                 ["prod.gro", "md.gro", "step7.gro", "system.gro"]
        struct = traj = None
        for d in (work, os.path.join(root, "d%d" % delta), root):
            if not os.path.isdir(d):
                continue
            for name in tnames:
                p = os.path.join(d, name)
                if traj is None and os.path.exists(p):
                    traj = p
            for name in snames:
                p = os.path.join(d, name)
                if struct is None and os.path.exists(p):
                    struct = p
            if struct and traj:
                break
        out.append((delta, rep, ns, js, struct, traj))
    return out


def settled_chol(js, work):
    """Cholesterol mol% of each leaflet after the run.

    The phospholipid numbers come from what the build recorded and never from
    counting the packed system. Counting misassigns a handful of molecules
    wherever the protein bends the membrane, and the error is not small: the
    system at delta = -20 was built with 120 phospholipid molecules above and
    153 below, and counting returns 132 and 141 whichever rule is used. The
    build wrote the numbers down, so they are read.

    built_counts in the relax json is preferred; memble_build.json in the work
    directory is the fallback. A system with neither is skipped rather than
    guessed at.
    """
    rec = json.load(open(js))
    c = rec["counts"][STEROL]
    n = len(c["upper"])
    k = int(0.6 * n)
    u = float(np.mean(c["upper"][k:]))
    l = float(np.mean(c["lower"][k:]))

    b = rec.get("built_counts")
    src = "relax json built_counts"
    if not b and work:
        mb = os.path.join(work, "memble_build.json")
        if os.path.exists(mb):
            d = json.load(open(mb))
            for key in ("composition", "leaflets", "lipids", None):
                cand = d.get(key) if key else d
                if isinstance(cand, dict) and "upper" in cand and "lower" in cand:
                    b = cand
                    src = "memble_build.json"
                    break
    if not b:
        return None
    pu = sum(v for kk, v in b["upper"].items() if kk != STEROL)
    pl = sum(v for kk, v in b["lower"].items() if kk != STEROL)
    return u, l, int(pu), int(pl), src


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--last-fraction", type=float, default=0.4)
    ap.add_argument("--tm", default="65:88")
    ap.add_argument("--jm", default="89:117")
    ap.add_argument("--replicates", action="store_true",
                    help="also read relax_d<delta>_s<n>.json, which is what "
                         "gives the spread between runs of one build")
    a = ap.parse_args()
    TM_A, TM_B = parse_range(a.tm)
    JM_A, JM_B = parse_range(a.jm)

    try:
        import MDAnalysis as mda
        from MDAnalysis.lib.distances import distance_array
    except ImportError:
        raise SystemExit("protein_vs_chol: MDAnalysis is not importable here. "
                         "Run this where version 1 ran, which reported using "
                         "MDAnalysis.")
    print("protein_vs_chol v2: MDAnalysis %s" % mda.__version__)

    systems = find_systems(a.dir, replicates=a.replicates)
    if not systems:
        raise SystemExit("no relax_d*.json under %s" % a.dir)

    rows = []
    first = True
    for delta, rep, ns, js, struct, traj in systems:
        tag = ("d%+d" % delta + ("_s%s" % rep if rep else "")
               + ("_%sns" % ns if ns else ""))
        if not struct or not traj:
            print("%-9s structure or trajectory not found" % tag)
            continue
        u = mda.Universe(struct, traj)
        tm = u.select_atoms("name BB and resid %d:%d" % (TM_A, TM_B))
        jm = u.select_atoms("resname %s and resid %d:%d and not name BB"
                            % (" ".join(BASIC), JM_A, JM_B))
        an = u.select_atoms("resname %s and name %s"
                            % (" ".join(ANIONIC), " ".join(HEADNAMES)))
        po4 = u.select_atoms("name PO4")
        lipid = u.select_atoms("not protein and not resname W ION NA CL")
        tails = u.select_atoms("name %s"
                               % " ".join(sorted({at.name for at in lipid
                                                  if TAIL.match(at.name)})
                                          or ["NOSUCHBEAD"]))
        if first:
            print("")
            print("selections, on the first system read:")
            print("  TM backbone beads, resid %d:%d      %d" % (TM_A, TM_B, len(tm)))
            print("  basic side chains, resid %d:%d      %d" % (JM_A, JM_B, len(jm)))
            print("  anionic head group beads             %d" % len(an))
            print("  PO4 beads                            %d" % len(po4))
            print("  tail beads                           %d" % len(tails))
            if len(jm) == 0:
                print("  ** no basic side chain matched: the contact count "
                      "below is meaningless. Check the residue numbering with "
                      "--jm, or the construct has no Arg or Lys there.")
            if len(an) == 0:
                print("  ** no anionic head group matched: check the bead "
                      "names of DOPS and the phosphoinositide.")
            print("")
            print("%-9s %-8s %-8s %-7s %-8s %-8s %-8s %-8s %-7s"
                  % ("system", "CHOLup%", "CHOLlo%", "tilt", "contacts",
                     "d_up nm", "d_lo nm", "S_up-lo", "frames"))
            first = False

        n = len(u.trajectory)
        start = int((1.0 - a.last_fraction) * n)
        tilt, con, dup, dlo, sup, slo = [], [], [], [], [], []
        for ts in u.trajectory[start:]:
            mid = lipid.positions[:, 2].mean()
            if len(tm) >= 4:
                x = tm.positions - tm.positions.mean(axis=0)
                ax = np.linalg.svd(x, full_matrices=False)[2][0]
                tilt.append(np.degrees(np.arccos(abs(ax[2]) / np.linalg.norm(ax))))
            if len(jm) and len(an):
                d = distance_array(jm.positions, an.positions, box=ts.dimensions)
                con.append(int((d < CUTOFF_NM * 10.0).sum()))
            if len(po4):
                z = po4.positions[:, 2]
                dup.append((z[z > mid].mean() - mid) / 10.0)
                dlo.append((mid - z[z < mid].mean()) / 10.0)
            if len(tails):
                for store, keep in ((sup, True), (slo, False)):
                    vals = []
                    for res in tails.residues:
                        b = res.atoms.select_atoms("name %s"
                                                   % " ".join(sorted(
                                                       {at.name for at in res.atoms
                                                        if TAIL.match(at.name)})))
                        if len(b) < 2:
                            continue
                        if (b.positions[:, 2].mean() > mid) != keep:
                            continue
                        v = np.diff(b.positions, axis=0)
                        nrm = np.linalg.norm(v, axis=1)
                        ok = nrm > 1e-6
                        if not ok.any():
                            continue
                        cz = v[ok, 2] / nrm[ok]
                        vals.append(np.mean(0.5 * (3 * cz ** 2 - 1)))
                    store.append(np.mean(vals) if vals else np.nan)

        got = settled_chol(js, os.path.dirname(traj))
        if got is None:
            print("%-11s built leaflet numbers not recorded; skipped" % tag)
            continue
        cu, cl, pu, pl, src = got
        if False:
            pu = pl = 0
            for r in lipid.residues:
                if r.resname == STEROL:
                    continue
                h = r.atoms.select_atoms("name %s" % " ".join(HEADNAMES))
                t = [at.index for at in r.atoms if TAIL.match(at.name)]
                if not len(h) or not t:
                    continue
                if h.positions[0, 2] > u.atoms.positions[t, 2].mean():
                    pu += 1
                else:
                    pl += 1
            print("           %s: phospholipid numbers counted by the head "
                  "against the tails of the same molecule, %d upper %d lower"
                  % (tag, pu, pl))
        up_pct = 100.0 * cu / (pu + cu)
        lo_pct = 100.0 * cl / (pl + cl)
        ds = (np.nanmean(sup) - np.nanmean(slo)) if sup and slo else np.nan
        rows.append(dict(tag=tag, delta=delta, rep=rep, up=up_pct, lo=lo_pct,
                         tilt=np.mean(tilt) if tilt else np.nan,
                         con=np.mean(con) if con else np.nan,
                         dup=np.mean(dup) if dup else np.nan,
                         dlo=np.mean(dlo) if dlo else np.nan, ds=ds))
        print("%-9s %-8.2f %-8.2f %-7.2f %-8.2f %-8.3f %-8.3f %-8.4f %-7d"
              % (tag, up_pct, lo_pct, rows[-1]["tilt"], rows[-1]["con"],
                 rows[-1]["dup"], rows[-1]["dlo"], ds, n - start))

    singles = [r for r in rows if r["rep"] is None]
    if len(singles) >= 3:
        x = np.array([r["up"] for r in singles])
        print("")
        print("regressed on the cholesterol of the leaflet holding the "
              "sphingomyelin, %.1f to %.1f mol%%, %d systems"
              % (x.min(), x.max(), len(singles)))
        for name, key, unit in (("tilt", "tilt", "deg per mol%"),
                                ("contacts", "con", "beads per mol%"),
                                ("d_upper", "dup", "nm per mol%"),
                                ("d_lower", "dlo", "nm per mol%"),
                                ("S_up-lo", "ds", "per mol%")):
            y = np.array([r[key] for r in singles], float)
            if np.isnan(y).any():
                print("  %-9s not measured" % name)
                continue
            b, a0 = np.polyfit(x, y, 1)
            res = y - (a0 + b * x)
            s = float(np.sqrt(np.sum(res ** 2) / (len(x) - 2)))
            sb = s / np.sqrt(np.sum((x - x.mean()) ** 2))
            print("  %-9s slope %+10.5f +- %-9.5f %-15s residual sd %.4f  %s"
                  % (name, b, sb, unit, s,
                     "separated from zero" if abs(b) > 2 * sb else "not separated"))

    reps = {}
    for r in rows:
        if r["rep"] is not None:
            reps.setdefault(r["delta"], []).append(r)
    if reps:
        print("")
        print("spread between runs of one build, which is the yardstick")
        for delta, rs in sorted(reps.items()):
            if len(rs) < 2:
                continue
            for name, key in (("tilt", "tilt"), ("contacts", "con"),
                              ("d_upper", "dup"), ("S_up-lo", "ds")):
                v = np.array([r[key] for r in rs], float)
                if np.isnan(v).any():
                    continue
                print("  d%+d  %-9s %s   sd %.4f over %d runs"
                      % (delta, name, " ".join("%.3f" % t for t in v),
                         v.std(ddof=1), len(v)))
    else:
        print("")
        print("No replicate json was read. Pass --replicates to include")
        print("relax_d<delta>_s<n>.json, which is what makes the slopes above")
        print("readable: a slope matters only beside the spread between runs")
        print("of one build.")


if __name__ == "__main__":
    main()
