#!/usr/bin/env python3
"""
fix_vsites.py

Rebuild every virtual site in a built system.gro EXACTLY from its itp definition.

COBY (and most membrane builders) place coarse-grained lipids by their real
beads and leave virtual sites at an approximate, often in-plane, position. For
sterols the out-of-plane virtual_sites3 funct 4 sites (CHOL ROH, R3) then sit
~0.1 nm away from where GROMACS reconstructs them every step; in a packed
membrane that reconstructed position can land on a neighbouring bead and produce
an infinite Lennard-Jones force at minimisation. This script walks the system,
and for every molecule whose moleculetype declares virtual sites, overwrites the
vsite coordinates with the exact funct-based reconstruction from that molecule's
own real beads. Run it right after the build (before declashing and minimising).

Usage:
  fix_vsites.py --gro system.gro --top system.top [--itp-dir DIR]
"""

import argparse
import os
import re
import sys
import numpy as np

# reuse the exact constructor from itp_to_struct
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itp_to_struct import construct_vsite, _isnum   # noqa: E402


def parse_all_moleculetypes(itp_path):
    """Return two dicts from one itp:
      napm[molname]  = number of atoms in [ atoms ]
      vdefs[molname] = [(site, kind, funct, cons, params), ...] (may be empty)
    """
    napm = {}
    vdefs = {}
    mol = None
    n = 0
    vlist = []
    section = None

    def flush():
        if mol is not None:
            napm[mol] = n
            if vlist:
                vdefs[mol] = list(vlist)

    for raw in open(itp_path):
        line = raw.split(";")[0].rstrip()
        s = line.strip()
        if not s:
            continue
        m = re.match(r"\[\s*([a-zA-Z_0-9]+)\s*\]", s)
        if m:
            section = m.group(1).lower()
            if section == "moleculetype":
                flush()
                mol = None
                n = 0
                vlist = []
            continue
        f = s.split()
        if section == "moleculetype":
            mol = f[0]
        elif section == "atoms":
            n += 1
        elif section and section.startswith("virtual_sites"):
            site = int(f[0]) - 1
            if section == "virtual_sites2":
                ii = [int(t) - 1 for t in f[1:3]]
                funct = int(f[3]) if len(f) > 3 else 1
                params = [float(t) for t in f[4:] if _isnum(t)]
                vlist.append((site, "vs2", funct, ii, params))
            elif section == "virtual_sitesn":
                funct = int(f[1]) if len(f) > 1 and f[1].isdigit() else 1
                cons = [int(t) - 1 for t in f[2:] if t.lstrip("-").isdigit()]
                vlist.append((site, "vsn", funct, cons, []))
            else:
                ii = [int(t) - 1 for t in f[1:4]]
                funct = int(f[4]) if len(f) > 4 and f[4].lstrip("-").isdigit() else 1
                params = [float(t) for t in f[5:] if _isnum(t)]
                vlist.append((site, "vs3", funct, ii, params))
    flush()
    return napm, vdefs


def collect_itps(top_path, itp_dir):
    """All itp files reachable from the top (its includes) plus any in itp_dir."""
    paths = []
    base = os.path.dirname(os.path.abspath(top_path))
    for ln in open(top_path):
        m = re.match(r'\s*#include\s+"([^"]+)"', ln)
        if m:
            p = m.group(1)
            cand = p if os.path.isabs(p) else os.path.join(base, p)
            if os.path.exists(cand):
                paths.append(cand)
    if itp_dir and os.path.isdir(itp_dir):
        for fn in os.listdir(itp_dir):
            if fn.endswith(".itp"):
                paths.append(os.path.join(itp_dir, fn))
    return paths


def read_gro(path):
    L = open(path).read().splitlines()
    return L[0], int(L[1]), L[2:2 + int(L[1])], L[2 + int(L[1])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--itp-dir", default="")
    args = ap.parse_args()

    # gather atoms-per-molecule and vsite defs from all reachable itps
    napm = {}
    vdefs = {}
    for itp in collect_itps(args.top, args.itp_dir):
        try:
            a, v = parse_all_moleculetypes(itp)
            for k, val in a.items():
                napm.setdefault(k, val)
            for k, val in v.items():
                vdefs.setdefault(k, val)
        except Exception:
            continue
    if not vdefs:
        print("fix_vsites: no virtual-site moleculetypes found; nothing to do")
        return

    # molecule order/counts from [ molecules ]
    order = []
    in_mol = False
    for ln in open(args.top):
        s = ln.split(";")[0].strip()
        if s.lower().startswith("[ molecules"):
            in_mol = True
            continue
        if in_mol and s:
            if s.startswith("["):
                break
            p = s.split()
            order.append((p[0], int(p[1])))

    title, n, body, box = read_gro(args.gro)
    coords = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                       for b in body])

    idx = 0
    nfix = 0
    nmol_fixed = 0
    for molname, count in order:
        a = napm.get(molname)
        if a is None:                       # last resort: infer from gro resid run
            a = guess_napm(body, idx)
        if molname in vdefs:
            vlist = vdefs[molname]
            for _ in range(count):
                X = {k: coords[idx + k] for k in range(a)}
                for site, kind, funct, cons, params in vlist:
                    if site >= a or any(c >= a for c in cons):
                        continue            # malformed; skip safely
                    newpos = construct_vsite(kind, funct, cons, params, X)
                    coords[idx + site] = newpos
                    X[site] = newpos
                    nfix += 1
                idx += a
                nmol_fixed += 1
        else:
            idx += a * count

    if idx != n:
        sys.exit("fix_vsites: ERROR atom count mismatch (walked %d, gro has %d); "
                 "aborting without writing to avoid corruption" % (idx, n))

    # write back
    out = [title, str(n)]
    for i, b in enumerate(body):
        out.append(b[:20] + "%8.3f%8.3f%8.3f" % (coords[i, 0], coords[i, 1],
                                                 coords[i, 2]) + b[44:])
    out.append(box)
    open(args.gro, "w").write("\n".join(out) + "\n")
    print("fix_vsites: rebuilt %d virtual sites across %d molecules"
          % (nfix, nmol_fixed))


def guess_napm(body, start):
    """Atoms in the residue block beginning at index `start` (same resid run)."""
    if start >= len(body):
        return 1
    rid0 = body[start][0:5]
    k = start
    while k < len(body) and body[k][0:5] == rid0:
        k += 1
    return max(k - start, 1)


if __name__ == "__main__":
    main()
