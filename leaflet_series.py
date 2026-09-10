#!/usr/bin/env python3
"""How many molecules of each species each leaflet holds, frame by frame.

leaflet_relax.py reads the whole trajectory into memory, which a trajectory of a
few microseconds does not fit into. This script streams it in chunks and reads
only the lipid beads, so the memory it needs does not grow with the length of
the run.

The leaflet of a lipid is read from the lipid: the head bead of a lipid in the
upper leaflet lies above the rest of that same molecule, and the head bead of a
lipid in the lower leaflet lies below it. No plane is drawn through the bilayer.

Usage:
  leaflet_series.py system.gro traj.xtc out.json [stride]
"""
import json
import sys

import numpy as np
import mdtraj as md

gro, xtc, out = sys.argv[1], sys.argv[2], sys.argv[3]
stride = int(sys.argv[4]) if len(sys.argv) > 4 else 1

STEROL_HEAD = {"ROH"}
PHOSPHATE = {"PO4"}

top = md.load(gro).topology
res = [r for r in top.residues
       if set(a.name for a in r.atoms) & (PHOSPHATE | STEROL_HEAD)]
if not res:
    sys.exit("no lipid found in %s" % gro)

lipid_atoms = sorted(a.index for r in res for a in r.atoms)
pos = {g: i for i, g in enumerate(lipid_atoms)}

head_mol, head_at, body_mol, body_at, po4_mol, po4_at = [], [], [], [], [], []
names = []
for m, r in enumerate(res):
    names.append(r.name.strip())
    h = [pos[a.index] for a in r.atoms if a.name in PHOSPHATE | STEROL_HEAD]
    b = [pos[a.index] for a in r.atoms if a.name not in PHOSPHATE | STEROL_HEAD]
    if not b:
        b = h
    head_mol += [m] * len(h); head_at += h
    body_mol += [m] * len(b); body_at += b
    p = [pos[a.index] for a in r.atoms if a.name in PHOSPHATE]
    po4_mol += [m] * len(p); po4_at += p

head_mol = np.array(head_mol); head_at = np.array(head_at)
body_mol = np.array(body_mol); body_at = np.array(body_at)
po4_mol = np.array(po4_mol); po4_at = np.array(po4_at)
nmol = len(res)
hn = np.bincount(head_mol, minlength=nmol).astype(float)
bn = np.bincount(body_mol, minlength=nmol).astype(float)
pn = np.bincount(po4_mol, minlength=nmol).astype(float)
has_po4 = pn > 0

namearr = np.array(names)
masks = {s: namearr == s for s in sorted(set(names))}

counts = {s: {"upper": [], "lower": []} for s in masks}
times, area, thick_u, thick_l = [], [], [], []

done = 0
for ch in md.iterload(xtc, top=gro, chunk=100, stride=stride,
                      atom_indices=lipid_atoms):
    z = ch.xyz[:, :, 2]
    for f in range(len(ch)):
        zf = z[f]
        zh = np.bincount(head_mol, weights=zf[head_at], minlength=nmol) / hn
        zb = np.bincount(body_mol, weights=zf[body_at], minlength=nmol) / bn
        up = zh > zb
        for s, m in masks.items():
            counts[s]["upper"].append(int(up[m].sum()))
            counts[s]["lower"].append(int((~up[m]).sum()))
        mid = 0.5 * (zh[up].mean() + zh[~up].mean())
        zp = np.bincount(po4_mol, weights=zf[po4_at], minlength=nmol)
        zp = np.divide(zp, pn, out=np.zeros(nmol), where=has_po4)
        thick_u.append(float(np.abs(zp[up & has_po4] - mid).mean()))
        thick_l.append(float(np.abs(zp[(~up) & has_po4] - mid).mean()))
    times += ch.time.tolist()
    area += (ch.unitcell_lengths[:, 0] * ch.unitcell_lengths[:, 1]).tolist()
    done += len(ch)
    print("  %d frames" % done, flush=True)

json.dump({"time_ps": times, "counts": counts, "area_nm2": area,
           "thickness_upper_nm": thick_u, "thickness_lower_nm": thick_l,
           "n_lipids": nmol, "stride": stride}, open(out, "w"))
print("written to", out)
