#!/usr/bin/env python3
"""Cholesterol around the protein, in each leaflet separately."""
import json, sys
import numpy as np
import mdtraj as md

gro, xtc, out = sys.argv[1], sys.argv[2], sys.argv[3]
stride = int(sys.argv[4]) if len(sys.argv) > 4 else 4
t0 = float(sys.argv[5]) if len(sys.argv) > 5 else 1e6

HEAD = {"PO4", "ROH"}
SOL = {"W", "WF", "NA", "CL", "ION"}
top = md.load(gro).topology
lip = [r for r in top.residues if set(a.name for a in r.atoms) & HEAD]
pro = [r for r in top.residues if r.name.strip() not in SOL
       and not r.name.strip().startswith("W")
       and not (set(a.name for a in r.atoms) & HEAD)]
lip_atoms = sorted(a.index for r in lip for a in r.atoms)
pro_atoms = sorted(a.index for r in pro for a in r.atoms)
keep = sorted(set(lip_atoms) | set(pro_atoms))
pos = {g: i for i, g in enumerate(keep)}
pro_sub = np.array([pos[i] for i in pro_atoms])

names, hm, ha, bm, ba, allm, alla, first = [], [], [], [], [], [], [], []
for m, r in enumerate(lip):
    names.append(r.name.strip())
    idx = [pos[a.index] for a in r.atoms]
    h = [pos[a.index] for a in r.atoms if a.name in HEAD]
    b = [pos[a.index] for a in r.atoms if a.name not in HEAD] or h
    hm += [m]*len(h); ha += h; bm += [m]*len(b); ba += b
    allm += [m]*len(idx); alla += idx; first.append(idx[0])
hm=np.array(hm); ha=np.array(ha); bm=np.array(bm); ba=np.array(ba)
allm=np.array(allm); alla=np.array(alla); first=np.array(first)
n=len(lip)
hn=np.bincount(hm,minlength=n).astype(float); bn=np.bincount(bm,minlength=n).astype(float)
an=np.bincount(allm,minlength=n).astype(float)
namearr=np.array(names)
species=sorted(set(names))

E=np.array([0,0.8,1.2,1.6,2.0,2.5,3.0,4.0,5.0,7.0,100.0]); nb=len(E)-1
cnt={"upper":np.zeros(nb),"lower":np.zeros(nb)}
sp={x:{"upper":np.zeros(nb),"lower":np.zeros(nb)} for x in species}
masks={x:(namearr==x) for x in species}
nfr=0
for ch in md.iterload(xtc, top=gro, chunk=50, stride=stride, atom_indices=keep):
    for f in range(len(ch)):
        if ch.time[f] < t0: continue
        xyz=ch.xyz[f]; box=ch.unitcell_lengths[f][:2]; z=xyz[:,2]
        zh=np.bincount(hm,weights=z[ha],minlength=n)/hn
        zb=np.bincount(bm,weights=z[ba],minlength=n)/bn
        up=zh>zb
        mid=0.5*(zh[up].mean()+zh[~up].mean())
        o=xyz[first,:2]; d=xyz[alla,:2]-o[allm]; d-=box*np.round(d/box)
        L=np.stack([o[:,0]+np.bincount(allm,weights=d[:,0],minlength=n)/an,
                    o[:,1]+np.bincount(allm,weights=d[:,1],minlength=n)/an],axis=1)
        tm=pro_sub[np.abs(z[pro_sub]-mid)<1.5]
        if len(tm)==0: continue
        P=xyz[tm,:2]
        dd=L[:,None,:]-P[None,:,:]; dd-=box*np.round(dd/box)
        r=np.sqrt((dd*dd).sum(axis=2)).min(axis=1)
        i=np.clip(np.digitize(r,E)-1,0,nb-1)
        for side,sel in (("upper",up),("lower",~up)):
            np.add.at(cnt[side],i[sel],1.0)
            for x in species:
                np.add.at(sp[x][side],i[sel&masks[x]],1.0)
        nfr+=1
json.dump({"edges":E.tolist(),"n_frames":nfr,"n_tm_beads":int(len(tm)),
           "count":{k:v.tolist() for k,v in cnt.items()},
           "species":{x:{k:v.tolist() for k,v in sp[x].items()} for x in species}},open(out,"w"))
print("frames",nfr,"->",out)
