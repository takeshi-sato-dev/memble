"""Is the upper leaflet demixed on its own, or only around the protein?

For every lipid of one leaflet, the composition of the neighbours within 1.5 nm
is taken. A leaflet that is mixed gives every lipid the composition of the
leaflet. A leaflet that has separated gives a DLPC molecule more DLPC around it
than the leaflet holds. The lipids are split by distance from the protein, so a
perturbation that belongs to the protein is separated from one that does not.
"""
import json, sys
import numpy as np
import mdtraj as md

gro, xtc, out = sys.argv[1], sys.argv[2], sys.argv[3]
stride = int(sys.argv[4]); t0 = float(sys.argv[5]); RC = 1.5

HEAD={"PO4","ROH"}; SOL={"W","WF","NA","CL","ION"}
top=md.load(gro).topology
lip=[r for r in top.residues if set(a.name for a in r.atoms)&HEAD]
pro=[r for r in top.residues if r.name.strip() not in SOL and not r.name.strip().startswith("W")
     and not (set(a.name for a in r.atoms)&HEAD)]
lip_atoms=sorted(a.index for r in lip for a in r.atoms)
pro_atoms=sorted(a.index for r in pro for a in r.atoms)
keep=sorted(set(lip_atoms)|set(pro_atoms)); pos={g:i for i,g in enumerate(keep)}
pro_sub=np.array([pos[i] for i in pro_atoms])
names=[];hm=[];ha=[];bm=[];ba=[];allm=[];alla=[];first=[]
for m,r in enumerate(lip):
    names.append(r.name.strip()); idx=[pos[a.index] for a in r.atoms]
    h=[pos[a.index] for a in r.atoms if a.name in HEAD]
    b=[pos[a.index] for a in r.atoms if a.name not in HEAD] or h
    hm+=[m]*len(h); ha+=h; bm+=[m]*len(b); ba+=b
    allm+=[m]*len(idx); alla+=idx; first.append(idx[0])
hm,ha,bm,ba=map(np.array,(hm,ha,bm,ba)); allm=np.array(allm); alla=np.array(alla); first=np.array(first)
n=len(lip); hn=np.bincount(hm,minlength=n).astype(float); bn=np.bincount(bm,minlength=n).astype(float)
an=np.bincount(allm,minlength=n).astype(float)
namearr=np.array(names); species=sorted(set(names))

BANDS=[(0,3.0),(3.0,7.0),(7.0,1e3)]
acc={s:{b:{"self":0.0,"nb":0.0,"n":0.0} for b in range(3)} for s in species}
leafmol={s:0.0 for s in species}; nlip=0.0; nfr=0
snap=None
for ch in md.iterload(xtc, top=gro, chunk=50, stride=stride, atom_indices=keep):
    for f in range(len(ch)):
        if ch.time[f]<t0: continue
        xyz=ch.xyz[f]; box=ch.unitcell_lengths[f][:2]; z=xyz[:,2]
        zh=np.bincount(hm,weights=z[ha],minlength=n)/hn
        zb=np.bincount(bm,weights=z[ba],minlength=n)/bn
        up=zh>zb; mid=0.5*(zh[up].mean()+zh[~up].mean())
        o=xyz[first,:2]; d=xyz[alla,:2]-o[allm]; d-=box*np.round(d/box)
        L=np.stack([o[:,0]+np.bincount(allm,weights=d[:,0],minlength=n)/an,
                    o[:,1]+np.bincount(allm,weights=d[:,1],minlength=n)/an],axis=1)
        tm=pro_sub[np.abs(z[pro_sub]-mid)<1.5]
        if len(tm)==0: continue
        P=xyz[tm,:2]
        dd=L[:,None,:]-P[None,:,:]; dd-=box*np.round(dd/box)
        rp=np.sqrt((dd*dd).sum(axis=2)).min(axis=1)
        U=np.where(up)[0]
        XY=L[U]; nm=namearr[U]
        D=XY[:,None,:]-XY[None,:,:]; D-=box*np.round(D/box)
        R=np.sqrt((D*D).sum(axis=2)); np.fill_diagonal(R,1e9)
        near=R<RC
        for s in species:
            leafmol[s]+= (nm==s).sum()
        nlip+=len(U); nfr+=1
        for bi,(a,b) in enumerate(BANDS):
            sel=(rp[U]>=a)&(rp[U]<b)
            for s in species:
                m=sel&(nm==s)
                if not m.any(): continue
                cnt=near[m].sum(axis=1).astype(float)
                same=(near[m]&(nm==s)[None,:]).sum(axis=1).astype(float)
                ok=cnt>0
                acc[s][bi]["self"]+=same[ok].sum(); acc[s][bi]["nb"]+=cnt[ok].sum()
                acc[s][bi]["n"]+=ok.sum()
        if snap is None:
            snap={"xy":XY.tolist(),"name":nm.tolist(),"box":box.tolist(),
                  "rp":rp[U].tolist(),"time_ps":float(ch.time[f])}
res={"n_frames":nfr,"leaflet_molpct":{s:100*leafmol[s]/nlip for s in species},
     "bands":[list(b) for b in BANDS],
     "local":{s:[(100*acc[s][b]["self"]/acc[s][b]["nb"] if acc[s][b]["nb"] else None) for b in range(3)] for s in species},
     "n_scored":{s:[acc[s][b]["n"] for b in range(3)] for s in species},
     "snapshot":snap}
json.dump(res,open(out,"w")); print("frames",nfr,"->",out)
