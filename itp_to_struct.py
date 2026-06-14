#!/usr/bin/env python3
"""
itp_to_struct.py

Build a minimization-stable single-molecule coarse-grained structure (GRO) for a
lipid (or any single-residue Martini molecule) directly from its GROMACS .itp,
using:

  * the [atoms] order                         -> the bead list and names
  * [bonds] + [constraints] connectivity      -> the graph edges
  * the Martini ffbonded file                 -> the real length of each edge
                                                 (named bondtypes like
                                                 b_NC3_PO4_def carry no length in
                                                 the lipid itp; lengths live in
                                                 #define entries in ffbonded)

The bead coordinates are obtained by classical multidimensional scaling (MDS) of
the graph's all-pairs shortest-path distance matrix (weighted by the real bond
lengths). Unlike a straight-line layout, MDS reproduces rings (e.g. cholesterol's
fused sterol core) and bent chains, so the resulting structure already satisfies
the bonded distances closely and does NOT explode under LINCS when COBY packs it
into a membrane and GROMACS minimises it.

This keeps the tool fully general: any new lipidome is usable from its itp plus
the ffbonded definitions alone, with no pre-built membrane or external structure.

Usage:
  itp_to_struct.py --itp lipid.itp --mol DLPC --out DLPC.gro [--ffbonded ffb.itp]
  # prints: UPDOWN <mol> up:<headbead> upidx:<1-based> down:<tail1,tail2>
"""

import re
import argparse
import sys
from collections import defaultdict

import numpy as np

HEAD_TOKENS = {
    "NC3", "NH3", "PO4", "PO1", "PO2", "CNO", "GL0", "ROH",
    "AM1", "AM2", "OH1", "OH2", "TP1",
    "GM", "B1", "B2", "B3", "S1", "S2", "S3",
}
DEFAULT_LEN = 0.47   # nm, fallback when an edge length cannot be resolved


def parse_ffbonded(path):
    """Map named bond/constraint type -> length (nm) from #define entries."""
    lengths = {}
    if not path:
        return lengths
    try:
        with open(path) as fh:
            for raw in fh:
                line = raw.split(";", 1)[0].split()
                if len(line) >= 4 and line[0] == "#define":
                    try:
                        lengths[line[1]] = float(line[3])
                    except ValueError:
                        pass
    except OSError:
        pass
    return lengths


def _isnum(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def construct_vsite(kind, funct, cons, params, X):
    """Exact virtual-site position from its GROMACS definition.

    virtual_sites3 funct 1 (3):    r = ri + a(rj-ri) + b(rk-ri)
    virtual_sites3 funct 4 (3out): r = ri + a(rj-ri) + b(rk-ri)
                                       + c[(rj-ri) x (rk-ri)]
    virtual_sites2 funct 1:        r = ri + a(rj-ri)
    virtual_sitesn / unknown:      centroid of constructing atoms
    """
    P = [X[c] for c in cons]
    if kind == "vs3" and len(P) == 3:
        ri, rj, rk = P
        rij, rik = rj - ri, rk - ri
        if funct == 4 and len(params) >= 3:
            a, b, c = params[:3]
            return ri + a * rij + b * rik + c * np.cross(rij, rik)
        if funct == 1 and len(params) >= 2:
            a, b = params[:2]
            return ri + a * rij + b * rik
        if funct == 2 and len(params) >= 2:            # 3fd
            a, d = params[:2]
            base = ri + a * rij
            dirv = rk - base
            nrm = np.linalg.norm(dirv)
            return base + d * dirv / nrm if nrm > 1e-9 else base
    if kind == "vs2" and len(P) == 2 and params:
        ri, rj = P
        return ri + params[0] * (rj - ri)
    return np.mean(P, axis=0)


def parse_itp(path, mol):
    """Return (bead_names, edges); edges = (i, j, typename|None, inlen|None).
    Bonds and constraints become real edges (their length drives the geometry).
    Virtual sites become short nominal edges to their constructing atoms so the
    vsite beads cluster sensibly; GROMACS reconstructs vsites every step, so only
    their rough initial placement matters."""
    atoms, edges = [], []
    vsites = {}
    vsite_defs = {}
    section = None
    in_target = False
    with open(path) as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].rstrip()
            if not line.strip():
                continue
            s = line.strip()
            if s.startswith("["):
                section = s.strip("[] ").lower()
                continue
            if s.startswith("#"):
                continue  # skip #ifdef/#else/#endif/#define inside molecule
            if section == "moleculetype":
                in_target = (s.split()[0] == mol)
                continue
            if not in_target:
                continue
            if section == "atoms":
                atoms.append(s.split()[4])
            elif section in ("bonds", "constraints"):
                f = s.split()
                if len(f) >= 2 and f[0].lstrip("-").isdigit() \
                        and f[1].lstrip("-").isdigit():
                    i, j = int(f[0]) - 1, int(f[1]) - 1
                    tname = None
                    if len(f) >= 3 and not f[2].lstrip("-").replace(".", "", 1).isdigit():
                        tname = f[2]
                    inlen = None
                    if len(f) >= 4:
                        try:
                            inlen = float(f[3])
                        except ValueError:
                            inlen = None
                    edges.append((i, j, tname, inlen))
            elif section in ("virtual_sites2", "virtual_sites3",
                             "virtual_sitesn", "virtual_sites4"):
                f = s.split()
                if not f or not f[0].isdigit():
                    continue
                site = int(f[0]) - 1
                if section == "virtual_sitesn":
                    # site funct a1 a2 ... (COG/COM of listed atoms)
                    funct = int(f[1]) if len(f) > 1 and f[1].isdigit() else 1
                    cons = [int(t) - 1 for t in f[2:] if t.lstrip("-").isdigit()]
                    vsite_defs[site] = ("vsn", funct, cons, [])
                    vsites[site] = [c for c in cons if c != site]
                elif section == "virtual_sites2":
                    # site i j funct a
                    ii = [int(t) - 1 for t in f[1:3]]
                    funct = int(f[3]) if len(f) > 3 else 1
                    params = [float(t) for t in f[4:] if _isnum(t)]
                    vsite_defs[site] = ("vs2", funct, ii, params)
                    vsites[site] = [c for c in ii if c != site]
                else:  # virtual_sites3 / 4
                    ii = [int(t) - 1 for t in f[1:4]]
                    funct = int(f[4]) if len(f) > 4 and f[4].lstrip("-").isdigit() else 1
                    params = [float(t) for t in f[5:] if _isnum(t)]
                    vsite_defs[site] = ("vs3", funct, ii, params)
                    vsites[site] = [c for c in ii if c != site]
    if not atoms:
        sys.exit("ERROR: moleculetype '%s' not found (or has no [atoms]) in %s"
                 % (mol, path))
    return atoms, edges, vsites, vsite_defs


def edge_length(tname, inlen, ffb):
    if inlen is not None:
        return inlen
    if tname is not None and tname in ffb:
        return ffb[tname]
    return DEFAULT_LEN


def all_pairs_dist(n, adj):
    """Floyd-Warshall shortest paths on the weighted bonded graph (small n)."""
    INF = 1e9
    D = np.full((n, n), INF)
    np.fill_diagonal(D, 0.0)
    for i, nbrs in adj.items():
        for j, w in nbrs:
            if w < D[i, j]:
                D[i, j] = D[j, i] = w
    for k in range(n):
        D = np.minimum(D, D[:, k][:, None] + D[k, :][None, :])
    finite_max = D[D < INF].max() if np.any(D < INF) else 1.0
    D[D >= INF] = finite_max * 2.0 + 1.0
    return D


def adj_full(edges, n):
    """Adjacency (unweighted) over all beads for head/tail detection."""
    a = defaultdict(list)
    for (i, j, _t, _l) in edges:
        if i != j:
            a[i].append(j)
            a[j].append(i)
    for i in range(n):
        a.setdefault(i, [])
    return a


def relax_bonds(X, bonds, iters=400):
    """Jakobsen-style iterative constraint relaxation: nudge every bonded pair
    toward its exact target length. Starting from the MDS layout (good global
    shape, rings preserved) this drives all bonded distances to the correct
    values, giving a minimization-stable structure."""
    if not bonds:
        return X
    rng = np.random.default_rng(0)
    for _ in range(iters):
        for (i, j, L) in bonds:
            d = X[j] - X[i]
            dist = float(np.linalg.norm(d))
            if dist < 1e-6:
                d = rng.standard_normal(3) * 1e-3
                dist = float(np.linalg.norm(d))
            corr = 0.5 * (dist - L) / dist
            X[i] = X[i] + corr * d
            X[j] = X[j] - corr * d
    return X


def mds_3d(D):
    n = D.shape[0]
    if n == 1:
        return np.zeros((1, 3))
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    B = (B + B.T) / 2.0
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    k = min(3, n)
    L = np.clip(vals[:k], 0, None)
    X = vecs[:, :k] * np.sqrt(L)[None, :]
    if X.shape[1] < 3:
        X = np.hstack([X, np.zeros((n, 3 - X.shape[1]))])
    return X


def pick_root(names, adj):
    for i, nm in enumerate(names):
        if nm.upper() in HEAD_TOKENS:
            return i
    for i in range(len(names)):
        if len(adj[i]) == 1:
            return i
    return 0


def _is_tail_name(nm):
    """A bead name that denotes an aliphatic acyl or sphingoid chain bead, e.g.
    C1A, C2A, C4B, D2A (one double bond), T1A. Sugar and head beads do not match.
    """
    return bool(re.match(r"^[CDT]\d+[AB]?$", nm.upper()))


def pick_head_root(names, adj):
    """Root the orientation at the head bead farthest from the aliphatic tails.

    Tail terminals are the degree-one beads whose names are aliphatic. Head
    candidates are the head-token beads, or, if none are named, the non-aliphatic
    degree-one beads. The chosen root is the head candidate whose nearest tail
    terminal is the most bonds away, which is the outermost head bead (a sugar
    tip for a ganglioside, the choline or sterol hydroxyl for a simple lipid).
    """
    n = len(names)
    tail_terms = [i for i in range(n) if len(adj[i]) == 1 and _is_tail_name(names[i])]
    heads = [i for i, nm in enumerate(names) if nm.upper() in HEAD_TOKENS]
    if not heads:
        heads = [i for i in range(n) if len(adj[i]) == 1 and not _is_tail_name(names[i])]
    if not heads or not tail_terms:
        return pick_root(names, adj)

    def bfs(src):
        d = [-1] * n
        d[src] = 0
        q = [src]
        while q:
            cur = q.pop(0)
            for nb in adj[cur]:
                if d[nb] < 0:
                    d[nb] = d[cur] + 1
                    q.append(nb)
        return d
    dist_to_tail = [min(bfs(t)[h] for t in tail_terms) for h in heads]
    return heads[int(np.argmax(dist_to_tail))]


def orient_head_up(X, root):
    X = X - X.mean(axis=0)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    axis = Vt[0]
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z)
    s = np.linalg.norm(v)
    c = float(np.dot(axis, z))
    if s > 1e-8:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))
        X = X @ R.T
    if X[root, 2] < X[:, 2].mean():
        X[:, 2] = -X[:, 2]
    return X - X.mean(axis=0)


def tail_ends(names, X, root, adj):
    terms = [(X[i, 2], names[i]) for i in range(len(names))
             if len(adj[i]) == 1 and i != root]
    terms.sort()
    return [nm for _, nm in terms[:2]] if terms else [names[-1]]


def write_gro(path, mol, names, X):
    n = len(names)
    with open(path, "w") as fh:
        fh.write("single-molecule CG structure from itp connectivity (MDS)\n")
        fh.write("%5d\n" % n)
        for i, (nm, xyz) in enumerate(zip(names, X), start=1):
            fh.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"
                     % (1, mol[:5], nm[:5], i, xyz[0], xyz[1], xyz[2]))
        fh.write("%10.5f%10.5f%10.5f\n" % (5.0, 5.0, 5.0))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--itp", required=True)
    p.add_argument("--mol", required=True, help="moleculetype name (e.g. DLPC)")
    p.add_argument("--out", required=True)
    p.add_argument("--ffbonded", default="",
                   help="Martini ffbonded itp providing #define bond lengths")
    args = p.parse_args()

    ffb = parse_ffbonded(args.ffbonded)
    names, edges, vsites, vsite_defs = parse_itp(args.itp, args.mol)
    n = len(names)
    vset = set(vsites)
    real = [i for i in range(n) if i not in vset]
    rindex = {g: r for r, g in enumerate(real)}

    # MDS over the real (non-virtual) beads, driven by real bond/constraint
    # lengths. Virtual sites are excluded here and positioned afterwards, since
    # their "distance" to constructors is not a real bonded length.
    adj = defaultdict(list)
    for (i, j, tname, inlen) in edges:
        if i == j or i in vset or j in vset:
            continue
        w = edge_length(tname, inlen, ffb)
        adj[rindex[i]].append((rindex[j], w))
        adj[rindex[j]].append((rindex[i], w))
    for r in range(len(real)):
        adj.setdefault(r, [])

    X = np.zeros((n, 3))
    if real:
        D = all_pairs_dist(len(real), adj)
        Xr = mds_3d(D)
        # refine: drive every real bonded distance to its exact target length
        rbonds = []
        for (i, j, tname, inlen) in edges:
            if i == j or i in vset or j in vset:
                continue
            rbonds.append((rindex[i], rindex[j], edge_length(tname, inlen, ffb)))
        Xr = relax_bonds(Xr, rbonds)
        for r, g in enumerate(real):
            X[g] = Xr[r]

    # Orient as a proper lipid: head at the top, ALL tails hanging below it.
    # The MDS/graph embedding alone tends to put the two acyl tails at opposite
    # ends of the long axis with the head branching off the middle; aligning that
    # principal axis to z then leaves one tail pointing up into the water, so COBY
    # builds an inverted, broken bilayer. Instead we set z from the graph-distance
    # (number of bonds) to the head bead: the head sits at the top and every chain
    # descends from it, guaranteeing both tails point the same way (down). The
    # lateral (x,y) spread is kept from the embedding so the two tails do not
    # overlap; equilibration relaxes the rest.
    adj = adj_full(edges, n)
    # virtual sites (e.g. sterol ROH) carry no bonds, so add links from each
    # vsite to its constructing atoms; otherwise BFS from a vsite head reaches
    # nothing and the depth (hence the orientation) is degenerate.
    for site, cons in vsites.items():
        for c in cons:
            if c not in adj[site]:
                adj[site].append(c)
            if site not in adj[c]:
                adj[c].append(site)
    # Choose the head root as the head bead farthest from the aliphatic tails.
    # Picking the first head-token bead fails for a ganglioside: its branched
    # sugar head is as many bonds deep as the acyl tails, so the depth from an
    # inner sugar makes a sialic acid tip rank among the deepest beads and that
    # tip is then placed at the bottom, into the membrane. Rooting at the
    # outermost head bead (the sugar tip) makes the tails the deepest beads, so
    # they go down and the whole sugar head goes up. For a simple lipid the
    # outermost head bead is still the choline or the sterol hydroxyl, so this
    # does not change those cases.
    root = pick_head_root(names, adj)
    depth = [-1] * n
    depth[root] = 0
    queue = [root]
    while queue:
        cur = queue.pop(0)
        for nb in adj[cur]:
            if depth[nb] < 0:
                depth[nb] = depth[cur] + 1
                queue.append(nb)
    maxd = max(d for d in depth if d >= 0)
    for i in range(n):
        if depth[i] < 0:          # disconnected bead: drop it to the bottom
            depth[i] = maxd + 1
    X = X - X.mean(axis=0)
    # Make the head->tail direction the LONG (principal) axis so COBY's
    # principal-axis alignment stands the lipid upright. Two acyl tails make the
    # lateral spread the largest-variance axis otherwise, and COBY then lays the
    # lipid on its side (tails poke into the water). We stretch z by graph depth
    # and compress the lateral coordinates so the z-extent dominates.
    layer = 0.40                   # nm per bond along the membrane normal
    z_from_depth = -np.array(depth, dtype=float) * layer
    z_from_depth -= z_from_depth.mean()
    # compress x,y so their total spread is well under the z spread
    xy = X[:, :2].copy()
    xy -= xy.mean(axis=0)
    z_extent = max(z_from_depth.max() - z_from_depth.min(), 1e-3)
    xy_extent = max(np.linalg.norm(xy, axis=1).max() * 2.0, 1e-3)
    target_xy = 0.35 * z_extent    # keep tails separated but laterally narrow
    if xy_extent > target_xy:
        xy *= target_xy / xy_extent
    X[:, :2] = xy
    X[:, 2] = z_from_depth

    # the depth layout fixes orientation but slightly distorts fused rings
    # (e.g. the sterol ring); relax all real bond/constraint lengths back to
    # target on top of the head-up layout (local nudging keeps the orientation).
    fbonds = []
    for (i, j, tname, inlen) in edges:
        if i == j or i in vset or j in vset:
            continue
        fbonds.append((i, j, edge_length(tname, inlen, ffb)))
    if fbonds:
        X = relax_bonds(X, fbonds)

    # Re-assert the upright orientation with a rigid rotation (preserves the
    # relaxed geometry, rings and bond lengths exactly): align the vector from
    # the head to the deepest tail beads with -z, so the head-to-tail direction
    # is the membrane normal. This makes COBY's principal-axis alignment stand
    # every lipid up (phospholipids and the flat sterol ring alike) instead of
    # laying it on its side.
    head_ids = [i for i in range(n) if names[i].upper() in HEAD_TOKENS] or [root]
    deepest = max(depth)
    tail_ids = [i for i in range(n) if depth[i] >= deepest - 0]
    if not tail_ids:
        tail_ids = [int(np.argmax(depth))]
    head_ref = X[head_ids].mean(axis=0)
    tail_ref = X[tail_ids].mean(axis=0)
    v = head_ref - tail_ref            # should point head-up
    nv = np.linalg.norm(v)
    if nv > 1e-6:
        v = v / nv
        zaxis = np.array([0.0, 0.0, 1.0])
        axis = np.cross(v, zaxis)
        s = np.linalg.norm(axis)
        c = float(np.dot(v, zaxis))
        if s > 1e-8:
            axis = axis / s
            K = np.array([[0, -axis[2], axis[1]],
                          [axis[2], 0, -axis[0]],
                          [-axis[1], axis[0], 0]])
            R = np.eye(3) + s * K + (1 - c) * (K @ K)
            X = X @ R.T
        elif c < 0:                    # already anti-aligned: flip
            X[:, 2] = -X[:, 2]
    X = X - X.mean(axis=0)

    # construct each virtual site EXACTLY from its itp definition (now from the
    # final, possibly-reflected real coordinates), so the initial coordinates
    # already match what GROMACS reconstructs every step. This is what keeps
    # sterol ring vsites from landing on a real bead and producing an infinite
    # Lennard-Jones force at minimisation.
    # Enforce a minimum spacing between non-bonded real beads. The upright layout
    # plus lateral compression can place beads from different branches (e.g. the
    # two acyl-tail roots C1A/C1B) on top of each other; COBY then copies that
    # overlap into every lipid and minimisation sees an infinite force. Push any
    # too-close non-bonded pair apart and restore bond lengths, a few rounds.
    bonded_pairs = set()
    for (i, j, _t, _l) in edges:
        if i in vset or j in vset:
            continue
        bonded_pairs.add((min(i, j), max(i, j)))
    dmin = 0.30
    for _ in range(60):
        moved = False
        for a in range(len(real)):
            for b in range(a + 1, len(real)):
                i, j = real[a], real[b]
                if (i, j) in bonded_pairs:
                    continue
                d = X[j] - X[i]
                r = float(np.linalg.norm(d))
                if r < dmin:
                    if r < 1e-6:
                        d = np.random.default_rng(i * 1000 + j).normal(size=3)
                        r = float(np.linalg.norm(d))
                    push = 0.5 * (dmin - r) / r * d
                    X[i] -= push
                    X[j] += push
                    moved = True
        if fbonds:
            X = relax_bonds(X, fbonds)
        if not moved:
            break

    placed = set(real)
    for _ in range(6):
        progressed = False
        for site, (kind, funct, cons, params) in vsite_defs.items():
            if site in placed:
                continue
            if not all(c in placed for c in cons):
                continue
            X[site] = construct_vsite(kind, funct, cons, params, X)
            placed.add(site)
            progressed = True
        if not progressed:
            break

    out = args.out if args.out.endswith(".gro") else args.out + ".gro"
    write_gro(out, mol=args.mol, names=names, X=X)

    downs = tail_ends(names, X, root, adj_full(edges, n))
    n_unresolved = sum(1 for (i, j, t, il) in edges
                       if i not in vset and j not in vset
                       and il is None and (t is None or t not in ffb))
    tail0 = int(np.argmax(depth))
    print("UPDOWN %s up:%s upidx:%d down:%s head0:%d tail0:%d"
          % (args.mol, names[root], root + 1, ",".join(downs), root, tail0))
    print("itp_to_struct: %s -> %s (%d beads, %d vsites placed, %d unresolved"
          "->default) via MDS"
          % (args.mol, out, n, len(vsites), n_unresolved))


if __name__ == "__main__":
    main()
