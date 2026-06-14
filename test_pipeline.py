import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
"""
Automated correctness tests for the helper scripts.

Each helper is exercised as it is used in the pipeline (as a command-line tool)
on small synthetic inputs, and the output or exit status is checked. Run with:

    python3 -m pytest tests/ -v

The helper directory defaults to the repository root (the parent of this file)
and can be overridden with the HELPERS_DIR environment variable.
"""

import os
import subprocess
import sys
import math

HELPERS = os.environ.get(
    "HELPERS_DIR", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)


def run(script, *args, expect_zero=True):
    cmd = [sys.executable, os.path.join(HELPERS, script), *map(str, args)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if expect_zero:
        assert p.returncode == 0, "%s failed: %s" % (script, p.stderr)
    return p


# ---------------------------------------------------------------- orient_tm
def test_orient_tm_aligns_to_z(tmp_path):
    """A tilted hydrophobic helix is oriented to near-zero tilt from z."""
    import numpy as np
    axis = np.array([1.0, 1.0, 2.0]); axis /= np.linalg.norm(axis)
    p = tmp_path / "aa.pdb"
    lines = []
    for r in range(1, 41):
        resn = "LEU" if 10 <= r <= 30 else "LYS"
        pos = np.array([5.0, 5.0, 5.0]) + axis * 1.5 * (r - 10)
        lines.append("ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C"
                     % (r, resn, r, pos[0], pos[1], pos[2]))
    p.write_text("\n".join(lines) + "\nEND\n")
    out = tmp_path / "ori.pdb"
    res = run("orient_tm.py", "--in", p, "--out", out)
    assert "residual tilt" in res.stdout
    tilt = float(res.stdout.split("residual tilt")[1].split("deg")[0])
    assert tilt < 1.0, "TM not aligned to z (tilt %.2f)" % tilt


# ------------------------------------------------------- replicate (grid)
def test_replicate_grid_count_and_separation(tmp_path):
    """Grid placement yields N copies, well separated inside the box."""
    src = tmp_path / "one.pdb"
    src.write_text(
        "ATOM      1 BB   ALA A   1      10.000  10.000  20.000  1.00  0.00\n"
        "ATOM      2 BB   ALA A   2      12.000  10.000  20.000  1.00  0.00\nEND\n")
    out = tmp_path / "four.pdb"
    run("replicate_and_fix_top.py", "replicate", "--in", src, "--out", out,
        "--grid", 4, "--box", "40,40", "--margin", 8)
    chains = {}
    for ln in out.read_text().splitlines():
        if ln.startswith("ATOM"):
            chains.setdefault(ln[21], []).append((float(ln[30:38]), float(ln[38:46])))
    assert len(chains) == 4, "expected 4 copies, got %d" % len(chains)
    centers = [(_avg(c, 0), _avg(c, 1)) for c in chains.values()]
    dmin = min(math.dist(a, b) for i, a in enumerate(centers)
              for b in centers[i + 1:])
    assert dmin > 80.0, "copies too close (min %.1f A)" % dmin  # >8 nm


def _avg(pts, i):
    return sum(p[i] for p in pts) / len(pts)


# ---------------------------------------------------- inject_posres (auto)
def test_inject_posres_auto_head(tmp_path):
    """--beads auto picks ROH for a sterol and PO4 for a phospholipid."""
    itp = tmp_path / "lip.itp"
    itp.write_text(
        "[ moleculetype ]\nCHOL 1\n[ atoms ]\n1 SP1 1 CHOL ROH 1 0\n2 SC1 1 CHOL R1 2 0\n\n"
        "[ moleculetype ]\nDPSM 1\n[ atoms ]\n1 Q1 1 DPSM NC3 1 1\n2 Q5 1 DPSM PO4 2 -1\n")
    run("inject_posres.py", "--itp", itp, "--mol", "CHOL", "--beads", "auto",
        "--stage", "POSRES_STEP1:400")
    run("inject_posres.py", "--itp", itp, "--mol", "DPSM", "--beads", "auto",
        "--stage", "POSRES_STEP1:400")
    text = itp.read_text()
    # two guarded blocks injected: CHOL on its ROH (local idx 1),
    # DPSM on its PO4 (local idx 2)
    assert text.count("#ifdef POSRES_STEP1") == 2
    assert "     1 1 400 400 400" in text   # ROH is atom 1 of CHOL
    assert "     2 1 400 400 400" in text   # PO4 is atom 2 of DPSM


def test_inject_posres_missing_molecule_errors(tmp_path):
    itp = tmp_path / "x.itp"
    itp.write_text("[ moleculetype ]\nCHOL 1\n[ atoms ]\n1 SP1 1 CHOL ROH 1 0\n")
    p = run("inject_posres.py", "--itp", itp, "--mol", "NOPE", "--beads", "auto",
            "--stage", "POSRES_STEP1:400", expect_zero=False)
    assert p.returncode != 0


# ------------------------------------------------ leaflet_area_check
def _write_gro(path, recs, box="10 10 10"):
    lines = ["title", str(len(recs))]
    for (resid, resname, aname, z) in recs:
        lines.append("%5s%-5s%5s%5d%8.3f%8.3f%8.3f"
                     % ("%5d" % resid, resname, aname, 1, 1.0, 1.0, z))
    lines.append(box)
    path.write_text("\n".join(lines) + "\n")


def test_leaflet_area_balanced_passes(tmp_path):
    gro = tmp_path / "sym.gro"
    recs = [(r, "POPC", "PO4", +2.0) for r in range(1, 101)] + \
           [(r, "POPC", "PO4", -2.0) for r in range(101, 201)]
    _write_gro(gro, recs)
    run("leaflet_area_check.py", "--gro", gro, "--lipids", "POPC", "--asym", 1)


def test_leaflet_area_mismatch_aborts(tmp_path):
    # severe mismatch (> hard-tol): upper all POPC (0.64), lower all CHOL (0.40)
    gro = tmp_path / "asym.gro"
    recs = [(r, "POPC", "PO4", +2.0) for r in range(1, 101)] + \
           [(r, "CHOL", "ROH", -2.0) for r in range(101, 201)]
    _write_gro(gro, recs)
    p = run("leaflet_area_check.py", "--gro", gro, "--lipids", "POPC CHOL",
            "--asym", 1, expect_zero=False)
    assert p.returncode != 0
    assert "MISMATCH" in (p.stdout + p.stderr)


def test_leaflet_area_moderate_warns_not_abort(tmp_path):
    # moderate mismatch (tol < dev < hard-tol): warn but continue so the build
    # still produces run files; the user balances or lets equilibration absorb it
    gro = tmp_path / "asym2.gro"
    recs = [(r, "POPC", "PO4", +2.0) for r in range(1, 101)] + \
           [(r, "POPC", "PO4", -2.0) for r in range(101, 151)] + \
           [(r, "CHOL", "ROH", -2.0) for r in range(151, 201)]
    _write_gro(gro, recs)
    p = run("leaflet_area_check.py", "--gro", gro, "--lipids", "POPC CHOL",
            "--asym", 1, expect_zero=False)
    assert p.returncode == 0
    assert "WARN" in (p.stdout + p.stderr)


# ------------------------------------------------ place_partner (post-COBY)
def _mini_membrane_system(tmp_path):
    """A tiny built system: TM protein, two DLPC leaflets, water, ions."""
    import numpy as np
    atoms = []
    rid = 1
    for k in range(10):
        atoms.append((rid, "PRO", "BB", 3.0, 3.0, 0.5 + 8.0 * k / 9)); 
    rid += 1
    for zc in (3.8, 5.2):
        for i in range(6):
            for j in range(6):
                atoms.append((rid, "DLPC", "PO4", 0.5 + i, 0.5 + j, zc)); rid += 1
    nw = 0
    for x in np.arange(0.3, 6, 0.7):
        for y in np.arange(0.3, 6, 0.7):
            for z in np.arange(2.0, 7.0, 0.7):
                atoms.append((rid, "W", "W", float(x), float(y), float(z)))
                rid += 1; nw += 1
    N = len(atoms)
    L = ["syn", str(N)]
    for s, (r, rn, an, x, y, z) in enumerate(atoms, 1):
        L.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                 % (r % 100000, rn, an, s % 100000, x, y, z))
    L.append("%10.5f%10.5f%10.5f" % (6, 6, 9))
    gro = tmp_path / "system.gro"; gro.write_text("\n".join(L) + "\n")
    top = tmp_path / "system.top"
    top.write_text('#include "ff.itp"\n[ system ]\nt\n[ molecules ]\n'
                   'PRO 1\nDLPC 72\nW %d\n' % nw)
    # CG partner: 12 BB beads
    rng = np.random.default_rng(3)
    pa = rng.uniform(0, 20, (12, 3))
    pcg = tmp_path / "partner_cg.pdb"
    with open(pcg, "w") as fh:
        for i, p in enumerate(pa, 1):
            fh.write("ATOM  %5d  BB  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
                     % (i, i, p[0], p[1], p[2]))
        fh.write("END\n")
    pitp = tmp_path / "PARTNER0.itp"
    pitp.write_text("[ moleculetype ]\nPARTNER0 1\n[ atoms ]\n")
    return gro, top, pcg, pitp, nw


def _gro_counts(gro):
    from collections import Counter
    L = gro.read_text().splitlines()
    n = int(L[1]); body = L[2:2 + n]
    box = [float(v) for v in L[2 + n].split()]
    return n, Counter(b[5:10].strip() for b in body), box


def test_place_partner_side_gap_and_consistency(tmp_path):
    """Partner lands on the chosen leaflet at the requested gap, box grows
    asymmetrically, and gro atom count matches the top molecule sum."""
    gro, top, pcg, pitp, nw = _mini_membrane_system(tmp_path)
    run("place_partner.py", "--system-gro", gro, "--system-top", top,
        "--partner-cg", pcg, "--partner-name", "PARTNER0",
        "--partner-itp", "PARTNER0.itp", "--side", "upper",
        "--gap", 1.5, "--water-nm", 3.0, "--lipids", "DLPC")
    n, c, box = _gro_counts(gro)
    assert c["PARTN"] == 12                    # 5-char-truncated resname present
    assert box[2] > 9.0                        # box grew
    # partner sits above the upper head (~5.2 nm) by ~gap
    L = gro.read_text().splitlines(); body = L[2:2 + n]
    pz = [float(b[36:44]) for b in body if b[5:10].strip() == "PARTN"]
    headz = [float(b[36:44]) for b in body
             if b[5:10].strip() == "DLPC" and float(b[36:44]) > 4.5]
    assert min(pz) > max(headz)                # partner is above the upper leaflet
    # gro vs top consistency
    apm = {"PRO": 10, "DLPC": 1, "W": 1, "NA": 1, "CL": 1, "PARTNER0": 12}
    tl = top.read_text().splitlines(); seen = False; tot = 0
    for x in tl:
        if x.strip().lower().startswith("[ molecules"):
            seen = True; continue
        if seen and x.strip() and not x.startswith("["):
            p = x.split(); tot += apm[p[0]] * int(p[1])
    assert tot == n                            # exact match (grompp-safe)


def test_place_partner_lower_side(tmp_path):
    """Lower-leaflet placement puts the partner below the lower head."""
    gro, top, pcg, pitp, nw = _mini_membrane_system(tmp_path)
    run("place_partner.py", "--system-gro", gro, "--system-top", top,
        "--partner-cg", pcg, "--partner-name", "PARTNER0",
        "--partner-itp", "PARTNER0.itp", "--side", "lower",
        "--gap", 1.5, "--water-nm", 3.0, "--lipids", "DLPC")
    n, c, box = _gro_counts(gro)
    L = gro.read_text().splitlines(); body = L[2:2 + n]
    pz = [float(b[36:44]) for b in body if b[5:10].strip() == "PARTN"]
    lowerz = [float(b[36:44]) for b in body
              if b[5:10].strip() == "DLPC" and float(b[36:44]) < 4.5]
    # after recenter everything shifted; check partner is on the OTHER side of
    # the membrane midplane from the upper leaflet
    midz = sum(float(b[36:44]) for b in body
               if b[5:10].strip() == "DLPC") / 72.0
    assert max(pz) < midz                      # partner below the bilayer middle


def test_add_partner_pull_flatbottom(tmp_path):
    """add_partner_pull adds a PARTNER index group and a one-sided flat-bottom
    pull (free to associate, blocked from wrapping) to the production mdp."""
    gro, top, pcg, pitp, nw = _mini_membrane_system(tmp_path)
    run("place_partner.py", "--system-gro", gro, "--system-top", top,
        "--partner-cg", pcg, "--partner-name", "PARTNER0",
        "--partner-itp", "PARTNER0.itp", "--side", "upper",
        "--gap", 1.5, "--water-nm", 3.0, "--lipids", "DLPC")
    ndx = tmp_path / "index.ndx"; ndx.write_text("[ SOLV ]\n1\n\n")
    mdp = tmp_path / "step7.mdp"; mdp.write_text("integrator = md\n")
    run("add_partner_pull.py", "--gro", gro, "--ndx", ndx, "--mdp", mdp,
        "--partner-name", "PARTNER0", "--lipids", "DLPC",
        "--margin", 1.0, "--k", 1000)
    m = mdp.read_text()
    assert "[ PARTNER ]" in ndx.read_text()
    assert "pull-coord1-type     = flat-bottom" in m
    assert "pull-coord1-dim      = N N Y" in m   # z-only confinement
    assert "pull-group2-name     = PARTNER" in m


# ------------------------------------------------ GUI helpers (app.py)
def test_app_param_builders():
    """The GUI form-to-parameter helpers produce the expected spec strings."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "app", os.path.join(HELPERS, "app.py"))
    app = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(app)  # streamlit is imported only inside main()

    assert app.build_lipid_spec([
        {"name": "CHOL", "ratio": 1, "head": "auto"},
        {"name": "DPG3", "ratio": 1, "head": "GM1"},
    ]) == "CHOL:1 DPG3:1:GM1"

    # asymmetric, no partner
    form = dict(pep_aa="p.pdb", m3_dir="/m3", gmx="gmx", asym=True,
                lipids="", upper="CHOL:1", lower="POPC:1", n_copy=4,
                box_x=None, box_y=None, water_nm=2.5, memb_thick_nm=4.0,
                temp=310, salt_m=0.15, tm_range="", partner_pdb="", n_rep=8)
    env = app.build_env(form, "/repo")
    assert env["UPPER"] == "CHOL:1" and "LIPIDS" not in env
    assert env["N_REP"] == "8"
    assert "PARTNER_PDB" not in env
    assert app.build_command(form, "/repo")[1].endswith("build_ensemble.sh")

    # prebuilt-multi suppresses N_COPY; partner maps to post-COBY env vars
    form2 = dict(pep_aa="p.pdb", m3_dir="/m3", gmx="gmx", asym=False,
                 lipids="CHOL:1 DIPC:1", upper="", lower="", n_copy=4,
                 prebuilt_multi=True, res_keep="A:1-10", tm_core="A:3-8",
                 box_x=None, box_y=None, water_nm=2.5, memb_thick_nm=4.0,
                 temp=310, salt_m=0.15, tm_range="", n_rep=1,
                 partner_pdb="/w/partner_aa.pdb", partner_side="lower",
                 partner_gap=2.0, partner_water=3.5)
    env2 = app.build_env(form2, "/repo")
    assert "N_COPY" not in env2                 # count comes from the PDB
    assert env2["PREBUILT_MULTI"] == "1"
    assert env2["PARTNER_PDB"] == "/w/partner_aa.pdb"
    assert env2["PARTNER_SIDE"] == "lower"
    assert env2["PARTNER_GAP"] == "2.0"
    assert app.build_command(form2, "/repo")[1].endswith("memble.sh")


# ------------------------------------------------ prebuild_orient
def test_prebuild_orient_trim_center_rigid(tmp_path):
    """Trim to per-chain ranges, center TM core at z=0, keep resSeq and rigid arrangement."""
    import numpy as np
    pdb = tmp_path / "multi.pdb"
    lines = []
    k = 0
    for ch, x0 in (("A", 0.0), ("B", 30.0)):
        for r in range(1, 16):              # resSeq 1..15 (11..15 will be trimmed)
            t = r * 3.0
            k += 1
            lines.append("ATOM  %5d  CA  GLY %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                         % (k, ch, r, x0 + t * 0.7, 0.0, 50 + t * 0.7))
    pdb.write_text("\n".join(lines) + "\nEND\n")
    out = tmp_path / "pre.pdb"
    run("prebuild_orient.py", "--in", pdb, "--out", out,
        "--keep", "A:1-10;B:1-10", "--core", "A:3-8;B:3-8")

    def rows(f):
        return [l for l in open(f) if l.startswith("ATOM")]

    def rc(l):
        return l[21], int(l[22:26]), (float(l[30:38]), float(l[38:46]), float(l[46:54]))

    r = rows(out)
    res = [rc(l)[1] for l in r]
    assert max(res) == 10 and min(res) == 1            # trimmed, resSeq preserved
    core_z = [rc(l)[2][2] for l in r if 3 <= rc(l)[1] <= 8]
    assert abs(np.mean(core_z)) < 1e-2                  # TM core centered at z=0

    def cen(f, ch):
        return np.mean([rc(l)[2] for l in rows(f) if rc(l)[0] == ch], axis=0)

    d_out = np.linalg.norm(cen(out, "A") - cen(out, "B"))
    # input centroids restricted to kept residues
    def cen_in(ch):
        return np.mean([rc(l)[2] for l in rows(pdb)
                        if rc(l)[0] == ch and rc(l)[1] <= 10], axis=0)
    d_in = np.linalg.norm(cen_in("A") - cen_in("B"))
    assert abs(d_out - d_in) < 1e-2                     # rigid: inter-chain distance preserved


# ------------------------------------------------ itp_to_struct
def test_itp_to_struct_from_connectivity(tmp_path):
    """A single-molecule CG structure is generated from the itp graph with
    correct bead count, head-bead detection, and bonded distances matching the
    ffbonded lengths (MDS init + constraint relaxation)."""
    import math
    ffb = tmp_path / "ffb.itp"
    ffb.write_text(
        "#define b_NC3_PO4_def 1 0.404 7000\n"
        "#define b_PO4_GL 1 0.40 7000\n"
        "#define b_GL_GL 1 0.37 7000\n"
        "#define b_GL_C1 1 0.47 5000\n"
        "#define b_CC 1 0.48 3800\n"
    )
    itp = tmp_path / "dlpc.itp"
    itp.write_text(
        "[moleculetype]\n DLPC 1\n[atoms]\n"
        "1 Q1 1 DLPC NC3 1 1.0\n2 Q5 1 DLPC PO4 2 -1.0\n"
        "3 SN4a 1 DLPC GL1 3 0\n4 SN4a 1 DLPC GL2 4 0\n"
        "5 C1 1 DLPC C1A 5 0\n6 C1 1 DLPC C2A 6 0\n7 C1 1 DLPC C3A 7 0\n"
        "8 C1 1 DLPC C1B 8 0\n9 C1 1 DLPC C2B 9 0\n10 C1 1 DLPC C3B 10 0\n"
        "[bonds]\n1 2 b_NC3_PO4_def\n2 3 b_PO4_GL\n2 4 b_PO4_GL\n3 4 b_GL_GL\n"
        "3 5 b_GL_C1\n5 6 b_CC\n6 7 b_CC\n4 8 b_GL_C1\n8 9 b_CC\n9 10 b_CC\n"
    )
    out = tmp_path / "DLPC.gro"
    res = run("itp_to_struct.py", "--itp", itp, "--mol", "DLPC",
              "--out", out, "--ffbonded", ffb)
    lines = out.read_text().splitlines()
    n = int(lines[1])
    body = lines[2:2 + n]
    assert n == 10
    pos = {l[10:15].strip(): (float(l[20:28]), float(l[28:36]), float(l[36:44]))
           for l in body}

    def dist(a, b):
        return math.dist(pos[a], pos[b])

    # bonded distances must match the ffbonded targets within 0.02 nm
    assert abs(dist("NC3", "PO4") - 0.404) < 0.02
    assert abs(dist("C1A", "C2A") - 0.48) < 0.02
    assert abs(dist("GL1", "C1A") - 0.47) < 0.02
    assert "UPDOWN DLPC up:NC3 upidx:1" in res.stdout


# ------------------------------------------------ declash
def test_declash_preserves_intramolecular(tmp_path):
    """Rigid-molecule declash separates overlapping molecules while keeping each
    molecule's internal distances exactly unchanged."""
    import math
    gro = tmp_path / "clash.gro"
    gro.write_text(
        "t\n    4\n"
        "    1AAA     B1    1   0.000   0.000   0.000\n"
        "    1AAA     B2    2   0.300   0.000   0.000\n"
        "    2BBB     B1    3   0.310   0.010   0.000\n"
        "    2BBB     B2    4   0.610   0.010   0.000\n"
        "   2.00000   2.00000   2.00000\n"
    )
    run("declash_gro.py", "--gro", gro, "--lipids", "AAA BBB",
        "--target", "0.30", "--iters", "300")
    body = gro.read_text().splitlines()[2:6]
    p = [(float(l[20:28]), float(l[28:36]), float(l[36:44])) for l in body]
    assert abs(math.dist(p[0], p[1]) - 0.30) < 0.01   # mol1 internal preserved
    assert abs(math.dist(p[2], p[3]) - 0.30) < 0.01   # mol2 internal preserved
    inter = min(math.dist(p[a], p[b]) for a in (0, 1) for b in (2, 3))
    assert inter > 0.27                               # molecules separated


# ------------------------------------------------ prebuild_orient bundle normal
def test_prebuild_orient_stands_helices_up(tmp_path):
    """A side-by-side multi-helix bundle (each helix axis = x, stacked along y)
    must be oriented so each helix points along z (membrane normal), not left
    lying flat. Guards against using the pooled principal axis."""
    import numpy as np
    from collections import defaultdict
    pdb = tmp_path / "bundle.pdb"
    lines = []
    s = 1
    for ci, ch in enumerate("ABCD"):
        y0 = ci * 20.0
        for k, res in enumerate(range(60, 83)):
            lines.append("ATOM  %5d  BB  ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                         % (s, ch, res, k * 1.5, y0, 0.0))
            s += 1
    pdb.write_text("\n".join(lines) + "\nEND\n")
    out = tmp_path / "oriented.pdb"
    run("prebuild_orient.py", "--in", pdb, "--out", out,
        "--keep", "A:60-82;B:60-82;C:60-82;D:60-82",
        "--core", "A:60-82;B:60-82;C:60-82;D:60-82")
    ch = defaultdict(list)
    for ln in out.read_text().splitlines():
        if ln.startswith("ATOM"):
            ch[ln[21]].append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    for c, pts in ch.items():
        A = np.array(pts) - np.mean(pts, axis=0)
        _, _, Vt = np.linalg.svd(A)
        ang = np.degrees(np.arccos(abs(Vt[0][2])))
        assert ang < 20, "chain %s still lying flat (z-angle %.0f)" % (c, ang)


# ------------------------------------------------ add_water box expansion
def test_add_water_expands_and_matches_top(tmp_path):
    """add_water grows box_z and solvates the new slabs; gro atom count must stay
    consistent with the per-moleculetype counts in the topology, and molecule
    blocks must remain in a single contiguous order."""
    import numpy as np
    from collections import Counter
    gro = tmp_path / "sys.gro"
    top = tmp_path / "sys.top"
    atoms = []
    rid = 1
    for k in range(40):
        atoms.append((rid, "PRO", "BB", 3.0, 3.0, 0.5 + 8.0 * k / 39)); 
    rid += 1
    nw = 0
    for x in np.arange(0.3, 6, 0.7):
        for y in np.arange(0.3, 6, 0.7):
            for z in np.arange(2.0, 7.0, 0.7):
                atoms.append((rid, "W", "W", float(x), float(y), float(z))); rid += 1; nw += 1
    L = ["syn", str(len(atoms))]
    for s, (r, rn, an, x, y, z) in enumerate(atoms, 1):
        L.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (r % 100000, rn, an, s % 100000, x, y, z))
    L.append("%10.5f%10.5f%10.5f" % (6, 6, 9))
    gro.write_text("\n".join(L) + "\n")
    top.write_text("[ molecules ]\nPRO 1\nW %d\n" % nw)
    run("add_water.py", "--gro", gro, "--top", top, "--water-nm", "3.5", "--salt", "0.15")
    body = gro.read_text().splitlines()
    n = int(body[1]); box = body[2 + n].split()
    assert float(box[2]) > 9.0                       # box grew
    apm = {"PRO": 40, "W": 1, "NA": 1, "CL": 1}
    tot = 0
    for line in top.read_text().splitlines():
        p = line.split()
        if len(p) == 2 and p[0] in apm:
            tot += apm[p[0]] * int(p[1])
    assert tot == n                                   # gro atoms == top atom-sum


# ------------------------------------------------ exact vsite construction
def test_itp_to_struct_exact_vsites_no_collapse(tmp_path):
    """Cholesterol ring virtual sites (virtual_sites3 funct 4) must be built from
    their exact definition so no bead lands on top of another (which caused
    infinite forces at minimisation). All intra-molecule bead pairs should be a
    sane distance apart."""
    import numpy as np, os, glob
    itp = None
    for c in ["/tmp/martini_v3.0.0_sterols_v1.itp"]:
        if os.path.exists(c):
            itp = c
    if itp is None:
        import pytest; pytest.skip("sterol itp not available in this environment")
    out = tmp_path / "CHOL.gro"
    args = ["--itp", itp, "--mol", "CHOL", "--out", str(out)]
    if os.path.exists("/tmp/ffb.itp"):
        args += ["--ffbonded", "/tmp/ffb.itp"]
    run("itp_to_struct.py", *args)
    L = out.read_text().splitlines(); n = int(L[1])
    P = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])]
                  for l in L[2:2 + n]])
    m = min(float(np.linalg.norm(P[i] - P[j]))
            for i in range(n) for j in range(i + 1, n))
    assert m > 0.15, "vsite collapse: min intra-pair %.3f nm" % m


def test_itp_to_struct_funct4_vsite_survives_reflection(tmp_path):
    """A sterol-like virtual_sites3 funct 4 site with a large out-of-plane
    coefficient must be reconstructed to match GROMACS even when orient_head_up
    reflects the molecule (a reflection flips the cross-product term). Regression
    for the infinite-force bug: vsites are built AFTER orientation."""
    import numpy as np
    itp = tmp_path / "sterol.itp"
    itp.write_text(
        "[ moleculetype ]\nCHOL 1\n\n[ atoms ]\n"
        " 1 TC4 1 CHOL ROH 1 0\n 2 SP1 1 CHOL R1 2 0\n 3 TC4 1 CHOL R2 3 0\n"
        " 4 TC4 1 CHOL R3 4 0\n 5 SC1 1 CHOL R4 5 0\n 6 SC3 1 CHOL R5 6 0\n"
        " 7 TC4 1 CHOL R6 7 0\n 8 SC1 1 CHOL C1 8 0\n 9 C1 1 CHOL C2 9 0\n\n"
        "[ bonds ]\n 8 9 1 0.425 1250.0\n\n"
        "[ constraints ]\n 8 3 1 0.75012\n 8 2 1 0.78504\n 3 2 1 0.43500\n"
        " 5 6 1 0.50000\n 6 7 1 0.46000\n 5 7 1 0.46000\n 2 5 1 0.45000\n"
        " 3 7 1 0.45000\n\n"
        "[ virtual_sites3 ]\n 1 8 3 2 4 1.08999 0.35891 0.22947\n"
        " 4 8 3 2 4 -0.41067 0.83597 -0.05311\n\n"
        "[ exclusions ]\n 1 2 3 4 5 6 7 8 9\n 2 3 4 5 6 7 8 9\n 3 4 5 6 7 8 9\n"
        " 4 5 6 7 8 9\n 5 6 7 8 9\n 6 7 8 9\n 7 8 9\n 8 9\n")
    out = tmp_path / "CHOL.gro"
    run("itp_to_struct.py", "--itp", itp, "--mol", "CHOL", "--out", out)
    L = out.read_text().splitlines(); n = int(L[1]); body = L[2:2 + n]
    nm = [b[10:15].strip() for b in body]
    P = {nm[i]: np.array([float(body[i][20:28]), float(body[i][28:36]),
                          float(body[i][36:44])]) for i in range(n)}
    order = ["ROH", "R1", "R2", "R3", "R4", "R5", "R6", "C1", "C2"]

    def vs(i, j, k, a, b, c):
        ri, rj, rk = P[order[i-1]], P[order[j-1]], P[order[k-1]]
        return ri + a*(rj-ri) + b*(rk-ri) + c*np.cross(rj-ri, rk-ri)

    roh = vs(8, 3, 2, 1.08999, 0.35891, 0.22947)
    r3 = vs(8, 3, 2, -0.41067, 0.83597, -0.05311)
    # must match the GROMACS funct-4 reconstruction to well under a bead radius
    assert np.linalg.norm(P["ROH"] - roh) < 0.01
    assert np.linalg.norm(P["R3"] - r3) < 0.01


def test_fix_vsites_rebuilds_flat_sterol_in_system(tmp_path):
    """fix_vsites overwrites COBY-style flat sterol vsites with the exact funct-4
    reconstruction, leaves a multi-residue protein and water untouched, and
    preserves the atom count (no index desync)."""
    import numpy as np
    st = tmp_path / "st.itp"
    st.write_text(
        "[ moleculetype ]\nCHOL 1\n\n[ atoms ]\n"
        " 1 TC4 1 CHOL ROH 1 0\n 2 SP1 1 CHOL R1 2 0\n 3 TC4 1 CHOL R2 3 0\n"
        " 4 TC4 1 CHOL R3 4 0\n 5 SC1 1 CHOL R4 5 0\n 6 SC3 1 CHOL R5 6 0\n"
        " 7 TC4 1 CHOL R6 7 0\n 8 SC1 1 CHOL C1 8 0\n 9 C1 1 CHOL C2 9 0\n\n"
        "[ virtual_sites3 ]\n 1 8 3 2 4 1.08999 0.35891 0.22947\n"
        " 4 8 3 2 4 -0.41067 0.83597 -0.05311\n")
    pro = tmp_path / "pro.itp"
    pro.write_text("[ moleculetype ]\nPRO 1\n[ atoms ]\n" +
                   "".join(" %d P5 1 PRO B%d %d 0\n" % (i, i, i)
                           for i in range(1, 6)))
    w = tmp_path / "w.itp"
    w.write_text("[ moleculetype ]\nW 1\n[ atoms ]\n 1 W 1 W W 1 0\n")
    top = tmp_path / "sys.top"
    top.write_text('#include "pro.itp"\n#include "st.itp"\n#include "w.itp"\n'
                   '[ system ]\nx\n[ molecules ]\nPRO 1\nCHOL 2\nW 3\n')
    order = ["ROH", "R1", "R2", "R3", "R4", "R5", "R6", "C1", "C2"]
    lines = []
    serial = 1

    def add(rid, rn, an, p):
        nonlocal serial
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                     % (rid % 100000, rn, an, serial % 100000, p[0], p[1], p[2]))
        serial += 1
    for i in range(5):
        add(1, "PRO", "B%d" % (i + 1), (3 + 0.1 * i, 3, 3))
    for c in range(2):
        cen = np.array([5.0 + c * 0.6, 5.0, 4.0])
        real = {"R1": cen + [0.1, 0, 0.2], "R2": cen + [0.3, 0, 0.1],
                "R4": cen + [0.2, 0.1, -0.2], "R5": cen + [0.0, 0.2, 0.1],
                "R6": cen + [0.15, 0.2, -0.1], "C1": cen + [0.2, 0, -0.3],
                "C2": cen + [0.2, 0, -0.7]}
        for an in order:
            add(100 + c, "CHOL", an, cen if an in ("ROH", "R3") else real[an])
    for i in range(3):
        add(900 + i, "W", "W", (i, i, i))
    N = len(lines)
    gro = tmp_path / "sys.gro"
    gro.write_text("t\n%d\n" % N + "\n".join(lines) +
                   "\n%10.5f%10.5f%10.5f\n" % (8, 8, 8))

    def roh_diff(rid):
        L = gro.read_text().splitlines(); n = int(L[1]); body = L[2:2 + n]
        b = {ln[10:15].strip(): np.array([float(ln[20:28]), float(ln[28:36]),
                                          float(ln[36:44])])
             for ln in body if ln[0:5].strip() == str(rid)}
        ri, rj, rk = b["C1"], b["R2"], b["R1"]
        roh = ri + 1.08999*(rj-ri) + 0.35891*(rk-ri) + 0.22947*np.cross(rj-ri, rk-ri)
        return np.linalg.norm(b["ROH"] - roh)

    assert roh_diff(100) > 0.1                      # COBY-style flat: wrong
    run("fix_vsites.py", "--gro", gro, "--top", top)
    assert roh_diff(100) < 0.01                      # rebuilt exactly
    assert roh_diff(101) < 0.01
    assert int(gro.read_text().splitlines()[1]) == N  # atom count preserved


def test_declash_resolves_pbc_boundary_clash(tmp_path):
    """A bead just outside the box clashing with the opposite face under PBC must
    be detected, separated, and wrapped back inside (regression for the
    infinite-force-on-a-water bug)."""
    import numpy as np
    box = 8.0
    rows = [(1, "W", "W", -0.088, 4.0, 4.0), (2, "W", "W", box - 0.05, 4.0, 4.0)]
    for k in range(3, 8):
        rows.append((k, "W", "W", 2.0 + 0.5 * k, 1.0, 1.0))
    gro = tmp_path / "pbc.gro"
    L = ["t", str(len(rows))]
    for s, (r, rn, an, x, y, z) in enumerate(rows, 1):
        L.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (r, rn, an, s, x, y, z))
    L.append("%10.5f%10.5f%10.5f" % (box, box, box))
    gro.write_text("\n".join(L) + "\n")
    run("declash_gro.py", "--gro", gro, "--lipids", "", "--target", 0.21,
        "--iters", 100)
    Lo = gro.read_text().splitlines(); n = int(Lo[1]); body = Lo[2:2 + n]
    P = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                  for b in body])
    rid = [int(b[0:5]) for b in body]
    bx = np.array([box] * 3)
    dm = 9.0
    for i in range(n):
        for j in range(i + 1, n):
            if rid[i] != rid[j]:
                d = P[i] - P[j]; d -= np.round(d / bx) * bx
                dm = min(dm, float(np.linalg.norm(d)))
    assert dm > 0.20                                   # boundary clash resolved
    assert ((P >= 0) & (P < box)).all()                # everything wrapped inside


def test_itp_to_struct_lipid_head_above_tails(tmp_path):
    """A two-tailed phospholipid template must come out as a proper lipid: head
    group at the top, BOTH acyl tails below it (regression for the inverted/
    splayed-tail bilayer where MDS put the two tails at opposite ends)."""
    import numpy as np
    itp = tmp_path / "dlpc.itp"
    itp.write_text(
        "[ moleculetype ]\nDLPC 1\n[ atoms ]\n"
        " 1 Q1 1 DLPC NC3 1 1.0\n 2 Q5 1 DLPC PO4 2 -1.0\n"
        " 3 SN4a 1 DLPC GL1 3 0\n 4 N4a 1 DLPC GL2 4 0\n"
        " 5 C1 1 DLPC C1A 5 0\n 6 C4h 1 DLPC D2A 6 0\n"
        " 7 C4h 1 DLPC D3A 7 0\n 8 C1 1 DLPC C4A 8 0\n"
        " 9 C1 1 DLPC C1B 9 0\n 10 C4h 1 DLPC D2B 10 0\n"
        " 11 C4h 1 DLPC D3B 11 0\n 12 C1 1 DLPC C4B 12 0\n"
        "[ bonds ]\n 1 2 1 0.47 1250\n 2 3 1 0.47 1250\n 3 4 1 0.37 1250\n"
        " 3 5 1 0.47 1250\n 5 6 1 0.47 1250\n 6 7 1 0.47 1250\n 7 8 1 0.47 1250\n"
        " 4 9 1 0.47 1250\n 9 10 1 0.47 1250\n 10 11 1 0.47 1250\n 11 12 1 0.47 1250\n")
    out = tmp_path / "DLPC.gro"
    run("itp_to_struct.py", "--itp", itp, "--mol", "DLPC", "--out", out)
    L = out.read_text().splitlines(); n = int(L[1]); body = L[2:2 + n]
    P = {b[10:15].strip(): np.array([float(b[20:28]), float(b[28:36]),
                                     float(b[36:44])]) for b in body}
    # head above both tail ends
    assert P["PO4"][2] > P["C4A"][2]
    assert P["PO4"][2] > P["C4B"][2]
    assert P["NC3"][2] > P["C4A"][2] and P["NC3"][2] > P["C4B"][2]
    # bond lengths preserved (~0.47)
    assert abs(np.linalg.norm(P["NC3"] - P["PO4"]) - 0.47) < 0.05
    assert abs(np.linalg.norm(P["C1A"] - P["D2A"]) - 0.47) < 0.05


def test_write_conect_pdb_bonds(tmp_path):
    """CONECT-PDB carries the right number of bonds (per-molecule bonds x count),
    so viewers show the system and trajectory WITH bonds by default."""
    pro = tmp_path / "pro.itp"
    pro.write_text("[ moleculetype ]\nPRO 1\n[ atoms ]\n" +
                   "".join(" %d P5 1 PRO B%d %d 0\n" % (i, i, i)
                           for i in range(1, 6)) +
                   "[ bonds ]\n 1 2 1\n 2 3 1\n 3 4 1\n 4 5 1\n")
    dlpc = tmp_path / "dlpc.itp"
    dlpc.write_text(
        "[ moleculetype ]\nDLPC 1\n[ atoms ]\n" +
        "".join(" %d C 1 DLPC X%d %d 0\n" % (i, i, i) for i in range(1, 13)) +
        "[ bonds ]\n" + "".join(" %d %d 1\n" % (i, i + 1) for i in range(1, 12)))
    w = tmp_path / "w.itp"
    w.write_text("[ moleculetype ]\nW 1\n[ atoms ]\n 1 W 1 W W 1 0\n")
    top = tmp_path / "v.top"
    top.write_text('#include "pro.itp"\n#include "dlpc.itp"\n#include "w.itp"\n'
                   '[ molecules ]\nPRO 1\nDLPC 2\nW 3\n')
    rows = []
    s = 1
    def add(rid, rn, an):
        nonlocal s
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                    % (rid % 100000, rn, an, s % 100000, 0.1*s, 0.1*s, 0.1*s))
        s += 1
    for i in range(5):
        add(1, "PRO", "B%d" % (i + 1))
    for c in range(2):
        for i in range(12):
            add(10 + c, "DLPC", "X%d" % (i + 1))
    for i in range(3):
        add(90 + i, "W", "W")
    N = len(rows)
    gro = tmp_path / "v.gro"
    gro.write_text("t\n%d\n" % N + "\n".join(rows) + "\n5 5 5\n")
    out = tmp_path / "v.pdb"
    run("write_conect_pdb.py", "--gro", gro, "--top", top, "--out", out)
    txt = out.read_text()
    # count unique bonds from CONECT (each bond appears twice, i->j and j->i)
    pairs = set()
    for ln in txt.splitlines():
        if ln.startswith("CONECT"):
            nums = [int(ln[6 + 5*k:6 + 5*(k+1)]) for k in range((len(ln)-6)//5)]
            a = nums[0]
            for b in nums[1:]:
                pairs.add(tuple(sorted((a, b))))
    assert len(pairs) == 4 + 11 * 2          # PRO 4 + DLPC 11 x2
    assert "CRYST1" in txt                    # box record present for viewers


def test_write_psf_atoms_and_bonds(tmp_path):
    """PSF carries NATOM matching the gro and NBOND = per-molecule bonds x count,
    so loading system.psf then a trajectory shows bonds on every frame."""
    import re
    pro = tmp_path / "pro.itp"
    pro.write_text("[ moleculetype ]\nPRO 1\n[ atoms ]\n" +
                   "".join(" %d P5 1 PRO B%d %d 0\n" % (i, i, i)
                           for i in range(1, 6)) +
                   "[ bonds ]\n 1 2 1\n 2 3 1\n 3 4 1\n 4 5 1\n")
    dlpc = tmp_path / "dlpc.itp"
    dlpc.write_text(
        "[ moleculetype ]\nDLPC 1\n[ atoms ]\n" +
        "".join(" %d C 1 DLPC X%d %d 0\n" % (i, i, i) for i in range(1, 13)) +
        "[ bonds ]\n" + "".join(" %d %d 1\n" % (i, i + 1) for i in range(1, 12)))
    w = tmp_path / "w.itp"
    w.write_text("[ moleculetype ]\nW 1\n[ atoms ]\n 1 W 1 W W 1 0\n")
    top = tmp_path / "v.top"
    top.write_text('#include "pro.itp"\n#include "dlpc.itp"\n#include "w.itp"\n'
                   '[ molecules ]\nPRO 1\nDLPC 2\nW 3\n')
    rows = []
    s = 1
    def add(rid, rn, an):
        nonlocal s
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                    % (rid % 100000, rn, an, s % 100000, 0.1*s, 0.1*s, 0.1*s))
        s += 1
    for i in range(5):
        add(1, "PRO", "B%d" % (i + 1))
    for c in range(2):
        for i in range(12):
            add(10 + c, "DLPC", "X%d" % (i + 1))
    for i in range(3):
        add(90 + i, "W", "W")
    N = len(rows)
    gro = tmp_path / "v.gro"
    gro.write_text("t\n%d\n" % N + "\n".join(rows) + "\n5 5 5\n")
    out = tmp_path / "system.psf"
    run("write_psf.py", "--gro", gro, "--top", top, "--out", out)
    psf = out.read_text()
    assert "PSF" in psf.splitlines()[0]
    na = int(re.search(r"(\d+) !NATOM", psf).group(1))
    nb = int(re.search(r"(\d+) !NBOND", psf).group(1))
    assert na == N
    assert nb == 4 + 11 * 2


def test_fix_protein_resid_restores_per_chain_numbering(tmp_path):
    """Protein residues renumbered from 1 by martinize are restored to the input
    PDB numbering (e.g. 54..), per chain, from oriented_aa.pdb. Not hardcoded."""
    ori = tmp_path / "oriented_aa.pdb"
    ori.write_text(
        "ATOM      1  N   GLU A  54       0.000   0.000   0.000\n"
        "ATOM      2  N   GLY A  55       0.000   0.000   0.000\n"
        "ATOM      3  N   CYS A  56       0.000   0.000   0.000\n"
        "ATOM      4  N   GLU B  54       0.000   0.000   0.000\n"
        "ATOM      5  N   GLY B  55       0.000   0.000   0.000\n"
        "ATOM      6  N   CYS B  56       0.000   0.000   0.000\n")
    for mol in ("molecule_0", "molecule_1"):
        (tmp_path / (mol + ".itp")).write_text(
            "[ moleculetype ]\n%s 1\n[ atoms ]\n" % mol +
            " 1 P5 1 GLU BB 1 0\n 2 P5 1 GLU SC1 1 0\n 3 P5 2 GLY BB 2 0\n"
            " 4 P5 3 CYS BB 3 0\n 5 P5 3 CYS SC1 3 0\n")
    top = tmp_path / "system.top"
    top.write_text('#include "molecule_0.itp"\n#include "molecule_1.itp"\n'
                   '[ molecules ]\nmolecule_0 1\nmolecule_1 1\n')
    rows = []
    s = 1
    def add(rid, rn, an):
        nonlocal s
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                    % (rid, rn, an, s, 0.1*s, 0.1*s, 0.1*s))
        s += 1
    for _ch in range(2):
        for rid, rn, an in [(1, "GLU", "BB"), (1, "GLU", "SC1"), (2, "GLY", "BB"),
                            (3, "CYS", "BB"), (3, "CYS", "SC1")]:
            add(rid, rn, an)
    gro = tmp_path / "system.gro"
    gro.write_text("t\n%d\n" % len(rows) + "\n".join(rows) + "\n5 5 5\n")
    run("fix_protein_resid.py", "--gro", gro, "--oriented", ori,
        "--top", top, "--itp-dir", str(tmp_path))
    L = gro.read_text().splitlines()
    body = L[2:2 + int(L[1])]
    resids = [int(b[0:5]) for b in body]
    assert resids[:5] == [54, 54, 55, 56, 56]      # chain A restored
    assert resids[5:10] == [54, 54, 55, 56, 56]    # chain B restored


def test_make_protein_whole_unsplits_chain(tmp_path):
    """A protein chain split across the z boundary is made contiguous and the
    system recentered, so consecutive backbone beads are bonded-distance apart
    (regression for the LINCS/infinite-force blowup at minimisation)."""
    (tmp_path / "molecule_0.itp").write_text(
        "[ moleculetype ]\nmolecule_0 1\n[ atoms ]\n" +
        "".join(" %d P5 %d ALA BB %d 0\n" % (i, i, i) for i in range(1, 7)) +
        "[ bonds ]\n" + "".join(" %d %d 1\n" % (i, i + 1) for i in range(1, 6)))
    (tmp_path / "w.itp").write_text("[ moleculetype ]\nW 1\n[ atoms ]\n 1 W 1 W W 1 0\n")
    top = tmp_path / "w.top"
    top.write_text('#include "molecule_0.itp"\n#include "w.itp"\n'
                   '[ molecules ]\nmolecule_0 1\nW 2\n')
    zs = [6.0, 6.35, 6.7, 7.05 + 14, 7.4, 7.75]   # bead 4 wrapped a full box
    rows = []
    s = 1
    for z in zs:
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (1, "ALA", "BB", s, 1.0, 1.0, z))
        s += 1
    for z in (1.0, 15.0):
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (2, "W", "W", s, 1.0, 1.0, z))
        s += 1
    gro = tmp_path / "w.gro"
    gro.write_text("t\n%d\n" % len(rows) + "\n".join(rows) +
                   "\n  14.00000  14.00000  14.00000\n")
    run("make_protein_whole.py", "--gro", gro, "--top", top, "--itp-dir", str(tmp_path))
    L = gro.read_text().splitlines()
    body = L[2:2 + int(L[1])]
    zc = [float(b[36:44]) for b in body[:6]]
    # consecutive backbone beads within a bond length (no box jump)
    for a, b in zip(zc, zc[1:]):
        assert abs(b - a) < 1.0
    # water wrapped inside the box
    for b in body[6:8]:
        assert 0.0 <= float(b[36:44]) < 14.0


def test_declash_freeze_protein_keeps_protein_fixed():
    """With frozen={-1} the protein never moves; only the clashing solvent does,
    so the protein cannot be distorted or re-split during the final declash."""
    import numpy as np, importlib
    d = importlib.import_module("declash_gro")
    coords = np.array([[1.0, 1.0, 1.0], [1.05, 1.0, 1.0], [1.02, 1.0, 1.0]])
    mol = np.array([-1, -1, 5])           # two protein beads + one water
    out, rem = d.declash(coords.copy(), mol, target=0.3, iters=300,
                         box=[10, 10, 10], frozen={-1})
    assert np.allclose(out[0], [1.0, 1.0, 1.0])     # protein frozen
    assert np.allclose(out[1], [1.05, 1.0, 1.0])
    assert np.linalg.norm(out[2] - out[0]) > 0.29   # water pushed clear


def test_itp_to_struct_no_intramolecular_overlap(tmp_path):
    """Template beads must keep a minimum non-bonded spacing; the two acyl-tail
    roots must not collapse onto each other (regression for the DLPC C1A/C1B
    overlap that gave an infinite force at minimisation)."""
    import numpy as np, itertools
    itp = tmp_path / "dlpc.itp"
    itp.write_text(
        "[ moleculetype ]\nDLPC 1\n[ atoms ]\n"
        " 1 Q1 1 DLPC NC3 1 1.0\n 2 Q5 1 DLPC PO4 2 -1.0\n"
        " 3 SN4a 1 DLPC GL1 3 0\n 4 N4a 1 DLPC GL2 4 0\n"
        " 5 C1 1 DLPC C1A 5 0\n 6 C4h 1 DLPC D2A 6 0\n"
        " 7 C4h 1 DLPC D3A 7 0\n 8 C1 1 DLPC C4A 8 0\n"
        " 9 C1 1 DLPC C1B 9 0\n 10 C4h 1 DLPC D2B 10 0\n"
        " 11 C4h 1 DLPC D3B 11 0\n 12 C1 1 DLPC C4B 12 0\n"
        "[ bonds ]\n 1 2 1 0.47 1250\n 2 3 1 0.47 1250\n 3 4 1 0.37 1250\n"
        " 3 5 1 0.47 1250\n 5 6 1 0.47 1250\n 6 7 1 0.47 1250\n 7 8 1 0.47 1250\n"
        " 4 9 1 0.47 1250\n 9 10 1 0.47 1250\n 10 11 1 0.47 1250\n 11 12 1 0.47 1250\n")
    out = tmp_path / "DLPC.gro"
    run("itp_to_struct.py", "--itp", itp, "--mol", "DLPC", "--out", out)
    L = out.read_text().splitlines(); n = int(L[1]); body = L[2:2 + n]
    P = [np.array([float(b[20:28]), float(b[28:36]), float(b[36:44])])
         for b in body]
    mind = min(np.linalg.norm(P[a] - P[b])
               for a, b in itertools.combinations(range(n), 2))
    assert mind > 0.20            # no two beads collapsed together


def test_make_protein_whole_centers_membrane(tmp_path):
    """The bilayer is centered at the box middle in z (membrane circular-mean),
    not on the protein, so the hydrophobic core does not land on the box edge
    and the bilayer is not displayed split."""
    import numpy as np
    (tmp_path / "p.itp").write_text(
        "[ moleculetype ]\nmolecule_0 1\n[ atoms ]\n" +
        "".join(" %d P5 %d ALA BB %d 0\n" % (i, i, i) for i in range(1, 6)) +
        "[ bonds ]\n" + "".join(" %d %d 1\n" % (i, i + 1) for i in range(1, 5)))
    (tmp_path / "d.itp").write_text(
        "[ moleculetype ]\nDLPC 1\n[ atoms ]\n" +
        "".join(" %d C 1 DLPC %s %d 0\n" % (i, nm, i)
                for i, nm in enumerate(["PO4", "GL1", "C1A", "C4A"], 1)) +
        "[ bonds ]\n 1 2 1\n 2 3 1\n 3 4 1\n")
    (tmp_path / "w.itp").write_text("[ moleculetype ]\nW 1\n[ atoms ]\n 1 W 1 W W 1 0\n")
    top = tmp_path / "c.top"
    top.write_text('#include "p.itp"\n#include "d.itp"\n#include "w.itp"\n'
                   '[ molecules ]\nmolecule_0 1\nDLPC 4\nW 3\n')
    rows = []
    s = 1
    def add(rid, rn, an, z):
        nonlocal s
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (rid, rn, an, s, 1.0, 1.0, z))
        s += 1
    for z in (13.0, 13.5, 0.2, 0.7, 1.2):
        add(1, "ALA", "BB", z)
    for rid, base in ((10, 2.0), (11, 2.2)):
        add(rid, "DLPC", "PO4", base); add(rid, "DLPC", "GL1", base - 0.7)
        add(rid, "DLPC", "C1A", base - 1.3); add(rid, "DLPC", "C4A", base - 1.8)
    for rid, base in ((12, 11.5), (13, 11.3)):
        add(rid, "DLPC", "PO4", base); add(rid, "DLPC", "GL1", base + 0.8)
        add(rid, "DLPC", "C1A", base + 1.5); add(rid, "DLPC", "C4A", base + 2.1)
    for z in (6.0, 7.0, 8.0):
        add(20, "W", "W", z)
    gro = tmp_path / "c.gro"
    gro.write_text("t\n%d\n" % len(rows) + "\n".join(rows) +
                   "\n  14.00000  14.00000  14.00000\n")
    run("make_protein_whole.py", "--gro", gro, "--top", top, "--itp-dir", str(tmp_path))
    L = gro.read_text().splitlines()
    body = L[2:2 + int(L[1])]
    c4a = [float(b[36:44]) for b in body if b[10:15].strip() == "C4A"]
    # hydrophobic core (tail ends) near the box center, not at the edge
    assert abs(np.mean(c4a) - 7.0) < 1.5


def test_add_water_sizes_box_to_protein_not_membrane(tmp_path):
    """box_z must be sized from the true protein z-extent plus the water layer,
    not from the lipid-filled thin box, or a protruding TM-JM protein ends up
    with almost no water beyond its ends (regression for the protein leaving the
    membrane when pressure coupling shrinks an over-tight box)."""
    import numpy as np
    # thin box_z=8; protein (amino acids) spans z 1..7 (=6 nm); lipids fill 2..6
    rows = []
    s = 1
    def add(rn, an, z):
        nonlocal s
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (1, rn, an, s, 1.0, 1.0, z))
        s += 1
    for z in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0):
        add("ALA", "BB", z)
    for z in (2.5, 3.5, 4.5, 5.5):
        add("DLPC", "C1A", z)
    rng = __import__("random").Random(0)
    for _ in range(60):
        add("W", "W", rng.uniform(0.3, 7.7))
    gro = tmp_path / "s.gro"
    gro.write_text("t\n%d\n" % len(rows) + "\n".join(rows) +
                   "\n   6.00000   6.00000   8.00000\n")
    top = tmp_path / "s.top"
    top.write_text("[ molecules ]\n")     # not needed by add_water
    run("add_water.py", "--gro", gro, "--top", top, "--water-nm", "2.5", "--salt", "0")
    L = gro.read_text().splitlines()
    box_z = float(L[-1].split()[2])
    # protein span 6 nm + 2*2.5 = 11 nm, not 8 (thin) and not 8+5 measured on lipids
    assert box_z > 10.0


def test_balance_apl_matches_leaflet_areas():
    """Per-leaflet apl from composition keeps the average at the packing apl but
    sets the leaflet ratio so the two leaflets come out area-matched."""
    import subprocess, sys, os
    here = os.path.join(os.path.dirname(__file__), "..")
    out = subprocess.run(
        [sys.executable, os.path.join(here, "balance_apl.py"),
         "--upper", "CHOL:1 DLPC:1 PSM:1", "--lower", "CHOL:1 DLPC:1 DOPS:1",
         "--base-apl", "0.70"], capture_output=True, text=True)
    up, lo = (float(x) for x in out.stdout.split())
    # average preserved
    assert abs((up + lo) / 2 - 0.70) < 1e-3
    # lower leaflet (larger DOPS) gets the larger apl -> fewer lipids -> balance
    assert lo > up
    # leaflet mean-area ratio reproduced (PSM 0.50 vs DOPS 0.60 in the table)
    a_up = (0.40 + 0.60 + 0.50) / 3
    a_lo = (0.40 + 0.60 + 0.60) / 3
    assert abs(up / lo - a_up / a_lo) < 1e-3


def test_inject_posres_resid_range_tm_core_only(tmp_path):
    """Restraints can be limited to a residue range (the TM core), leaving the
    juxtamembrane backbone free so it is not frozen in an extended pose during
    equilibration."""
    itp = tmp_path / "molecule_0.itp"
    itp.write_text("[ moleculetype ]\nmolecule_0 1\n[ atoms ]\n" +
                   "".join(" %d P5 %d ALA BB %d 0\n" % (i, i, i)
                           for i in range(1, 21)) + "[ bonds ]\n 1 2 1\n")
    run("inject_posres.py", "--itp", itp, "--mol", "molecule_0", "--beads", "BB",
        "--resid-min", "12", "--resid-max", "15", "--stage", "POSRES_STEP1:1000")
    txt = itp.read_text()
    blk = txt.split("#ifdef POSRES_STEP1")[1].split("#endif")[0]
    ids = [int(l.split()[0]) for l in blk.splitlines()
           if l.split() and l.split()[0].isdigit()]
    assert ids == [12, 13, 14, 15]


def test_make_protein_whole_unifies_chains_one_image(tmp_path):
    """All protein chains are brought into the same periodic image (preserving
    their input xy positions), so a multi-chain bundle is not split into separate
    z clusters that leave the membrane two at a time."""
    import numpy as np
    for m in range(4):
        (tmp_path / ("molecule_%d.itp" % m)).write_text(
            "[ moleculetype ]\nmolecule_%d 1\n[ atoms ]\n" % m +
            "".join(" %d P5 %d ALA BB %d 0\n" % (i, i, i) for i in range(1, 6)) +
            "[ bonds ]\n" + "".join(" %d %d 1\n" % (i, i + 1) for i in range(1, 5)))
    (tmp_path / "w.itp").write_text("[ moleculetype ]\nW 1\n[ atoms ]\n 1 W 1 W W 1 0\n")
    top = tmp_path / "q.top"
    top.write_text("".join('#include "molecule_%d.itp"\n' % m for m in range(4)) +
                   '#include "w.itp"\n[ molecules ]\n' +
                   "".join("molecule_%d 1\n" % m for m in range(4)) + "W 4\n")
    rows = []
    s = 1
    box = 15.0
    def add(rn, an, x, y, z):
        nonlocal s
        rows.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                    % (1, rn, an, s, x, y, z % box))
        s += 1
    # chains 2,3 shifted by +box in z (different image)
    for x, y, z0 in [(5, 5, 5), (7, 5, 5), (5, 8, 5 + 15), (7, 8, 5 + 15)]:
        for k in range(5):
            add("ALA", "BB", x, y, z0 + k)
    for i in range(4):
        add("W", "W", 2 + i, 2, 2.0)
    gro = tmp_path / "q.gro"
    gro.write_text("t\n%d\n" % len(rows) + "\n".join(rows) +
                   "\n  15.00000  15.00000  15.00000\n")
    run("make_protein_whole.py", "--gro", gro, "--top", top, "--itp-dir", str(tmp_path))
    L = gro.read_text().splitlines()
    body = L[2:2 + int(L[1])]
    X = np.array([[float(b[20:28]), float(b[28:36]), float(b[36:44])]
                  for b in body if b[5:10].strip() == "ALA"])
    # all chains in one z cluster (spread ~ protein thickness, not ~box)
    assert X[:, 2].max() - X[:, 2].min() < 6.0
    # xy positions preserved (chains still at their distinct spots)
    xs = [X[m*5:(m+1)*5, 0].mean() for m in range(4)]
    assert max(xs) - min(xs) > 1.5


def test_prebuild_orient_keeps_cores_coplanar(tmp_path):
    """A multi-chain TM bundle with tilted helices must come out with the TM-core
    centroids coplanar (membrane-plane normal aligned to z), not spread along z
    (regression for the assembly being tilted so chains sit at staircase heights
    and leave the membrane two at a time)."""
    import numpy as np
    from collections import defaultdict
    tilt = np.deg2rad(20)
    axis = np.array([np.sin(tilt), 0, np.cos(tilt)])
    lines = []
    serial = 1
    for ch, (gx, gy) in {"A": (0, 0), "B": (0, 40), "C": (100, 0),
                         "D": (100, 40)}.items():
        for r in range(54, 104):
            p = np.array([gx, gy, 45.0]) + axis * ((r - 76) * 1.5)
            lines.append("ATOM  %5d  CA  ALA %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                         % (serial, ch, r, p[0], p[1], p[2]))
            serial += 1
    pdb = tmp_path / "synth4.pdb"
    pdb.write_text("\n".join(lines) + "\nEND\n")
    out = tmp_path / "synth4_ori.pdb"
    run("prebuild_orient.py", "--in", pdb, "--out", out,
        "--keep", "A:54-103;B:54-103;C:54-103;D:54-103",
        "--core", "A:65-88;B:65-88;C:65-88;D:65-88")
    tm = defaultdict(list)
    for ln in out.read_text().splitlines():
        if ln.startswith("ATOM"):
            r = int(ln[22:26])
            if 65 <= r <= 88:
                tm[ln[21]].append(float(ln[46:54]))
    cz = [np.mean(v) for v in tm.values()]
    assert max(cz) - min(cz) < 0.5     # cores coplanar (nm)
