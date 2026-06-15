#!/usr/bin/env python3
"""
Streamlit GUI for the Martini 3 membrane-protein system builder.

This front end is Martini 3 only. It collects the same parameters the
command-line pipeline accepts, builds the command, runs it, and shows the
leaflet area check and a structure preview. The backend tools (COBY,
martinize2, GROMACS, DSSP) must be installed; the GUI builds and runs the
command, it does not replace those tools.

Run:  streamlit run app.py
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path

# ----------------------------------------------------------------------
# Pure helpers (no Streamlit) so they can be unit tested.
# ----------------------------------------------------------------------
def build_lipid_spec(rows):
    """rows: list of dicts {name, ratio, head}. -> 'NAME:ratio[:head] ...'."""
    toks = []
    for r in rows:
        name = (r.get("name") or "").strip()
        if not name:
            continue
        ratio = str(r.get("ratio") or 1).strip()
        head = (r.get("head") or "").strip()
        tok = "%s:%s" % (name, ratio)
        if head and head.lower() != "auto":
            tok += ":" + head
        toks.append(tok)
    return " ".join(toks)


def build_partner_spec(partners):
    """partners: list of dicts. -> COBY-free PARTNER string with '|' keys, ';' entries."""
    entries = []
    for p in partners:
        pdb = (p.get("pdb") or "").strip()
        if not pdb:
            continue
        kv = ["pdb=%s" % pdb, "offset=%s" % p.get("offset", 3)]
        if p.get("x") not in (None, ""):
            kv.append("x=%s" % p["x"])
        if p.get("y") not in (None, ""):
            kv.append("y=%s" % p["y"])
        if p.get("below"):
            kv.append("below=1")
        rot = p.get("rotate", "random")
        kv.append("rotate=%s" % rot)
        if rot == "euler":
            kv.append("euler=%s" % (p.get("euler") or "0,0,0"))
        entries.append("|".join(kv))
    return ";".join(entries)


def build_env(form, helpers_dir):
    """Assemble the environment dict for the build command."""
    env = dict(os.environ)
    env["PY"] = sys.executable
    env["MARTINIZE2"] = str(Path(sys.executable).parent / "martinize2")
    env["M3_DIR"] = form["m3_dir"]
    env["GMX"] = form["gmx"]
    # DSSP 4.x is rejected by vermouth; 'mdtraj' uses mdtraj's internal DSSP
    env["DSSP"] = form.get("dssp") or "mdtraj"
    env["SS_MODE"] = form.get("ss_mode") or "dssp"
    if form.get("multi_tm"):
        env["MULTI_TM"] = "1"
        if form.get("multi_tm_minlen"):
            env["MULTI_TM_MINLEN"] = str(form["multi_tm_minlen"])
    if form.get("nterm_side") and form["nterm_side"] != "auto":
        env["NTERM_SIDE"] = form["nterm_side"]
    if form.get("ss_override"):
        env["SS_OVERRIDE"] = form["ss_override"]
    if form.get("z_shift"):
        env["Z_SHIFT"] = str(form["z_shift"])
    if form.get("prebuilt_multi"):
        env["PREBUILT_MULTI"] = "1"
        env["RES_KEEP"] = form["res_keep"]
        if form.get("tm_core"):
            env["TM_CORE"] = form["tm_core"]
    env["HELPER_REP"] = str(Path(helpers_dir) / "replicate_and_fix_top.py")
    env["HELPER_POS"] = str(Path(helpers_dir) / "inject_posres.py")
    env["HELPER_ORI"] = str(Path(helpers_dir) / "orient_tm.py")
    env["HELPER_SSDSSP"] = str(Path(helpers_dir) / "ss_from_dssp.py")
    env["HELPER_AREA"] = str(Path(helpers_dir) / "leaflet_area_check.py")
    env["HELPER_PART"] = str(Path(helpers_dir) / "place_partner.py")
    env["HELPER_PRE"] = str(Path(helpers_dir) / "prebuild_orient.py")
    env["HELPER_ITP2STRUCT"] = str(Path(helpers_dir) / "itp_to_struct.py")
    env["HELPER_DECLASH"] = str(Path(helpers_dir) / "declash_gro.py")
    env["HELPER_ADDWATER"] = str(Path(helpers_dir) / "add_water.py")
    env["HELPER_PARTPULL"] = str(Path(helpers_dir) / "add_partner_pull.py")
    env["HELPER_FIXVS"] = str(Path(helpers_dir) / "fix_vsites.py")
    env["HELPER_WHOLE"] = str(Path(helpers_dir) / "make_protein_whole.py")
    env["HELPER_ZSHIFT"] = str(Path(helpers_dir) / "shift_protein_z.py")
    env["HELPER_FIXRESID"] = str(Path(helpers_dir) / "fix_protein_resid.py")
    env["HELPER_MINDIST"] = str(Path(helpers_dir) / "check_min_distance.py")
    env["HELPER_GENSS"] = str(Path(helpers_dir) / "gen_ss.py")
    if form.get("asym"):
        env["UPPER"] = form["upper"]
        env["LOWER"] = form["lower"]
    else:
        env["LIPIDS"] = form["lipids"]
    # protein count comes from the PDB itself in prebuilt-multi mode
    if not form.get("prebuilt_multi"):
        env["N_COPY"] = str(form["n_copy"])
    env["WATER_NM"] = str(form["water_nm"])
    env["MEMB_THICK_NM"] = str(form["memb_thick_nm"])
    env["TEMP"] = str(form["temp"])
    env["SALT_M"] = str(form["salt_m"])
    if form.get("coby_apl"):
        env["COBY_APL"] = str(form["coby_apl"])
    if form.get("box_x"):
        env["BOX_X"] = str(form["box_x"])
    if form.get("box_y"):
        env["BOX_Y"] = str(form["box_y"])
    if form.get("tm_range"):
        env["TM_RANGE"] = form["tm_range"]
    # post-COBY peripheral protein (atomistic PDB -> martinize -> place)
    if form.get("partner_pdb"):
        env["PARTNER_PDB"] = form["partner_pdb"]
        env["PARTNER_SIDE"] = form.get("partner_side", "upper")
        env["PARTNER_GAP"] = str(form.get("partner_gap", 1.5))
        env["PARTNER_WATER"] = str(form.get("partner_water", 3.0))
    if form.get("n_rep") and int(form["n_rep"]) > 1:
        env["N_REP"] = str(form["n_rep"])
    return env


def build_command(form, helpers_dir):
    """Choose the single build or the ensemble driver."""
    ensemble = form.get("n_rep") and int(form["n_rep"]) > 1
    script = "build_ensemble.sh" if ensemble else "memble.sh"
    return ["bash", str(Path(helpers_dir) / script), form["pep_aa"]]


# ----------------------------------------------------------------------
# Streamlit UI
# ----------------------------------------------------------------------
def main():
    import streamlit as st

    st.set_page_config(page_title="Martini 3 membrane builder", layout="wide")
    st.title("Martini 3 membrane-protein system builder")
    st.caption("Martini 3 only. Builds CHARMM-GUI-equivalent inputs (GROMACS "
               "topology, PSF, CRD, staged equilibration) for arbitrary "
               "compositions, leaflet asymmetry, and peripheral partners.")

    here = Path(__file__).resolve().parent

    with st.sidebar:
        st.header("Environment")
        m3_dir = st.text_input("Martini 3 lipidome directory (M3_DIR)", "")
        gmx = st.text_input("GROMACS 2023 binary (GMX)", "gmx")
        out_dir = st.text_input("Output directory (where the built system is "
                                "saved)", str(Path(__file__).resolve().parent))
        dssp = st.text_input("DSSP (path, or 'mdtraj' to use mdtraj's DSSP)",
                             "mdtraj")
        st.caption("martinize2 needs DSSP 2.2.1 or 3.x; 4.x is not compatible. "
                   "Give a 3.x binary path, 'mdtraj', or use SS_MODE=tm (no DSSP).")
        ss_override = st.text_input("Secondary structure override (SS_OVERRIDE): "
                                    "H for all-helix, blank to use DSSP", "")
        ss_mode = st.selectbox(
            "Secondary structure source (SS_MODE)",
            ["dssp", "tm"],
            help="dssp: DSSP assigns secondary structure (use for GPCRs and other "
                 "multi-helix proteins). tm: the TM range is set to helix and "
                 "everything else to coil, so the juxtamembrane stays flexible "
                 "(single-pass TM and TM-JM peptides). SS_MODE=tm needs the TM "
                 "range or TM-core ranges set below, and uses no DSSP.")
        multi_tm = st.toggle(
            "Multi-pass TM bundle (e.g. 7-TM GPCR)",
            help="Orient on the helix bundle instead of one principal axis. "
                 "Detects every membrane-length helix from DSSP and sets the "
                 "membrane normal to the mean helix axis. Use for GPCRs and "
                 "other multi-pass receptors. Single-pass TM peptides leave "
                 "this off.")
        multi_tm_minlen = ""
        if multi_tm:
            multi_tm_minlen = st.text_input(
                "Minimum helix length counted as a TM helix (residues)", "12")
        nterm_side = st.selectbox(
            "N-terminus side (membrane orientation)",
            ["auto", "up", "down"],
            help="The orientation axis sign is arbitrary, so a multi-pass "
                 "receptor can come out inverted. Set up to put the N-terminus "
                 "on the upper leaflet, down for the lower leaflet. For a GPCR "
                 "the N-terminus is extracellular, so up places the "
                 "extracellular side on the upper leaflet. auto leaves the sign "
                 "as computed.")
        st.caption("Backend tools required: COBY, martinize2, GROMACS, DSSP.")

    workdir = Path(out_dir).expanduser() if out_dir.strip() \
        else Path(tempfile.gettempdir()) / "memble-gui_work"
    workdir.mkdir(parents=True, exist_ok=True)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Protein")
        pep = st.file_uploader("All-atom protein PDB", type=["pdb"])
        tm_range = st.text_input("TM range override (resSeq, e.g. 619:641)", "")
        prebuilt_multi = st.toggle("Use multi-chain assembly as-is "
                                   "(PREBUILT_MULTI; no replication)")
        if prebuilt_multi:
            res_keep = st.text_input("Residue ranges to keep per chain "
                                     "(RES_KEEP, e.g. A:54-103;B:54-103;"
                                     "C:54-103;D:54-103)", "")
            tm_core = st.text_input("TM-core ranges per chain for z-orient "
                                    "(TM_CORE, e.g. A:65-88;B:65-88;C:65-88;"
                                    "D:65-88)", "")
        else:
            res_keep = tm_core = ""

        st.subheader("Composition")
        asym = st.toggle("Asymmetric leaflets")
        if not asym:
            lip = st.text_area("Lipids  name:ratio[:head] (space separated)",
                               "CHOL:1 DLPC:1 PSM:1")
            upper = lower = ""
        else:
            upper = st.text_area("Upper leaflet", "CHOL:1 DLPC:1 PSM:1")
            lower = st.text_area("Lower leaflet", "CHOL:1 DLPC:1 DOPS:1")
            lip = ""

    with col2:
        st.subheader("Geometry and conditions")
        if prebuilt_multi:
            st.caption("Protein count is taken from the multi-chain PDB "
                       "(no N_COPY needed).")
            n_copy = 1
        else:
            n_copy = st.number_input("Protein copies (N_COPY)", 1, 64, 4)
        box_x = st.number_input("Box X nm (0 = auto)", 0.0, 200.0, 0.0)
        box_y = st.number_input("Box Y nm (0 = auto)", 0.0, 200.0, 0.0)
        water_nm = st.number_input("Water per side nm (WATER_NM)", 1.0, 20.0, 2.5)
        memb_thick = st.number_input("Bilayer thickness nm (MEMB_THICK_NM)",
                                     2.0, 8.0, 4.0)
        z_shift = st.number_input("Protein z shift nm (Z_SHIFT): fine tune how "
                                  "deep the protein sits; 0 keeps the TM core at "
                                  "the midplane", -2.0, 2.0, 0.0, 0.1)
        temp = st.number_input("Temperature K (TEMP)", 250, 400, 310)
        salt = st.number_input("NaCl molarity (SALT_M)", 0.0, 1.0, 0.15)
        coby_apl = st.number_input("Area per lipid nm^2 (COBY packing; raise if "
                                   "lipids overlap)", 0.60, 0.85, 0.70, 0.01)

    st.subheader("Peripheral protein (optional)")
    st.caption("Atomistic PDB is coarse-grained and placed AFTER the membrane "
               "is built, on the leaflet you choose, a set gap away. It is "
               "frozen during equilibration and, in production, kept from "
               "wrapping to the other (asymmetric) leaflet by a one-sided "
               "flat-bottom restraint while remaining free to associate.")
    pf = st.file_uploader("Peripheral protein PDB (all-atom)", type=["pdb"],
                          key="partner_pdb")
    partner_pdb = ""
    partner_side = "upper"
    partner_gap = 1.5
    partner_water = 3.0
    if pf is not None:
        ppath = workdir / "partner_aa.pdb"
        ppath.write_bytes(pf.getbuffer())
        partner_pdb = str(ppath)
        partner_side = st.radio("Which leaflet to associate with",
                                ["upper", "lower"], horizontal=True)
        partner_gap = st.number_input("Initial gap: leaflet head to partner "
                                      "edge (nm)", 0.5, 6.0, 1.5, 0.1)
        partner_water = st.number_input("Bulk water beyond partner (nm)",
                                        1.0, 10.0, 3.0, 0.5)

    n_rep = st.number_input("Replicates (N_REP, >1 = orientation-randomized "
                            "ensemble)", 1, 200, 1)

    if st.button("Build system", type="primary"):
        if pep is None or not m3_dir:
            st.error("Upload a protein PDB and set the Martini 3 directory.")
            return
        if prebuilt_multi and not res_keep.strip():
            st.error("PREBUILT_MULTI needs RES_KEEP (residue ranges per chain).")
            return
        pep_path = workdir / "protein_aa.pdb"
        pep_path.write_bytes(pep.getbuffer())

        form = dict(
            pep_aa=str(pep_path), m3_dir=m3_dir, gmx=gmx, dssp=dssp,
            ss_override=ss_override, ss_mode=ss_mode, z_shift=z_shift,
            multi_tm=multi_tm, multi_tm_minlen=multi_tm_minlen,
            nterm_side=nterm_side,
            prebuilt_multi=prebuilt_multi, res_keep=res_keep, tm_core=tm_core,
            asym=asym, lipids=lip, upper=upper, lower=lower,
            n_copy=int(n_copy), box_x=(box_x or None), box_y=(box_y or None),
            water_nm=water_nm, memb_thick_nm=memb_thick, temp=int(temp),
            salt_m=salt, tm_range=tm_range, coby_apl=coby_apl,
            partner_pdb=partner_pdb, partner_side=partner_side,
            partner_gap=partner_gap, partner_water=partner_water,
            n_rep=int(n_rep),
        )
        env = build_env(form, here)
        env["PYTHONUNBUFFERED"] = "1"   # stream helper output live to the GUI
        cmd = build_command(form, here)

        st.code(" ".join("%s=%s" % (k, env[k]) for k in
                ("M3_DIR", "GMX", "LIPIDS", "UPPER", "LOWER", "PARTNER", "N_REP")
                if k in env) + "\n" + " ".join(cmd), language="bash")

        import time
        t0 = time.time()
        log = st.empty()
        status = st.empty()
        buf = []
        proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1,
                                cwd=workdir)
        for line in proc.stdout:
            buf.append(line.rstrip())
            log.code("\n".join(buf[-40:]))
            status.caption("building... %.0f s elapsed" % (time.time() - t0))
        rc = proc.wait()
        status.caption("finished in %.0f s" % (time.time() - t0))
        gro = next(workdir.glob("**/system.gro"), None)
        if rc == 0:
            outloc = gro.parent if gro is not None else workdir
            st.success("Build complete. System saved in: %s" % outloc)
            st.code(
                "cd %s\n"
                "vmd -e view.vmd                 # bonded view: system.psf + trajectory\n"
                "#   chain P = the proteins, protein = all protein, chain M = membrane\n"
                "bash run_md.sh 6.0              # minimize, then 6.1 .. 6.6, then prod\n"
                "grep 'water cushion' build.log  # check the protein has enough water"
                % outloc, language="bash")
        else:
            st.error("Build stopped (exit %d). See the log, including any "
                     "leaflet area mismatch or COBY token issue." % rc)

        if gro is not None:
            _preview(st, gro)


def _preview(st, gro):
    """3D preview of the built coarse-grained system. The structure is wrapped
    into the box per molecule and centered on the membrane, so a correct system
    shows as a continuous bilayer rather than one split across the periodic box.
    Beads are drawn as spheres because a plain structure carries no bonds, and
    water and ions are hidden. For a full bonded view, open view.vmd in VMD."""
    try:
        import py3Dmol
        import streamlit.components.v1 as components
    except Exception as e:
        st.info("3D preview unavailable (%s). Open view.vmd in VMD." % e)
        return
    try:
        pdb, natoms = _gro_to_wrapped_pdb(gro)
    except Exception as e:
        st.info("3D preview skipped (%s). Open view.vmd in VMD." % e)
        return
    if natoms == 0:
        st.info("nothing to preview. Open view.vmd in VMD.")
        return
    view = py3Dmol.view(width=700, height=500)
    view.addModel(pdb, "pdb")
    view.setStyle({"chain": "P"}, {"sphere": {"radius": 1.6, "color": "orange"}})
    view.setStyle({"chain": "U"}, {"sphere": {"radius": 1.1, "color": "lightblue"}})
    view.setStyle({"chain": "L"}, {"sphere": {"radius": 1.1, "color": "palegreen"}})
    view.zoomTo()
    components.html(view._make_html(), height=520)
    st.caption("orange: protein   light blue: outer leaflet   green: inner "
               "leaflet   (water and ions hidden). Wrapped for display; for a "
               "bonded view open view.vmd in VMD.")


_AA3 = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
        "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
        "HSD", "HSE", "HSP", "HID", "HIE", "HIP"}
_SKIP = {"W", "WF", "WAT", "SOL", "ION", "NA", "CL", "NA+", "CL-", "K", "CA",
         "MG", "POT", "SOD", "CLA"}


def _gro_to_wrapped_pdb(gro):
    """Read a gro file, wrap each residue into the box by its centroid, center
    the membrane in z, drop water and ions, and return a PDB string plus the
    atom count. Protein gets chain P, outer leaflet U, inner leaflet L."""
    import numpy as np
    lines = gro.read_text().splitlines()
    n = int(lines[1])
    body = lines[2:2 + n]
    box = [float(v) for v in lines[2 + n].split()[:3]]
    bx, by, bz = box[0], box[1], box[2]
    recs = []
    for ln in body:
        resid = ln[0:5].strip()
        resn = ln[5:10].strip().upper()
        name = ln[10:15].strip()
        x = float(ln[20:28]); y = float(ln[28:36]); z = float(ln[36:44])
        recs.append([resid, resn, name, x, y, z])
    from collections import OrderedDict
    groups = OrderedDict()
    for i, r in enumerate(recs):
        groups.setdefault((r[0], r[1]), []).append(i)
    # membrane midplane from lipid z, by circular mean so a bilayer that
    # straddles the periodic boundary is located correctly (a plain mean of z
    # would point at the gap between the leaflets, not at the membrane).
    lip_z = []
    for (resid, resn), idx in groups.items():
        if resn in _AA3 or resn in _SKIP:
            continue
        for i in idx:
            lip_z.append(recs[i][5])
    if lip_z and bz > 0:
        ang = 2.0 * np.pi * np.array(lip_z) / bz
        m = np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang)))
        if m < 0:
            m += 2.0 * np.pi
        midz = m / (2.0 * np.pi) * bz
    else:
        midz = bz / 2.0
    shift_z = bz / 2.0 - midz
    for (resid, resn), idx in groups.items():
        for i in idx:
            recs[i][5] += shift_z
        cx = np.mean([recs[i][3] for i in idx])
        cy = np.mean([recs[i][4] for i in idx])
        cz = np.mean([recs[i][5] for i in idx])
        sx = -np.floor(cx / bx) * bx if bx > 0 else 0.0
        sy = -np.floor(cy / by) * by if by > 0 else 0.0
        sz = -np.floor(cz / bz) * bz if bz > 0 else 0.0
        for i in idx:
            recs[i][3] += sx; recs[i][4] += sy; recs[i][5] += sz
    out = []
    serial = 1
    for (resid, resn), idx in groups.items():
        if resn in _SKIP:
            continue
        if resn in _AA3:
            chain = "P"
        else:
            cz = np.mean([recs[i][5] for i in idx])
            chain = "U" if cz >= bz / 2.0 else "L"
        try:
            rid = int(resid) % 10000
        except ValueError:
            rid = 1
        for i in idx:
            name = recs[i][2]; x = recs[i][3]; y = recs[i][4]; z = recs[i][5]
            out.append("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00" %
                       (serial % 100000, name[:4], resn[:3], chain, rid,
                        x * 10.0, y * 10.0, z * 10.0))
            serial += 1
    out.append("END")
    return "\n".join(out), serial - 1


if __name__ == "__main__":
    main()
