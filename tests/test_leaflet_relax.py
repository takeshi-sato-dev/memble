"""The leaflet a lipid belongs to, and the plane its thickness is measured from.

leaflet_relax.py reports how many molecules of each species changed leaflet
during a run, and how thick each leaflet is. These tests build bilayers whose
leaflets are known by construction, and hand the script the cases that a
bilayer with a protein in it produces: leaflets of unequal composition, a lipid
whose head beads outnumber its tail beads, a sterol that carries no phosphate,
and a leaflet that is not flat.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

md = pytest.importorskip("mdtraj")

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "leaflet_relax.py"

# A short lipid, five beads, and a long one, ten. The head bead of each is PO4.
SHORT = [("PO4", 0.0), ("C1A", 0.5), ("C2A", 1.0), ("C1B", 0.5), ("C2B", 1.0)]
LONG = [("PO4", 0.0), ("C1A", 0.4), ("C2A", 0.8), ("C3A", 1.2), ("C1B", 0.4),
        ("C2B", 0.8), ("C3B", 1.2), ("C4B", 1.5), ("C5B", 1.7), ("C6B", 1.9)]
# A sterol carries no phosphate: its head bead is the hydroxyl.
STEROL = [("ROH", 0.0), ("R1", 0.3), ("R2", 0.6), ("C1", 0.9), ("C2", 1.2)]


def _bilayer(tmp_path, upper, lower, head_z=2.0, n_frames=3):
    """A bilayer.

    upper and lower are lists of (residue name, bead list), or of (residue
    name, bead list, offset) where offset moves that one molecule along z, so
    that a leaflet can be given a dimple.
    """
    top = md.Topology()
    ch = top.add_chain()
    xyz = []
    # dz runs from the head bead toward the tail, so the upper leaflet has its
    # head at +head_z and its tails below it, and the lower leaflet the reverse.
    for side, sign, leaflet in (("u", +1.0, upper), ("l", -1.0, lower)):
        for i, entry in enumerate(leaflet):
            name, beads = entry[0], entry[1]
            off = entry[2] if len(entry) > 2 else 0.0
            r = top.add_residue(name, ch)
            for bead, dz in beads:
                top.add_atom(bead, md.element.carbon, r)
                xyz.append([0.5 * i, 0.0 if side == "u" else 1.0,
                            sign * head_z - sign * dz + off])
    frames = np.repeat(np.array(xyz, dtype=np.float32)[None], n_frames, axis=0)
    t = md.Trajectory(frames, top)
    t.unitcell_lengths = np.tile([5.0, 5.0, 10.0], (n_frames, 1)).astype("f4")
    t.unitcell_angles = np.tile([90.0, 90.0, 90.0], (n_frames, 1)).astype("f4")
    t.time = np.arange(n_frames) * 100.0
    gro, xtc = tmp_path / "s.gro", tmp_path / "s.xtc"
    t[0].save_gro(str(gro))
    t.save_xtc(str(xtc))
    return gro, xtc


def _run(gro, xtc, tmp_path):
    out = tmp_path / "relax.json"
    p = subprocess.run([sys.executable, str(SCRIPT), "--gro", str(gro),
                        "--xtc", str(xtc), "--json", str(out)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    import json
    return json.load(open(out))


def test_symmetric_bilayer_is_split_in_half(tmp_path):
    up = [("DLPC", SHORT)] * 10
    lo = [("DLPC", SHORT)] * 10
    r = _run(*_bilayer(tmp_path, up, lo), tmp_path=tmp_path)
    assert r["counts"]["DLPC"]["upper"][0] == 10
    assert r["counts"]["DLPC"]["lower"][0] == 10


def test_the_plane_does_not_follow_the_heavier_leaflet(tmp_path):
    """The two leaflets hold different numbers of beads.

    The mean z of every lipid bead then sits toward the leaflet that holds more
    of them, and a leaflet measured from that mean is reported thicker than it
    is by exactly the offset. Both leaflets here are 2.0 nm thick.
    """
    up = [("DLPC", SHORT)] * 10          # 50 beads
    lo = [("DPPC", LONG)] * 10           # 100 beads
    r = _run(*_bilayer(tmp_path, up, lo, head_z=2.0), tmp_path=tmp_path)
    assert r["counts"]["DLPC"]["upper"][0] == 10
    assert r["counts"]["DPPC"]["upper"][0] == 0
    th = r["leaflet_thickness_nm"]
    assert th["upper"]["settled"] == pytest.approx(2.0, abs=0.02)
    assert th["lower"]["settled"] == pytest.approx(2.0, abs=0.02)


def test_a_sterol_is_placed_by_its_hydroxyl(tmp_path):
    """A sterol carries no phosphate, and its head bead is ROH."""
    up = [("DLPC", SHORT)] * 8 + [("CHOL", STEROL)] * 4
    lo = [("DLPC", SHORT)] * 8 + [("CHOL", STEROL)] * 4
    r = _run(*_bilayer(tmp_path, up, lo), tmp_path=tmp_path)
    assert r["counts"]["CHOL"]["upper"][0] == 4
    assert r["counts"]["CHOL"]["lower"][0] == 4


def test_a_lipid_whose_beads_straddle_the_plane_is_placed_by_its_head(tmp_path):
    """A long headgroup over short tails.

    Every bead of this lipid but the head sits within 0.3 nm of the head, and
    the tails reach 1.6 nm past it, so the mean over the whole molecule falls on
    the far side of the plane while the head bead does not. The molecule belongs
    to the leaflet its head bead is in.
    """
    weird = [("PO4", 0.0), ("P1", -0.1), ("P2", -0.2),
             ("C1A", 1.6), ("C2A", 1.7), ("C3A", 1.8), ("C1B", 1.6),
             ("C2B", 1.7)]
    up = [("DLPC", SHORT)] * 10
    lo = [("DLPC", SHORT)] * 9 + [("POP2_", weird)]
    # the mean over the whole molecule lands above the plane, its head bead below
    r = _run(*_bilayer(tmp_path, up, lo, head_z=0.9), tmp_path=tmp_path)
    assert r["counts"]["POP2_"]["lower"][0] == 1
    assert r["counts"]["POP2_"]["upper"][0] == 0
    # and it stays there for every frame, so nothing is reported as moved
    assert r["between_leaflets"]["POP2_"]["moved"] == pytest.approx(0.0)


def test_a_lipid_in_a_dimple_stays_in_its_own_leaflet(tmp_path):
    """A bilayer that holds a protein is not flat.

    Two lipids of the upper leaflet are pushed 2.4 nm down, far enough that
    their head beads fall below the midpoint of the two leaflets, and their
    tails go down with them. Each still points out of the leaflet it belongs
    to, and belongs to the upper leaflet.
    """
    up = [("DLPC", SHORT)] * 8 + [("DLPC", SHORT, -2.4)] * 2
    lo = [("DPPC", LONG)] * 10
    r = _run(*_bilayer(tmp_path, up, lo, head_z=2.0), tmp_path=tmp_path)
    assert r["counts"]["DLPC"]["upper"][0] == 10
    assert r["counts"]["DLPC"]["lower"][0] == 0
    assert r["counts"]["DPPC"]["upper"][0] == 0


# --- the baseline the moved numbers are measured against --------------------
#
# A trajectory's first frame is not the composition the system was built with:
# the equilibration has already moved cholesterol by the time that frame is
# written. The build is read from system.top, whose [ molecules ] section holds
# the upper leaflet first and the lower leaflet second, and from
# memble_build.json, which names the species of each leaflet.

import importlib.util

_spec = importlib.util.spec_from_file_location("leaflet_relax", SCRIPT)
_lr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lr)


def _write_build(tmp_path, molecules, upper="", lower=""):
    top = tmp_path / "system.top"
    top.write_text("#include \"martini.itp\"\n\n[ molecules ]\n; name number\n"
                   + "".join("%s %d\n" % (n, k) for n, k in molecules))
    bj = tmp_path / "memble_build.json"
    bj.write_text(json.dumps({"upper": upper, "lower": lower}))
    return str(top), str(bj)


def test_the_build_is_read_from_the_topology_leaflet_by_leaflet(tmp_path):
    top, bj = _write_build(
        tmp_path,
        [("molecule_0", 1), ("CHOL", 71), ("DLPC", 50), ("PSM", 70),
         ("CHOL", 63), ("DLPC", 83), ("DOPS", 62), ("POP2_45", 8),
         ("W", 16023), ("NA", 264), ("CL", 174)],
        upper="CHOL:1.014 DLPC:0.714 PSM:1",
        lower="CHOL:1.016 DLPC:1.339 DOPS:1 POP2_45:0.113")
    got = _lr.built_counts(top, bj)
    assert got["upper"] == {"CHOL": 71, "DLPC": 50, "PSM": 70}
    assert got["lower"] == {"CHOL": 63, "DLPC": 83, "DOPS": 62, "POP2_45": 8}


def test_a_build_given_one_composition_divides_its_lipids_in_the_middle(tmp_path):
    top, bj = _write_build(
        tmp_path,
        [("molecule_0", 1), ("CHOL", 47), ("DLPC", 47), ("PSM", 47),
         ("CHOL", 47), ("DLPC", 46), ("PSM", 46), ("W", 11065)])
    got = _lr.built_counts(top, bj)
    assert got["upper"] == {"CHOL": 47, "DLPC": 47, "PSM": 47}
    assert got["lower"] == {"CHOL": 47, "DLPC": 46, "PSM": 46}


def test_a_missing_topology_returns_nothing_rather_than_a_wrong_baseline(tmp_path):
    assert _lr.built_counts(str(tmp_path / "no.top"),
                            str(tmp_path / "no.json")) is None
