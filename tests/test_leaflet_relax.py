"""The leaflet a lipid belongs to, and the plane the two leaflets are measured from.

leaflet_relax.py reports how many molecules of each species changed leaflet
during a run, and how thick each leaflet is. Both rest on one plane. These tests
build bilayers whose leaflets are known by construction and check that the plane
the script finds is the plane that separates them.
"""
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
    """A flat bilayer. upper and lower are lists of (residue name, bead list)."""
    top = md.Topology()
    ch = top.add_chain()
    xyz = []
    # dz runs from the head bead toward the tail, so the upper leaflet has its
    # head at +head_z and its tails below it, and the lower leaflet the reverse.
    for side, sign, leaflet in (("u", +1.0, upper), ("l", -1.0, lower)):
        for i, (name, beads) in enumerate(leaflet):
            r = top.add_residue(name, ch)
            for bead, dz in beads:
                top.add_atom(bead, md.element.carbon, r)
                xyz.append([0.5 * i, 0.0 if side == "u" else 1.0,
                            sign * head_z - sign * dz])
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
