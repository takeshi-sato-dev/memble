"""How often a lipid changes leaflet, against bilayers whose answer is known.

count_crossings.py counts the crossings themselves, where leaflet_relax.py
counts how many molecules each leaflet holds. The difference matters for a
species that keeps moving after the composition has settled: its net stops
changing while its crossings do not.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

md = pytest.importorskip("mdtraj")

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "count_crossings.py"

STEROL = [("ROH", 0.0), ("R1", 0.3), ("R2", 0.6), ("C1", 0.9), ("C2", 1.2)]
LIPID = [("PO4", 0.0), ("C1A", 0.5), ("C2A", 1.0), ("C1B", 0.5), ("C2B", 1.0)]


def _traj(tmp_path, mols, sides, head_z=2.0):
    """mols: [(resname, beads)].  sides: [n_frames][n_mols] of +1 (upper)/-1.

    A molecule is placed with its head bead at +head_z and its tails below when
    it is in the upper leaflet, and the mirror image when it is in the lower
    one, which is what the assignment rule reads.
    """
    top = md.Topology()
    ch = top.add_chain()
    for name, beads in mols:
        r = top.add_residue(name, ch)
        for bead, _ in beads:
            top.add_atom(bead, md.element.carbon, r)
    frames = []
    for row in sides:
        xyz = []
        for i, ((name, beads), s) in enumerate(zip(mols, row)):
            for _, dz in beads:
                xyz.append([0.5 * i, 0.0, s * head_z - s * dz])
        frames.append(xyz)
    t = md.Trajectory(np.array(frames, dtype=np.float32), top)
    n = len(sides)
    t.unitcell_lengths = np.tile([5.0, 5.0, 10.0], (n, 1)).astype("f4")
    t.unitcell_angles = np.tile([90.0, 90.0, 90.0], (n, 1)).astype("f4")
    t.time = (np.arange(n) * 1000.0).astype("f4")        # 1 ns per frame
    gro, xtc = tmp_path / "s.gro", tmp_path / "s.xtc"
    t[0].save_gro(str(gro))
    t.save_xtc(str(xtc))
    return gro, xtc


def _run(gro, xtc, tmp_path, *extra):
    out = tmp_path / "c.json"
    p = subprocess.run([sys.executable, str(SCRIPT), "--gro", str(gro),
                        "--xtc", str(xtc), "--json", str(out), *extra],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return json.load(open(out))


def test_a_molecule_that_never_crosses_is_counted_as_never_crossing(tmp_path):
    mols = [("CHOL", STEROL), ("DLPC", LIPID)]
    sides = [[+1, -1]] * 8
    d = _run(*_traj(tmp_path, mols, sides), tmp_path)
    assert d["species"]["CHOL"]["crossings_total"] == 0
    assert d["species"]["DLPC"]["crossings_total"] == 0
    assert d["species"]["CHOL"]["net_change"] == 0


def test_one_crossing_that_stays_is_one_crossing(tmp_path):
    mols = [("CHOL", STEROL)]
    sides = [[+1]] * 4 + [[-1]] * 4
    d = _run(*_traj(tmp_path, mols, sides), tmp_path)
    c = d["species"]["CHOL"]
    assert c["crossings_total"] == 1
    assert c["crossings_down"] == 1 and c["crossings_up"] == 0
    assert c["net_change"] == -1


def test_a_molecule_that_goes_and_returns_moves_nothing_and_crosses_twice(tmp_path):
    """This is the case the net cannot see. Two crossings, no change of
    composition, and it is what a species at equilibrium does all run."""
    mols = [("CHOL", STEROL)]
    sides = [[+1]] * 3 + [[-1]] * 3 + [[+1]] * 3
    d = _run(*_traj(tmp_path, mols, sides), tmp_path)
    c = d["species"]["CHOL"]
    assert c["crossings_total"] == 2
    assert c["crossings_up"] == 1 and c["crossings_down"] == 1
    assert c["net_change"] == 0


def test_a_one_frame_flicker_is_not_a_crossing(tmp_path):
    """A molecule whose head and body sit almost level answers differently from
    one frame to the next. Those flickers are not crossings, and counting them
    would turn a still membrane into a busy one."""
    mols = [("CHOL", STEROL)]
    sides = [[+1], [+1], [-1], [+1], [+1], [+1], [-1], [+1]]
    d = _run(*_traj(tmp_path, mols, sides), tmp_path)
    assert d["species"]["CHOL"]["crossings_total"] == 0
    # with the filter switched off every change is counted
    d1 = _run(*_traj(tmp_path, mols, sides), tmp_path, "--dwell", "1")
    assert d1["species"]["CHOL"]["crossings_total"] == 4


def test_the_rate_is_reported_before_and_after_the_composition_settles(tmp_path):
    """The claim of the article is that the sterol keeps moving after the net
    has stopped changing, so the two windows are counted apart."""
    mols = [("CHOL", STEROL)]
    # ns 0-3: one crossing that changes the composition.  ns 4-13: it goes back
    # and forth, so the net is flat while the crossings continue.
    sides = [[+1]] * 2 + [[-1]] * 2 + [[+1]] * 2 + [[-1]] * 2 + [[+1]] * 2 + [[-1]] * 2
    gro, xtc = _traj(tmp_path, mols, sides)
    d = _run(gro, xtc, tmp_path, "--settled-from-ns", "4")
    c = d["species"]["CHOL"]
    assert c["before_settled"]["up"] + c["before_settled"]["down"] == 1
    assert c["after_settled"]["up"] + c["after_settled"]["down"] == 4
    assert c["after_settled"]["crossings_per_us"] > 0


def test_the_imbalance_of_a_curve_point_comes_from_its_counts(tmp_path):
    """The imbalance of a point is the difference between the two phospholipid
    numbers over their sum, and those two numbers are read from the topology of
    that build. Writing the imbalance as a function of delta gives the wrong
    answer at two of the six points, because the packing rounds the ratio of the
    minor species and puts one more molecule in the lower leaflet there."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("bigread", str(ROOT / "bigread.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for delta, f in ((-20, -11.76), (-12, -6.23), (-6, -1.83),
                     (0, 2.94), (6, 7.35), (12, 11.76)):
        up, lo = m.COUNTS[delta]
        got = m.imbalance_from_counts(up, lo)
        assert abs(got - f) < 0.01, (delta, got, f)
    # the packing rounds the minor species, so the total is not the same at
    # every point and the imbalance cannot be written as a function of delta
    totals = sorted(set(sum(v) for v in m.COUNTS.values()))
    assert totals == [272, 273], totals


def test_the_curve_in_bigread_matches_the_fit(tmp_path):
    """bigread.py carries the quadratic as three constants, and those three have
    to be the ones curve_stats.py returns from the six points."""
    import importlib.util
    import numpy as np
    spec = importlib.util.spec_from_file_location("bigread", str(ROOT / "bigread.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    spec2 = importlib.util.spec_from_file_location("cs", str(ROOT / "curve_stats.py"))
    cs = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(cs)
    q = np.polyfit(cs.F, cs.C, 2)
    assert abs(m.A2 - q[0]) < 1e-5, (m.A2, q[0])
    assert abs(m.A1 - q[1]) < 1e-4, (m.A1, q[1])
    assert abs(m.A0 - q[2]) < 1e-3, (m.A0, q[2])
