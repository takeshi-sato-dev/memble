"""The baseline the moved numbers of Table 3 and of the response curve use.

summarise_counts.py collects the runs of compare_counts.sh, and
curve_report.py collects the points of run_curve.sh. Both report how far
cholesterol moved, and both take that move against the numbers the build gave
each leaflet.

The first frame of a production run is not those numbers. The equilibration
has already moved cholesterol by the time that frame is written, and in a
bilayer that holds a protein a few molecules are assigned to the leaflet they
did not come from, which puts a different error on each seed of one build.
The tests below build a directory whose builds and whose runs are known, and
give the first frame of every run an offset from the build.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARISE = ROOT / "summarise_counts.py"
CURVE = ROOT / "curve_report.py"

BUILD_JSON = json.dumps({"upper": "CHOL:1 DLPC:1 PSM:1",
                         "lower": "CHOL:1 DLPC:1 DOPS:1 POP2_45:0.1"})


def _build(work, chol_upper, dlpc_upper, chol_lower, dlpc_lower):
    os.makedirs(work, exist_ok=True)
    Path(work, "memble_build.json").write_text(BUILD_JSON)
    Path(work, "system.top").write_text(
        "[ molecules ]\n; name number\nmolecule_0 1\n"
        "CHOL %d\nDLPC %d\nPSM 70\n"
        "CHOL %d\nDLPC %d\nDOPS 62\nPOP2_45 7\n"
        "W 16023\nNA 264\nCL 174\n"
        % (chol_upper, dlpc_upper, chol_lower, dlpc_lower))


def _run_json(chol_upper_built, first_frame_offset, settled_move, n=50):
    """A run whose first frame is off the build and whose sterol then moves."""
    cut = int(0.6 * n)
    upper = [chol_upper_built + first_frame_offset] * cut \
        + [chol_upper_built + settled_move] * (n - cut)
    rec = {"time_ps": [float(i * 1000) for i in range(n)],
           "counts": {"CHOL": {"upper": upper,
                               "lower": [134 - x for x in upper]}},
           "between_leaflets": {
               "CHOL": {"start": upper[0],
                        "settled": float(chol_upper_built + settled_move),
                        "moved": float(settled_move - first_frame_offset)}},
           "leaflet_thickness_nm": {"upper": {"settled": 2.0},
                                    "lower": {"settled": 2.1}},
           "tail_order": {"upper": {"settled": 0.30},
                          "lower": {"settled": 0.25}}}
    return rec


def test_the_move_of_each_arm_is_taken_against_its_build(tmp_path):
    d = tmp_path / "counts"
    d.mkdir()
    _build(str(d / "table" / "table_work"), 71, 70, 63, 63)
    # Three seeds of one build. Each starts its production run at a different
    # place, and all three settle at 80. The move is 9 for every seed.
    for seed, off in enumerate((+2, -1, +3), 1):
        (d / ("relax_table_s%d.json" % seed)).write_text(
            json.dumps(_run_json(71, off, 9)))

    p = subprocess.run([sys.executable, str(SUMMARISE), str(d)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    got = json.load(open(d / "counts_summary.json"))
    assert got["table"]["moved"][0] == 9.0        # the mean over the seeds
    assert got["table"]["moved"][1] == 0.0        # and no scatter between them
    assert got["table"]["baseline_from"] == ["system.top"]


def test_an_arm_whose_build_is_gone_says_so_rather_than_reporting_the_frame(tmp_path):
    d = tmp_path / "counts"
    d.mkdir()
    for seed, off in enumerate((+2, -1, +3), 1):
        (d / ("relax_eqn_s%d.json" % seed)).write_text(
            json.dumps(_run_json(67, off, 9)))
    p = subprocess.run([sys.executable, str(SUMMARISE), str(d)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert "the build of the eqn arm was not found" in p.stdout
    got = json.load(open(d / "counts_summary.json"))
    assert got["eqn"]["moved"][0] != 9.0          # the first frame, and wrong
    assert got["eqn"]["baseline_from"] == []


def test_the_curve_is_drawn_against_the_numbers_the_builds_were_given(tmp_path):
    d = tmp_path / "curve"
    d.mkdir()
    # DLPC is moved between the leaflets, so the two leaflets hold
    # 50+70 against 83+62+7, 70+70 against 63+62+7, and 90+70 against 43+62+7,
    # which are differences of -32, +8 and +48.
    for tag, dl_u, dl_l in (("d-20", 50, 83), ("d0", 70, 63), ("d20", 90, 43)):
        _build(str(d / tag / ("%s_work" % tag)), 71, dl_u, 63, dl_l)
        (d / ("relax_%s.json" % tag)).write_text(json.dumps(_run_json(71, +2, 8)))

    out = d / "curve.csv"
    p = subprocess.run([sys.executable, str(CURVE), str(d), "--csv", str(out)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr

    import csv
    rows = list(csv.DictReader(open(out)))
    assert [r["point"] for r in rows] == ["d-20", "d0", "d20"]
    assert [int(r["pl_difference"]) for r in rows] == [-32, +8, +48]
    # every point moved 8 molecules against its build, and started its run 2
    # molecules away from that build
    assert all(float(r["moved"]) == 8.0 for r in rows)
    assert all(int(r["chol_upper_first_frame"]) == 73 for r in rows)
    assert all(r["baseline_from"] == "system.top" for r in rows)
