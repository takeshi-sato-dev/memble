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
import re
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


def test_the_repeats_of_a_point_are_not_read_as_points_of_the_curve(tmp_path):
    """A point of the curve is named for its delta and nothing else. Seeds of a
    point carry the seed and the run length as well (relax_d0_s2_250ns.json),
    and they share the directory. Leave them out and say so, rather than
    stopping on the first name that does not parse."""
    d = tmp_path / "curve"
    d.mkdir()
    for tag, dl_u, dl_l in (("d-20", 50, 83), ("d0", 70, 63)):
        _build(str(d / tag / ("%s_work" % tag)), 71, dl_u, 63, dl_l)
        (d / ("relax_%s.json" % tag)).write_text(json.dumps(_run_json(71, +2, 8)))
    for name in ("relax_d0_s2_250ns.json", "relax_d0_s2_500ns.json",
                 "relax_d-20_s3.json"):
        (d / name).write_text(json.dumps(_run_json(71, +2, 12)))

    out = d / "curve.csv"
    p = subprocess.run([sys.executable, str(CURVE), str(d), "--csv", str(out)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert "not points of the curve" in p.stdout
    assert "relax_d0_s2_250ns.json" in p.stdout

    import csv
    rows = list(csv.DictReader(open(out)))
    assert [r["point"] for r in rows] == ["d-20", "d0"]
    assert all(float(r["moved"]) == 8.0 for r in rows)


SOLVE = ROOT / "solve_target.py"


def _curve_dir(tmp_path, rows):
    """rows: (delta, pl_upper, pl_lower, settled share %, runs)."""
    d = tmp_path / "curve"
    d.mkdir()
    for delta, pu, pl, share, runs in rows:
        for k in range(runs):
            tot = 134.0
            up = share / 100.0 * tot
            rec = {"counts": {"CHOL": {"upper": [up] * 50,
                                       "lower": [tot - up] * 50}},
                   "time_ps": [i * 1000.0 for i in range(50)],
                   "between_leaflets": {"CHOL": {"built": 71, "first_frame": 71,
                                                 "settled": up, "moved": up - 71,
                                                 "baseline": "system.top"}},
                   "built_counts": {"upper": {"CHOL": 71, "DLPC": pu - 70,
                                              "PSM": 70},
                                    "lower": {"CHOL": 63, "DLPC": pl - 69,
                                              "DOPS": 62, "POP2_45": 7}}}
            name = ("relax_d%d.json" % delta if k == 0
                    else "relax_d%d_s%d_250ns.json" % (delta, k + 1))
            (d / name).write_text(json.dumps(rec))
    return d


CURVE_ROWS = [(-20, 120, 152, 64.78, 3), (-12, 128, 145, 63.20, 1),
              (-6, 134, 139, 61.17, 1), (0, 140, 132, 56.86, 3),
              (6, 146, 126, 50.26, 1), (12, 152, 120, 46.36, 1)]


def test_a_target_returns_the_two_phospholipid_numbers(tmp_path):
    """The point of the curve is to be read backwards. A share the run must
    hold goes in, and the numbers a build has to be given come out."""
    d = _curve_dir(tmp_path, CURVE_ROWS)
    p = subprocess.run([sys.executable, str(SOLVE), str(d), "--target", "60",
                        "--emit-env"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr + p.stdout
    env = dict(line.split("=", 1) for line in p.stdout.splitlines()
               if re.match(r"^[A-Z_]+=", line))
    up, lo = int(env["PL_UPPER"]), int(env["PL_LOWER"])
    assert up + lo == 272                      # the curve's own total
    assert lo > up                             # 60% needs the upper leaflet short
    assert abs(float(env["IMBALANCE_PCT"])) < 3.0


def test_the_numbers_scale_to_the_box_the_system_will_use(tmp_path):
    """The axis is a fraction, so a target asked of a larger system returns
    the same fraction and a larger pair of numbers."""
    d = _curve_dir(tmp_path, CURVE_ROWS)
    out = {}
    for total in (272, 1937):
        p = subprocess.run([sys.executable, str(SOLVE), str(d), "--target", "60",
                            "--total-pl", str(total), "--emit-env"],
                           capture_output=True, text=True)
        assert p.returncode == 0, p.stderr + p.stdout
        env = dict(line.split("=", 1) for line in p.stdout.splitlines()
                   if re.match(r"^[A-Z_]+=", line))
        out[total] = (int(env["PL_UPPER"]), int(env["PL_LOWER"]),
                      float(env["IMBALANCE_PCT"]))
    assert out[272][2] == out[1937][2]
    assert sum(out[1937][:2]) == 1937
    assert out[1937][0] > out[272][0]


def test_a_target_the_curve_does_not_reach_is_refused(tmp_path):
    """Extrapolating a fitted quadratic past its data returns a number with no
    measurement behind it. Refuse, and say what the curve does reach."""
    d = _curve_dir(tmp_path, CURVE_ROWS)
    p = subprocess.run([sys.executable, str(SOLVE), str(d), "--target", "85"],
                       capture_output=True, text=True)
    assert p.returncode == 2
    assert "does not reach" in p.stdout
    assert "46" in p.stdout and "64" in p.stdout


def test_the_repeats_of_a_point_are_averaged_into_that_point(tmp_path):
    """Three runs of one build are one point of the curve, not three."""
    d = _curve_dir(tmp_path, CURVE_ROWS)
    p = subprocess.run([sys.executable, str(SOLVE), str(d), "--target", "60"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr + p.stdout
    assert "6 points" in p.stdout


def _app():
    import importlib.util
    spec = importlib.util.spec_from_file_location("app", str(ROOT / "app.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # streamlit is imported inside pages
    return mod


def test_the_gui_reads_the_curve_backwards_like_the_command_line(tmp_path):
    """Everything the command line does is on a page of the GUI, and the two
    return the same numbers, because the page calls the same code."""
    d = _curve_dir(tmp_path, CURVE_ROWS)
    app = _app()
    pts = app.curve_points(str(d))
    assert len(pts) == 6
    gui = app.curve_solution(pts, 60.0, 1937)
    assert gui["ok"]

    p = subprocess.run([sys.executable, str(SOLVE), str(d), "--target", "60",
                        "--total-pl", "1937", "--emit-env"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr + p.stdout
    env = dict(line.split("=", 1) for line in p.stdout.splitlines()
               if re.match(r"^[A-Z_]+=", line))
    assert gui["pl_upper"] == int(env["PL_UPPER"])
    assert gui["pl_lower"] == int(env["PL_LOWER"])
    assert abs(gui["imbalance"] - float(env["IMBALANCE_PCT"])) < 1e-3


def test_the_gui_refuses_a_target_the_curve_does_not_reach(tmp_path):
    """The refusal is a property of the measurement, not of the front end."""
    d = _curve_dir(tmp_path, CURVE_ROWS)
    app = _app()
    sol = app.curve_solution(app.curve_points(str(d)), 85.0)
    assert sol["ok"] is False
    assert "does not reach" in sol["reason"]
    assert "pl_upper" not in sol


def test_the_gui_states_both_protocols_as_commands():
    """A protocol whose runs take hours is handed over as a command line, and
    the command has to name the steps in the order the protocol gives them."""
    app = _app()
    a = "\n".join(app.protocol_a_commands("/w", "prod.xtc", "/repo"))
    assert "leaflet_relax.py" in a and "run_settled.sh" in a
    assert a.index("leaflet_relax.py") < a.index("run_settled.sh")
    b = "\n".join(app.protocol_b_commands("/out", "/repo"))
    assert "run_curve.sh" in b and "solve_target.py" in b
    assert b.index("run_curve.sh") < b.index("solve_target.py")


def test_the_gui_shows_what_a_run_did_against_the_build(tmp_path):
    """The settled table is read against system.top, not against the first
    frame, and the page says which."""
    d = tmp_path / "r.json"
    d.write_text(json.dumps({
        "between_leaflets": {
            "CHOL": {"built": 71, "first_frame": 73, "settled": 79.0,
                     "moved": 8.0, "baseline": "system.top"},
            "DLPC": {"built": 70, "first_frame": 70, "settled": 70.0,
                     "moved": 0.0, "baseline": "system.top"}}}))
    app = _app()
    rows = app.settled_rows(str(d))
    chol = [r for r in rows if r["species"] == "CHOL"][0]
    assert chol["built"] == 71
    assert chol["moved"] == 8.0
    assert chol["baseline"] == "system.top"
