"""Tests for the verification stage.

Each of the eight properties that memble measures is exercised twice: on a
system that holds it and on a system built to break it. A check that cannot
fail is not a check, so the failing case matters as much as the passing one.

The systems are written here rather than built, so the tests run in a second
and need neither GROMACS nor COBY.
"""

import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VERIFY = os.path.join(ROOT, "verify_system.py")
EQUIL = os.path.join(ROOT, "check_equilibration.py")

LIPIDS = "CHOL POPC"


# --------------------------------------------------------------- builders

def _gro_line(resid, resname, atom, n, x, y, z):
    return "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (
        resid % 100000, resname, atom, n % 100000, x, y, z)


def write_system(d, n_upper=40, n_lower=40, tm_offset=0.0, box=(6.0, 6.0, 12.0),
                 water_span=(8.6, 11.5), overlap=False, na=6, cl=6,
                 n_res=20, seed=3):
    """Write a small bilayer with a transmembrane helix through it.

    Every argument is a way of breaking one of the eight properties.
    """
    import random
    rng = random.Random(seed)
    rows, rid = [], 0
    mid = 0.5 * box[2]

    for i in range(n_res):
        rid += 1
        z = mid - 1.0 + i * (2.0 / max(n_res - 1, 1)) + tm_offset
        rows.append(_gro_line(rid, "ALA", "BB", len(rows) + 1, 0.08, 0.08, z))

    # A grid with a small jitter, so that no two lipids land on one point and
    # the overlap check has nothing to find in a sound system.
    for side, count in (("up", n_upper), ("lo", n_lower)):
        head = mid + 1.9 if side == "up" else mid - 1.9
        tail = mid + 0.7 if side == "up" else mid - 0.7
        ncol = int(count ** 0.5) + 1
        step_x, step_y = box[0] / ncol, box[1] / ncol
        for i in range(count):
            rid += 1
            name = "CHOL" if i % 2 else "POPC"
            x = (i % ncol + 0.5) * step_x + rng.uniform(-0.05, 0.05) * step_x
            y = (i // ncol + 0.5) * step_y + rng.uniform(-0.05, 0.05) * step_y
            hb = "ROH" if name == "CHOL" else "PO4"
            rows.append(_gro_line(rid, name, hb, len(rows) + 1, x, y, head))
            rows.append(_gro_line(rid, name, "C1A", len(rows) + 1, x, y, tail))

    if overlap:
        rid += 1
        x0 = float(rows[-2][20:28])
        y0 = float(rows[-2][28:36])
        z0 = float(rows[-2][36:44])
        rows.append(_gro_line(rid, "POPC", "PO4", len(rows) + 1,
                              x0 + 0.01, y0, z0))
        rows.append(_gro_line(rid, "POPC", "C1A", len(rows) + 1,
                              x0 + 0.01, y0, z0 - 1.2))
        n_upper += 1

    lo, hi = water_span
    for i in range(120):
        rid += 1
        z = rng.uniform(hi, box[2] - 0.1) if i % 2 else rng.uniform(0.1, lo - 8.6 + 1.0)
        rows.append(_gro_line(rid, "W", "W", len(rows) + 1,
                              (i % 11 + 0.5) * box[0] / 11,
                              (i // 11 % 11 + 0.5) * box[1] / 11, z))
    k = 0
    for name, count in (("NA", na), ("CL", cl)):
        for _ in range(count):
            rid += 1; k += 1
            rows.append(_gro_line(rid, name, name, len(rows) + 1,
                                  (k % 5 + 0.5) * box[0] / 5,
                                  (k // 5 % 5 + 0.5) * box[1] / 5,
                                  1.2 + 0.06 * (k % 6)))

    gro = os.path.join(d, "system.gro")
    with open(gro, "w") as fh:
        fh.write("test\n%d\n" % len(rows))
        fh.write("\n".join(rows))
        fh.write("\n%.5f %.5f %.5f\n" % box)

    n_chol = sum(1 for i in range(n_upper) if i % 2) + \
             sum(1 for i in range(n_lower) if i % 2)
    n_popc = (n_upper + n_lower) - n_chol
    top = os.path.join(d, "system.top")
    with open(top, "w") as fh:
        fh.write("[ molecules ]\nProtein 1\nCHOL %d\nPOPC %d\nW 120\nNA %d\nCL %d\n"
                 % (n_chol, n_popc, na, cl))

    def itp(name, beads, charges):
        s = "[ moleculetype ]\n%s 1\n\n[ atoms ]\n" % name
        for i, (b, q) in enumerate(zip(beads, charges), 1):
            s += "%d P4 1 %s %s %d %.3f 72.0\n" % (i, name, b, i, q)
        return s + "\n"

    with open(os.path.join(d, "all.itp"), "w") as fh:
        fh.write(itp("Protein", ["BB"] * n_res, [0.0] * n_res))
        fh.write(itp("CHOL", ["ROH", "C1A"], [0, 0]))
        fh.write(itp("POPC", ["PO4", "C1A"], [0, 0]))
        fh.write(itp("W", ["W"], [0]))
        fh.write(itp("NA", ["NA"], [1.0]))
        fh.write(itp("CL", ["CL"], [-1.0]))
    return gro, top


def run_verify(d, **kw):
    """Run the verification stage and return the report it wrote."""
    cmd = [sys.executable, VERIFY,
           "--gro", os.path.join(d, "system.gro"),
           "--top", os.path.join(d, "system.top"),
           "--itp-dir", d, "--lipids", LIPIDS,
           "--out-json", os.path.join(d, "r.json"),
           "--out-txt", os.path.join(d, "r.txt")]
    for k, v in kw.items():
        cmd += ["--" + k.replace("_", "-"), str(v)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    with open(os.path.join(d, "r.json")) as fh:
        rep = json.load(fh)
    return p.returncode, rep


def status_of(rep, name):
    for r in rep["checks"]:
        if r["check"] == name:
            return r
    raise AssertionError("no check named %s" % name)


# ----------------------------------------------------------------- tests

def test_a_sound_system_passes_every_check(tmp_path):
    d = str(tmp_path)
    write_system(d)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05)
    assert rep["result"] == "PASS", rep
    assert code == 0


def test_secondary_structure_string_shorter_than_the_protein_fails(tmp_path):
    d = str(tmp_path)
    write_system(d)
    code, rep = run_verify(d, ss_string="C" * 15, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05)
    r = status_of(rep, "secondary_structure")
    assert r["pass"] is False
    assert r["measured"] == 15 and r["expected"] == 20
    assert code != 0


def test_a_build_with_no_recorded_string_fails(tmp_path):
    d = str(tmp_path)
    write_system(d)
    code, rep = run_verify(d, tm_resids="1-20", water_nm=1.0, min_dist=0.05)
    assert status_of(rep, "secondary_structure")["pass"] is False
    assert code != 0


def test_composition_reports_what_was_placed(tmp_path):
    d = str(tmp_path)
    write_system(d, n_upper=40, n_lower=40)
    _, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                        water_nm=1.0, min_dist=0.05)
    m = status_of(rep, "composition")["measured"]
    assert m["upper"]["CHOL"] + m["upper"]["POPC"] == 40
    assert m["lower"]["CHOL"] + m["lower"]["POPC"] == 40


def test_a_ratio_that_was_not_placed_fails(tmp_path):
    d = str(tmp_path)
    write_system(d)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05,
                           expect_upper="CHOL:9 POPC:1")
    assert status_of(rep, "composition")["pass"] is False
    assert code != 0


def test_a_protein_at_the_wrong_depth_fails(tmp_path):
    d = str(tmp_path)
    write_system(d, tm_offset=1.4)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05)
    r = status_of(rep, "protein_placement")
    assert r["pass"] is False
    assert r["measured"] > 1.0
    assert code != 0


def test_a_protein_at_the_midplane_passes(tmp_path):
    d = str(tmp_path)
    write_system(d, tm_offset=0.0)
    _, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                        water_nm=1.0, min_dist=0.05)
    assert status_of(rep, "protein_placement")["pass"] is True


def test_an_unbalanced_charge_fails(tmp_path):
    d = str(tmp_path)
    write_system(d, na=9, cl=6)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05)
    r = status_of(rep, "charge")
    assert r["pass"] is False
    assert abs(r["measured"] - 3.0) < 1e-6
    assert code != 0


def test_a_box_shorter_than_the_protein_fails(tmp_path):
    d = str(tmp_path)
    write_system(d, box=(6.0, 6.0, 3.0), water_span=(2.0, 2.6))
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05, cutoff=1.1)
    r = status_of(rep, "periodic_image")
    assert r["pass"] is False
    assert code != 0


def test_a_thin_water_layer_fails(tmp_path):
    d = str(tmp_path)
    write_system(d)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=5.0, min_dist=0.05)
    assert status_of(rep, "water_layer")["pass"] is False
    assert code != 0


def test_two_beads_on_one_point_fail(tmp_path):
    d = str(tmp_path)
    write_system(d, overlap=True)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.12)
    r = status_of(rep, "overlap")
    assert r["pass"] is False
    assert r["measured"] < 0.12
    assert code != 0


def test_a_named_check_is_accepted_with_allow(tmp_path):
    d = str(tmp_path)
    write_system(d, tm_offset=1.4)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05, allow="protein_placement")
    assert status_of(rep, "protein_placement")["skipped"] is True
    assert rep["result"] == "PASS"
    assert code == 0


def test_a_mismatched_leaflet_measurement_fails(tmp_path):
    d = str(tmp_path)
    write_system(d)
    la = {"result": "FAIL", "shared_fraction": 0.9,
          "shared_species": [{"lipid": "POPC", "upper_nm2": 0.80,
                              "lower_nm2": 0.60,
                              "relative_difference": 0.286,
                              "relative_standard_error": 0.03}]}
    path = os.path.join(d, "leaflet_area.json")
    with open(path, "w") as fh:
        json.dump(la, fh)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05, leaflet_json=path)
    assert status_of(rep, "leaflet_area")["pass"] is False
    assert code != 0


def test_leaflets_that_share_little_are_reported_and_not_failed(tmp_path):
    d = str(tmp_path)
    write_system(d)
    la = {"result": "REPORTED", "shared_fraction": 0.31,
          "shared_species": [{"lipid": "POPC", "upper_nm2": 0.80,
                              "lower_nm2": 0.60,
                              "relative_difference": 0.286,
                              "relative_standard_error": 0.05}]}
    path = os.path.join(d, "leaflet_area.json")
    with open(path, "w") as fh:
        json.dump(la, fh)
    code, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                           water_nm=1.0, min_dist=0.05, leaflet_json=path)
    r = status_of(rep, "leaflet_area")
    assert r["skipped"] is True
    assert "in common" in r["detail"]
    assert code == 0


def test_the_report_carries_the_string_and_the_range(tmp_path):
    d = str(tmp_path)
    write_system(d)
    _, rep = run_verify(d, ss_string="C" * 20, tm_resids="1-20",
                        water_nm=1.0, min_dist=0.05, ss_mode="tm")
    assert rep["build"]["ss_string"] == "C" * 20
    assert rep["build"]["ss_mode"] == "tm"
    assert rep["build"]["n_protein_residues"] == 20


def test_every_failing_check_carries_what_to_do(tmp_path):
    d = str(tmp_path)
    write_system(d, tm_offset=1.4, na=9, cl=6)
    _, rep = run_verify(d, ss_string="C" * 15, tm_resids="1-20",
                        water_nm=5.0, min_dist=0.05)
    failed = [r for r in rep["checks"] if not r["pass"] and not r["skipped"]]
    assert len(failed) >= 3
    for r in failed:
        assert r["remedy"].strip(), r["check"]


# ------------------------------------------------- the equilibration read

def _xvg(path, drift_per_ns):
    import math
    lines = []
    for k in range(0, 1001):
        t = k * 10.0
        g = 1.0 + drift_per_ns * (t / 1000.0)
        w = 0.02 * math.sin(k * 1.7)
        lines.append("%g %g %g" % (t, (8.944 + w) * g, (8.944 - w) * g))
    open(path, "w").write("\n".join(lines) + "\n")


def test_a_settled_membrane_passes(tmp_path):
    p = str(tmp_path / "flat.xvg")
    _xvg(p, 0.0)
    r = subprocess.run([sys.executable, EQUIL, "--xvg", p],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert "settled" in r.stdout


def test_a_drifting_membrane_is_held_back(tmp_path):
    p = str(tmp_path / "drift.xvg")
    _xvg(p, 0.02)
    r = subprocess.run([sys.executable, EQUIL, "--xvg", p],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "still changing" in r.stdout
    assert "What to do" in r.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
