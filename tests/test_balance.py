"""The balance pass corrects a leaflet mismatch only when the check would fail on it.

The pass lives in memble.sh, in the heredoc that ends with BALEOF. These tests
read that block out of the script itself, so the code under test is the code
that runs, and feed it the leaflet measurements that a build writes.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "memble.sh"


def _extract_block():
    text = SCRIPT.read_text()
    m = re.search(r"<<'BALEOF'\n(.*?)\nBALEOF\n", text, re.S)
    assert m, "the balance block was not found in memble.sh"
    return m.group(1)


@pytest.fixture(scope="module")
def block(tmp_path_factory):
    p = tmp_path_factory.mktemp("balance") / "balance_block.py"
    p.write_text(_extract_block())
    return p


def run(block, tmp_path, species, shared_fraction=1.0, up=0.65, lo=0.65, tol=0.08,
        damp=None):
    """Run the block on one measurement and return worst, error, apls and the flag."""
    (tmp_path / "leaflet_area.json").write_text(json.dumps(
        {"shared_fraction": shared_fraction, "shared_species": species}))
    env = dict(os.environ)
    if damp is not None:
        env["MEMBLE_BALANCE_DAMP"] = str(damp)
    out = subprocess.run([sys.executable, str(block), str(up), str(lo), str(tol)],
                         cwd=tmp_path, env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    f = out.stdout.split()
    assert len(f) == 5, out.stdout
    return dict(worst=float(f[0]), error=float(f[1]), upper=float(f[2]),
                lower=float(f[3]), correct=int(f[4]))


def species(name, upper, lower, error):
    d = abs(upper - lower) / (0.5 * (upper + lower))
    return {"lipid": name, "upper_nm2": upper, "lower_nm2": lower,
            "relative_difference": round(d, 4), "relative_standard_error": error}


def test_two_matched_leaflets_are_not_corrected(block, tmp_path):
    r = run(block, tmp_path, [species("DLPC", 0.700, 0.706, 0.01)])
    assert r["correct"] == 0


def test_a_difference_the_size_of_its_own_error_is_not_corrected(block, tmp_path):
    """A difference of 8.2% carrying a standard error of 5% says nothing."""
    r = run(block, tmp_path, [species("CHOL", 0.594, 0.645, 0.05)])
    assert r["worst"] > 0.08
    assert r["correct"] == 0


def test_a_difference_the_check_would_fail_is_corrected(block, tmp_path):
    r = run(block, tmp_path, [species("CHOL", 0.560, 0.700, 0.01)])
    assert r["correct"] == 1


def test_the_correction_moves_the_two_leaflets_towards_each_other(block, tmp_path):
    r = run(block, tmp_path, [species("CHOL", 0.560, 0.700, 0.01)], up=0.65, lo=0.65)
    # The upper leaflet holds the smaller lipids, so it takes more of them.
    assert r["upper"] > 0.65 > r["lower"]


def test_the_correction_is_damped(block, tmp_path):
    """Half a step, so one pass cannot carry the two leaflets past each other."""
    half = run(block, tmp_path, [species("CHOL", 0.560, 0.700, 0.01)], damp=0.5)
    full = run(block, tmp_path, [species("CHOL", 0.560, 0.700, 0.01)], damp=1.0)
    assert 0.65 < half["upper"] < full["upper"]
    assert 0.65 > half["lower"] > full["lower"]


def test_leaflets_that_share_few_lipids_are_not_corrected(block, tmp_path):
    """A lipid takes the area its neighbours leave it, and the neighbours differ."""
    r = run(block, tmp_path, [species("CHOL", 0.560, 0.700, 0.01)],
            shared_fraction=0.2)
    assert r["correct"] == 0


def test_the_worst_lipid_sets_the_difference(block, tmp_path):
    r = run(block, tmp_path, [species("DLPC", 0.700, 0.706, 0.01),
                              species("CHOL", 0.560, 0.700, 0.01)])
    assert r["worst"] == pytest.approx(0.2222, abs=1e-3)
    assert r["error"] == pytest.approx(0.01)
    assert r["correct"] == 1


def test_no_shared_lipid_writes_nothing(block, tmp_path):
    (tmp_path / "leaflet_area.json").write_text(json.dumps({"shared_species": []}))
    out = subprocess.run([sys.executable, str(block), "0.65", "0.65", "0.08"],
                         cwd=tmp_path, capture_output=True, text=True)
    assert out.stdout.strip() == ""
