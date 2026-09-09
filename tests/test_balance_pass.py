"""A balance pass keeps the build it started from until the next one is measured.

The pass rebuilds the whole system, and a rebuild can come back worse than the
build that asked for it. These tests run the pass itself, taken out of memble.sh,
around a stand-in for the build step whose measurements are scripted, so that the
sequence of leaflet differences is under control. What is tested is which build
the pass leaves behind.
"""
import json
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "memble.sh"

HARNESS = r"""
_SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
_START_DIR=$(pwd)
_ARG1=${1:-}
PY=python3
OUTTAG=memble
ASYM=1
AREA_TOL=0.08
WORK=$(pwd)/${OUTTAG}_work
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
APL_UP=${APL_UPPER:-0.6500}; APL_LO=${APL_LOWER:-0.6500}
# the stand-in for the build: the measurement of this pass is scripted
_I=${_MEMBLE_ITER:-0}
cp "$MEASUREMENTS/pass${_I}.json" leaflet_area.json
echo "$_I" > pass.txt
echo "BUILD pass=$_I upper=$APL_UP lower=$APL_LO"
__BLOCK__
echo "KEPT pass=$(cat pass.txt) log=${_BAL_LOG:-none}"
"""


def _block():
    text = SCRIPT.read_text()
    m = re.search(r"(# 4c\. BALANCE PASS.*?)\n# =+\n# 5\. per-lipid", text, re.S)
    assert m, "the balance pass was not found in memble.sh"
    return "# " + m.group(1)


def measurement(diff, error=0.01, shared_fraction=1.0):
    """One leaflet measurement: DLPC differing by diff, with that standard error."""
    upper = 0.70
    lower = upper * (2.0 + diff) / (2.0 - diff)   # so that the relative difference is diff
    return {"shared_fraction": shared_fraction,
            "shared_species": [{"lipid": "DLPC", "upper_nm2": round(upper, 4),
                                "lower_nm2": round(lower, 4),
                                "relative_difference": round(diff, 4),
                                "relative_standard_error": error}]}


def run(tmp_path, diffs, iters=2):
    home = tmp_path / "build"
    meas = tmp_path / "meas"
    home.mkdir()
    meas.mkdir()
    for i, d in enumerate(diffs):
        (meas / ("pass%d.json" % i)).write_text(json.dumps(measurement(d)))
    # the last measurement repeats, so a pass that asks for one more build gets one
    for i in range(len(diffs), len(diffs) + 4):
        (meas / ("pass%d.json" % i)).write_text(json.dumps(measurement(diffs[-1])))
    script = home / "fake_memble.sh"
    script.write_text(HARNESS.replace("__BLOCK__", _block()))
    out = subprocess.run(["bash", str(script), "protein.pdb"], cwd=home,
                         capture_output=True, text=True,
                         env={"PATH": "/usr/bin:/bin:/usr/local/bin",
                              "MEASUREMENTS": str(meas),
                              "MEMBLE_BALANCE_ITER": str(iters)})
    assert out.returncode == 0, out.stdout + out.stderr
    return out.stdout, home


def kept_pass(home):
    return int((home / "memble_work" / "pass.txt").read_text().strip())


def test_a_matched_membrane_is_built_once(tmp_path):
    out, home = run(tmp_path, [0.02])
    assert "balance pass" not in out
    assert kept_pass(home) == 0


def test_a_correction_that_works_is_kept(tmp_path):
    out, home = run(tmp_path, [0.20, 0.02])
    assert "balance pass 1" in out
    assert kept_pass(home) == 1


def test_a_correction_that_makes_the_membrane_worse_is_undone(tmp_path):
    out, home = run(tmp_path, [0.20, 0.30])
    assert "balance pass 1" in out
    assert "restores the build" in out
    assert kept_pass(home) == 0


def test_the_pass_stops_at_the_first_build_that_is_not_better(tmp_path):
    """Three builds at most, and the smallest difference of the three survives."""
    out, home = run(tmp_path, [0.30, 0.20, 0.25], iters=5)
    assert kept_pass(home) == 1
    assert "restores the build" in out


def test_the_kept_build_is_not_left_behind(tmp_path):
    out, home = run(tmp_path, [0.20, 0.30])
    assert not (home / "memble_work.balance_keep").exists()
    assert sorted(p.name for p in home.iterdir() if p.is_dir()) == ["memble_work"]


def test_the_build_record_carries_every_difference(tmp_path):
    out, home = run(tmp_path, [0.30, 0.20, 0.25], iters=5)
    log = out.strip().splitlines()[-1]
    assert "pass0:30.0%" in log and "pass1:20.0%" in log and "pass2:25.0%" in log
    assert "restored:pass1" in log


def test_the_iteration_limit_holds(tmp_path):
    out, home = run(tmp_path, [0.30, 0.25, 0.20, 0.15], iters=1)
    assert out.count("balance pass") == 1
    assert kept_pass(home) == 1
