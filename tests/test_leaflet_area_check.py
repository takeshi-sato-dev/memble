"""The leaflet a lipid is given, and the species the two leaflets are compared on.

leaflet_area_check.py measures the area of every lipid in each leaflet of a
freshly built system, and memble stops or corrects the build on the largest
difference it reports between the two leaflets. Both of those depend on two
decisions the script makes before it measures anything: which leaflet each
molecule belongs to, and which species the two leaflets hold in enough numbers
to be compared on.

These tests build systems whose leaflets are known by construction and hand the
script the cases a real membrane produces: leaflets of unequal composition and
unequal bead count, a leaflet that is not flat, a sterol that carries no
phosphate, and a species that one leaflet holds and the other holds a handful of.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "leaflet_area_check.py"

BOX = 12.0
# A phospholipid, five beads, and a longer one, ten. The head bead of each is
# PO4, and the numbers are the distance of each bead from the head.
SHORT = [("PO4", 0.0), ("C1A", 0.5), ("C2A", 1.0), ("C1B", 0.5), ("C2B", 1.0)]
LONG = [("PO4", 0.0), ("C1A", 0.4), ("C2A", 0.8), ("C3A", 1.2), ("C4A", 1.5),
        ("C1B", 0.4), ("C2B", 0.8), ("C3B", 1.2), ("C4B", 1.5), ("C5B", 1.8)]
# A sterol carries no phosphate: its head bead is the hydroxyl.
STEROL = [("ROH", 0.0), ("R1", 0.3), ("R2", 0.6), ("C1", 0.9), ("C2", 1.2)]


def _grid(n, jitter=0.0, seed=0):
    """n positions in the box, on a square grid."""
    side = int(np.ceil(np.sqrt(n)))
    step = BOX / side
    rng = np.random.default_rng(seed)
    pts = []
    for i in range(n):
        x = (i % side + 0.5) * step
        y = (i // side + 0.5) * step
        if jitter:
            x += rng.uniform(-jitter, jitter)
            y += rng.uniform(-jitter, jitter)
        pts.append((x % BOX, y % BOX))
    return pts


def _write_gro(path, mols):
    """mols is a list of (resname, beads, x, y, z_head, sign).

    sign is +1 for a molecule of the upper leaflet, whose beads run downward
    from its head bead, and -1 for one of the lower leaflet.
    """
    lines = []
    n = 0
    for rid, (rn, beads, x, y, zh, sign) in enumerate(mols, start=1):
        for name, dz in beads:
            n += 1
            lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f"
                         % (rid % 100000, rn, name, n % 100000,
                            x, y, zh - sign * dz))
    with open(path, "w") as fh:
        fh.write("test\n%d\n" % n)
        fh.write("\n".join(lines))
        fh.write("\n%10.5f%10.5f%10.5f\n" % (BOX, BOX, 12.0))


def _leaflet(rn, beads, positions, zhead, sign):
    return [(rn, beads, x, y, zhead, sign) for x, y in positions]


def _run(tmp_path, mols, lipids, extra=()):
    gro = tmp_path / "system.gro"
    _write_gro(gro, mols)
    out = tmp_path / "leaflet_area.json"
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--gro", str(gro), "--lipids", lipids,
         "--points", "20000", "--json", str(out)] + list(extra),
        capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return json.loads(out.read_text()), r.stdout


def _naive_misassigned(mols):
    """How many molecules the midplane of the mean lipid bead z would misplace.

    This is the assignment the script used before, and the tests below use it to
    show which of these systems it fails on.
    """
    allz = [zh - sign * dz for rn, beads, x, y, zh, sign in mols
            for _n, dz in beads]
    mid = float(np.mean(allz))
    bad = 0
    for rn, beads, x, y, zh, sign in mols:
        z = float(np.mean([zh - sign * dz for _n, dz in beads]))
        if ("upper" if z >= mid else "lower") != ("upper" if sign > 0 else "lower"):
            bad += 1
    return bad


def _counts(d, leaf):
    return {s: v["n"] for s, v in d["leaflets"][leaf]["species"].items()}


def test_flat_symmetric_bilayer_is_read_correctly(tmp_path):
    """Every molecule of a flat bilayer is given the leaflet it was built into."""
    up = _leaflet("DLPC", SHORT, _grid(200, 0.1, 1), +2.0, +1.0)
    lo = _leaflet("DLPC", SHORT, _grid(200, 0.1, 2), -2.0, -1.0)
    d, _ = _run(tmp_path, up + lo, "DLPC")
    assert _counts(d, "upper") == {"DLPC": 200}
    assert _counts(d, "lower") == {"DLPC": 200}
    assert d["n_placed_by_midplane"] == 0


def test_a_sterol_is_read_from_its_hydroxyl(tmp_path):
    """A sterol carries no phosphate, and its head bead sits below the phosphates."""
    up = (_leaflet("DLPC", SHORT, _grid(150, 0.1, 1), +2.0, +1.0)
          + _leaflet("CHOL", STEROL, _grid(60, 0.1, 3), +1.4, +1.0))
    lo = (_leaflet("DLPC", SHORT, _grid(150, 0.1, 2), -2.0, -1.0)
          + _leaflet("CHOL", STEROL, _grid(60, 0.1, 4), -1.4, -1.0))
    d, _ = _run(tmp_path, up + lo, "DLPC CHOL")
    assert _counts(d, "upper") == {"CHOL": 60, "DLPC": 150}
    assert _counts(d, "lower") == {"CHOL": 60, "DLPC": 150}


def test_leaflets_of_unequal_composition_are_read_correctly(tmp_path):
    """The two leaflets hold different species, and different numbers of beads."""
    up = _leaflet("DLPC", SHORT, _grid(220, 0.1, 1), +2.0, +1.0)
    lo = _leaflet("DOPS", LONG, _grid(180, 0.1, 2), -2.2, -1.0)
    d, _ = _run(tmp_path, up + lo, "DLPC DOPS")
    assert _counts(d, "upper") == {"DLPC": 220}
    assert _counts(d, "lower") == {"DOPS": 180}


def test_a_leaflet_that_is_not_flat_is_read_correctly(tmp_path):
    """A protein depresses the bilayer, and a plane carries lipids across itself.

    A patch of the lower leaflet is raised above the midplane of the mean lipid
    bead z. Those molecules still point downward, so the leaflet they belong to
    is not in doubt, and reading it from the molecule returns it.
    """
    up = _leaflet("DLPC", SHORT, _grid(200, 0.1, 1), +2.0, +1.0)
    pos = _grid(200, 0.1, 2)
    lo = (_leaflet("DLPC", SHORT, pos[:170], -2.0, -1.0)
          + _leaflet("DLPC", SHORT, pos[170:], +1.2, -1.0))
    mols = up + lo
    assert _naive_misassigned(mols) == 30, "the plane should fail on this system"
    d, _ = _run(tmp_path, mols, "DLPC")
    assert _counts(d, "upper") == {"DLPC": 200}
    assert _counts(d, "lower") == {"DLPC": 200}


def test_a_species_one_leaflet_barely_holds_is_not_compared(tmp_path):
    """Four molecules of a species do not make it a species both leaflets hold."""
    up = (_leaflet("DLPC", SHORT, _grid(200, 0.1, 1), +2.0, +1.0)
          + _leaflet("DOPS", LONG, [(0.3, 0.3), (0.5, 0.5), (0.7, 0.3), (0.4, 0.7)],
                     +2.2, +1.0))
    lo = (_leaflet("DLPC", SHORT, _grid(60, 0.1, 2), -2.0, -1.0)
          + _leaflet("DOPS", LONG, _grid(140, 0.1, 4), -2.2, -1.0))
    d, stdout = _run(tmp_path, up + lo, "DLPC DOPS")

    assert _counts(d, "upper")["DOPS"] == 4
    compared = [x["lipid"] for x in d["shared_species"]]
    aside = [x["lipid"] for x in d["not_compared"]]
    assert compared == ["DLPC"]
    assert aside == ["DOPS"]
    assert "not compared" in stdout


def test_a_species_one_leaflet_barely_holds_does_not_drive_the_correction(tmp_path):
    """memble corrects on the largest difference among the compared species.

    The four molecules above are packed into one corner, so the area they are
    given is a fraction of the area of the same species in the other leaflet.
    Left in the comparison, that difference would be the largest one and it
    would set the correction on its own.
    """
    up = (_leaflet("DLPC", SHORT, _grid(200, 0.1, 1), +2.0, +1.0)
          + _leaflet("DOPS", LONG, [(0.3, 0.3), (0.5, 0.5), (0.7, 0.3), (0.4, 0.7)],
                     +2.2, +1.0))
    lo = (_leaflet("DLPC", SHORT, _grid(60, 0.1, 2), -2.0, -1.0)
          + _leaflet("DOPS", LONG, _grid(140, 0.1, 4), -2.2, -1.0))
    d, _ = _run(tmp_path, up + lo, "DLPC DOPS")

    ghost = [x for x in d["not_compared"] if x["lipid"] == "DOPS"][0]
    worst = max(d["shared_species"], key=lambda x: x["relative_difference"])
    assert ghost["relative_difference"] > worst["relative_difference"]
    assert worst["lipid"] == "DLPC"


def test_the_minimum_is_settable(tmp_path):
    """A build that wants those four molecules compared can ask for it."""
    up = (_leaflet("DLPC", SHORT, _grid(200, 0.1, 1), +2.0, +1.0)
          + _leaflet("DOPS", LONG, _grid(4, 0.1, 3), +2.2, +1.0))
    lo = (_leaflet("DLPC", SHORT, _grid(60, 0.1, 2), -2.0, -1.0)
          + _leaflet("DOPS", LONG, _grid(140, 0.1, 4), -2.2, -1.0))
    d, _ = _run(tmp_path, up + lo, "DLPC DOPS",
                extra=("--min-shared", "1", "--min-shared-ratio", "0.0"))
    assert sorted(x["lipid"] for x in d["shared_species"]) == ["DLPC", "DOPS"]
    assert d["not_compared"] == []
