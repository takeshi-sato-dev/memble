"""The text of the build, as one string.

memble.sh carries the defaults and sources the stages under lib/ in the order
they run. A test that looks for a line of the build has to look in all of them,
so this returns memble.sh followed by every stage file in the order memble.sh
sources them. A test that reads memble.sh alone finds only the defaults.
"""
import io
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(p):
    return io.open(p, encoding="utf-8").read()


def stage_files(root=None):
    """The stage files, in the order memble.sh sources them."""
    root = root or ROOT
    head = _read(os.path.join(root, "memble.sh"))
    m = re.search(r'MEMBLE_STAGES="\n(.*?)\n"', head, re.S)
    if not m:
        return []
    return [os.path.join(root, "lib", n.strip())
            for n in m.group(1).split("\n") if n.strip()]


def source(root=None):
    root = root or ROOT
    parts = [_read(os.path.join(root, "memble.sh"))]
    parts += [_read(p) for p in stage_files(root)]
    return "\n".join(parts)
