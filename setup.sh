#!/usr/bin/env bash
#
# setup.sh : create (or recreate) the dedicated virtual environment and verify it.
# Safe to run any time: it deletes and rebuilds the venv from scratch, so a
# botched install is fixed by rerunning this. Your system Python and other
# project environments are never touched.
#
set -euo pipefail
VENV=${MEMBLE_VENV:-$HOME/memble-venv}
HERE=$(cd "$(dirname "$0")" && pwd)

echo ">>> (re)creating venv at $VENV"
rm -rf "$VENV"
python3 -m venv "$VENV"
"$VENV/bin/pip" install --upgrade pip >/dev/null
"$VENV/bin/pip" install -r "$HERE/requirements.txt"

echo ">>> verifying imports"
"$VENV/bin/python" - <<'PYCHECK'
import importlib, sys
mods = ["numpy", "scipy", "matplotlib", "networkx", "alphashape", "shapely",
        "parmed", "COBY", "vermouth", "mdtraj", "streamlit"]
missing = []
for m in mods:
    try:
        importlib.import_module(m)
    except Exception as e:
        missing.append((m, str(e).splitlines()[0] if str(e) else e.__class__.__name__))
if missing:
    print("FAILED imports:")
    for m, why in missing:
        print("  %-12s %s" % (m, why))
    sys.exit(1)
print("all imports OK")
PYCHECK

echo ">>> checking the martinize2 command"
[ -x "$VENV/bin/martinize2" ] && echo "  martinize2: $VENV/bin/martinize2" || echo "  WARNING: martinize2 not found in venv (vermouth should provide it)"

echo ">>> checking GROMACS (not pip-managed)"
command -v "${GMX:-gmx}" >/dev/null && echo "  gmx: $(command -v ${GMX:-gmx})" || echo "  note: gmx not on PATH; install GROMACS separately and pass GMX=gmx"

echo ">>> done. Use ./memble and ./memble-gui (no activate needed)."
echo ">>> reset anytime: rm -rf $VENV  (or just rerun ./setup.sh)"
