#!/usr/bin/env bash
#
# run_replicate.sh : run one built system again with independent velocities.
#
# A point of the curve carries an uncertainty only if the same system was run
# more than once. Two runs of one composition are separated by their velocities
# and by nothing else: SEED seeds the rotation of a peripheral protein and does
# not reach the packing, so two builds that differ only in SEED are the same
# build and return the same numbers.
#
# This script takes a system that is already built and already minimized, copies
# the state as it stood after step6.1, gives the equilibration a different
# gen-seed, and runs stages 6.2 to 6.6 and the production again in a directory of
# its own. The build is not repeated, so a repeat costs the run and not the
# packing.
#
# Usage:
#   WORK=$HOME/bigbox_d-20/curve/d-20/d-20_work VEL_SEED=2 bash run_replicate.sh
#
# Environment:
#   WORK       the work directory of the built point (holds step6.1.gro)   required
#   VEL_SEED   the gen-seed of this repeat, an integer                     required
#   OUT        where the json is written (default: two levels above WORK)
#   TAG        the name of the point (default: read from the WORK directory name)
#   PROD_NS    production length in ns (default: read from step7_production.mdp)
#   NTMPI, NT, GPU, MDRUN_EXTRA, LAST   as in run_curve.sh
#
set -uo pipefail

R=${R:-$HOME/memble}
VENV=${MEMBLE_VENV:-$HOME/memble-venv}
[ -x "$VENV/bin/python" ] && PY=${PY:-$VENV/bin/python}
PY=${PY:-python3}
GMX=${GMX:-$(command -v gmx)}
NTMPI=${NTMPI:-1}; NT=${NT:-8}; GPU=${GPU:-0}; MDRUN_EXTRA=${MDRUN_EXTRA:-}
LAST=${LAST:-0.4}

WORK=${WORK:?set WORK to the work directory of a point that is already built}
VEL_SEED=${VEL_SEED:?set VEL_SEED to an integer, the gen-seed of this repeat}
[ -d "$WORK" ] || { echo "run_replicate: $WORK is not a directory"; exit 1; }
WORK=$(cd "$WORK" && pwd)
TAG=${TAG:-$(basename "$WORK")}; TAG=${TAG%_work}
OUT=${OUT:-$(cd "$WORK/../.." && pwd)}

if [ ! -f "$WORK/step6.1.gro" ]; then
  echo "run_replicate: $WORK holds no step6.1.gro, so the point is not built and"
  echo "minimized yet. Build it first with run_curve.sh, then repeat it here."
  exit 1
fi

REP="$OUT/${TAG}_v${VEL_SEED}"
mkdir -p "$REP" || exit 1
cd "$REP" || exit 1

# Everything the stages read. system.gro is copied because leaflet_relax.py
# takes the topology of the trajectory from it.
for f in system.gro system.top index.ndx step6.1.gro \
         step6.2_equilibration.mdp step6.3_equilibration.mdp \
         step6.4_equilibration.mdp step6.5_equilibration.mdp \
         step6.6_equilibration.mdp step7_production.mdp; do
  cp -f "$WORK/$f" . || { echo "run_replicate: $WORK holds no $f"; exit 1; }
done
cp -f "$WORK"/*.itp . 2>/dev/null

# The one difference between this run and the run it repeats. A line that is
# already there is replaced rather than added twice, so the file can be read
# and the seed it carries is the seed that was used.
for k in 2 3 4 5 6; do
  m=step6.${k}_equilibration.mdp
  if grep -q "^gen-seed" "$m"; then
    "$PY" - "$m" "$VEL_SEED" <<'PYEOF'
import io, re, sys
p, seed = sys.argv[1], sys.argv[2]
s = io.open(p, encoding="utf-8").read()
io.open(p, "w", encoding="utf-8").write(
    re.sub(r"(?m)^gen-seed\s*=.*$", "gen-seed = %s" % seed, s))
PYEOF
  else
    echo "gen-seed = $VEL_SEED" >> "$m"
  fi
done
echo ">>> $TAG repeat with gen-seed $VEL_SEED, under $REP"

[ -n "${PROD_NS:-}" ] && {
  "$PY" - step7_production.mdp "$PROD_NS" <<'PYEOF'
import io, re, sys
p, ns = sys.argv[1], sys.argv[2]
s = io.open(p, encoding="utf-8").read()
io.open(p, "w", encoding="utf-8").write(
    re.sub(r"(?m)^nsteps\s*=.*$", "nsteps = %d" % (int(ns) * 50000), s))
PYEOF
}

run_stage(){
  local mdp=$1 out=$2 start=$3 restr=${4:-} rc=0 cont=""
  [ -f "$out.gro" ] && return 0
  if [ ! -f "$out.tpr" ]; then
    if [ -n "$restr" ]; then
      "$GMX" grompp -f "$mdp" -c "$start" -r "$restr" -p system.top -n index.ndx \
          -o "$out.tpr" -maxwarn 10 > "grompp_$out.log" 2>&1 || rc=1
    else
      "$GMX" grompp -f "$mdp" -c "$start" -p system.top -n index.ndx \
          -o "$out.tpr" -maxwarn 10 > "grompp_$out.log" 2>&1 || rc=1
    fi
    [ $rc -eq 0 ] || { echo "    grompp failed for $out; see $REP/grompp_$out.log"; return 1; }
  fi
  [ -f "$out.cpt" ] && cont="-cpi $out.cpt -append"
  CUDA_VISIBLE_DEVICES=$GPU "$GMX" mdrun -deffnm "$out" -v -ntmpi "$NTMPI" -ntomp "$NT" \
      $MDRUN_EXTRA $cont >> "mdrun_$out.log" 2>&1 \
      || { echo "    mdrun failed for $out; see $REP/mdrun_$out.log"; return 1; }
  return 0
}

# The restraint reference is the minimized structure of the build that is being
# repeated, exactly as run_curve.sh uses step6.0.gro for every stage.
cp -f "$WORK/step6.0.gro" . 2>/dev/null || cp -f step6.1.gro step6.0.gro
prev=step6.1
for k in 2 3 4 5 6; do
  run_stage step6.${k}_equilibration.mdp step6.${k} "${prev}.gro" step6.0.gro || exit 1
  prev=step6.${k}
done
run_stage step7_production.mdp prod "${prev}.gro" || exit 1

NS=${PROD_NS:-$(awk -F= '/^nsteps/{printf "%d", $2/50000}' step7_production.mdp)}
JSON="$OUT/relax_${TAG}_v${VEL_SEED}_${NS}ns.json"
"$PY" "$R/leaflet_relax.py" --gro system.gro --xtc prod.xtc \
    --last-fraction "$LAST" --json "$JSON" | sed 's/^/    /'
echo ">>> written $JSON"
echo ">>> the name carries the seed and the run length, so curve_report.py"
echo ">>> leaves it out of the curve and reads it as a repeat of $TAG."
