#!/usr/bin/env bash
#
# compare_counts.sh : what the lipid numbers of each leaflet cost the membrane.
#
# Three ways of deciding how many lipid molecules each leaflet of an asymmetric
# bilayer receives, built from one composition and run under one protocol:
#
#   eqn       both leaflets get the same area per lipid, so both get the same
#             number of molecules (AUTO_BALANCE=0)
#   table     the area per lipid of each leaflet comes from a table of areas
#             per lipid, which is what a builder does today
#   measured  the area per lipid of each leaflet comes from a measurement on
#             the packed system, and memble builds the system again
#
# Every arm is built once and run from the same protocol with N different sets
# of starting velocities. The membrane then relieves whatever mismatch the
# numbers left, and leaflet_relax.py measures what it did.
#
# Usage:
#   R=/path/to/memble bash compare_counts.sh /path/to/output/dir
#
# The environment memble needs (M3_DIR, GMX, PY, MARTINIZE2, HELPER_*) must be
# set, exactly as for a normal build.
set -uo pipefail

OUT=${1:-$PWD/counts_compare}
R=${R:?set R to the memble directory}
PEP=${PEP:?set PEP to the all-atom protein PDB}
GMX=${GMX:-$(command -v gmx || command -v gmx_mpi)}
PY=${PY:-python3}
NT=${NT:-8}                     # threads for mdrun
SEEDS=${SEEDS:-"1 2 3"}
PROD_NS=${PROD_NS:-100}         # length of each production run, in ns
ARMS=${ARMS:-"eqn table measured"}

# the composition of the comparison: asymmetric, five components, one sterol
export UPPER=${UPPER:-"CHOL:1 DLPC:1 PSM:1"}
export LOWER=${LOWER:-"CHOL:1 DLPC:1 DOPS:1 POP2_45:0.1"}
export BOX_X=${BOX_X:-12} BOX_Y=${BOX_Y:-12} BOX_Z=${BOX_Z:-16}
export TM_RANGE=${TM_RANGE:-65:88} SS_MODE=${SS_MODE:-tm}
export WATER_NM=${WATER_NM:-2.5} SALT_M=${SALT_M:-0.15} N_COPY=${N_COPY:-1}
export NPROD_STEPS=$(( PROD_NS * 50000 ))     # dt = 0.02 ps

mkdir -p "$OUT"
echo "output          $OUT"
echo "arms            $ARMS"
echo "seeds           $SEEDS"
echo "production      $PROD_NS ns per seed ($NPROD_STEPS steps)"
echo ""

run_stage(){   # run_stage <mdp> <deffnm> <start.gro> [restraint.gro]
  local mdp=$1 out=$2 start=$3 restr=${4:-}
  local rflag=(); [ -n "$restr" ] && rflag=(-r "$restr")
  "$GMX" grompp -f "$mdp" -c "$start" "${rflag[@]}" -p system.top -n index.ndx \
      -o "$out.tpr" -maxwarn 10 > "grompp_$out.log" 2>&1 || {
      echo "  grompp failed for $out; see grompp_$out.log"; return 1; }
  "$GMX" mdrun -deffnm "$out" -nt "$NT" > "mdrun_$out.log" 2>&1 || {
      echo "  mdrun failed for $out; see mdrun_$out.log"; return 1; }
  if ls step*[0-9]b.pdb > /dev/null 2>&1; then
      echo "  INSTABILITY during $out: GROMACS wrote step*b.pdb"; return 1; fi
  return 0
}

for arm in $ARMS; do
  echo "=== $arm ==="
  D="$OUT/$arm"
  rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1

  case $arm in
    eqn)      export AUTO_BALANCE=0 MEMBLE_BALANCE_ITER=0 ;;
    table)    export AUTO_BALANCE=1 MEMBLE_BALANCE_ITER=0 ;;
    measured) export AUTO_BALANCE=1 MEMBLE_BALANCE_ITER=2 ;;
  esac
  export OUTTAG=$arm

  bash "$R/memble.sh" "$PEP" > build.log 2>&1
  rc=$?
  W="$D/${arm}_work"
  if [ $rc -ne 0 ] || [ ! -f "$W/system.gro" ]; then
    echo "  the build did not finish (exit $rc); see $D/build.log"
    # a build that memble refuses to hand over is itself a result, and the
    # arm is carried no further
    continue
  fi
  cd "$W" || exit 1
  awk '/^[A-Z0-9_]+ +[0-9]+$/{print "  "$0}' system.top | tail -12

  echo "  minimization and stages 6.1 to 6.5"
  run_stage step6.0_minimization.mdp step6.0 system.gro || continue
  prev=step6.0
  ok=1
  for k in 1 2 3 4 5; do
    run_stage step6.${k}_equilibration.mdp step6.${k} "${prev}.gro" step6.0.gro \
      || { ok=0; break; }
    prev=step6.${k}
  done
  [ $ok -eq 1 ] || continue

  for s in $SEEDS; do
    echo "  seed $s: stage 6.6 and $PROD_NS ns"
    sed '/gen-seed/d' step6.6_equilibration.mdp > s${s}_6.6.mdp
    echo "gen-seed = $s" >> s${s}_6.6.mdp
    run_stage s${s}_6.6.mdp s${s}_step6.6 "${prev}.gro" step6.0.gro || continue
    run_stage step7_production.mdp s${s}_step7 "s${s}_step6.6.gro" || continue
    "$PY" "$R/leaflet_relax.py" --gro system.gro --xtc s${s}_step7.xtc \
        --json "$OUT/relax_${arm}_s${s}.json" | sed 's/^/    /'
  done
  cd "$OUT" || exit 1
done

echo ""
echo "the measurements are in $OUT/relax_<arm>_s<seed>.json"
echo "summarise them with:  python3 $R/summarise_counts.py $OUT"
