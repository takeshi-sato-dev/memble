#!/usr/bin/env bash
#
# run_curve.sh : where cholesterol settles, as a function of how many
#                phospholipids each leaflet was given.
#
# The number of phospholipids in each leaflet is fixed when the system is built
# and no molecule of them changes leaflet afterwards. Cholesterol does change
# leaflet, and where it settles is not the number it was given. This script asks
# what that settled number depends on.
#
# DLPC is the one phospholipid both leaflets hold, so it is the one moved: the
# upper leaflet is given 70+d of it and the lower 63-d, with every other count
# left alone and with cholesterol started at 71/63 in every case. Each system is
# built, equilibrated and run, and the cholesterol of each leaflet is measured
# over the settled part of the run.
#
# Usage:
#   bash ~/run_curve.sh                       # every point, one GPU
#   GPU=0 DELTAS="-20 -12 -6"  bash ~/run_curve.sh &     # split over two GPUs
#   GPU=1 DELTAS="0 6 12 20"   bash ~/run_curve.sh &
#   bash ~/run_curve.sh --report              # print the curve, run nothing
#
set -uo pipefail
export R=${R:-$HOME/memble}
export PEP=${PEP:-$HOME/chainA.pdb}
export M3_DIR=${M3_DIR:-$HOME/m3lipidome}
export GMX=${GMX:-/usr/local/gromacs-2023.3/bin/gmx}
VENV=${MEMBLE_VENV:-$HOME/memble-venv}
[ -x "$VENV/bin/python" ] && PY=${PY:-$VENV/bin/python}
[ -x "$VENV/bin/martinize2" ] && export MARTINIZE2=${MARTINIZE2:-$VENV/bin/martinize2}
PY=${PY:-python3}; export PY GMX
export DSSP=${DSSP:-mdtraj}
for h in REP:replicate_and_fix_top POS:inject_posres ORI:orient_tm \
         SSDSSP:ss_from_dssp AREA:leaflet_area_check PART:place_partner \
         PRE:prebuild_orient ITP2STRUCT:itp_to_struct DECLASH:declash_gro \
         ADDWATER:add_water PARTPULL:add_partner_pull FIXVS:fix_vsites \
         WHOLE:make_protein_whole ZSHIFT:shift_protein_z FIXRESID:fix_protein_resid \
         MINDIST:check_min_distance GENSS:gen_ss; do
  v=HELPER_${h%%:*}; [ -n "${!v:-}" ] || export "$v=$R/${h##*:}.py"
done

OUT=${OUT:-$HOME/counts_run}/curve
DELTAS=${DELTAS:-"-20 -12 -6 0 6 12 20"}
PROD_NS=${PROD_NS:-250}
LAST=${LAST:-0.4}                # the fraction of the run that is averaged
NT=${NT:-8}                       # OpenMP threads per rank
NTMPI=${NTMPI:-1}                # thread-MPI ranks; raise on a machine
                                 # with many cores and no GPU, for example
                                 # NTMPI=8 NT=8 on 64 cores
MDRUN_EXTRA=${MDRUN_EXTRA:-}     # extra mdrun arguments
GPU=${GPU:-0}
mkdir -p "$OUT"

# the build the curve is measured against: the composition and the areas per
# lipid of the arm built from the table, at delta = 0
CH_U=${CH_U:-71}; DL_U=${DL_U:-70}; SM_U=${SM_U:-70}
CH_L=${CH_L:-63}; DL_L=${DL_L:-63}; PS_L=${PS_L:-62}; P2_L=${P2_L:-7}
APL_U0=${APL_U0:-0.6760}; APL_L0=${APL_L0:-0.7240}
N_U0=$((CH_U + DL_U + SM_U)); N_L0=$((CH_L + DL_L + PS_L + P2_L))

report(){
  # Both axes are taken from the build of each point and not from the first
  # frame of its run. curve_report.py reads system.top of each build.
  "$PY" "$R/curve_report.py" "$OUT" --prefix relax_d --last-fraction "$LAST" \
      || echo "  (nothing measured yet)"
}

[ "${1:-}" = "--report" ] && { report; exit 0; }

echo "output     $OUT"
echo "deltas     $DELTAS"
echo "production $PROD_NS ns per point, GPU $GPU, $NT threads"

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
    [ $rc -eq 0 ] || { echo "    grompp failed for $out"; return 1; }
  fi
  [ -f "$out.cpt" ] && cont="-cpi $out.cpt -append"
  CUDA_VISIBLE_DEVICES=$GPU "$GMX" mdrun -deffnm "$out" -v -ntmpi "$NTMPI" -ntomp "$NT" \
      $MDRUN_EXTRA $cont >> "mdrun_$out.log" 2>&1 || { echo "    mdrun failed for $out"; return 1; }
  return 0
}

for d in $DELTAS; do
  echo ""
  echo "=== delta $d ==="
  if [ -f "$OUT/relax_d${d}.json" ]; then
    echo "  already measured"; continue
  fi
  du=$((DL_U + d)); dl=$((DL_L - d))
  if [ "$du" -lt 5 ] || [ "$dl" -lt 5 ]; then
    echo "  skipped: a leaflet would hold fewer than 5 DLPC"; continue
  fi
  nu=$((CH_U + du + SM_U)); nl=$((CH_L + dl + PS_L + P2_L))
  # the two leaflets still cover the same box, so the area per lipid of each
  # follows the number of molecules that leaflet now holds
  read -r au al <<< "$("$PY" -c "print('%.5f %.5f' % ($APL_U0*$N_U0/$nu, $APL_L0*$N_L0/$nl))")"
  echo "  upper CHOL $CH_U DLPC $du PSM $SM_U   ($nu lipids, apl $au)"
  echo "  lower CHOL $CH_L DLPC $dl DOPS $PS_L POP2_45 $P2_L   ($nl lipids, apl $al)"
  echo "  phospholipids upper $((du + SM_U))  lower $((dl + PS_L + P2_L))  difference $((du + SM_U - dl - PS_L - P2_L))"

  export UPPER="CHOL:$("$PY" -c "print('%.5f'%($CH_U/$SM_U))") DLPC:$("$PY" -c "print('%.5f'%($du/$SM_U))") PSM:1"
  export LOWER="CHOL:$("$PY" -c "print('%.5f'%($CH_L/$PS_L))") DLPC:$("$PY" -c "print('%.5f'%($dl/$PS_L))") DOPS:1 POP2_45:$("$PY" -c "print('%.5f'%($P2_L/$PS_L))")"
  export APL_UPPER="$au" APL_LOWER="$al"
  # box_z is pinned rather than derived. add_water.py would size it from the
  # z span of the protein plus 2*WATER_NM, and every point of the curve must
  # carry the same box so that the only difference between points is the
  # number of phospholipid molecules. 16 nm holds the 10.53 nm protein with
  # 2.74 nm of water each side.
  export BOX_X=${BOX_X:-12} BOX_Y=${BOX_Y:-12} BOX_Z=${BOX_Z:-16}
  export TM_RANGE=${TM_RANGE:-65:88} SS_MODE=${SS_MODE:-tm}
  export WATER_NM=${WATER_NM:-2.5} SALT_M=${SALT_M:-0.15} N_COPY=${N_COPY:-1}
  export NPROD_STEPS=$(( PROD_NS * 50000 ))
  # The area check is measured and written to the report, and it decides
  # nothing here: the point of the curve is what the membrane does with the
  # numbers it was given, not what the area check says about them.
  export AUTO_BALANCE=0 MEMBLE_BALANCE_ITER=0 MEMBLE_ALLOW=leaflet_area
  export OUTTAG=d${d}

  D="$OUT/d${d}"; W="$D/d${d}_work"
  mkdir -p "$D"; cd "$D" || exit 1
  if [ ! -f "$W/system.gro" ]; then
    bash "$R/memble.sh" "$PEP" > build.log 2>&1
    [ -f "$W/system.gro" ] || { echo "  the build did not finish; see $D/build.log"; continue; }
  fi
  cd "$W" || exit 1
  awk '/^[A-Z0-9_]+ +[0-9]+$/{print "    "$0}' system.top

  run_stage step6.0_minimization.mdp step6.0 system.gro || continue
  prev=step6.0; ok=1
  for k in 1 2 3 4 5 6; do
    run_stage step6.${k}_equilibration.mdp step6.${k} "${prev}.gro" step6.0.gro \
      || { ok=0; break; }
    prev=step6.${k}
  done
  [ $ok -eq 1 ] || continue
  echo "  $PROD_NS ns"
  run_stage step7_production.mdp prod "${prev}.gro" || continue
  "$PY" "$R/leaflet_relax.py" --gro system.gro --xtc prod.xtc \
      --last-fraction "$LAST" --json "$OUT/relax_d${d}.json" | sed 's/^/    /'
  cd "$OUT" || exit 1
done

report
echo "print the curve again at any time with:  bash ~/run_curve.sh --report"
