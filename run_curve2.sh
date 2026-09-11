#!/usr/bin/env bash
#
# run_curve2.sh : the same curve as run_curve.sh, along a second path.
#
# run_curve.sh moves DLPC between the leaflets, which changes the number of
# phospholipids each leaflet holds and changes the ratio of sphingomyelin to
# DLPC in the upper leaflet at the same time. This script changes the numbers
# and leaves both ratios alone: the two phospholipids of the upper leaflet are
# scaled together, the three of the lower leaflet are scaled together, and the
# cholesterol of each leaflet is left at the number run_curve.sh starts from.
#
# The two scripts meet at a phospholipid number difference of +8, which is the
# delta = 0 point of run_curve.sh, so that point is not repeated here. Two
# curves that lie on top of each other say the numbers carry the effect. Two
# that separate say the ratios carry part of it.
#
# Usage:
#   GPU=0 ASYM="-32 -16 -4" nohup bash ~/run_curve2.sh > ~/curve2a.log 2>&1 &
#   GPU=1 ASYM="20 32 48"   nohup bash ~/run_curve2.sh > ~/curve2b.log 2>&1 &
#   bash ~/run_curve2.sh --report
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

OUT=${OUT:-$HOME/counts_run}/curve2
ASYM=${ASYM:-"-32 -16 -4 20 32 48"}     # upper minus lower phospholipids
PROD_NS=${PROD_NS:-250}
LAST=${LAST:-0.4}
NT=${NT:-8}                       # OpenMP threads per rank
NTMPI=${NTMPI:-1}                # thread-MPI ranks; raise on a machine
                                 # with many cores and no GPU, for example
                                 # NTMPI=8 NT=8 on 64 cores
MDRUN_EXTRA=${MDRUN_EXTRA:-}     # extra mdrun arguments
GPU=${GPU:-0}
mkdir -p "$OUT"

# the reference build, at a phospholipid number difference of +8
CH_U=${CH_U:-71}; CH_L=${CH_L:-63}
PL_TOT=${PL_TOT:-272}                   # phospholipids of both leaflets together
DL_L0=63; PS_L0=62; P2_L0=7             # the ratio the lower leaflet keeps
APL_U0=${APL_U0:-0.6760}; APL_L0=${APL_L0:-0.7240}
N_U0=211; N_L0=195

report(){
  # Both axes are taken from the build of each point and not from the first
  # frame of its run. curve_report.py reads system.top of each build.
  "$PY" "$R/curve_report.py" "$OUT" --prefix relax_a --last-fraction "$LAST" \
      || echo "  (nothing measured yet)"
}
[ "${1:-}" = "--report" ] && { report; exit 0; }

echo "output     $OUT"
echo "asymmetry  $ASYM   (upper minus lower phospholipids)"
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

for A in $ASYM; do
  echo ""
  echo "=== phospholipid number difference $A ==="
  if [ -f "$OUT/relax_a${A}.json" ]; then echo "  already measured"; continue; fi
  read -r plu pll du sl ok <<< "$("$PY" -c "
A=$A; T=$PL_TOT
plu=(T+A)//2; pll=T-plu
du=plu//2
ok = 1 if (plu%2==0 and plu>=20 and pll>=20) else 0
print(plu, pll, du, pll/132.0, ok)")"
  [ "$ok" = "1" ] || { echo "  skipped: a leaflet would hold too few phospholipids"; continue; }
  nu=$((CH_U + plu)); nl=$((CH_L + pll))
  read -r au al <<< "$("$PY" -c "print('%.5f %.5f'%($APL_U0*$N_U0/$nu, $APL_L0*$N_L0/$nl))")"
  echo "  upper CHOL $CH_U  DLPC $du  PSM $du            ($nu lipids, apl $au)"
  echo "  lower CHOL $CH_L  DLPC:DOPS:POP2 = 63:62:7 scaled to $pll   ($nl lipids, apl $al)"

  export UPPER="CHOL:$("$PY" -c "print('%.5f'%($CH_U/$du))") DLPC:1 PSM:1"
  export LOWER="CHOL:$("$PY" -c "print('%.5f'%($CH_L/($PS_L0*$sl)))") DLPC:$("$PY" -c "print('%.5f'%($DL_L0/$PS_L0))") DOPS:1 POP2_45:$("$PY" -c "print('%.5f'%($P2_L0/$PS_L0))")"
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
  export AUTO_BALANCE=0 MEMBLE_BALANCE_ITER=0 MEMBLE_ALLOW=leaflet_area
  export OUTTAG=a${A}

  D="$OUT/a${A}"; W="$D/a${A}_work"
  mkdir -p "$D"; cd "$D" || exit 1
  if [ ! -f "$W/system.gro" ]; then
    bash "$R/memble.sh" "$PEP" > build.log 2>&1
    [ -f "$W/system.gro" ] || { echo "  the build did not finish; see $D/build.log"; continue; }
  fi
  cd "$W" || exit 1
  awk '/^[A-Z0-9_]+ +[0-9]+$/{print "    "$0}' system.top

  run_stage step6.0_minimization.mdp step6.0 system.gro || continue
  prev=step6.0; ok2=1
  for k in 1 2 3 4 5 6; do
    run_stage step6.${k}_equilibration.mdp step6.${k} "${prev}.gro" step6.0.gro || { ok2=0; break; }
    prev=step6.${k}
  done
  [ $ok2 -eq 1 ] || continue
  echo "  $PROD_NS ns"
  run_stage step7_production.mdp prod "${prev}.gro" || continue
  "$PY" "$R/leaflet_relax.py" --gro system.gro --xtc prod.xtc \
      --last-fraction "$LAST" --json "$OUT/relax_a${A}.json" | sed 's/^/    /'
  cd "$OUT" || exit 1
done

report
echo "print the curve again at any time with:  bash ~/run_curve2.sh --report"
