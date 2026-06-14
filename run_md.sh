#!/usr/bin/env bash
#
# run_md.sh : run the staged Martini 3 equilibration + production for a system
# built by memble, CHARMM-GUI style (step6.0 min -> step6.1..6.6 equilibration
# -> step7 production), one stage at a time with a clear success check after each.
#
# Run from inside the build's memble_work directory:
#     cd .../memble_work
#     bash run_md.sh                # min + all equilibration, then ask before production
#     bash run_md.sh 6.0            # run a single stage (6.0,6.1,...,6.6,7)
#     bash run_md.sh prod           # run only production (step7)
#     NT=8 GMX=gmx bash run_md.sh   # override threads / gmx binary
#
# DD note: the elastic network on the multi-chain protein forces a large minimum
# domain-decomposition cell, so we run single-rank (-ntmpi 1) with OpenMP threads
# (-ntomp NT); this avoids "no domain decomposition compatible with the box".

set -eo pipefail
GMX=${GMX:-gmx}
NT=${NT:-8}
MAXWARN=${MAXWARN:-20}

# stage -> (mdp, previous coordinate basename)
mdp_for() {
  case "$1" in
    6.0) echo "step6.0_minimization.mdp system" ;;
    6.1) echo "step6.1_equilibration.mdp step6.0_minimization" ;;
    6.2) echo "step6.2_equilibration.mdp step6.1_equilibration" ;;
    6.3) echo "step6.3_equilibration.mdp step6.2_equilibration" ;;
    6.4) echo "step6.4_equilibration.mdp step6.3_equilibration" ;;
    6.5) echo "step6.5_equilibration.mdp step6.4_equilibration" ;;
    6.6) echo "step6.6_equilibration.mdp step6.5_equilibration" ;;
    7)   echo "step7_production.mdp step6.6_equilibration" ;;
    *)   echo ""; ;;
  esac
}

out_for() {
  case "$1" in
    6.0) echo "step6.0_minimization" ;;
    7)   echo "step7_production" ;;
    *)   echo "step${1}_equilibration" ;;
  esac
}

run_stage() {
  local S="$1"
  local info mdp prev out
  info=$(mdp_for "$S"); [ -n "$info" ] || { echo "unknown stage: $S"; exit 1; }
  mdp=${info%% *}; prev=${info##* }; out=$(out_for "$S")

  # restraint reference: the position restraints must pull toward a FIXED target
  # (the minimized structure), not the rolling previous output, or the reference
  # drifts each stage and the protein is gradually pulled out of place. The start
  # coordinates (-c) continue from the previous stage; the restraint reference
  # (-r) is always step6.0_minimization.gro for the equilibration stages.
  local ref="$prev"
  case "$S" in 6.1|6.2|6.3|6.4|6.5|6.6) ref="step6.0_minimization" ;; esac

  [ -f "$mdp" ]        || { echo "!! missing mdp: $mdp"; exit 1; }
  [ -f "${prev}.gro" ] || { echo "!! missing input coords: ${prev}.gro (run the previous stage first)"; exit 1; }

  echo ""
  echo "==================================================================="
  echo "  STAGE $S   ($mdp  <-  ${prev}.gro)"
  echo "==================================================================="
  $GMX grompp -f "$mdp" -o "${out}.tpr" -c "${prev}.gro" -r "${ref}.gro" \
       -p system.top -n index.ndx -maxwarn "$MAXWARN"
  $GMX mdrun -deffnm "$out" -v -ntmpi 1 -ntomp "$NT"

  if [ -f "${out}.gro" ]; then
    local warn=0
    [ -f "${out}.log" ] && warn=$(grep -c -i "LINCS warning\|INSTABILITY" "${out}.log" || true)
    echo ""
    echo ">>>>>> STAGE $S DONE: ${out}.gro written  (LINCS/instability lines: ${warn}) <<<<<<"
    # centered, whole copy for viewing: GROMACS wraps molecules into the box, so
    # a membrane sitting near the z edge looks split. Center on the protein+
    # membrane and make molecules whole (the CHARMM-GUI viewing convention). This
    # does NOT touch ${out}.gro, which carries the velocities for the next stage.
    if [ -f index.ndx ] && grep -q "SOLU_MEMB" index.ndx; then
      printf "SOLU_MEMB\nSystem\n" | $GMX trjconv -s "${out}.tpr" -f "${out}.gro" \
        -o "${out}_view.gro" -pbc mol -center -n index.ndx >/dev/null 2>&1 \
        && echo "        wrote ${out}_view.gro (centered, whole; for VMD)"
      if [ -f "${out}.xtc" ]; then
        printf "SOLU_MEMB\nSystem\n" | $GMX trjconv -s "${out}.tpr" -f "${out}.xtc" \
          -o "${out}_view.xtc" -pbc mol -center -n index.ndx >/dev/null 2>&1 \
          && echo "        wrote ${out}_view.xtc (centered, whole trajectory; for VMD)"
      fi
    fi
  else
    echo "!!!!!! STAGE $S FAILED: ${out}.gro not produced !!!!!!"
    exit 1
  fi
}

# ---- single-stage mode ---------------------------------------------------
if [ -n "$1" ]; then
  if [ "$1" = "prod" ]; then run_stage 7; exit 0; fi
  run_stage "$1"
  exit 0
fi

# ---- full run: min + equilibration -------------------------------------
for S in 6.0 6.1 6.2 6.3 6.4 6.5 6.6; do
  run_stage "$S"
done

echo ""
echo "==================================================================="
echo "  Equilibration complete (step6.0 -> step6.6)."
echo "  Production (step7, long) is NOT started automatically."
echo "  To run it in the background:"
echo "      nohup bash run_md.sh prod > log_step7.txt 2>&1 &"
echo "      tail -f log_step7.txt"
echo "==================================================================="
