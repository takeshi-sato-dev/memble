#!/usr/bin/env bash
# verify_refactor.sh : does a change to the build change what the build returns?
#
# Builds one small system twice, once with each of two clones, at the same COBY
# seed, and reports every file that differs. A change that only moves code
# leaves this silent. A change that moves a molecule does not.
#
# Usage:
#   OLD=/path/to/clone-before NEW=/path/to/clone-after \
#   M3_DIR=/path/to/lipidome PEP=/path/to/protein.pdb \
#   bash verify_refactor.sh
set -uo pipefail
OLD=${OLD:?set OLD to a clone of the code the systems were built with}
NEW=${NEW:?set NEW to the clone being changed}
SEED=${SEED:-777}
BASE=${BASE:-$(pwd)/verify_refactor}
M3=${M3_DIR:?set M3_DIR to the Martini 3 lipidome}
PEPIN=${PEP:?set PEP to an all-atom protein pdb}
build(){
  local src=$1 out=$2
  rm -rf "$out"; mkdir -p "$out"
  ( cd "$out"
    export R=$src M3_DIR=$M3 PEP=$PEPIN
    export GMX=${GMX:-$(command -v gmx)} PY=${PY:-python3} MARTINIZE2=martinize2 DSSP=mdtraj
    for h in REP:replicate_and_fix_top POS:inject_posres ORI:orient_tm \
             SSDSSP:ss_from_dssp AREA:leaflet_area_check PART:place_partner \
             PRE:prebuild_orient ITP2STRUCT:itp_to_struct DECLASH:declash_gro \
             ADDWATER:add_water PARTPULL:add_partner_pull FIXVS:fix_vsites \
             WHOLE:make_protein_whole ZSHIFT:shift_protein_z FIXRESID:fix_protein_resid \
             MINDIST:check_min_distance GENSS:gen_ss; do
      v=HELPER_${h%%:*}; export "$v=$src/${h##*:}.py"
    done
    export UPPER="CHOL:1 DLPC:1 PSM:1" LOWER="CHOL:1 DLPC:1 PSM:1"
    export BOX_X=10 BOX_Y=10 BOX_Z=16 N_COPY=1 TM_RANGE=65:88 SS_MODE=tm
    export WATER_NM=1.5 SALT_M=0.15 COBY_SEED=$SEED
    bash "$src/memble.sh" "$PEP" > build.log 2>&1 )
  return $?
}
build "$OLD" "$BASE/old"; rc_old=$?
build "$NEW" "$BASE/new"; rc_new=$?
echo "old exit $rc_old   new exit $rc_new"
[ $rc_old -eq 0 ] && [ $rc_new -eq 0 ] || { echo "FAIL: a build did not finish"
  tail -5 "$BASE/new/build.log"; exit 1; }
bad=0
for f in system.gro system.top index.ndx memble_report.txt \
         step6.0_minimization.mdp step6.1_equilibration.mdp \
         step6.2_equilibration.mdp step6.3_equilibration.mdp \
         step6.4_equilibration.mdp step6.5_equilibration.mdp \
         step6.6_equilibration.mdp step7_production.mdp; do
  a=$(md5sum "$BASE/old/memble_work/$f" 2>/dev/null | cut -d' ' -f1)
  b=$(md5sum "$BASE/new/memble_work/$f" 2>/dev/null | cut -d' ' -f1)
  if [ "$a" != "$b" ]; then echo "DIFFER  $f"; bad=1; fi
done
# Every itp the build wrote, restraints included. vermouth writes the citations
# of the force field out of a set, so their order changes between two runs of one
# command and says nothing about the topology. The comment lines are therefore
# compared as a set and every other line byte for byte.
itp_key(){ grep -v '^;' "$1" | md5sum | cut -d' ' -f1; }
itp_cit(){ grep '^;' "$1" | sort | md5sum | cut -d' ' -f1; }
for p in "$BASE/old/memble_work"/*.itp; do
  f=$(basename "$p")
  q="$BASE/new/memble_work/$f"
  [ -f "$q" ] || { echo "MISSING $f"; bad=1; continue; }
  if [ "$(itp_key "$p")" != "$(itp_key "$q")" ]; then echo "DIFFER  $f"; bad=1
  elif [ "$(itp_cit "$p")" != "$(itp_cit "$q")" ]; then echo "DIFFER  $f (comments)"; bad=1; fi
done
[ $bad -eq 0 ] && echo "PASS: the two builds are identical" || echo "FAIL"
exit $bad
