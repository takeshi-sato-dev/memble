#!/usr/bin/env bash
#
# run_settled.sh : is the leaflet composition the membrane settled to a fixed point?
#
# The table arm was built with one number of cholesterol in each leaflet and the
# run moved that number. This script builds the same membrane again with the
# number the run settled to, leaving every phospholipid count as it was, and runs
# it. Cholesterol that stays where it was put means the composition reproduces
# itself, and a composition that reproduces itself is the answer the counts were
# asked for.
#
# Usage:  bash ~/run_settled.sh
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

OUT=${OUT:-$HOME/counts_run}
SRC=${SRC:-table}                 # the arm whose settled composition is used
SEED=${SEED:-1}
PROD_NS=${PROD_NS:-200}
NT=${NT:-8}
GPU=${GPU:-1}
D=$OUT/settled; W=$D/settled_work
REL=$OUT/relax_${SRC}_s${SEED}.json
LOG=$OUT/${SRC}/${SRC}_work/coby.log

for f in "$REL" "$LOG" "$PEP"; do
  [ -f "$f" ] || { echo "not found: $f"; exit 1; }
done

# --- what to build ------------------------------------------------------
# The cholesterol of each leaflet is the number the run settled to. Every other
# species keeps the number it was built with, because no molecule of it changed
# leaflet. The area per lipid of each leaflet is scaled by the number of
# molecules that leaflet now holds, so the two leaflets still cover the box.
EV=$("$PY" - "$REL" "$LOG" <<'PYEOF'
import json, re, sys
d = json.load(open(sys.argv[1]))
c = d["counts"]
n = len(d["time_ps"]); cut = n // 2
start = {s: (v["upper"][0], v["lower"][0]) for s, v in c.items()}
tot_chol = sum(start["CHOL"])
up = int(round(sum(c["CHOL"]["upper"][cut:]) / len(c["CHOL"]["upper"][cut:])))
new = {s: list(v) for s, v in start.items()}
new["CHOL"] = [up, tot_chol - up]

txt = open(sys.argv[2]).read()
m = re.search(r'membrane = "(.*?)"', txt, re.S)
apl = [float(x) for x in re.findall(r"apl:([0-9.]+)", m.group(1))]
order = re.findall(r"leaflet:(\w+)(.*?)(?=leaflet:|$)", m.group(1), re.S)
side_species = {k: re.findall(r"lipid:([A-Za-z0-9_]+):", v) for k, v in order}

# The residue field of a GRO file is five characters, so POP2_45 is read back
# as POP2_. The build needs the full name, and the counts are held under the
# key the trajectory reports.
def key(s):
    return s[:5]

out = []
for k, (side, apl0) in enumerate((("upper", apl[0]), ("lower", apl[1]))):
    sp = side_species[side]
    for s in sp:
        if key(s) not in new:
            sys.exit("%s is in the build but not in the trajectory" % s)
    old = sum(start[key(s)][k] for s in sp)
    now = sum(new[key(s)][k] for s in sp)
    ref = float(new[key(sp[1])][k]) if len(sp) > 1 else float(new[key(sp[0])][k])
    spec = " ".join("%s:%.5f" % (s, new[key(s)][k] / ref) for s in sp)
    out.append(('UPPER' if k == 0 else 'LOWER', spec))
    out.append(('APL_UPPER' if k == 0 else 'APL_LOWER', "%.5f" % (apl0 * old / now)))
    print("# %s: %s  (%d -> %d lipids, apl %.4f -> %.4f)"
          % (side, " ".join("%s %d" % (s, new[key(s)][k]) for s in sp),
             old, now, apl0, apl0 * old / now), file=sys.stderr)
for k, v in out:
    print('%s="%s"' % (k, v))
PYEOF
) || { echo "could not read the settled composition"; exit 1; }
echo "$EV" | sed 's/^/  /'
eval "export $(echo "$EV" | tr '\n' ' ')"

export BOX_X=${BOX_X:-12} BOX_Y=${BOX_Y:-12} BOX_Z=${BOX_Z:-16}
export TM_RANGE=${TM_RANGE:-65:88} SS_MODE=${SS_MODE:-tm}
export WATER_NM=${WATER_NM:-2.5} SALT_M=${SALT_M:-0.15} N_COPY=${N_COPY:-1}
export NPROD_STEPS=$(( PROD_NS * 50000 ))
# The area check is measured and reported, and it does not decide this build:
# the composition being tested is the one the membrane chose, and the question
# is what the membrane does with it, not what the area check says about it.
export AUTO_BALANCE=0 MEMBLE_BALANCE_ITER=0 MEMBLE_ALLOW=leaflet_area
export OUTTAG=settled

mkdir -p "$D"; cd "$D" || exit 1
if [ -f "$W/system.gro" ]; then
  echo ">>> the build was already there, and is used as it stands"
else
  echo ">>> building"
  bash "$R/memble.sh" "$PEP" > build.log 2>&1
  [ -f "$W/system.gro" ] || { echo "the build did not finish; see $D/build.log"; exit 1; }
fi
cd "$W" || exit 1
awk '/^[A-Z0-9_]+ +[0-9]+$/{print "  "$0}' system.top
grep -E "LEAFLET AREA MISMATCH|The leaflets are matched" "$D/build.log" | tail -1 | sed 's/^/  check: /'

run_stage(){
  local mdp=$1 out=$2 start=$3 restr=${4:-} rc=0 cont=""
  [ -f "$out.gro" ] && { echo "  $out was already finished"; return 0; }
  if [ ! -f "$out.tpr" ]; then
    if [ -n "$restr" ]; then
      "$GMX" grompp -f "$mdp" -c "$start" -r "$restr" -p system.top -n index.ndx \
          -o "$out.tpr" -maxwarn 10 > "grompp_$out.log" 2>&1 || rc=1
    else
      "$GMX" grompp -f "$mdp" -c "$start" -p system.top -n index.ndx \
          -o "$out.tpr" -maxwarn 10 > "grompp_$out.log" 2>&1 || rc=1
    fi
    [ $rc -eq 0 ] || { echo "  grompp failed for $out"; return 1; }
  fi
  [ -f "$out.cpt" ] && cont="-cpi $out.cpt -append"
  CUDA_VISIBLE_DEVICES=$GPU "$GMX" mdrun -deffnm "$out" -v -ntmpi 1 -ntomp "$NT" \
      $cont >> "mdrun_$out.log" 2>&1 || { echo "  mdrun failed for $out"; return 1; }
  return 0
}

echo ">>> minimization and stages 6.0 to 6.6"
run_stage step6.0_minimization.mdp step6.0 system.gro || exit 1
prev=step6.0
for k in 1 2 3 4 5 6; do
  run_stage step6.${k}_equilibration.mdp step6.${k} "${prev}.gro" step6.0.gro || exit 1
  prev=step6.${k}
done

echo ">>> $PROD_NS ns"
run_stage step7_production.mdp prod "${prev}.gro" || exit 1
"$PY" "$R/leaflet_relax.py" --gro system.gro --xtc prod.xtc \
    --json "$OUT/relax_settled.json" | sed 's/^/  /'
echo ""
echo "the composition reproduces itself when CHOL moved is zero."
