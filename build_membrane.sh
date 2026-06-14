#!/usr/bin/env bash
# build_membrane.sh : build a protein free Martini 3 bilayer for validation.
# Same lipid template, auto balance, water, and repair path as the main build,
# without any protein step, so the area per lipid and thickness of a new
# lipidome can be measured against the values reported for its parameters.
#
# Usage (symmetric single composition):
#   M3_DIR=/path/m3lipidome GMX=gmx PY=python3 \
#   HELPER_ITP2STRUCT=$PWD/itp_to_struct.py \
#   BOX_X=12 BOX_Y=12 LIPIDS="POPC:1" \
#   bash build_membrane.sh
#
# Usage (asymmetric, two leaflets):
#   ... UPPER="CHOL:1 DLPC:1 PSM:1" LOWER="CHOL:1 DLPC:1 DOPS:1" \
#   bash build_membrane.sh
set -euo pipefail

M3_DIR=${M3_DIR:?set M3_DIR=/path/to/martini3 lipid itp dir}
GMX=${GMX:-gmx}; PY=${PY:-python3}
BOX_X=${BOX_X:-12}; BOX_Y=${BOX_Y:-12}
WATER_NM=${WATER_NM:-2.0}; MEMB_THICK_NM=${MEMB_THICK_NM:-4.0}
TEMP=${TEMP:-310}; SALT_M=${SALT_M:-0.15}
COBY_APL=${COBY_APL:-0.70}; COBY_OPT_STEPS=${COBY_OPT_STEPS:-30}; COBY_PUSH=${COBY_PUSH:-1.0}
AUTO_BALANCE=${AUTO_BALANCE:-1}; APL_TABLE=${APL_TABLE:-}
UPPER=${UPPER:-}; LOWER=${LOWER:-}; LIPIDS=${LIPIDS:-}
OUTTAG=${OUTTAG:-membrane}
WORK=${WORK:-$PWD/membrane_work}

HELPER_ITP2STRUCT=${HELPER_ITP2STRUCT:?set HELPER_ITP2STRUCT=/path/itp_to_struct.py}
HD=$(dirname "$HELPER_ITP2STRUCT")
HELPER_ADDWATER=${HELPER_ADDWATER:-$HD/add_water.py}
HELPER_FIXVS=${HELPER_FIXVS:-$HD/fix_vsites.py}
HELPER_DECLASH=${HELPER_DECLASH:-$HD/declash_gro.py}
HELPER_AREA=${HELPER_AREA:-$HD/leaflet_area_check.py}
HELPER_BAL=${HELPER_BAL:-$HD/balance_apl.py}

M3_CORE_ITP=${M3_CORE_ITP:-$M3_DIR/martini_v3.0.0.itp}
M3_SOLV_ITP=${M3_SOLV_ITP:-$M3_DIR/martini_v3.0.0_solvents_v1.itp}
M3_ION_ITP=${M3_ION_ITP:-$M3_DIR/martini_v3.0.0_ions_v1.itp}
M3_FFBONDED_ITP=${M3_FFBONDED_ITP:-}
if [ -z "$M3_FFBONDED_ITP" ]; then
  for _f in "$M3_DIR"/*ffbonded*.itp; do [ -f "$_f" ] && { M3_FFBONDED_ITP="$_f"; break; }; done
fi
[ -n "$M3_FFBONDED_ITP" ] && echo ">>> using bonded parameters from $M3_FFBONDED_ITP"

rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"

# symmetric vs asymmetric
if [ -n "$UPPER" ] && [ -n "$LOWER" ]; then ASYM=1
elif [ -n "$LIPIDS" ]; then ASYM=0; UPPER="$LIPIDS"; LOWER="$LIPIDS"
else echo "ERROR: set LIPIDS=... (symmetric) or UPPER=... LOWER=... (asymmetric)"; exit 1; fi

read -r -a UN <<< "$(echo "$UPPER" | tr ' ' '\n' | cut -d: -f1 | tr '\n' ' ')"
read -r -a UR <<< "$(echo "$UPPER" | tr ' ' '\n' | cut -d: -f2 | tr '\n' ' ')"
read -r -a DN <<< "$(echo "$LOWER" | tr ' ' '\n' | cut -d: -f1 | tr '\n' ' ')"
read -r -a DR <<< "$(echo "$LOWER" | tr ' ' '\n' | cut -d: -f2 | tr '\n' ' ')"
ALL=($(printf "%s\n" "${UN[@]}" "${DN[@]}" | sort -u))

find_itp_for_mol(){ local nm="$1" f; for f in "$M3_DIR"/*.itp; do
  awk -v mol="$nm" '/^\[ *moleculetype *\]/{g=1;next} g&&NF&&$1!~/^;/{if($1==mol)fd=1;g=0} END{exit !fd}' "$f" && { echo "$f"; return 0; }; done; return 0; }
UNIQ_SRC=(); for nm in "${ALL[@]}"; do s=$(find_itp_for_mol "$nm") || true; [ -n "$s" ] || { echo "ERROR: no itp for $nm in $M3_DIR (each lipid must be the first moleculetype in an itp there)"; exit 1; }; UNIQ_SRC+=("$s"); done
UNIQ_SRC=($(printf "%s\n" "${UNIQ_SRC[@]}" | sort -u))

# 1. lipid templates + molecule_import (orientation aware, like the main build)
: > mol_import.txt
for nm in "${ALL[@]}"; do
  src=$(find_itp_for_mol "$nm")
  HT=$("$PY" "$HELPER_ITP2STRUCT" --itp "$src" --mol "$nm" --out "${nm}.gro" ${M3_FFBONDED_ITP:+--ffbonded "$M3_FFBONDED_ITP"} | sed -n 's/^UPDOWN .* head0:\([0-9][0-9]*\) tail0:\([0-9][0-9]*\)$/\1 \2/p')
  head0=${HT% *}; tail0=${HT#* }
  echo "file:${nm}.gro moleculetype:${nm} params:m3lib alignment:manual upbead:bead:${head0} downbead:bead:${tail0} library_types:lipid" >> mol_import.txt
done

# 2. leaflet auto balance (per leaflet apl from composition)
PACK="optimize_run:yes optimize_max_steps:${COBY_OPT_STEPS} optimize_lipid_push_multiplier:${COBY_PUSH}"
APL_UP="$COBY_APL"; APL_LO="$COBY_APL"
if [ "$ASYM" -eq 1 ] && [ "$AUTO_BALANCE" != 0 ] && [ -f "$HELPER_BAL" ]; then
  BAL=$("$PY" "$HELPER_BAL" --upper "$UPPER" --lower "$LOWER" --base-apl "$COBY_APL" ${APL_TABLE:+--apl "$APL_TABLE"} 2>/dev/null) || BAL=""
  [ -n "$BAL" ] && { APL_UP=${BAL%% *}; APL_LO=${BAL##* }; echo ">>> leaflet auto-balance: upper apl=$APL_UP lower apl=$APL_LO"; }
fi
if [ "$ASYM" -eq 1 ]; then
  MEMB="$PACK leaflet:upper apl:${APL_UP}"; for i in "${!UN[@]}"; do MEMB+=" lipid:${UN[$i]}:${UR[$i]}:params:m3lib"; done
  MEMB+=" leaflet:lower apl:${APL_LO}"; for i in "${!DN[@]}"; do MEMB+=" lipid:${DN[$i]}:${DR[$i]}:params:m3lib"; done
else
  MEMB="$PACK apl:${COBY_APL} "; for i in "${!UN[@]}"; do MEMB+="lipid:${UN[$i]}:${UR[$i]}:params:m3lib "; done; MEMB=${MEMB% }
fi
echo ">>> COBY membrane = $MEMB"

BOX_Z=$(awk -v m="$MEMB_THICK_NM" -v w="$WATER_NM" 'BEGIN{print m+2*w}')
ITP_LIST=("$M3_CORE_ITP" ${M3_FFBONDED_ITP:+"$M3_FFBONDED_ITP"} "${UNIQ_SRC[@]}" "$M3_SOLV_ITP" "$M3_ION_ITP")

# 3. COBY build (membrane + solvent, no protein)
"$PY" - "$BOX_X" "$BOX_Y" "$BOX_Z" "$SALT_M" "$MEMB" "${ITP_LIST[@]}" <<'PYEOF'
import sys, os, COBY
bx, by, bz, salt, memb = sys.argv[1:6]
itps = sys.argv[6:]
mol_import = [l.strip() for l in open("mol_import.txt")] if os.path.exists("mol_import.txt") else []
COBY.COBY(
    box=[float(bx), float(by), float(bz)], box_type="rectangular",
    membrane=memb, molecule_import=mol_import,
    solvation="solv:W pos:NA neg:CL salt_molarity:%s" % salt,
    itp_input=["include:%s" % p for p in itps],
    sn="memble", out_sys="system.gro", out_top="system.top", out_log="coby.log",
)
print("COBY membrane build done")
PYEOF
sed -i 's/NA+/NA /g; s/CL-/CL /g' system.gro system.top 2>/dev/null || true

# 4. repair: rebuild vsites, then declash inter-molecular overlaps so the
# equilibration does not hit non-finite forces. The declash arguments match the
# protein build (target 0.21 nm, sterol ROH and R3 vsites excluded), since the
# earlier call passed flags declash does not accept, so it silently did nothing.
"$PY" "$HELPER_FIXVS" --gro system.gro --top system.top --itp-dir "$M3_DIR" || true
"$PY" "$HELPER_DECLASH" --gro system.gro --lipids "${ALL[*]}" --target 0.21 --iters 200 --exclude-beads "ROH R3" || true
"$PY" "$HELPER_AREA" --gro system.gro --lipids "${ALL[*]}" --asym "$ASYM" --tol 0.08 --hard-tol 0.25 || true

# 5. mdp: minimization then NPT equilibration (semiisotropic, no posres)
cat > min.mdp <<MDP
integrator = steep
nsteps = 10000
emtol = 10
emstep = 0.01
define = -DFLEXIBLE
nstlist = 20
cutoff-scheme = Verlet
coulombtype = reaction-field
coulomb-modifier = potential-shift
epsilon_r = 15
epsilon_rf = 0
rcoulomb = 1.1
vdwtype = cutoff
vdw-modifier = potential-shift-verlet
rvdw = 1.1
MDP
cat > eq_soft.mdp <<MDP
integrator = md
dt = 0.005
nsteps = 50000
nstlist = 20
cutoff-scheme = Verlet
coulombtype = reaction-field
coulomb-modifier = potential-shift
epsilon_r = 15
epsilon_rf = 0
rcoulomb = 1.1
vdwtype = cutoff
vdw-modifier = potential-shift-verlet
rvdw = 1.1
tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = ${TEMP}
pcoupl = c-rescale
pcoupltype = semiisotropic
tau-p = 12.0
compressibility = 3e-4 3e-4
ref-p = 1.0 1.0
gen-vel = yes
gen-temp = ${TEMP}
gen-seed = -1
MDP

cat > eq.mdp <<MDP
integrator = md
dt = 0.02
nsteps = ${EQ_STEPS:-250000}
nstlist = 20
cutoff-scheme = Verlet
coulombtype = reaction-field
coulomb-modifier = potential-shift
epsilon_r = 15
epsilon_rf = 0
rcoulomb = 1.1
vdwtype = cutoff
vdw-modifier = potential-shift-verlet
rvdw = 1.1
tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = ${TEMP}
pcoupl = c-rescale
pcoupltype = semiisotropic
tau-p = 12.0
compressibility = 3e-4 3e-4
ref-p = 1.0 1.0
nstxout-compressed = 5000
nstlog = 5000
nstenergy = 5000
gen-vel = yes
gen-temp = ${TEMP}
gen-seed = -1
MDP

cat > run_membrane.sh <<RUN
#!/usr/bin/env bash
set -euo pipefail
${GMX} grompp -f min.mdp -c system.gro -p system.top -o min.tpr -maxwarn 10
${GMX} mdrun -deffnm min -v -ntmpi 1 -ntomp ${NTOMP:-8}
${GMX} grompp -f eq_soft.mdp -c min.gro -p system.top -o eq_soft.tpr -maxwarn 10
${GMX} mdrun -deffnm eq_soft -v -ntmpi 1 -ntomp ${NTOMP:-8}
${GMX} grompp -f eq.mdp -c eq_soft.gro -p system.top -o eq.tpr -maxwarn 10
${GMX} mdrun -deffnm eq -v -ntmpi 1 -ntomp ${NTOMP:-8}
echo "DONE: eq.gro eq.xtc  ->  measure with membrane_props.py"
RUN
chmod +x run_membrane.sh

echo ""
echo ">>> membrane build complete in $WORK"
echo ">>> next: bash run_membrane.sh   then   membrane_props.py --gro eq.gro --traj eq.xtc --top system.gro"
