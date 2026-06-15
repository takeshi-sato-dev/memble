#!/usr/bin/env bash
#
# memble.sh   (Martini 3, full CHARMM-GUI-equivalent; GROMACS 2023.x)
#
# AA protein PDB in -> oriented, coarse-grained, membrane-embedded, solvated,
# equilibrated-ready system out. Composition, leaflet asymmetry, box xy, water
# thickness and temperature are all parameters; nothing protein- or
# composition-specific is hardcoded.
#
# COMPOSITION
#   Symmetric:   LIPIDS="CHOL:1 DIPC:1 DPSM:1"        # name:ratio[:head]
#   Asymmetric:  UPPER="CHOL:1 DPSM:2 POPC:1"
#                LOWER="POPC:3 POPE:1"                # both must be set
#   Escape hatch: COBY_MEMBRANE="<verbatim COBY membrane string>"
#
# GEOMETRY / THERMODYNAMICS
#   BOX_X, BOX_Y   nm (default: auto square from N_COPY x SPACING_NM)
#   WATER_NM       water layer per side, nm (default 2.5); BOX_Z is derived as
#                  MEMB_THICK_NM + 2*WATER_NM unless BOX_Z is set directly
#   MEMB_THICK_NM  approx bilayer thickness, nm (default 4.0; composition-dep.)
#   TEMP           kelvin (default 310)
#
# COBY membrane-string grammar (lipid ratio token, leaflet token) is COBY-
# version-dependent and was not verifiable here (README rate-limited). Defaults
# are lipid:NAME:RATIO and leaflet:upper/lower; confirm with `COBY -h`, change
# COBY_LEAFLET if needed, or override the whole string with COBY_MEMBRANE.
# Everything else (itp routing, restraints, assert, box, water, temp) is
# independent of that grammar.
#
# "Any composition" still requires: the lipid exists in your M3 lipidome;
# multi-headgroup lipids may need an explicit head bead in the spec. No global
# elastic network (JM stays flexible).
#
set -eo pipefail

PEP_AA=${1:?usage: memble.sh <all_atom_protein.pdb>}
M3_DIR=${M3_DIR:?set M3_DIR to the Martini 3 lipidome itp directory}
GMX=${GMX:?set GMX to your GROMACS 2023 binary}
HELPER_REP=${HELPER_REP:?set HELPER_REP=/path/replicate_and_fix_top.py}
HELPER_POS=${HELPER_POS:?set HELPER_POS=/path/inject_posres.py}
HELPER_ORI=${HELPER_ORI:?set HELPER_ORI=/path/orient_tm.py}
HELPER_SSDSSP=${HELPER_SSDSSP:-}       # path to ss_from_dssp.py (MULTI_TM orient only)
HELPER_AREA=${HELPER_AREA:?set HELPER_AREA=/path/leaflet_area_check.py}
HELPER_PART=${HELPER_PART:-}   # path to place_partner.py (required only if PARTNER set)
HELPER_ITP2STRUCT=${HELPER_ITP2STRUCT:?set HELPER_ITP2STRUCT=/path/itp_to_struct.py}
HELPER_DECLASH=${HELPER_DECLASH:-$(dirname "$HELPER_ITP2STRUCT")/declash_gro.py}
HELPER_ADDWATER=${HELPER_ADDWATER:-$(dirname "$HELPER_ITP2STRUCT")/add_water.py}
HELPER_PARTPULL=${HELPER_PARTPULL:-$(dirname "$HELPER_ITP2STRUCT")/add_partner_pull.py}
HELPER_FIXVS=${HELPER_FIXVS:-$(dirname "$HELPER_ITP2STRUCT")/fix_vsites.py}
HELPER_FIXRESID=${HELPER_FIXRESID:-$(dirname "$HELPER_ITP2STRUCT")/fix_protein_resid.py}
HELPER_WHOLE=${HELPER_WHOLE:-$(dirname "$HELPER_ITP2STRUCT")/make_protein_whole.py}
HELPER_ZSHIFT=${HELPER_ZSHIFT:-$(dirname "$HELPER_ITP2STRUCT")/shift_protein_z.py}
HELPER_MINDIST=${HELPER_MINDIST:-$(dirname "$HELPER_ITP2STRUCT")/check_min_distance.py}
HELPER_GENSS=${HELPER_GENSS:-$(dirname "$HELPER_ITP2STRUCT")/gen_ss.py}
PY=${PY:-python3}; MARTINIZE2=${MARTINIZE2:-martinize2}

# composition
LIPIDS=${LIPIDS:-"CHOL:1 DIPC:1 DPSM:1"}
UPPER=${UPPER:-}; LOWER=${LOWER:-}
COBY_MEMBRANE=${COBY_MEMBRANE:-}; COBY_LEAFLET=${COBY_LEAFLET:-leaflet}
COBY_APL=${COBY_APL:-0.70}   # nm^2 area per lipid for COBY packing (0.6 default is too dense -> overlaps)
COBY_OPT_STEPS=${COBY_OPT_STEPS:-30}   # COBY overlap-optimizer max steps; kept SMALL because
                                       # declash (fast post-step) resolves overlaps. Large values
                                       # (100s) can hang for hours on dense/large systems.
COBY_PUSH=${COBY_PUSH:-1.0}             # COBY lipid-lipid push multiplier (default 1.0)
# geometry / thermodynamics
N_COPY=${N_COPY:-4}; SPACING_NM=${SPACING_NM:-20}; MARGIN_NM=${MARGIN_NM:-8}
WATER_NM=${WATER_NM:-3.0}; MEMB_THICK_NM=${MEMB_THICK_NM:-4.0}; TEMP=${TEMP:-310}
AREA_TOL=${AREA_TOL:-0.08}; APL_OVERRIDE=${APL_OVERRIDE:-}   # leaflet area pre-check
AUTO_BALANCE=${AUTO_BALANCE:-1}   # auto per-leaflet apl so asymmetric leaflets match in area
AREA_HARD_TOL=${AREA_HARD_TOL:-0.25}   # abort an asymmetric build only above this mismatch
APL_TABLE=${APL_TABLE:-}          # optional "NAME:area ..." overrides for balancing
SALT_M=${SALT_M:-0.15}
# protein / sim
WATER_BIAS=${WATER_BIAS:-0}; SS_OVERRIDE=${SS_OVERRIDE:-}; TM_RANGE=${TM_RANGE:-}
SS_MODE=${SS_MODE:-dssp}   # dssp = let DSSP assign SS (GPCR/multi-helix); tm = TM ranges helix, rest coil (TM-JM peptides); string = use SS_OVERRIDE
OUTTAG=${OUTTAG:-memble}; NPROD_STEPS=${NPROD_STEPS:-400000000}
PARTNER=${PARTNER:-}      # legacy peripheral partner spec; empty = none
SEED=${SEED:-0}           # seed for partner random rotation (use replicate index)
# post-COBY peripheral protein (recommended): give an atomistic PDB + side
PARTNER_PDB=${PARTNER_PDB:-}                 # atomistic peripheral protein PDB
PARTNER_SIDE=${PARTNER_SIDE:-upper}          # which leaflet: upper | lower
PARTNER_GAP=${PARTNER_GAP:-1.5}              # nm head-to-partner-edge initial gap
PARTNER_WATER=${PARTNER_WATER:-3.0}          # nm bulk water beyond the partner
PARTNER_FREEZE_FC=${PARTNER_FREEZE_FC:-1000} # posres fc holding partner in eq.
PARTNER_MARGIN=${PARTNER_MARGIN:-1.0}        # nm outward play before flat-bottom
PARTNER_K=${PARTNER_K:-1000}                 # flat-bottom force constant
PARTNER_ROTATE=${PARTNER_ROTATE:-none}       # none | random
PREBUILT_MULTI=${PREBUILT_MULTI:-0}   # 1 = input PDB already holds the assembled chains
MULTI_TM=${MULTI_TM:-0}               # 1 = multi-pass TM bundle (e.g. 7-TM GPCR): orient on the helix bundle
MULTI_TM_MINLEN=${MULTI_TM_MINLEN:-}  # min helix length (res) counted as a TM helix in MULTI_TM orient
NTERM_SIDE=${NTERM_SIDE:-}            # up|down: force the N-terminus to face +z (up) or -z (down) after orienting
RES_KEEP=${RES_KEEP:-}                 # per-chain residues to keep, e.g. "A:54-103;B:54-103"
TM_CORE=${TM_CORE:-}                   # per-chain TM core to center at z=0, e.g. "A:65-88;..."
HELPER_PRE=${HELPER_PRE:-}             # path to prebuild_orient.py (PREBUILT_MULTI only)
PROT_FC=(1000 500 200 100 50 10); LIP_FC=(400 200 100 50 20 0)
STAGE_DT=(0.002 0.005 0.010 0.015 0.020 0.020); STAGE_NS=(${STAGE_NS:-0.5 0.5 1 1 2 5})

WORK=$(pwd)/${OUTTAG}_work
rm -rf "$WORK"                     # start clean: remove any previous build output
mkdir -p "$WORK"; cd "$WORK"; cp "$PEP_AA" input_aa.pdb

# --- parse composition (symmetric or asymmetric) ---
sed_inplace(){ local e=$1; shift; local f; for f in "$@"; do sed "$e" "$f" > "$f.__si__" && mv "$f.__si__" "$f"; done; }
parse_into(){ local tok nm ra hd
  for tok in $1; do IFS=: read -r nm ra hd <<<"$tok"
    eval "$2+=(\"\$nm\")"; eval "$3+=(\"\${ra:-1}\")"; eval "$4+=(\"\${hd:-auto}\")"
  done; }
declare -a UN UR UH DN DR DH
if [ -n "$UPPER" ] && [ -n "$LOWER" ]; then ASYM=1; parse_into "$UPPER" UN UR UH; parse_into "$LOWER" DN DR DH
else ASYM=0; parse_into "$LIPIDS" UN UR UH; DN=("${UN[@]}"); DR=("${UR[@]}"); DH=("${UH[@]}"); fi
ALL=(); HEADS=()
add_all(){ local NN HH i n h
  eval "NN=(\"\${$1[@]}\")"; eval "HH=(\"\${$2[@]}\")"
  for i in "${!NN[@]}"; do n="${NN[$i]}"; h="${HH[$i]}"
    case " ${ALL[*]} " in *" $n "*) : ;; *) ALL+=("$n"); HEADS+=("$h");; esac
  done; }
add_all UN UH; add_all DN DH
head_of(){ local i; for i in "${!ALL[@]}"; do [ "${ALL[$i]}" = "$1" ] && { echo "${HEADS[$i]}"; return 0; }; done; }
echo ">>> composition ASYM=$ASYM ; lipids to route/restrain: ${ALL[*]}"

# --- auto: dssp + core/solvent/ion itps ---
# DSSP selection. martinize2 needs DSSP 2.2.1/3.0.x; the 4.x CLI is incompatible
# (martinize2 cannot parse it). Prefer a real 3.x binary, then mdtraj's DSSP.
# SS_MODE=tm and SS_OVERRIDE need no DSSP at all.
if [ -z "$DSSP" ]; then
  _DSSP_BIN=$(command -v mkdssp || command -v dssp || true)
  if [ -n "$_DSSP_BIN" ]; then
    _DSSP_VER=$("$_DSSP_BIN" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -1)
    case "$_DSSP_VER" in
      4.*) echo ">>> WARNING: $_DSSP_BIN is DSSP $_DSSP_VER, which martinize2 cannot parse."
           echo ">>>          Install 3.x: micromamba create -n dssp3 -c bioconda 'dssp=3.1.4'"
           echo ">>>          then pass DSSP=/path/to/that/mkdssp, or use SS_MODE=tm. Falling back to mdtraj."
           "$PY" -c 'import mdtraj' 2>/dev/null && DSSP=mdtraj || DSSP="" ;;
      *)   DSSP="$_DSSP_BIN"; echo ">>> using DSSP ${_DSSP_VER:-?} at $_DSSP_BIN" ;;
    esac
  else
    "$PY" -c 'import mdtraj' 2>/dev/null && DSSP=mdtraj || DSSP=""
  fi
fi
[ "$SS_MODE" = tm ] || [ -n "$DSSP" ] || [ -n "$SS_OVERRIDE" ] || { echo "ERROR: no usable DSSP (4.x is not compatible). Install DSSP 3.x and set DSSP=, or use SS_MODE=tm, or set SS_OVERRIDE."; exit 1; }
pick_itp(){ ls $1 2>/dev/null | head -1; }
M3_CORE_ITP=${M3_CORE_ITP:-$(pick_itp "$M3_DIR/martini_v3.0.0.itp")}
M3_SOLV_ITP=${M3_SOLV_ITP:-$(pick_itp "$M3_DIR/*solvent*.itp")}
M3_ION_ITP=${M3_ION_ITP:-$(pick_itp "$M3_DIR/*ion*.itp")}
M3_FFBONDED_ITP=${M3_FFBONDED_ITP:-$(pick_itp "$M3_DIR/*ffbonded*.itp")}   # v2 lipids only; optional
for v in M3_CORE_ITP M3_SOLV_ITP M3_ION_ITP; do [ -n "${!v}" ] || { echo "ERROR: locate $v in $M3_DIR"; exit 1; }; done

# --- route each lipid to the itp that defines it ---
find_itp_for_mol(){ local m=$1 f; for f in "$M3_DIR"/*.itp; do
  awk -v mol="$m" '/^\[ *moleculetype *\]/{g=1;next} g&&NF&&$1!~/^;/{if($1==mol)fd=1;g=0} END{exit !fd}' "$f" && { echo "$f"; return 0; }; done; return 1; }
UNIQ_SRC=()
in_list(){ local x=$1; shift; case " $* " in *" $x "*) return 0;; *) return 1;; esac; }
local_for(){ echo "local_$(basename "$1")"; }
for nm in "${ALL[@]}"; do
  src=$(find_itp_for_mol "$nm") || { echo "ERROR: lipid '$nm' not in any itp under $M3_DIR (not in your M3 lipidome?)"; exit 1; }
  in_list "$src" "${UNIQ_SRC[@]}" || UNIQ_SRC+=("$src")
done

# ====================================================================
# 1. orient TM along z
# ====================================================================
if [ "$PREBUILT_MULTI" = 1 ]; then
  [ -n "$HELPER_PRE" ] || { echo "ERROR: PREBUILT_MULTI=1 needs HELPER_PRE=/path/prebuild_orient.py"; exit 1; }
  [ -n "$RES_KEEP" ]   || { echo "ERROR: PREBUILT_MULTI=1 needs RES_KEEP (e.g. 'A:54-103;B:54-103')"; exit 1; }
  PRE=(--in input_aa.pdb --out oriented_aa.pdb --keep "$RES_KEEP"); [ -n "$TM_CORE" ] && PRE+=(--core "$TM_CORE")
  "$PY" "$HELPER_PRE" "${PRE[@]}"
else
  ORI=(--in input_aa.pdb --out oriented_aa.pdb)
  if [ "$MULTI_TM" = 1 ]; then
    # multi-pass TM (e.g. 7-TM GPCR): orient on the helix bundle, not a single
    # principal axis. Need a per-residue SS string; reuse SS_OVERRIDE if given,
    # else derive helices from DSSP on the input structure.
    if [ -n "$SS_OVERRIDE" ]; then
      ORI_SS="$SS_OVERRIDE"
    elif [ "$DSSP" = mdtraj ] || [ -z "$DSSP" ]; then
      ORI_SS=$("$PY" "$HELPER_SSDSSP" --pdb input_aa.pdb 2>/dev/null) || ORI_SS=""
    else
      ORI_SS=$("$PY" "$HELPER_SSDSSP" --pdb input_aa.pdb --dssp "$DSSP" 2>/dev/null) || ORI_SS=""
    fi
    [ -n "$ORI_SS" ] || { echo "ERROR: MULTI_TM=1 needs a secondary-structure string; set SS_OVERRIDE or provide a working DSSP"; exit 1; }
    ORI+=(--multi-tm --ss "$ORI_SS")
    [ -n "$MULTI_TM_MINLEN" ] && ORI+=(--multi-tm-minlen "$MULTI_TM_MINLEN")
    [ -n "$NTERM_SIDE" ] && ORI+=(--nterm-side "$NTERM_SIDE")
    echo ">>> MULTI_TM=1: orienting on TM helix bundle (SS length ${#ORI_SS})"
  else
    [ -n "$TM_RANGE" ] && ORI+=(--tm-range "$TM_RANGE")
    [ -n "$NTERM_SIDE" ] && ORI+=(--nterm-side "$NTERM_SIDE")
  fi
  "$PY" "$HELPER_ORI" "${ORI[@]}"
fi

# ====================================================================
# 2. martinize2 (auto SS; no global EN) + detect protein itp/molname
# ====================================================================
MZ=(-ff martini3001 -f oriented_aa.pdb -x cg_peptide.pdb -o protein_only.top -p backbone -cys auto -maxwarn 10)
# Secondary structure source (SS_OVERRIDE > SS_MODE=tm > SS_MODE=dssp):
if [ -n "$SS_OVERRIDE" ]; then
  MZ+=(-ss "$SS_OVERRIDE")
elif [ "$SS_MODE" = tm ]; then
  SS_TM=${TM_CORE:-}; [ -n "$SS_TM" ] || { [ -n "$TM_RANGE" ] && SS_TM="ALL:$TM_RANGE"; }
  [ -n "$SS_TM" ] || { echo "ERROR: SS_MODE=tm needs TM_CORE (e.g. 'A:65-88;B:65-88') or TM_RANGE"; exit 1; }
  SS_STR=$("$PY" "$HELPER_GENSS" --pdb oriented_aa.pdb --tm "$SS_TM") || { echo "ERROR: SS generation failed"; exit 1; }
  [ -n "$SS_STR" ] || { echo "ERROR: SS_MODE=tm produced an empty SS string"; exit 1; }
  echo ">>> SS_MODE=tm: TM=$SS_TM -> SS length ${#SS_STR} (TM helix, rest coil)"
  MZ+=(-ss "$SS_STR")
elif [ "$DSSP" = mdtraj ]; then
  MZ+=(-dssp)
else
  MZ+=(-dssp "$DSSP")
fi
[ "$WATER_BIAS" -eq 1 ] && MZ+=(-water-bias -water-bias-eps E:-0.5 C:1.0 H:-1.0)
"$MARTINIZE2" "${MZ[@]}"
res_count_itp(){ awk '/^\[/{a=($2=="atoms")?1:0;next} a&&NF>0&&$1!~/^;/{print $3}' "$1" | sort -un | wc -l | tr -d ' '; }
if [ "$PREBUILT_MULTI" = 1 ]; then
  PROT_ITPS=(molecule_*.itp)
  PROT_MOLMAP=$(awk '/\[ *molecules *\]/{m=1;next} m&&/^\[/{m=0} m&&NF&&$1!~/^;/{printf (n++?":":"") $1} END{printf "\n"}' protein_only.top)
  [ -n "$PROT_MOLMAP" ] || { echo "ERROR: could not read [molecules] from protein_only.top"; exit 1; }
  PROT_BLOCKS=""; for mt in $(echo "$PROT_MOLMAP" | tr ':' ' '); do PROT_BLOCKS="$PROT_BLOCKS $(res_count_itp "${mt}.itp")"; done
  NCHAINS=$(echo "$PROT_MOLMAP" | tr ':' '\n' | grep -c .)
  echo ">>> PREBUILT_MULTI: chains=$NCHAINS molmap=$PROT_MOLMAP blocks=[$PROT_BLOCKS ] itps=${PROT_ITPS[*]}"
else
  PROT_ITP=$(ls molecule_*.itp 2>/dev/null | head -1); [ -n "$PROT_ITP" ] || PROT_ITP=$(grep -l moleculetype ./*.itp | head -1)
  PROT_NAME=$(awk '/^\[ *moleculetype *\]/{f=1;next} f&&NF&&$1!~/^;/{print $1; exit}' "$PROT_ITP")
  PROT_ITPS=("$PROT_ITP")
  N_TMJM=$(res_count_itp "$PROT_ITP")
  PROT_MOLMAP=""; for ((i=0;i<N_COPY;i++)); do PROT_MOLMAP="$PROT_MOLMAP:$PROT_NAME"; done; PROT_MOLMAP=${PROT_MOLMAP#:}
  PROT_BLOCKS=""; for ((i=0;i<N_COPY;i++)); do PROT_BLOCKS="$PROT_BLOCKS $N_TMJM"; done
  NCHAINS=$N_COPY
  echo ">>> protein itp=$PROT_ITP moleculetype=$PROT_NAME ; N_TMJM=$N_TMJM"
fi

# ====================================================================
# 3. box (xy explicit or auto) + water-derived z + grid placement
# ====================================================================
# box_z: membrane thickness + water per side. NOTE: COBY's lipid-grid optimiser
# can fail to converge (hang) when box_z is much larger than the membrane, so we
# size box_z from the membrane, not from a long protruding protein. A protein
# whose hydrophilic ends stick out slightly is tolerated; for a much taller water
# box, build at this size then expand+re-solvate downstream.
ZSPAN=$("$PY" - cg_peptide.pdb <<'PZ'
import sys
zs=[float(l[46:54]) for l in open("cg_peptide.pdb") if l.startswith(("ATOM","HETATM"))]
print("%.3f" % ((max(zs)-min(zs))/10.0) if zs else "0.0")
PZ
)
if [ "$PREBUILT_MULTI" = 1 ]; then
  ASSEMBLY=cg_peptide.pdb
  EX=$("$PY" - cg_peptide.pdb <<'PZ'
import sys
xs=[];ys=[]
for l in open(sys.argv[1]):
    if l.startswith(("ATOM","HETATM")): xs.append(float(l[30:38])); ys.append(float(l[38:46]))
print("%.3f %.3f"%((max(xs)-min(xs))/10.0,(max(ys)-min(ys))/10.0))
PZ
)
  EXX=${EX%% *}; EXY=${EX##* }
  BOX_X=${BOX_X:-$(awk -v e="$EXX" -v m="$MARGIN_NM" 'BEGIN{print e+2*m}')}
  BOX_Y=${BOX_Y:-$(awk -v e="$EXY" -v m="$MARGIN_NM" 'BEGIN{print e+2*m}')}
  BOX_Z=${BOX_Z:-$(awk -v m="$MEMB_THICK_NM" -v w="$WATER_NM" 'BEGIN{print m+2*w}')}
  echo ">>> box=${BOX_X}x${BOX_Y}x${BOX_Z} nm (prebuilt assembly; protein z-span ${ZSPAN}, water ${WATER_NM}/side); T=${TEMP}K"
else
  NCOLS=$(awk -v n="$N_COPY" 'BEGIN{c=sqrt(n);ci=int(c);if(ci<c)ci++;print ci}')
  BOX_X=${BOX_X:-$(awk -v c="$NCOLS" -v s="$SPACING_NM" 'BEGIN{print c*s}')}
  BOX_Y=${BOX_Y:-$BOX_X}
  BOX_Z=${BOX_Z:-$(awk -v m="$MEMB_THICK_NM" -v w="$WATER_NM" 'BEGIN{print m+2*w}')}
  echo ">>> box=${BOX_X}x${BOX_Y}x${BOX_Z} nm (protein z-span ${ZSPAN}, water ${WATER_NM}/side); N_COPY=$N_COPY; T=${TEMP}K"
  "$PY" "$HELPER_REP" replicate --in cg_peptide.pdb --out cg_xN.pdb --grid "$N_COPY" --box "${BOX_X},${BOX_Y}" --margin "$MARGIN_NM"
  ASSEMBLY=cg_xN.pdb
fi
declare -a PART_NAMES=() PART_ITPS=() PART_COUNTS=()
# Peripheral protein is now placed AFTER the membrane is built (post-COBY), so it
# never inflates the box COBY sees. Here we only coarse-grain it and remember the
# files; placement + box growth + restraints happen further down.
PARTNER_CG=""; PARTNER_NAME=""; PARTNER_ITP=""
if [ -n "$PARTNER_PDB" ]; then
  [ -n "$HELPER_PART" ] || { echo "ERROR: PARTNER_PDB set but HELPER_PART not given"; exit 1; }
  pd=partner_cg; mkdir -p "$pd"
  "$MARTINIZE2" -ff martini3001 -f "$PARTNER_PDB" -x "$pd/cg.pdb" -o "$pd/top.top" \
      -elastic -p backbone -cys auto -maxwarn 10
  pit=$(ls "$pd"/molecule_*.itp 2>/dev/null | head -1); [ -n "$pit" ] || pit=$(grep -l moleculetype "$pd"/*.itp | head -1)
  PARTNER_NAME="PARTNER0"
  sed "s/\bmolecule_0\b/$PARTNER_NAME/g" "$pit" > "${PARTNER_NAME}.itp"
  PARTNER_CG="$pd/cg.pdb"; PARTNER_ITP="${PARTNER_NAME}.itp"
  echo ">>> peripheral protein coarse-grained: $PARTNER_NAME (placed post-COBY on $PARTNER_SIDE leaflet)"
fi

# ====================================================================
# 4. COBY build (membrane string from composition / leaflets)
# ====================================================================
PACK="optimize_run:yes optimize_max_steps:${COBY_OPT_STEPS} optimize_lipid_push_multiplier:${COBY_PUSH}"
# Auto-balance the two leaflets of an asymmetric membrane: lipids occupy
# different areas, so one apl for both leaflets leaves them area-mismatched and
# the bilayer under stress. Size a per-leaflet apl from each composition (the
# average stays at COBY_APL to preserve packing density), unless the user pinned
# the apl with COBY_APL explicitly or gave a full COBY_MEMBRANE string.
APL_UP="$COBY_APL"; APL_LO="$COBY_APL"
_HELPDIR_B=$(dirname "$HELPER_ITP2STRUCT")
if [ "$ASYM" -eq 1 ] && [ -z "$COBY_MEMBRANE" ] && [ "$AUTO_BALANCE" != 0 ] \
   && [ -f "$_HELPDIR_B/balance_apl.py" ]; then
  BAL=$("$PY" "$_HELPDIR_B/balance_apl.py" --upper "$UPPER" --lower "$LOWER" \
        --base-apl "$COBY_APL" ${APL_TABLE:+--apl "$APL_TABLE"} 2>/dev/null) || BAL=""
  if [ -n "$BAL" ]; then
    APL_UP=${BAL%% *}; APL_LO=${BAL##* }
    echo ">>> leaflet auto-balance: upper apl=${APL_UP} lower apl=${APL_LO} (from composition; set AUTO_BALANCE=0 to disable)"
  fi
fi
if [ -n "$COBY_MEMBRANE" ]; then MEMB="$COBY_MEMBRANE"
elif [ "$ASYM" -eq 1 ]; then
  MEMB="$PACK ${COBY_LEAFLET}:upper apl:${APL_UP}"; for i in "${!UN[@]}"; do MEMB+=" lipid:${UN[$i]}:${UR[$i]}:params:m3lib"; done
  MEMB+=" ${COBY_LEAFLET}:lower apl:${APL_LO}"; for i in "${!DN[@]}"; do MEMB+=" lipid:${DN[$i]}:${DR[$i]}:params:m3lib"; done
else MEMB="$PACK apl:${COBY_APL} "; for i in "${!UN[@]}"; do MEMB+="lipid:${UN[$i]}:${UR[$i]}:params:m3lib "; done; MEMB=${MEMB% }; fi
echo ">>> COBY membrane = $MEMB"
# Center the protein in COBY on its TM-core residues, so the TM core lands at
# the box center where COBY lays the bilayer midplane. The default centering is
# the mean of all beads, which for a transmembrane-juxtamembrane protein is
# pulled off the membrane by the many juxtamembrane beads on one side; the TM
# core then sits away from the bilayer and the long side overhangs the box, and
# COBY wraps the overhanging beads to the other side (which split chains off two
# at a time). The TM-core residue list is the itp-local numbering (martinize
# renumbers each chain from 1), i.e. original minus RES_KEEP start plus 1.
CEN_RES=""
if [ -n "$TM_CORE" ] && [ -n "$RES_KEEP" ]; then
  # Center COBY on the combined centroid of EVERY chain TM core, not just the
  # first chain. COBY indexes residues by position in the assembly, so chain A
  # cores are positions 0..(lenA-1), chain B continues after that, and so on.
  # Centering on one chain core would leave a laterally spread assembly (for
  # example a 2x2 receptor grid) off center, so a smaller box would clip the
  # far chains. Building every chain core range, each shifted by the running
  # residue offset, centers the whole group.
  _offset=0; _ranges=""
  IFS=';' read -ra _KEEPS <<< "$RES_KEEP"
  IFS=';' read -ra _CORES <<< "$TM_CORE"
  for _kspec in "${_KEEPS[@]}"; do
    _ch="${_kspec%%:*}"; _kr="${_kspec#*:}"
    _kl="${_kr%-*}"; _kh="${_kr#*-}"
    [ -n "$_kl" ] && [ -n "$_kh" ] || continue
    _clen=$((_kh - _kl + 1))
    _tspec=""
    for _c in "${_CORES[@]}"; do
      [ "${_c%%:*}" = "$_ch" ] && _tspec="${_c#*:}"
    done
    if [ -n "$_tspec" ]; then
      _tl="${_tspec%-*}"; _th="${_tspec#*-}"
      # COBY indexes residues 0-based by position in the assembly, so the first
      # kept residue is position 0 and the TM core start is (core - keep) with
      # no offset of one. _offset is the running count of residues in earlier
      # chains. These numbers are computed from the TM_CORE and RES_KEEP that
      # the user passes, for example A:65-88 with A:54-103 gives 11-34.
      _a=$((_offset + _tl - _kl)); _b=$((_offset + _th - _kl))
      _ranges="${_ranges:+${_ranges}:}${_a}-${_b}"
    fi
    _offset=$((_offset + _clen))
  done
  CEN_RES="$_ranges"
  [ -n "$CEN_RES" ] && echo ">>> COBY will center the protein on all TM cores (assembly resid ${CEN_RES})"
fi
MOLMAP="$PROT_MOLMAP"
for pn in "${PART_NAMES[@]}"; do MOLMAP="$MOLMAP:$pn"; done
# generate one single-molecule structure per lipid from its itp connectivity and
# register them with COBY via molecule_import (this is what lets a brand-new
# lipidome be used from its itp alone, with no pre-built membrane needed)
: > mol_import.txt
for nm in "${ALL[@]}"; do
  src=$(find_itp_for_mol "$nm")
  HT=$("$PY" "$HELPER_ITP2STRUCT" --itp "$src" --mol "$nm" --out "${nm}.gro" ${M3_FFBONDED_ITP:+--ffbonded "$M3_FFBONDED_ITP"} | sed -n 's/^UPDOWN .* head0:\([0-9][0-9]*\) tail0:\([0-9][0-9]*\)$/\1 \2/p')
  head0=${HT% *}; tail0=${HT#* }
  # COBY 'manual' alignment aligns the head->tail line to the membrane normal,
  # which stands every lipid up regardless of ring shape (the 'principal' method
  # lays flat sterols on their side). Bead indices are 0-based in COBY.
  echo "file:${nm}.gro moleculetype:${nm} params:m3lib alignment:manual upbead:bead:${head0} downbead:bead:${tail0} library_types:lipid" >> mol_import.txt
done
echo ">>> molecule_import:"; cat mol_import.txt
ITP_LIST=("$M3_CORE_ITP" ${M3_FFBONDED_ITP:+"$M3_FFBONDED_ITP"} "${UNIQ_SRC[@]}" "$M3_SOLV_ITP" "$M3_ION_ITP" "${PROT_ITPS[@]}" "${PART_ITPS[@]}")
"$PY" - "$BOX_X" "$BOX_Y" "$BOX_Z" "$SALT_M" "$MOLMAP" "$MEMB" "$ASSEMBLY" "$CEN_RES" "${ITP_LIST[@]}" <<'PYEOF'
import sys, os, COBY
bx, by, bz, salt, molmap, memb, assembly, cen_res = sys.argv[1:9]
itps = sys.argv[9:]
mol_import = [l.strip() for l in open("mol_import.txt")] if os.path.exists("mol_import.txt") else []
# center the protein on its TM-core residues so the core sits at the bilayer
# midplane (box center) and the juxtamembrane side does not overhang the box
prot = "file:%s moleculetypes:%s" % (assembly, molmap)
if cen_res:
    prot += " cen_method:res:%s" % cen_res
COBY.COBY(
    box=[float(bx), float(by), float(bz)], box_type="rectangular",
    membrane=memb,                                       # grammar: confirm via COBY -h
    molecule_import=mol_import,                          # 12-bead lipidome structures from itp
    protein=prot,                                        # mapping + TM-core centering
    solvation="solv:W pos:NA neg:CL salt_molarity:%s" % salt,
    itp_input=["include:%s" % p for p in itps],
    sn="memble", out_sys="system.gro", out_top="system.top", out_log="coby.log",
)
print("COBY build done")
PYEOF
sed_inplace 's/NA+/NA /g; s/CL-/CL /g' system.gro system.top || true

# Optional fine tuning of how deep the protein sits in the membrane. Z_SHIFT
# (nm, default 0) moves only the protein beads in z, found by trying a few
# values and rebuilding. The declash below relaxes the lipids around the
# shifted protein, and the position restraint reference is read from these
# coordinates, so the offset is held through equilibration.
if [ "${Z_SHIFT:-0}" != "0" ] && [ "${Z_SHIFT:-0}" != "0.0" ]; then
  echo ">>> shifting the protein in z by Z_SHIFT=${Z_SHIFT} nm"
  "$PY" "$HELPER_ZSHIFT" --gro system.gro --dz "$Z_SHIFT" || true
fi

# ====================================================================
# 4a. DECLASH: push apart inter-molecular bead overlaps from dense packing
#     (rigid-molecule moves; intramolecular geometry preserved) so that
#     minimisation does not hit infinite Lennard-Jones forces.
# ====================================================================
# COBY places coarse-grained lipids by their real beads but leaves out-of-plane
# virtual sites (sterol ROH/R3, funct 4) flat; GROMACS then reconstructs them
# ~0.1 nm away and they can explode at minimisation. Rebuild every vsite exactly
# from its itp definition before declashing/minimising.
"$PY" "$HELPER_FIXVS" --gro system.gro --top system.top --itp-dir "$M3_DIR" || true
"$PY" "$HELPER_DECLASH" --gro system.gro --lipids "${ALL[*]}" --target 0.21 --iters 200 --exclude-beads "ROH R3" || true

# ====================================================================
# 4b. ADD WATER: COBY builds the membrane in a thin box (a large box_z hangs its
#     lipid-grid optimiser). Now grow box_z so a protruding TM-JM protein gets a
#     real bulk-water cushion, and solvate the new slabs (+ salt). No-op if the
#     protein already fits with the requested water on each side.
# ====================================================================
if [ -n "$PARTNER_CG" ]; then
  # peripheral protein: place on the chosen leaflet, expand box, then fill water
  "$PY" "$HELPER_PART" --system-gro system.gro --system-top system.top \
      --partner-cg "$PARTNER_CG" --partner-name "$PARTNER_NAME" \
      --partner-itp "$PARTNER_ITP" --side "$PARTNER_SIDE" --gap "$PARTNER_GAP" \
      --water-nm "$PARTNER_WATER" --lipids "${ALL[*]}" \
      --rotate "$PARTNER_ROTATE" --seed "$SEED"
  if [ -n "$HELPER_ADDWATER" ]; then
    "$PY" "$HELPER_ADDWATER" --gro system.gro --top system.top \
        --water-nm "$PARTNER_WATER" --salt "$SALT_M" --keep-box || true
  fi
elif [ -n "$HELPER_ADDWATER" ]; then
  "$PY" "$HELPER_ADDWATER" --gro system.gro --top system.top \
      --water-nm "$WATER_NM" --salt "$SALT_M" || true
fi

# final PBC-aware declash: catches any bead sitting just outside the box that
# clashes with the opposite face under periodic boundaries (a common COBY edge
# effect that produces an infinite force on a water at minimisation).
# Make the protein contiguous across PBC and center the system on it, BEFORE the
# final declash. A tall TM-JM protein can straddle the z boundary after
# assembly/solvation; per-atom wrapping upstream then splits a chain across the
# box (consecutive backbone beads a full box apart -> LINCS blowup -> infinite
# force). Recentering can push some lipids/water across the edge, so the final
# declash must run AFTER this, as the last coordinate step, to clean up.
if [ -n "$HELPER_WHOLE" ] && [ -f "$HELPER_WHOLE" ]; then
  "$PY" "$HELPER_WHOLE" --gro system.gro --top system.top --itp-dir . || true
fi

# FINAL declash and the last coordinate-modifying step: guarantees no two beads
# from different molecules overlap. The protein is frozen (never moved), so it
# is not distorted or re-split; only solvent and lipids are pushed apart.
"$PY" "$HELPER_DECLASH" --gro system.gro --lipids "${ALL[*]}" --target 0.21 --iters 200 --exclude-beads "ROH R3" --freeze-protein || true

# Confirm there is no residual overlap that would give an infinite force, so a
# broken system is never shipped (the user sees this before running gmx).
if [ -n "$HELPER_MINDIST" ] && [ -f "$HELPER_MINDIST" ]; then
  "$PY" "$HELPER_MINDIST" --gro system.gro --lipids "${ALL[*]}" --min 0.12 || \
    echo ">>> WARNING: residual overlap detected; see message above."
fi

# Restore the protein's original per-chain residue numbers (martinize renumbers
# every chain from 1, so the assembled chains overlap). The true numbers come
# from oriented_aa.pdb (which preserves the input PDB resSeq); nothing hardcoded.
if [ -f oriented_aa.pdb ] && [ -n "$HELPER_FIXRESID" ]; then
  "$PY" "$HELPER_FIXRESID" --gro system.gro --oriented oriented_aa.pdb \
        --top system.top --itp-dir . || true
fi

# ====================================================================
# 4b. LEAFLET AREA PRE-CHECK (before any gmx MD; abort the build if asymmetric
#     leaflets are area-mismatched, so it is fixed now, not after a melted run)
# ====================================================================
"$PY" "$HELPER_AREA" --gro system.gro --lipids "${ALL[*]}" --asym "$ASYM" --tol "$AREA_TOL" --hard-tol "$AREA_HARD_TOL" ${APL_OVERRIDE:+--apl "$APL_OVERRIDE"}

# ====================================================================
# 5. per-lipid bead-count sanity assert
# ====================================================================
for nm in "${ALL[@]}"; do src=$(find_itp_for_mol "$nm")
  bg=$(awk -v m="$nm" 'BEGIN{mk=substr(m,1,5)} NR>2{rid=substr($0,1,5);rn=substr($0,6,5);gsub(/ /,"",rid);gsub(/ /,"",rn);
    if(rn==mk){if(first==""){first=rid};if(rid==first)c++}}END{print c+0}' system.gro)
  bi=$(awk -v mol="$nm" '
    function secname(l,t){t=l;gsub(/[][ \t\r]/,"",t);return t}
    /^\[/{sec=secname($0); if(sec=="moleculetype"){inmol=0;expectname=1}; next}
    expectname && NF && $1!~/^;/ {inmol=($1==mol);expectname=0;next}
    sec=="atoms" && inmol && NF && $1!~/^;/ {n++}
    END{print n+0}' "$src")
  echo ">>> $nm beads gro=$bg itp=$bi"
  [ "$bg" -eq "$bi" ] && [ "$bg" -ne 0 ] || { echo "ASSERT FAILED for $nm"; exit 1; }
done

# ====================================================================
# 6. staged restraints: protein BB + every lipid head
# ====================================================================
for s in "${UNIQ_SRC[@]}"; do cp "$s" "$(local_for "$s")"; done
stage_args(){ local o="" k=1 v; for v in "$@"; do o="$o --stage POSRES_STEP${k}:${v}"; k=$((k+1)); done; echo "$o"; }
# Map TM_CORE (original residue numbers) to the martinize itp-local numbering
# (each chain renumbered from 1), so equilibration restrains only the membrane-
# spanning core and leaves the juxtamembrane free to move. tm_local_range <idx>
# returns "--resid-min L --resid-max H" for the idx-th protein chain, or empty.
tm_local_range(){
  local idx="$1" ci=0 tok keep_lo tm_lo tm_hi
  [ -z "$TM_CORE" ] && return
  local IFS=';'
  for tok in $RES_KEEP; do
    if [ "$ci" -eq "$idx" ]; then keep_lo=$(echo "${tok#*:}" | cut -d- -f1); fi
    ci=$((ci+1))
  done
  ci=0
  for tok in $TM_CORE; do
    if [ "$ci" -eq "$idx" ]; then tm_lo=$(echo "${tok#*:}" | cut -d- -f1); tm_hi=$(echo "${tok#*:}" | cut -d- -f2); fi
    ci=$((ci+1))
  done
  [ -z "$keep_lo" ] || [ -z "$tm_lo" ] || [ -z "$tm_hi" ] && return
  echo "--resid-min $((tm_lo - keep_lo + 1)) --resid-max $((tm_hi - keep_lo + 1))"
}
# shellcheck disable=SC2046
_pi=0
for itp in "${PROT_ITPS[@]}"; do
  mt=$(awk '/^\[ *moleculetype *\]/{f=1;next} f&&NF&&$1!~/^;/{print $1; exit}' "$itp")
  "$PY" "$HELPER_POS" --itp "$itp" --mol "$mt" --beads BB $(tm_local_range "$_pi") $(stage_args "${PROT_FC[@]}")
  _pi=$((_pi+1))
done
for nm in "${ALL[@]}"; do src=$(find_itp_for_mol "$nm"); loc=$(local_for "$src")
  "$PY" "$HELPER_POS" --itp "$loc" --mol "$nm" --beads "$(head_of "$nm")" $(stage_args "${LIP_FC[@]}")
done
# partner backbone restraints (held during equilibration, free in production)
for pi in "${!PART_ITPS[@]}"; do
  "$PY" "$HELPER_POS" --itp "${PART_ITPS[$pi]}" --mol "${PART_NAMES[$pi]}" --beads BB $(stage_args "${PROT_FC[@]}")
done
# post-COBY peripheral protein: freeze it strongly through ALL equilibration
# stages (so it never touches the relaxing membrane); production has no posres
# (only the flat-bottom pull guard added later).
if [ -n "$PARTNER_ITP" ]; then
  "$PY" "$HELPER_POS" --itp "$PARTNER_ITP" --mol "$PARTNER_NAME" --beads BB \
      $(stage_args "$PARTNER_FREEZE_FC" "$PARTNER_FREEZE_FC" "$PARTNER_FREEZE_FC" \
                   "$PARTNER_FREEZE_FC" "$PARTNER_FREEZE_FC" "$PARTNER_FREEZE_FC")
fi
for s in "${UNIQ_SRC[@]}"; do sed_inplace "s#[^\"]*$(basename "$s")#$(local_for "$s")#g" system.top; done

# ====================================================================
# 7. PSF / CRD / PDB (ParmEd)
# ====================================================================
"$PY" - "${ALL[*]}" "$PROT_BLOCKS" "${PART_COUNTS[*]}" <<'PYEOF'
import sys, os, re, string, parmed as pmd
lipids = set(sys.argv[1].split())
prot_blocks = [int(x) for x in sys.argv[2].split()] if sys.argv[2].strip() else []
pcounts = [int(x) for x in sys.argv[3].split()] if len(sys.argv) > 3 and sys.argv[3].strip() else []
water  = {"W", "WF"}
ions   = {"NA", "CL", "ION", "NA+", "CL-"}
SEG_MEMB = "MEMB"; SEG_SOLV = "SOLV"; SEG_ION = "ION"  # CHARMM-GUI conventions

# ParmEd's GROMACS reader only supports 3-point vsite type 1, but Martini 3
# cholesterol (and some other lipids) use other virtual-site constructions.
# Those sections are irrelevant for PSF/segid generation, so build a flattened,
# vsite-stripped copy of the topology just for ParmEd. The real system.top
# (with vsites intact) is what GROMACS uses; this copy never touches the run.
def _inline(path, seen):
    ap = os.path.abspath(path); base = os.path.dirname(ap); out = []
    if ap in seen: return out
    seen.add(ap)
    with open(ap) as fh:
        for ln in fh:
            m = re.match(r'\s*#include\s+"([^"]+)"', ln)
            if m:
                inc = m.group(1); cand = None
                for d in ([''] if os.path.isabs(inc) else [os.getcwd(), base]):
                    p = inc if os.path.isabs(inc) else os.path.join(d, inc)
                    if os.path.exists(p): cand = p; break
                if cand: out += _inline(cand, seen); continue
            out.append(ln)
    return out
flat = _inline('system.top', set())
keep = []; skip = False
for ln in flat:
    st = ln.strip()
    if st.startswith('['):
        skip = st.strip('[] ').lower().startswith('virtual_sites')
    # Drop the data lines of virtual_sites sections (parmed cannot parse them),
    # but always keep preprocessor directives (#ifdef/#ifndef/#else/#endif/
    # #define) so their pairing stays balanced. A virtual_sites section can sit
    # between a #ifdef and its #endif; dropping the directive lines too would
    # leave an orphan #endif and break the parmed read.
    if skip and not st.startswith('#'):
        continue
    keep.append(ln)
open('system_parmed.top', 'w').write(''.join(keep))
top = pmd.load_file('system_parmed.top', xyz='system.gro', parametrize=False)
# protein residues, in order: ncopy TM-JM copies then each partner; one PRO chain per block
prot = [r for r in top.residues if r.name not in lipids
        and r.name not in water and r.name not in ions]
block_sizes = prot_blocks + pcounts
idx = 0; chain = 0
for b in block_sizes:
    for _ in range(b):
        if idx < len(prot):
            prot[idx].segid = "PRO" + string.ascii_uppercase[min(chain, 25)]; idx += 1
    chain += 1
for r in prot[idx:]:
    r.segid = "PRO" + string.ascii_uppercase[min(chain, 25)]
for r in top.residues:
    if r.name in lipids: r.segid = SEG_MEMB
    elif r.name in water: r.segid = SEG_SOLV
    elif r.name in ions: r.segid = SEG_ION
top.save('system.psf', overwrite=True); top.save('system_vmd.psf', vmd=True, overwrite=True)
top.save('system.crd', format='charmmcrd', overwrite=True); top.save('system_parmed.pdb', overwrite=True)
segs = sorted({r.segid for r in top.residues})
print('ParmEd wrote psf/crd/pdb; segids =', segs)
PYEOF
"$GMX" editconf -f system.gro -o system.pdb >/dev/null 2>&1 || true

# ====================================================================
# 8. mdp (GROMACS 2023.x; TEMP-parameterized) + run.sh
# ====================================================================
RF='nstlist = 20
cutoff-scheme = Verlet
verlet-buffer-tolerance = 0.005
coulombtype = reaction-field
coulomb-modifier = Potential-shift
rcoulomb = 1.1
epsilon_r = 15
epsilon_rf = 0
vdwtype = cutoff
vdw-modifier = Potential-shift-verlet
rvdw = 1.1'
{ echo "integrator = steep"; echo "nsteps = 50000"; echo "emtol = 100.0"; echo "emstep = 0.001"; echo "define = -DFLEXIBLE"; echo "$RF"; } > step6.0_minimization.mdp
# step6.1 is a SECOND minimization with the protein position-restrained (no md
# yet), the CHARMM-GUI Martini Maker approach: it lets the restrained system
# relax from the rough packed start before any dynamics, so the protein cannot
# be kicked out of the membrane by a finite-dt integration of a stiff structure.
{ echo "integrator = steep"; echo "nsteps = 50000"; echo "emtol = 200.0"; echo "emstep = 0.001"
  echo "define = -DPOSRES_STEP1"; echo "refcoord-scaling = all"; echo "$RF"
} > step6.1_equilibration.mdp
# md equilibration stages 6.2..6.6: ramp the timestep, keep the protein
# restrained (weakening per stage via the POSRES_STEP blocks), refcoord-scaling
# = all so the restraint reference follows the box under semiisotropic pressure.
for k in 2 3 4 5 6; do i=$((k-1)); nst=$(awk -v ns="${STAGE_NS[$i]}" -v dt="${STAGE_DT[$i]}" 'BEGIN{printf "%d",(ns/dt)*1000}')
  { echo "integrator = md"; echo "dt = ${STAGE_DT[$i]}"; echo "nsteps = $nst"; echo "$RF"
    echo "tcoupl = v-rescale"; echo "tc-grps = SOLU_MEMB SOLV"; echo "tau-t = 1.0 1.0"; echo "ref-t = $TEMP $TEMP"
    echo "pcoupl = c-rescale"; echo "pcoupltype = semiisotropic"; echo "tau-p = 4.0"
    echo "compressibility = 3e-4 3e-4"; echo "ref-p = 1.0 1.0"; echo "refcoord-scaling = all"
    echo "gen-vel = yes"; echo "gen-temp = $TEMP"; echo "constraints = none"; echo "define = -DPOSRES_STEP$k"
  } > step6.${k}_equilibration.mdp
done
{ echo "integrator = md"; echo "dt = 0.02"; echo "nsteps = $NPROD_STEPS"; echo "$RF"
  echo "tcoupl = v-rescale"; echo "tc-grps = SOLU_MEMB SOLV"; echo "tau-t = 1.0 1.0"; echo "ref-t = $TEMP $TEMP"
  echo "pcoupl = parrinello-rahman"; echo "pcoupltype = semiisotropic"; echo "tau-p = 12.0"
  echo "compressibility = 3e-4 3e-4"; echo "ref-p = 1.0 1.0"
  echo "nstxout-compressed = 5000"; echo "compressed-x-precision = 1000"
  echo "nstlog = 5000"; echo "nstenergy = 5000"
} > step7_production.mdp
LIPRESN="${ALL[*]}"
# index.ndx generated here (not via interactive gmx select in run.sh): groups
# SOLV / MEMB / SOLU / SOLU_MEMB by residue name, matching system.gro atom order.
"$PY" - "$LIPRESN" <<'NDXEOF'
import sys
lip = set(sys.argv[1].split())
solv = {"W", "WF", "NA", "CL", "ION", "NA+", "CL-"}
L = open("system.gro").read().splitlines(); n = int(L[1]); body = L[2:2+n]
g = {"SOLV": [], "MEMB": [], "SOLU": [], "SOLU_MEMB": []}
for i, ln in enumerate(body, 1):
    rn = ln[5:10].strip()
    if rn in solv:
        g["SOLV"].append(i)
    elif rn in lip:
        g["MEMB"].append(i); g["SOLU_MEMB"].append(i)
    else:
        g["SOLU"].append(i); g["SOLU_MEMB"].append(i)
with open("index.ndx", "w") as fh:
    for name in ("SOLV", "MEMB", "SOLU", "SOLU_MEMB"):
        fh.write("[ %s ]\n" % name)
        idx = g[name]
        for k in range(0, len(idx), 15):
            fh.write(" ".join("%d" % x for x in idx[k:k+15]) + "\n")
        fh.write("\n")
print(">>> index.ndx: SOLV=%d MEMB=%d SOLU=%d SOLU_MEMB=%d"
      % (len(g["SOLV"]), len(g["MEMB"]), len(g["SOLU"]), len(g["SOLU_MEMB"])))
NDXEOF
# peripheral protein: add a PARTNER index group and a one-sided flat-bottom pull
# (MEMB vs PARTNER COM along z) to production, so it can associate with its
# leaflet but can never wrap to the other leaflet through PBC.
if [ -n "$PARTNER_NAME" ] && [ -n "$HELPER_PARTPULL" ]; then
  "$PY" "$HELPER_PARTPULL" --gro system.gro --ndx index.ndx \
      --mdp step7_production.mdp --partner-name "$PARTNER_NAME" \
      --lipids "${ALL[*]}" --margin "$PARTNER_MARGIN" --k "$PARTNER_K"
fi
cat > run.sh <<RUNEOF
#!/usr/bin/env bash
set -eo pipefail
GMX="$GMX"
# GROMACS dumps step<N>b.pdb when atoms move too far (system blowing up).
# set -e misses a run that "succeeds" while melting, so check explicitly.
check_blowup(){ if ls step*[0-9]b.pdb >/dev/null 2>&1; then
  echo "INSTABILITY during \$1: GROMACS wrote step*b.pdb (atoms moving too far)."
  echo "  The membrane is blowing up. Inspect: clashes (EM max force), box too small,"
  echo "  or asymmetric leaflet area mismatch. Do not continue."; exit 1; fi; }
\$GMX grompp -f step6.0_minimization.mdp -c system.gro -p system.top -o step6.0.tpr -maxwarn 10
\$GMX mdrun -deffnm step6.0 -v
grep -i "Maximum force" step6.0.log | tail -1 || true   # should be finite, not astronomical
check_blowup step6.0
prev=step6.0
for k in 1 2 3 4 5 6; do
  \$GMX grompp -f step6.\${k}_equilibration.mdp -c \${prev}.gro -r step6.0.gro -p system.top -n index.ndx -o step6.\${k}.tpr -maxwarn 10
  \$GMX mdrun -deffnm step6.\${k} -v
  check_blowup step6.\${k}
  prev=step6.\${k}
done
\$GMX grompp -f step7_production.mdp -c \${prev}.gro -p system.top -n index.ndx -o step7.tpr -maxwarn 10
\$GMX mdrun -deffnm step7 -v
check_blowup step7
echo "DONE: step7.xtc"
RUNEOF
chmod +x run.sh
# also drop the staged, stop-on-failure MD runner (CHARMM-GUI style: one stage at
# a time with a success check) next to the system, copied from the helper folder.
_HELPDIR=$(dirname "$HELPER_ITP2STRUCT")
if [ -f "$_HELPDIR/run_md.sh" ]; then
  cp "$_HELPDIR/run_md.sh" run_md.sh && chmod +x run_md.sh
fi
# Connectivity for viewers, CHARMM-GUI style: write a PSF with real bonds from
# the topology (ParmEd drops Martini bonds), then load system.psf and read
# system.gro or a trajectory on top of it; bonds show on every frame.
if [ -f "$_HELPDIR/write_psf.py" ]; then
  "$PY" "$_HELPDIR/write_psf.py" --gro system.gro --top system.top \
        --out system.psf --itp-dir "$M3_DIR" || true
fi
if [ -f "$_HELPDIR/write_conect_pdb.py" ]; then
  "$PY" "$_HELPDIR/write_conect_pdb.py" --gro system.gro --top system.top \
        --out system_view.pdb --itp-dir "$M3_DIR" || true
fi
cat > view.vmd <<'VMDEOF'
# CHARMM-GUI-style bonded view:
#   vmd -e view.vmd
# system.psf supplies BONDS, system_view.pdb supplies chain IDs (protein = P).
# run_md.sh writes centered, whole *_view.gro/_view.xtc (gmx trjconv -pbc mol
# -center) so the bilayer shows contiguous and mid-box; this loads those when
# present, else the raw files. Select the proteins with:  chain P   (or) protein
mol new system.psf type psf waitfor all
mol addfile system_view.pdb type pdb waitfor all
if { [file exists step7_production_view.xtc] } {
    mol addfile step7_production_view.xtc type xtc waitfor all
} elseif { [file exists step6.6_equilibration_view.gro] } {
    mol addfile step6.6_equilibration_view.gro type gro waitfor all
} elseif { [file exists step7_production.xtc] } {
    mol addfile step7_production.xtc type xtc waitfor all
}
mol delrep 0 top
mol representation VDW 0.6 12.0
mol addrep top
mol representation Bonds 0.3 12.0
mol addrep top
display resetview
VMDEOF
# and a plain-text sheet of the 8 MD stages to copy-paste one at a time
cat > md_steps.txt <<'MDEOF'
==============================================================================
 memble : run equilibration + production manually, one stage at a time
==============================================================================
 How to use:
   - Copy-paste the blocks below one at a time, from top to bottom.
   - After each block, the  ls -l ...gro  line shows whether the .gro was
     produced. If it exists, the stage succeeded -> go to the next block.
   - If the .gro is missing or you see an error, read the end of the mdrun output.
   - mdrun is run single-rank (-ntmpi 1 -ntomp 8) to avoid domain-decomposition
     errors from the protein elastic network. Change the -ntomp number to match
     your core count (find it with:  sysctl -n hw.ncpu  or  nproc).
 Layout: step6.0 = minimization, step6.1..6.6 = equilibration (restraints are
         released in stages), step7 = production.
 Run everything inside the build output directory  memble_work/ .
==============================================================================


# ---- 0) go to the build directory (once) -----------------------------------
cd memble_work


# ---- 1) step6.0 : energy minimization --------------------------------------
gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr \
  -c system.gro -r system.gro -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.0_minimization -v -ntmpi 1 -ntomp 8
ls -l step6.0_minimization.gro
# expect: "Potential Energy" negative, "Maximum force" finite (not inf)


# ---- 2) step6.1 : equilibration 1 ------------------------------------------
gmx grompp -f step6.1_equilibration.mdp -o step6.1_equilibration.tpr \
  -c step6.0_minimization.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.1_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.1_equilibration.gro


# ---- 3) step6.2 : equilibration 2 ------------------------------------------
gmx grompp -f step6.2_equilibration.mdp -o step6.2_equilibration.tpr \
  -c step6.1_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.2_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.2_equilibration.gro


# ---- 4) step6.3 : equilibration 3 ------------------------------------------
gmx grompp -f step6.3_equilibration.mdp -o step6.3_equilibration.tpr \
  -c step6.2_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.3_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.3_equilibration.gro


# ---- 5) step6.4 : equilibration 4 ------------------------------------------
gmx grompp -f step6.4_equilibration.mdp -o step6.4_equilibration.tpr \
  -c step6.3_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.4_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.4_equilibration.gro


# ---- 6) step6.5 : equilibration 5 ------------------------------------------
gmx grompp -f step6.5_equilibration.mdp -o step6.5_equilibration.tpr \
  -c step6.4_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.5_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.5_equilibration.gro


# ---- 7) step6.6 : equilibration 6 (final) ----------------------------------
gmx grompp -f step6.6_equilibration.mdp -o step6.6_equilibration.tpr \
  -c step6.5_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.6_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.6_equilibration.gro
# equilibration is complete once all of these .gro files exist


# ---- 8) step7 : production (long, ~us scale; run in the background) ---------
gmx grompp -f step7_production.mdp -o step7_production.tpr \
  -c step6.6_equilibration.gro -r step6.6_equilibration.gro \
  -p system.top -n index.ndx -maxwarn 20
nohup gmx mdrun -deffnm step7_production -v -ntmpi 1 -ntomp 8 > log_step7.txt 2>&1 &
sleep 20; tail -15 log_step7.txt
# watch progress:  tail -f log_step7.txt
# stop it:         pkill -f step7_production


==============================================================================
 Troubleshooting
==============================================================================
 - "no domain decomposition ..."  -> -ntmpi 1 is not in effect. Re-paste the
   command exactly (it must contain  -ntmpi 1 -ntomp 8 ).
 - mdrun stops with inf           -> problem in system.gro; rebuild from scratch.
 - a few LINCS warnings           -> normal. Many warnings with no .gro produced
   -> read the end of that stage's output.
 - threads: change the -ntomp number to your core count
   (sysctl -n hw.ncpu on macOS, nproc on Linux).
==============================================================================
MDEOF
echo ">>> Build complete in $WORK."
echo ">>> Equilibrate + produce stage-by-stage:  cd $WORK && bash run_md.sh"
echo ">>> Or copy-paste stages manually from:    $WORK/md_steps.txt"
echo ">>> View WITH BONDS (CHARMM-GUI style):     cd $WORK && vmd -e view.vmd"
echo ">>>   (loads system.psf for connectivity, then the trajectory/gro)"
echo ">>> (legacy all-in-one script:             ./run.sh)"
