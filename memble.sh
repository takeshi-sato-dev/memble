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

MEMBLE_VERSION=1.2.2
_SELF=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")
_START_DIR=$(pwd)
_ARG1=${1:-}
case "${1:-}" in
  -v|--version) echo "memble $MEMBLE_VERSION"; exit 0 ;;
esac

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
HELPER_VERIFY=${HELPER_VERIFY:-$(dirname "$HELPER_ITP2STRUCT")/verify_system.py}
_HELPDIR=$(dirname "$HELPER_ITP2STRUCT")
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
# The leaflet area measurement is reported for every build and stops a build only
# where the difference is too large to come from the packing. The mean area per
# lipid of a leaflet is the area of the box, less the area the protein occupies
# in that leaflet, divided by the number of molecules assigned to that leaflet,
# so a difference of a few percent returns the numbers that were assigned rather
# than a fault. A difference of tens of percent does report a fault: a species
# placed in the wrong leaflet, or a packing that failed. AREA_TOL is set to catch
# the second and to leave the first alone.
AREA_TOL=${AREA_TOL:-0.25}; APL_OVERRIDE=${APL_OVERRIDE:-}   # leaflet area pre-check
AUTO_BALANCE=${AUTO_BALANCE:-1}   # auto per-leaflet apl so asymmetric leaflets match in area
AREA_HARD_TOL=${AREA_HARD_TOL:-0.25}   # abort an asymmetric build only above this mismatch
# The pass that rebuilds from the area measurement is off. Correcting the counts
# from that measurement moves them toward equal areas per lipid, which is the
# build the area measurement scores best and the membrane agrees with least.
# MEMBLE_BALANCE_ITER=2 turns the pass back on.
MEMBLE_BALANCE_ITER=${MEMBLE_BALANCE_ITER:-0}
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

# Resolve every user-supplied path before the cd below. The build changes
# directory into $WORK, so a relative path given on the command line would no
# longer resolve: "./memble protein_AA.pdb" used to die at the first cp with
# "No such file or directory". Resolve against the directory the command was
# run from, and fail early with a clear message if the input is missing.
_abs(){ case "$1" in ""|/*) printf '%s' "$1" ;; *) printf '%s/%s' "$(pwd)" "$1" ;; esac; }
PEP_AA=$(_abs "$PEP_AA")
M3_DIR=$(_abs "$M3_DIR")
[ -z "$PARTNER_PDB" ] || PARTNER_PDB=$(_abs "$PARTNER_PDB")
case "$DSSP" in ""|mdtraj|/*) ;; */*) DSSP=$(_abs "$DSSP") ;; esac
[ -f "$PEP_AA" ] || stop "the protein PDB was not found: $PEP_AA" \
  "Give the path to an all-atom PDB as the first argument:\n  memble.sh /path/to/protein.pdb\nA relative path is resolved against the directory you started memble in."
[ -d "$M3_DIR" ] || stop "M3_DIR is not a directory: $M3_DIR" \
  "M3_DIR points at the Martini 3 lipidome, the directory that holds the lipid\nitp files. Set it to the directory you unpacked the lipidome into:\n  export M3_DIR=/path/to/martini3-lipidome\n  ls \$M3_DIR/*.itp | head"
[ -z "$PARTNER_PDB" ] || [ -f "$PARTNER_PDB" ] || stop "PARTNER_PDB was not found: $PARTNER_PDB" \
  "PARTNER_PDB is the peripheral protein that memble places on one leaflet.\nGive the path to its all-atom PDB, or unset PARTNER_PDB to build the membrane\nprotein alone."

WORK=$(pwd)/${OUTTAG}_work
rm -rf "$WORK"                     # start clean: remove any previous build output
mkdir -p "$WORK"; cd "$WORK"; cp "$PEP_AA" input_aa.pdb

# --- parse composition (symmetric or asymmetric) ---
sed_inplace(){ local e=$1; shift; local f; for f in "$@"; do sed "$e" "$f" > "$f.__si__" && mv "$f.__si__" "$f"; done; }

# step <what memble is doing now>
# The terminal shows where a build is in the pipeline, and how long it has taken
# so far. A build that is rerun by the balance pass carries the same clock.
_MEMBLE_T0=${_MEMBLE_T0:-$(date +%s)}; export _MEMBLE_T0
_STEP=0
_NSTEP=7
step(){ _STEP=$((_STEP + 1))
  printf '>>> [%d/%d] %s  (%ds)\n' "$_STEP" "$_NSTEP" "$1" "$(( $(date +%s) - _MEMBLE_T0 ))"; }
elapsed(){ echo $(( $(date +%s) - _MEMBLE_T0 )); }

# must <what this step does> <what to do about it> -- <command...>
# A step that changes the system stops the build when it fails, and it says what
# to do next. A step that only writes a viewer file keeps "|| true". The
# distinction matters because a system whose vsites, water, ions or residue
# numbers were never fixed still starts, still minimizes, and still returns
# numbers that look ordinary.
must(){
  local what=$1 remedy=$2; shift 2
  [ "$1" = "--" ] && shift
  if ! "$@"; then
    stop "$what failed." "$remedy" "$*"
    MEMBLE_DEGRADED="${MEMBLE_DEGRADED}${MEMBLE_DEGRADED:+; }$what"
    export MEMBLE_DEGRADED
  fi
}

# stop <what happened> <what to do> [<command that failed>]
# Every exit of memble goes through here, so every stop carries a remedy.
stop(){
  local what=$1 remedy=$2 cmd=${3:-}
  echo "" >&2
  echo "ERROR: $what" >&2
  [ -z "$cmd" ] || echo "  command: $cmd" >&2
  echo "" >&2
  echo "  What to do:" >&2
  printf '%b\n' "$remedy" | while IFS= read -r line; do echo "    $line" >&2; done
  echo "" >&2
  if [ "${MEMBLE_IGNORE_ERRORS:-0}" = "1" ]; then
    echo "  MEMBLE_IGNORE_ERRORS=1 is set. memble continues with a system that" >&2
    echo "  did not pass this step, and records the step in memble_build.json." >&2
    return 0
  fi
  echo "  Nothing was deleted. The working directory holds the files as they" >&2
  echo "  stood when the step failed, so the state can be inspected." >&2
  echo "" >&2
  echo "MEMBLE STOPPED after $(elapsed) s: $what" >&2
  # The build a balance pass set aside is removed here. A build that ended early
  # leaves one working directory, so nothing that reads the directory back finds
  # a second system.gro from a pass that was abandoned.
  case "${_MEMBLE_BEST_DIR:-}" in
    *_work.balance_keep) rm -rf "$_MEMBLE_BEST_DIR" ;;
  esac
  exit 1
}
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
[ "$SS_MODE" = tm ] || [ -n "$DSSP" ] || [ -n "$SS_OVERRIDE" ] || stop "no usable DSSP was found (DSSP 4.x is not compatible)" \
  "The secondary structure assignment sets the backbone bonded parameters, and\nMartini holds it for the whole run, so memble needs one. Three routes:\n  1. install DSSP 3.x and point at it:   export DSSP=/path/to/mkdssp\n  2. install mdtraj:                     pip install mdtraj\n  3. assign the transmembrane range by hand, with no DSSP at all:\n       export SS_MODE=tm TM_RANGE=619:641\n  4. give the whole string yourself:     export SS_OVERRIDE=CCCHHHH..."
pick_itp(){ ls $1 2>/dev/null | head -1; }
M3_CORE_ITP=${M3_CORE_ITP:-$(pick_itp "$M3_DIR/martini_v3.0.0.itp")}
M3_SOLV_ITP=${M3_SOLV_ITP:-$(pick_itp "$M3_DIR/*solvent*.itp")}
M3_ION_ITP=${M3_ION_ITP:-$(pick_itp "$M3_DIR/*ion*.itp")}
M3_FFBONDED_ITP=${M3_FFBONDED_ITP:-$(pick_itp "$M3_DIR/*ffbonded*.itp")}   # v2 lipids only; optional
for v in M3_CORE_ITP M3_SOLV_ITP M3_ION_ITP; do [ -n "${!v}" ] || stop "$v was not found under $M3_DIR" \
  "memble looks for the Martini 3 core, solvent and ion itp files in M3_DIR.\nCheck what is there:\n  ls \$M3_DIR/*.itp\nIf the files sit under a different name, set the variable directly, for example\n  export M3_SOLV_ITP=\$M3_DIR/martini_v3.0.0_solvents_v1.itp"; done

# --- route each lipid to the itp that defines it ---
find_itp_for_mol(){ local m=$1 f; for f in "$M3_DIR"/*.itp; do
  awk -v mol="$m" '/^\[ *moleculetype *\]/{g=1;next} g&&NF&&$1!~/^;/{if($1==mol)fd=1;g=0} END{exit !fd}' "$f" && { echo "$f"; return 0; }; done; return 1; }
UNIQ_SRC=()
in_list(){ local x=$1; shift; case " $* " in *" $x "*) return 0;; *) return 1;; esac; }
local_for(){ echo "local_$(basename "$1")"; }
for nm in "${ALL[@]}"; do
  src=$(find_itp_for_mol "$nm") || stop "the lipid '$nm' is in no itp under $M3_DIR" \
  "Every lipid named in LIPIDS, UPPER or LOWER needs a Martini 3 topology.\nCheck the spelling against what the lipidome holds:\n  grep -h moleculetype -A2 \$M3_DIR/*.itp | grep -i $nm\nA lipid that exists only as a published itp is used by copying that itp into\nM3_DIR; COBY imports its structure from the topology."
  in_list "$src" "${UNIQ_SRC[@]}" || UNIQ_SRC+=("$src")
done

# Two files of a lipidome can define the same moleculetype. Including both gives
# a topology that GROMACS refuses, and the refusal names a molecule the build
# never asked for. memble chooses the smallest set of files that covers the
# requested lipids and that defines nothing twice.
_CONS=$("$PY" - "$M3_DIR" "${ALL[*]}" <<'CONSEOF'
import glob, os, sys
d, want = sys.argv[1], sys.argv[2].split()
defs = {}
for f in sorted(glob.glob(os.path.join(d, "*.itp"))):
    names, sec = set(), None
    try:
        fh = open(f, errors="replace")
    except OSError:
        continue
    with fh:
        for raw in fh:
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("["):
                sec = line.strip("[] ").lower()
                continue
            if sec == "moleculetype":
                names.add(line.split()[0]); sec = None
    if names:
        defs[f] = names
need = set(want)
chosen = []
while need:
    best, cover = None, 0
    for f, names in defs.items():
        if f in chosen:
            continue
        if any(defs[c] & names for c in chosen):
            continue
        n = len(names & need)
        if n > cover:
            best, cover = f, n
    if best is None:
        left = sorted(need)
        blocked = [f for f, n in defs.items() if n & need]
        print("COLLIDE|%s|%s" % (",".join(left), ",".join(os.path.basename(x) for x in blocked[:4])))
        sys.exit(0)
    chosen.append(best); need -= defs[best]
for f in chosen:
    print(f)
CONSEOF
)
case "$_CONS" in
  COLLIDE\|*)
    _left=$(printf '%s' "$_CONS" | cut -d'|' -f2)
    _files=$(printf '%s' "$_CONS" | cut -d'|' -f3)
    stop "the lipids $_left cannot be taken from $M3_DIR without defining a molecule twice" \
      "Two files of the lipidome define the same moleculetype, and a topology that\nincludes both is refused by GROMACS.\nThe files that carry these lipids are: $_files\n  1. choose lipids that one of those files supplies on its own\n  2. or move the file you do not need out of \$M3_DIR\nList what each file defines:\n  for f in \$M3_DIR/*.itp; do echo \"== \$f\"; awk '/^\\[ *moleculetype/{g=1;next} g&&NF&&\$1!~/^;/{print \"   \"\$1;g=0}' \"\$f\"; done"
    ;;
  "") ;;
  *)
    UNIQ_SRC=()
    while IFS= read -r _f; do [ -n "$_f" ] && UNIQ_SRC+=("$_f"); done <<< "$_CONS"
    echo ">>> itp sources: ${#UNIQ_SRC[@]} file(s)"
    for _f in "${UNIQ_SRC[@]}"; do echo "      $(basename "$_f")"; done
    ;;
esac
# Every later lookup has to resolve to the files that were chosen, and not to
# the first file of the lipidome that happens to define the molecule.
find_itp_for_mol(){ local m=$1 f; for f in "${UNIQ_SRC[@]}"; do
  awk -v mol="$m" '/^\[ *moleculetype *\]/{g=1;next} g&&NF&&$1!~/^;/{if($1==mol)fd=1;g=0} END{exit !fd}' "$f" && { echo "$f"; return 0; }; done; return 1; }

# ====================================================================
# 1. orient TM along z
# ====================================================================
step "orienting the protein on the transmembrane range"
if [ "$PREBUILT_MULTI" = 1 ]; then
  [ -n "$HELPER_PRE" ] || stop "PREBUILT_MULTI=1 needs HELPER_PRE" \
  "Point HELPER_PRE at prebuild_orient.py in the memble directory:\n  export HELPER_PRE=/path/to/memble/prebuild_orient.py\nsetup.sh sets every HELPER_ variable at once; source it instead of setting them\none at a time."
  [ -n "$RES_KEEP" ]   || stop "PREBUILT_MULTI=1 needs RES_KEEP" \
  "RES_KEEP names the residues to keep from a prebuilt multimer, one range per\nchain, for example\n  export RES_KEEP='A:54-103;B:54-103'\nThe chain letters are the ones in the input PDB."
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
    [ -n "$ORI_SS" ] || stop "MULTI_TM=1 needs a secondary structure string and DSSP returned none" \
  "MULTI_TM=1 orients a multi-pass protein on its helix bundle, so it needs to\nknow which residues are helical. Either\n  export SS_OVERRIDE=CCCHHHHH...    (one character per residue)\nor install a working DSSP (3.x) or mdtraj and try again."
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
step "assigning the secondary structure and coarse-graining the protein"
MZ=(-ff martini3001 -f oriented_aa.pdb -x cg_peptide.pdb -o protein_only.top -p backbone -cys auto -maxwarn 10)
# Secondary structure source (SS_OVERRIDE > SS_MODE=tm > SS_MODE=dssp):
if [ -n "$SS_OVERRIDE" ]; then
  MZ+=(-ss "$SS_OVERRIDE")
elif [ "$SS_MODE" = tm ]; then
  # gen_ss.py wants CHAIN:start-end, but TM_RANGE is written start:end on the
  # command line, so ALL:$TM_RANGE would hand it "ALL:65:88" and it would
  # reject the range. Convert the separator instead of failing.
  SS_TM=${TM_CORE:-}; [ -n "$SS_TM" ] || { [ -n "$TM_RANGE" ] && SS_TM="ALL:${TM_RANGE/:/-}"; }
  [ -n "$SS_TM" ] || stop "SS_MODE=tm needs a transmembrane range" \
  "Give the residue numbers of the transmembrane helix, in the numbering of the\ninput PDB:\n  export TM_RANGE=619:641                 (one chain)\n  export TM_CORE='A:619-641;B:619-641'    (several chains)\nTM_RANGE uses a colon, TM_CORE uses a hyphen inside each chain."
  SS_STR=$("$PY" "$HELPER_GENSS" --pdb oriented_aa.pdb --tm "$SS_TM") || stop "the secondary structure string could not be generated from $SS_TM" \
  "gen_ss.py wants CHAIN:start-end. Check that the residue numbers exist in the\ninput PDB and that the chain letter matches:\n  grep '^ATOM' oriented_aa.pdb | cut -c22-26 | sort -un | head\n  grep '^ATOM' oriented_aa.pdb | cut -c22 | sort -u"
  [ -n "$SS_STR" ] || stop "SS_MODE=tm produced an empty secondary structure string" \
  "No residue matched the range $SS_TM. Check the residue numbering of the input\nPDB, which memble does not renumber:\n  grep '^ATOM' oriented_aa.pdb | cut -c22-26 | sort -un | head"
  echo ">>> SS_MODE=tm: TM=$SS_TM -> SS length ${#SS_STR} (TM helix, rest coil)"
  MZ+=(-ss "$SS_STR")
elif [ "$SS_MODE" = dssp-internal ]; then
  # Old behavior, kept as an escape hatch: martinize2 runs DSSP itself and the
  # string it used is never written down. A build made this way cannot be
  # repeated from the output alone, and verify_system.py reports that.
  echo ">>> SS_MODE=dssp-internal: martinize2 runs DSSP; the assignment is not recorded"
  if [ "$DSSP" = mdtraj ]; then MZ+=(-dssp); else MZ+=(-dssp "$DSSP"); fi
else
  # Default: run DSSP here, keep the string, and hand it to martinize2. The
  # assignment sets the backbone bonded parameters and Martini holds it for the
  # whole run, so the string belongs in the output next to the composition.
  if [ "$DSSP" = mdtraj ] || [ -z "$DSSP" ]; then
    SS_STR=$("$PY" "$HELPER_SSDSSP" --pdb oriented_aa.pdb) \
      || stop "DSSP failed on oriented_aa.pdb" \
        "memble runs DSSP itself so that the assignment is recorded. Three ways on:\n  1. install one:            pip install mdtraj      (or install DSSP 3.x)\n  2. assign by hand:         export SS_MODE=tm TM_RANGE=619:641\n  3. give the whole string:  export SS_OVERRIDE=CCCHHHH...\nTo let martinize2 run DSSP the old way, without recording the string:\n  export SS_MODE=dssp-internal"
  else
    SS_STR=$("$PY" "$HELPER_SSDSSP" --pdb oriented_aa.pdb --dssp "$DSSP") \
      || stop "DSSP ($DSSP) failed on oriented_aa.pdb" \
        "Check that the binary runs and is version 3.x; 4.x writes a different format:\n  \$DSSP --version\nOtherwise assign the transmembrane range by hand, with no DSSP at all:\n  export SS_MODE=tm TM_RANGE=619:641"
  fi
  [ -n "$SS_STR" ] || stop "DSSP returned an empty secondary structure string" \
  "DSSP ran and assigned nothing, so it read no protein residue in oriented_aa.pdb.\nCheck that the file holds ATOM records with backbone atoms:\n  grep -c '^ATOM' oriented_aa.pdb\n  grep '^ATOM' oriented_aa.pdb | cut -c13-16 | sort -u | head"
  echo ">>> SS_MODE=dssp: SS length ${#SS_STR}, helical residues $(printf '%s' "$SS_STR" | tr -cd 'H' | wc -c | tr -d ' ')"
  MZ+=(-ss "$SS_STR")
fi
# Record the assignment that martinize2 actually received, whatever produced it.
SS_USED=""; SS_SOURCE=""
if [ -n "$SS_OVERRIDE" ]; then SS_USED="$SS_OVERRIDE"; SS_SOURCE="override"
elif [ "$SS_MODE" = tm ]; then SS_USED="$SS_STR"; SS_SOURCE="tm"
elif [ "$SS_MODE" = dssp-internal ]; then SS_USED=""; SS_SOURCE="dssp-internal"
else SS_USED="$SS_STR"; SS_SOURCE="dssp"
fi
export SS_USED SS_SOURCE
[ "$WATER_BIAS" -eq 1 ] && MZ+=(-water-bias -water-bias-eps E:-0.5 C:1.0 H:-1.0)
"$MARTINIZE2" "${MZ[@]}"
res_count_itp(){ awk '/^\[/{a=($2=="atoms")?1:0;next} a&&NF>0&&$1!~/^;/{print $3}' "$1" | sort -un | wc -l | tr -d ' '; }
if [ "$PREBUILT_MULTI" = 1 ]; then
  PROT_ITPS=(molecule_*.itp)
  PROT_MOLMAP=$(awk '/\[ *molecules *\]/{m=1;next} m&&/^\[/{m=0} m&&NF&&$1!~/^;/{printf (n++?":":"") $1} END{printf "\n"}' protein_only.top)
  [ -n "$PROT_MOLMAP" ] || stop "the [ molecules ] section of protein_only.top could not be read" \
  "martinize2 wrote protein_only.top and memble found no molecule in it, so the\ncoarse-graining produced nothing. Read the martinize2 output above for the\nreason, and check the input:\n  head -30 protein_only.top\nA PDB with alternate locations or with no CA atoms is the usual cause."
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
# box_z: membrane thickness + water per side. NOTE: COBY's lipid-grid optimizer
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
  [ -n "$HELPER_PART" ] || stop "PARTNER_PDB is set and HELPER_PART is not" \
  "Point HELPER_PART at place_partner.py in the memble directory:\n  export HELPER_PART=/path/to/memble/place_partner.py\nsetup.sh sets every HELPER_ variable at once."
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
step "packing the lipids around the protein"
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
# A measured value replaces the tabulated one. The balance pass below sets these
# from the areas it measured on a first build, so the counts stop depending on a
# table of areas per lipid.
[ -z "${APL_UPPER:-}" ] || { APL_UP="$APL_UPPER"; echo ">>> upper apl set to ${APL_UP} from a measurement"; }
[ -z "${APL_LOWER:-}" ] || { APL_LO="$APL_LOWER"; echo ">>> lower apl set to ${APL_LO} from a measurement"; }
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
if [ "${MEMBLE_NO_TM_CENTER:-0}" = "1" ]; then
  # Reproduces the placement a build gives when COBY is not told which residues
  # cross the membrane. Kept so the comparison of Section 3.1 can be repeated.
  CEN_RES=""
  echo ">>> MEMBLE_NO_TM_CENTER=1: COBY centers the protein on the whole molecule"
elif [ -z "$CEN_RES" ] && [ -n "$TM_RANGE" ] && [ "$PREBUILT_MULTI" != 1 ]; then
  # The single-chain path also has to tell COBY where the membrane-spanning part
  # is. Without it COBY centers the protein on the centroid of the whole
  # molecule, and a construct whose extramembrane parts differ in length between
  # the two sides then sits with its transmembrane helix off the bilayer
  # midplane by half that difference. The system builds, minimizes and runs.
  _first=$(awk '/^ATOM/{r=substr($0,23,4)+0; print r; exit}' oriented_aa.pdb 2>/dev/null)
  _tl=${TM_RANGE%%:*}; _th=${TM_RANGE##*:}
  if [ -n "$_first" ] && [ -n "$_tl" ] && [ -n "$_th" ] && [ -n "$N_TMJM" ]; then
    _ranges=""; _off=0; _i=0
    while [ "$_i" -lt "$N_COPY" ]; do
      _a=$((_off + _tl - _first)); _b=$((_off + _th - _first))
      _ranges="${_ranges:+${_ranges}:}${_a}-${_b}"
      _off=$((_off + N_TMJM)); _i=$((_i + 1))
    done
    CEN_RES="$_ranges"
    echo ">>> COBY will center the protein on the transmembrane range $TM_RANGE (assembly resid ${CEN_RES})"
  fi
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
must "renaming the ion beads to the Martini names" \
  "COBY wrote NA+ and CL-, and the Martini itp files call them NA and CL.\nsed could not write system.gro or system.top.\nCheck that the output directory is writable and that the disk is not full." -- \
  sed_inplace 's/NA+/NA /g; s/CL-/CL /g' system.gro system.top

# Optional fine tuning of how deep the protein sits in the membrane. Z_SHIFT
# (nm, default 0) moves only the protein beads in z, found by trying a few
# values and rebuilding. The declash below relaxes the lipids around the
# shifted protein, and the position restraint reference is read from these
# coordinates, so the offset is held through equilibration.
if [ "${Z_SHIFT:-0}" != "0" ] && [ "${Z_SHIFT:-0}" != "0.0" ]; then
  echo ">>> shifting the protein in z by Z_SHIFT=${Z_SHIFT} nm"
  must "shifting the protein in z by Z_SHIFT=${Z_SHIFT} nm" \
  "Z_SHIFT moves only the protein beads in z.\nCheck that Z_SHIFT is a number in nm, for example Z_SHIFT=0.4 or Z_SHIFT=-0.4.\nSet Z_SHIFT=0 to skip the shift and place the protein where COBY put it." -- \
    "$PY" "$HELPER_ZSHIFT" --gro system.gro --dz "$Z_SHIFT"
fi

# ====================================================================
# 4a. DECLASH: push apart inter-molecular bead overlaps from dense packing
#     (rigid-molecule moves; intramolecular geometry preserved) so that
#     minimization does not hit infinite Lennard-Jones forces.
# ====================================================================
# COBY places coarse-grained lipids by their real beads but leaves out-of-plane
# virtual sites (sterol ROH/R3, funct 4) flat; GROMACS then reconstructs them
# ~0.1 nm away and they can explode at minimization. Rebuild every vsite exactly
# from its itp definition before declashing/minimizing.
must "rebuilding the sterol virtual sites from their itp definitions" \
  "COBY leaves out-of-plane virtual sites (sterol ROH and R3, funct 4) flat, and\nGROMACS rebuilds them about 0.1 nm away, which explodes at minimization.\nCheck that M3_DIR points at the Martini 3 lipidome directory and that it holds\nthe itp of every sterol in the composition:  ls $M3_DIR/*.itp | grep -i chol\nA composition without a sterol does not need this step; remove the sterol or\nadd its itp to M3_DIR." -- \
  "$PY" "$HELPER_FIXVS" --gro system.gro --top system.top --itp-dir "$M3_DIR"
must "separating overlapping beads left by the packing" \
  "The packing left two beads of different molecules on top of each other.\nRaise BOX_X and BOX_Y by 1 nm and build again, which gives the packing room.\nOr lower the lipid density with a larger COBY_APL.\nTo keep the system and look at it, set MEMBLE_ALLOW_OVERLAP=1." -- \
  "$PY" "$HELPER_DECLASH" --gro system.gro --lipids "${ALL[*]}" --target 0.21 --iters 200 --exclude-beads "ROH R3"

# ====================================================================
# 4b. ADD WATER: COBY builds the membrane in a thin box (a large box_z hangs its
#     lipid-grid optimizer). Now grow box_z so a protruding TM-JM protein gets a
#     real bulk-water cushion, and solvate the new slabs (+ salt). No-op if the
#     protein already fits with the requested water on each side.
# ====================================================================
step "sizing the box and adding water and ions"
if [ -n "$PARTNER_CG" ]; then
  # peripheral protein: place on the chosen leaflet, expand box, then fill water
  "$PY" "$HELPER_PART" --system-gro system.gro --system-top system.top \
      --partner-cg "$PARTNER_CG" --partner-name "$PARTNER_NAME" \
      --partner-itp "$PARTNER_ITP" --side "$PARTNER_SIDE" --gap "$PARTNER_GAP" \
      --water-nm "$PARTNER_WATER" --lipids "${ALL[*]}" \
      --rotate "$PARTNER_ROTATE" --seed "$SEED"
  if [ -n "$HELPER_ADDWATER" ]; then
    must "adding water and ions around the peripheral protein" \
  "The peripheral protein was placed, and the solvation of the new box failed.\nCheck PARTNER_WATER (nm, for example 2.5) and SALT_M (molar, for example 0.15).\nCheck that the Martini water itp is in M3_DIR:  ls $M3_DIR | grep -i water" -- \
      "$PY" "$HELPER_ADDWATER" --gro system.gro --top system.top \
        --water-nm "$PARTNER_WATER" --salt "$SALT_M" --keep-box
  fi
elif [ -n "$HELPER_ADDWATER" ]; then
  # Make the protein contiguous BEFORE sizing the box. add_water derives
  # box_z from the protein z-extent it can see; while a chain is still split
  # across the z boundary that extent is the wrapped one, which for a tall
  # TM-JM assembly is far smaller than the real span (6.95 vs 13.73 nm in the
  # four-copy EGFR build), so the box came out ~7 nm too short and the protein
  # ended up overlapping its own periodic image. Unwrap first, then size.
  # The call after the solvation step below stays: it re-centers in the final
  # box and reports the resulting water cushion.
  if [ -n "$HELPER_WHOLE" ] && [ -f "$HELPER_WHOLE" ]; then
    must "making the protein contiguous before the box is sized" \
  "add_water measures the z span of the protein to size the box. While a chain is\nstill split across the z boundary that span is the wrapped one, which is far\nshorter than the protein, and the box comes out too short.\nCheck that the itp files of the protein are in the working directory:  ls molecule_*.itp\nIf the protein has one chain and does not cross the boundary, this step can be\nskipped with HELPER_WHOLE= (empty)." -- \
      "$PY" "$HELPER_WHOLE" --gro system.gro --top system.top --itp-dir .
  fi
  must "sizing the box and adding water and ions" \
  "Check WATER_NM (nm per side, for example 2.5) and SALT_M (molar, for example 0.15).\nCheck that the Martini water and ion itp files are in M3_DIR.\nA very tall protein in a small xy box can leave no room for water; raise BOX_X\nand BOX_Y, or lower WATER_NM." -- \
    "$PY" "$HELPER_ADDWATER" --gro system.gro --top system.top \
      --water-nm "$WATER_NM" --salt "$SALT_M"
fi

# final PBC-aware declash: catches any bead sitting just outside the box that
# clashes with the opposite face under periodic boundaries (a common COBY edge
# effect that produces an infinite force on a water at minimization).
# Make the protein contiguous across PBC and center the system on it, BEFORE the
# final declash. A tall TM-JM protein can straddle the z boundary after
# assembly/solvation; per-atom wrapping upstream then splits a chain across the
# box (consecutive backbone beads a full box apart -> LINCS blowup -> infinite
# force). Recentering can push some lipids/water across the edge, so the final
# declash must run AFTER this, as the last coordinate step, to clean up.
if [ -n "$HELPER_WHOLE" ] && [ -f "$HELPER_WHOLE" ]; then
  must "making the protein contiguous in the final box" \
  "Recentering can split a chain across the box edge, and two consecutive backbone\nbeads a full box apart give an infinite force at minimization.\nCheck that the protein itp files are in the working directory:  ls molecule_*.itp" -- \
    "$PY" "$HELPER_WHOLE" --gro system.gro --top system.top --itp-dir . --box-is-final
fi

# FINAL declash and the last coordinate-modifying step: guarantees no two beads
# from different molecules overlap. The protein is frozen (never moved), so it
# is not distorted or re-split; only solvent and lipids are pushed apart.
must "the final separation of overlapping beads" \
  "Recentering pushed some lipids or water across the box edge and they now overlap.\nRaise BOX_X and BOX_Y by 1 nm and build again.\nTo keep the system and look at it, set MEMBLE_ALLOW_OVERLAP=1." -- \
  "$PY" "$HELPER_DECLASH" --gro system.gro --lipids "${ALL[*]}" --target 0.21 --iters 200 --exclude-beads "ROH R3" --freeze-protein

# Confirm there is no residual overlap that would give an infinite force. A
# warning here was read past and the build shipped, so this now stops the build.
# MEMBLE_ALLOW_OVERLAP=1 keeps the old behavior for a system that is being
# inspected rather than run.
if [ -n "$HELPER_MINDIST" ] && [ -f "$HELPER_MINDIST" ]; then
  if ! "$PY" "$HELPER_MINDIST" --gro system.gro --lipids "${ALL[*]}" --min 0.12; then
    if [ "${MEMBLE_ALLOW_OVERLAP:-0}" = "1" ]; then
      echo ">>> MEMBLE_ALLOW_OVERLAP=1: continuing with a residual overlap."
    else
      stop "beads of different molecules are closer than 0.12 nm" \
        "Minimization of this system reports an infinite force, or it moves the two\nmolecules far enough apart to distort them. The pair is named just above.\n  1. give the packing room:   raise BOX_X and BOX_Y by 1 nm\n  2. lower the lipid density: raise COBY_APL\n  3. if the pair involves the protein, raise SPACING_NM or lower N_COPY\n  4. to keep this system and look at it:  export MEMBLE_ALLOW_OVERLAP=1"
    fi
  fi
fi

# Restore the protein's original per-chain residue numbers (martinize renumbers
# every chain from 1, so the assembled chains overlap). The true numbers come
# from oriented_aa.pdb (which preserves the input PDB resSeq); nothing hardcoded.
if [ -f oriented_aa.pdb ] && [ -n "$HELPER_FIXRESID" ]; then
  must "restoring the per-chain residue numbers of the input PDB" \
  "martinize2 renumbers every chain from 1, so the chains of a multi-chain protein\noverlap. The true numbers come from oriented_aa.pdb.\nCheck that oriented_aa.pdb exists and holds the input residue numbers.\nThe system is correct without this step; only the numbering in system.gro is\nmartinize2 numbering. Set HELPER_FIXRESID= (empty) to accept that." -- \
    "$PY" "$HELPER_FIXRESID" --gro system.gro --oriented oriented_aa.pdb \
        --top system.top --itp-dir .
fi

# ====================================================================
# 4b. LEAFLET AREA PRE-CHECK (before any gmx MD; abort the build if asymmetric
#     leaflets are area-mismatched, so it is fixed now, not after a melted run)
# ====================================================================
step "measuring the area of every lipid in each leaflet"
"$PY" "$HELPER_AREA" --gro system.gro --lipids "${ALL[*]}" --asym "$ASYM" --tol "$AREA_TOL" --hard-tol "$AREA_HARD_TOL" --json leaflet_area.json ${APL_OVERRIDE:+--apl "$APL_OVERRIDE"}

# ====================================================================
# 4c. BALANCE PASS
#     The number of lipids a leaflet receives came from a table of areas per
#     lipid, and a table cannot know what a mixture with a sterol and a protein
#     in it will do. This pass reads the areas that were just measured, corrects
#     the area per lipid of each leaflet by what the measurement says, and builds
#     the system again. MEMBLE_BALANCE_ITER=0 keeps the first build.
#
#     Three rules keep the pass from making the system worse. It corrects only a
#     difference that the check itself would fail, so a difference of the size of
#     its own standard error is left alone. It moves the areas per lipid half way
#     to the correction, so one pass cannot overshoot. It keeps every build it
#     makes, compares each new difference against the smallest one so far, and
#     restores the earlier build when a pass does not lower the difference.
# ====================================================================
_ITER=${_MEMBLE_ITER:-0}
_KEEP="$_START_DIR/${OUTTAG}_work.balance_keep"
_BAL_LOG=${_MEMBLE_BAL_LOG:-}
if [ "$ASYM" = 1 ] && [ -f leaflet_area.json ]; then
  _FIX=$("$PY" - "$APL_UP" "$APL_LO" "${MEMBLE_BALANCE_TOL:-$AREA_TOL}" <<'BALEOF'
import json, os, sys
up, lo, tol = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
damp = float(os.environ.get("MEMBLE_BALANCE_DAMP", "0.5"))
try:
    d = json.load(open("leaflet_area.json"))
except Exception:
    sys.exit(0)
sh = d.get("shared_species", [])
if not sh:
    sys.exit(0)
# The lipid whose two areas differ by the most, and the standard error of that
# difference. A difference no larger than twice its own standard error is what
# the measurement returns on a balanced membrane, and it is not corrected.
w = max(sh, key=lambda x: x["relative_difference"])
worst, worst_e = w["relative_difference"], w["relative_standard_error"]
# The area a shared lipid takes in each leaflet, averaged over the shared
# lipids. A leaflet whose lipids are larger than the other holds too few of
# them, and its area per lipid has to come down by that ratio.
au = sum(x["upper_nm2"] for x in sh) / len(sh)
al = sum(x["lower_nm2"] for x in sh) / len(sh)
ok = 1
if d.get("shared_fraction", 0.0) < 0.5: ok = 0
if worst <= tol or worst <= 2.0 * worst_e: ok = 0
if au <= 0 or al <= 0: ok = 0
mid = 0.5 * (au + al)
nu = up * (1.0 + damp * (mid / au - 1.0)) if au > 0 else up
nl = lo * (1.0 + damp * (mid / al - 1.0)) if al > 0 else lo
print("%.4f %.4f %.4f %.4f %d" % (worst, worst_e, nu, nl, ok))
BALEOF
)
  if [ -n "$_FIX" ]; then
    _W=$( printf '%s' "$_FIX" | awk '{print $1}')
    _WE=$(printf '%s' "$_FIX" | awk '{print $2}')
    _NU=$(printf '%s' "$_FIX" | awk '{print $3}')
    _NL=$(printf '%s' "$_FIX" | awk '{print $4}')
    _OK=$(printf '%s' "$_FIX" | awk '{print $5}')
    _WS=$(awk -v x="$_W" 'BEGIN{printf "%.1f", 100*x}')
    _BEST=${_MEMBLE_BEST_DIFF:-}
    _BAL_LOG="${_BAL_LOG:+$_BAL_LOG }pass${_ITER}:${_WS}%"
    # Is this build the best one so far?
    if [ -z "$_BEST" ] || awk -v a="$_W" -v b="$_BEST" 'BEGIN{exit !(a<b)}'; then
      _IS_BEST=1
    else
      _IS_BEST=0
    fi
    if [ "$_IS_BEST" = 1 ] && [ "$_OK" = 1 ] && [ "$_ITER" -lt "${MEMBLE_BALANCE_ITER:-0}" ]; then
      echo ""
      echo ">>> balance pass $((_ITER + 1)): the leaflets differ by ${_WS}%, above the"
      echo "    tolerance of $(awk -v x="${MEMBLE_BALANCE_TOL:-$AREA_TOL}" 'BEGIN{printf "%.1f", 100*x}')% and above twice the standard error of $(awk -v x="$_WE" 'BEGIN{printf "%.1f", 100*x}')%."
      echo "    The areas per lipid are corrected half way to the measurement:"
      echo "      upper ${APL_UP} -> ${_NU}"
      echo "      lower ${APL_LO} -> ${_NL}"
      echo "    This build is kept, and memble builds the system again."
      cd "$_START_DIR" || cd ..
      rm -rf "$_KEEP"; mv "$WORK" "$_KEEP"
      _MEMBLE_ITER=$((_ITER + 1)) _MEMBLE_BEST_DIFF="$_W" _MEMBLE_BEST_DIR="$_KEEP" \
      _MEMBLE_BAL_LOG="$_BAL_LOG" APL_UPPER="$_NU" APL_LOWER="$_NL" \
        exec bash "$_SELF" "$_ARG1"
    fi
    if [ "$_IS_BEST" = 0 ] && [ -d "${_MEMBLE_BEST_DIR:-/nonexistent}" ]; then
      _BS=$(awk -v x="$_BEST" 'BEGIN{printf "%.1f", 100*x}')
      echo ""
      echo ">>> the correction of pass ${_ITER} left the leaflets ${_WS}% apart, against"
      echo "    ${_BS}% before it, so the correction did not improve the membrane."
      echo "    memble restores the build that gave ${_BS}% and stops correcting."
      cd "$_START_DIR" || cd ..
      rm -rf "$WORK"; mv "$_MEMBLE_BEST_DIR" "$WORK"; cd "$WORK"
      _BAL_LOG="${_BAL_LOG} restored:pass$((_ITER - 1))"
    fi
  fi
fi
rm -rf "$_KEEP"

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
  [ "$bg" -eq "$bi" ] && [ "$bg" -ne 0 ] || stop "the lipid $nm has $bg beads in system.gro and $bi beads in its topology" \
  "The structure that COBY placed and the itp that the topology includes describe\ndifferent molecules, so GROMACS would read the wrong beads.\nUsually two itp files in M3_DIR define the same moleculetype with different bead\ncounts. Find them:\n  grep -l \"^ *$nm \" \$M3_DIR/*.itp\nKeep one and move the other out of M3_DIR."
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
# ParmEd writes the PSF, the CRD and a PDB for a viewer. The system runs without
# them, so a missing ParmEd is reported and the build carries on.
if ! "$PY" -c 'import parmed' >/dev/null 2>&1; then
  echo ">>> ParmEd is not installed, so system.psf, system.crd and system_parmed.pdb"
  echo "    were not written. The system itself is complete and runs without them."
  echo "    To get the viewer files:  pip install ParmEd"
else
if ! "$PY" - "${ALL[*]}" "$PROT_BLOCKS" "${PART_COUNTS[*]}" <<'PYEOF'
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
then :; else
  echo ">>> ParmEd could not write system.psf, system.crd and system_parmed.pdb."
  echo "    Those files are for a viewer. The system itself is complete and runs"
  echo "    without them, so the build carries on."
fi
fi
"$GMX" editconf -f system.gro -o system.pdb >/dev/null 2>&1 || true

# ====================================================================
# 8. mdp (GROMACS 2023.x; TEMP-parameterized) + run.sh
# ====================================================================
step "writing the staged equilibration"
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
# ====================================================================
# 5b. BUILD RECORD
#     One file that says what this system is and how it was made, so the
#     directory alone answers the question six months later.
# ====================================================================
{
  printf '{\n'
  printf '  "memble_version": "%s",\n' "$MEMBLE_VERSION"
  printf '  "built_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "input_pdb": "%s",\n' "$PEP_AA"
  printf '  "input_pdb_sha256": "%s",\n' "$( { sha256sum "$PEP_AA" 2>/dev/null || shasum -a 256 "$PEP_AA" 2>/dev/null; } | awk '{print $1}')"
  printf '  "ss_source": "%s",\n' "$SS_SOURCE"
  printf '  "ss_string": "%s",\n' "$SS_USED"
  printf '  "tm_range": "%s",\n' "$TM_RANGE"
  printf '  "tm_core": "%s",\n' "${TM_CORE:-}"
  printf '  "n_copy": "%s",\n' "$N_COPY"
  printf '  "lipids": "%s",\n' "$LIPIDS"
  printf '  "upper": "%s",\n' "${UPPER:-}"
  printf '  "lower": "%s",\n' "${LOWER:-}"
  printf '  "box_x_nm": "%s",\n' "${BOX_X:-}"
  printf '  "box_y_nm": "%s",\n' "${BOX_Y:-}"
  printf '  "water_nm": "%s",\n' "$WATER_NM"
  printf '  "salt_M": "%s",\n' "$SALT_M"
  printf '  "temperature_K": "%s",\n' "$TEMP"
  printf '  "m3_dir": "%s",\n' "$M3_DIR"
  printf '  "gromacs": "%s",\n' "$("$GMX" --version 2>/dev/null | awk -F': *' '/GROMACS version/{print $2; exit}')"
  printf '  "martinize2": "%s",\n' "$("$MARTINIZE2" --version 2>&1 | head -1 | tr -d '"')"
  printf '  "coby": "%s",\n' "$("$PY" -c 'import COBY;print(getattr(COBY,"__version__","unknown"))' 2>/dev/null)"
  printf '  "balance_passes": "%s",\n' "${_BAL_LOG:-}"
  printf '  "degraded_steps": "%s"\n' "${MEMBLE_DEGRADED:-}"
  printf '}\n'
} > memble_build.json
echo ">>> build record written to memble_build.json"

# ====================================================================
# 6. VERIFICATION GATE
#    Every post-condition that a broken system still survives is measured here.
#    memble writes run.sh only after the gate passes, so a build that did not
#    pass is never handed over as one that is ready to run.
# ====================================================================
step "measuring the eight properties of the finished system"
if [ -f "$HELPER_VERIFY" ]; then
  VERIFY=(--gro system.gro --top system.top --itp-dir "$M3_DIR" --itp-dir .
          --lipids "${ALL[*]}" --water-nm "$WATER_NM"
          --ss-mode "$SS_SOURCE" --meta memble_build.json
          --leaflet-json leaflet_area.json)
  [ -n "$SS_USED" ]  && VERIFY+=(--ss-string "$SS_USED")
  [ -n "$TM_RANGE" ] && VERIFY+=(--tm-resids "$TM_RANGE")
  [ -n "$UPPER" ]    && VERIFY+=(--expect-upper "$UPPER")
  [ -n "$LOWER" ]    && VERIFY+=(--expect-lower "$LOWER")
  if [ -z "$UPPER$LOWER" ] && [ -n "$LIPIDS" ]; then
    VERIFY+=(--expect-upper "$LIPIDS" --expect-lower "$LIPIDS")
  fi
  for c in ${MEMBLE_ALLOW:-}; do VERIFY+=(--allow "$c"); done
  if ! "$PY" "$HELPER_VERIFY" "${VERIFY[@]}"; then
    stop "the finished system did not pass verification" \
      "Every check that failed is printed above and in memble_report.txt, and each\none carries what to do about it. system.gro and system.top are kept so the\nsystem can be looked at; run.sh was not written, so nothing runs by accident.\nTo accept one named check and continue:\n  export MEMBLE_ALLOW=\"water_layer overlap\""
  fi
  echo ">>> verification passed; memble_report.txt and memble_report.json written"
fi

cat > run.sh <<RUNEOF
#!/usr/bin/env bash
set -eo pipefail
GMX="$GMX"
# One thread-MPI rank and NT OpenMP threads. A system of this size needs no
# domain decomposition, and GROMACS otherwise takes one rank per GPU it can see,
# which stops the run on a machine that carries more than one. Set NT to the
# number of cores to use:  NT=16 bash run.sh
NT=\${NT:-8}
# GROMACS dumps step<N>b.pdb when atoms move too far (system blowing up).
# set -e misses a run that "succeeds" while melting, so check explicitly.
check_blowup(){ if ls step*[0-9]b.pdb >/dev/null 2>&1; then
  echo "INSTABILITY during \$1: GROMACS wrote step*b.pdb (atoms moving too far)."
  echo "  The membrane is blowing up. Inspect: clashes (EM max force), box too small,"
  echo "  or asymmetric leaflet area mismatch. Do not continue."; exit 1; fi; }
# A stage whose coordinates are already written is left as it stands, and a
# stage that was stopped partway continues from the checkpoint GROMACS wrote,
# so a run that was interrupted is picked up by starting this script again.
# REDO=1 runs every stage from the beginning:  REDO=1 bash run.sh
stage(){   # stage <name> <mdp> <start.gro> [restraint.gro]
  out=\$1; mdp=\$2; start=\$3; restr=\${4:-}; nx=""
  [ -f index.ndx ] && nx="-n index.ndx"
  if [ -f "\$out.gro" ] && [ "\${REDO:-0}" != "1" ]; then
    echo ">>> \$out was already finished"; return 0; fi
  if [ ! -f "\$out.tpr" ] || [ "\${REDO:-0}" = "1" ]; then
    if [ -n "\$restr" ]; then
      \$GMX grompp -f "\$mdp" -c "\$start" -r "\$restr" -p system.top \$nx -o "\$out.tpr" -maxwarn 10
    else
      \$GMX grompp -f "\$mdp" -c "\$start" -p system.top \$nx -o "\$out.tpr" -maxwarn 10
    fi
  fi
  if [ -f "\$out.cpt" ] && [ "\${REDO:-0}" != "1" ]; then
    echo ">>> \$out continues from the checkpoint it left"
    \$GMX mdrun -deffnm "\$out" -v -ntmpi 1 -ntomp \$NT -cpi "\$out.cpt" -append
  else
    \$GMX mdrun -deffnm "\$out" -v -ntmpi 1 -ntomp \$NT
  fi
  check_blowup "\$out"
}
stage step6.0 step6.0_minimization.mdp system.gro
grep -i "Maximum force" step6.0.log | tail -1 || true   # should be finite, not astronomical
prev=step6.0
for k in 1 2 3 4 5 6; do
  stage step6.\${k} step6.\${k}_equilibration.mdp \${prev}.gro step6.0.gro
  prev=step6.\${k}
done
# The last equilibration stage is read before the production run starts. A
# membrane whose area is still drifting is not equilibrated, whatever the length
# of the stage was, and the production run inherits the drift.
if [ -f "$_HELPDIR/check_equilibration.py" ]; then
  if ! "$PY" "$_HELPDIR/check_equilibration.py" --edr \${prev}.edr --gmx "\$GMX" \
        --json equilibration.json; then
    if [ "\${MEMBLE_ALLOW_DRIFT:-0}" = "1" ]; then
      echo ">>> MEMBLE_ALLOW_DRIFT=1: starting the production run anyway."
    else
      echo "The production run was not started. Extend \${prev} and read it again," >&2
      echo "or set MEMBLE_ALLOW_DRIFT=1 to start regardless." >&2
      exit 1
    fi
  fi
fi
stage step7 step7_production.mdp \${prev}.gro
echo "DONE: step7.xtc"
RUNEOF
chmod +x run.sh
# also drop the staged, stop-on-failure MD runner (CHARMM-GUI style: one stage at
# a time with a success check) next to the system, copied from the helper folder.
if [ -f "$_HELPDIR/run_md.sh" ]; then
  cp "$_HELPDIR/run_md.sh" run_md.sh && chmod +x run_md.sh
fi
for _h in check_equilibration.py leaflet_area_check.py verify_system.py; do
  if [ -f "$_HELPDIR/$_h" ] && [ ! -f "./$_h" ]; then cp "$_HELPDIR/$_h" "./$_h"; fi
done
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
   - A stage that was stopped partway leaves a checkpoint file, <name>.cpt.
     Continue it by adding  -cpi <name>.cpt -append  to that stage's mdrun line
     and running the line again. The grompp line is not run a second time.
     Example, for a production run that was stopped:
       gmx mdrun -deffnm step7_production -v -ntmpi 1 -ntomp 8 \
         -cpi step7_production.cpt -append
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
# One last line, so a build that ran in the background is read from its tail.
_NPASS=$(awk '/^RESULT:/{print $2}' memble_report.txt 2>/dev/null)
echo ""
echo "MEMBLE ${_NPASS:-DONE}: ${OUTTAG} built in $(elapsed) s, all eight properties measured."
echo ">>> Build complete in $WORK."
echo ">>> Equilibrate + produce stage-by-stage:  cd $WORK && bash run_md.sh"
echo ">>> Or copy-paste stages manually from:    $WORK/md_steps.txt"
echo ">>> View WITH BONDS (CHARMM-GUI style):     cd $WORK && vmd -e view.vmd"
echo ">>>   (loads system.psf for connectivity, then the trajectory/gro)"
echo ">>> (legacy all-in-one script:             ./run.sh)"
