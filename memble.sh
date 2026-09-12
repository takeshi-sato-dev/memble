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
# The helpers sit beside memble.sh in the clone, so each one defaults to the
# copy next to this script. Sixteen exported paths were required before, which
# every caller had to repeat and every new machine got wrong once. Setting a
# HELPER_ variable still overrides the file it names, one at a time.
_HERE=$(dirname "$_SELF")
HELPER_REP=${HELPER_REP:-$_HERE/replicate_and_fix_top.py}
HELPER_POS=${HELPER_POS:-$_HERE/inject_posres.py}
HELPER_ORI=${HELPER_ORI:-$_HERE/orient_tm.py}
HELPER_SSDSSP=${HELPER_SSDSSP:-$_HERE/ss_from_dssp.py}
HELPER_AREA=${HELPER_AREA:-$_HERE/leaflet_area_check.py}
HELPER_PART=${HELPER_PART:-$_HERE/place_partner.py}
HELPER_ITP2STRUCT=${HELPER_ITP2STRUCT:-$_HERE/itp_to_struct.py}
for _h in REP POS ORI AREA ITP2STRUCT; do
  _v=HELPER_$_h
  [ -f "${!_v}" ] || { echo "memble: $_v points at a file that is not there: ${!_v}" >&2
    echo "The helpers live beside memble.sh. Check the clone:  ls $_HERE/*.py" >&2; exit 1; }
done
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
# COBY seeds its own random number generator from the wall clock, so two builds
# of one composition return two different packings and a built system cannot be
# built again. COBY_SEED hands COBY an integer instead, and the build is then
# reproducible bead for bead. Left empty, COBY keeps its own default, which is
# what every system built before this option carries; the seed COBY chose is
# written to coby.log, so an earlier build can be reproduced by reading the seed
# out of that file and giving it back here.
COBY_SEED=${COBY_SEED:-}
# The declash pass runs twice: once at a wide target to open the packing, and
# once at a target just above the gate to remove the few pairs that are left.
# The wide pass spreads its effort over every contact below its target, and in a
# dense box that is tens of thousands of pairs; a molecule whose pushes cancel
# is then displaced at random and the count stops falling. The tight pass sees
# only the pairs that would stop the build and moves those alone.
DECLASH_TARGET=${DECLASH_TARGET:-0.21}        # nm, the wide pass
# The distance at which a pair of bonded molecules stops the build. It is
# written here once: the gate reads it, the message that names it reads it, and
# the tight declash pass takes its target from it, so the two can never drift
# apart. A tight pass below the gate would leave the build stopping on pairs it
# had just been asked to ignore.
MIN_DIST=${MIN_DIST:-0.12}                    # nm, the gate
DECLASH_TIGHT=${DECLASH_TIGHT:-$(awk -v m="${MIN_DIST:-0.12}" 'BEGIN{printf "%.3f", m+0.01}')}
DECLASH_TIGHT_ITERS=${DECLASH_TIGHT_ITERS:-150}
# Every declash pass goes through this, so the beads it ignores and the lipids
# it is told about are written once. Two hand-written calls is how one of them
# kept a target the other had moved away from.
declash_pass(){   # declash_pass <what it is doing> <target nm> <iters> [extra args]
  local what=$1 target=$2 iters=$3; shift 3
  must "$what" \
  "The packing left two beads of different molecules on top of each other.\nRaise BOX_X and BOX_Y by 1 nm and build again, which gives the packing room.\nOr lower the lipid density with a larger COBY_APL.\nOr build again with another packing: set COBY_SEED to another integer." -- \
    "$PY" "$HELPER_DECLASH" --gro system.gro --lipids "${ALL[*]}" \
        --target "$target" --iters "$iters" --exclude-beads "$DECLASH_EXCLUDE" "$@"
}
# The sterol out-of-plane virtual sites. GROMACS rebuilds them from the itp every
# step, so a contact that involves one of them says nothing about the packing.
DECLASH_EXCLUDE=${DECLASH_EXCLUDE:-"ROH R3"}
# geometry / thermodynamics
N_COPY=${N_COPY:-4}; SPACING_NM=${SPACING_NM:-20}; MARGIN_NM=${MARGIN_NM:-8}
# WATER_NM is the water each side of the protein, in nm, and it is what sets
# box_z: add_water.py measures the z span of the protein and grows the box to
# that span plus 2*WATER_NM. 1.5 nm is three Martini water beads, which puts
# 3.0 nm between one protein end and its periodic image, above twice the 1.1 nm
# cutoff. A taller cushion costs water beads and nothing else: at 3.0 nm the
# 32 by 32 nm build of Section 3.4 carries 114,034 water beads against 90,562
# at 1.5 nm, which is 15% of the whole system.
WATER_NM=${WATER_NM:-1.5}; MEMB_THICK_NM=${MEMB_THICK_NM:-4.0}; TEMP=${TEMP:-310}
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
VEL_SEED=${VEL_SEED:-}    # gen-seed for the equilibration; set it to repeat one
                          # system with independent velocities (see run_replicate.sh)
SEED=${SEED:-0}           # seeds the random rotation of a peripheral protein ONLY.
                          # It does not reach COBY, so it does not change how the
                          # lipids are packed. Two builds that differ only in SEED
                          # and carry no peripheral protein are the same build.
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

# The build is carried out by the files below, sourced in the order they run.
# Each one is a stage of the build and reads the variables the stages before it
# set. The list is the order of the build, and it is the only place that order
# is written down.
_LIB=$(dirname "$_SELF")/lib
MEMBLE_STAGES="
  10_base.sh
  20_composition.sh
  30_orient.sh
  40_coarse_grain.sh
  50_box.sh
  60_pack.sh
  65_declash.sh
  70_solvate.sh
  80_area.sh
  85_restraints.sh
  90_formats.sh
  95_mdp.sh
  96_record.sh
  97_finish.sh
"
# MEMBLE_PLAN=1 prints what the build will do and what it will do it with, and
# writes nothing. It is the answer to "what is this command going to build", and
# it is what a methods section is written from.
if [ "${MEMBLE_PLAN:-0}" = "1" ]; then
  echo "memble $MEMBLE_VERSION would build, in this order:"
  for _s in $MEMBLE_STAGES; do
    printf '  %-20s %s\n' "$_s" "$(sed -n '1s/^# //p' "$_LIB/$_s" 2>/dev/null)"
  done
  echo ""
  echo "with:"
  for _v in PEP_AA M3_DIR BOX_X BOX_Y BOX_Z N_COPY WATER_NM SALT_M TEMP \
            UPPER LOWER LIPIDS COBY_APL COBY_SEED COBY_OPT_STEPS \
            DECLASH_TARGET DECLASH_TIGHT DECLASH_TIGHT_ITERS MIN_DIST \
            TM_RANGE TM_CORE RES_KEEP SS_MODE VEL_SEED SEED; do
    printf '  %-20s %s\n' "$_v" "${!_v:-}"
  done
  exit 0
fi

# The stages are sourced, not run as scripts: each one reads the variables the
# stages before it set. stop() is defined by the first of them, so a missing
# file here is reported without it.
for _s in $MEMBLE_STAGES; do
  if [ ! -f "$_LIB/$_s" ]; then
    echo "memble: the stage file $_s is missing from $_LIB" >&2
    echo "memble.sh carries out the build by sourcing the files under lib/, and a" >&2
    echo "clone that is missing one of them cannot build. Check what the clone holds:" >&2
    echo "  ls $_LIB" >&2
    exit 1
  fi
  . "$_LIB/$_s"
done
