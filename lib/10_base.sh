# shell utilities: paths, step, must, stop
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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
