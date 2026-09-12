# the composition and the topology of every lipid
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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

