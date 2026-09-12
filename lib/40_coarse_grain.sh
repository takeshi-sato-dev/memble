# secondary structure and martinize2
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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

