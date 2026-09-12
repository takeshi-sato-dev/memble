# orienting the protein on its transmembrane range
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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

