# COBY: packing the lipids around the protein
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 4. COBY build (membrane string from composition / leaflets)
# ====================================================================
step "packing the lipids around the protein"
export COBY_SEED
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
    # COBY seeds itself from the wall clock unless it is given a seed
    **({"randseed": int(os.environ["COBY_SEED"])} if os.environ.get("COBY_SEED") else {}),
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

