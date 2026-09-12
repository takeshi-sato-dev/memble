# the box and the placement grid
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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

