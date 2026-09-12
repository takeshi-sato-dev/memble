# the staged position restraints
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 6. staged restraints: protein BB + every lipid head
# ====================================================================
for s in "${UNIQ_SRC[@]}"; do _l=$(local_for "$s"); cp "$s" "$_l"; chmod u+w "$_l"; done
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

