# the leaflet areas and the balance pass
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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

