# the build record
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

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
  printf '  "coby_seed": "%s",\n' "$(awk -F': *' '/Setting random seed to/{print $2; exit}' coby.log 2>/dev/null)"
  printf '  "gromacs": "%s",\n' "$("$GMX" --version 2>/dev/null | awk -F': *' '/GROMACS version/{print $2; exit}')"
  printf '  "martinize2": "%s",\n' "$("$MARTINIZE2" --version 2>&1 | head -1 | tr -d '"')"
  printf '  "coby": "%s",\n' "$("$PY" -c 'import COBY;print(getattr(COBY,"__version__","unknown"))' 2>/dev/null)"
  printf '  "balance_passes": "%s",\n' "${_BAL_LOG:-}"
  printf '  "degraded_steps": "%s"\n' "${MEMBLE_DEGRADED:-}"
  printf '}\n'
} > memble_build.json
echo ">>> build record written to memble_build.json"

