# sizing the box, water and ions
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 4b. ADD WATER: COBY builds the membrane in a thin box (a large box_z hangs its
#     lipid-grid optimizer). Now grow box_z so a protruding TM-JM protein gets a
#     real bulk-water cushion, and solvate the new slabs (+ salt). No-op if the
#     protein already fits with the requested water on each side.
# ====================================================================
step "sizing the box and adding water and ions"
if [ -n "$PARTNER_CG" ]; then
  # peripheral protein: place on the chosen leaflet, expand box, then fill water
  "$PY" "$HELPER_PART" --system-gro system.gro --system-top system.top \
      --partner-cg "$PARTNER_CG" --partner-name "$PARTNER_NAME" \
      --partner-itp "$PARTNER_ITP" --side "$PARTNER_SIDE" --gap "$PARTNER_GAP" \
      --water-nm "$PARTNER_WATER" --lipids "${ALL[*]}" \
      --rotate "$PARTNER_ROTATE" --seed "$SEED"
  if [ -n "$HELPER_ADDWATER" ]; then
    must "adding water and ions around the peripheral protein" \
  "The peripheral protein was placed, and the solvation of the new box failed.\nCheck PARTNER_WATER (nm, for example 2.5) and SALT_M (molar, for example 0.15).\nCheck that the Martini water itp is in M3_DIR:  ls $M3_DIR | grep -i water" -- \
      "$PY" "$HELPER_ADDWATER" --gro system.gro --top system.top \
        --water-nm "$PARTNER_WATER" --salt "$SALT_M" --keep-box
  fi
elif [ -n "$HELPER_ADDWATER" ]; then
  # Make the protein contiguous BEFORE sizing the box. add_water derives
  # box_z from the protein z-extent it can see; while a chain is still split
  # across the z boundary that extent is the wrapped one, which for a tall
  # TM-JM assembly is far smaller than the real span (6.95 vs 13.73 nm in the
  # four-copy EGFR build), so the box came out ~7 nm too short and the protein
  # ended up overlapping its own periodic image. Unwrap first, then size.
  # The call after the solvation step below stays: it re-centers in the final
  # box and reports the resulting water cushion.
  if [ -n "$HELPER_WHOLE" ] && [ -f "$HELPER_WHOLE" ]; then
    must "making the protein contiguous before the box is sized" \
  "add_water measures the z span of the protein to size the box. While a chain is\nstill split across the z boundary that span is the wrapped one, which is far\nshorter than the protein, and the box comes out too short.\nCheck that the itp files of the protein are in the working directory:  ls molecule_*.itp\nIf the protein has one chain and does not cross the boundary, this step can be\nskipped with HELPER_WHOLE= (empty)." -- \
      "$PY" "$HELPER_WHOLE" --gro system.gro --top system.top --itp-dir .
  fi
  must "sizing the box and adding water and ions" \
  "Check WATER_NM (nm per side, for example 2.5) and SALT_M (molar, for example 0.15).\nCheck that the Martini water and ion itp files are in M3_DIR.\nA very tall protein in a small xy box can leave no room for water; raise BOX_X\nand BOX_Y, or lower WATER_NM." -- \
    "$PY" "$HELPER_ADDWATER" --gro system.gro --top system.top \
      --water-nm "$WATER_NM" --salt "$SALT_M"
fi

# final PBC-aware declash: catches any bead sitting just outside the box that
# clashes with the opposite face under periodic boundaries (a common COBY edge
# effect that produces an infinite force on a water at minimization).
# Make the protein contiguous across PBC and center the system on it, BEFORE the
# final declash. A tall TM-JM protein can straddle the z boundary after
# assembly/solvation; per-atom wrapping upstream then splits a chain across the
# box (consecutive backbone beads a full box apart -> LINCS blowup -> infinite
# force). Recentering can push some lipids/water across the edge, so the final
# declash must run AFTER this, as the last coordinate step, to clean up.
if [ -n "$HELPER_WHOLE" ] && [ -f "$HELPER_WHOLE" ]; then
  must "making the protein contiguous in the final box" \
  "Recentering can split a chain across the box edge, and two consecutive backbone\nbeads a full box apart give an infinite force at minimization.\nCheck that the protein itp files are in the working directory:  ls molecule_*.itp" -- \
    "$PY" "$HELPER_WHOLE" --gro system.gro --top system.top --itp-dir . --box-is-final
fi

# FINAL declash and the last coordinate-modifying step: guarantees no two beads
# from different molecules overlap. The protein is frozen (never moved), so it
# is not distorted or re-split; only solvent and lipids are pushed apart.
declash_pass "the final separation of overlapping beads" "$DECLASH_TARGET" 200 --freeze-protein

# The tight pass. A 32 by 32 nm box holds seven times the molecules of a 12 by
# 12 nm box, and the wide pass above leaves several thousand contacts below its
# own target in a box that size. Those contacts are almost all water, and the
# handful that hold two bonded molecules are the ones that stop the build. This
# pass takes a target just above the gate, so it sees those few pairs and moves
# only the molecules that carry them.
declash_pass "separating the last pairs that would stop the build" "$DECLASH_TIGHT" "$DECLASH_TIGHT_ITERS" --freeze-protein

# Confirm there is no residual overlap that would give an infinite force. A
# warning here was read past and the build shipped, so this now stops the build.
# MEMBLE_ALLOW_OVERLAP=1 keeps the old behavior for a system that is being
# inspected rather than run.
if [ -n "$HELPER_MINDIST" ] && [ -f "$HELPER_MINDIST" ]; then
  if ! "$PY" "$HELPER_MINDIST" --gro system.gro --lipids "${ALL[*]}" \
          --min "$MIN_DIST" --exclude-beads "$DECLASH_EXCLUDE"; then
    if [ "${MEMBLE_ALLOW_OVERLAP:-0}" = "1" ]; then
      echo ">>> MEMBLE_ALLOW_OVERLAP=1: continuing with a residual overlap"
      echo ">>> between two bonded molecules. Watch stage 6.3: a molecule that is"
      echo ">>> still overlapped when the restraints are eased is torn apart, and"
      echo ">>> LINCS then reports a constraint deviation of millions."
    else
      stop "two bonded molecules are closer than $MIN_DIST nm" \
        "Minimization does not always separate such a pair. A molecule that is still\noverlapped when the restraints of stage 6.3 are eased is torn apart, and the run\nthen makes no progress. The pair is named just above. A pair holding a water\nbead or an ion does not stop a build: both are single free particles and the\nminimization moves them apart in its first steps.\nSEED does not change the packing: it seeds the rotation of a peripheral protein\nand nothing else, so building again with another SEED returns the same pair.\n  1. build again with another packing:  set COBY_SEED to another integer\n  2. give the packing room:             raise BOX_X and BOX_Y by 1 nm\n  3. lower the lipid density:           raise COBY_APL\n  4. let the packing optimizer work:    raise COBY_OPT_STEPS\n  5. push the last pairs harder:        raise DECLASH_TIGHT_ITERS\n  6. if the pair holds the protein:     raise SPACING_NM or lower N_COPY\n  7. to keep this system and look at it: export MEMBLE_ALLOW_OVERLAP=1"
    fi
  fi
fi

# Restore the protein's original per-chain residue numbers (martinize renumbers
# every chain from 1, so the assembled chains overlap). The true numbers come
# from oriented_aa.pdb (which preserves the input PDB resSeq); nothing hardcoded.
if [ -f oriented_aa.pdb ] && [ -n "$HELPER_FIXRESID" ]; then
  must "restoring the per-chain residue numbers of the input PDB" \
  "martinize2 renumbers every chain from 1, so the chains of a multi-chain protein\noverlap. The true numbers come from oriented_aa.pdb.\nCheck that oriented_aa.pdb exists and holds the input residue numbers.\nThe system is correct without this step; only the numbering in system.gro is\nmartinize2 numbering. Set HELPER_FIXRESID= (empty) to accept that." -- \
    "$PY" "$HELPER_FIXRESID" --gro system.gro --oriented oriented_aa.pdb \
        --top system.top --itp-dir .
fi

