# separating the beads the packing left overlapping
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 4a. DECLASH: push apart inter-molecular bead overlaps from dense packing
#     (rigid-molecule moves; intramolecular geometry preserved) so that
#     minimization does not hit infinite Lennard-Jones forces.
# ====================================================================
# COBY places coarse-grained lipids by their real beads but leaves out-of-plane
# virtual sites (sterol ROH/R3, funct 4) flat; GROMACS then reconstructs them
# ~0.1 nm away and they can explode at minimization. Rebuild every vsite exactly
# from its itp definition before declashing/minimizing.
must "rebuilding the sterol virtual sites from their itp definitions" \
  "COBY leaves out-of-plane virtual sites (sterol ROH and R3, funct 4) flat, and\nGROMACS rebuilds them about 0.1 nm away, which explodes at minimization.\nCheck that M3_DIR points at the Martini 3 lipidome directory and that it holds\nthe itp of every sterol in the composition:  ls $M3_DIR/*.itp | grep -i chol\nA composition without a sterol does not need this step; remove the sterol or\nadd its itp to M3_DIR." -- \
  "$PY" "$HELPER_FIXVS" --gro system.gro --top system.top --itp-dir "$M3_DIR"
declash_pass "separating overlapping beads left by the packing" "$DECLASH_TARGET" 200

