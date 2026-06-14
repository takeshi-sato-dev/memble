#!/usr/bin/env bash
#
# build_ensemble.sh <all_atom_protein.pdb>
#
# Build N_REP independent systems, each with the peripheral partner placed in a
# different random orientation (rotation seed = replicate index). This is the
# orientation-randomized encounter ensemble for unbiased membrane-association
# sampling. All other parameters (LIPIDS, PARTNER, M3_DIR, GMX, helpers, box,
# water, temperature, salt) are inherited from the environment.
#
set -euo pipefail
PEP=${1:?usage: build_ensemble.sh <all_atom_protein.pdb>}
: "${N_REP:?set N_REP to the number of replicates}"
here=$(cd "$(dirname "$0")" && pwd)
for i in $(seq 1 "$N_REP"); do
  echo ">>> replicate $i / $N_REP (partner orientation seed=$i)"
  SEED="$i" OUTTAG="rep_$i" "$here/memble.sh" "$PEP"
done
echo "ensemble built: rep_1 .. rep_${N_REP} (run each rep_i/run.sh)"
