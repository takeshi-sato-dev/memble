# the mdp files and index.ndx
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 8. mdp (GROMACS 2023.x; TEMP-parameterized) + run.sh
# ====================================================================
step "writing the staged equilibration"
RF='nstlist = 20
cutoff-scheme = Verlet
verlet-buffer-tolerance = 0.005
coulombtype = reaction-field
coulomb-modifier = Potential-shift
rcoulomb = 1.1
epsilon_r = 15
epsilon_rf = 0
vdwtype = cutoff
vdw-modifier = Potential-shift-verlet
rvdw = 1.1'
{ echo "integrator = steep"; echo "nsteps = 50000"; echo "emtol = 100.0"; echo "emstep = 0.001"; echo "define = -DFLEXIBLE"; echo "$RF"; } > step6.0_minimization.mdp
# step6.1 is a SECOND minimization with the protein position-restrained (no md
# yet), the CHARMM-GUI Martini Maker approach: it lets the restrained system
# relax from the rough packed start before any dynamics, so the protein cannot
# be kicked out of the membrane by a finite-dt integration of a stiff structure.
{ echo "integrator = steep"; echo "nsteps = 50000"; echo "emtol = 200.0"; echo "emstep = 0.001"
  echo "define = -DPOSRES_STEP1"; echo "refcoord-scaling = all"; echo "$RF"
} > step6.1_equilibration.mdp
# md equilibration stages 6.2..6.6: ramp the timestep, keep the protein
# restrained (weakening per stage via the POSRES_STEP blocks), refcoord-scaling
# = all so the restraint reference follows the box under semiisotropic pressure.
for k in 2 3 4 5 6; do i=$((k-1)); nst=$(awk -v ns="${STAGE_NS[$i]}" -v dt="${STAGE_DT[$i]}" 'BEGIN{printf "%d",(ns/dt)*1000}')
  { echo "integrator = md"; echo "dt = ${STAGE_DT[$i]}"; echo "nsteps = $nst"; echo "$RF"
    echo "tcoupl = v-rescale"; echo "tc-grps = SOLU_MEMB SOLV"; echo "tau-t = 1.0 1.0"; echo "ref-t = $TEMP $TEMP"
    echo "pcoupl = c-rescale"; echo "pcoupltype = semiisotropic"; echo "tau-p = 4.0"
    echo "compressibility = 3e-4 3e-4"; echo "ref-p = 1.0 1.0"; echo "refcoord-scaling = all"
    echo "gen-vel = yes"; echo "gen-temp = $TEMP"; echo "constraints = none"; echo "define = -DPOSRES_STEP$k"
    # gen-seed is left at the GROMACS default (-1, drawn from the process) unless
    # VEL_SEED is set. A repeat of one system is made by running the same built
    # and minimized structure again from a different VEL_SEED, which is the only
    # thing that separates two runs of one composition: SEED seeds the rotation
    # of a peripheral protein and does not reach the packing.
    [ -n "${VEL_SEED:-}" ] && echo "gen-seed = $VEL_SEED"
  } > step6.${k}_equilibration.mdp
done
{ echo "integrator = md"; echo "dt = 0.02"; echo "nsteps = $NPROD_STEPS"; echo "$RF"
  echo "tcoupl = v-rescale"; echo "tc-grps = SOLU_MEMB SOLV"; echo "tau-t = 1.0 1.0"; echo "ref-t = $TEMP $TEMP"
  echo "pcoupl = parrinello-rahman"; echo "pcoupltype = semiisotropic"; echo "tau-p = 12.0"
  echo "compressibility = 3e-4 3e-4"; echo "ref-p = 1.0 1.0"
  echo "nstxout-compressed = 5000"; echo "compressed-x-precision = 1000"
  echo "nstlog = 5000"; echo "nstenergy = 5000"
} > step7_production.mdp
LIPRESN="${ALL[*]}"
# index.ndx generated here (not via interactive gmx select in run.sh): groups
# SOLV / MEMB / SOLU / SOLU_MEMB by residue name, matching system.gro atom order.
"$PY" - "$LIPRESN" <<'NDXEOF'
import sys
lip = set(sys.argv[1].split())
solv = {"W", "WF", "NA", "CL", "ION", "NA+", "CL-"}
L = open("system.gro").read().splitlines(); n = int(L[1]); body = L[2:2+n]
g = {"SOLV": [], "MEMB": [], "SOLU": [], "SOLU_MEMB": []}
for i, ln in enumerate(body, 1):
    rn = ln[5:10].strip()
    if rn in solv:
        g["SOLV"].append(i)
    elif rn in lip:
        g["MEMB"].append(i); g["SOLU_MEMB"].append(i)
    else:
        g["SOLU"].append(i); g["SOLU_MEMB"].append(i)
with open("index.ndx", "w") as fh:
    for name in ("SOLV", "MEMB", "SOLU", "SOLU_MEMB"):
        fh.write("[ %s ]\n" % name)
        idx = g[name]
        for k in range(0, len(idx), 15):
            fh.write(" ".join("%d" % x for x in idx[k:k+15]) + "\n")
        fh.write("\n")
print(">>> index.ndx: SOLV=%d MEMB=%d SOLU=%d SOLU_MEMB=%d"
      % (len(g["SOLV"]), len(g["MEMB"]), len(g["SOLU"]), len(g["SOLU_MEMB"])))
NDXEOF
# peripheral protein: add a PARTNER index group and a one-sided flat-bottom pull
# (MEMB vs PARTNER COM along z) to production, so it can associate with its
# leaflet but can never wrap to the other leaflet through PBC.
if [ -n "$PARTNER_NAME" ] && [ -n "$HELPER_PARTPULL" ]; then
  "$PY" "$HELPER_PARTPULL" --gro system.gro --ndx index.ndx \
      --mdp step7_production.mdp --partner-name "$PARTNER_NAME" \
      --lipids "${ALL[*]}" --margin "$PARTNER_MARGIN" --k "$PARTNER_K"
fi
