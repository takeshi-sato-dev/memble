# the verification gate and the run scripts
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 6. VERIFICATION GATE
#    Every post-condition that a broken system still survives is measured here.
#    memble writes run.sh only after the gate passes, so a build that did not
#    pass is never handed over as one that is ready to run.
# ====================================================================
step "measuring the eight properties of the finished system"
if [ -f "$HELPER_VERIFY" ]; then
  VERIFY=(--gro system.gro --top system.top --itp-dir "$M3_DIR" --itp-dir .
          --lipids "${ALL[*]}" --water-nm "$WATER_NM"
          --ss-mode "$SS_SOURCE" --meta memble_build.json
          --leaflet-json leaflet_area.json)
  [ -n "$SS_USED" ]  && VERIFY+=(--ss-string "$SS_USED")
  [ -n "$TM_RANGE" ] && VERIFY+=(--tm-resids "$TM_RANGE")
  [ -n "$UPPER" ]    && VERIFY+=(--expect-upper "$UPPER")
  [ -n "$LOWER" ]    && VERIFY+=(--expect-lower "$LOWER")
  if [ -z "$UPPER$LOWER" ] && [ -n "$LIPIDS" ]; then
    VERIFY+=(--expect-upper "$LIPIDS" --expect-lower "$LIPIDS")
  fi
  for c in ${MEMBLE_ALLOW:-}; do VERIFY+=(--allow "$c"); done
  if ! "$PY" "$HELPER_VERIFY" "${VERIFY[@]}"; then
    stop "the finished system did not pass verification" \
      "Every check that failed is printed above and in memble_report.txt, and each\none carries what to do about it. system.gro and system.top are kept so the\nsystem can be looked at; run.sh was not written, so nothing runs by accident.\nTo accept one named check and continue:\n  export MEMBLE_ALLOW=\"water_layer overlap\""
  fi
  echo ">>> verification passed; memble_report.txt and memble_report.json written"
fi

cat > run.sh <<RUNEOF
#!/usr/bin/env bash
set -eo pipefail
GMX="$GMX"
# One thread-MPI rank and NT OpenMP threads. A system of this size needs no
# domain decomposition, and GROMACS otherwise takes one rank per GPU it can see,
# which stops the run on a machine that carries more than one. Set NT to the
# number of cores to use:  NT=16 bash run.sh
NT=\${NT:-8}
# GROMACS dumps step<N>b.pdb when atoms move too far (system blowing up).
# set -e misses a run that "succeeds" while melting, so check explicitly.
check_blowup(){ if ls step*[0-9]b.pdb >/dev/null 2>&1; then
  echo "INSTABILITY during \$1: GROMACS wrote step*b.pdb (atoms moving too far)."
  echo "  The membrane is blowing up. Inspect: clashes (EM max force), box too small,"
  echo "  or asymmetric leaflet area mismatch. Do not continue."; exit 1; fi; }
# A stage whose coordinates are already written is left as it stands, and a
# stage that was stopped partway continues from the checkpoint GROMACS wrote,
# so a run that was interrupted is picked up by starting this script again.
# REDO=1 runs every stage from the beginning:  REDO=1 bash run.sh
stage(){   # stage <name> <mdp> <start.gro> [restraint.gro]
  out=\$1; mdp=\$2; start=\$3; restr=\${4:-}; nx=""
  [ -f index.ndx ] && nx="-n index.ndx"
  if [ -f "\$out.gro" ] && [ "\${REDO:-0}" != "1" ]; then
    echo ">>> \$out was already finished"; return 0; fi
  if [ ! -f "\$out.tpr" ] || [ "\${REDO:-0}" = "1" ]; then
    if [ -n "\$restr" ]; then
      \$GMX grompp -f "\$mdp" -c "\$start" -r "\$restr" -p system.top \$nx -o "\$out.tpr" -maxwarn 10
    else
      \$GMX grompp -f "\$mdp" -c "\$start" -p system.top \$nx -o "\$out.tpr" -maxwarn 10
    fi
  fi
  if [ -f "\$out.cpt" ] && [ "\${REDO:-0}" != "1" ]; then
    echo ">>> \$out continues from the checkpoint it left"
    \$GMX mdrun -deffnm "\$out" -v -ntmpi 1 -ntomp \$NT -cpi "\$out.cpt" -append
  else
    \$GMX mdrun -deffnm "\$out" -v -ntmpi 1 -ntomp \$NT
  fi
  check_blowup "\$out"
}
stage step6.0 step6.0_minimization.mdp system.gro
grep -i "Maximum force" step6.0.log | tail -1 || true   # should be finite, not astronomical
prev=step6.0
for k in 1 2 3 4 5 6; do
  stage step6.\${k} step6.\${k}_equilibration.mdp \${prev}.gro step6.0.gro
  prev=step6.\${k}
done
# The last equilibration stage is read before the production run starts. A
# membrane whose area is still drifting is not equilibrated, whatever the length
# of the stage was, and the production run inherits the drift.
if [ -f "$_HELPDIR/check_equilibration.py" ]; then
  if ! "$PY" "$_HELPDIR/check_equilibration.py" --edr \${prev}.edr --gmx "\$GMX" \
        --json equilibration.json; then
    if [ "\${MEMBLE_ALLOW_DRIFT:-0}" = "1" ]; then
      echo ">>> MEMBLE_ALLOW_DRIFT=1: starting the production run anyway."
    else
      echo "The production run was not started. Extend \${prev} and read it again," >&2
      echo "or set MEMBLE_ALLOW_DRIFT=1 to start regardless." >&2
      exit 1
    fi
  fi
fi
stage step7 step7_production.mdp \${prev}.gro
echo "DONE: step7.xtc"
RUNEOF
chmod +x run.sh
# also drop the staged, stop-on-failure MD runner (CHARMM-GUI style: one stage at
# a time with a success check) next to the system, copied from the helper folder.
if [ -f "$_HELPDIR/run_md.sh" ]; then
  cp "$_HELPDIR/run_md.sh" run_md.sh && chmod +x run_md.sh
fi
for _h in check_equilibration.py leaflet_area_check.py verify_system.py; do
  if [ -f "$_HELPDIR/$_h" ] && [ ! -f "./$_h" ]; then cp "$_HELPDIR/$_h" "./$_h"; fi
done
# Connectivity for viewers, CHARMM-GUI style: write a PSF with real bonds from
# the topology (ParmEd drops Martini bonds), then load system.psf and read
# system.gro or a trajectory on top of it; bonds show on every frame.
if [ -f "$_HELPDIR/write_psf.py" ]; then
  "$PY" "$_HELPDIR/write_psf.py" --gro system.gro --top system.top \
        --out system.psf --itp-dir "$M3_DIR" || true
fi
if [ -f "$_HELPDIR/write_conect_pdb.py" ]; then
  "$PY" "$_HELPDIR/write_conect_pdb.py" --gro system.gro --top system.top \
        --out system_view.pdb --itp-dir "$M3_DIR" || true
fi
cat > view.vmd <<'VMDEOF'
# CHARMM-GUI-style bonded view:
#   vmd -e view.vmd
# system.psf supplies BONDS, system_view.pdb supplies chain IDs (protein = P).
# run_md.sh writes centered, whole *_view.gro/_view.xtc (gmx trjconv -pbc mol
# -center) so the bilayer shows contiguous and mid-box; this loads those when
# present, else the raw files. Select the proteins with:  chain P   (or) protein
mol new system.psf type psf waitfor all
mol addfile system_view.pdb type pdb waitfor all
if { [file exists step7_production_view.xtc] } {
    mol addfile step7_production_view.xtc type xtc waitfor all
} elseif { [file exists step6.6_equilibration_view.gro] } {
    mol addfile step6.6_equilibration_view.gro type gro waitfor all
} elseif { [file exists step7_production.xtc] } {
    mol addfile step7_production.xtc type xtc waitfor all
}
mol delrep 0 top
mol representation VDW 0.6 12.0
mol addrep top
mol representation Bonds 0.3 12.0
mol addrep top
display resetview
VMDEOF
# and a plain-text sheet of the 8 MD stages to copy-paste one at a time
cat > md_steps.txt <<'MDEOF'
==============================================================================
 memble : run equilibration + production manually, one stage at a time
==============================================================================
 How to use:
   - Copy-paste the blocks below one at a time, from top to bottom.
   - After each block, the  ls -l ...gro  line shows whether the .gro was
     produced. If it exists, the stage succeeded -> go to the next block.
   - If the .gro is missing or you see an error, read the end of the mdrun output.
   - mdrun is run single-rank (-ntmpi 1 -ntomp 8) to avoid domain-decomposition
     errors from the protein elastic network. Change the -ntomp number to match
     your core count (find it with:  sysctl -n hw.ncpu  or  nproc).
   - A stage that was stopped partway leaves a checkpoint file, <name>.cpt.
     Continue it by adding  -cpi <name>.cpt -append  to that stage's mdrun line
     and running the line again. The grompp line is not run a second time.
     Example, for a production run that was stopped:
       gmx mdrun -deffnm step7_production -v -ntmpi 1 -ntomp 8 \
         -cpi step7_production.cpt -append
 Layout: step6.0 = minimization, step6.1..6.6 = equilibration (restraints are
         released in stages), step7 = production.
 Run everything inside the build output directory  memble_work/ .
==============================================================================


# ---- 0) go to the build directory (once) -----------------------------------
cd memble_work


# ---- 1) step6.0 : energy minimization --------------------------------------
gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr \
  -c system.gro -r system.gro -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.0_minimization -v -ntmpi 1 -ntomp 8
ls -l step6.0_minimization.gro
# expect: "Potential Energy" negative, "Maximum force" finite (not inf)


# ---- 2) step6.1 : equilibration 1 ------------------------------------------
gmx grompp -f step6.1_equilibration.mdp -o step6.1_equilibration.tpr \
  -c step6.0_minimization.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.1_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.1_equilibration.gro


# ---- 3) step6.2 : equilibration 2 ------------------------------------------
gmx grompp -f step6.2_equilibration.mdp -o step6.2_equilibration.tpr \
  -c step6.1_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.2_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.2_equilibration.gro


# ---- 4) step6.3 : equilibration 3 ------------------------------------------
gmx grompp -f step6.3_equilibration.mdp -o step6.3_equilibration.tpr \
  -c step6.2_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.3_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.3_equilibration.gro


# ---- 5) step6.4 : equilibration 4 ------------------------------------------
gmx grompp -f step6.4_equilibration.mdp -o step6.4_equilibration.tpr \
  -c step6.3_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.4_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.4_equilibration.gro


# ---- 6) step6.5 : equilibration 5 ------------------------------------------
gmx grompp -f step6.5_equilibration.mdp -o step6.5_equilibration.tpr \
  -c step6.4_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.5_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.5_equilibration.gro


# ---- 7) step6.6 : equilibration 6 (final) ----------------------------------
gmx grompp -f step6.6_equilibration.mdp -o step6.6_equilibration.tpr \
  -c step6.5_equilibration.gro -r step6.0_minimization.gro \
  -p system.top -n index.ndx -maxwarn 20
gmx mdrun -deffnm step6.6_equilibration -v -ntmpi 1 -ntomp 8
ls -l step6.6_equilibration.gro
# equilibration is complete once all of these .gro files exist


# ---- 8) step7 : production (long, ~us scale; run in the background) ---------
gmx grompp -f step7_production.mdp -o step7_production.tpr \
  -c step6.6_equilibration.gro -r step6.6_equilibration.gro \
  -p system.top -n index.ndx -maxwarn 20
nohup gmx mdrun -deffnm step7_production -v -ntmpi 1 -ntomp 8 > log_step7.txt 2>&1 &
sleep 20; tail -15 log_step7.txt
# watch progress:  tail -f log_step7.txt
# stop it:         pkill -f step7_production


==============================================================================
 Troubleshooting
==============================================================================
 - "no domain decomposition ..."  -> -ntmpi 1 is not in effect. Re-paste the
   command exactly (it must contain  -ntmpi 1 -ntomp 8 ).
 - mdrun stops with inf           -> problem in system.gro; rebuild from scratch.
 - a few LINCS warnings           -> normal. Many warnings with no .gro produced
   -> read the end of that stage's output.
 - threads: change the -ntomp number to your core count
   (sysctl -n hw.ncpu on macOS, nproc on Linux).
==============================================================================
MDEOF
# One last line, so a build that ran in the background is read from its tail.
_NPASS=$(awk '/^RESULT:/{print $2}' memble_report.txt 2>/dev/null)
echo ""
echo "MEMBLE ${_NPASS:-DONE}: ${OUTTAG} built in $(elapsed) s, all eight properties measured."
echo ">>> Build complete in $WORK."
echo ">>> Equilibrate + produce stage-by-stage:  cd $WORK && bash run_md.sh"
echo ">>> Or copy-paste stages manually from:    $WORK/md_steps.txt"
echo ">>> View WITH BONDS (CHARMM-GUI style):     cd $WORK && vmd -e view.vmd"
echo ">>>   (loads system.psf for connectivity, then the trajectory/gro)"
echo ">>> (legacy all-in-one script:             ./run.sh)"

