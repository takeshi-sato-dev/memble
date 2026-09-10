# memble

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20684812.svg)](https://doi.org/10.5281/zenodo.20684812)

A composition-general builder for Martini 3 membrane-protein systems, producing
CHARMM-GUI-equivalent inputs (GROMACS topology, CHARMM PSF and CRD, a staged
equilibration protocol, and an index file) without CHARMM-GUI. The tool takes an
all-atom protein structure and a lipid composition and returns a system that is
ready to minimize, equilibrate, and run in GROMACS 2023.

## What it does

Given one all-atom protein PDB, the pipeline:

1. orients the protein so its transmembrane helix lies along z (automatic
   detection, with an optional residue-range override),
2. coarse-grains it with martinize2 using secondary structure from DSSP and no
   global elastic network, so juxtamembrane segments stay flexible,
3. places a chosen number of copies on a lateral grid,
4. builds the bilayer, inserts the protein copies, solvates, and adds ions with
   COBY (Martini 3 native, so cholesterol and other sterols are placed with the
   correct Martini 3 mapping by construction),
5. checks leaflet area balance before any molecular dynamics is run,
6. injects a staged set of position restraints on the protein backbone and on
   every lipid head bead,
7. writes the CHARMM PSF, CRD, and PDB with ParmEd,
8. writes a minimization, a six-stage equilibration, and a production parameter
   set for GROMACS 2023, plus a single run script.

Lipid composition, leaflet asymmetry, lateral box size, water-layer thickness,
and temperature are all parameters. Nothing protein-specific or
composition-specific is hardcoded.

## Requirements

- Python 3.9 or newer with `numpy`, `parmed`, and `COBY`
- `martinize2` (Vermouth) and `dssp` (or `mkdssp`) on the PATH
- GROMACS 2023.x
- The Martini 3 lipidome: the itp files that define the lipids, the water and
  the ions

`setup.sh` below installs the Python tools. DSSP and GROMACS are installed
separately through your package manager or module system, and the lipidome is
downloaded once, as follows.

### Getting the Martini 3 lipidome

memble places a lipid only if the supplied lipidome defines it, so the itp files
have to be on disk before the first build. The parameters are the refined
Martini 3 lipidome of Pedersen et al., ACS Cent. Sci. 2025, 11, 1598, and the
Martini Force Field Initiative distributes them:

```
git clone https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters
```

Collect every itp file into one directory, including `martini_v3.0.0.itp`, the
ions and the solvents, and give that directory as `M3_DIR`. A build reads the
directory where it stands, and copies or installs nothing:

```
M3_DIR=~/m3lipidome
ls $M3_DIR
  martini_v3.0.0.itp                        martini_v3.0.0_phospholipids_PC_v2.itp
  martini_v3.0.0_ions_v1.itp                martini_v3.0.0_sterols_v1.itp
  martini_v3.0.0_solvents_v1.itp            martini_v3.0_phosphoinositides_v1.0.itp
  ...
```

## Installation

Two steps, in this order. A build started before `setup.sh` has run stops with
`command not found`, and that message names neither the tool that is missing nor
the step that was skipped.

```
git clone https://github.com/takeshi-sato-dev/memble.git
cd memble
bash setup.sh
```

`setup.sh` creates a virtual environment at `~/memble-venv` (override with
`MEMBLE_VENV`), installs numpy, parmed, COBY, vermouth, mdtraj and streamlit
into it, and then reports every import that failed, where martinize2 is, and
whether GROMACS is on the PATH. Nothing in an existing Python setup changes.
Rerunning it deletes the environment and builds it again, so a botched install
is repaired by running it a second time, and `rm -rf ~/memble-venv` resets it.

Then build, without activating anything:

```
M3_DIR=~/m3lipidome GMX=gmx LIPIDS="CHOL:1 DIPC:1 DPSM:1" \
  ./memble protein_AA.pdb        # one build, or an ensemble if N_REP > 1

./memble-gui                     # the same build, from a browser
```

`memble` and `memble-gui` call the interpreter and the tools of that virtual
environment and set the helper paths, so `source .../activate` is never needed.
The helper scripts themselves (`orient_tm.py`, `replicate_and_fix_top.py`,
`inject_posres.py`, `leaflet_area_check.py` and the rest) live at the repository
root, and their paths are given by hand only when `memble.sh` is called directly,
as in the last example below.

## Quick start

Symmetric ternary raft mixture, four protein copies:

```
LIPIDS="CHOL:1 DIPC:1 DPSM:1" \
M3_DIR=~/m3lipidome \
GMX=gmx \
./memble protein_AA.pdb

cd memble_work
NT=8 bash run.sh
```

Asymmetric bilayer, explicit box, thicker water, lower temperature:

```
UPPER="CHOL:1 DPSM:2 POPC:1" LOWER="POPC:3 POPE:1" \
BOX_X=60 BOX_Y=50 WATER_NM=3 TEMP=300 \
M3_DIR=~/m3lipidome GMX=gmx \
./memble protein_AA.pdb
```

`memble.sh` can be called in place of `memble`, and it then reads the interpreter
and the helper paths from the command line:

```
PY=~/memble-venv/bin/python MARTINIZE2=~/memble-venv/bin/martinize2 \
HELPER_REP=$PWD/replicate_and_fix_top.py HELPER_POS=$PWD/inject_posres.py \
HELPER_ORI=$PWD/orient_tm.py HELPER_AREA=$PWD/leaflet_area_check.py \
LIPIDS="CHOL:1 DIPC:1 DPSM:1" M3_DIR=~/m3lipidome GMX=gmx \
bash memble.sh protein_AA.pdb
```

If the two leaflets are not area-balanced, the build stops before any molecular
dynamics is run and reports the per-leaflet areas, so the composition or box can
be adjusted and the system rebuilt.

Before it stops, memble tries to correct the counts itself. It measures the areas,
corrects the area per lipid of each leaflet, and builds again. The pass runs only
on a difference that exceeds both `MEMBLE_BALANCE_TOL` and twice its own standard
error, so a difference the size of the measurement scatter is left alone. Each
pass moves the areas per lipid `MEMBLE_BALANCE_DAMP` of the way to the
measurement, and memble restores the earlier build when a pass does not lower the
difference. `memble_build.json` records the difference of every pass in
`balance_passes`.

## Parameters

| Variable | Default | Meaning |
| --- | --- | --- |
| `LIPIDS` | `CHOL:1 DIPC:1 DPSM:1` | symmetric composition, `name:ratio[:head]` |
| `UPPER`, `LOWER` | unset | per-leaflet compositions; set both for asymmetry |
| `COBY_MEMBRANE` | unset | full override of the COBY membrane string |
| `N_COPY` | 4 | number of protein copies |
| `SPACING_NM` | 20 | lateral spacing used for automatic box sizing |
| `BOX_X`, `BOX_Y` | auto | explicit lateral box dimensions in nm |
| `WATER_NM` | 2.5 | water layer per side in nm; sets the box height |
| `MEMB_THICK_NM` | 4.0 | approximate bilayer thickness for box height |
| `TEMP` | 310 | temperature in kelvin |
| `SALT_M` | 0.15 | NaCl molarity |
| `TM_RANGE` | auto | residue range for orientation, for example `619:641` |
| `SS_OVERRIDE` | unset | explicit secondary-structure string (highest priority) |
| `SS_MODE` | `dssp` | `dssp` runs DSSP in memble and hands the resulting string to martinize2 with `-ss`, so the assignment is recorded in the output; `dssp-internal` lets martinize2 run DSSP itself, and the string is then not recorded; `tm` sets the TM range to helix (GPCRs, multi-helix proteins); `tm` sets the TM range to helix and everything else to coil, so the juxtamembrane stays flexible (single-pass TM and TM-JM peptides) |
| `DSSP` | auto | path to a DSSP 2.2.1/3.x binary, or `mdtraj`; 4.x is not compatible with martinize2 |
| `TM_CORE` | unset | per-chain TM core, e.g. `A:65-88;B:65-88`; used by `SS_MODE=tm` and by prebuilt orientation |
| `PREBUILT_MULTI` | 0 | `1` uses the chains already assembled in the input PDB instead of replicating one chain |
| `MULTI_TM` | 0 | `1` orients a multi-pass TM bundle (e.g. a 7-TM GPCR) on its helix bundle instead of one principal axis: every membrane-length helix is detected from DSSP and the membrane normal is set to the mean helix axis |
| `MULTI_TM_MINLEN` | 12 | minimum helix length (residues) counted as a TM helix in `MULTI_TM` orientation |
| `NTERM_SIDE` | unset | `up` or `down`: force the N-terminus to face +z (upper leaflet) or -z (lower leaflet). The orientation axis sign is otherwise arbitrary, so a multi-pass receptor can come out inverted. For a GPCR the N-terminus is extracellular, so `up` puts the extracellular side on the upper leaflet |
| `RES_KEEP` | unset | per-chain residues to keep with `PREBUILT_MULTI`, e.g. `A:54-103;B:54-103` |
| `AREA_TOL` | 0.08 | largest leaflet area difference accepted |
| `AREA_HARD_TOL` | 0.25 | an asymmetric build stops above this difference |
| `MEMBLE_BALANCE_ITER` | 2 | extra builds allowed for the leaflet balance pass; `0` keeps the first build |
| `MEMBLE_BALANCE_TOL` | `AREA_TOL` | leaflet difference above which the balance pass runs |
| `MEMBLE_BALANCE_DAMP` | 0.5 | share of the correction that one balance pass applies |
| `APL_UPPER`, `APL_LOWER` | unset | area per lipid handed to COBY for each leaflet, bypassing the composition estimate |
| `MEMBLE_ALLOW` | unset | names of checks to accept, e.g. `"water_layer overlap"` |
| `MEMBLE_ALLOW_OVERLAP` | 0 | `1` keeps a system that still holds a bead overlap |
| `MEMBLE_ALLOW_DRIFT` | 0 | `1` starts the production run while the membrane area is still drifting |
| `MEMBLE_IGNORE_ERRORS` | 0 | `1` continues past a failed build step and records it in `memble_build.json` |
| `MEMBLE_NO_TM_CENTER` | 0 | `1` leaves COBY to center the protein on the whole molecule instead of the TM range |
| `NPROD_STEPS` | 4e8 | production steps at 20 fs (8 microseconds) |

### Choosing the secondary structure source

DSSP reads secondary structure from the input coordinates. If the input has a
helical juxtamembrane, DSSP marks it helix, and the backbone dihedral restraints
then hold that region rigid. For a single-pass TM or TM-JM peptide where the
juxtamembrane should stay flexible, set `SS_MODE=tm`: the TM range becomes helix
and everything else coil, independent of the input conformation and of DSSP. For
multi-helix proteins such as GPCRs, keep `SS_MODE=dssp` so the loops, the
membrane helices, and any peripheral helix are each assigned correctly.

### Orienting a multi-pass TM bundle (GPCRs)

A single-pass TM peptide has one membrane-spanning helix, so the long axis of
that helix is the membrane normal. A multi-pass receptor (for example a 7-TM
GPCR) is one chain folded into a bundle of helices, and the first principal
component of the pooled coordinates points along the widest in-plane spread, not
along the normal. Orienting on that axis lays the receptor on its side.

Set `MULTI_TM=1` for these proteins. Every membrane-length helix is detected
(from `SS_OVERRIDE` if given, otherwise from DSSP on the input), each helix axis
is flipped into a common hemisphere, and the mean of those axes is used as the
membrane normal. `PREBUILT_MULTI` and `MULTI_TM` address different things:
`PREBUILT_MULTI` keeps several chains already assembled in the input, while
`MULTI_TM` orients one chain that crosses the membrane several times. A GPCR is
a single chain, so it uses `MULTI_TM=1` without `PREBUILT_MULTI`. Note that
input structures with non-standard residues (bound ligands, lipids, sugars,
ions) must have those removed before martinize2, which reads standard residues
only.

martinize2 cannot parse DSSP 4.x. Install a 3.x build, for example

```
micromamba create -n dssp3 -c bioconda 'dssp=3.1.4'
```

and pass `DSSP=$(micromamba run -n dssp3 which mkdssp)`, or avoid DSSP entirely
with `SS_MODE=tm`.

## Outputs

In the working directory: `system.gro`, `system.pdb`, `system.top` with its itp
includes, `system.psf` and `system_vmd.psf`, `system.crd`, the minimization,
six equilibration, and production parameter files, and `run.sh`.

Four files record what the system is and how it was built:

| File | Holds |
| --- | --- |
| `memble_report.txt`, `memble_report.json` | the eight properties measured on the finished system, and what to do about any that failed |
| `memble_build.json` | the checksum of the input structure, the secondary structure string, the composition, the box, the temperature, the salt concentration, and the versions of COBY, martinize2 and GROMACS, and the leaflet difference of every balance pass |
| `leaflet_area.json` | the measured area of every lipid in each leaflet, the area the protein occupies, and how much of each leaflet is made of lipids that both leaflets hold |
| `equilibration.json` | the area drift and the surface tension over the second half of the last equilibration stage |

## What memble checks before it writes the run files

`run.sh` is written only when all eight of these hold. A check that fails prints
what failed and what to do about it, and `MEMBLE_ALLOW` accepts one by name.

1. The secondary structure string covers every residue of the protein.
2. The lipids placed match the composition requested.
3. The membrane-spanning range sits at the bilayer midplane.
4. The system is neutral.
5. The protein and its image in z are further apart than the non-bonded cutoff.
6. The water above and below the bilayer reaches the requested depth.
7. No two beads of different molecules are closer than 0.12 nm.
8. A lipid present in both leaflets occupies the same area in both.

The leaflet areas are measured, not looked up. Each leaflet is tessellated in
the membrane plane, with every lipid and every protein bead inside the leaflet
given the region closer to it than to anything else, and the regions are
integrated by Monte Carlo under the periodic boundary conditions. memble builds
the system once, measures, corrects the area per lipid of each leaflet from the
measurement, and builds again.

## How the equilibration avoids early instability

The six-stage equilibration releases the protein-backbone and lipid-head
restraints gradually while the time step rises from 2 to 20 fs, the pressure is
controlled with the c-rescale barostat until the bilayer settles, and the
Parrinello-Rahman barostat is used only for production. The run script checks
for the step coordinate files that GROMACS writes when atoms move too far and
stops if the system is destabilizing.

### Running on a machine that carries more than one GPU

GROMACS takes one thread-MPI rank for every GPU it can see, and it stops the run
when the number of threads is not divisible by the number of ranks. A system of
this size needs no domain decomposition, so `run.sh`, `run_md.sh` and
`md_steps.txt` all place one rank and give it `NT` OpenMP threads. `NT` defaults
to 8, and the number of cores to use is set on the command line:

```
NT=16 bash run.sh
```

To give one run one card, name the card as well:

```
CUDA_VISIBLE_DEVICES=0 NT=16 bash run.sh
```

## Tests

```
python3 -m pytest tests/ -v
```

The suite exercises each helper on synthetic inputs: orientation tilt, grid
placement and separation, automatic and explicit head-bead restraint injection,
the leaflet area check for both a balanced and a mismatched bilayer, the eight
properties of the verification gate on a passing and a failing system, and the
balance pass on a sequence of measurements it improves and one it does not.

## Limitations

- A lipid can only be placed if it exists in the supplied Martini 3 lipidome.
- The leaflet areas are measured on the system as built, before minimization,
  and a packed bilayer is not an equilibrated one. The check is a build-time
  guard, not a measurement of the equilibrated membrane.
- The bilayer is treated as flat and the areas are projected on the xy plane.
- The area of a region is set by the neighbors of the molecule, so comparing a
  lipid against itself in the other leaflet reports the state of the two
  leaflets only while both present a similar neighborhood. Where the leaflets
  hold less than half of their lipids in common, the difference is reported and
  does not stop the build.
- Multi-headgroup lipids (cardiolipin, phosphoinositides, gangliosides) may need
  an explicit head bead given in the composition string.
- The COBY membrane-string grammar (lipid ratio token, leaflet token) can depend
  on the COBY version. Confirm it with `COBY -h`, set `COBY_LEAFLET` if needed,
  or override the whole string with `COBY_MEMBRANE`.

## License

Released under the Apache License 2.0. See `LICENSE`.

## Contributing

Issues and pull requests are welcome.

## Citation

If this tool is useful in your work, please cite the accompanying paper and the
underlying tools: Martini 3, the Martini 3 lipidome, martinize2 and Vermouth,
COBY, GROMACS, and ParmEd.

The software is archived on Zenodo: https://doi.org/10.5281/zenodo.20684812
