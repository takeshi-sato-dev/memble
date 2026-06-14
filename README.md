# memble

A composition-general builder for Martini 3 membrane-protein systems, producing
CHARMM-GUI-equivalent inputs (GROMACS topology, CHARMM PSF and CRD, a staged
equilibration protocol, and an index file) without CHARMM-GUI. The tool takes an
all-atom protein structure and a lipid composition and returns a system that is
ready to minimize, equilibrate, and run in GROMACS 2023.

> Working name. Rename the repository, the command, and the paper title to your
> final project name before release.

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
- A Martini 3 lipidome distribution (the itp files defining your lipids, water,
  and ions)

Install the Python tools, for example:

```
pip install numpy parmed COBY vermouth
```

DSSP and GROMACS are installed separately through your package manager or module
system.

## Installation

```
git clone https://github.com/XXX/memble
cd memble
chmod +x memble.sh
```

The four helper scripts (`orient_tm.py`, `replicate_and_fix_top.py`,
`inject_posres.py`, `leaflet_area_check.py`) live at the repository root.

## Setup (one command) and no-activate use

Create a dedicated virtual environment so nothing in your existing setup
changes. This is safe to rerun: it rebuilds the venv from scratch.

```
bash setup.sh
```

Then run without activating the venv each time:

```
M3_DIR=/path/to/martini3_lipidome GMX=gmx LIPIDS="CHOL:1 DIPC:1 DPSM:1" \
  ./memble protein_AA.pdb        # single build (or ensemble if N_REP > 1)

./memble-gui                           # launch the GUI
```

`memble` and `memble-gui` use the venv interpreter and tools directly and wire the
helper paths for you, so no `source .../activate` is needed. The venv lives at
`~/memble-venv` (override with `MEMBLE_VENV`); delete that folder to reset.

## Quick start

Symmetric ternary raft mixture, four protein copies:

```
LIPIDS="CHOL:1 DIPC:1 DPSM:1" \
M3_DIR=/path/to/martini3_lipidome \
GMX=gmx \
HELPER_REP=$PWD/replicate_and_fix_top.py \
HELPER_POS=$PWD/inject_posres.py \
HELPER_ORI=$PWD/orient_tm.py \
HELPER_AREA=$PWD/leaflet_area_check.py \
./memble.sh protein_AA.pdb

cd memble_work
./run.sh
```

Asymmetric bilayer, explicit box, thicker water, lower temperature:

```
UPPER="CHOL:1 DPSM:2 POPC:1" LOWER="POPC:3 POPE:1" \
BOX_X=60 BOX_Y=50 WATER_NM=3 TEMP=300 \
M3_DIR=... GMX=gmx HELPER_REP=... HELPER_POS=... HELPER_ORI=... HELPER_AREA=... \
./memble.sh protein_AA.pdb
```

If the two leaflets are not area-balanced, the build stops before any molecular
dynamics is run and reports the per-leaflet areas, so the composition or box can
be adjusted and the system rebuilt.

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
| `SS_MODE` | `dssp` | `dssp` lets DSSP assign secondary structure (GPCRs, multi-helix proteins); `tm` sets the TM range to helix and everything else to coil, so the juxtamembrane stays flexible (single-pass TM and TM-JM peptides) |
| `DSSP` | auto | path to a DSSP 2.2.1/3.x binary, or `mdtraj`; 4.x is not compatible with martinize2 |
| `TM_CORE` | unset | per-chain TM core, e.g. `A:65-88;B:65-88`; used by `SS_MODE=tm` and by prebuilt orientation |
| `PREBUILT_MULTI` | 0 | `1` uses the chains already assembled in the input PDB instead of replicating one chain |
| `RES_KEEP` | unset | per-chain residues to keep with `PREBUILT_MULTI`, e.g. `A:54-103;B:54-103` |
| `AREA_TOL` | 0.08 | maximum leaflet area mismatch before abort |
| `NPROD_STEPS` | 4e8 | production steps at 20 fs (8 microseconds) |

### Choosing the secondary structure source

DSSP reads secondary structure from the input coordinates. If the input has a
helical juxtamembrane, DSSP marks it helix, and the backbone dihedral restraints
then hold that region rigid. For a single-pass TM or TM-JM peptide where the
juxtamembrane should stay flexible, set `SS_MODE=tm`: the TM range becomes helix
and everything else coil, independent of the input conformation and of DSSP. For
multi-helix proteins such as GPCRs, keep `SS_MODE=dssp` so the loops, the
membrane helices, and any peripheral helix are each assigned correctly.

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

## How the equilibration avoids early instability

The six-stage equilibration releases the protein-backbone and lipid-head
restraints gradually while the time step rises from 2 to 20 fs, the pressure is
controlled with the c-rescale barostat until the bilayer settles, and the
Parrinello-Rahman barostat is used only for production. The run script checks
for the step coordinate files that GROMACS writes when atoms move too far and
stops if the system is destabilizing.

## Tests

```
python3 -m pytest tests/ -v
```

The suite exercises each helper on synthetic inputs: orientation tilt, grid
placement and separation, automatic and explicit head-bead restraint injection,
and the leaflet area check for both a balanced and a mismatched bilayer.

## Limitations

- A lipid can only be placed if it exists in the supplied Martini 3 lipidome.
- The leaflet area check uses approximate area-per-lipid values and is a coarse
  pre-run guard, not an exact measurement. Override values with `APL_OVERRIDE`.
- Multi-headgroup lipids (cardiolipin, phosphoinositides, gangliosides) may need
  an explicit head bead given in the composition string.
- The COBY membrane-string grammar (lipid ratio token, leaflet token) can depend
  on the COBY version. Confirm it with `COBY -h`, set `COBY_LEAFLET` if needed,
  or override the whole string with `COBY_MEMBRANE`.

## License

Released under the Apache License 2.0. See `LICENSE`.

## Contributing

Issues and pull requests are welcome. See `CONTRIBUTING.md`.

## Citation

If this tool is useful in your work, please cite the accompanying paper (see
`paper.md`) and the underlying tools: Martini 3, martinize2 and Vermouth, COBY,
GROMACS, and ParmEd.
