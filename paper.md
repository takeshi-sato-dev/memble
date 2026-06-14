---
title: 'memble: a CHARMM-GUI-free builder for Martini 3 membrane-protein systems with arbitrary lipidomes'
tags:
  - molecular dynamics
  - coarse-grained
  - Martini 3
  - membrane proteins
  - lipid bilayers
  - GROMACS
authors:
  - name: Takeshi Sato
    orcid: XXX
    affiliation: 1
affiliations:
  - name: Center for International Exchange and Department of Basic Sciences, Kyoto Pharmaceutical University, Kyoto, Japan
    index: 1
date: XXX
bibliography: paper.bib
---

# Summary

`memble` assembles solvated, ready-to-run Martini 3 coarse-grained
membrane-protein systems from an all-atom protein structure and a
user-supplied set of lipid topologies. It targets the case that general
assembly tools handle awkwardly: an asymmetric bilayer of arbitrary
composition built from the most recent Martini 3 lipid parameters, with one or
more transmembrane proteins and an optional peripheral partner protein, prepared
entirely from the command line or a small graphical interface and without a web
service. The output is a complete GROMACS input set (coordinates, topology,
index groups, position restraints, and staged minimization, equilibration, and
production run parameters), together with connectivity files that let standard
viewers display the coarse-grained system and trajectory with bonds.

# Statement of need

Coarse-grained Martini 3 simulations of membrane proteins depend on a starting
structure in which lipids are correctly oriented, virtual sites match the
reconstruction the engine performs at run time, the protein is contiguous across
periodic boundaries, and no two beads overlap. Producing such a structure is the
part of the workflow that published methods sections rarely describe. CHARMM-GUI
[@Jo2008; @Qi2015] solves it behind a web interface, and command-line packers
such as `insane` [@Wassenaar2015] embed a fixed internal lipid library. Neither
path makes it straightforward to build a system from a freshly released lipidome
that the tool does not already know, for example the refined Martini 3 lipid set
of @Pedersen2025, while keeping an asymmetric composition and a membrane-spanning
receptor.

`memble` fills that gap. It reads any Martini 3 lipid topology the user
provides, generates a three-dimensional template for each species directly from
the bonded definitions in the topology, and hands those templates to the COBY
membrane packer [@Andreasen2025], which places the bilayer and the
pre-coarse-grained protein. The protein is converted with `martinize2`
[@Kroon2024]. Because the lipid templates come from the topology rather than a
fixed library, a user can simulate a new or unusual lipid the same day its
parameters are published, which is the situation that motivated the tool:
building an epidermal growth factor receptor transmembrane-juxtamembrane system
in a refined Martini 3 lipidome with and without the ganglioside GM3.

# Functionality

A single command produces a build directory containing `system.gro`,
`system.top`, `index.ndx`, position-restraint files, minimization through
production run-parameter files, a staged run script that stops on the first
failed stage, and viewer files (`system.psf` and a connectivity PDB) that
display the coarse-grained system and trajectory with bonds in the manner of a
CHARMM-GUI session. Supported features include symmetric or asymmetric bilayers
of arbitrary composition, multiple transmembrane copies from a single
pre-oriented structure, an optional peripheral partner protein placed against a
chosen leaflet with a one-sided flat-bottom restraint, and orientation control
of the transmembrane core. A small Streamlit interface exposes the same options.

# Starting-structure failure modes and how the build avoids them

The contribution that is most useful to other practitioners is the set of
starting-structure errors the build detects and corrects, each of which
otherwise produces an infinite force at the first minimization step and is
difficult to diagnose from the engine output alone. They fall into a small
number of categories.

*Lipid orientation.* A template generated naively from bonded connectivity can
place the two acyl tails of a phospholipid on opposite sides of the head group,
so that the largest-variance axis is not the head-to-tail axis. A packer that
aligns the principal axis to the membrane normal then lays the lipid on its
side. `memble` orients each template by graph distance from the head bead and
aligns it in COBY through an explicit head-to-tail vector rather than the
principal axis. In a standalone test this produced bilayers with every lipid
head pointing outward (DLPC 196 of 196, cholesterol 256 of 256, with no flips).

*Virtual sites.* Sterol virtual sites with a large out-of-plane coefficient must
match the position the engine reconstructs each step, or minimization diverges.
The build constructs all virtual sites after the final real-bead geometry is
fixed, and reconstructs every virtual site in the assembled system from its
topology definition.

*Periodic splitting of the protein.* A transmembrane-juxtamembrane protein is
tall along the membrane normal and can straddle the box edge after assembly and
solvation. Per-atom wrapping upstream then places consecutive backbone beads a
full box apart, the bond constraint cannot be satisfied, and the constraint
solver diverges. The build makes each protein chain contiguous across the
boundary and centers the system on the protein before any further step.

*Operation order and residual overlap.* Any coordinate operation that runs after
clash removal can reintroduce overlap. The build performs all coordinate and
periodic-boundary operations first and runs clash removal last, with the protein
frozen so that only solvent and lipids move. A final check reports the smallest
intra-molecular and inter-molecular bead-bead distances under the minimum image
convention, so a system that still contains an overlap is identified before the
engine runs.

# Quality control

`memble` ships with a unit-test suite (26 tests) covering template
generation, virtual-site reconstruction, residue-number restoration, protein
unwrapping, clash removal with a frozen protein, and the connectivity writers.
The lipidome parameters themselves are not distributed; the documentation
directs users to obtain the refined Martini 3 lipid set [@Pedersen2025] and pass
it at build time.

# Acknowledgements

XXX

# References
