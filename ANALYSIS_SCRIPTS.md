# Analysis scripts of the article

Every script below reads a finished system and returns a number the article
quotes. They are listed here so that a reader who wants one number knows which
script returns it. The builder itself is `memble.sh`; these are the scripts that
measure what the builder produced.

| script | needs | what it returns |
|---|---|---|
| `count_leaflets_gro.py` | the standard library only | the phospholipid number of each leaflet, read from a `.gro` file. Run on `step6.0.gro` and `step6.6.gro` it gives the imbalance the build assigned and the imbalance the production run starts from |
| `compare_step6.py` | the standard library only | the molecules that change leaflet during the equilibration, one by one, with the head-to-tail separation and the lateral distance from the protein of each |
| `f_from_run.py` | MDAnalysis | the phospholipid number of each leaflet over the settled part of a production run, taken as the leaflet each molecule occupies in most frames |
| `settle2.py` | MDAnalysis | where the cholesterol of one finished run divided between the two leaflets, in molecules and in mol% of each leaflet |
| `leaflet_identity2.py` | MDAnalysis | whether a phospholipid that reads as crossed has crossed, or is misread |
| `ident_mdtraj2.py` | mdtraj | the same, for a machine that carries mdtraj and not MDAnalysis |
| `fit_curve_f.py` | numpy | the fit of Section 3.5 and the tests of it, with every system entered in the file |
| `protein_vs_chol.py` | MDAnalysis | the tilt of the helix, the juxtamembrane contacts, the leaflet thickness and the tail order across the systems of the curve |
| `curve2_run.sh` | GROMACS | runs a curve of six systems on two cards, a pair at a time, over as many days as it takes |

`solve_target.py`, `run_settled.sh`, `curve_stats.py`, `leaflet_relax.py`,
`leaflet_area_check.py` and the figure scripts are described in the article and
in Section S1.

## The two that need nothing installed

`count_leaflets_gro.py` and `compare_step6.py` read the `.gro` format directly
and import nothing outside the standard library, so they run on a machine that
carries neither MDAnalysis nor mdtraj.

```
python3 count_leaflets_gro.py --walk '~/counts_run/curve/*/*_work'
python3 compare_step6.py <work>/step6.0.gro <work>/step6.6.gro
```
