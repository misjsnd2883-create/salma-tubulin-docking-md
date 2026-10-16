# Salma series: docking, MD and analysis code

Code for the manuscript on the Salma series (Salma-12) at the colchicine site of tubulin
(PDB 1SA0). Manuscript in preparation. Preprint: https://doi.org/10.26434/chemrxiv.15005087/v2

## Contents

| Folder | Purpose |
|---|---|
| `docking/dock.py` | Receptor preparation (Meeko), grid maps, docking with AutoDock Vina 1.2.7, CN2 re-dock, results table (`results.json`) |
| `docking/export_pose.py` | Docked pose to ligand SDF and protein-ligand complex PDB (input for CHARMM-GUI) |
| `md/run_md.py` | Resumable equilibration and production run of a CHARMM-GUI/OpenMM system (segments of nstep × dt) |
| `md/check_system.py` | Sanity checks of a CHARMM-GUI system (atom counts, box, net charge, HMR, time step) |
| `analysis/md_analysis.py`, `analysis/plot_md.py` | Ligand RMSD, protein RMSD/RMSF, contact occupancy, polar contacts, figures |
| `docs/` | `reference_results.md` (re-run docking table), `charmm_gui_settings.md` (system build and MD protocol) |
| `tests/` | pytest suite |

## Protocol in brief

- Receptor: 1SA0 chains A+B with GTP/GDP and Mg2+; reference ligand CN2 (DAMA-colchicine); box centre = CN2 centroid (117.219, 90.180, 6.290), 22 × 22 × 22 Å.
- Docking: AutoDock Vina 1.2.7, exhaustiveness 32, seeds 1–5, Meeko 0.7.1 (`-a`). Reported spread is the **sample SD** (n−1) over seeds.
- Stereochemistry: undefined centres are set explicitly (S) before embedding; RDKit would otherwise pick an enantiomer at random.
- Colchicine biaryl twist: the default embeds the twist observed in the crystal CN2 (+55°). `--colchicine-twist free` reproduces the unconstrained embedding used in the first submission (−59°); the scores are statistically indistinguishable.
- MD: CHARMM-GUI + OpenMM, CHARMM36m / CGenFF, TIP3P, HMR, 310.15 K, 4 fs. Production = 50 segments × 1 ns (250,000 steps) = 50 ns, frames every 0.1 ns.
- Analysis: heavy-atom ligand RMSD after Cα superposition; contacts ≤ 4.0 Å; polar contacts ≤ 3.5 Å; Cα RMSF.

## Install and run

```bash
# Colab / Ubuntu, Python 3.13 (needed only to build the vina wheel)
apt-get install -y swig libboost-all-dev
pip install -r requirements-docking.txt

python docking/dock.py --workdir work --pdb 1SA0.pdb          # all ligands, seeds 1-5
python docking/dock.py --workdir work --only Salma-12 --seeds 1 2 3
python docking/export_pose.py --pose work/Salma-12_s2.pdbqt --site work/1SA0_site.pdb --out work/S12 --anion

python md/check_system.py --help
python md/run_md.py --system-dir <charmm-gui>/openmm --rep-dir reps/rep1 --n-cycles 50

pip install -r requirements-analysis.txt
python analysis/md_analysis.py --help
```

`dock.py` is resumable (`results.json`). The input PDB must be downloaded from RCSB (not included).

## Tests

```bash
pip install -r requirements-dev.txt
SALMA_TEST_PDB=/path/to/1SA0.pdb python -m pytest tests -q
```

Tests that need 1SA0 are skipped without `SALMA_TEST_PDB`. The MD scripts are tested with a stub `openmm_run.py` and synthetic CHARMM-GUI files only; the production simulations themselves were run on a Tesla T4 (OpenCL, OpenMM 8.6 / 8.6.1) and are not reproduced by the tests.

## Reproducibility notes

- The trajectories in the paper were produced with the author's Colab notebooks (same `.inp` files). `md/run_md.py` is a cleaned, restartable version of that loop; it has been tested with a stub only, not on a full simulation. The analysis scripts were run on the real trajectories and reproduce the reported values.

- Docking scores depend on the starting conformer of flexible ligands (up to ~0.25 kcal/mol, see `docs/reference_results.md`) and can differ by a few hundredths between CPU types and Vina builds. The reported SD covers Vina seeds only.
- Meeko removes 27 residues that fail template matching; the nearest (B:322) is 1.4 Å from the box edge and ≥ 14 Å from CN2 and all poses, beyond Vina's 8 Å cutoff. `dock.py` checks this and records the distance for every pose.
- Trajectories are not stored in this repository (size). Link to the archive (e.g. Zenodo): TODO.

## Use of AI assistance

The scientific design, choice of targets and ligands, protocol, and interpretation are the author's. AI assistants (Claude Opus 5 and Claude Sonnet 5.5, Anthropic) were used to review, debug and test the scripts and to help write the documentation. All code was run and checked by the author; the results in the manuscript were reproduced with these scripts.

## Licence

MIT (see `LICENSE`). Third-party components: AutoDock Vina (Apache-2.0), Meeko (LGPL-2.1), RDKit (BSD), MDTraj (LGPL-2.1), OpenMM (MIT/LGPL).

## Citation

TODO
