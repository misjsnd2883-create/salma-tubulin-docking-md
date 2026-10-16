# CHARMM-GUI system set-up and MD protocol

Systems were built with CHARMM-GUI (Solution Builder / Input Generator) and run with OpenMM.
Values below are taken from the generated `step3_input.psf`, `sysinfo.dat` and the `.inp` files of the two runs.

## Systems

| | Salma-12 | CN2 (control) |
|---|---|---|
| Atoms | 77,701 | 77,732 |
| Protein | α-tubulin 436 res (2–437), β-tubulin 427 res (2–438) | same |
| Cofactors | GTP, GDP, 2 Mg2+ | same |
| Ligand | LIG, 44 atoms (20 H), carboxylate anion, S | CN2, 55 atoms (30 heavy) |
| Water (TIP3) | 21,365 | 21,372 |
| Ions | 90 Na+, 59 Cl− | 89 Na+, 59 Cl− |
| Box (Å) | 124 × 82 × 82 | 124 × 82 × 82 |

- Force field: CHARMM36m (protein), CGenFF (ligands; Salma-12 penalties: parameter 7.2, charge 4.309), TIP3P water.
- Hydrogen mass repartitioning (hydrogen mass 3.024), constraints on bonds to H.
- Chain termini: NH3+ at the start, N-methylamide cap at the end.
- Input complex: best docking pose (Salma-12: seed 2), built with `docking/export_pose.py`.

## Protocol (identical in both runs; `step4_equilibration.inp`, `step5_production.inp`)

- Equilibration: 5,000 minimisation steps, then 125,000 steps × 1 fs without pressure coupling; positional restraints (backbone 400, side chain 40 kJ/mol/nm²); velocities at 310.15 K.
- Production: 50 segments × 250,000 steps × 4 fs = 50 ns; no restraints; Langevin 310.15 K, friction 1 ps⁻¹; Monte Carlo barostat 1 bar every 100 steps; PME (tolerance 0.0005); force-switch vdW 1.0–1.2 nm; frames every 100 ps (500 per run).
- Hardware/software: Tesla T4, OpenCL; OpenMM 8.6 (Salma-12), 8.6.1 (CN2).
- Single replicate per system.

## Not recorded here (fill in from the CHARMM-GUI job page if you want it in the repository)

- Salt concentration and ion-placing method, water-box padding and box type, CHARMM-GUI version/job IDs.
