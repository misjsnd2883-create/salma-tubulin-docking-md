# Reference results of the repository code

Clean re-run of `docking/dock.py` (Vina 1.2.7, Meeko 0.7.1, RDKit 2026.9.1, Python 3.13, CPU, exhaustiveness 32, **seeds 1–3**).
Mean ± sample SD. The manuscript values use seeds 1–5; small third-decimal differences are expected.

| Ligand | Heavy atoms | Score (kcal/mol) | LE |
|---|---|---|---|
| Salma-1 | 16 | −6.196 ± 0.009 | 0.387 |
| Salma-3 | 18 | −6.504 ± 0.022 | 0.361 |
| Salma-4 | 19 | −6.653 ± 0.009 | 0.350 |
| Salma-5 | 17 | −6.327 ± 0.005 | 0.372 |
| Salma-6 | 18 | −6.474 ± 0.009 | 0.360 |
| Salma-12 | 24 | −7.696 ± 0.005 | 0.321 |
| Salma-20 | 21 | −6.991 ± 0.002 | 0.333 |
| Colchicine | 29 | −6.496 ± 0.028 | 0.224 |
| Combretastatin A-4 | 23 | −6.969 ± 0.003 | 0.303 |
| Podophyllotoxin | 30 | −8.337 ± 0.000 | 0.278 |
| Nocodazole | 21 | −7.602 ± 0.013 | 0.362 |
| CN2 re-dock (crystal geometry) | 30 | −6.706 ± 0.011 | 0.224 |

CN2 re-dock heavy-atom RMSD to the crystal pose: 1.14 ± 0.01 Å (seeds 1–3).
Colchicine here uses the crystal biaryl twist (+55°); the first-submission value (free embedding, −59°) was −6.494.

Tests: 33 passed. Each test module was also checked by mutation (22 deliberate bugs in `dock.py` and the MD scripts, all caught).

## Sensitivity to the starting conformer

The SD in the table above comes from Vina seeds only; every seed starts from the same RDKit/MMFF conformer.
Docking with Vina seed 1 from further starting conformers (embedding seed offsets 1000–4000, plus 5000 for Salma-3; all S) gave:

| Ligand | Default start (offset 0) | Other starts (range) |
|---|---|---|
| Salma-1 | −6.187 | −6.171 to −6.205 (4 starts) |
| Salma-3 | −6.523 | −6.271 to −6.287 (5 starts) |
| Salma-12 | −7.692 | −7.656 to −7.795 (4 starts) |

For flexible ligands the starting geometry changes the score by up to ~0.25 kcal/mol (Salma-3), more than the seed-to-seed SD (≤ 0.03).
The Salma-3 value in the first submission (−6.446 ± 0.005, Colab) is consistent with this: it was obtained from a different starting conformer than the one this repository generates (−6.49 ± 0.02).
Scores should therefore be read as "Vina seed SD from one starting conformer", and differences below ~0.25 kcal/mol between ligands are not interpreted.
