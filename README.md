# CellSim

**Open-source, non-AI, physics-first drug-discovery triage platform.**

CellSim is the in-silico front of a wet-lab shortlist. Drop in a
receptor PDB and a list of SMILES, walk away for coffee, come back
to a ranked CSV with ΔG ± CI, pocket-fit flags, ADMET descriptors,
and one-page drug profile dashboards for your top hits.

Every method is physics-grounded (classical force fields, semi-
empirical QM, alchemical free energy, flux-balance analysis) and
every rate constant cites a PMID or a cached physics calculation.
No neural scoring, no learned potentials, no black-box surrogates —
see [`MISSION.md`](MISSION.md) for the discipline and
[`GOAL`](GOAL) for the 5-campaign roadmap.

> **Status:** Campaign 1 (Atomic → Molecular Foundation) is in
> progress. The pre-restart HeLa / p53 / cisplatin cell prototype
> is frozen under [`OLD/`](OLD/) as a regression snapshot (still
> builds; still passes its 8 headless validators).

## Quickstart — first docking run in 5 minutes

```bash
# 1. Install the environment (one-time).
conda env create -f environment.yml      # or mamba
conda activate cellsim

# 2. Run the end-to-end biologist workflow on a bundled cocrystal.
python -m src.dock.batch \
    --smi benchmarks/dock/1stp_batch_5.smi \
    --receptor benchmarks/dock/1stp.pdb \
    --out-csv /tmp/run/report.csv \
    --mc 4 --profile-top-k 3 \
    --crystal-pdb benchmarks/dock/1stp.pdb \
    --crystal-resname BTN
```

Output:

```
RANK  NAME                    ΔG(kcal)   K_d       POCKET  RMSD    Ro5  QED   logS
   1  biotin_TRUE_BINDER        -7.44    3.5 µM    ✓       1.96 Å  ✓    0.49  -1.53
   2  ibuprofen_negative        -7.36    4.1 µM    ?       -       ✓    0.82  -3.09
   3  aspirin_negative          -6.66   13.0 µM    ✓       -       ✓    0.55  -1.99
   4  acetaminophen_negative    -6.31   23.8 µM    ?       -       ✓    0.59  -1.97
   5  caffeine_negative         -5.48   95.7 µM    ?       -       ✓    0.54  -0.87
```

Plus `profile_01_biotin_TRUE_BINDER.png`, `profile_02_...png`,
`profile_03_...png` — one-page dashboards showing 3D + charges,
predicted CYP3A4 sites-of-metabolism, HOMO/LUMO, and the full
Lipinski / QED / logS datasheet.

If you do **not** know the binding-site coordinates of your target,
omit `--center` / `--box` and CellSim auto-detects pockets via fpocket
(the canonical non-ML geometric pocket finder).

## What CellSim can do today

| Layer | What it does | Status |
|---|---|---|
| **1.1 Chem** | SMILES → OpenFF-parametrised system (AM1-BCC charges) | ✅ 9/10 full tier, 10/10 RDKit tier |
| **1.2 MD** | Classical Langevin MD, solvated protein loader (AMBER14 + TIP3P) | ✅ 1 ps ubiquitin Cα RMSD 0.74 Å |
| **1.3 Docking** | Vina + Meeko + PoseBusters + fpocket auto-site | ✅ mini-bench 2/3 canonical gate |
| **1.4 Quantum** | xTB GFN2 single-point + CYP3A4 SoM predictor (BDE) | ✅ 10/10 sane + 3/3 SoM smoke |
| **1.5 Coarse-grained** | Martini 3 membrane / bilayer MD | ⏳ scaffold pending |
| **1.6 UQ** | Monte-Carlo over Vina seeds → ΔG ± 95 % CI | ✅ MVP shipped; Sobol pending |
| **1.7 Blind harness** | PDBBind scale gate + red-team slot | ⏳ scale harness pending |

**Cross-cutting UX:**

- `src/dock/batch.py` — one-command ranked screen with MC error bars
  and optional `--profile-top-k` auto-dashboards.
- `src/chem/profile.py` — six-panel per-compound profile combining
  3D charges, SoM predictions, HOMO/LUMO, BDE chart, and
  Ro5 / QED / logS callouts.
- `src/chem/admet.py` — Lipinski / TPSA / QED / ESOL solubility
  (all published formulae; no ML).
- `src/dock/pocket_detect.py` — auto-binding-site detection so any
  receptor PDB works out-of-box.
- `src/uq/dock_mc.py` — honest ΔG ± CI from N-seed Monte-Carlo.

## Non-AI discipline

CellSim is strictly physics-first and non-AI by design. This is
load-bearing, not a preference. Every prediction must trace to a
physics calculation or a literature-cited empirical formula. See
[`MISSION.md`](MISSION.md) §"No black-box / no AI surrogates" for
the five ground rules and the narrow "ML as accelerator only"
exception clause.

Explicitly excluded:

- ML potentials (MACE / NequIP / OrbNet / Allegro) as the force path.
- GNINA CNN-scored docking as the primary evidence (Vina only; GNINA
  may ship as an explicitly labeled fast-guess alongside).
- "Deep ensembles" for UQ. Sobol + Monte-Carlo + MAPIE conformal
  (post-hoc, non-parametric) only.

## Validation that runs on every PR

[`.github/workflows/smoke.yml`](.github/workflows/smoke.yml) provisions
the `cellsim` conda env and runs the full gate stack in ~15 min:

1. `tests/chem/test_parametrize_smoke.py` — 10 canonical drugs →
   OpenMM system (Layer 1.1).
2. `tests/chem/test_admet_smoke.py` — Lipinski / QED / logS sanity.
3. `tests/md/test_ligand_vacuum.py --max 3 --gate 3` — 10-ps vacuum
   Langevin (Layer 1.2).
4. `tests/md/test_protein_load.py` — 1UBQ load + solvate + minimise +
   1-ps MD.
5. `tests/dock/test_prep_smoke.py` — Meeko PDBQT prep.
6. `tests/dock/test_redocking.py` — 1STP biotin top-1 ≤ 2.5 Å AND
   top-3 best < 2.0 Å.
7. `tests/dock/test_pocket_detect_smoke.py` — fpocket finds biotin
   pocket within 3 Å.
8. `tests/dock/test_mini_bench.py` — 3-cocrystal aggregate pose
   recovery ≥ 66 %.
9. `tests/quantum/test_xtb_smoke.py` — 10-drug GFN2 sanity.
10. `tests/quantum/test_som_smoke.py` — 3-drug CYP3A4 BDE smoke.
11. `tests/uq/test_mc_dock_smoke.py` — 4-seed MC dock sanity.

A regression in any gate blocks merge.

## Layout

```
src/
  chem/      Layer 1.1  parametrise, ADMET, profile dashboard
  md/        Layer 1.2  OpenMM MD driver, PDBFixer protein loader
  dock/      Layer 1.3  Vina + Meeko + PoseBusters + fpocket + batch
  quantum/   Layer 1.4  xTB GFN2, CYP3A4 SoM predictor
  cg/        Layer 1.5  Martini 3 (scaffold; not populated yet)
  cache/     cross-cut  SQLite + HDF5 physics-prior cache (scaffold)
  uq/        Layer 1.6  Monte-Carlo dock, Sobol (pending)
  bridge/    cross-cut  Layer-1 → Layer-2 rate-law emitter (future)
  core/      physics-neutral RNG, telemetry, constants
benchmarks/
  chembl/    10 canonical drugs
  md/        1UBQ ubiquitin
  dock/      1STP / 1M17 / 3PTB cocrystals + mini_bench.yaml
tests/       one folder per src/ module, same smoke pattern
OLD/         frozen pre-restart Campaign-2 prototype (builds; passes
             its own 8 headless validators as regression snapshot)
docs/        strategic plan, professor debriefs, campaign scope
scripts/     CLI utilities (fetch_pdb, fetch_chembl_sample, ...)
```

## How to cite

Every prediction in CellSim carries its method provenance (tool
version, force field, seed, search parameters). To cite a run,
include the commit SHA and the `provenance` block from the relevant
result envelope in your publication's methods.

## License

MIT. See [`LICENSE`](LICENSE).
