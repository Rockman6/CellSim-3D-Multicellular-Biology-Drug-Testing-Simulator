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
RANK  NAME                     TRIAGE         ΔG(kcal)   K_d        POCKET  STRAIN       Ro5  QED   logS
   1  biotin_TRUE_BINDER       follow_up         -7.44    3.5 µM    ✓       acceptable  ✓    0.49  -1.53
   2  ibuprofen_negative       deprioritise      -7.36    4.1 µM    ?       good        ✓    0.82  -3.09
   3  aspirin_negative         drop              -6.66   13.0 µM    ✓       acceptable  ✓    0.55  -1.99
   4  acetaminophen_negative   drop              -6.31   23.8 µM    ?       acceptable  ✓    0.59  -1.97
   5  caffeine_negative        drop              -5.48   95.7 µM    ?       good        ✓    0.54  -0.87
```

The `TRIAGE` column (`follow_up` / `review` / `deprioritise` /
`drop`) is the one-decision column: CellSim synthesises ΔG +
pose-strain + PoseBusters + ADMET into one verdict with a
paste-ready reason string. Wet-lab users read one column, not
five booleans.

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
| **1.5 Coarse-grained** | Martini 3 membrane / bilayer MD | ⏳ scaffold only |
| **1.6 UQ** | Monte-Carlo / Sobol / split-conformal for ΔG bounds | ✅ triad shipped |
| **1.7 Blind harness** | PDBBind scale gate + red-team slot | ⏳ 3-cocrystal mini-bench shipped; PDBBind scale pending |
| **x-cut cache** | SQLite physics-prior memoisation (AM1-BCC, Vina, xTB) | ✅ shipped; wired into Layers 1.1 / 1.3 / 1.4 |

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
- `src/cache/` — content-addressed SQLite cache for every
  expensive physics call; pass `--cache /path/to/c.sqlite` to
  `cellsim dock` and repeat runs short-circuit identical
  (ligand, receptor, box, seed) tuples (~1000× per-compound
  Vina speedup, ~19× end-to-end on a small screen).

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
the `cellsim` conda env and runs **34 gates** in ~15 min, grouped
by layer:

- **Layer 1.1 chem** — 10-drug parametrise, ADMET, AM1-BCC cache.
- **Layer 1.2 MD** — 3-drug vacuum MD, 1UBQ solvate + 1 ps MD.
- **Layer 1.3 docking + triage** — Meeko prep, 1STP re-dock,
  fpocket detect, 3-cocrystal mini-bench, refine, Vina cache,
  pose SDF/PDB export, off-target selectivity, strain diagnostic
  (UFF-ensemble), strain-gate top-pose promotion, triage rule
  table, shortlist filter, triage-PNG dashboard, kinase-receptor
  heads-up, CYP3A4 DDI-risk strain-downgrade.
- **Layer 1.4 quantum** — 10-drug xTB, 3-drug CYP3A4 BDE SoM,
  heme-accessibility SoM, PySCF DFT single-point.
- **Layer 1.6 UQ** — MC-dock, conformal, streptavidin + EGFR
  calibrations.
- **Cross-cut** — cache round-trip.

A regression in any gate blocks merge. See
[`BENCHMARKS.md`](BENCHMARKS.md) for the numeric targets each
gate enforces.

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

## How good is it?

See [`BENCHMARKS.md`](BENCHMARKS.md) for the full scorecard —
every current number (pose-recovery %, Spearman on two receptor
families, cache speed-ups, CI gate list) with honest caveats and
per-row reproducers.

Headline numbers today:

- **3-cocrystal mini-bench:** 2/3 = 67 % pose recovery at canonical
  top-3 < 2 Å gate.
- **Streptavidin calibration:** Spearman ρ = **1.00** across
  14 orders of magnitude of K_d (but Pearson r = 0.98 hides the
  fact that Vina's absolute ΔG saturates on tight binders — MAE
  4.99 kcal/mol).
- **Trypsin calibration:** MAE = **0.91 kcal/mol** on benzamidine
  analogs (Vina's absolute is accurate here), but Spearman
  ρ = 0.09 within the narrow 4 kcal/mol window (noise floor).
- **End-to-end cache speed-up:** 19× on a 5-compound batch rerun.

## How to cite

Every prediction in CellSim carries its method provenance (tool
version, force field, seed, search parameters). To cite a run,
include the commit SHA and the `provenance` block from the relevant
result envelope in your publication's methods.

## License

MIT. See [`LICENSE`](LICENSE).
