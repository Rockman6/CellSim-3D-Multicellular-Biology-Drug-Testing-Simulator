# CellSim

**A physics-first, uncertainty-honest in-silico cell and drug-response
simulator in Python. Open source, Apple-silicon and Linux, no GPU
required.**

**[Open a simulated spheroid in your browser →](https://rockman6.github.io/CellSim-3D-Multicellular-Biology-Drug-Testing-Simulator/)**
(no install; the demo runs are simulated fresh on every deploy, so they
show what the current code produces)

```bash
pip install cellsim
```

CellSim is converging on one product: an interactive dish. Pick a cell
line and a drug, set a dose schedule, watch the colony respond, and get a
dose-response curve whose error bar has been checked against public
cell-line data. The molecular layer (docking, ADMET, quantum descriptors)
feeds that engine as an input with provenance and a measured accuracy
flag; it is not the product on its own.

> **Status (October 2026).** Phases 0–2 of the [plan](docs/PLAN.md) are
> done and Phase 3 is in review. The smoke suite is green, the repository
> has one identity, and `cellsim/cell/` holds a validated drug-response
> engine that runs well-mixed or in space.
>
> **Phase 1 closed with a negative result that set the direction.** The
> engine reproduces each drug's potency scale, but predicting *which
> cell line* is more sensitive turned out to be unreachable: eleven
> canonical markers, fitted optimally across 156–683 lines, explain only
> 9–20 % of the difference between lines. So per-line ranking is out of
> scope, and Phase 2 aimed at what mechanism is actually good for —
> **within-line dynamics**: dose schedule, wash-out, resistance under
> selection, combination timing.
>
> **Phase 2 built the dish.** `cellsim dish` grows those cells in space —
> a monolayer with contact inhibition, or a spheroid with oxygen and drug
> diffusing in, hypoxic quiescence and a necrotic core — and
> [`web/viewer/`](web/viewer/) renders the per-cell stream in a browser.
> Every Phase-2 claim was re-measured over seeds and doses before the
> dish was built on it; one was retracted and one qualified.
>
> **Phase 3 found the exception to the Phase-1 result.** Across ten drugs
> and 91 held-out predictions the conclusion holds — except for the one
> drug whose target is part of the modelled mechanism. An MDM2 inhibitor
> works by stabilising p53, and the engine calls all eleven lines
> correctly (a constant gets five), including HeLa, whose TP53 is
> wild-type by sequence but whose p53 is degraded by HPV E6 rather than
> by MDM2. So: no ranking lines for cytotoxic drugs, but yes for a drug
> whose mechanism is modelled. Every number is in
> [`docs/VALIDATION.md`](docs/VALIDATION.md), misses included.
>
> The 2026 C++/Metal prototype lives in [`OLD/`](OLD/) as a biology
> reference; it still builds and passes its 8 headless benchmarks, but
> its UI is retired and a new interface is planned on top of the Python
> engine.

## What works today

| Layer | What you can run | How it is validated |
|---|---|---|
| Docking triage (`cellsim dock`) | SMILES list + receptor PDB → ranked CSV with ΔG ± CI, pose-trust (UFF strain), PoseBusters, ADMET, one triage column | 15-cocrystal blind re-dock, structures and SMILES fetched from RCSB: top-3 pose recovery 87 %, top-1 73 %, PoseBusters validity 100 % |
| Target-class reliability (`cellsim/uq`) | Measured docking error per receptor family, so every ΔG carries an **accuracy** flag, not just seed scatter | trypsin-like 0.9, kinase ATP-site 2.2, ultra-tight binders 5.0 kcal/mol MAE (n = 4–6 each) |
| ADMET descriptors (`cellsim admet`, `profile`) | Lipinski, TPSA, QED, ESOL logS, BBB / hERG / Ames rule flags, one-page profile PNG | Published formulae, cited at the point of use in `cellsim/chem/admet.py` |
| CYP3A4 site of metabolism (`cellsim som`) | xTB C–H bond-dissociation ranking with a heme-accessibility re-rank | 2/3 on the bundled literature set; blind to N-dealkylation. **Advisory only.** |
| Knock a gene out (`api.knockout`) | 13 pathway genes — knockout, knockdown or overexpression — and what the rest of the network does without them | TP53 loss lands where the p53-mutant lines are; PUMA/BAX/CASP3 knockouts are identical because they are one pathway; ABCB1 protects doxorubicin and **not** cisplatin, which it does not transport |
| Starve a tumour (`api.starve`) | Glucose as a second diffusing nutrient; whichever of it and oxygen is scarcer sets the viable rim | **Fails its validation and is labelled so**: across the published 0.8–16.5 mM range the model gives identical answers, because the centre never falls below its viability threshold at reachable spheroid sizes ([details](docs/VALIDATION.md)) |
| Cells that crawl (`DishParams.migration_*`) | A random walk on the lattice, optionally up the oxygen gradient | Matches D = k·h²/4 to within 1 % at 2 hops/h |
| Immune killing (`api.killing`) | Effector cells that crawl, engage, kill and tire, on a tumour monolayer | Specific lysis 0 → 37 → 88 → 100 % across E:T ratios; kill counts fall at high E:T as targets run out. Kinetics are literature values — calibrated, not validated |
| Your own plate (`cellsim plate`) | Read a plate-reader export (grid or tidy), normalise to the plate's own controls, fit a dose-response curve, and get GR metrics comparable across lines | Coverage measured against known IC50s: a nominal 95 % interval contains the truth ~90 % of the time on ten concentrations — and only **68 % on five**, which is flagged. Flat plates are refused; ranges that clip the response are marked extrapolated |
| Calibrated uncertainty (`scripts/calibrate_cell_uncertainty.py`) | Every predicted IC50 comes with an interval whose coverage has been scored | Leave-one-drug-out: ±2.6× covers 64 % (68 % nominal), ±4.8× covers 85 % (90 %), ±6.5× covers 95 % (95 %). Wide on purpose — GDSC's own replicate screens disagree by a median 5.3× |
| Cell drug-response engine (`cellsim cell-sim`, `cellsim cell-stream`) | Pick a cell line, a drug and a **dose schedule**; get the colony over time, viability, and a per-cell state stream (JSON Lines) for a UI. Built for *within-line* questions: schedule, wash-out, resistance, combinations | Each drug's potency scale is right (A549 cisplatin 10.9 µM vs 9.8 measured; 20/25 held-out predictions inside GDSC's replicate span). **Ranking one cell line against another is explicitly out of scope** — measured as unreachable from canonical markers ([why](docs/VALIDATION.md)). |
| Spatial dish (`cellsim dish`) and web viewer (`web/viewer/`) | Grow cells as a monolayer or a spheroid under a dose schedule; oxygen and drug fields, contact inhibition, hypoxic quiescence, necrosis, clones. Open the JSON Lines stream in the viewer: 3-D cells coloured by phase, p53, caspase-3, oxygen or lineage, time scrubbing, charts, CSV export | Field solvers match closed forms (Grimes's 233 µm oxygen limit; the analytic slab profile). Re-running the Cell Tracking Challenge HeLa movie: colony size within 3–4 % at 46 h. Spheroid growth **calibrated** on DLD-1 (one constant). **Open misses:** the spheroid's necrotic core is too large, and real HeLa splits into fast cyclers and non-dividers, which the engine does not ([details](docs/VALIDATION.md)). |
| Cell-level PK/PD modules (`cellsim/cell`) | Occupancy, permeation, pH trapping, efflux, binding sink, tissue penetration, fate, resistance, clearance, cell cycle, lattice agents | Each module is checked against its analytic limit (22 tests). |
| Alchemical FEP (`cellsim/fep`) | Hydration and binding ΔG scaffolds on openmmtools + MBAR | **Experimental, not a product path.** FreeSolv-12 MAE 1.42 kcal/mol on 10/12 with a size-dependent bias (no barostat). Binding ΔG has never produced a number on a real binder. |

The full set of numbers, with reproducers and caveats, is in
[`docs/VALIDATION.md`](docs/VALIDATION.md).

## Quickstart: simulate cells in one minute

The cell simulator needs only NumPy and SciPy, so it installs with plain
pip — no conda, no compilers, no chemistry stack. (Checked in CI on a bare
virtual environment: `tests/test_pip_install.py`.)

```bash
pip install cellsim                 # or: pip install -e . from a checkout

# A dose-response curve and an IC50 for one line and drug.
cellsim cell-sim --line A549 --drug cisplatin --hours 72

# Grow a spheroid for five days and watch oxygen run out in its core.
cellsim dish --line DLD-1 --geometry spheroid --hours 120 \
    --cells 3000 --grid 48 --every 8 --out spheroid.jsonl

# Open spheroid.jsonl in web/viewer/index.html (no server needed for a
# local file: use the "Open .jsonl" button).
```

In a notebook, `cellsim.api` returns tidy tables (a `pandas.DataFrame`
when pandas is installed, a list of dicts otherwise):

```python
from cellsim.api import curve, exposure, washout, combination, spheroid

curve("A549", "paclitaxel")               # viability vs concentration, with IC50
exposure("A549", "paclitaxel")            # concentration needed for half kill at
                                          # each exposure time — the schedule question
washout("A549", "cisplatin", 20.0, 12.0)  # a 12 h pulse, then recovery
combination("A549", "cisplatin", "paclitaxel", 10.0, 0.03)   # vs independent action
spheroid("DLD-1", days=6)                 # size, hypoxia, necrotic core over time
```

[`docs/cell_tutorial.md`](docs/cell_tutorial.md) walks through a real
question end to end — *should I pulse this drug or leave it on?* — with
every command run and its actual output.

`ic50_spread()` is worth running before quoting any IC50 from this
engine: the number moves 10–20 % between seeds, because a population's
IC50 is set by which resistant lineages were drawn rather than by how
many cells were counted ([why](docs/VALIDATION.md)). Two numbers that
differ by less than that do not differ.

`exposure()` is the one worth knowing about: it returns `inf` where no
concentration reaches half kill, which is the honest answer for a taxane
given for three hours, and the reason AUC is the wrong exposure metric
for that class.

The molecular layers (docking, MD, FEP, quantum) need the conda stack
below.

## Quickstart: first docking screen in five minutes

```bash
# 1. One-time environment (conda or mamba).
mamba env create -f environment.yml
conda activate cellsim
pip install -e .                    # registers the `cellsim` command
cellsim doctor                      # expect "42/42 checks passed"

# 2. Screen the bundled 5-compound batch against streptavidin (PDB 1STP).
cellsim dock \
    --smi benchmarks/dock/1stp_batch_5.smi \
    --receptor benchmarks/dock/1stp.pdb \
    --out-csv /tmp/run/report.csv \
    --mc 4 --profile-top-k 3 \
    --crystal-pdb benchmarks/dock/1stp.pdb --crystal-resname BTN
```

You get a ranked CSV (`triage`, `dG_kcalmol`, `Kd_human`, `strain_band`,
`pocket_ok`, ADMET columns) plus one-page profile PNGs for the top hits.
Omit `--center` / `--box` for any other receptor and fpocket finds the
site. `cellsim help` lists every subcommand;
[`docs/docking_tutorial.md`](docs/docking_tutorial.md) walks through a
full screen and, in §8, says which target classes to trust.

## Where this is going

| Phase | Deliverable | Target |
|---|---|---|
| 0 — done | July fixes merged, CI green, one identity, honest validation page | Oct 2026 |
| 1 — closed, with a finding | Single-cell drug-response engine with one state object and one integrator: cell cycle and p53 axis ported from `OLD/`, death that actually executes, `cellsim/cell` PK on top. Calibrated against GDSC for cisplatin, doxorubicin and paclitaxel (ten drugs by Phase 3). Headless state stream for any UI. Closed early: the per-line prediction it aimed at was measured to be unreachable | Oct 2026 |
| 2 — code done | The dish: spatial colony, diffusion, drug penetration, contact inhibition, validated on the CTC HeLa movie and calibrated on DLD-1 spheroids; web viewer reading the engine stream. Outside users still to find | Apr 2027 |
| 3 — in review | Ten drugs × ten lines, combinations and evolved resistance, calibrated uncertainty with coverage, JOSS paper, conda-forge | Oct 2027 |
| 4 — started | **Your own cells**: `cellsim plate` reads a plate-reader export and fits it with a measured error bar; next, fit this line's own simulator constants and give predictions as bands rather than lines ([plan](docs/PLAN.md)) | 2028 |

Details, the keep/kill list and success metrics: [`docs/PLAN.md`](docs/PLAN.md).

## Layout

```
cellsim/
  dock/      Vina + Meeko + PoseBusters + fpocket, batch screen, triage, strain, off-target
  chem/      SMILES → OpenFF Sage + AM1-BCC; ADMET descriptors; profile dashboards
  uq/        Monte-Carlo, Sobol, split-conformal, target-class reliability table
  quantum/   xTB GFN2, CYP3A4 site-of-metabolism (advisory), PySCF hook
  cell/      cell-level PK/PD modules (Phase 1 engine grows here)
  bridge/    binding ΔG → rate-law priors with provenance and accuracy
  cache/     content-addressed SQLite cache for every expensive physics call
  fep/       alchemical ΔG scaffolds (experimental; see VALIDATION.md)
  md/        OpenMM ligand / protein MD helpers
benchmarks/  every input needed to reproduce a number in docs/VALIDATION.md
tests/       one folder per src/ module; each file is a standalone gate
scripts/     the `cellsim` CLI dispatcher, fetchers, benchmark runners, doctor
docs/        plan, validation, docking tutorial, archive of the 2026 campaign docs
OLD/         frozen 2026 C++/Metal prototype, its 6 headless validators and scripts (see OLD/README.md)
```

## Validation and CI

Every push and pull request runs [`smoke.yml`](.github/workflows/smoke.yml):
54 gates across chem, MD, docking, FEP scaffolds, quantum, UQ, bridge,
cell and cache, about 15 minutes on a 2-core runner. A regression in any
gate blocks merge. Numbers that take hours (FEP sampling, scale tests)
run only on `workflow_dispatch` and are reported in
[`docs/VALIDATION.md`](docs/VALIDATION.md) with the commit that produced
them.

## History, briefly

v1.0 (March 2026) was a real-time C++/Metal 3D cell simulator with drug
PK/PD. In April 2026 the project restarted around a non-AI molecular
stack after a critique of circular validation; that stack is the
`src/` tree above. The July 2026 fixes corrected the docking numbers and
the FEP binding path. In October 2026 the project was refocused on the
single product described at the top. The complete record, including the
original roadmap, the professor debrief, and every benchmark post-mortem,
is preserved in [`docs/archive/campaign1/`](docs/archive/campaign1/).

## License

MIT, see [`LICENSE`](LICENSE). Third-party attributions are in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
