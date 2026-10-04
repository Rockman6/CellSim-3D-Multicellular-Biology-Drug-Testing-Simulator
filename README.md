# CellSim

**A physics-first, uncertainty-honest in-silico cell and drug-response
simulator in Python. Open source, Apple-silicon and Linux, no GPU
required.**

CellSim is converging on one product: an interactive dish. Pick a cell
line and a drug, set a dose schedule, watch the colony respond, and get a
dose-response curve whose error bar has been checked against public
cell-line data. The molecular layer (docking, ADMET, quantum descriptors)
feeds that engine as an input with provenance and a measured accuracy
flag; it is not the product on its own.

> **Status (October 2026).** Phases 0 and 1 of the [plan](docs/PLAN.md)
> are closed. The July fixes are merged, the smoke suite is green, the
> repository has one identity, and `cellsim/cell/` holds a validated
> single-cell drug-response engine.
>
> **Phase 1 closed with a negative result that set the direction.** The
> engine reproduces each drug's potency scale, but predicting *which
> cell line* is more sensitive turned out to be unreachable: eleven
> canonical markers, fitted optimally across 156–683 lines, explain only
> 9–20 % of the difference between lines. So per-line ranking is out of
> scope, and Phase 2 aims at what mechanism is actually good for —
> **within-line dynamics**: dose schedule, wash-out, resistance under
> selection, combination timing. The reasoning and numbers are in
> [`docs/VALIDATION.md`](docs/VALIDATION.md).
>
> **Phase 2 has its dish.** `cellsim dish` grows the engine's cells in
> space — a monolayer with contact inhibition, or a spheroid with
> oxygen and drug diffusing in, hypoxic quiescence and a necrotic core —
> and [`web/viewer/`](web/viewer/) renders the per-cell stream in a
> browser. Every Phase-2 claim was re-measured over seeds before the dish
> was built on it, and one was retracted.
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

In Python:

```python
import numpy as np
from cellsim.cell.engine import dose_response, ic50
from cellsim.cell.library import get_drug, get_line

line, drug = get_line("A549"), get_drug("paclitaxel")
print(ic50(line, drug, guess_uM=0.03))              # 72 h continuous exposure
print(ic50(line, drug, guess_uM=0.03, exposure_h=6))  # 6 h pulse, read at 72 h
```

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
| 1 | Single-cell drug-response engine with one state object and one integrator: cell cycle and p53 axis ported from `OLD/`, death that actually executes, `cellsim/cell` PK on top. Calibrated against GDSC dose-response for cisplatin, doxorubicin and paclitaxel. Headless state stream for any UI. | Jan 2027 |
| 2 | The dish: spatial colony, diffusion, drug penetration, contact inhibition, validated on spheroid growth curves; new web UI reading the engine stream; first outside users | Apr 2027 |
| 3 | Ten drugs × five lines, combinations and resistance, calibrated uncertainty, JOSS paper, conda-forge | Oct 2027 |

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
