# src/chem — Cheminformatics + force-field foundation

**Campaign 1, Layer 1.1.**

## Scope
Any SMILES / SDF input → fully parameterised OpenMM system with
AM1-BCC (fast) or xTB (accurate) partial charges.

## Upstream tools (build-vs-buy: BUY)
- **RDKit** (BSD-3) — SMILES/SDF I/O, 3D embedding, conformer generation.
- **OpenFF Toolkit** (MIT) — small-molecule force-field parameter
  assignment (Sage / Parsley).
- **AmberTools** or **RDKit built-in** for AM1-BCC partial charges.
- **PDBFixer** (MIT) — protein protonation, missing atoms.
- **PROPKA** / **DimorphiteDL** — pKa prediction.

## What we write
- Python orchestration glue.
- C++ bindings where RDKit/OpenFF outputs feed the OpenMM driver in
  `src/md/`.

## Modules in this layer
- `parametrize.py` — `parametrize_smiles(smi) → ParametrizeResult`.
  Every result carries method + tool versions + failure reason;
  never raises on recoverable errors.
- `viewer.py` — matplotlib 3D ball-and-stick with per-atom charge
  colouring (blue = negative, red = positive). Stopgap until the
  Metal viewer is rewired in `src/render/` + `src/viewer/`.
- `__init__.py` — package marker re-exporting the public API.

## Quickstart

### One-time env setup
```bash
conda env create -f environment.yml       # or `mamba env create ...`
conda activate cellsim
```

If you want just Layer 1.1 (minimum footprint):
```bash
pip install rdkit openff-toolkit openff-forcefields openmm \
            matplotlib numpy pytest
```

### Parameterise one SMILES
```bash
python -m src.chem.parametrize "CC(=O)OC1=CC=CC=C1C(=O)O" --json
```

### Render a ligand in 3D
```bash
python -m src.chem.viewer aspirin \
    --smiles "CC(=O)OC1=CC=CC=C1C(=O)O"
# or save to PNG (no window needed):
python -m src.chem.viewer aspirin \
    --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --save aspirin.png
```

### Run the smoke test
```bash
pytest tests/chem/test_parametrize_smoke.py -v -s
# or directly:
python tests/chem/test_parametrize_smoke.py
```

Expected: at least 8/10 canonical drugs parameterise cleanly
(`benchmarks/chembl/smoke_10.smi`). The full Layer-1.1 exit
criterion is 10 k ChEMBL compounds with ≥ 99 % success; the
10-compound gate keeps CI runnable in seconds.

## Exit criteria (two-tier gate)
- **RDKit-only tier** (no conda on the box): ≥ 9/10 RDKit-embed
  + 3D coords pass on `smoke_10.smi`. Runs on any dev laptop.
- **Full tier** (OpenFF + AM1-BCC present): ≥ 8/10 end-to-end to
  an OpenMM System on `smoke_10.smi`; ≥ 99 % on 10 k ChEMBL in a
  follow-up PR that adds the dataset.
- Every result carries `charge_method`, `ff_version`, and
  `tool_versions`.

The `test_parametrize_smoke.py` gate auto-detects which tier
applies based on whether `openff.toolkit` imports — no flags
needed.

## Current status on this machine
- ✓ RDKit 2022.09.5 present (via `.molgen_venv/`).
- ✓ OpenMM 8.5.1 present.
- ✗ OpenFF Toolkit absent (PyPI has no wheels for arm64 macOS 26;
  use `conda env create -f environment.yml` to install it).

Until OpenFF lands, `cellsim-chem` runs in RDKit-only mode:
viewer works, charges unavailable.

## Viewer (`cellsim-chem --view`)
3D ball-and-stick of the parameterised ligand; atoms coloured by
partial charge.
