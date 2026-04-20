# src/quantum — Semi-empirical quantum (Campaign 1, Layer 1.4, non-AI)

xTB GFN2 (and PySCF DFT on demand) for reactive fragments, HOMO /
LUMO, electrostatic properties, and CYP-family site-of-metabolism
prediction. All methods are deterministic and published; no ML.

## Why this layer matters for wet-lab replacement
- HOMO-LUMO gap + Fukui indices predict where a metabolising enzyme
  (CYP3A4 etc.) will oxidise a drug — informs half-life / clearance.
- ESP + dipole inform binding-site complementarity.
- Provides electronic features for reactive-metabolite risk
  assessment (drug-induced liver injury screen).

## Upstream tools (BUY, non-AI)
- **xtb** (LGPL) — Grimme-group semi-empirical QM. Primary engine.
  CLI + optional Python bindings.
- **PySCF** (Apache-2) — DFT for cases xtb can't handle. Hook only;
  bulk work stays in xtb.
- **RDKit** — SMILES → 3D conformer input.

## Modules (Layer 1.4.1 shipped)
- `xtb.py` — `xtb_single_point(smiles, charge, multiplicity,
  solvent, method, random_seed, timeout_s) → XtbResult`.
  Subprocess the `xtb` CLI; parse stdout for total energy (Hartree
  + eV), HOMO/LUMO (eV), HL-gap (eV), molecular dipole (Debye),
  Mulliken charges, and CM5 charges via the xtb JSON output.
  ALPB implicit solvent optional (water, DMSO, …).
  Deterministic: ETKDGv3 conformer seed pinned; xtb dispersion is
  deterministic given the same input geometry.

## Planned (pending)
- `metabolism.py` — CYP3A4 site-of-metabolism prediction via Fukui
  indices (global and local softness). Benchmarked against ~20
  marketed CYP3A4 substrates from the literature.
- `viewer.py` — ESP isosurface + HOMO/LUMO density clouds overlaid
  on the RDKit-embedded molecule.
- `dft.py` — PySCF wrapper for cases xtb fails (transition metals
  beyond xtb's coverage, open-shell, excited states).

## Biologist quickstart
```bash
conda activate cellsim
python -m src.quantum.xtb "CC(=O)OC1=CC=CC=C1C(=O)O"
# [OK]   xtb GFN2-xTB  CC(=O)...  atoms=21  E=-1078.20 eV
#        HOMO=-11.38  LUMO=-7.87  gap=3.51 eV  µ=2.39 D  (0.31 s)
```

## Exit tests
- **Smoke (shipped):** `tests/quantum/test_xtb_smoke.py` — 10
  canonical drugs from `benchmarks/chembl/smoke_10.smi`; sanity
  gates on energy sign, HOMO<LUMO, gap in (0, 12) eV, dipole
  finite. 10/10 pass in < 1 s total.
- **CYP3A4 SoM gate (planned):** match literature CYP3A4 primary
  site-of-metabolism on ≥ 15/20 marketed drugs.

## Non-AI disclaimer
All numbers here come from semi-empirical or ab-initio methods with
known parametrisation. There is no learned model in this layer;
SoM predictions use classical Fukui / softness analysis of the
xtb-derived electron density, not a neural classifier. See
`MISSION.md` for the repository-wide no-AI commitment.
