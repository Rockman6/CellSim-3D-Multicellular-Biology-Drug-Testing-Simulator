# benchmarks/ — Blind-validation data + harness

## Subdirectories (shipping)

- `chembl/`     — small ChEMBL subset for SMILES round-trip smoke
                  + held-out kinase IC50 panels.
- `dock/`       — cocrystal PDBs + per-system calibration YAMLs
                  (streptavidin, trypsin, EGFR kinase, CYP3A4).
- `md/`         — 1UBQ + loader smoke inputs.
- `fep/`        — FEP benchmark YAMLs:
                  - `freesolv_12.yaml`           — Milestone A hydration (12 cpds)
                  - `binding_streptavidin.yaml`  — Milestone B absolute ΔG_bind (4 cpds)
                  - `binding_egfr.yaml`          — Milestone B EGFR 6-cpd kinase series
                  Each is dry-run validated by `cellsim doctor`
                  step 13 on every PR.
- `quantum/`    — CYP3A4 SoM validation (3 drugs).

## Subdirectories (Campaign 1 Layer 1.8, reserved)

- `pdbbind/`    — PDBBind refined-set pose recovery blind test.
- `casf/`       — CASF-2016 scoring / ranking / docking / screening.
- `posebusters/`— PoseBusters physical-validity suite (already
                  wired via `src/dock/validity.py`).
- `redteam/`    — quarterly adversarial compound drops from external
                  contributors; each one a new regression gate.

## Rule

Every benchmark is **public**, **blind** (the relevant splits must
never have been seen by the model under test), and **runs on every
PR** via GitHub Actions. Pass-rate regressions block merge.
