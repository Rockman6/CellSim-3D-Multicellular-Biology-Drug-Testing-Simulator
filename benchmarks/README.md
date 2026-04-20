# benchmarks/ — Blind-validation data + harness

Populated during Campaign 1, Layer 1.8.

## Subdirectories

- `pdbbind/`    — PDBBind refined-set pose recovery blind test.
- `casf/`       — CASF-2016 scoring / ranking / docking / screening.
- `chembl/`     — ChEMBL held-out kinase IC50 panels.
- `posebusters/`— PoseBusters physical-validity suite.
- `redteam/`    — quarterly adversarial compound drops from external
                  contributors; each one a new regression gate.

## Rule

Every benchmark is **public**, **blind** (the relevant splits must
never have been seen by the model under test), and **runs on every
PR** via GitHub Actions. Pass-rate regressions block merge.
