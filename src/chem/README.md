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

## Exit test
Round-trip 10 000 ChEMBL compounds → parameterised OpenMM system;
every one produces a valid, energy-minimizable system.

## Viewer (`cellsim-chem --view`)
3D ball-and-stick of the parameterised ligand; atoms coloured by
partial charge.
