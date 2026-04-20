# src/dock — Docking + binding free energy

**Campaign 1, Layer 1.3.** Non-AI.

## Scope
Given a receptor (from `src/md/protein.load_protein_pdb`) and one or
more ligands (from `src/chem/parametrize_smiles`), produce a
`DockingResult` per pair with:

- Top-N 3D poses (seed-pinned, reproducible).
- AutoDock Vina ΔG in both **kcal/mol** (biologist default) and
  **kJ/mol** (SI).
- PoseBusters physical-validity flags (no clashes, sane bond lengths,
  correct chirality).
- Pose RMSD vs a supplied reference pose when available (canonical
  re-docking validation metric).
- Provenance: Vina version, receptor PDB, ligand SMILES + InChI key,
  exhaustiveness, num modes, seed.

The `src/cache/` layer keys on `(ligand_hash, receptor_hash,
exhaustiveness, seed)` so re-docking an identical pair never re-runs
the search.

## Upstream tools (build-vs-buy: BUY, non-AI)
- **AutoDock Vina 1.2+** (Apache-2) — primary docking engine.
  Empirical scoring function, auditable, widely validated.
- **Meeko** (Apache-2) — PDBQT conversion for ligand and receptor
  (the file format Vina consumes).
- **PoseBusters** (BSD-3) — physical-validity test suite (clashes,
  bonds, chirality, ring flattening).
- **RDKit** — SMILES/SDF handling, pose RMSD calculation.
- **PDBFixer** (via `src/md/protein.py`) — receptor preparation.

## Explicitly not used
- GNINA CNN-scored mode as the primary path. Per `MISSION.md` it may
  ship only as an explicitly labeled "fast-guess" option, shown
  alongside the Vina score, never as the sole value.
- Any learned scoring function.

## Modules
- `vina.py` — `dock_ligand(receptor_result, parametrize_result, *,
  center, box_size, exhaustiveness, num_modes, seed) ->
  DockingResult`. Pure Python, subprocess-driven Vina CLI.
- `prep.py` — receptor + ligand → PDBQT via Meeko. Hot path, cached.
- `pose_rmsd.py` — RMSD vs reference pose helper (re-docking gate).
- `validity.py` — PoseBusters wrapper.
- `viewer.py` — matplotlib static view of top pose overlaid on
  receptor Cα ribbon + ΔG bar with Monte-Carlo CI.

## Exit tests
- **Re-docking gate (MVP):** take a well-characterised cocrystal
  (e.g. PDB 1M17 EGFR + erlotinib, or 1STP streptavidin + biotin),
  strip the ligand, re-dock, require top pose RMSD < 2.0 Å vs the
  crystal pose. This is the canonical validation metric every
  docking method in the field is judged by.
- **PDBBind refined gate (full):** ≥ 75 % pose recovery within
  2 Å on the refined set. Deferred to a follow-up PR.

## Biologist-UX notes
Outputs that matter to a biologist (and must be first-class in the
API):
- ΔG in kcal/mol with ± CI. Secondary: kJ/mol.
- K_d estimate (implied from ΔG), reported in nM/µM as a biologist
  would use it.
- Top-1 pose RMSD vs reference (when available). Green if < 2 Å,
  yellow 2–3 Å, red > 3 Å.
- PoseBusters flags per pose: "ready for wet-lab?" tag.
- A single HTML/CSV report the biologist can hand to wet-lab.

## Quickstart (planned; coming in subsequent PRs)
```bash
conda activate cellsim
python -m src.dock.vina \
    --receptor benchmarks/md/1ubq.pdb \
    --ligand-smiles "CC(=O)OC1=CC=CC=C1C(=O)O" \
    --center 20,20,20 --box 20,20,20 \
    --save dock_report.html
```
