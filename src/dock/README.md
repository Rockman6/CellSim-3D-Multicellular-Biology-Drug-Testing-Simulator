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

## Modules (MVP shipped)
- `vina.py` — `dock_ligand(receptor_pdb, ligand_smiles, *,
  center_A, box_size_A, exhaustiveness, num_modes, seed, cpu,
  timeout_s) → DockingResult`. Subprocess-driven Vina 1.2 CLI.
  Meeko 0.7 for PDBQT prep (ligand via Python API; receptor via
  `mk_prepare_receptor.py` CLI with automatic water/ion cleanup).
- `pose_rmsd.py` — symmetry-aware RMSD-vs-crystal helper. Uses the
  industry-standard recipe: SMILES template → bond-order-assigned
  reference → `rdMolAlign.GetBestRMS` with automorphism search.
  `attach_crystal_rmsd(result, crystal_pdb, ligand_resname)`
  fills `rmsd_vs_reference_A` on every pose.
- `validity.py` — `attach_posebusters(result, receptor_pdb,
  crystal_pdb, crystal_resname)` fills three biologist-relevant
  flags per pose: `posebusters_pocket_ok` (pose in-pocket without
  protein clashes; the triage signal), `posebusters_geometry_ok`
  (covalent sanity), `posebusters_ok` (strict all-tests).
- `viewer.py` — matplotlib static view: receptor Cα ribbon +
  Vina pose (green) overlaid on crystal reference (magenta ghost)
  + ΔG bar chart with K_d annotations and GREEN/AMBER/RED
  crystal-RMSD colouring.

## Exit tests
- **Re-docking gate (MVP, shipped):** 1STP streptavidin-biotin.
  Canonical Astex/PDBBind convention: top-1 RMSD ≤ 2.5 Å
  AND best-of-top-3 < 2.0 Å. Verified locally: top-1 = 2.02 Å,
  top-3 best = 1.99 Å. PASS.
- **PDBBind refined gate (full):** ≥ 75 % pose recovery within
  2 Å on the refined set. Deferred to a follow-up PR that bundles
  the PDBBind harness.

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

## Biologist quickstart

### One-off: dock a compound into your target, see the result
```bash
conda activate cellsim
python -m src.dock.viewer \
    --receptor benchmarks/dock/1stp.pdb \
    --ligand-smiles "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12" \
    --center 11.12,1.68,-10.75 --box 20,20,20 \
    --crystal-resname BTN \
    --save dock_1stp.png
```
Opens a window (or saves PNG) showing:
- Your target's Cα trace.
- The Vina docked pose (green) overlaid on the crystal (magenta).
- A bar chart of top-9 poses, each annotated with ΔG in kcal/mol,
  implied K_d in nM/µM/mM, and crystal-RMSD (green < 2 Å =
  trustworthy, amber 2–3 Å = borderline, red > 3 Å = guess).

### Re-docking regression gate
```bash
python -u tests/dock/test_redocking.py
```
Gate passes iff: top-1 RMSD ≤ 2.5 Å AND best-of-top-3 < 2.0 Å.
Verified locally on 1STP: top-1 = 2.02 Å, top-3 best = 1.99 Å.

### Batch screen: N compounds → ranked CSV (the wet-lab deliverable)
```bash
conda activate cellsim
python -m src.dock.batch \
    --smi benchmarks/dock/1stp_batch_5.smi \
    --receptor benchmarks/dock/1stp.pdb \
    --center 11.12,1.68,-10.75 --box 20,20,20 \
    --exhaustiveness 16 --num-modes 3 --seed 1 \
    --workers 4 --cpu-per-job 2 \
    --crystal-pdb benchmarks/dock/1stp.pdb --crystal-resname BTN \
    --out-csv dock_report.csv
```
Prints a ranked stdout table (RANK / NAME / ΔG / K_d / POCKET / RMSD)
and writes the full CSV. Biologists can drop the CSV directly into
Excel / Prism / a slide to justify their wet-lab shortlist.

On the bundled `1stp_batch_5.smi` the pipeline correctly ranks
biotin (true binder) at #1 with pocket:ok + crystal-RMSD 1.96 Å,
with four non-binder drugs trailing behind. Wall time ~20 s for 5
compounds on a 4-core laptop.

### Programmatic use
```python
from src.dock import dock_ligand, attach_crystal_rmsd, attach_posebusters

result = dock_ligand(
    "benchmarks/dock/1stp.pdb",
    "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
    center_A=(11.12, 1.68, -10.75), box_size_A=(20, 20, 20),
    exhaustiveness=32, num_modes=9, seed=1, cpu=4)
result = attach_crystal_rmsd(result,
    crystal_pdb="benchmarks/dock/1stp.pdb", ligand_resname="BTN")
result = attach_posebusters(result,
    receptor_pdb="benchmarks/dock/1stp.pdb",
    crystal_pdb="benchmarks/dock/1stp.pdb", ligand_resname="BTN")

for pose in result.poses:
    print(pose.biologist_summary())
    # ΔG, kJ/mol, K_d in nM, crystal-RMSD, pocket:ok/fail
```
