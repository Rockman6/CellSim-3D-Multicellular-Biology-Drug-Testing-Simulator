# benchmarks/dock — Docking validation cocrystals

Bundled cocrystal PDBs used by `tests/dock/test_redocking.py` as
canonical pose-recovery benchmarks. Every entry ships with a
published crystal pose so `rmsd_vs_reference_A` can be computed
exactly.

## 1stp.pdb — streptavidin + biotin

The "ubiquitin of docking benchmarks." Small (~120 residues), tight
binder, well-behaved, used in every docking-method validation paper
since AutoDock 1.0.

| Property | Value |
|---|---|
| PDB ID | 1STP |
| Resolution | 2.6 Å |
| Ligand residue name | `BTN` (biotin) |
| Ligand SMILES | `OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12` |
| Crystal ligand centroid (Å) | (11.12, 1.68, −10.75) |
| Recommended Vina search box (Å) | 20 × 20 × 20 around centroid |
| Experimental ΔG | ≈ −18 kcal/mol (K_d ≈ 10⁻¹⁴ M) |

### Vina reality check
When biotin is redocked into 1STP with AutoDock Vina 1.2.5 at
`exhaustiveness=16, num_modes=9, seed=1`, the top pose reports
**ΔG ≈ −7.5 kcal/mol**. This is ~10 kcal/mol weaker than the
experimental affinity — a **known limitation** of Vina's empirical
scoring function for extreme tight binders, not a bug.

Implication for our gates: on 1STP we validate **pose geometry**
(RMSD vs crystal BTN < 2 Å), not absolute ΔG magnitude. The ΔG gate
is reserved for compounds whose experimental affinity falls in
Vina's calibrated range (−6 to −12 kcal/mol, roughly nM to pM).

### Full-range ΔG gate target (next PR)
The canonical ΔG validation is **PDBBind refined set Pearson r**,
where Vina has been extensively benchmarked and is expected to
achieve r ≥ 0.5 against experimental ΔG. That harness lands with
Layer 1.3's PDBBind gate.

## How to add a new cocrystal

1. Fetch the PDB: `python scripts/fetch_pdb.py <ID> --out benchmarks/dock/<id>.pdb`
2. Identify the ligand residue name and compute centroid:
   ```bash
   grep "^HETATM" benchmarks/dock/<id>.pdb | awk '{print $4}' | sort -u
   ```
3. Note the experimental ΔG (from BindingMOAD / PDBBind / the
   original paper) and the Vina-expected ΔG range.
4. Add a section in this README + a case in
   `tests/dock/test_redocking.py`.
