# src/cache — Physics-prior cache

**Campaign 1, cross-cutting.**

## Scope
Memoise everything expensive. Every docking pose, FEP ΔG, xTB
energy, MACE force, partial-charge table is written here keyed by
content hashes and force-field / method versions. Invalidates on FF
bump.

No existing open-source tool provides this exact schema, so this is
a BUILD item — but it's thin glue over mature storage.

## Upstream storage
- **SQLite** — metadata + query index.
- **HDF5** (hdf5, h5py) — trajectories, force tensors, ESP grids.
- **Parquet** (pyarrow) — tabular ΔG tables for downstream analysis.

## Schema sketch
```
cache.sqlite
  compounds(hash, smiles, inchi, formula, added_at, …)
  receptors(hash, pdb_id, chains, resolution, added_at, …)
  poses(ligand_hash, receptor_hash, method, score, rmsd, …)
  deltaG(ligand_hash, receptor_hash, method, value, ci_lo, ci_hi, …)
  xtb_energies(molecule_hash, method, energy, homo, lumo, …)
  ff_version(tool, version, hash) — invalidates derived rows on bump
```

## Exit test
Cache hit on re-running an already-computed PDBBind entry: ≥ 1000×
speed-up vs cold run. Cache invalidation on FF version bump drops all
affected rows and re-computes.
