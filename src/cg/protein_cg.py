"""Atomistic PDB → Martini 3 coarse-grained topology.

**Scaffold only.** Wraps the `martinize2` / `vermouth` command-line
tool (already installed via environment.yml) to convert an
all-atom protein PDB into Martini 3 bead + elastic-network topology
that can be dropped into a bilayer by `bilayer.py`.

Planned API:

    atomistic_to_martini3(
        pdb_path: Path,
        *,
        ff: str = "martini3001",
        elastic_network: bool = True,
        elastic_force_kjmol: float = 500.0,
        elastic_lower_nm: float = 0.5,
        elastic_upper_nm: float = 0.9,
    ) -> MartiniProtein

MartiniProtein bundles the CG .gro + .itp files so `bilayer.py`
can splice them into a composite topology.

Exit gate (Layer 1.5, protein-in-membrane):
  - 10 µs GPCR (β2-adrenergic receptor, PDB 2RH1) in POPC
  - backbone CG-RMSD vs starting structure stable within 3 Å
  - no helix-bundle unwinding (visual / Cα secondary-structure
    preserved)

Non-AI: Martini force field with published parametrisation,
deterministic CG mapping rules.
"""

from __future__ import annotations


def atomistic_to_martini3(*args, **kwargs):  # pragma: no cover
    raise NotImplementedError(
        "src.cg.protein_cg.atomistic_to_martini3 is planned but "
        "not yet implemented. See src/cg/README.md for scope.")
