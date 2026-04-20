"""CellSim MD layer (Campaign 1, Layer 1.2).

OpenMM-backed MD driver. Consumes Layer 1.1's `ParametrizeResult`
and produces a `TrajectoryResult` with temperature / energy / RMSD
telemetry on every stored frame.

First PR: ligand-in-vacuum dynamics. Protein loading (PDBFixer +
ff14SB) and the ubiquitin 100 ns RMSD < 3 Å gate follow.
"""

from .simulate import (
    TrajectoryResult,
    simulate_ligand,
)
from .protein import (
    ProteinSystemResult,
    ProteinTrajectoryResult,
    load_protein_pdb,
    short_protein_md,
)
from .protein_viewer import render_protein_trajectory

__all__ = [
    "TrajectoryResult",
    "simulate_ligand",
    "ProteinSystemResult",
    "ProteinTrajectoryResult",
    "load_protein_pdb",
    "short_protein_md",
    "render_protein_trajectory",
]
