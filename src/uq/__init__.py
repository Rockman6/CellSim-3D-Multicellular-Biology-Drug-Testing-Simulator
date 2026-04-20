"""CellSim UQ layer (Campaign 1, Layer 1.6, non-AI).

Monte-Carlo + Sobol sensitivity over docking / MD / SoM outputs
using SALib and simple seed ensembles. No deep ensembles; no
learned surrogates.
"""

from .dock_mc import (
    DockingMCResult,
    monte_carlo_dock,
)

__all__ = [
    "DockingMCResult",
    "monte_carlo_dock",
]
