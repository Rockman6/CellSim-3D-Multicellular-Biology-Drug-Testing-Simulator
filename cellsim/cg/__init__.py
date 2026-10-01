"""CellSim coarse-grained layer (Campaign 1, Layer 1.5, non-AI).

Martini 3 lipid bilayer + protein elastic-network MD via OpenMM,
for membrane-embedded drug targets (GPCRs, ion channels, kinases
with membrane anchors).

**Scaffold only.** The actual bilayer builder, atomistic→Martini
converter, and long-run MD driver are pending. Planned API lives
in `bilayer.py` and `protein_cg.py`.
"""

from .bilayer import build_martini3_bilayer, run_cg_md
from .protein_cg import atomistic_to_martini3

__all__ = [
    "build_martini3_bilayer",
    "run_cg_md",
    "atomistic_to_martini3",
]
