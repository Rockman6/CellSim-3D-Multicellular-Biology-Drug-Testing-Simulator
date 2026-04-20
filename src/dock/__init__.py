"""CellSim dock layer (Campaign 1, Layer 1.3, non-AI).

AutoDock Vina primary docking engine with Meeko PDBQT prep, PoseBusters
validity, and canonical re-docking pose-RMSD validation.
"""

from .vina import (
    DockingResult,
    DockingPose,
    dock_ligand,
    attach_reference_rmsd,
)
from .pose_rmsd import (
    attach_crystal_rmsd,
    extract_hetatm_ligand,
    pose_rmsd_symmetry_aware,
)
from .validity import attach_posebusters
from .batch import BatchConfig, run_batch
from .pocket_detect import PocketCandidate, detect_pockets
from .refine import refine_pose_openff
from .export import export_poses_sdf, export_poses_pdb
from .off_target import (
    OffTargetEntry,
    OffTargetResult,
    off_target_screen,
)
from .off_target_viewer import render_off_target_result

__all__ = [
    "DockingResult",
    "DockingPose",
    "dock_ligand",
    "attach_reference_rmsd",
    "attach_crystal_rmsd",
    "extract_hetatm_ligand",
    "pose_rmsd_symmetry_aware",
    "attach_posebusters",
    "BatchConfig",
    "run_batch",
    "PocketCandidate",
    "detect_pockets",
    "refine_pose_openff",
    "export_poses_sdf",
    "export_poses_pdb",
    "OffTargetEntry",
    "OffTargetResult",
    "off_target_screen",
    "render_off_target_result",
]
