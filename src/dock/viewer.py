"""Layer 1.3 dock viewer — receptor + top pose + crystal overlay.

What a biologist actually wants to look at after a docking run:

  LEFT-3D  Receptor Cα ribbon (sequence-coloured viridis to give
           chain orientation), the top-ranked Vina pose as ligand
           ball-and-stick (carbons green), and — if a crystal
           reference was supplied — the crystal ligand pose as a
           semi-transparent ghost (carbons magenta) for direct
           eye comparison.

  RIGHT    Pose-ranked bar chart: ΔG per pose in kcal/mol (primary,
           biologist default) with the implied K_d annotated on each
           bar. Crystal-RMSD coloured green / amber / red tag per
           pose, so "which poses are trustworthy" is obvious in one
           glance.

This is the Layer-1.3 stopgap using matplotlib — the Metal viewer
lands when src/render/ is wired up.

Run:
    conda activate cellsim
    python -m src.dock.viewer \\
        --receptor benchmarks/dock/1stp.pdb \\
        --ligand-smiles "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12" \\
        --center 11.12,1.68,-10.75 --box 20,20,20 \\
        --crystal-resname BTN \\
        --save dock_1stp.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import (  # noqa: E402
    DockingResult,
    attach_crystal_rmsd,
    dock_ligand,
    extract_hetatm_ligand,
)


def _rmsd_tag_color(r: Optional[float]) -> str:
    if r is None:
        return "#888888"
    if r < 2.0:
        return "#2a9d3a"     # green
    if r < 3.0:
        return "#e6a82a"     # amber
    return "#c62828"         # red


def _extract_ca_trace(pdb_path: Path):
    """Return (positions_A, chains) for protein Cα atoms."""
    import numpy as np

    positions: list = []
    chains: list = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        chain = line[21]
        positions.append([x, y, z])
        chains.append(chain)
    return np.array(positions), chains


def _radii_map(sym: str) -> float:
    return {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
        "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    }.get(sym.title(), 0.75)


def render_dock_result(
    result: DockingResult,
    *,
    receptor_pdb: str | Path,
    crystal_resname: Optional[str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    if not result.ok or not result.poses:
        raise SystemExit(f"cannot render: {result.reason or 'no poses'}")

    receptor_pdb = Path(receptor_pdb)
    ca, chains = _extract_ca_trace(receptor_pdb)
    if len(ca) == 0:
        raise SystemExit(f"no Cα atoms in {receptor_pdb}")

    # Chain segmentation for the Cα line (don't connect C-term to
    # next N-term).
    segments: list[tuple[int, int]] = []
    start = 0
    for i in range(1, len(chains)):
        if chains[i] != chains[start]:
            segments.append((start, i))
            start = i
    segments.append((start, len(chains)))

    # Crystal ligand (if present).
    crystal_heavy = []
    if crystal_resname:
        crystal_heavy = extract_hetatm_ligand(receptor_pdb, crystal_resname)

    # Top pose ligand coords + elements.
    top = result.poses[0]
    pose_xyz = np.array(top.positions_A)
    pose_elems = top.elements

    fig = plt.figure(figsize=(14, 6.5))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_bar = fig.add_subplot(1, 2, 2)

    # Cα trace.
    for a, b in segments:
        t = np.linspace(0.0, 1.0, max(1, b - a))
        colors = cm.viridis(t)[:, :3]
        ax3d.plot(ca[a:b, 0], ca[a:b, 1], ca[a:b, 2],
                  color=colors[(len(t)-1)//2], linewidth=1.6,
                  alpha=0.85, zorder=1)

    # Search-box wireframe (faint).
    if result.center_A and result.box_size_A:
        cx, cy, cz = result.center_A
        dx, dy, dz = result.box_size_A
        corners = np.array([
            [cx-dx/2, cy-dy/2, cz-dz/2], [cx+dx/2, cy-dy/2, cz-dz/2],
            [cx+dx/2, cy+dy/2, cz-dz/2], [cx-dx/2, cy+dy/2, cz-dz/2],
            [cx-dx/2, cy-dy/2, cz+dz/2], [cx+dx/2, cy-dy/2, cz+dz/2],
            [cx+dx/2, cy+dy/2, cz+dz/2], [cx-dx/2, cy+dy/2, cz+dz/2],
        ])
        edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),
                 (0,4),(1,5),(2,6),(3,7)]
        for i, j in edges:
            ax3d.plot(*zip(corners[i], corners[j]),
                      color="#aaaaaa", linestyle="--",
                      linewidth=0.6, alpha=0.5, zorder=0)

    # Crystal ghost (behind) — magenta carbons.
    if crystal_heavy:
        xs = [r["x"] for r in crystal_heavy]
        ys = [r["y"] for r in crystal_heavy]
        zs = [r["z"] for r in crystal_heavy]
        elems = [r["elem"] for r in crystal_heavy]
        sizes = [_radii_map(e) * 250.0 for e in elems]
        colors_xtal = [("#d43caf" if e.title() == "C" else
                        ("#2f7fd4" if e.title() == "N" else
                         "#d43a3a" if e.title() == "O" else
                         "#ddaa00" if e.title() == "S" else "#bbbbbb"))
                       for e in elems]
        ax3d.scatter(xs, ys, zs, s=sizes, c=colors_xtal,
                     alpha=0.35, edgecolors="none", zorder=2,
                     label="crystal")

    # Top pose (in front) — green carbons.
    heavy_idx = [i for i, e in enumerate(pose_elems) if e.upper() != "H"]
    pose_heavy_xyz = pose_xyz[heavy_idx]
    pose_heavy_elems = [pose_elems[i].title() for i in heavy_idx]
    pose_sizes = [_radii_map(e) * 280.0 for e in pose_heavy_elems]
    pose_colors = [("#2a9d3a" if e == "C" else
                    ("#2f7fd4" if e == "N" else
                     "#d43a3a" if e == "O" else
                     "#ddaa00" if e == "S" else "#222222"))
                   for e in pose_heavy_elems]
    ax3d.scatter(pose_heavy_xyz[:, 0], pose_heavy_xyz[:, 1],
                 pose_heavy_xyz[:, 2],
                 s=pose_sizes, c=pose_colors, edgecolors="black",
                 linewidths=0.5, zorder=3, label="top pose")

    # Fit view to the receptor + pose.
    everything = np.vstack([ca, pose_xyz] +
                           ([np.array([[r["x"], r["y"], r["z"]]
                                       for r in crystal_heavy])]
                            if crystal_heavy else []))
    mins = everything.min(axis=0)
    maxs = everything.max(axis=0)
    c = (mins + maxs) / 2.0
    span = (maxs - mins).max() * 0.6 + 2.0
    ax3d.set_xlim(c[0]-span, c[0]+span)
    ax3d.set_ylim(c[1]-span, c[1]+span)
    ax3d.set_zlim(c[2]-span, c[2]+span)
    ax3d.set_xlabel("x (Å)")
    ax3d.set_ylabel("y (Å)")
    ax3d.set_zlabel("z (Å)")
    ax3d.legend(loc="upper left", fontsize=8)

    # Pose bar chart
    poses = result.poses
    xs = list(range(1, len(poses) + 1))
    dGs = [p.affinity_kcalmol for p in poses]
    rmsds = [p.rmsd_vs_reference_A for p in poses]
    bar_colors = [_rmsd_tag_color(r) for r in rmsds]
    bars = ax_bar.bar(xs, dGs, color=bar_colors, edgecolor="black",
                      linewidth=0.6)
    for i, (p, b) in enumerate(zip(poses, bars)):
        h = b.get_height()
        kd = p.kd_implied_nM
        kd_str = (f"{kd:.1f} nM" if kd < 1e3
                  else f"{kd/1e3:.1f} µM" if kd < 1e6
                  else f"{kd/1e6:.1f} mM")
        label = f"{h:.2f}\n{kd_str}"
        if p.rmsd_vs_reference_A is not None:
            label += f"\n{p.rmsd_vs_reference_A:.2f} Å"
        ax_bar.text(i + 1, h - 0.2, label, ha="center", va="top",
                    fontsize=8, color="white"
                    if h < -3.0 else "black")
    ax_bar.set_xticks(xs)
    ax_bar.set_xlabel("pose rank")
    ax_bar.set_ylabel("ΔG (kcal/mol)")
    ax_bar.axhline(0, color="black", linewidth=0.3)
    ax_bar.invert_yaxis()  # more-negative ΔG = better, push up
    ax_bar.set_title(
        f"bar colour: RMSD-vs-crystal "
        f"(green < 2 Å, amber < 3 Å, red ≥ 3 Å, grey = n/a)",
        fontsize=9)

    # Top-strip summary.
    best = poses[0]
    title_lines = [
        (title or f"{result.ligand_formula or 'ligand'}  →  "
         f"{receptor_pdb.name}"),
        f"top-1: ΔG = {best.affinity_kcalmol:.2f} kcal/mol  "
        f"({best.affinity_kJmol:.1f} kJ/mol)  "
        f"K_d ≈ {best.kd_implied_nM:.2g} nM  "
        f"seed={result.seed}  exh={result.exhaustiveness}"
    ]
    if best.rmsd_vs_reference_A is not None:
        title_lines.append(
            f"top-1 crystal RMSD = {best.rmsd_vs_reference_A:.2f} Å "
            f"  |  best-of-top-3 = "
            f"{min(p.rmsd_vs_reference_A for p in poses[:3] if p.rmsd_vs_reference_A is not None):.2f} Å")

    # Strain badge. Skips silently if the diagnostic fails (e.g.
    # RDKit connectivity inference mismatches on this ligand).
    try:
        from src.dock.strain import ligand_strain
        s = ligand_strain(best.elements, best.positions_A,
                           result.ligand_smiles, ensemble_n=20)
        if s.ok:
            title_lines.append(
                f"top-1 strain = {s.strain_kcalmol:+.1f} kcal/mol  "
                f"(UFF ratio {s.energy_ratio:.2f}, {s.band}-band)")
    except Exception:
        pass
    fig.suptitle("\n".join(title_lines), fontsize=10)

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150)
        print(f"[dock-viewer] wrote {save}")
    if show:
        plt.show()


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand-smiles", required=True)
    ap.add_argument("--center", required=True,
                    help="comma-separated x,y,z in Å")
    ap.add_argument("--box", default="22,22,22")
    ap.add_argument("--crystal-resname",
                    help="HETATM residue name to compare against "
                         "(e.g. BTN for 1STP biotin)")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--save", help="save PNG (no window)")
    args = ap.parse_args(argv)

    center = tuple(float(x) for x in args.center.split(","))
    box = tuple(float(x) for x in args.box.split(","))

    r = dock_ligand(
        args.receptor, args.ligand_smiles,
        center_A=center, box_size_A=box,
        exhaustiveness=args.exhaustiveness, num_modes=args.num_modes,
        seed=args.seed, cpu=args.cpu)
    print(r.summary())
    if not r.ok:
        return 1

    if args.crystal_resname:
        r = attach_crystal_rmsd(
            r, crystal_pdb=args.receptor,
            ligand_resname=args.crystal_resname)
        if r.best_rmsd_vs_reference_A is not None:
            print(f"top-1 crystal RMSD = "
                  f"{r.best_rmsd_vs_reference_A:.3f} Å")

    render_dock_result(
        r, receptor_pdb=args.receptor,
        crystal_resname=args.crystal_resname,
        save=args.save, show=args.save is None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
