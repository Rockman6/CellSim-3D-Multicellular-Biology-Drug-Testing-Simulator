"""Layer 1.2 protein viewer — animated Cα trace + telemetry panel.

LEFT   3D Cα-trace of the protein, one line segment per residue pair,
       animated through the trajectory frames. Residues are coloured
       by secondary-structure proxy (B-factor-free substitute: a
       smooth gradient over sequence index, distinguishing the
       chains).

RIGHT  Temperature / potential energy / Cα-RMSD time series with a
       vertical marker tracking the current frame (same layout as
       the ligand viewer).

Run:
    conda activate cellsim
    python -m src.md.protein_viewer benchmarks/md/1ubq.pdb \\
        --md-steps 1000 --save 1ubq_md.gif
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.md.protein import (  # noqa: E402
    ProteinTrajectoryResult,
    load_protein_pdb,
    short_protein_md,
)


def render_protein_trajectory(
    tr: ProteinTrajectoryResult,
    *,
    title: str | None = None,
    save: str | None = None,
    show: bool = True,
    fps: int = 20,
    stride: int = 1,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    if tr.n_frames == 0 or not tr.ca_positions_A_frames:
        raise SystemExit(f"nothing to render: {tr.reason}")

    ca = np.array(tr.ca_positions_A_frames)[::stride]      # (nf, nca, 3)
    times = np.array(tr.times_ps)[::stride]
    temps = np.array(tr.temperatures_K)[::stride]
    energies = np.array(tr.potential_energies_kJmol)[::stride]
    rmsds = np.array(tr.rmsd_ca_A)[::stride]
    n_frames, n_ca, _ = ca.shape

    # Group Cα indices by chain so we don't draw a bond between the
    # C-terminus of one chain and the N-terminus of the next.
    chains_of = [r.get("chain", "A") for r in (tr.ca_residues or [])]
    if len(chains_of) != n_ca:
        chains_of = ["A"] * n_ca
    chain_segments: list[tuple[int, int]] = []
    start = 0
    for i in range(1, n_ca):
        if chains_of[i] != chains_of[start]:
            chain_segments.append((start, i))
            start = i
    chain_segments.append((start, n_ca))

    # Per-residue colour: smooth rainbow along sequence per chain.
    colors = np.zeros((n_ca, 3))
    import matplotlib.cm as cm
    for (a, b) in chain_segments:
        t = np.linspace(0.0, 1.0, max(1, b - a))
        colors[a:b] = cm.viridis(t)[:, :3]

    fig = plt.figure(figsize=(13, 6))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_tel = fig.add_subplot(1, 2, 2)

    # Fixed world bounds across the trajectory.
    flat = ca.reshape(-1, 3)
    mins = flat.min(axis=0)
    maxs = flat.max(axis=0)
    center = (mins + maxs) / 2.0
    span = (maxs - mins).max() * 0.6 + 2.0

    # One Line3D per chain + scatter of Cα atoms.
    lines = []
    for (a, b) in chain_segments:
        ln, = ax3d.plot(
            ca[0, a:b, 0], ca[0, a:b, 1], ca[0, a:b, 2],
            color=colors[(a + b) // 2], linewidth=2.2, zorder=1)
        lines.append((ln, a, b))
    scat = ax3d.scatter(
        ca[0, :, 0], ca[0, :, 1], ca[0, :, 2],
        s=28.0, c=colors, edgecolors="black", linewidths=0.3, zorder=2)

    ax3d.set_xlim(center[0] - span, center[0] + span)
    ax3d.set_ylim(center[1] - span, center[1] + span)
    ax3d.set_zlim(center[2] - span, center[2] + span)
    ax3d.set_xlabel("x (Å)")
    ax3d.set_ylabel("y (Å)")
    ax3d.set_zlabel("z (Å)")

    # Telemetry — same layout as the ligand viewer.
    ax_tel.plot(times, temps, color="tab:red", label="T (K)")
    if tr.temperature_setpoint_K is not None:
        ax_tel.axhline(tr.temperature_setpoint_K, color="tab:red",
                       linestyle=":", alpha=0.4)
    ax_tel.set_xlabel("time (ps)")
    ax_tel.set_ylabel("T (K)", color="tab:red")
    ax_tel.tick_params(axis="y", labelcolor="tab:red")

    ax_e = ax_tel.twinx()
    ax_e.plot(times, energies, color="tab:blue")
    ax_e.set_ylabel("E_pot (kJ/mol)", color="tab:blue")
    ax_e.tick_params(axis="y", labelcolor="tab:blue")

    ax_r = ax_tel.twinx()
    ax_r.spines["right"].set_position(("outward", 55))
    ax_r.plot(times, rmsds, color="tab:green")
    ax_r.set_ylabel("Cα RMSD (Å)", color="tab:green")
    ax_r.tick_params(axis="y", labelcolor="tab:green")

    marker = ax_tel.axvline(times[0], color="black", linestyle="--",
                            alpha=0.4)

    def title_text(i: int) -> str:
        return (f"{title or tr.pdb_id or tr.source}   "
                f"residues={tr.n_residues_protein}   "
                f"atoms={tr.n_atoms}   platform={tr.platform}   "
                f"t={times[i]:.2f} ps   T={temps[i]:.1f} K   "
                f"Cα RMSD={rmsds[i]:.3f} Å")

    fig.suptitle(title_text(0), fontsize=10)

    def update(frame_idx: int):
        pts = ca[frame_idx]
        scat._offsets3d = (pts[:, 0], pts[:, 1], pts[:, 2])
        for ln, a, b in lines:
            ln.set_data(pts[a:b, 0], pts[a:b, 1])
            ln.set_3d_properties(pts[a:b, 2])
        marker.set_xdata([times[frame_idx], times[frame_idx]])
        fig.suptitle(title_text(frame_idx), fontsize=10)
        return scat, *[ln for ln, _, _ in lines], marker

    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000 / fps, blit=False)

    if save:
        save_path = Path(save)
        if save_path.suffix.lower() in (".mp4", ".mov"):
            anim.save(save, fps=fps, dpi=110)
        else:
            anim.save(save, fps=fps, dpi=110, writer="pillow")
        print(f"[protein-viewer] wrote {save}")
    if show:
        plt.tight_layout()
        plt.show()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb")
    ap.add_argument("--md-steps", type=int, default=1000,
                    help="N Langevin steps after minimise (default 1000)")
    ap.add_argument("--dt", type=float, default=2.0, help="fs")
    ap.add_argument("--save-every", type=int, default=50,
                    help="record a frame every N steps (default 50)")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--padding", type=float, default=1.0)
    ap.add_argument("--platform", default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--save", help="save MP4/GIF instead of showing")
    args = ap.parse_args(argv)

    sys_res = load_protein_pdb(
        args.pdb, padding_nm=args.padding, platform=args.platform)
    print(sys_res.summary(), flush=True)
    if not sys_res.ok:
        return 1
    tr = short_protein_md(sys_res, n_steps=args.md_steps, dt_fs=args.dt,
                          save_every=args.save_every,
                          temperature_K=args.temp)
    print(tr.summary(), flush=True)
    if tr.n_frames == 0:
        return 1
    render_protein_trajectory(tr, save=args.save,
                              show=args.save is None,
                              stride=args.stride)
    return 0 if tr.ok else 1


if __name__ == "__main__":
    sys.exit(main())
