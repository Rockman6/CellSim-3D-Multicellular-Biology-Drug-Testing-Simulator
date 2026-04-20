"""Layer 1.2 viewer — animated ligand trajectory + telemetry panel.

Two panels side-by-side:

  LEFT   3D ball-and-stick animation of the trajectory (frame slider
         + auto-play). Atoms coloured by partial charge if Layer 1.1
         produced them; otherwise element-grey.

  RIGHT  Telemetry panel: temperature, potential energy, RMSD
         time-series with a vertical marker that tracks the
         current frame.

Run:
    python -m src.md.viewer aspirin \\
        --smiles "CC(=O)OC1=CC=CC=C1C(=O)O"
    python -m src.md.viewer aspirin \\
        --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --save aspirin_md.mp4
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.md.simulate import TrajectoryResult, simulate_ligand  # noqa: E402

_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Pt": 1.36, "Na": 1.66, "K": 2.03, "Ca": 1.76, "Mg": 1.41,
}


def _radius(sym: str) -> float:
    return _RADII.get(sym, 0.75)


def render_trajectory(
    result: TrajectoryResult,
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

    if result.n_frames == 0 or not result.positions_A:
        raise SystemExit(f"nothing to render: {result.reason}")

    frames = np.array(result.positions_A)[::stride]
    times = np.array(result.times_ps)[::stride]
    temps = np.array(result.temperatures_K)[::stride]
    energies = np.array(result.potential_energies_kJmol)[::stride]
    rmsds = np.array(result.rmsd_A)[::stride]
    n_frames = len(frames)

    elements = result.elements
    bonds = result.bonds or []
    sizes = np.array([_radius(s) for s in elements]) * 220.0

    fig = plt.figure(figsize=(13, 6))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_tel = fig.add_subplot(1, 2, 2)

    # Fixed world bounds across trajectory so animation doesn't jitter.
    mins = frames.reshape(-1, 3).min(axis=0)
    maxs = frames.reshape(-1, 3).max(axis=0)
    center = (mins + maxs) / 2.0
    span = (maxs - mins).max() * 0.6 + 1.5

    # Scatter placeholder (positions updated per frame)
    scat = ax3d.scatter(
        frames[0, :, 0], frames[0, :, 1], frames[0, :, 2],
        s=sizes, c="lightgray", edgecolors="black", linewidths=0.4)
    bond_lines = []
    for (i, j, order) in bonds:
        lw = 1.0 + 1.2 * (order - 1)
        ln, = ax3d.plot(
            [frames[0, i, 0], frames[0, j, 0]],
            [frames[0, i, 1], frames[0, j, 1]],
            [frames[0, i, 2], frames[0, j, 2]],
            color=(0.3, 0.3, 0.3), linewidth=lw, zorder=1)
        bond_lines.append(ln)

    ax3d.set_xlim(center[0] - span, center[0] + span)
    ax3d.set_ylim(center[1] - span, center[1] + span)
    ax3d.set_zlim(center[2] - span, center[2] + span)
    ax3d.set_xlabel("x (Å)")
    ax3d.set_ylabel("y (Å)")
    ax3d.set_zlabel("z (Å)")

    # Telemetry: three traces on shared-x axes.
    ax_tel.plot(times, temps, color="tab:red", label="T (K)")
    ax_tel.axhline(result.temperature_setpoint_K, color="tab:red",
                   linestyle=":", alpha=0.4, label="T setpoint")
    ax_tel.set_xlabel("time (ps)")
    ax_tel.set_ylabel("T (K)", color="tab:red")
    ax_tel.tick_params(axis="y", labelcolor="tab:red")

    ax_e = ax_tel.twinx()
    ax_e.plot(times, energies, color="tab:blue",
              label="E_pot (kJ/mol)")
    ax_e.set_ylabel("E_pot (kJ/mol)", color="tab:blue")
    ax_e.tick_params(axis="y", labelcolor="tab:blue")

    ax_r = ax_tel.twinx()
    ax_r.spines["right"].set_position(("outward", 55))
    ax_r.plot(times, rmsds, color="tab:green", label="RMSD (Å)")
    ax_r.set_ylabel("RMSD (Å)", color="tab:green")
    ax_r.tick_params(axis="y", labelcolor="tab:green")

    marker = ax_tel.axvline(times[0], color="black", linestyle="--",
                            alpha=0.4)

    def title_text(i: int) -> str:
        return (f"{title or result.smiles}   {result.formula or ''}   "
                f"atoms={result.n_atoms}   "
                f"platform={result.platform}   "
                f"t={times[i]:.2f} ps   "
                f"T={temps[i]:.1f} K   RMSD={rmsds[i]:.3f} Å")

    fig.suptitle(title_text(0), fontsize=10)

    def update(frame_idx: int):
        pts = frames[frame_idx]
        scat._offsets3d = (pts[:, 0], pts[:, 1], pts[:, 2])
        for ln, (i, j, _) in zip(bond_lines, bonds):
            ln.set_data([pts[i, 0], pts[j, 0]], [pts[i, 1], pts[j, 1]])
            ln.set_3d_properties([pts[i, 2], pts[j, 2]])
        marker.set_xdata([times[frame_idx], times[frame_idx]])
        fig.suptitle(title_text(frame_idx), fontsize=10)
        return scat, *bond_lines, marker

    anim = FuncAnimation(fig, update, frames=n_frames,
                         interval=1000 / fps, blit=False)

    if save:
        save_path = Path(save)
        if save_path.suffix.lower() in (".mp4", ".mov"):
            anim.save(save, fps=fps, dpi=110)
        else:
            anim.save(save, fps=fps, dpi=110, writer="pillow")
        print(f"[md-viewer] wrote {save}")
    if show:
        plt.tight_layout()
        plt.show()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("name")
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--dt", type=float, default=2.0, help="fs")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--platform", default=None)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--stride", type=int, default=1,
                    help="render every Nth frame (default 1)")
    ap.add_argument("--save", help="save MP4/GIF instead of showing")
    args = ap.parse_args(argv)

    r = simulate_ligand(
        args.smiles, n_steps=args.steps, dt_fs=args.dt,
        temperature_K=args.temp, platform=args.platform,
        random_seed=args.seed)
    print(r.summary())
    if r.n_frames == 0:
        return 1
    render_trajectory(r, title=args.name, save=args.save,
                      show=args.save is None, stride=args.stride)
    return 0 if r.ok else 1


if __name__ == "__main__":
    sys.exit(main())
