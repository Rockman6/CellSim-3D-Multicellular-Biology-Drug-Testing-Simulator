"""Calibration scatter viewer — predicted ΔG vs experimental ΔG.

Canonical biologist plot: each calibration compound is a labelled
point at (predicted, experimental); the y = x line is the perfect-
calibration reference; the ±q95 band around y = x shows CellSim's
honest uncertainty interval. Pearson / Spearman / MAE / RMSE
annotated in the corner.

Reading:
  - Points on the diagonal → CellSim predicts the experimental
    affinity accurately.
  - Points above the diagonal → CellSim under-estimates binding
    (Vina's classic behaviour on tight binders).
  - Points inside the q95 band → CellSim's conformal CI covers the
    truth.
  - Spearman = 1 → ranking is perfect (even if absolute values are
    off): CellSim is usable for triage.

Usage:
    from src.uq import run_calibration, render_calibration_result
    r = run_calibration("benchmarks/dock/streptavidin_calibration.yaml")
    render_calibration_result(r, save="streptavidin_calib.png")
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.uq.calibration import CalibrationResult  # noqa: E402


def render_calibration_result(
    result: CalibrationResult,
    *,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    if not result.ok:
        raise SystemExit(f"cannot render: {result.reason}")

    ok_pts = [p for p in result.points if p.dG_pred_kcalmol is not None]
    if len(ok_pts) < 2:
        raise SystemExit("need >= 2 docked points to plot")

    preds = np.array([p.dG_pred_kcalmol for p in ok_pts])
    truths = np.array([p.dG_expt_kcalmol for p in ok_pts])
    names = [p.name for p in ok_pts]

    lo = float(min(preds.min(), truths.min())) - 2.0
    hi = float(max(preds.max(), truths.max())) + 2.0

    fig, ax = plt.subplots(figsize=(8, 7))

    # Conformal band
    if result.conformal_q95_kcalmol is not None:
        q = result.conformal_q95_kcalmol
        xs = np.array([lo, hi])
        ax.fill_between(xs, xs - q, xs + q,
                        color="#e6eef9", alpha=0.9,
                        label=f"±q95 = ±{q:.2f} kcal/mol")

    # y = x reference
    ax.plot([lo, hi], [lo, hi],
            color=(0.45, 0.45, 0.45), linestyle="--",
            linewidth=1.0, label="y = x  (perfect)")

    # Points + labels
    ax.scatter(preds, truths, s=120, c="#2a9d3a",
               edgecolors="black", linewidths=0.7, zorder=5)
    for x, y, n in zip(preds, truths, names):
        ax.annotate(n, (x, y),
                    textcoords="offset points", xytext=(8, 4),
                    fontsize=10)

    ax.set_xlabel("CellSim predicted ΔG (kcal/mol)")
    ax.set_ylabel("experimental ΔG (kcal/mol)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", fontsize=9)

    # Stats annotation in top-left
    lines = []
    if result.pearson_r is not None:
        lines.append(f"Pearson r  = {result.pearson_r:+.3f}")
    if result.spearman_rho is not None:
        lines.append(f"Spearman ρ = {result.spearman_rho:+.3f}")
    if result.mae_kcalmol is not None:
        lines.append(f"MAE        = {result.mae_kcalmol:.2f} kcal/mol")
    if result.rmse_kcalmol is not None:
        lines.append(f"RMSE       = {result.rmse_kcalmol:.2f} kcal/mol")
    lines.append(f"n = {result.n_ok}/{result.n_points}")
    ax.text(0.03, 0.97, "\n".join(lines),
            transform=ax.transAxes,
            va="top", ha="left", fontsize=10,
            family="monospace",
            bbox=dict(facecolor="white", edgecolor=(0.5, 0.5, 0.5),
                       alpha=0.9, pad=6))

    fig.suptitle(
        title or f"CellSim calibration  —  {Path(result.receptor_pdb).name}",
        fontsize=12)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
        print(f"[calib-viewer] wrote {save}")
    if show:
        plt.show()
