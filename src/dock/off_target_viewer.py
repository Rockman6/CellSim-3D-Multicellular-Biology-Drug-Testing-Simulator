"""Off-target panel viewer — per-receptor ΔG bars.

Biologist UX:
  LEFT: horizontal bar chart of ΔG per receptor, sorted ascending
        (strongest binding at top). Green = best, amber = middle,
        grey = weak. Optional 'intended target' overlay (bold frame)
        if the caller marks one.
  RIGHT: per-receptor K_d annotations (nM/µM/mM auto-formatted) +
         selectivity ΔΔG callout at the bottom.

Usage:
    from src.dock import off_target_screen
    from src.dock.off_target_viewer import render_off_target_result
    r = off_target_screen(smiles, receptors, …)
    render_off_target_result(r, intended_target="streptavidin",
                              save="off_target.png")
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock.off_target import OffTargetResult  # noqa: E402


def _kd_label(kd_nM: Optional[float]) -> str:
    if kd_nM is None:
        return "-"
    if kd_nM < 1.0:
        return f"{kd_nM * 1000:.1f} pM"
    if kd_nM < 1e3:
        return f"{kd_nM:.1f} nM"
    if kd_nM < 1e6:
        return f"{kd_nM / 1e3:.1f} µM"
    return f"{kd_nM / 1e6:.1f} mM"


def render_off_target_result(
    result: OffTargetResult,
    *,
    intended_target: Optional[str] = None,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
) -> None:
    import matplotlib.pyplot as plt

    ranked = [e for e in result.sorted_by_affinity() if e.ok
              and e.dG_kcalmol is not None]
    failed = [e for e in result.entries
              if not (e.ok and e.dG_kcalmol is not None)]
    if not ranked:
        raise SystemExit("no successful dockings to render")

    fig, ax = plt.subplots(figsize=(10, 0.6 + 0.5 * len(ranked)))

    names = [e.name for e in ranked]
    dG = [e.dG_kcalmol for e in ranked]

    # Colour per bar by rank:
    #   #1 strongest → green
    #   middle → amber
    #   weakest → grey
    #   intended target (if named and present) → bold border
    n = len(ranked)
    colors = []
    for i in range(n):
        frac = i / max(1, n - 1)
        if frac < 0.34:
            colors.append("#2a9d3a")       # green
        elif frac < 0.67:
            colors.append("#e6a82a")       # amber
        else:
            colors.append("#888888")       # grey

    y = list(range(n))
    bars = ax.barh(y, dG, color=colors, edgecolor="black",
                    linewidth=0.6)

    if intended_target:
        for i, e in enumerate(ranked):
            if e.name == intended_target:
                bars[i].set_edgecolor("#1b4f87")
                bars[i].set_linewidth(2.5)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=10)
    ax.invert_yaxis()             # rank-1 at top
    ax.set_xlabel("ΔG (kcal/mol)   —   more negative = tighter binding")
    # Room for K_d annotations on the right.
    xmin = min(dG) - 2.5
    xmax = 0.0
    ax.set_xlim(xmin, xmax)

    # Annotate each bar with ΔG + K_d + drug_score
    for i, e in enumerate(ranked):
        kd = _kd_label(e.kd_implied_nM)
        drug = (f"  drug={e.pocket_drug_score:.2f}"
                if e.pocket_drug_score is not None else "")
        ax.text(e.dG_kcalmol + 0.1, i,
                f"{e.dG_kcalmol:+.2f}  |  K_d ≈ {kd}{drug}",
                va="center", fontsize=9, color="#1c1c1c")

    # Selectivity annotation
    sel = result.selectivity_kcalmol()
    lines = []
    if sel is not None:
        tag = ("GOOD" if sel >= 3.0
               else "WEAK" if sel >= 1.0
               else "FLAT")
        lines.append(f"selectivity ΔΔG (best vs 2nd-best): "
                     f"{sel:+.2f} kcal/mol  [{tag}]")
    if intended_target:
        lines.append(f"intended target highlighted: {intended_target}")
    if failed:
        lines.append("FAILED: " + ", ".join(
            f"{e.name} ({(e.reason or '?')[:30]})" for e in failed))

    if lines:
        ax.text(0.02, 0.02, "\n".join(lines),
                transform=ax.transAxes,
                fontsize=9, family="monospace",
                bbox=dict(facecolor="white",
                           edgecolor=(0.5, 0.5, 0.5), alpha=0.9,
                           pad=5),
                va="bottom", ha="left")

    fig.suptitle(
        title or (f"Off-target panel  —  "
                  f"{result.ligand_formula or result.ligand_smiles}"),
        fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
        print(f"[off-target-viewer] wrote {save}")
    if show:
        plt.show()
