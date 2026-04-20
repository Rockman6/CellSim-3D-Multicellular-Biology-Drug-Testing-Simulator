"""Layer 1.4 SoM viewer — molecule coloured by predicted CYP3A4 hotspots.

LEFT-3D  Ligand ball-and-stick. Heavy atoms are coloured by the
         BDE of their weakest X-H bond: **red** = lowest BDE (= primary
         site of metabolism), **orange** = top-3, **grey** = not an
         SoM candidate (no X-H or high BDE). A biologist looking at
         the view immediately sees where the drug is likely to be
         oxidised.

RIGHT    Ranked BDE bar chart (top-10 candidates), with per-bar
         annotations showing BDE in kcal/mol and the parent atom
         symbol + index. Green/yellow/grey colour key matches the
         3D panel.

Run:
    conda activate cellsim
    python -m src.quantum.som_viewer testosterone \\
        --smiles "CC12CCC3C(C1CCC2=O)CCC4=CC(=O)CCC34C" \\
        --save som_testosterone.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.quantum import predict_cyp_som_bde, SoMResult  # noqa: E402


_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}


def _atom_size(sym: str) -> float:
    return _RADII.get(sym, 0.75) * 240.0


def render_som_result(
    result: SoMResult,
    *,
    top_k_highlight: int = 3,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    if not result.ok:
        raise SystemExit(f"cannot render: {result.reason}")

    pos = np.array(result.positions_A)
    elements = result.elements
    n_atoms = len(elements)

    # Per-atom "hottest BDE" lookup: for each heavy atom, take the
    # minimum BDE among its X-H candidates.
    hottest: dict[int, tuple[float, int]] = {}
    for c in result.candidates:
        p = c.parent_atom_idx
        if c.bde_kcalmol is None:
            continue
        cur = hottest.get(p)
        if cur is None or c.bde_kcalmol < cur[0]:
            hottest[p] = (c.bde_kcalmol, c.rank or 999)

    # Colour scheme: primary SoM = red, top_k = orange, others =
    # grey (if X-H candidate) / light-grey (not a candidate).
    def _atom_color(i: int):
        if i not in hottest:
            return (0.78, 0.78, 0.78)
        bde, rank = hottest[i]
        if rank == 1:
            return (0.85, 0.15, 0.15)     # red — primary SoM
        if rank <= top_k_highlight:
            return (0.95, 0.55, 0.20)     # orange — top-3
        return (0.50, 0.50, 0.50)         # grey — low-priority X-H

    fig = plt.figure(figsize=(13, 6.5))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_bar = fig.add_subplot(1, 2, 2)

    # Bonds (distance-based for organic small molecules).
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            d = float(np.linalg.norm(pos[i] - pos[j]))
            r_cut = (_RADII.get(elements[i], 0.75)
                     + _RADII.get(elements[j], 0.75)) * 1.25
            if d < r_cut:
                ax3d.plot(
                    [pos[i, 0], pos[j, 0]],
                    [pos[i, 1], pos[j, 1]],
                    [pos[i, 2], pos[j, 2]],
                    color=(0.3, 0.3, 0.3), linewidth=1.0, zorder=1)

    sizes = np.array([_atom_size(s) for s in elements])
    colors = [_atom_color(i) for i in range(n_atoms)]
    ax3d.scatter(pos[:, 0], pos[:, 1], pos[:, 2],
                 s=sizes, c=colors, edgecolors="black",
                 linewidths=0.5, zorder=3)

    # Label top-3 SoM candidates on the 3D view.
    for c in result.top_k(top_k_highlight):
        x, y, z = pos[c.parent_atom_idx]
        ax3d.text(x, y, z + 0.3,
                  f"#{c.rank}\n{c.bde_kcalmol:.1f}",
                  fontsize=8, ha="center",
                  color=("red" if c.rank == 1 else "darkorange"),
                  weight="bold", zorder=5)

    mins = pos.min(axis=0)
    maxs = pos.max(axis=0)
    ctr = (mins + maxs) / 2.0
    span = (maxs - mins).max() * 0.6 + 1.5
    ax3d.set_xlim(ctr[0]-span, ctr[0]+span)
    ax3d.set_ylim(ctr[1]-span, ctr[1]+span)
    ax3d.set_zlim(ctr[2]-span, ctr[2]+span)
    ax3d.set_xlabel("x (Å)")
    ax3d.set_ylabel("y (Å)")
    ax3d.set_zlabel("z (Å)")
    ax3d.set_title(
        "Predicted CYP3A4 SoM  —  red = primary, orange = top-3",
        fontsize=9)

    # BDE bar chart (top-10).
    top = result.candidates[:10]
    ranks = [c.rank for c in top]
    bdes = [c.bde_kcalmol for c in top]
    labels = [f"{c.parent_element}(idx={c.parent_atom_idx})" for c in top]
    bar_colors = []
    for c in top:
        if c.rank == 1:
            bar_colors.append((0.85, 0.15, 0.15))
        elif c.rank and c.rank <= top_k_highlight:
            bar_colors.append((0.95, 0.55, 0.20))
        else:
            bar_colors.append((0.60, 0.60, 0.60))
    ax_bar.barh(ranks, bdes, color=bar_colors, edgecolor="black",
                linewidth=0.5)
    ax_bar.set_yticks(ranks)
    ax_bar.set_yticklabels(labels, fontsize=9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("BDE (kcal/mol)  — lower = easier to abstract")
    ax_bar.set_title(f"Top-10 X-H candidates  "
                     f"({result.method},  {result.enzyme})",
                     fontsize=10)
    # Annotate each bar with the BDE value.
    for r, b in zip(ranks, bdes):
        ax_bar.text(b + 0.3, r, f"{b:.1f}", va="center", fontsize=8)

    fig.suptitle(
        f"{title or result.smiles}  "
        f"({len(result.candidates)} X-H candidates, "
        f"wall {result.wall_seconds:.1f} s)",
        fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    if save:
        plt.savefig(save, dpi=150)
        print(f"[som-viewer] wrote {save}")
    if show:
        plt.show()


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("name")
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--top-k", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", help="save PNG")
    args = ap.parse_args(argv)

    r = predict_cyp_som_bde(args.smiles, seed=args.seed)
    print(r.summary())
    if not r.ok:
        return 1
    render_som_result(r, top_k_highlight=args.top_k,
                      title=args.name, save=args.save,
                      show=args.save is None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
