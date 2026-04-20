"""Tornado-plot viewer for Sobol sensitivity results.

Classic biology/engineering UX for "which input matters most?" —
horizontal bar chart with inputs sorted by total-effect magnitude,
first-order (S1) and total (ST) indices plotted side-by-side with
Saltelli 95 % CI error bars.

Green = "this knob matters a lot for ΔG variance"
Grey  = "this knob barely moves the needle"
Red   = "this knob has a large total-effect but small first-order
         → its effect is entirely through interactions with others"

Usage:
    from src.uq import sobol_dock
    from src.uq.sobol_viewer import render_sobol_result
    r = sobol_dock(..., n_base=32)
    render_sobol_result(r, save="sobol_1stp.png")
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.uq.sobol import SobolResult  # noqa: E402


def render_sobol_result(
    result: SobolResult,
    *,
    title: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = True,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    if not result.ok:
        raise SystemExit(f"cannot render: {result.reason}")

    names = list(result.input_names)
    s1 = [result.s1.get(n, 0.0) for n in names]
    s1c = [result.s1_conf.get(n, 0.0) for n in names]
    st = [result.st.get(n, 0.0) for n in names]
    stc = [result.st_conf.get(n, 0.0) for n in names]

    # Sort by ST descending (biggest total-effect on top).
    order = np.argsort([-abs(x) for x in st])
    names = [names[i] for i in order]
    s1 = [s1[i] for i in order]
    s1c = [s1c[i] for i in order]
    st = [st[i] for i in order]
    stc = [stc[i] for i in order]

    y = np.arange(len(names))

    fig, (ax_bar, ax_samples) = plt.subplots(
        1, 2, figsize=(13, 5),
        gridspec_kw={"width_ratios": [3, 2]})

    ax_bar.barh(y - 0.2, s1, height=0.4, xerr=s1c,
                color="#2a5fb4", ecolor="#2a5fb4",
                alpha=0.85, label="S1 (first-order)",
                capsize=3, edgecolor="black", linewidth=0.4)
    ax_bar.barh(y + 0.2, st, height=0.4, xerr=stc,
                color="#b43a3a", ecolor="#b43a3a",
                alpha=0.85, label="ST (total-effect)",
                capsize=3, edgecolor="black", linewidth=0.4)
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(names, fontsize=10)
    ax_bar.axvline(0, color="black", linewidth=0.3)
    ax_bar.set_xlabel("Sobol index")
    ax_bar.set_title("Sobol sensitivity of ΔG   (S1 = independent, "
                     "ST = total incl. interactions)",
                     fontsize=10)
    ax_bar.legend(loc="lower right", fontsize=9)
    ax_bar.invert_yaxis()

    # ΔG sample histogram on the right.
    samples = np.array(result.dG_samples) if result.dG_samples else np.array([])
    if samples.size:
        ax_samples.hist(samples, bins=min(20, max(5, len(samples) // 3)),
                         color="#2a9d3a", edgecolor="black", linewidth=0.4,
                         alpha=0.85)
        ax_samples.axvline(samples.mean(), color="black", linestyle="--",
                           linewidth=0.8, label=f"mean {samples.mean():+.2f}")
        ax_samples.axvline(np.median(samples), color="black",
                           linestyle=":", linewidth=0.8,
                           label=f"median {np.median(samples):+.2f}")
        ax_samples.set_xlabel("Vina top-pose ΔG (kcal/mol)")
        ax_samples.set_ylabel("count")
        ax_samples.set_title(
            f"ΔG distribution over {len(samples)} Saltelli samples  "
            f"(range {samples.ptp():.2f} kcal/mol)", fontsize=10)
        ax_samples.legend(fontsize=9)

    fig.suptitle(
        title or (f"Sobol dock sensitivity  —  "
                  f"{Path(result.receptor_pdb).name}  "
                  f"×  {result.ligand_smiles[:30]}" +
                  ("…" if len(result.ligand_smiles) > 30 else "")),
        fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.94))
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
        print(f"[sobol-viewer] wrote {save}", flush=True)
    if show:
        plt.show()
