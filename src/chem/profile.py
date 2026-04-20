"""Per-compound drug profile dashboard — one PNG per compound.

Combines Campaign-1 Layer 1.1 + 1.4 output for a single SMILES
into a single biologist-facing dashboard:

  TOP-LEFT    3D ball-and-stick (xTB Mulliken charges colour-coded).
  TOP-MID     SoM predictions (primary SoM = red, top-3 = orange).
  TOP-RIGHT   Electronic + ADMET datasheet (text).
  BOTTOM-LEFT  HOMO/LUMO bar.
  BOTTOM-MID  SoM BDE bar chart (top-10).
  BOTTOM-RIGHT Ro5 / QED / logS callouts.

Usage:
    conda activate cellsim
    python -m src.chem.profile aspirin \\
        --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" \\
        --save profile_aspirin.png

Intended workflow: after the batch screen lands the top-10 hits,
run this for each hit to generate a 1-page biologist summary that
goes into the wet-lab shortlist meeting deck.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.chem import compute_admet, parametrize_smiles  # noqa: E402
from src.quantum import predict_cyp_som_bde, xtb_single_point  # noqa: E402


_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}


def _atom_size(sym: str, scale: float = 200.0) -> float:
    return _RADII.get(sym, 0.75) * scale


def _charge_color(q: float, vmax: float = 0.6):
    t = max(-1.0, min(1.0, q / vmax))
    if t >= 0:
        return (1.0, 1.0 - 0.85 * t, 1.0 - 0.85 * t)
    return (1.0 + 0.85 * t, 1.0 + 0.85 * t, 1.0)


def _bonds_by_distance(pos, elements):
    import numpy as np
    bonds = []
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(pos[i] - pos[j]))
            r_cut = (_RADII.get(elements[i], 0.75)
                     + _RADII.get(elements[j], 0.75)) * 1.25
            if d < r_cut:
                bonds.append((i, j))
    return bonds


def render_profile(
    name: str, smiles: str, *,
    save: Optional[str] = None, show: bool = True,
    seed: int = 1,
) -> int:
    """Render a one-page drug profile dashboard."""
    import numpy as np
    import matplotlib.pyplot as plt

    # --- Run all three analyses.
    print(f"[profile] {name}  {smiles}", flush=True)
    admet = compute_admet(smiles)
    if not admet.ok:
        print(f"[profile] ADMET failed: {admet.reason}")
        return 1
    xtb = xtb_single_point(smiles, random_seed=seed)
    som = predict_cyp_som_bde(smiles, seed=seed)
    parametrise = parametrize_smiles(smiles)

    # --- Figure layout (3-column x 2-row grid + side panels).
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.25)

    # ---- TOP ROW ----
    # (0,0) 3D structure with xTB charges (or AM1-BCC fallback).
    ax3d = fig.add_subplot(gs[0, 0], projection="3d")
    if xtb.ok and xtb.positions_A and xtb.elements:
        pos = np.array(xtb.positions_A)
        elems = xtb.elements
        charges = xtb.mulliken_charges or [0.0] * len(elems)
    elif parametrise.positions_nm and parametrise.elements:
        # fall back to AM1-BCC charges from parametrize
        pos = np.array(parametrise.positions_nm) * 10.0  # nm → Å
        elems = parametrise.elements
        charges = parametrise.partial_charges_e or [0.0] * len(elems)
    else:
        pos = None

    if pos is not None:
        for (i, j) in _bonds_by_distance(pos, elems):
            ax3d.plot([pos[i, 0], pos[j, 0]],
                      [pos[i, 1], pos[j, 1]],
                      [pos[i, 2], pos[j, 2]],
                      color=(0.3, 0.3, 0.3), linewidth=1.0, zorder=1)
        sizes = [_atom_size(s) for s in elems]
        colors = [_charge_color(q) for q in charges]
        ax3d.scatter(pos[:, 0], pos[:, 1], pos[:, 2],
                     s=sizes, c=colors,
                     edgecolors="black", linewidths=0.4, zorder=3)
        mn, mx = pos.min(0), pos.max(0)
        ctr = (mn + mx) / 2.0
        span = (mx - mn).max() * 0.6 + 1.5
        for a, c in zip(
                (ax3d.set_xlim, ax3d.set_ylim, ax3d.set_zlim), ctr):
            a(c - span, c + span)
        ax3d.set_xlabel("x (Å)"); ax3d.set_ylabel("y (Å)")
        ax3d.set_zlabel("z (Å)")
    method = (xtb.method if xtb.ok else
              parametrise.charge_method if parametrise.ok else "none")
    ax3d.set_title(f"3D + charges ({method})", fontsize=10)

    # (0,1) SoM view (primary red / top-3 orange).
    ax_som = fig.add_subplot(gs[0, 1], projection="3d")
    if som.ok and som.positions_A:
        pos_s = np.array(som.positions_A)
        elems_s = som.elements
        hottest: dict[int, tuple[float, int]] = {}
        for c in som.candidates:
            p = c.parent_atom_idx
            if c.bde_kcalmol is None:
                continue
            cur = hottest.get(p)
            if cur is None or c.bde_kcalmol < cur[0]:
                hottest[p] = (c.bde_kcalmol, c.rank or 999)

        def _som_color(i):
            if i not in hottest:
                return (0.78, 0.78, 0.78)
            _, rank = hottest[i]
            if rank == 1:
                return (0.85, 0.15, 0.15)
            if rank <= 3:
                return (0.95, 0.55, 0.20)
            return (0.55, 0.55, 0.55)
        for (i, j) in _bonds_by_distance(pos_s, elems_s):
            ax_som.plot([pos_s[i, 0], pos_s[j, 0]],
                        [pos_s[i, 1], pos_s[j, 1]],
                        [pos_s[i, 2], pos_s[j, 2]],
                        color=(0.3, 0.3, 0.3), linewidth=1.0, zorder=1)
        sizes_s = [_atom_size(s) for s in elems_s]
        colors_s = [_som_color(i) for i in range(len(elems_s))]
        ax_som.scatter(pos_s[:, 0], pos_s[:, 1], pos_s[:, 2],
                       s=sizes_s, c=colors_s,
                       edgecolors="black", linewidths=0.4, zorder=3)
        for c in som.top_k(3):
            if c.parent_atom_idx < len(pos_s):
                x, y, z = pos_s[c.parent_atom_idx]
                ax_som.text(x, y, z + 0.2,
                            f"#{c.rank}",
                            fontsize=8, ha="center",
                            color="red" if c.rank == 1 else "darkorange",
                            weight="bold", zorder=5)
        mn, mx = pos_s.min(0), pos_s.max(0)
        ctr = (mn + mx) / 2.0
        span = (mx - mn).max() * 0.6 + 1.5
        for a, c in zip(
                (ax_som.set_xlim, ax_som.set_ylim, ax_som.set_zlim), ctr):
            a(c - span, c + span)
        ax_som.set_xlabel("x (Å)"); ax_som.set_ylabel("y (Å)")
        ax_som.set_zlabel("z (Å)")
    ax_som.set_title("predicted CYP3A4 site-of-metabolism", fontsize=10)

    # (0,2) text datasheet.
    ax_text = fig.add_subplot(gs[0, 2])
    ax_text.axis("off")
    lines = []
    lines.append(f"  SMILES:  {smiles}")
    if admet.ok:
        lines.append(f"  formula: {admet.formula}   InChIKey: "
                     f"{admet.inchi_key or '-'}")
        lines.append("")
        lines.append("  PHYSICOCHEMICAL (non-AI):")
        lines.append(f"    MW           {admet.MW:>7.2f} Da")
        lines.append(f"    logP         {admet.logP:>+7.2f}")
        lines.append(f"    TPSA         {admet.tpsa:>7.1f} Å²")
        lines.append(f"    HBA / HBD    {admet.hba} / {admet.hbd}")
        lines.append(f"    rotb / aromatic rings   {admet.rotb} / "
                     f"{admet.aromatic_rings}")
        lines.append(f"    heavy atoms  {admet.heavy_atoms}")
        lines.append("")
        lines.append("  DRUG-LIKENESS:")
        ro5 = ("✓  Ro5 pass" if admet.ro5_pass else
               f"✗  {admet.ro5_violations} Ro5 violation(s)")
        lines.append(f"    Lipinski    {ro5}")
        lines.append(f"    QED         {admet.qed:.3f} (Bickerton 2012)")
        lines.append(f"    logS (ESOL) {admet.logS_ESOL:+.2f} → "
                     f"{admet.solubility_class}")
    if xtb.ok:
        lines.append("")
        lines.append("  ELECTRONIC (xTB GFN2):")
        lines.append(f"    HOMO  {xtb.homo_eV:>+7.2f} eV")
        lines.append(f"    LUMO  {xtb.lumo_eV:>+7.2f} eV")
        lines.append(f"    gap   {xtb.homo_lumo_gap_eV:>+7.2f} eV")
        lines.append(f"    dipole {xtb.dipole_Debye:>6.2f} D")
    if som.ok:
        lines.append("")
        lines.append("  CYP3A4 SoM (xTB BDE, top-3):")
        for c in som.top_k(3):
            lines.append(f"    #{c.rank} {c.parent_element}"
                         f"(idx={c.parent_atom_idx})  "
                         f"BDE = {c.bde_kcalmol:.1f} kcal/mol")
    ax_text.text(0.0, 1.0, "\n".join(lines), fontsize=9,
                 family="monospace", va="top", ha="left")

    # ---- BOTTOM ROW ----
    # (1,0) HOMO/LUMO horizontal bar.
    ax_hl = fig.add_subplot(gs[1, 0])
    if xtb.ok and xtb.homo_eV is not None and xtb.lumo_eV is not None:
        ax_hl.barh([0.4], [abs(xtb.homo_eV)], left=xtb.homo_eV,
                   color="#2a5fb4", alpha=0.85)
        ax_hl.barh([1.0], [abs(xtb.lumo_eV)], left=xtb.lumo_eV,
                   color="#b43a3a", alpha=0.85)
        ax_hl.set_yticks([0.4, 1.0])
        ax_hl.set_yticklabels(
            [f"HOMO\n{xtb.homo_eV:+.2f}", f"LUMO\n{xtb.lumo_eV:+.2f}"],
            fontsize=9)
        ax_hl.set_xlim(-15, 0)
        ax_hl.axvline(0, color="black", linewidth=0.3)
        ax_hl.set_xlabel("orbital energy (eV)")
        ax_hl.set_title(
            f"HL-gap {xtb.homo_lumo_gap_eV:.2f} eV", fontsize=9)
    else:
        ax_hl.axis("off")
        ax_hl.text(0.5, 0.5, "xTB unavailable", ha="center")

    # (1,1) SoM BDE bar chart (top-10).
    ax_bde = fig.add_subplot(gs[1, 1])
    if som.ok and som.candidates:
        top = som.candidates[:10]
        bdes = [c.bde_kcalmol for c in top]
        labels = [f"{c.parent_element}(idx={c.parent_atom_idx})" for c in top]
        colors = [
            (0.85, 0.15, 0.15) if c.rank == 1 else
            (0.95, 0.55, 0.20) if c.rank and c.rank <= 3 else
            (0.60, 0.60, 0.60)
            for c in top]
        ranks = list(range(1, len(top) + 1))
        ax_bde.barh(ranks, bdes, color=colors, edgecolor="black",
                    linewidth=0.5)
        ax_bde.set_yticks(ranks)
        ax_bde.set_yticklabels(labels, fontsize=8)
        ax_bde.invert_yaxis()
        ax_bde.set_xlabel("BDE (kcal/mol)")
        ax_bde.set_title("top-10 SoM candidates", fontsize=9)
    else:
        ax_bde.axis("off")
        ax_bde.text(0.5, 0.5, "SoM unavailable", ha="center")

    # (1,2) Ro5 / QED callout.
    ax_call = fig.add_subplot(gs[1, 2])
    ax_call.axis("off")
    if admet.ok:
        # Big coloured badge.
        y = 0.92
        if admet.ro5_pass:
            ax_call.text(0.5, y, "Ro5 ✓", ha="center",
                         fontsize=30, color="#2a9d3a", weight="bold")
        else:
            ax_call.text(0.5, y, f"Ro5 ✗ ×{admet.ro5_violations}",
                         ha="center", fontsize=26, color="#c62828",
                         weight="bold")
        ax_call.text(0.5, 0.76,
                     f"QED = {admet.qed:.2f}",
                     ha="center", fontsize=18,
                     color=("#2a9d3a" if admet.qed >= 0.6
                            else "#e6a82a" if admet.qed >= 0.4
                            else "#c62828"))
        ax_call.text(0.5, 0.60,
                     f"logS = {admet.logS_ESOL:+.2f}",
                     ha="center", fontsize=16)
        ax_call.text(0.5, 0.50, admet.solubility_class, ha="center",
                     fontsize=13, style="italic")
        ax_call.text(0.5, 0.32,
                     f"MW {admet.MW:.0f}  logP {admet.logP:+.2f}",
                     ha="center", fontsize=14)
        ax_call.text(0.5, 0.20,
                     f"TPSA {admet.tpsa:.1f} Å²",
                     ha="center", fontsize=12)
        ax_call.text(0.5, 0.08,
                     f"HBA {admet.hba}  HBD {admet.hbd}  "
                     f"rotb {admet.rotb}",
                     ha="center", fontsize=11)

    fig.suptitle(f"CellSim drug profile  —  {name}   "
                 f"{admet.formula or smiles}",
                 fontsize=13)
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
        print(f"[profile] wrote {save}", flush=True)
    if show:
        plt.show()
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("name")
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", help="save PNG instead of showing")
    args = ap.parse_args(argv)
    return render_profile(args.name, args.smiles, seed=args.seed,
                           save=args.save, show=args.save is None)


if __name__ == "__main__":
    sys.exit(main())
