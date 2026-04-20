"""Layer 1.4 quantum viewer — xTB charges + dipole + electronic panel.

LEFT-3D  Ligand ball-and-stick, atoms colour-mapped by xtb Mulliken
         partial charge (blue = electron-rich, red = electron-poor).
         Dipole drawn as a black arrow from the charge centroid in
         the dipole's direction (magnitude proportional to Debye).

RIGHT    Compact "electronic datasheet" panel:
           - Molecular formula + SMILES
           - HOMO / LUMO / gap (eV bar chart in a style biologists
             know from UV-vis literature)
           - Dipole magnitude
           - Total energy (eV)
           - Per-atom charge histogram
           - Tool provenance (xtb version, method, solvent)

This is the eye-check a medicinal chemist wants before running an
xtb calculation through a metabolism / ADMET predictor: "are the
charges sensible? is this molecule polar or apolar? where's the
nucleophile?"

Run:
    conda activate cellsim
    python -m src.quantum.viewer aspirin \\
        --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --save aspirin_qm.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.quantum import xtb_single_point, XtbResult  # noqa: E402


_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Pt": 1.36, "Na": 1.66, "K": 2.03, "Ca": 1.76, "Mg": 1.41,
}


def _atom_size(sym: str) -> float:
    return _RADII.get(sym, 0.75) * 240.0


def _charge_color(q: float, vmax: float = 0.6):
    """Blue (q<0) → white (q=0) → red (q>0)."""
    t = max(-1.0, min(1.0, q / vmax))
    if t >= 0:
        return (1.0, 1.0 - 0.85 * t, 1.0 - 0.85 * t)
    return (1.0 + 0.85 * t, 1.0 + 0.85 * t, 1.0)


def render_xtb_result(
    result: XtbResult,
    *,
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
    charges = result.mulliken_charges or [0.0] * len(elements)
    if len(charges) != len(elements):
        # fall back: pad or truncate
        charges = (charges + [0.0] * len(elements))[:len(elements)]

    sizes = np.array([_atom_size(s) for s in elements])
    colors = [_charge_color(q) for q in charges]

    fig = plt.figure(figsize=(13, 6.5))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax_panel = fig.add_subplot(1, 2, 2)

    # Atoms
    ax3d.scatter(pos[:, 0], pos[:, 1], pos[:, 2],
                 s=sizes, c=colors, edgecolors="black",
                 linewidths=0.5, zorder=3)

    # Bonds — guess from distances (no topology here; good enough for
    # organic small molecules where covalent cut is ~1.8 Å).
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            d = float(np.linalg.norm(pos[i] - pos[j]))
            # Sum of covalent radii × 1.25 as the cutoff.
            r_cut = (_RADII.get(elements[i], 0.75)
                     + _RADII.get(elements[j], 0.75)) * 1.25
            if d < r_cut:
                ax3d.plot(
                    [pos[i, 0], pos[j, 0]],
                    [pos[i, 1], pos[j, 1]],
                    [pos[i, 2], pos[j, 2]],
                    color=(0.3, 0.3, 0.3), linewidth=1.0, zorder=1)

    # Atom labels (heavy atoms only)
    for (x, y, z), sym, q in zip(pos, elements, charges):
        if sym != "H":
            ax3d.text(x, y, z + 0.2, f"{sym}\n{q:+.2f}",
                      fontsize=7, ha="center", zorder=4)

    # Dipole arrow: from centroid along (x,y,z) — we only have the
    # magnitude, so use charge-weighted-centroid direction as proxy.
    centroid = pos.mean(axis=0)
    if charges:
        weighted = (pos * np.array(charges).reshape(-1, 1)).sum(axis=0)
        if np.linalg.norm(weighted) > 1e-6:
            direction = weighted / np.linalg.norm(weighted)
            scale = result.dipole_Debye or 1.0
            tip = centroid + direction * min(scale * 0.8, 4.0)
            ax3d.plot([centroid[0], tip[0]],
                      [centroid[1], tip[1]],
                      [centroid[2], tip[2]],
                      color="black", linewidth=2.2, zorder=5)

    mins = pos.min(axis=0)
    maxs = pos.max(axis=0)
    c = (mins + maxs) / 2.0
    span = (maxs - mins).max() * 0.6 + 1.5
    ax3d.set_xlim(c[0]-span, c[0]+span)
    ax3d.set_ylim(c[1]-span, c[1]+span)
    ax3d.set_zlim(c[2]-span, c[2]+span)
    ax3d.set_xlabel("x (Å)")
    ax3d.set_ylabel("y (Å)")
    ax3d.set_zlabel("z (Å)")
    ax3d.set_title("xTB Mulliken charges (blue < 0, red > 0)  +  "
                   "dipole arrow", fontsize=9)

    # Electronic datasheet panel
    ax_panel.axis("off")
    lines = []
    lines.append(f"SMILES:   {result.smiles}")
    lines.append(f"method:   {result.method}"
                 + (f"   solvent: {result.solvent}" if result.solvent
                    else ""))
    if result.total_energy_eV is not None:
        lines.append(f"E_total:  {result.total_energy_eV:>11.2f} eV  "
                     f"({result.total_energy_Hartree:.6f} Eh)")
    if result.homo_eV is not None:
        lines.append(f"HOMO:     {result.homo_eV:>11.2f} eV")
    if result.lumo_eV is not None:
        lines.append(f"LUMO:     {result.lumo_eV:>11.2f} eV")
    if result.homo_lumo_gap_eV is not None:
        lines.append(f"gap:      {result.homo_lumo_gap_eV:>11.2f} eV")
    if result.dipole_Debye is not None:
        lines.append(f"dipole:   {result.dipole_Debye:>11.2f} D")
    lines.append(f"atoms:    {len(elements):>11d}")
    lines.append(f"seed:     {result.random_seed:>11d}")
    lines.append(f"wall:     {(result.wall_seconds or 0.0):>11.2f} s")
    lines.append("")
    lines.append("tool versions:")
    for k, v in (result.tool_versions or {}).items():
        lines.append(f"  {k}: {v}")

    ax_panel.text(0.02, 0.98, "\n".join(lines), fontsize=10,
                  family="monospace", va="top", ha="left")

    # Tiny HOMO/LUMO bar (x-axis = eV).
    if result.homo_eV is not None and result.lumo_eV is not None:
        inset = fig.add_axes([0.55, 0.08, 0.35, 0.18])
        inset.set_xlim(-15, 0)
        inset.barh([0], [abs(result.homo_eV)], left=result.homo_eV,
                   color="#2a5fb4", alpha=0.8,
                   label=f"HOMO {result.homo_eV:.2f} eV")
        inset.barh([0.6], [abs(result.lumo_eV)], left=result.lumo_eV,
                   color="#b43a3a", alpha=0.8,
                   label=f"LUMO {result.lumo_eV:.2f} eV")
        inset.set_yticks([0, 0.6])
        inset.set_yticklabels(["HOMO", "LUMO"], fontsize=8)
        inset.set_xlabel("orbital energy (eV)", fontsize=8)
        inset.axvline(0, color="black", linewidth=0.3)
        inset.set_title(
            f"HL-gap = {(result.homo_lumo_gap_eV or 0.0):.2f} eV",
            fontsize=8)
        inset.spines["right"].set_visible(False)
        inset.spines["top"].set_visible(False)

    fig.suptitle(title or result.smiles, fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    if save:
        plt.savefig(save, dpi=150)
        print(f"[quantum-viewer] wrote {save}")
    if show:
        plt.show()


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("name")
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--solvent", default=None)
    ap.add_argument("--method", default="gfn2")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", help="save PNG instead of showing")
    args = ap.parse_args(argv)

    r = xtb_single_point(
        args.smiles, charge=args.charge, multiplicity=args.mult,
        solvent=args.solvent, method=args.method,
        random_seed=args.seed)
    print(r.summary())
    if not r.ok:
        return 1

    render_xtb_result(r, title=args.name, save=args.save,
                      show=args.save is None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
