"""Layer 1.1 viewer — ball-and-stick of a parameterised ligand.

Minimal matplotlib 3D visualisation. Atom spheres are scaled to
covalent radius (roughly), atoms are coloured by AM1-BCC partial
charge (blue = negative, red = positive, white = neutral), and
bonds draw as line segments with thickness ∝ bond order.

This is the stopgap until src/render/ + src/viewer/ wire up the
Metal pipeline from the OLD/ tree. The point is to have a visual
gate NOW so Layer 1.1's output can be eye-checked.

Run:
    python -m src.chem.viewer "aspirin" \\
        --smiles "CC(=O)OC1=CC=CC=C1C(=O)O"
    python -m src.chem.viewer aspirin \\
        --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --save out.png
"""

from __future__ import annotations

import argparse
import sys

from .parametrize import ParametrizeResult, parametrize_smiles


# Pauling-ish covalent radii in Å (only used for sphere size; approx).
_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Pt": 1.36, "Na": 1.66, "K": 2.03, "Ca": 1.76, "Mg": 1.41,
    "Fe": 1.32, "Zn": 1.22, "Cu": 1.32,
}

_DEFAULT_RADIUS = 0.75


def _element_radius(sym: str) -> float:
    return _RADII.get(sym, _DEFAULT_RADIUS)


def _charge_color(q: float, vmax: float = 0.6):
    """Blue (negative) -> white (0) -> red (positive)."""
    import numpy as np

    t = max(-1.0, min(1.0, q / vmax))
    if t >= 0:
        return (1.0, 1.0 - 0.8 * t, 1.0 - 0.8 * t)
    return (1.0 + 0.8 * t, 1.0 + 0.8 * t, 1.0)


def render_ligand(
    result: ParametrizeResult,
    *,
    title: str | None = None,
    save: str | None = None,
    show: bool = True,
) -> None:
    """Render a parameterised ligand as 3D ball-and-stick."""
    try:
        import matplotlib.pyplot as plt  # noqa: F401
        import numpy as np
    except ImportError as e:
        raise SystemExit(
            "matplotlib + numpy required for viewer: "
            f"{e}\npip install matplotlib numpy")

    if not result.positions_nm or not result.elements:
        raise SystemExit(
            f"nothing to render — parametrise failed: {result.reason}")

    pos_nm = np.array(result.positions_nm, dtype=float)
    pos_A = pos_nm * 10.0  # nm -> Å for display
    elements = result.elements
    bonds = result.bonds or []
    charges = result.partial_charges_e

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Bonds
    for i, j, order in bonds:
        xi, yi, zi = pos_A[i]
        xj, yj, zj = pos_A[j]
        lw = 1.2 + 1.4 * (order - 1)
        ax.plot([xi, xj], [yi, yj], [zi, zj],
                color=(0.3, 0.3, 0.3), linewidth=lw, zorder=1)

    # Atoms
    sizes = np.array([_element_radius(s) for s in elements]) * 200.0
    if charges is not None:
        colors = [_charge_color(q) for q in charges]
    else:
        colors = ["lightgray"] * len(elements)

    ax.scatter(pos_A[:, 0], pos_A[:, 1], pos_A[:, 2],
               s=sizes, c=colors, edgecolors="black",
               linewidths=0.6, zorder=2)

    # Element labels for non-H
    for (x, y, z), sym in zip(pos_A, elements):
        if sym != "H":
            ax.text(x, y, z + 0.2, sym, fontsize=8, zorder=3)

    # Equal aspect
    mins = pos_A.min(axis=0)
    maxs = pos_A.max(axis=0)
    center = (mins + maxs) / 2.0
    span = (maxs - mins).max() * 0.6 + 1.5
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")

    charge_note = (f"charges: {result.charge_method}"
                   if result.charge_method else
                   f"charges: unavailable ({result.reason})")
    ff_note = (f"ff: {result.ff_version}"
               if result.ff_version else "ff: n/a (rdkit-only)")
    ax.set_title(
        f"{title or result.smiles}\n"
        f"{result.formula or ''}   "
        f"{result.n_atoms or 0} atoms   {charge_note}   {ff_note}",
        fontsize=10)

    if save:
        plt.tight_layout()
        plt.savefig(save, dpi=150)
        print(f"[viewer] wrote {save}")
    if show:
        plt.show()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("name", help="display name")
    ap.add_argument("--smiles", required=True, help="SMILES input")
    ap.add_argument("--charge", default="am1bcc")
    ap.add_argument("--ff", default="openff-2.1.0.offxml")
    ap.add_argument("--save", help="save PNG instead of showing")
    args = ap.parse_args(argv)

    result = parametrize_smiles(
        args.smiles, charge_method=args.charge, ff_name=args.ff)

    tag = "OK" if result.ok else "PARTIAL"
    print(f"[{tag}] {args.name}  {args.smiles}")
    if result.ok:
        print(f"  formula={result.formula}  atoms={result.n_atoms}  "
              f"ff={result.ff_version}  charges={result.charge_method}")
        print(f"  hash={result.hash_key()}")
    else:
        print(f"  reason: {result.reason}")

    render_ligand(result, title=args.name, save=args.save,
                  show=args.save is None)
    return 0 if result.ok else 1


if __name__ == "__main__":
    sys.exit(main())
