"""Reference-scene viewer for the single-cell drug-disposition layer.

Renders one figure that shows, entirely from the REAL validated model
(`src/cell`), the four behaviours the layer captures:

  A. Dose–response + resistance — target occupancy vs extracellular dose,
     passive vs a P-gp efflux pump (with a Monte-Carlo uncertainty band
     from the K_d prior).
  B. Time course — the full composed transient converging to steady state,
     and the binding sink acting as a capacitor (same endpoint, slower).
  C. pH ion trapping — steady-state accumulation of a weak base vs a weak
     acid across cellular compartments.
  D. Composed waterfall — how occupancy changes as each transport effect
     is layered on (permeation → +trapping → +efflux → +sink).

This is the Layer viewer for criterion 8: it renders the reference scene
(`save=...`) from the validated numeric engine, not a re-implementation.

Palette: Okabe–Ito (published colourblind-safe). Categorical hues are
assigned in fixed order, never cycled.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

# Okabe-Ito, fixed order.
_OI = {
    "blue": "#0072B2", "vermillion": "#D55E00", "green": "#009E73",
    "orange": "#E69F00", "purple": "#CC79A7", "sky": "#56B4E9",
    "grey": "#8A8A8A",
}


def _recessive_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, which="both", color="#DDDDDD", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def render_disposition_scene(
    *,
    dG_kcalmol: float = -9.0,
    receptor: str = "benchmarks/dock/3ptb.pdb",
    save: Optional[str] = None,
    show: bool = True,
    mc_samples: int = 1500,
) -> None:
    import numpy as np
    import matplotlib
    if save is not None and not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.bridge import binding_to_hill
    from src.cell import (
        Permeability, spherical_cell_geometry,
        EffluxPump, BindingSink, COMPARTMENT_PH,
        occupancy_from_prior, hill_occupancy,
        efflux_steady_state, accumulation_ratio,
        solve_steady_state, simulate_disposition,
        monte_carlo_propagate,
    )

    prior = binding_to_hill(dG_kcalmol, uncertainty_kcalmol=0.3,
                            receptor=receptor, method="vina")
    perm = Permeability(P_cm_per_s=1e-6, source="literature")
    geom = spherical_cell_geometry(10.0)
    pump = EffluxPump(Vmax_M_per_s=3e-9, Km_M=1e-7, name="P-gp")

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))
    fig.suptitle(
        f"Single-cell drug disposition — reference scene "
        f"(ΔG={dG_kcalmol:g} kcal/mol, K_d={prior.parameters['Kd_M']:.1e} M)",
        fontsize=13, fontweight="bold")

    # ---- A. Dose-response + resistance + MC band ------------------------
    ax = axes[0, 0]
    doses = np.logspace(-10, -4, 40)
    occ_passive, occ_pgp = [], []
    band_lo, band_hi = [], []
    for d in doses:
        occ_passive.append(hill_occupancy(d, prior.parameters["Kd_M"], 1.0))
        ss = efflux_steady_state(d, perm, geom, pump)
        occ_pgp.append(hill_occupancy(ss.C_in_M, prior.parameters["Kd_M"], 1.0))
        mc = monte_carlo_propagate(
            {"A": prior}, lambda kd: hill_occupancy(d, kd["A"], 1.0),
            n=mc_samples, seed=0)
        band_lo.append(mc.ci95[0])
        band_hi.append(mc.ci95[1])
    occ_passive = np.array(occ_passive) * 100
    occ_pgp = np.array(occ_pgp) * 100
    ax.fill_between(doses, np.array(band_lo) * 100, np.array(band_hi) * 100,
                    color=_OI["blue"], alpha=0.15, lw=0,
                    label="passive 95% CI (K_d unc.)")
    ax.plot(doses, occ_passive, color=_OI["blue"], lw=2, label="passive")
    ax.plot(doses, occ_pgp, color=_OI["vermillion"], lw=2,
            label="with P-gp efflux")
    ax.set_xscale("log")
    ax.set_xlabel("extracellular dose (M)")
    ax.set_ylabel("target occupancy (%)")
    ax.set_title("A · Dose–response and efflux resistance", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    _recessive_axes(ax)

    # ---- B. Time course: sink as capacitor (log time) -------------------
    # The three curves equilibrate on timescales spanning ~10^3, so a log
    # time axis is the only way to show all three; the capacitor effect is
    # the rightward shift with sink capacity (same plateau).
    ax = axes[0, 1]
    for label, sink, col in [
        ("no sink", None, _OI["blue"]),
        ("modest sink", BindingSink(1e-5, 1e-6), _OI["green"]),
        ("large sink", BindingSink(1e-3, 1e-6), _OI["purple"]),
    ]:
        # Each curve gets its OWN horizon (the module's default scales with
        # sink capacity); a shared grid would render the fast curves as a
        # single flat point on a log axis.
        tc = simulate_disposition(1e-7, permeability=perm, geometry=geom,
                                  prior=prior, sink=sink, n_points=400)
        t = np.array(tc.t_s)
        occ = np.array(tc.occupancy) * 100
        ax.plot(t[1:], occ[1:], color=col, lw=2, label=label)
        th = tc.time_to_fraction_of_steady_state(0.5)
        if th:
            ax.plot([th], [np.interp(th, t, occ)], "o", color=col, ms=7,
                    zorder=4, markeredgecolor="white", markeredgewidth=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("time (s, log)   •   dot = t½")
    ax.set_ylabel("target occupancy (%)")
    ax.set_title("B · Transient — binding sink is a capacitor", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    _recessive_axes(ax)

    # ---- C. pH ion trapping across compartments -------------------------
    ax = axes[1, 0]
    comps = ["lysosome", "early_endosome", "cytosol", "mitochondrion"]
    x = np.arange(len(comps))
    base_R = [accumulation_ratio(9.0, COMPARTMENT_PH[c], 7.4, "base")
              for c in comps]
    acid_R = [accumulation_ratio(4.0, COMPARTMENT_PH[c], 7.4, "acid")
              for c in comps]
    w = 0.38
    ax.bar(x - w / 2, base_R, w, color=_OI["green"], label="weak base (pKa 9)",
           zorder=3)
    ax.bar(x + w / 2, acid_R, w, color=_OI["orange"], label="weak acid (pKa 4)",
           zorder=3)
    ax.axhline(1.0, color=_OI["grey"], lw=1, ls="--", zorder=2)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", "\n") for c in comps], fontsize=9)
    ax.set_ylabel("accumulation  C_in / C_out")
    ax.set_title("C · pH ion trapping (vs blood pH 7.4)", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    _recessive_axes(ax)

    # ---- D. Composed waterfall ------------------------------------------
    ax = axes[1, 1]
    C_out = 1e-7
    stages = [
        ("permeation", dict()),
        ("+ lysosomal\ntrapping", dict(pKa=9.0, ion_type="base",
                                       pH_in=COMPARTMENT_PH["lysosome"])),
        ("+ P-gp\nefflux", dict(pKa=9.0, ion_type="base",
                                pH_in=COMPARTMENT_PH["lysosome"], pump=pump)),
    ]
    labels, occs = [], []
    for name, kw in stages:
        ss = solve_steady_state(C_out, permeability=perm, geometry=geom, **kw)
        labels.append(name)
        occs.append(ss.occupancy(prior) * 100)
    xs = np.arange(len(labels))
    bars = ax.bar(xs, occs, 0.6, color=[_OI["blue"], _OI["green"],
                                        _OI["vermillion"]], zorder=3)
    for b, v in zip(bars, occs):
        # Label small values inside-above so a ~0 % bar still reports its
        # number (the efflux collapse is the point of the panel).
        txt = f"{v:.1f}%" if v < 10 else f"{v:.0f}%"
        ax.text(b.get_x() + b.get_width() / 2, max(v, 0) + 1.5, txt,
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("target occupancy (%)")
    ax.set_ylim(0, 112)
    ax.set_title(f"D · Composed steady state (dose {C_out:.0e} M)", loc="left",
                 fontsize=11, fontweight="bold")
    ax.text(0.99, 0.96, "trust: " + (
        "calibrated" if solve_steady_state(
            C_out, permeability=perm, geometry=geom).inputs_calibrated
        else "UNCALIBRATED inputs"),
        transform=ax.transAxes, ha="right", va="top", fontsize=8,
        color=_OI["grey"])
    _recessive_axes(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if save is not None:
        fig.savefig(save, dpi=130, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "cell_disposition_scene.png"
    render_disposition_scene(save=out, show=False)
    print(f"wrote {out}")
