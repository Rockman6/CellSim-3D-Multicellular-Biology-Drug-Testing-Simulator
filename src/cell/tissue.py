"""Spatial drug penetration into tissue — cells at depth see less drug.

Every module so far describes ONE cell in a well-stirred bath. Real tissue
is not well-stirred: drug leaves a blood vessel, diffuses through the
interstitium, and is consumed by the cells it passes. Cells deep in the
tissue therefore see a much lower concentration than cells next to the
vessel — the reason drugs that work in a dish fail against solid tumours.

Steady-state reaction-diffusion in a slab of tissue (x = 0 at the vessel
wall, x = L at the midpoint between vessels, where symmetry gives zero
flux):

    D·d²C/dx²  =  uptake(C)
    C(0) = C_vessel        dC/dx(L) = 0

For FIRST-ORDER uptake (uptake = k·C) this has the closed form

    C(x) = C_vessel · cosh((L − x)/λ) / cosh(L/λ),      λ = √(D/k)

where λ is the PENETRATION DEPTH: the length scale over which the drug is
consumed. If λ ≪ L, the far half of the tissue is effectively untreated no
matter how much drug is in the blood.

For SATURABLE (Michaelis-Menten) uptake the equation is nonlinear and is
solved numerically. Rigor anchor (tests): at C ≪ K_m the saturable solver
must reproduce the analytic first-order profile with k = V_max/K_m.

Provenance: D and the uptake kinetics are inputs carrying source + trust.
Typical small-molecule interstitial diffusivity is ~1e-6 cm²/s (Nugent &
Jain 1984, Cancer Res 44:238); intercapillary half-distance in tumours is
~50-200 µm (Krogh-cylinder geometry; Thomlinson & Gray 1955, Br J Cancer
9:539 — the classic 150 µm oxygen-diffusion limit).

Non-AI: closed-form + a standard BVP solver. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

# Interstitial diffusivity of a typical small molecule (cm²/s).
# Nugent & Jain 1984, Cancer Res 44:238.
_DEFAULT_D_CM2_PER_S = 1.0e-6


def penetration_depth_um(D_cm2_per_s: float, k_per_s: float) -> float:
    """λ = √(D/k), the depth over which the drug is consumed (µm)."""
    if D_cm2_per_s <= 0 or k_per_s <= 0:
        raise ValueError("D and k must be > 0")
    lam_cm = math.sqrt(D_cm2_per_s / k_per_s)
    return lam_cm * 1e4          # cm → µm


@dataclass
class TissueProfile:
    """Steady-state drug concentration vs depth into tissue."""

    x_um: List[float]
    C_M: List[float]
    C_vessel_M: float
    thickness_um: float
    penetration_depth_um: Optional[float]   # None for saturable uptake
    method: str                              # 'analytic' | 'numerical'
    trust: str = "uncalibrated"

    @property
    def C_far_M(self) -> float:
        """Concentration at the far edge (the least-treated cells)."""
        return self.C_M[-1]

    @property
    def far_fraction(self) -> float:
        """C(far) / C(vessel) — how much drug reaches the deepest cells."""
        return self.C_far_M / self.C_vessel_M if self.C_vessel_M > 0 else 0.0

    def occupancy_profile(self, prior) -> List[float]:
        """Target occupancy of the cells at each depth."""
        from src.cell.occupancy import hill_occupancy
        Kd = float(prior.parameters["Kd_M"])
        n = float(prior.parameters.get("n_hill", 1.0))
        return [hill_occupancy(c, Kd, n) for c in self.C_M]

    def treated_fraction(self, prior, threshold: float = 0.5) -> float:
        """Fraction of the tissue depth whose occupancy exceeds `threshold`.

        The therapeutically meaningful number: a drug that saturates the
        first 20 µm and nothing beyond has 'worked' on 10 % of the tissue.
        """
        occ = self.occupancy_profile(prior)
        n_ok = sum(1 for o in occ if o >= threshold)
        return n_ok / len(occ) if occ else 0.0

    def summary(self) -> str:
        lam = (f"λ={self.penetration_depth_um:.0f}µm"
               if self.penetration_depth_um is not None else "λ n/a (saturable)")
        return (f"{lam}, tissue {self.thickness_um:.0f}µm: "
                f"far-edge {100*self.far_fraction:.1f}% of vessel conc "
                f"[{self.method}, {self.trust}]")


def penetration_profile_first_order(
    C_vessel_M: float, *,
    thickness_um: float,
    k_per_s: float,
    D_cm2_per_s: float = _DEFAULT_D_CM2_PER_S,
    n_points: int = 200,
    trust: str = "uncalibrated",
) -> TissueProfile:
    """Analytic steady-state profile for first-order uptake (uptake = k·C)."""
    if C_vessel_M < 0 or thickness_um <= 0:
        raise ValueError("C_vessel_M ≥ 0 and thickness_um > 0 required")
    lam_um = penetration_depth_um(D_cm2_per_s, k_per_s)
    L = thickness_um
    xs = [L * i / (n_points - 1) for i in range(n_points)]
    denom = math.cosh(L / lam_um)
    cs = [C_vessel_M * math.cosh((L - x) / lam_um) / denom for x in xs]
    return TissueProfile(x_um=xs, C_M=cs, C_vessel_M=C_vessel_M,
                         thickness_um=L, penetration_depth_um=lam_um,
                         method="analytic", trust=trust)


def penetration_profile_saturable(
    C_vessel_M: float, *,
    thickness_um: float,
    Vmax_M_per_s: float,
    Km_M: float,
    D_cm2_per_s: float = _DEFAULT_D_CM2_PER_S,
    n_points: int = 200,
    trust: str = "uncalibrated",
) -> TissueProfile:
    """Numerical steady-state profile for Michaelis-Menten uptake.

    Solves D·C'' = V_max·C/(K_m + C) with C(0) = C_vessel, C'(L) = 0.
    Reduces to the first-order profile when C ≪ K_m (k = V_max/K_m).
    """
    import numpy as np
    from scipy.integrate import solve_bvp

    if C_vessel_M < 0 or thickness_um <= 0:
        raise ValueError("C_vessel_M ≥ 0 and thickness_um > 0 required")
    if Vmax_M_per_s < 0 or Km_M <= 0:
        raise ValueError("Vmax ≥ 0 and Km > 0 required")

    # Work in µm so the geometry is well-scaled: D in µm²/s.
    D_um2 = D_cm2_per_s * 1e8
    L = thickness_um

    def rhs(x, y):
        C = np.maximum(y[0], 0.0)
        return np.vstack([y[1], (Vmax_M_per_s * C / (Km_M + C)) / D_um2])

    def bc(ya, yb):
        return np.array([ya[0] - C_vessel_M, yb[1]])

    x_init = np.linspace(0.0, L, 60)
    # Initial guess: the first-order profile (good starting point).
    k_eff = Vmax_M_per_s / Km_M if Km_M > 0 else 1.0
    lam0 = penetration_depth_um(D_cm2_per_s, k_eff) if k_eff > 0 else L
    y_init = np.vstack([
        C_vessel_M * np.cosh((L - x_init) / lam0) / np.cosh(L / lam0),
        np.zeros_like(x_init)])

    sol = solve_bvp(rhs, bc, x_init, y_init, tol=1e-8, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(f"BVP solve failed: {sol.message}")

    xs = np.linspace(0.0, L, n_points)
    cs = np.maximum(sol.sol(xs)[0], 0.0)
    return TissueProfile(x_um=[float(v) for v in xs],
                         C_M=[float(v) for v in cs],
                         C_vessel_M=C_vessel_M, thickness_um=L,
                         penetration_depth_um=None, method="numerical",
                         trust=trust)
