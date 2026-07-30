"""Joint Monte-Carlo propagation of prior uncertainty through a cell readout.

The other cell modules each give an EXACT interval CI, but only for their
own single source and only because their readout is monotone in one K_d.
When several priors feed one readout (competition, or a whole
permeation→binding→competition scenario), the honest uncertainty is the
JOINT distribution of the output given the joint distribution of the
inputs. This module samples it.

Fundamental sampled quantity: the binding FREE ENERGY, not K_d. A prior's
K_d CI comes from a Gaussian ΔG (ΔG ~ N(μ, σ), K_d = exp(ΔG/RT)), so we
sample ΔG and map to K_d — giving the correct log-normal K_d spread. The
σ recovered from the prior is its EFFECTIVE σ (max of precision and
measured accuracy; see `src/bridge`), so the propagation carries accuracy,
not just reproducibility — the theme of the whole stack.

Two things make this rigorous rather than decorative:
  * it is DETERMINISTIC given a seed (numpy default_rng), matching the
    project's reproducibility discipline;
  * for a single-prior monotone readout it must AGREE with the exact
    interval-arithmetic CI the analytic modules compute — a mutual
    cross-check (see tests). Where it DISAGREES (multiple correlated
    inputs, non-monotone readouts) is exactly where the analytic interval
    is not valid and the MC is the right tool.

Non-AI: Gaussian sampling of a physics quantity + the deterministic cell
model. No learned surrogate.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional

_KB_KCAL_PER_MOL_K = 0.0019872041

# Trust ordering, worst-first, for combining several priors' verdicts.
_TRUST_RANK = {
    "trustworthy_absolute": 3,
    "rank_order_only": 2,
    "do_not_trust_absolute": 1,
    "uncalibrated": 0,
}


def _prior_dG_mean_sigma(prior, T_K: float = 298.15):
    """Recover (ΔG mean, effective σ) in kcal/mol from a Hill prior.

    ΔG mean from the point K_d; σ backed out of the K_d 95 % CI
    (σ = RT·ln(K_hi/K_lo)/(2·1.96)). No CI → σ = 0 (a point mass — the
    prior carried no uncertainty to propagate).
    """
    if getattr(prior, "type", None) != "hill":
        raise ValueError("Monte-Carlo propagation needs Hill priors")
    RT = _KB_KCAL_PER_MOL_K * T_K
    Kd = float(prior.parameters["Kd_M"])
    dG_mean = RT * math.log(Kd)
    ci = (getattr(prior, "parameter_ci95", None) or {}).get("Kd_M")
    if ci is None:
        return dG_mean, 0.0
    lo, hi = ci
    sigma = RT * math.log(hi / lo) / (2.0 * 1.96)
    return dG_mean, sigma


@dataclass
class MCResult:
    """Distribution of a scalar cell readout under joint prior uncertainty."""

    mean: float
    std: float
    median: float
    ci95: tuple                       # (2.5th, 97.5th percentile)
    n: int
    seed: int
    trust: str                        # worst trust among the input priors
    all_calibrated: bool
    inputs: Dict[str, dict] = field(default_factory=dict)   # per-prior μ/σ
    samples: Optional[List[float]] = None                    # kept if requested

    @property
    def decision_grade(self) -> bool:
        return self.all_calibrated and self.trust == "trustworthy_absolute"

    def summary(self) -> str:
        lo, hi = self.ci95
        grade = "decision-grade" if self.decision_grade else \
            f"NOT decision-grade [{self.trust}]"
        return (f"{self.mean:.4g} ± {self.std:.2g} "
                f"(95% CI {lo:.4g}–{hi:.4g}, n={self.n}, seed={self.seed}) "
                f"— {grade}")


def monte_carlo_propagate(
    priors: Mapping[str, object],
    readout_fn: Callable[[Dict[str, float]], float],
    *,
    n: int = 4000,
    seed: int = 0,
    T_K: float = 298.15,
    keep_samples: bool = False,
) -> MCResult:
    """Propagate the priors' joint uncertainty through `readout_fn`.

    Args:
        priors: {name: RateLawPrior}. Each name is a knob the readout uses.
        readout_fn: maps {name: sampled K_d in M} → a scalar cell-level
            output (occupancy, complex fraction, …).
        n: sample count. seed: RNG seed (deterministic).
        keep_samples: retain the raw output samples on the result.

    Returns MCResult with mean/std/median, a percentile 95 % CI, and the
    worst trust verdict across the inputs (a readout is only as trustworthy
    as its least-trustworthy K_d).
    """
    import numpy as np

    if not priors:
        raise ValueError("need at least one prior")
    if n < 2:
        raise ValueError("n must be ≥ 2")

    rng = np.random.default_rng(seed)
    RT = _KB_KCAL_PER_MOL_K * T_K

    names = list(priors.keys())
    specs = {nm: _prior_dG_mean_sigma(priors[nm], T_K) for nm in names}

    # Pre-draw all Gaussian ΔG deviates: shape (n, n_priors).
    z = rng.standard_normal(size=(n, len(names)))
    outputs = np.empty(n, dtype=float)
    for i in range(n):
        kd_map = {}
        for j, nm in enumerate(names):
            mu, sig = specs[nm]
            dG = mu + sig * z[i, j] if sig > 0 else mu
            kd_map[nm] = math.exp(dG / RT)
        outputs[i] = float(readout_fn(kd_map))

    trusts = [getattr(priors[nm], "trust", "uncalibrated") for nm in names]
    worst_trust = min(trusts, key=lambda t: _TRUST_RANK.get(t, 0))
    all_cal = all(getattr(priors[nm], "accuracy_known", False) for nm in names)

    lo, hi = (float(x) for x in np.percentile(outputs, [2.5, 97.5]))
    return MCResult(
        mean=float(np.mean(outputs)),
        std=float(np.std(outputs, ddof=1)),
        median=float(np.median(outputs)),
        ci95=(lo, hi),
        n=n,
        seed=seed,
        trust=worst_trust,
        all_calibrated=all_cal,
        inputs={nm: {"dG_mean_kcalmol": specs[nm][0],
                     "dG_sigma_kcalmol": specs[nm][1]} for nm in names},
        samples=[float(x) for x in outputs] if keep_samples else None,
    )
