"""Split-conformal prediction for distribution-free ΔG bounds.

Post-hoc statistical wrapper (Vovk et al 2005, Shafer & Vovk 2008,
Angelopoulos & Bates 2023 tutorial). Given a calibration set of
(predicted, true) pairs, attach a calibrated CI to any new
prediction:

    y_hat ± q

where `q` is the (1-α) quantile of the absolute calibration
residuals. With enough calibration points this satisfies

    P(|y_pred - y_true| ≤ q) ≥ 1 - α

*without any parametric assumption* on the error distribution.
MAPIE is an industrial-strength implementation of this; we use a
simple explicit version because our Layer-1.6 wants deterministic
seedless math that doesn't pull in sklearn as a hard dependency.

Per MISSION.md §"No black-box / no AI surrogates":
Conformal prediction is acceptable as a *post-hoc non-parametric
statistical wrapper* around physics outputs. It provides
statistical bounds, **not** mechanistic insight — this is
documented on the result envelope so biologists do not mistake
the CI for a physics claim. The underlying calibration set's
quality caps everything.

Usage:
    from src.uq import ConformalBounds

    # Bundle of (Vina ΔG, experimental ΔG) from literature pairs.
    cb = ConformalBounds(alpha=0.05)
    cb.calibrate(predictions=[-7.4, -9.2, -6.1, -8.0, -5.5, -7.9],
                 truths=[-7.1, -9.8, -5.7, -8.5, -5.2, -8.1])
    # {'q': 0.60, 'alpha': 0.05, 'n_cal': 6, 'coverage_target': 0.95}

    lo, hi = cb.interval(-7.44)
    # (-8.04, -6.84)  — 95 % CI under the calibration distribution

Caveat: at n_cal < 20 the conformal bound is conservative and the
coverage guarantee weakens (Angelopoulos 2023). Biologists should
only trust the CI if the calibration set is representative of the
compound class they're predicting on.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional
import logging

logger = logging.getLogger(__name__)


@dataclass
class ConformalBounds:
    """Split-conformal wrapper; calibrate once, then
    `interval(prediction) → (lo, hi)` on every new prediction."""

    alpha: float = 0.05       # miscoverage rate (0.05 → 95 % CI)

    q: Optional[float] = None              # conformity quantile
    n_cal: int = 0
    cal_residuals: list = field(default_factory=list)
    cal_metadata: dict = field(default_factory=dict)
    # Set by calibrate(). True => n_cal was below (1-alpha)/alpha, so the
    # requested level was unreachable, q collapsed to the max residual,
    # and the "interval" has NO coverage guarantee.
    level_saturated: bool = False
    n_min_for_alpha: int = 0

    def calibrate(
        self, predictions: list, truths: list, *,
        metadata: Optional[dict] = None,
    ) -> dict:
        """Compute the conformity quantile q from calibration pairs.

        Uses the (1 - alpha)(1 + 1/n) finite-sample corrected
        quantile. Returns a summary dict for logging.
        """
        import numpy as np

        if len(predictions) != len(truths):
            raise ValueError(
                f"len(predictions)={len(predictions)} != "
                f"len(truths)={len(truths)}")
        if len(predictions) < 2:
            raise ValueError("need >= 2 calibration pairs")

        preds = np.asarray(predictions, dtype=float)
        truv = np.asarray(truths, dtype=float)
        residuals = np.abs(preds - truv).tolist()

        n = len(residuals)
        level = min(1.0, (1.0 - self.alpha) * (1.0 + 1.0 / n))
        # k-th smallest, 1-indexed. k = ceil(n * level).
        import math
        k = max(1, min(n, math.ceil(n * level)))
        self.q = float(sorted(residuals)[k - 1])

        # A split-conformal (1-alpha) quantile only EXISTS when
        # (1-alpha)(1 + 1/n) <= 1, i.e. n >= (1-alpha)/alpha
        # (n >= 19 for alpha=0.05). Below that the level saturates at
        # 1.0 and q silently degrades to "the largest residual we
        # happened to see" — which carries NO coverage guarantee and
        # is driven entirely by the worst outlier. Flag it loudly
        # instead of returning it dressed up as a 95 % interval.
        n_min = math.ceil((1.0 - self.alpha) / self.alpha)
        self.level_saturated = bool(level >= 1.0 or n < n_min)
        self.n_min_for_alpha = n_min
        if self.level_saturated:
            logger.warning(
                "conformal: n_cal=%d < %d needed for a genuine %.0f%% "
                "interval at alpha=%.2f — q=%.2f is just the MAXIMUM "
                "residual and has no coverage guarantee. Widen alpha or "
                "add calibration points before quoting this as a CI.",
                n, n_min, 100 * (1 - self.alpha), self.alpha, self.q)
        self.n_cal = n
        self.cal_residuals = residuals
        self.cal_metadata = metadata or {}
        return {
            "q_kcalmol": self.q, "alpha": self.alpha,
            "n_cal": n,
            "coverage_target": 1 - self.alpha,
            "cal_residuals_mean": float(np.mean(residuals)),
            "cal_residuals_max": float(np.max(residuals)),
            "k_quantile_index": k,
            # True => n_cal was too small for the requested alpha, so
            # q is the max residual and carries NO coverage guarantee.
            "level_saturated": self.level_saturated,
            "n_min_for_alpha": self.n_min_for_alpha,
        }

    def interval(self, prediction: float) -> tuple[float, float]:
        """Return (lo, hi) conformal interval around a prediction."""
        if self.q is None:
            raise RuntimeError(
                "ConformalBounds not calibrated — "
                "call .calibrate(predictions, truths) first")
        return (float(prediction) - self.q,
                float(prediction) + self.q)

    def check_coverage(self, predictions: list, truths: list) -> float:
        """Empirical coverage on a held-out test set.

        Returns the fraction of truths falling inside the conformal
        interval around each prediction. Should approach
        (1 - alpha) as the test set grows, by construction.
        """
        if self.q is None:
            raise RuntimeError("calibrate first")
        covered = 0
        for p, t in zip(predictions, truths):
            lo, hi = self.interval(p)
            if lo <= t <= hi:
                covered += 1
        return covered / max(1, len(predictions))

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if self.q is None:
            return "[NOT CALIBRATED]"
        return (f"ConformalBounds α={self.alpha}  n_cal={self.n_cal}  "
                f"q = {self.q:.3f} kcal/mol  "
                f"(95 % CI = y_pred ± {self.q:.3f})  "
                f"  NOTE: statistical bound only, not mechanistic")


if __name__ == "__main__":
    # Tiny demo on a synthetic calibration set.
    import numpy as np
    rng = np.random.default_rng(1)
    true_dG = rng.uniform(-10, -4, size=20)
    # Simulate Vina's typical systematic offset + noise.
    vina_dG = true_dG + rng.normal(0.0, 0.8, size=20)
    cb = ConformalBounds(alpha=0.05)
    info = cb.calibrate(predictions=vina_dG.tolist(),
                        truths=true_dG.tolist())
    print("calibration:", info)
    print(cb.summary())
    for pred in [-7.4, -9.0, -5.0]:
        lo, hi = cb.interval(pred)
        print(f"  pred = {pred:+.2f}  →  95 % CI [{lo:+.2f}, {hi:+.2f}]")
    # Coverage on a fresh test set.
    true_test = rng.uniform(-10, -4, size=100)
    vina_test = true_test + rng.normal(0.0, 0.8, size=100)
    cov = cb.check_coverage(vina_test.tolist(), true_test.tolist())
    print(f"empirical coverage on n=100 test: {cov:.3f} "
          f"(target {1 - cb.alpha:.3f})")
