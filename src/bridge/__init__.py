"""src/bridge — Campaign-1 → Campaign-2 rate-law emitter.

This module is the EXIT PATH from the professor's "closed
ontology / circular validation" critique. Campaign-2 pathway
models will not hand-tune binding rate constants; they will
cite this module's output with explicit method provenance and
calibrated uncertainty.

What ships here today (Campaign-1 scope)
----------------------------------------
- `binding_to_hill(dG_kcalmol, *, T_K=298.15, n_hill=1.0,
   uncertainty_kcalmol=None)` — convert an absolute binding ΔG
   (typically from `cellsim fep-binding bench --sample`) into a
   Hill-equation rate-law record with CI.
- `affinity_to_michaelis(kcat_per_s, KM_M, *,
   uncertainty_kcat=None, uncertainty_KM=None)` — same idea for
   enzyme kinetics. Pass-through for now since FEP doesn't
   produce kcat directly; the module exists so Campaign-2
   loaders have a single place to import from.

Provenance dataclass (`RateLawPrior`) carries:
  - the rate-law type ('hill' | 'michaelis')
  - the parameter values + CI (95%)
  - source ('FEP', 'phenomenological', 'literature')
  - method tag ('amber14-DDM-MBAR' or whatever produced the input)
  - ACCURACY: a measured target-class error + a `trust` verdict, so the
    CI reflects how RIGHT the calculation is (accuracy), not merely how
    REPRODUCIBLE it is (precision). The two diverge badly — a ±0.29
    kcal/mol biotin seed bar hid an ~11.8 kcal/mol true error — and a
    cell model must not inherit a tight bar on a wrong rate.
  - timestamp

Why a thin module: Campaign-2's signalling ODEs read these
priors at config time. Centralising the conversion + provenance
here means any Campaign-2 model can be audited back to a
specific FEP run; if we later switch from amber14 to ff19SB,
the audit trail tells us which rate-law records need re-emission.

Scope NOT in this commit (Campaign-2 work, gated by prof's
Milestone-A clearance):
  - Anything that runs a cell-level ODE.
  - Bulk emission to a Campaign-2 YAML config.
  - Cache integration (the `(ligand_hash, receptor_hash) → ΔG`
    lookup that lets Campaign 2 ingest 10⁴ compounds).

Non-AI: the conversions are closed-form thermodynamics
(Kd = exp(ΔG/RT)) + propagation of σ via the dG-uncertainty.
No learned surrogate.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Optional


# Boltzmann constant in kcal/mol/K (the unit FEP outputs use).
_KB_KCAL_PER_MOL_K = 0.0019872041


@dataclass
class RateLawPrior:
    """Provenance-tracked rate-law record consumed by Campaign-2.

    All numeric fields use SI molar units where applicable so a
    Campaign-2 ODE can plug them in without unit fiddling.
    """
    type: str                       # 'hill' or 'michaelis'
    parameters: dict                # see per-type docs below
    parameter_ci95: dict             # 2-tuples (lo, hi) per param
    source: str                     # 'FEP' | 'phenomenological' | 'literature'
    method: str                     # e.g. 'amber14-DDM-MBAR'
    temperature_K: float
    # ACCURACY, not just precision. The CI above is driven by the LARGER
    # of the input reproducibility σ and the measured target-class error
    # (see `_resolve_accuracy`), so a Campaign-2 cell model never inherits
    # a bar tighter than the calculation is actually known to be right —
    # the precision≠accuracy trap fixed one layer down in docking.
    accuracy_kcalmol: Optional[float] = None   # measured accuracy σ, or None
    accuracy_basis: str = "unknown"            # where accuracy came from
    trust: str = "uncalibrated"                # verdict for the consumer
    timestamp_iso: str = field(
        default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                time.gmtime()))
    notes: str = ""

    @property
    def accuracy_known(self) -> bool:
        """True only when a MEASURED accuracy backs the CI.

        When False the CI reflects reproducibility (or nothing) and the
        true error is unknown — a Campaign-2 consumer must treat the rate
        as rank-order-only, not as a calibrated absolute value.
        """
        return self.accuracy_kcalmol is not None

    def summary(self) -> str:
        """One-line summary suitable for Campaign-2 config logs."""
        params = ", ".join(
            f"{k}={v}" for k, v in self.parameters.items())
        acc = (f"±{self.accuracy_kcalmol:.1f}kcal/mol"
               if self.accuracy_known else "accuracy=UNKNOWN")
        return (f"[{self.type}] {params}  "
                f"src={self.source}/{self.method}  "
                f"T={self.temperature_K:.1f}K  "
                f"trust={self.trust} ({acc})")


def _resolve_accuracy(method: str, receptor):
    """Return (accuracy_kcalmol, basis, trust) for a binding ΔG.

    Accuracy is the MEASURED absolute error of the source calculation —
    distinct from the input reproducibility σ. Rules (NEVER GUESS, same
    as `src/uq/reliability`):

    * FEP-sourced ΔG uses a different, more accurate method than docking,
      so the docking reliability table does NOT apply to it. CellSim's
      absolute FEP binding ΔG is not yet GPU-validated, so an FEP rate
      law is currently `uncalibrated` for accuracy — flagged, not faked.
    * Docking-sourced ΔG uses the measured target-class error from
      `reliability_table.yaml` when the receptor is known and in the
      table; an unknown receptor stays `uncalibrated`.

    We never let the input precision masquerade as accuracy: if accuracy
    cannot be established we return None, and the caller's CI will carry
    a `trust` flag saying so.
    """
    m = (method or "").lower()
    if m.startswith("fep") or "ddm" in m or "mbar" in m:
        return None, "fep-unvalidated", "uncalibrated"
    if receptor is not None:
        try:
            from src.uq.reliability import reliability_for
            rel = reliability_for(receptor)
        except Exception:  # noqa: BLE001 — accuracy is optional, never fatal
            rel = None
        if rel is not None and rel.calibrated:
            return (rel.mean_abs_err_kcalmol,
                    f"docking-target-class:{rel.target_class}",
                    rel.verdict)
        return None, "docking-uncalibrated-target", "uncalibrated"
    return None, "unknown", "uncalibrated"


def binding_to_hill(
    dG_kcalmol: float,
    *,
    T_K: float = 298.15,
    n_hill: float = 1.0,
    uncertainty_kcalmol: Optional[float] = None,
    receptor=None,
    method: str = "FEP-DDM-MBAR",
) -> RateLawPrior:
    """Absolute binding ΔG → Hill-equation rate-law prior.

    Hill equation:
        fraction_bound([L]) = [L]^n / (K_d^n + [L]^n)

    where K_d = exp(ΔG / RT). Uncertainty in ΔG (1σ) propagates
    multiplicatively to K_d:
        σ(ln K_d) = σ(ΔG) / RT
        K_d_lo = K_d × exp(-1.96 × σ_ln)   (95% CI lower)
        K_d_hi = K_d × exp(+1.96 × σ_ln)

    Args:
        dG_kcalmol: absolute binding free energy. Negative for
            binders; sign convention matches `compute_absolute_
            binding_dg` output.
        T_K: temperature (default 298.15 K — physiological-ish,
            matches the FEP sampler default).
        n_hill: Hill cooperativity coefficient (default 1.0 —
            non-cooperative single-site binding; Campaign-2
            allosteric models override).
        uncertainty_kcalmol: 1σ uncertainty on ΔG (typically MBAR
            σ from the bench CSV's uncertainty_kcalmol column).
            None → no CI propagated; the resulting record will
            have parameter_ci95 = {} so the consumer knows the
            estimate is uncalibrated.
        method: provenance string (default 'FEP-DDM-MBAR' for
            absolute binding from compute_absolute_binding_dg;
            override to 'phenomenological-Lipinski' if the input
            came from the OLD/ Tier-0 BindingMatcher).

    Returns: RateLawPrior with parameters['Kd_M'] = K_d in M
             (NOT mM), parameters['n_hill'] = n_hill.
    """
    RT = _KB_KCAL_PER_MOL_K * T_K
    Kd_M = math.exp(dG_kcalmol / RT)

    parameters = {"Kd_M": Kd_M, "n_hill": float(n_hill)}

    # Precision (input reproducibility σ) vs accuracy (measured
    # target-class error). The CI must reflect the LARGER of the two so a
    # Campaign-2 model never gets a bar tighter than the calculation is
    # known to be correct — the biotin case where a ±0.29 seed bar hid an
    # ~11.8 kcal/mol true error must not propagate up as a tight K_d.
    precision_sigma = (uncertainty_kcalmol
                       if uncertainty_kcalmol and uncertainty_kcalmol > 0
                       else None)
    accuracy_sigma, accuracy_basis, trust = _resolve_accuracy(method,
                                                              receptor)
    candidate_sigmas = [s for s in (precision_sigma, accuracy_sigma) if s]
    sigma_eff = max(candidate_sigmas) if candidate_sigmas else None

    parameter_ci95: dict = {}
    if sigma_eff is not None:
        sigma_ln_Kd = sigma_eff / RT
        Kd_lo = Kd_M * math.exp(-1.96 * sigma_ln_Kd)
        Kd_hi = Kd_M * math.exp(+1.96 * sigma_ln_Kd)
        parameter_ci95 = {"Kd_M": (Kd_lo, Kd_hi)}

    notes = ""
    if trust == "uncalibrated":
        notes = ("accuracy UNKNOWN for this source — CI reflects "
                 "reproducibility only; treat as rank-order, not a "
                 "calibrated absolute K_d")
    elif accuracy_sigma and precision_sigma and accuracy_sigma > precision_sigma:
        notes = (f"CI widened from input σ={precision_sigma:.2f} to measured "
                 f"accuracy σ={accuracy_sigma:.2f} kcal/mol ({accuracy_basis})")

    return RateLawPrior(
        type="hill",
        parameters=parameters,
        parameter_ci95=parameter_ci95,
        source="FEP" if method.startswith("FEP") else "phenomenological",
        method=method,
        temperature_K=T_K,
        accuracy_kcalmol=accuracy_sigma,
        accuracy_basis=accuracy_basis,
        trust=trust,
        notes=notes,
    )


def affinity_to_michaelis(
    kcat_per_s: float,
    KM_M: float,
    *,
    T_K: float = 298.15,
    uncertainty_kcat: Optional[float] = None,
    uncertainty_KM: Optional[float] = None,
    method: str = "literature",
) -> RateLawPrior:
    """Enzyme kinetics → Michaelis-Menten rate-law prior.

    Pass-through for now: FEP doesn't directly produce kcat (that
    requires transition-state theory + a reactive sub-tier we
    haven't implemented). The function exists so Campaign-2
    loaders have a single import surface — when CellSim grows a
    reactive-FEP / Eyring TS module, the body changes here and
    every downstream consumer benefits.

    Args:
        kcat_per_s: turnover number in s⁻¹.
        KM_M: Michaelis constant in M.
        T_K: temperature.
        uncertainty_kcat: 1σ on kcat in s⁻¹ (default None → no CI).
        uncertainty_KM: 1σ on KM in M (default None → no CI).
        method: provenance ('literature' if the values are from
            BRENDA / SABIO; 'reactive-FEP' once that ships).

    Returns: RateLawPrior with type='michaelis'.
    """
    parameters = {"kcat_per_s": kcat_per_s, "KM_M": KM_M}
    parameter_ci95: dict = {}
    if uncertainty_kcat is not None and uncertainty_kcat > 0:
        parameter_ci95["kcat_per_s"] = (
            kcat_per_s - 1.96 * uncertainty_kcat,
            kcat_per_s + 1.96 * uncertainty_kcat)
    if uncertainty_KM is not None and uncertainty_KM > 0:
        parameter_ci95["KM_M"] = (
            KM_M - 1.96 * uncertainty_KM,
            KM_M + 1.96 * uncertainty_KM)

    return RateLawPrior(
        type="michaelis",
        parameters=parameters,
        parameter_ci95=parameter_ci95,
        source="FEP" if method.startswith("reactive-FEP") else
               "literature" if method == "literature" else
               "phenomenological",
        method=method,
        temperature_K=T_K,
    )


__all__ = [
    "RateLawPrior",
    "binding_to_hill",
    "affinity_to_michaelis",
]
