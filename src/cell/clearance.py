"""Metabolic clearance — the drug is consumed, so exposure decays.

Every module so far holds the extracellular concentration fixed. Real
dosing does not: the liver metabolises the drug, so exposure rises after a
dose and then decays. This matters for fate because killing needs the
occupancy held ABOVE the critical threshold θ* (see `fate`) — a very
potent drug that is cleared in an hour may never sustain θ* at all. Potency
without exposure is not efficacy.

Well-stirred hepatic model (Pang & Rowland 1977, J Pharmacokinet Biopharm
5:625 — the standard):

    CL_h = Q_h·fu·CL_int / (Q_h + fu·CL_int)
    E_h  = CL_h / Q_h                       (extraction ratio)
    F_h  = 1 − E_h                          (first-pass survival)

Two limits fall out and are worth stating because they are the whole
reason the model is shaped this way:
  * CL_int → ∞ gives CL_h → Q_h — you cannot clear a drug faster than
    blood delivers it to the liver (FLOW-limited);
  * CL_int small gives CL_h → fu·CL_int (CAPACITY-limited).

Elimination is then first-order: k_el = CL_h/V_d, t½ = ln2/k_el, and the
fundamental PK identity AUC = Dose/CL holds exactly (a test anchors it).

WHERE THE NUMBER COMES FROM — honesty note. `src/quantum/metabolism.py`
predicts the SITE of CYP metabolism (which C–H is abstracted), and per
`docs/campaign1_closeout.md` it is ADVISORY: it is not a validated rate
predictor, and it is blind to N-dealkylation. So CL_int here is a MEASURED
input (microsomal/hepatocyte assay) carrying its own trust flag. The SoM
prediction can accompany it as a qualitative flag ("a labile site exists")
but never substitutes for the rate. We do not turn a bond energy into a
clearance number.

Cited: human hepatic blood flow Q_h ≈ 90 L/h for a 70 kg adult (20.7
mL/min/kg — Davies & Morris 1993, Pharm Res 10:1093).

Non-AI: closed-form PK. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

# Human hepatic blood flow, 70 kg adult (Davies & Morris 1993,
# Pharm Res 10:1093 — 20.7 mL/min/kg).
_Q_HEPATIC_L_PER_H = 90.0
_LN2 = math.log(2.0)


@dataclass
class HepaticClearance:
    """Well-stirred hepatic clearance of a drug.

    CL_int is INTRINSIC clearance (L/h) from a microsomal or hepatocyte
    assay — a measured input, not something we derive from structure.
    """

    CLint_L_per_h: float
    fu_plasma: float = 1.0                 # unbound fraction in plasma
    Vd_L: float = 70.0                     # volume of distribution
    Q_hepatic_L_per_h: float = _Q_HEPATIC_L_PER_H
    source: str = "input:microsomal-assay"
    trust: str = "uncalibrated"
    som_flag: str = ""                     # advisory SoM note, never a rate

    def __post_init__(self):
        if self.CLint_L_per_h < 0 or self.Vd_L <= 0 or \
                self.Q_hepatic_L_per_h <= 0:
            raise ValueError("CLint ≥ 0, Vd > 0, Q_h > 0 required")
        if not (0.0 < self.fu_plasma <= 1.0):
            raise ValueError("fu_plasma must be in (0,1]")

    @property
    def CL_hepatic_L_per_h(self) -> float:
        """Well-stirred hepatic clearance — never exceeds hepatic flow."""
        fu_cl = self.fu_plasma * self.CLint_L_per_h
        return self.Q_hepatic_L_per_h * fu_cl / (self.Q_hepatic_L_per_h + fu_cl)

    @property
    def extraction_ratio(self) -> float:
        return self.CL_hepatic_L_per_h / self.Q_hepatic_L_per_h

    @property
    def bioavailability_first_pass(self) -> float:
        return 1.0 - self.extraction_ratio

    @property
    def k_elimination_per_h(self) -> float:
        return self.CL_hepatic_L_per_h / self.Vd_L

    @property
    def half_life_h(self) -> Optional[float]:
        k = self.k_elimination_per_h
        return (_LN2 / k) if k > 0 else None

    def summary(self) -> str:
        t = self.half_life_h
        ts = f"t½={t:.1f} h" if t is not None else "t½=∞ (no clearance)"
        return (f"CL_h={self.CL_hepatic_L_per_h:.1f} L/h "
                f"(E_h={self.extraction_ratio:.2f}, F={self.bioavailability_first_pass:.2f}), "
                f"{ts} [{self.trust}]")


@dataclass
class ExposureProfile:
    """Plasma concentration vs time under a dosing regimen."""

    t_h: List[float]
    C_M: List[float]
    dose_mol: float
    n_doses: int
    interval_h: Optional[float]
    clearance: HepaticClearance

    @property
    def Cmax_M(self) -> float:
        return max(self.C_M) if self.C_M else 0.0

    def auc_M_h(self) -> float:
        """Trapezoidal AUC over the simulated window (M·h)."""
        total = 0.0
        for (t0, c0), (t1, c1) in zip(zip(self.t_h, self.C_M),
                                      zip(self.t_h[1:], self.C_M[1:])):
            total += 0.5 * (c0 + c1) * (t1 - t0)
        return total

    def time_above_h(self, threshold_M: float) -> float:
        """Hours spent above a concentration threshold (trapezoid crossing)."""
        total = 0.0
        for (t0, c0), (t1, c1) in zip(zip(self.t_h, self.C_M),
                                      zip(self.t_h[1:], self.C_M[1:])):
            dt = t1 - t0
            if c0 >= threshold_M and c1 >= threshold_M:
                total += dt
            elif c0 >= threshold_M or c1 >= threshold_M:
                # Linear crossing inside the interval.
                frac = abs((max(c0, c1) - threshold_M) / (c1 - c0)) \
                    if c1 != c0 else 1.0
                total += dt * min(1.0, max(0.0, frac))
        return total

    def summary(self) -> str:
        return (f"Cmax={self.Cmax_M:.2e} M, AUC={self.auc_M_h():.2e} M·h "
                f"over {self.t_h[-1]:.0f} h, {self.n_doses} dose(s)")


def single_dose_exposure(dose_mol: float, clearance: HepaticClearance,
                         *, duration_h: Optional[float] = None,
                         n_points: int = 400) -> ExposureProfile:
    """IV-bolus single dose: C(t) = (Dose/V_d)·exp(−k_el·t)."""
    if dose_mol < 0:
        raise ValueError("dose_mol must be ≥ 0")
    k = clearance.k_elimination_per_h
    C0 = dose_mol / clearance.Vd_L
    if duration_h is None:
        duration_h = (7.0 * clearance.half_life_h) if clearance.half_life_h \
            else 24.0
    ts = [duration_h * i / (n_points - 1) for i in range(n_points)]
    cs = [C0 * math.exp(-k * t) for t in ts]
    return ExposureProfile(t_h=ts, C_M=cs, dose_mol=dose_mol, n_doses=1,
                           interval_h=None, clearance=clearance)


def repeat_dose_exposure(dose_mol: float, clearance: HepaticClearance,
                         *, interval_h: float, n_doses: int,
                         n_points_per_interval: int = 60) -> ExposureProfile:
    """Repeated IV boluses — superposition of decaying single doses."""
    if interval_h <= 0 or n_doses < 1:
        raise ValueError("interval_h > 0 and n_doses ≥ 1 required")
    k = clearance.k_elimination_per_h
    C0 = dose_mol / clearance.Vd_L
    duration = interval_h * n_doses
    n_pts = n_points_per_interval * n_doses
    ts = [duration * i / (n_pts - 1) for i in range(n_pts)]
    cs = []
    for t in ts:
        c = 0.0
        for d in range(n_doses):
            t_d = d * interval_h
            if t >= t_d:
                c += C0 * math.exp(-k * (t - t_d))
        cs.append(c)
    return ExposureProfile(t_h=ts, C_M=cs, dose_mol=dose_mol,
                           n_doses=n_doses, interval_h=interval_h,
                           clearance=clearance)


def accumulation_ratio(clearance: HepaticClearance, interval_h: float
                       ) -> float:
    """Steady-state accumulation on repeat dosing: 1/(1 − e^{−k·τ})."""
    k = clearance.k_elimination_per_h
    if k <= 0:
        return math.inf
    return 1.0 / (1.0 - math.exp(-k * interval_h))


@dataclass
class RegimenOutcome:
    """Does this dosing regimen actually kill? (PK/PD coupling)"""

    hours_above_threshold: float
    fraction_of_time_above: float
    log10_kill: float                  # net log-kill over the regimen
    net_effect: str                    # 'kills' | 'holds' | 'fails'
    critical_occupancy: Optional[float]
    threshold_conc_M: Optional[float]

    def summary(self) -> str:
        thr = (f"C*={self.threshold_conc_M:.2e} M" if self.threshold_conc_M
               else "no threshold (drug cannot kill)")
        return (f"{thr}; above it {self.hours_above_threshold:.1f} h "
                f"({100*self.fraction_of_time_above:.0f}% of regimen) → "
                f"{self.log10_kill:+.2f} log10 kill — {self.net_effect.upper()}")


def evaluate_regimen(exposure: ExposureProfile, prior,
                     params=None) -> RegimenOutcome:
    """Integrate cell fate over a time-varying exposure.

    The payoff of coupling PK to PD: computes the concentration C* that
    achieves the critical occupancy θ*, how long the regimen holds above
    it, and the net log-kill ∫k_net(θ(C(t)))·dt — the honest answer to
    "does this dosing schedule work?", which neither potency nor exposure
    alone can give.
    """
    from src.cell.fate import (CellFateParams, critical_occupancy,
                               fate_from_occupancy)
    from src.cell.occupancy import hill_occupancy

    p = params or CellFateParams()
    Kd = float(prior.parameters["Kd_M"])
    n_hill = float(prior.parameters.get("n_hill", 1.0))
    theta_star = critical_occupancy(p)

    # Concentration achieving θ*: invert the Hill equation.
    C_star = None
    if theta_star is not None and 0.0 < theta_star < 1.0:
        C_star = Kd * (theta_star / (1.0 - theta_star)) ** (1.0 / n_hill)

    hours_above = exposure.time_above_h(C_star) if C_star else 0.0
    total_h = exposure.t_h[-1] - exposure.t_h[0]

    # Net log-kill: integrate k_net over the trajectory.
    integral = 0.0
    for (t0, c0), (t1, c1) in zip(zip(exposure.t_h, exposure.C_M),
                                  zip(exposure.t_h[1:], exposure.C_M[1:])):
        k0 = fate_from_occupancy(hill_occupancy(c0, Kd, n_hill), p).k_net_per_s
        k1 = fate_from_occupancy(hill_occupancy(c1, Kd, n_hill), p).k_net_per_s
        integral += 0.5 * (k0 + k1) * (t1 - t0) * 3600.0     # h → s
    log10_kill = -integral / math.log(10.0)

    if log10_kill > 0.5:
        verdict = "kills"
    elif log10_kill > -0.5:
        verdict = "holds"
    else:
        verdict = "fails"

    return RegimenOutcome(
        hours_above_threshold=hours_above,
        fraction_of_time_above=(hours_above / total_h) if total_h > 0 else 0.0,
        log10_kill=log10_kill,
        net_effect=verdict,
        critical_occupancy=theta_star,
        threshold_conc_M=C_star,
    )
