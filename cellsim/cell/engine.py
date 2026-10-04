"""cellsim.cell.engine — Phase-1 single-cell drug-response engine.

One state array, one integrator, one `step`. A population is a set of
representative cells with weights (a weight doubles at division and is
zeroed at death), integrated together with NumPy, so a full dose-response
is a single simulation in which each cell carries its own extracellular
concentration.

Biology ported from the frozen C++ prototype (OLD/src/simulation) and
re-expressed in hours with the published constants it cited:

* Cell cycle — the Novak-Tyson-style 7-variable CDK/cyclin oscillator
  (CycD, Rb, E2F, CycE, CycA, CycB, p21) exactly as in
  OLD/src/simulation/Simulation.h `CDKState::step`, with a single rate
  scale calibrated so the untreated population doubles at the cell
  line's quoted doubling time.
* DNA-damage response — ATM sigmoid on the damage index (threshold 0.10,
  Rothkamm 2003), p53 with MDM2-mediated Michaelis degradation shielded
  by ATM (K_d 0.10, shield 0.70; Picksley 1994, Bottger 1997, Loewer
  2013), MDM2 mRNA as the explicit delay intermediate with Hill n=8
  transactivation (Geva-Zatorsky 2006; Riley 2008), MDM2 protein t1/2
  ~30 min (Momand 1992). The whole block is slowed by a factor
  3.03/5.5 so the sustained-damage pulse period matches Purvis & Lahav
  2012 (5.5 h) instead of the prototype's 3.0 h.
* p21 — p53-driven CDKN1A transcription, Hill n=8, K 0.18 (El-Deiry
  1993) plus a p53-independent CHK1/CHK2 S/G2 term.
* Apoptosis — the prototype's intrinsic-pathway topology (p53 -> PUMA
  -> BAX activation against a Bcl-2 reserve -> MOMP switch -> cyt c and
  Smac release -> caspase-9 -> caspase-3, XIAP neutralised by Smac),
  re-parametrised per hour so that MOMP is a sharp irreversible switch
  and caspase-3 crosses 0.5 about 20-30 min after MOMP (Albeck 2008).
  The prototype's XIAP/Smac stoichiometry (XIAP 1.0 vs Smac 0.5) kept
  caspase-3 clamped at 0.04 forever; here Smac is in excess over XIAP
  (Rehm 2006), which is what lets execution actually happen.
* Drug — intracellular concentration by first-order permeation; damage
  index driven by Pt-adduct / TopII mechanisms, or SAC-dependent mitotic
  arrest for tubulin binders with slippage after 24 h (Gascoigne &
  Taylor 2008).

Mutant p53 (``CellLine.p53_functional=False``) removes transactivation
(MDM2 feedback, p21, PUMA) while leaving p53 protein in place, which is
what loss-of-function hotspot mutants do.

Units: hours, µM, dimensionless normalised protein levels as in the
prototype (p53 basal 0.089 ≈ 1x). Non-AI: closed-form ODEs, RK4.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Union

import numpy as np

from cellsim.cell.library import CellLine, Drug

# ── state layout ──────────────────────────────────────────────────────
NAMES = [
    "CycD", "Rb", "E2F", "CycE", "CycA", "CycB", "p21",   # 0-6   cycle
    "ATM", "p53", "MDM2_mRNA", "MDM2",                    # 7-10  stress
    "D",                                                  # 11    DNA damage index
    "Puma", "Bax", "MOMP", "CytC", "Smac", "C9", "C3",    # 12-18 death
    "Exec",                                               # 19    execution progress
    "Cin",                                                # 20    intracellular drug (µM)
]
IX = {n: i for i, n in enumerate(NAMES)}
NV = len(NAMES)
# One intracellular-concentration slot per drug, appended after the base
# variables, so a combination carries each drug's own Cin. A single drug
# or none leaves the layout exactly as it was.
N_BASE = IX["Cin"]
_UPPER_BASE = [1.5, 1, 1, 1.5, 1.5, 1.5, 1,   # cycle
               1, 2, 2, 2,                     # stress
               2,                              # damage
               5, 1, 1, 1, 1, 1, 1,            # death
               1]                              # exec
_UPPER = np.array(_UPPER_BASE + [np.inf])      # + one Cin


def state_width(n_drugs: int) -> int:
    """Number of state variables for a run carrying `n_drugs` drugs."""
    return N_BASE + max(1, n_drugs)


def _upper_for(n_drugs: int) -> np.ndarray:
    return np.array(_UPPER_BASE + [np.inf] * max(1, n_drugs))


def _as_drug_tuple(drug):
    """Accept None, one drug, or a sequence; always return a tuple."""
    if drug is None:
        return ()
    if isinstance(drug, (list, tuple)):
        return tuple(drug)
    return (drug,)


@dataclass(frozen=True)
class Params:
    """Cited constants (see module docstring). Hours unless stated."""
    # DNA-damage response
    atm_threshold: float = 0.10      # Rothkamm 2003
    atm_gain: float = 20.0
    atm_tau_h: float = 0.1           # ATM autophosphorylation ~5 min
    p53_Kd: float = 0.10             # Picksley 1994 / Bottger 1997 (normalised)
    p53_k_basal: float = 1.0
    p53_k_deg: float = 10.0
    atm_shield: float = 0.70         # Loewer 2013 / Stewart-Ornstein 2017
    mdm2_K_hill: float = 0.18
    mdm2_k_basal: float = 0.315
    mdm2_k_max: float = 3.0
    mdm2_k_mrna_deg: float = 1.5
    mdm2_k_translate: float = 1.5
    mdm2_k_prot_deg: float = 1.5     # t1/2 ~30 min, Momand 1992
    p53_time_scale: float = 3.03 / 5.5  # prototype period -> Purvis 2012 period
    # p21
    p21_K_hill: float = 0.18         # El-Deiry 1993 response element; n=8
    p21_k_max: float = 0.15
    chk_damage_threshold: float = 0.2
    chk_gain: float = 0.008
    # apoptosis (per hour)
    p53_act_lo: float = 0.12         # linMap of p53 -> apoptotic p53 activity
    p53_act_hi: float = 0.35
    puma_k_on: float = 1.0
    puma_k_off: float = 0.3          # PUMA t1/2 ~2.3 h
    bax_k: float = 0.5
    bax_puma_gain: float = 0.3
    bcl2_block: float = 3.0
    bcl2_neutralised_by_puma: float = 0.2
    # Anti-apoptotic proteins SEQUESTER activated Bax stoichiometrically,
    # so MOMP needs activated Bax to exceed the reserve rather than merely
    # to accumulate more slowly against it (mitochondrial priming; Certo
    # et al. 2006 Cancer Cell 9:351, Letai). Modelling it as a rate
    # penalty instead made the reserve almost irrelevant: an 8-fold change
    # in CellLine.bcl2_level moved the IC50 by 1.1x, because a slower
    # climb still crosses a fixed threshold within 72 h. As a buffer it
    # sets whether the threshold is reachable at all. 0 restores the old
    # rate-only behaviour.
    bcl2_buffer: float = 0.30
    # Saturable P-glycoprotein efflux, applied only to substrate drugs and
    # scaled by the line's ABCB1 activity. Saturable rather than linear
    # because a pump has finite turnover: it defeats a low dose and is
    # overwhelmed by a high one, which is what makes efflux shift an IC50
    # instead of simply rescaling the concentration axis.
    #   d[Cin]/dt gains  - V_max * efflux_level * Cin / (K_m + Cin)
    # 1.0 uM/h reproduces the IC50 spread that the measured ABCB1 range
    # implies (1.5-2.1x); 0 disables efflux entirely.
    efflux_vmax_uM_per_h: float = 1.0
    efflux_km_uM: float = 0.5
    momp_k: float = 10.0
    momp_K: float = 0.5
    momp_n: float = 6.0
    release_k: float = 20.0          # cyt c / Smac release once MOMP (minutes)
    xiap_total: float = 0.3
    smac_total: float = 1.0          # Smac in excess over XIAP (Rehm 2006)
    c9_k: float = 5.0
    c9_xiap_block: float = 10.0
    c3_k: float = 6.0                # caspase-3 crosses 0.5 ~20-30 min after MOMP (Albeck 2008)
    exec_k: float = 2.0
    c3_commit: float = 0.5
    # mitosis
    m_duration_h: float = 1.0
    g2m_damage_block: float = 0.5
    slippage_h: float = 24.0         # Gascoigne & Taylor 2008
    slippage_damage: float = 0.3
    mitotic_death_lag_h: float = 2.0
    # Cell-to-cell variability. Each cell draws log-normal multipliers
    # (mean 1) on its anti-apoptotic Bcl-2 reserve and on its drug
    # accumulation. Pre-existing protein-level differences are the
    # documented origin of fractional killing (Spencer et al. 2009
    # Nature 459:428; Roux et al. 2015 Mol Syst Biol 11:803), and they are
    # what turns an all-or-none population into a graded dose-response.
    # 0 reproduces identical cells exactly.
    het_sigma: float = 0.30
    # p53-independent damage -> PUMA route: ATM / c-Abl activate p73, which
    # transactivates PUMA without p53 (Gong et al. 1999 Nature 399:806;
    # Agami et al. 1999 Nature 399:809, both with cisplatin). Engages only
    # once ATM is genuinely activated (> 0.5; basal ATM sits at 0.12-0.27).
    # 0 = off (p53-only). Default 0.1: chosen by a leave-one-mutant-out
    # test (scripts/experiment_p73.py) in which the gain is fitted on one
    # TP53-mutant line and scored on the other. Both folds independently
    # picked 0.1, and in both it improved the held-out line over the
    # p53-only model and over a constant-IC50 null.
    p73_gain: float = 0.1
    p73_atm_on: float = 0.5
    # Cell-to-cell variability of cycle duration: each cell's interphase
    # time is scaled by a log-normal factor of mean 1 and this coefficient
    # of variation. Real cycles vary: HeLa in the Cell Tracking Challenge
    # time-lapse (Fluo-N2DL-HeLa) completes cycles with CV 0.22-0.33, a
    # lower bound because a 46 h movie truncates long cycles
    # (scripts/ctc_reference.py). 0 keeps every cell on the same clock and
    # every earlier result unchanged.
    cycle_cv: float = 0.0


# ── small helpers ─────────────────────────────────────────────────────
def _hill8(x: np.ndarray, K: float) -> np.ndarray:
    x8 = x ** 8
    return x8 / (K ** 8 + x8)


def _phase(Y: np.ndarray) -> np.ndarray:
    """0=G1 1=S 2=G2 3=M from cyclin levels (prototype thresholds)."""
    CycE, CycA, CycB = Y[:, IX["CycE"]], Y[:, IX["CycA"]], Y[:, IX["CycB"]]
    ph = np.zeros(len(Y), dtype=np.int8)
    ph[(CycA > 0.08) | (CycE > 0.18)] = 1
    ph[CycA > 0.30] = 2
    ph[CycB > 0.25] = 3
    return ph


def fresh_cycle_state() -> np.ndarray:
    """Early-G1 cycle variables as the prototype's resetForNewCycle(gs=1)."""
    y = np.zeros(7)
    y[:] = [0.08, 0.90, 0.03, 0.01, 0.01, 0.01, 0.02]
    return y


def random_cycle_state(rng: np.random.Generator) -> np.ndarray:
    """Asynchronous-culture initial cycle state (Casciari 1992 phase
    fractions 46/33/17/4 %), as the prototype's CDKState::randomize."""
    u = rng.random()
    r = rng.random
    if u < 0.46:
        g = r()
        return np.array([0.10 + g * 0.50, 0.95 - g * 0.40, 0.02 + g * 0.15,
                         0.02 + g * 0.18, 0.01 + g * 0.05, 0.01 + r() * 0.02,
                         0.05 + r() * 0.05])
    if u < 0.79:
        s = r()
        return np.array([0.25 + r() * 0.20, 0.35 - s * 0.15, 0.40 + r() * 0.30,
                         0.20 + r() * 0.25, 0.15 + s * 0.25, 0.02 + r() * 0.05,
                         0.05 + r() * 0.05])
    if u < 0.96:
        g = r()
        return np.array([0.20 + r() * 0.20, 0.20, 0.55, 0.15 + r() * 0.10,
                         0.40 + r() * 0.15, 0.08 + g * 0.15, 0.05 + r() * 0.05])
    return np.array([0.30, 0.15, 0.50, 0.10, 0.55, 0.28 + r() * 0.10,
                     0.05 + r() * 0.05])


# ── the right-hand side ───────────────────────────────────────────────
def rhs(Y: np.ndarray, aux: dict, line: CellLine, drug: Optional[Drug],
        p: Params, k_cyc: float) -> np.ndarray:
    """dY/dt for every cell (vectorised). `aux` carries the per-cell
    non-ODE inputs: C_out (µM), arrest_h, forced_D (or NaN), the
    heterogeneity multipliers, and optionally a growth signal `gs`."""
    dY = np.zeros_like(Y)
    f53 = 1.0 if line.p53_functional else 0.0

    CycD, Rb, E2F, CycE, CycA, CycB, p21 = (Y[:, i] for i in range(7))
    ATM, p53, mRNA, MDM2 = (Y[:, i] for i in range(7, 11))
    D = Y[:, IX["D"]]
    Puma, Bax, MOMP, CytC, Smac, C9, C3 = (Y[:, i] for i in range(12, 19))
    Exec = Y[:, IX["Exec"]]
    drugs = _as_drug_tuple(drug)
    # Cin for drug i lives at N_BASE + i; with no drug the single
    # legacy slot is still present and simply stays at zero.
    Cin_all = [Y[:, N_BASE + i] for i in range(max(1, len(drugs)))]
    Cin = Cin_all[0]
    phase = _phase(Y)

    # ── cell cycle (prototype CDKState::step, growth signal = 1) ──
    gs = 1.0
    p21e = np.minimum(1.0, p21 * 1.5)
    CDK4act = CycD * (1 - p21e * 0.5)
    CDK2Eact = CycE * (1 - p21e)
    RbP = 1 - Rb
    APC_Cdh1 = np.maximum(0, 1 - (CDK2Eact + CycA) * 1.2)
    APC_Cdc20 = np.maximum(0, CycB - 0.25) * 2.0
    Cdc25 = np.maximum(0, CycA - 0.20) * 2.0
    dY[:, 0] = gs * 1.95 - CycD * (0.90 + (1 - gs) * 0.45)
    dY[:, 1] = 0.38 * (1 - Rb) - (CDK4act * 2.70 + CDK2Eact * 1.80) * Rb
    dY[:, 2] = RbP * 2.25 * (1 + E2F * 1.2) - (Rb * 1.80 + 0.45) * E2F
    dY[:, 3] = (E2F * 2.00 * (1 - CycA) - CycE * (0.65 + CycA * 2.45)) * (1 - p21e * 0.8)
    dY[:, 4] = (E2F * 1.80 * np.maximum(0, CycE - 0.12)
                - CycA * (0.22 + APC_Cdh1 * 1.60 + APC_Cdc20)) * (1 - p21e * 0.6)
    dY[:, 5] = (2.50 * Cdc25 * (1 + CycB * 0.8) - CycB * (0.18 + APC_Cdc20 * 2.5)) * (1 - p21e * 0.3)
    dY[:, 6] = -p21 * 0.18
    dY[:, :7] *= k_cyc
    # Growth signal per cell (0-1), set by the spatial dish for contact
    # inhibition and hypoxia; absent in a well-mixed population. It scales
    # progression through G1 only: S, G2 and M run their fixed course, so
    # a cell already past the restriction point finishes its cycle and its
    # daughters wait in G1 (quiescence) until the signal returns. Simply
    # lowering cyclin D synthesis would not do it, because this cycle's
    # Rb/E2F/cyclin E feedback restarts itself without cyclin D. PhysiCell
    # modulates the same G1 -> S transition with oxygen (Ghaffarizadeh et
    # al. 2018 PLoS Comput Biol 14:e1005991).
    g1_signal = aux.get("gs")
    if g1_signal is not None:
        dY[:, :7] *= np.where(phase == 0, g1_signal, 1.0)[:, None]
    cycle_rate = aux.get("cycle_rate")
    if cycle_rate is not None:
        dY[:, :7] *= cycle_rate[:, None]
    # p53 -> p21 (El-Deiry 1993) + p53-independent CHK1/2 S/G2 term
    dY[:, 6] += p.p21_k_max * _hill8(p53, p.p21_K_hill) * f53
    chk = (D > p.chk_damage_threshold) & ((phase == 1) | (phase == 2))
    dY[:, 6] += np.where(chk, D * p.chk_gain, 0.0)

    # ── DNA-damage response ──
    s = p.p53_time_scale
    atm_target = 1.0 / (1.0 + np.exp(-p.atm_gain * (D - p.atm_threshold)))
    dY[:, 7] = (atm_target - ATM) / p.atm_tau_h
    eff_mdm2 = MDM2 * (1.0 - ATM * p.atm_shield)
    dY[:, 8] = s * (p.p53_k_basal - p.p53_k_deg * eff_mdm2 / (p.p53_Kd + p53) * p53)
    dY[:, 9] = s * (p.mdm2_k_basal + p.mdm2_k_max * _hill8(p53, p.mdm2_K_hill) * f53
                    - p.mdm2_k_mrna_deg * mRNA)
    dY[:, 10] = s * (p.mdm2_k_translate * mRNA - p.mdm2_k_prot_deg * MDM2)

    # ── damage index ──
    dD = -line.repair_rate_per_h * D
    for i, dg in enumerate(drugs):
        ci = Cin_all[i]
        if dg.mechanism == "dna_adduct":
            dD = dD + dg.k_damage_per_uM_h * ci
        elif dg.mechanism == "topo2":
            dD = dD + dg.k_damage_per_uM_h * ci * np.where(phase == 1, dg.s_phase_factor, 1.0)
    forced = aux["forced_D"]
    dY[:, 11] = np.where(np.isnan(forced), dD, 0.0)

    # ── apoptosis ──
    p53_act = np.clip((p53 - p.p53_act_lo) / (p.p53_act_hi - p.p53_act_lo), 0, 1) * f53
    p73_act = np.clip((ATM - p.p73_atm_on) / (1.0 - p.p73_atm_on), 0, 1)
    dY[:, 12] = p.puma_k_on * (p53_act + p.p73_gain * p73_act) - p.puma_k_off * Puma
    mit_drive = np.zeros(len(Y))
    if any(dg.mechanism == "tubulin" for dg in drugs):
        k = max(dg.k_mitotic_death_per_h for dg in drugs if dg.mechanism == "tubulin")
        mit_drive = k * np.maximum(0.0, aux["arrest_h"] - p.mitotic_death_lag_h) / 10.0
    bcl2_eff = np.maximum(0.0, line.bcl2_level * aux["het_bcl2"]
                          - p.bcl2_neutralised_by_puma * Puma)
    dY[:, 13] = p.bax_k * (p.bax_puma_gain * Puma + mit_drive) / (1 + p.bcl2_block * bcl2_eff) * (1 - Bax)
    # Only Bax in excess of the anti-apoptotic buffer can form pores.
    bax_free = np.maximum(0.0, Bax - p.bcl2_buffer * bcl2_eff)
    bn = bax_free ** p.momp_n
    dY[:, 14] = p.momp_k * bn / (p.momp_K ** p.momp_n + bn) * (1 - MOMP)
    dY[:, 15] = p.release_k * MOMP * (1 - CytC)
    dY[:, 16] = p.release_k * MOMP * (1 - Smac)
    xiap_free = np.maximum(0.0, p.xiap_total - p.smac_total * Smac)
    dY[:, 17] = p.c9_k * CytC * (1 - C9) / (1 + p.c9_xiap_block * xiap_free)
    dY[:, 18] = p.c3_k * C9 * (1 - C3) * (1 - xiap_free / (0.1 + xiap_free))
    dY[:, 19] = p.exec_k * C3 * (1 - Exec)

    # ── drug uptake ──
    C_out_all = aux["C_out"]
    if C_out_all.ndim == 1:                      # single drug: (n_cells,)
        C_out_all = C_out_all[:, None]
    for i, dg in enumerate(drugs):
        ci = Cin_all[i]
        col = N_BASE + i
        dY[:, col] = ((dg.partition * aux["het_uptake"] * C_out_all[:, i] - ci)
                      / dg.tau_uptake_h)
        if p.efflux_vmax_uM_per_h > 0 and getattr(dg, "pgp_substrate", False):
            vmax = p.efflux_vmax_uM_per_h * line.efflux_level * aux["het_uptake"]
            dY[:, col] -= vmax * ci / (p.efflux_km_uM + ci)
    return dY


# ── population simulation ─────────────────────────────────────────────
@dataclass
class Population:
    Y: np.ndarray
    weight: np.ndarray
    alive: np.ndarray
    in_M: np.ndarray
    t_in_M: np.ndarray
    arrest_h: np.ndarray
    generation: np.ndarray
    C_out: np.ndarray
    forced_D: np.ndarray
    het_bcl2: np.ndarray             # per-cell Bcl-2 reserve multiplier (mean 1)
    het_uptake: np.ndarray           # per-cell drug-accumulation multiplier (mean 1)
    t_h: float = 0.0
    cycle_rate: Optional[np.ndarray] = None   # per-cell cycle speed (None = all 1)


@dataclass
class SimResult:
    t_h: np.ndarray                   # time axis
    alive_weight: np.ndarray          # (T, n_groups) summed alive weight per group
    groups: np.ndarray                # group id per cell (e.g. concentration index)
    group_C_out: np.ndarray           # C_out per group
    mean_p53: np.ndarray              # (T,) over alive cells (group 0)
    phase_frac: np.ndarray            # (T, 4) alive-weighted, group 0
    deaths: np.ndarray                # (T, n_groups) cumulative death weight
    final: Population

    def viability(self) -> np.ndarray:
        """Final alive weight of each group relative to group 0."""
        w = self.alive_weight[-1]
        return w / max(w[0], 1e-12)


def calibrate_cycle_scale(line: CellLine, p: Params = Params(), dt_h: float = 0.01,
                          max_h: float = 400.0) -> float:
    """Rate scale k_cyc so that one cycle (fresh G1 -> mitosis + M) takes
    the line's doubling time. One cell, no drug, RK4 at k_cyc=1."""
    Y = np.zeros((1, state_width(0)))
    Y[0, :7] = fresh_cycle_state()
    Y[0, IX["p53"]], Y[0, IX["MDM2_mRNA"]], Y[0, IX["MDM2"]] = 0.089, 0.21, 0.21
    aux = {"C_out": np.zeros(1), "arrest_h": np.zeros(1), "forced_D": np.full(1, np.nan),
           "het_bcl2": np.ones(1), "het_uptake": np.ones(1)}
    t = 0.0
    while t < max_h:
        _rk4(Y, aux, line, None, p, 1.0, dt_h)
        t += dt_h
        if Y[0, IX["CycB"]] > 0.25 and Y[0, IX["p21"]] < 0.35 and Y[0, IX["CycA"]] > 0.30:
            period_at_unit_scale = t + p.m_duration_h  # M itself is a fixed delay
            # k_cyc scales interphase only; solve k so interphase/k + M = T_d
            return t / max(line.doubling_time_h - p.m_duration_h, 1e-6)
    raise RuntimeError("cycle never reached mitosis at unit scale; check cycle ODE")


def _rk4(Y: np.ndarray, aux: dict, line, drug, p: Params, k_cyc: float, dt: float) -> None:
    k1 = rhs(Y, aux, line, drug, p, k_cyc)
    k2 = rhs(Y + 0.5 * dt * k1, aux, line, drug, p, k_cyc)
    k3 = rhs(Y + 0.5 * dt * k2, aux, line, drug, p, k_cyc)
    k4 = rhs(Y + dt * k3, aux, line, drug, p, k_cyc)
    Y += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    # Upper bounds must match the state width, which grows with the number
    # of drugs carried; every extra slot is an unbounded concentration.
    np.clip(Y, 0.0, _upper_for(Y.shape[1] - N_BASE), out=Y)


def _lognormal_mean_one(rng: np.random.Generator, sigma: float, n: int) -> np.ndarray:
    if sigma <= 0:
        return np.ones(n)
    return np.exp(rng.normal(-0.5 * sigma ** 2, sigma, n))


def make_population(n_cells: int, C_out: Union[float, np.ndarray], rng: np.random.Generator,
                    forced_D: Optional[float] = None, het_sigma: float = 0.0,
                    n_drugs: int = 1, cycle_cv: float = 0.0) -> Population:
    # C_out is (n_cells,) for one drug, or (n_cells, n_drugs) for several.
    C = np.asarray(C_out, dtype=float)
    C = (np.broadcast_to(C, (n_cells,)).copy() if n_drugs <= 1
         else np.broadcast_to(C, (n_cells, n_drugs)).copy())
    Y = np.zeros((n_cells, state_width(n_drugs)))
    for i in range(n_cells):
        Y[i, :7] = random_cycle_state(rng)
    Y[:, IX["p53"]], Y[:, IX["MDM2_mRNA"]], Y[:, IX["MDM2"]] = 0.089, 0.21, 0.21
    Y[:, IX["D"]] = rng.random(n_cells) * 0.05
    fD = np.full(n_cells, np.nan if forced_D is None else float(forced_D))
    if forced_D is not None:
        Y[:, IX["D"]] = forced_D
    # Drawn after the cycle states so het_sigma=0 leaves every earlier
    # random draw, and therefore every earlier result, unchanged.
    het_bcl2 = _lognormal_mean_one(rng, het_sigma, n_cells)
    het_uptake = _lognormal_mean_one(rng, het_sigma, n_cells)
    # Drawn last, and only when asked for, so cycle_cv=0 leaves every
    # earlier random draw unchanged. The duration factor has mean 1; the
    # rate is its reciprocal.
    cycle_rate = None
    if cycle_cv > 0:
        s = np.sqrt(np.log(1.0 + cycle_cv ** 2))
        cycle_rate = 1.0 / np.exp(rng.normal(-0.5 * s * s, s, n_cells))
    return Population(Y=Y, weight=np.ones(n_cells), alive=np.ones(n_cells, bool),
                      in_M=np.zeros(n_cells, bool), t_in_M=np.zeros(n_cells),
                      arrest_h=np.zeros(n_cells), generation=np.zeros(n_cells, int),
                      C_out=C, forced_D=fD, het_bcl2=het_bcl2, het_uptake=het_uptake,
                      cycle_rate=cycle_rate)


def _cell_events(pop: "Population", drug, p: Params, dt_h: float) -> tuple[np.ndarray, np.ndarray]:
    """Mitosis, mitotic arrest, slippage and death after one integration
    step, for any set of cells.

    Dividing cells are reset to early G1 (both daughters start identical)
    and returned in `divide`; what a division means for the population is
    the caller's business — `simulate` doubles a representative cell's
    weight, the spatial dish places a daughter on the lattice. Dying cells
    are marked dead and returned in `die`."""
    Y = pop.Y
    ready = (Y[:, IX["CycB"]] > 0.25) & (Y[:, IX["p21"]] < 0.35) & (Y[:, IX["CycA"]] > 0.30)
    enter = pop.alive & ~pop.in_M & ready & (Y[:, IX["D"]] < p.g2m_damage_block)
    pop.in_M[enter] = True
    pop.t_in_M[enter] = 0.0
    pop.t_in_M[pop.in_M] += dt_h
    arrested = np.zeros(len(Y), bool)
    for i, dg in enumerate(_as_drug_tuple(drug)):
        if dg.mechanism != "tubulin":
            continue
        ci = Y[:, N_BASE + i]
        theta = ci / (ci + dg.Kd_tubulin_uM)
        arrested = arrested | (pop.in_M & (theta > dg.theta_arrest))
    pop.arrest_h[arrested] += dt_h
    pop.arrest_h[~pop.in_M] = 0.0
    divide = pop.alive & pop.in_M & ~arrested & (pop.t_in_M >= p.m_duration_h)
    if divide.any():
        Y[divide, :7] = fresh_cycle_state()
        Y[divide, IX["p21"]] = np.maximum(0.02, Y[divide, IX["p21"]] * 0.2)
        pop.generation[divide] += 1
        pop.in_M[divide] = False
        pop.t_in_M[divide] = 0.0
    slip = pop.alive & arrested & (pop.t_in_M >= p.slippage_h)
    if slip.any():
        Y[slip, :7] = fresh_cycle_state()
        Y[slip, IX["D"]] = np.minimum(2.0, Y[slip, IX["D"]] + p.slippage_damage)
        pop.in_M[slip] = False
        pop.t_in_M[slip] = 0.0
        pop.arrest_h[slip] = 0.0
    die = pop.alive & (Y[:, IX["C3"]] > p.c3_commit)
    if die.any():
        pop.alive[die] = False
        pop.in_M[die] = False
    return divide, die


def simulate(line: CellLine, drug: Optional[Drug], C_out_uM: Union[float, np.ndarray],
             *, t_end_h: float = 72.0, n_cells: int = 128, dt_h: float = 0.02,
             seed: int = 1, groups: Optional[np.ndarray] = None,
             forced_D: Optional[float] = None, p: Params = Params(),
             k_cyc: Optional[float] = None, record_every_h: float = 1.0,
             dose_fn: Optional[Callable[[float], float]] = None,
             on_record: Optional[Callable[["Population", float], None]] = None,
             traits: Optional[dict] = None) -> SimResult:
    """Integrate a population for t_end_h. `C_out_uM` may be a per-cell
    vector; `groups` labels cells (e.g. by concentration) for readouts.

    `dose_fn(t_h)`, if given, sets every cell's extracellular
    concentration (µM) before each step, which is how dosing schedules
    and wash-outs are expressed. `on_record(pop, t_h)` is called at every
    record tick with the live population (read-only use), which is what
    `cellsim.cell.stream` turns into a per-tick state stream for a UI.

    `traits`, if given, sets each cell's heterogeneity multipliers
    instead of drawing them: {"bcl2": (n_cells,), "uptake": (n_cells,)}.
    Mapping outcome against chosen traits measures selection exactly,
    where random draws only sample it."""
    rng = np.random.default_rng(seed)
    drugs = _as_drug_tuple(drug)
    n_drugs = max(1, len(drugs))
    pop = make_population(n_cells, C_out_uM, rng, forced_D, het_sigma=p.het_sigma,
                          n_drugs=n_drugs, cycle_cv=p.cycle_cv)
    if traits is not None:
        pop.het_bcl2[:] = np.asarray(traits["bcl2"], float)
        pop.het_uptake[:] = np.asarray(traits["uptake"], float)
    if groups is None:
        groups = np.zeros(n_cells, int)
    n_groups = int(groups.max()) + 1
    _c1 = pop.C_out if pop.C_out.ndim == 1 else pop.C_out[:, 0]
    group_C = np.array([_c1[groups == g].mean() if np.any(groups == g) else np.nan
                        for g in range(n_groups)])
    if k_cyc is None:
        k_cyc = calibrate_cycle_scale(line, p)
    aux = {"C_out": pop.C_out, "arrest_h": pop.arrest_h, "forced_D": pop.forced_D,
           "het_bcl2": pop.het_bcl2, "het_uptake": pop.het_uptake}
    if pop.cycle_rate is not None:
        aux["cycle_rate"] = pop.cycle_rate
    n_steps = int(round(t_end_h / dt_h))
    rec_stride = max(1, int(round(record_every_h / dt_h)))
    T = n_steps // rec_stride + 1
    t_axis = np.zeros(T)
    alive_w = np.zeros((T, n_groups))
    deaths = np.zeros((T, n_groups))
    mean_p53 = np.zeros(T)
    phase_frac = np.zeros((T, 4))
    dead_weight = np.zeros(n_groups)

    def record(k: int) -> None:
        t_axis[k] = pop.t_h
        aw = pop.weight * pop.alive
        alive_w[k] = np.bincount(groups, weights=aw, minlength=n_groups)
        deaths[k] = dead_weight
        g0 = (groups == 0) & pop.alive
        mean_p53[k] = pop.Y[g0, IX["p53"]].mean() if g0.any() else np.nan
        ph = _phase(pop.Y)
        ph = np.where(pop.in_M, 3, ph)
        for q in range(4):
            phase_frac[k, q] = aw[g0 & (ph == q)].sum() / max(aw[g0].sum(), 1e-12)
        if on_record is not None:
            on_record(pop, pop.t_h)

    def _apply_dose(t_h: float) -> None:
        v = dose_fn(t_h)
        if pop.C_out.ndim == 1:
            pop.C_out[:] = v
        else:
            pop.C_out[:, :] = np.asarray(v, dtype=float)

    if dose_fn is not None:
        _apply_dose(0.0)
    record(0)
    for step in range(1, n_steps + 1):
        if dose_fn is not None:
            _apply_dose(pop.t_h)               # aux["C_out"] is this same array
        _rk4(pop.Y, aux, line, drug, p, k_cyc, dt_h)
        pop.t_h += dt_h
        divide, die = _cell_events(pop, drug, p, dt_h)
        if divide.any():
            pop.weight[divide] *= 2.0          # a representative cell stands for both daughters
        if die.any():
            dead_weight += np.bincount(groups[die], weights=pop.weight[die], minlength=n_groups)
        if step % rec_stride == 0:
            record(step // rec_stride)
    return SimResult(t_axis, alive_w, groups, group_C, mean_p53, phase_frac, deaths, pop)


def dose_response(line: CellLine, drug: Drug, conc_uM: np.ndarray, *, t_end_h: float = 72.0,
                  n_cells_per_conc: int = 96, dt_h: float = 0.02, seed: int = 1,
                  p: Params = Params(), k_cyc: Optional[float] = None,
                  exposure_h: Optional[float] = None) -> tuple[np.ndarray, SimResult]:
    """Viability at t_end_h for each concentration (index 0 is the
    untreated control, prepended automatically). One simulation.

    `exposure_h`, if given, washes the drug out at that time and reads
    the plate at t_end_h, which is how exposure-duration studies are run
    (e.g. a 3 h or 24 h pulse read at 72 h)."""
    conc = np.concatenate([[0.0], np.asarray(conc_uM, float)])
    groups = np.repeat(np.arange(len(conc)), n_cells_per_conc)
    C_out = conc[groups]
    dose_fn = None
    if exposure_h is not None:
        none = np.zeros_like(C_out)
        dose_fn = lambda t_h: C_out if t_h < exposure_h - 1e-9 else none  # noqa: E731
    res = simulate(line, drug, C_out, t_end_h=t_end_h, n_cells=len(groups), dt_h=dt_h,
                   seed=seed, groups=groups, p=p, k_cyc=k_cyc, dose_fn=dose_fn)
    return res.viability(), res


def ic50_from_curve(conc_uM: np.ndarray, viability: np.ndarray) -> float:
    """Concentration at which viability (relative to control) crosses 0.5,
    by log-linear interpolation; +inf if never, lowest conc if always."""
    c = np.asarray(conc_uM, float)
    v = np.asarray(viability, float)
    order = np.argsort(c)
    c, v = c[order], v[order]
    if v[0] <= 0.5:
        return float(c[0])
    for i in range(1, len(c)):
        if v[i] <= 0.5 < v[i - 1]:
            lc0, lc1 = np.log(c[i - 1]), np.log(c[i])
            f = (v[i - 1] - 0.5) / max(v[i - 1] - v[i], 1e-12)
            return float(np.exp(lc0 + f * (lc1 - lc0)))
    return float("inf")


def ic50(line: CellLine, drug: Drug, *, guess_uM: float = 1.0, t_end_h: float = 72.0,
         n_cells_per_conc: int = 48, seed: int = 1, p: Params = Params(),
         k_cyc: Optional[float] = None, exposure_h: Optional[float] = None) -> float:
    """The engine's continuous-exposure IC50 (µM) for one line and drug.

    Two passes: a 17-point grid over eight decades around `guess_uM`
    (3.2-fold steps) to bracket the crossing, then 13 points across that
    bracket (1.1-fold steps). Experiments use it to state doses as
    multiples of a drug's own potency, which keeps a design meaningful
    when the library's fitted constants are refitted. `exposure_h` gives
    the IC50 of a pulse washed out at that time and read at t_end_h; it
    is +inf when no concentration reaches half kill."""
    if k_cyc is None:
        k_cyc = calibrate_cycle_scale(line, p)
    coarse = guess_uM * np.logspace(-4, 4, 17)
    viab, _ = dose_response(line, drug, coarse, t_end_h=t_end_h, n_cells_per_conc=n_cells_per_conc,
                            seed=seed, p=p, k_cyc=k_cyc, exposure_h=exposure_h)
    first = ic50_from_curve(coarse, viab[1:])
    if not np.isfinite(first) or first <= coarse[0]:
        return first
    i = int(np.searchsorted(coarse, first))
    fine = np.geomspace(coarse[max(0, i - 1)], coarse[min(len(coarse) - 1, i)], 13)
    viab, _ = dose_response(line, drug, fine, t_end_h=t_end_h, n_cells_per_conc=n_cells_per_conc,
                            seed=seed, p=p, k_cyc=k_cyc, exposure_h=exposure_h)
    refined = ic50_from_curve(fine, viab[1:])
    return refined if np.isfinite(refined) else first


# ── CLI ───────────────────────────────────────────────────────────────
def main(argv: Optional[list[str]] = None) -> int:
    import argparse
    from cellsim.cell.library import get_drug, get_line

    ap = argparse.ArgumentParser(prog="cellsim cell-sim",
                                 description="Phase-1 single-cell drug-response engine")
    ap.add_argument("--line", default="A549")
    ap.add_argument("--drug", default="cisplatin")
    ap.add_argument("--conc", default="0.01,0.03,0.1,0.3,1,3,10,30,100",
                    help="µM, comma-separated (control added automatically)")
    ap.add_argument("--hours", type=float, default=72.0)
    ap.add_argument("--cells", type=int, default=96, help="representative cells per concentration")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(argv)
    line, drug = get_line(a.line), get_drug(a.drug)
    conc = np.array([float(x) for x in a.conc.split(",") if x.strip()])
    viab, res = dose_response(line, drug, conc, t_end_h=a.hours, n_cells_per_conc=a.cells, seed=a.seed)
    print(f"{line.name} + {drug.name}, {a.hours:.0f} h, {a.cells} cells/conc, seed {a.seed}")
    print(f"control growth: x{res.alive_weight[-1, 0] / res.alive_weight[0, 0]:.2f} "
          f"(doubling {line.doubling_time_h} h)")
    print(f"{'conc_uM':>10} {'viability':>10}")
    for c, v in zip(res.group_C_out, viab):
        print(f"{c:>10.4g} {v:>10.3f}")
    print(f"IC50 = {ic50_from_curve(res.group_C_out[1:], viab[1:]):.4g} µM")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
