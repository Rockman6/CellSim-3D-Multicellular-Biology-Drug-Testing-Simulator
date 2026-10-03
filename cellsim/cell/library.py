"""Curated cell lines and drugs for the Phase-1 engine.

Every number here is an *input with provenance*, not a fitted value.
The only fitted quantities in Phase 1 are the per-drug constants named
by `fit_target` (values below, marked FITTED). `scripts/validate_gdsc.py`
fits each on the clean TP53 wild-type lines A549 and MCF7 against GDSC
release 8.4 (with the engine's default parameters, including the p73
route) and applies it unchanged to every other line, so HeLa and
both TP53-mutant lines are genuine out-of-sample tests. Re-running the
script reproduces these values to within 1 %.

Doubling times are lab-dependent (±30 % between reports is common);
the values below are the ATCC / Cellosaurus figures and are used only
to scale the cycle ODE so that an untreated population doubles at the
quoted rate. TP53 status follows the IARC TP53 database as summarised
in Cellosaurus.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class CellLine:
    name: str
    doubling_time_h: float
    p53_functional: bool
    tissue: str = ""
    bcl2_level: float = 1.0          # relative anti-apoptotic reserve (1 = reference)
    # Relative P-glycoprotein (ABCB1) activity, 1 = reference. ABCB1
    # expression is one of only two per-line signals that survive
    # correction against GDSC at n = 427-683 (see docs/VALIDATION.md),
    # and it is the ONLY channel through which a taxane's potency can
    # vary by line, because paclitaxel's IC50 is set by a tubulin
    # occupancy threshold that belongs to the drug.
    efflux_level: float = 1.0
    repair_rate_per_h: float = 0.0866  # DNA-adduct/DSB repair, ln2 / 8 h
    source: str = ""


@dataclass(frozen=True)
class Drug:
    name: str
    mechanism: str                   # 'dna_adduct' | 'topo2' | 'tubulin'
    k_damage_per_uM_h: float = 0.0   # damage-index gain per µM intracellular per hour (FIT)
    s_phase_factor: float = 1.0      # extra damage in S phase (TopII poisons)
    Kd_tubulin_uM: float = 0.0       # tubulin-site dissociation constant
    theta_arrest: float = 0.3        # tubulin occupancy that blocks anaphase
    k_mitotic_death_per_h: float = 0.0  # Bax drive per hour of mitotic arrest (FIT)
    tau_uptake_h: float = 0.1        # passive permeation time constant
    partition: float = 1.0           # C_in / C_out at equilibrium
    # Is this drug a P-glycoprotein substrate? Doxorubicin and paclitaxel
    # are classical ones; cisplatin is not effluxed by ABCB1, which is
    # also what the GDSC correlation shows (ABCB1 significant for the
    # first two, absent for cisplatin).
    pgp_substrate: bool = False
    fit_target: str = ""             # which field is fitted, if any
    source: str = ""


# `efflux_level` is derived, not fitted: 2^(ABCB1 log2(TPM+1) − panel
# median) from DepMap 24Q4, i.e. each line's (TPM+1) ratio to the median
# of this panel. It is a fixed formula with no free parameter, so the
# held-out lines are genuinely out of sample. Using (TPM+1) rather than
# TPM is deliberate and conservative: it compresses the low end, where
# TPM is near zero and unreliable, and understates the true spread.
# Regenerate with scripts/experiment_expression_scale.py's data files.

# ── Cell lines ─────────────────────────────────────────────────────────
# Doubling times: ATCC product sheets / Cellosaurus "doubling time"
# field. TP53: Cellosaurus genome-ancestry + IARC TP53 somatic data.
CELL_LINES: dict[str, CellLine] = {
    "A549": CellLine("A549", 22.0, True, "lung (NSCLC)",
                     efflux_level=0.96,
                     source="ATCC CCL-185 ~22 h; TP53 wild-type (Cellosaurus CVCL_0023)"),
    "MCF7": CellLine("MCF7", 29.0, True, "breast (ER+)",
                     efflux_level=0.95,
                     source="ATCC HTB-22 ~29 h; TP53 wild-type (CVCL_0031)"),
    "HCT116": CellLine("HCT116", 18.0, True, "colon",
                       efflux_level=2.08,
                       source="ATCC CCL-247 ~18 h; TP53 wild-type (CVCL_0291)"),
    "MDA-MB-231": CellLine("MDA-MB-231", 26.0, False, "breast (TNBC)",
                           efflux_level=0.97,
                           source="ATCC HTB-26 ~25-38 h; TP53 p.R280K (CVCL_0062)"),
    "HT-29": CellLine("HT-29", 20.0, False, "colon",
                      efflux_level=1.03,
                      source="ATCC HTB-38 ~20-23 h; TP53 p.R273H (CVCL_0320)"),
    "SW480": CellLine("SW480", 26.0, False, "colon",
                      efflux_level=1.69,
                      source="ATCC CCL-228 ~26 h; TP53 p.R273H + p.P309S (CVCL_0546)"),
    "HeLa": CellLine("HeLa", 20.0, True, "cervix",
                     efflux_level=9.79,
                     source="CTC Fluo-N2DL-HeLa ~20 h; TP53 wild-type but HPV18 E6-degraded "
                            "(functionally hypomorphic; modelled as functional, flagged)"),

    # ── Added 2026-10 to test the p53-independent death route against
    # lines that played no part in choosing its parameter. Doubling
    # times are the ONLY line-specific input the engine has, so where
    # sources disagree the spread is recorded: it bounds how well any
    # line-to-line prediction can do.
    "U-2-OS": CellLine("U-2-OS", 27.0, True, "osteosarcoma",
                       efflux_level=5.18,
                       source="DSMZ ACC 785 25-30 h; TP53 wild-type, the standard "
                              "p53-functional positive control (CVCL_0042)"),
    "HT-1080": CellLine("HT-1080", 26.0, True, "fibrosarcoma",
                        efflux_level=0.96,
                        source="DSMZ ACC 315 ~26 h; TP53 wild-type (CVCL_0317)"),
    "T47D": CellLine("T47D", 32.0, False, "breast (ER+)",
                     efflux_level=0.97,
                     source="ATCC HTB-133 ~32 h; TP53 p.L194F (CVCL_0553)"),
    "MDA-MB-468": CellLine("MDA-MB-468", 38.0, False, "breast (TNBC)",
                           efflux_level=2.25,
                           source="TP53 p.R273H (CVCL_0419). Doubling time is poorly "
                                  "agreed: DSMZ 30-40 h, others 40.6 h, NCI-DCTD 62 h; "
                                  "38 h used, uncertainty ~1.6x"),
    "MIA-PaCa-2": CellLine("MIA-PaCa-2", 30.0, False, "pancreas",
                           efflux_level=0.94,
                           source="TP53 p.R248W homozygous (CVCL_0428). Doubling time "
                                  "reported 26 h to 40 h depending on conditions; 30 h "
                                  "used, uncertainty ~1.5x"),
}


# ── Drugs ──────────────────────────────────────────────────────────────
# Mechanism constants are literature inputs; the single potency gain per
# drug is fitted by scripts/validate_gdsc.py (see fit_target).
DRUGS: dict[str, Drug] = {
    "cisplatin": Drug(
        "cisplatin", "dna_adduct",
        k_damage_per_uM_h=0.000951,   # FITTED: A549 GDSC1 IC50 9.77 uM (scripts/validate_gdsc.py)
        tau_uptake_h=0.5,             # slow uptake (CTR1 + passive), hours
        partition=1.0,
        fit_target="k_damage_per_uM_h",
        source="Pt-DNA adducts form in proportion to intracellular Pt (Jamieson & Lippard 1999 "
               "Chem Rev 99:2467); adduct repair t1/2 of hours (NER) sets repair_rate"),
    "doxorubicin": Drug(
        "doxorubicin", "topo2",
        k_damage_per_uM_h=0.00794,    # FITTED: geo-mean of A549 + MCF7 fits to GDSC1 (validate_gdsc.py)
        s_phase_factor=3.0,           # TopII poison: DSBs mostly during replication
        pgp_substrate=True,           # classical P-gp substrate
        tau_uptake_h=0.5,
        partition=10.0,               # weak base + DNA intercalation: high intracellular accumulation
        fit_target="k_damage_per_uM_h",
        source="TopII poison, S-phase-dependent DSBs (Nitiss 2009 Nat Rev Cancer 9:338); "
               "nuclear accumulation 10-100x (Gigli 1988 Cancer Res 48:4100)"),
    "paclitaxel": Drug(
        "paclitaxel", "tubulin",
        Kd_tubulin_uM=0.010,          # ~10 nM for microtubule sites (Diaz & Andreu 1993 Biochemistry)
        theta_arrest=0.3,
        pgp_substrate=True,           # classical P-gp substrate
        k_mitotic_death_per_h=0.3,
        tau_uptake_h=0.2,
        partition=0.153,              # FITTED: geo-mean of A549 + MCF7 fits to GDSC2 (validate_gdsc.py)
        # The 72 h potency is set by the arrest threshold: once tubulin
        # occupancy passes theta_arrest every dividing cell stalls in M, so
        # the death rate barely moves the IC50. What does move it is how
        # much free drug the cell holds, i.e. uptake against efflux and
        # sequestration, so that is the one constant fitted.
        fit_target="partition",
        source="Microtubule stabiliser -> SAC-dependent mitotic arrest; death or slippage after "
               "many hours (Gascoigne & Taylor 2008 Cancer Cell 14:111)"),
}


def get_line(name: str) -> CellLine:
    try:
        return CELL_LINES[name]
    except KeyError:
        raise KeyError(f"unknown cell line {name!r}; known: {sorted(CELL_LINES)}") from None


def get_drug(name: str) -> Drug:
    try:
        return DRUGS[name.lower()]
    except KeyError:
        raise KeyError(f"unknown drug {name!r}; known: {sorted(DRUGS)}") from None
