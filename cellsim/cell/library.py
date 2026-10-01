"""Curated cell lines and drugs for the Phase-1 engine.

Every number here is an *input with provenance*, not a fitted value.
The only fitted quantities in Phase 1 are the per-drug potency gains in
`Drug` (marked `fit_target`), which `scripts/validate_gdsc.py` fits on
the p53-wild-type reference lines and then applies unchanged to every
other line, so the mutant-line predictions are genuine out-of-sample
tests.

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
    fit_target: str = ""             # which field is fitted, if any
    source: str = ""


# ── Cell lines ─────────────────────────────────────────────────────────
# Doubling times: ATCC product sheets / Cellosaurus "doubling time"
# field. TP53: Cellosaurus genome-ancestry + IARC TP53 somatic data.
CELL_LINES: dict[str, CellLine] = {
    "A549": CellLine("A549", 22.0, True, "lung (NSCLC)",
                     source="ATCC CCL-185 ~22 h; TP53 wild-type (Cellosaurus CVCL_0023)"),
    "MCF7": CellLine("MCF7", 29.0, True, "breast (ER+)",
                     source="ATCC HTB-22 ~29 h; TP53 wild-type (CVCL_0031)"),
    "HCT116": CellLine("HCT116", 18.0, True, "colon",
                       source="ATCC CCL-247 ~18 h; TP53 wild-type (CVCL_0291)"),
    "MDA-MB-231": CellLine("MDA-MB-231", 26.0, False, "breast (TNBC)",
                           source="ATCC HTB-26 ~25-38 h; TP53 p.R280K (CVCL_0062)"),
    "HT-29": CellLine("HT-29", 20.0, False, "colon",
                      source="ATCC HTB-38 ~20-23 h; TP53 p.R273H (CVCL_0320)"),
    "SW480": CellLine("SW480", 26.0, False, "colon",
                      source="ATCC CCL-228 ~26 h; TP53 p.R273H + p.P309S (CVCL_0546)"),
    "HeLa": CellLine("HeLa", 20.0, True, "cervix",
                     source="CTC Fluo-N2DL-HeLa ~20 h; TP53 wild-type but HPV18 E6-degraded "
                            "(functionally hypomorphic; modelled as functional, flagged)"),
}


# ── Drugs ──────────────────────────────────────────────────────────────
# Mechanism constants are literature inputs; the single potency gain per
# drug is fitted by scripts/validate_gdsc.py (see fit_target).
DRUGS: dict[str, Drug] = {
    "cisplatin": Drug(
        "cisplatin", "dna_adduct",
        k_damage_per_uM_h=0.05,       # placeholder; FIT on wt lines
        tau_uptake_h=0.5,             # slow uptake (CTR1 + passive), hours
        partition=1.0,
        fit_target="k_damage_per_uM_h",
        source="Pt-DNA adducts form in proportion to intracellular Pt (Jamieson & Lippard 1999 "
               "Chem Rev 99:2467); adduct repair t1/2 of hours (NER) sets repair_rate"),
    "doxorubicin": Drug(
        "doxorubicin", "topo2",
        k_damage_per_uM_h=0.5,        # placeholder; FIT
        s_phase_factor=3.0,           # TopII poison: DSBs mostly during replication
        tau_uptake_h=0.5,
        partition=10.0,               # weak base + DNA intercalation: high intracellular accumulation
        fit_target="k_damage_per_uM_h",
        source="TopII poison, S-phase-dependent DSBs (Nitiss 2009 Nat Rev Cancer 9:338); "
               "nuclear accumulation 10-100x (Gigli 1988 Cancer Res 48:4100)"),
    "paclitaxel": Drug(
        "paclitaxel", "tubulin",
        Kd_tubulin_uM=0.010,          # ~10 nM for microtubule sites (Diaz & Andreu 1993 Biochemistry)
        theta_arrest=0.3,
        k_mitotic_death_per_h=0.3,    # placeholder; FIT
        tau_uptake_h=0.2,
        partition=1.0,
        fit_target="k_mitotic_death_per_h",
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
