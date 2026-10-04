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
    # p53 degradation that MDM2 inhibitors cannot block, per hour in the
    # engine's normalised units. HPV E6 hands p53 to the E6AP ligase; in
    # HPV-positive cervical cancer cells degradation switches COMPLETELY
    # from MDM2 to E6 (Hengstermann et al. 2001 PNAS 98:1218), which is
    # why their wild-type TP53 neither rises under nutlin nor mounts a
    # full damage response. 0 for HPV-negative lines.
    p53_mdm2_independent_deg_per_h: float = 0.0
    gdsc_name: str = ""              # name in GDSC tables, when it differs
    source: str = ""


@dataclass(frozen=True)
class Drug:
    name: str
    # 'dna_adduct' | 'topo2' | 'tubulin' | 's_phase' (damage only while the
    # cell replicates: TOP1 poisons, antimetabolites) | 'mdm2' (blocks
    # MDM2-mediated p53 degradation)
    mechanism: str
    k_damage_per_uM_h: float = 0.0   # damage-index gain per µM intracellular per hour (FIT)
    s_phase_factor: float = 1.0      # extra damage in S phase (TopII poisons)
    Kd_tubulin_uM: float = 0.0       # tubulin-site dissociation constant
    theta_arrest: float = 0.3        # tubulin occupancy that blocks anaphase
    k_mitotic_death_per_h: float = 0.0  # Bax drive per hour of mitotic arrest (FIT)
    Kd_target_uM: float = 0.0        # target affinity for 'mdm2' (MDM2-p53 site)
    tau_uptake_h: float = 0.1        # passive permeation time constant
    partition: float = 1.0           # C_in / C_out at equilibrium
    # Is this drug a P-glycoprotein substrate? Doxorubicin and paclitaxel
    # are classical ones; cisplatin is not effluxed by ABCB1, which is
    # also what the GDSC correlation shows (ABCB1 significant for the
    # first two, absent for cisplatin).
    pgp_substrate: bool = False
    # Total drug a cell holds at equilibrium (free + bound) relative to the
    # medium. Used only by the spatial dish, where it is the sink that
    # slows penetration into tissue: every cell layer must fill before
    # drug passes deeper. Order-of-magnitude inputs, flagged per drug.
    accumulation_ratio: float = 1.0
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
                       efflux_level=2.08, gdsc_name="HCT-116",
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
    "HeLa": CellLine("HeLa", 31.0, True, "cervix",
                     efflux_level=9.79, p53_mdm2_independent_deg_per_h=10.0,
                     source="Doubling 1.3 d (Cellosaurus CVCL_0030, PubMed 29156801; DSMZ "
                            "~48 h); the Cell Tracking Challenge HeLa movie's own counts "
                            "double every 27-31 h (scripts/ctc_reference.py). An earlier "
                            "20 h here was the mean of its COMPLETE cycles, which a 46 h "
                            "movie biases short. TP53 wild-type but degraded by HPV18 "
                            "E6 independently of MDM2 (Hengstermann 2001); the E6 rate "
                            "holds p53 below its apoptotic threshold even with MDM2 "
                            "fully blocked"),

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

    # ── Added 2026-10 for the spatial dish: the line whose spheroids have
    # a measured oxygen consumption rate and diffusion limit (Grimes et al.
    # 2014 J R Soc Interface 11:20131124) and whose multicellular layers
    # are the standard drug-penetration model (Phillips et al. 1998 Br J
    # Cancer 77:2112). High ABCB1, like its sister line HCT-15.
    "DLD-1": CellLine("DLD-1", 25.0, False, "colon",
                      efflux_level=16.02,
                      source="TP53 p.S241F (CVCL_0248). Doubling time reported 15, 20, "
                             "25.3, 33 and 48 h (Cellosaurus); median 25 h used, "
                             "uncertainty ~1.8x"),
}


# ── Drugs ──────────────────────────────────────────────────────────────
# Mechanism constants are literature inputs; the single potency gain per
# drug is fitted by scripts/validate_gdsc.py (see fit_target).
DRUGS: dict[str, Drug] = {
    "cisplatin": Drug(
        "cisplatin", "dna_adduct",
        k_damage_per_uM_h=0.001181,   # FITTED: A549 GDSC1 IC50 9.77 uM (scripts/validate_gdsc.py)
        tau_uptake_h=0.5,             # slow uptake (CTR1 + passive), hours
        partition=1.0,
        accumulation_ratio=2.0,       # Pt accumulates only a few-fold over medium (order of magnitude)
        fit_target="k_damage_per_uM_h",
        source="Pt-DNA adducts form in proportion to intracellular Pt (Jamieson & Lippard 1999 "
               "Chem Rev 99:2467); adduct repair t1/2 of hours (NER) sets repair_rate"),
    "doxorubicin": Drug(
        "doxorubicin", "topo2",
        k_damage_per_uM_h=0.01749,    # FITTED: geo-mean of A549 + MCF7 fits to GDSC1 (validate_gdsc.py)
        s_phase_factor=3.0,           # TopII poison: DSBs mostly during replication
        pgp_substrate=True,           # classical P-gp substrate
        tau_uptake_h=0.5,
        partition=10.0,               # weak base + DNA intercalation: high intracellular accumulation
        accumulation_ratio=100.0,     # DNA-bound anthracycline, 10-100x (Gigli 1988); order of magnitude
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
        partition=0.212,              # FITTED: geo-mean of A549 + MCF7 fits to GDSC2 (validate_gdsc.py)
        # The 72 h potency is set by the arrest threshold: once tubulin
        # occupancy passes theta_arrest every dividing cell stalls in M, so
        # the death rate barely moves the IC50. What does move it is how
        # much free drug the cell holds, i.e. uptake against efflux and
        # sequestration, so that is the one constant fitted.
        accumulation_ratio=100.0,     # microtubule-bound, saturable (Kuh et al. 2000 JPET 293:761); order of magnitude
        fit_target="partition",
        source="Microtubule stabiliser -> SAC-dependent mitotic arrest; death or slippage after "
               "many hours (Gascoigne & Taylor 2008 Cancer Cell 14:111)"),

    # ── Phase 3 panel (October 2026): seven more drugs, chosen so that
    # every mechanism class the engine can represent is covered and each
    # has GDSC screens on our lines. Each carries ONE fitted constant,
    # fitted exactly as above; values marked PROVISIONAL are starting
    # points that scripts/validate_gdsc.py replaces.
    "etoposide": Drug(
        "etoposide", "topo2",
        k_damage_per_uM_h=0.001,      # PROVISIONAL, refitted on GDSC1
        s_phase_factor=3.0,           # TopII poison, as doxorubicin
        pgp_substrate=True,           # MDR1 substrate
        tau_uptake_h=0.5, partition=1.0, accumulation_ratio=2.0,
        fit_target="k_damage_per_uM_h",
        source="TopII poison without intercalation; S/G2 double-strand breaks (Nitiss 2009 "
               "Nat Rev Cancer 9:338); P-glycoprotein substrate"),
    "sn-38": Drug(
        "sn-38", "s_phase",
        k_damage_per_uM_h=0.01,       # PROVISIONAL, refitted on GDSC2
        tau_uptake_h=0.3, partition=1.0, accumulation_ratio=5.0,
        fit_target="k_damage_per_uM_h",
        source="Active metabolite of irinotecan; TOP1 cleavage complexes become DSBs when "
               "replication forks collide with them (Pommier 2006 Nat Rev Cancer 6:789). "
               "Effluxed mainly by ABCG2, not modelled"),
    "gemcitabine": Drug(
        "gemcitabine", "s_phase",
        k_damage_per_uM_h=0.01,       # PROVISIONAL, refitted on GDSC2
        tau_uptake_h=1.0,             # nucleoside transport + dCK phosphorylation
        partition=1.0, accumulation_ratio=10.0,
        fit_target="k_damage_per_uM_h",
        source="dFdCTP incorporation with masked chain termination stalls replication "
               "(Plunkett 1995 Semin Oncol 22:3). Metabolic activation and retention are not "
               "modelled; the fitted constant absorbs them"),
    "5-fluorouracil": Drug(
        "5-fluorouracil", "s_phase",
        k_damage_per_uM_h=0.0005,     # PROVISIONAL, refitted on GDSC2
        tau_uptake_h=0.5, partition=1.0, accumulation_ratio=1.0,
        fit_target="k_damage_per_uM_h",
        source="FdUMP inhibits thymidylate synthase, starving DNA synthesis (Longley 2003 Nat "
               "Rev Cancer 3:330). Its RNA-directed toxicity, which is not S-phase specific, "
               "is not modelled"),
    "docetaxel": Drug(
        "docetaxel", "tubulin",
        Kd_tubulin_uM=0.005,          # ~2x paclitaxel's microtubule affinity (Diaz & Andreu 1993)
        theta_arrest=0.3, k_mitotic_death_per_h=0.3,
        pgp_substrate=True, tau_uptake_h=0.2,
        partition=0.2,                # PROVISIONAL, refitted on GDSC2
        accumulation_ratio=100.0,
        fit_target="partition",
        source="Taxane microtubule stabiliser, mechanism as paclitaxel with higher affinity "
               "(Diaz & Andreu 1993 Biochemistry 32:2747); P-glycoprotein substrate"),
    "vinorelbine": Drug(
        "vinorelbine", "tubulin",
        Kd_tubulin_uM=0.01,           # high-affinity binding at microtubule ends, order of magnitude
        theta_arrest=0.3, k_mitotic_death_per_h=0.3,
        pgp_substrate=True, tau_uptake_h=0.2,
        partition=0.2,                # PROVISIONAL, refitted on GDSC2
        accumulation_ratio=50.0,
        fit_target="partition",
        source="Vinca alkaloid: suppresses microtubule dynamics at low nM and arrests cells "
               "in mitosis like a stabiliser (Jordan & Wilson 2004 Nat Rev Cancer 4:253); "
               "P-glycoprotein substrate"),
    "nutlin-3a": Drug(
        "nutlin-3a", "mdm2",
        Kd_target_uM=0.09,            # IC50 90 nM for the MDM2-p53 interaction (Vassilev 2004)
        tau_uptake_h=0.3,
        partition=1.0,                # PROVISIONAL, refitted on GDSC2
        accumulation_ratio=5.0,
        fit_target="partition",
        source="Occupies MDM2's p53 pocket, stabilising wild-type p53 without DNA damage "
               "(Vassilev et al. 2004 Science 303:844). No effect on mutant p53's targets, and "
               "none where p53 is degraded independently of MDM2 (HPV E6)"),
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
