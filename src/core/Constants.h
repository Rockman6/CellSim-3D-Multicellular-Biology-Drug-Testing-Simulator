#pragma once

// ══════════════════════════════════════════════════════════════════════════
//  Constants.h — Biology constants calibrated against published data
//  Target: HeLa-like cells, ~20h doubling, G1:46% S:33% G2:17% M:4%
//  References: Novak & Tyson (2008), Casciari (1992), Green & Kroemer (2004)
// ══════════════════════════════════════════════════════════════════════════

// ── Simulation Modes ────────────────────────────────────────────────────
enum SimMode {
    MODE_SINGLE_CELL = 0,  // Single cell, no division, detailed view
    MODE_COLONY      = 1,  // Multi-cell colony growth (original behavior)
};

// ── World Scale ─────────────────────────────────────────────────────────
// EVERY distance constant in the simulator derives from this single anchor.
// 1 world unit = 5 µm places CELL_RADIUS_BASE=2.0 at exactly 10 µm, which is
// HeLa's median cell radius (20 µm diameter — Alberts MBoC Ch 1; BioNumbers
// BNID 100434; data/reference/hela_reference.csv). Change this only if
// every other distance constant is rescaled in lockstep. SCENE_BOUND = 55 wu
// at this scale represents a 275 µm × 550 µm field of view — deliberately
// a magnified sub-region of a real 10 cm Petri dish, not the whole plate.
constexpr float MICROMETERS_PER_WORLD_UNIT = 5.0f;
constexpr float UM_TO_WU = 1.0f / MICROMETERS_PER_WORLD_UNIT;

// ── Molecule physical radii (in µm) ─────────────────────────────────────
// VALUES ARE RADII, not diameters. Sources:
//   Alberts MBoC 7e Ch 6, 8, 14
//   BioNumbers BNID 111385 (ribosome), 103861 (nuclear pore)
//   Adrio 2017 (spliceosome), Cramer 2004 (RNA Pol II)
//   Melnik 2012 (Hsp70 crystal structure), Lentz 2010 (vesicle sizes)
// Note these are the PHYSICAL radii. The actual rendered sphere size is
// radius_um * UM_TO_WU * VIS_BOOST because a 5 nm particle inside a
// 10 µm cell is 0.05% of the diameter — invisible at normal camera
// distance. The boost is cosmetic only; physics uses the real radius
// indirectly through particle spacing rules, not the rendered size.
namespace MoleculeRadiusUm {
    // Nuclear machinery
    constexpr float DNA_NODE         = 0.005f;   // nucleosome (147 bp wrapped on histone octamer, ~10 nm ⌀)
    constexpr float RNA_POL          = 0.007f;   // Pol II complex 14 nm diameter
    constexpr float DNA_POL          = 0.006f;   // replisome 12 nm diameter
    constexpr float SPLICEOSOME      = 0.015f;   // 30 nm diameter
    constexpr float NUCLEAR_PORE     = 0.060f;   // NPC 120 nm diameter
    // Cytoplasmic macromolecules
    constexpr float RIBOSOME_SMALL   = 0.010f;   // 40S ~20 nm diameter
    constexpr float RIBOSOME_LARGE   = 0.0125f;  // 60S ~25 nm diameter
    constexpr float TRNA             = 0.0035f;  // tRNA ~7 nm long
    constexpr float CHAPERONE        = 0.0025f;  // Hsp70 ~5 nm diameter (average folded protein)
    // Membrane carriers
    constexpr float COPII_VESICLE    = 0.030f;   // COPII 60 nm diameter
    constexpr float SECRETORY_VESICLE = 0.050f;  // 100 nm diameter
    // Small molecules — radii from compiled E. coli vdW table (Ouzounis 1998
    // / glycolysis compendium) plus Alberts Ch 2. Hydrodynamic radii used
    // for moving species where available.
    constexpr float ATP              = 0.0007f;  // ~0.7 nm hydrodynamic radius
    constexpr float GLUCOSE          = 0.0010f;  // ~1.0 nm vdW radius
    constexpr float NADH             = 0.0040f;  // ~4.0 nm vdW (large dinucleotide cofactor)
    constexpr float AMINO_ACID       = 0.00035f; // glycine ~0.7 nm long → ~0.35 nm half-extent
    constexpr float CAMP             = 0.0007f;  // nucleotide, matches ATP scale
    constexpr float WATER            = 0.00014f; // water vdW radius 1.4 Å (0.28 nm diameter)
    constexpr float CALCIUM_ION      = 0.0001f;  // Ca²⁺ bare ionic radius ~1 Å
}

// Render-only visibility boosts. Physical particle radius in a 10 µm cell
// is sub-pixel at normal camera distance, so we scale up for visibility.
// Lowered values so the sea of tiny metabolite particles (ATP/glucose/
// water/AA/NADH — ~1000 per cell in the sim) reads as soft texture
// rather than a visual storm that drowns the organelles.
constexpr float MOL_VIS_BOOST_MACROMOL = 3.0f;   // ribosome, pores, spliceosome
constexpr float MOL_VIS_BOOST_SMALL    = 15.0f;  // ATP, ions, small molecules

// Convert a physical molecular radius (µm) into render-space world units
// with the appropriate visibility boost. Use this at every particle-add
// site so the "why is this 0.012?" answer is a one-line constant.
static inline float molRadiusWu(float radius_um, float vis_boost) {
    return radius_um * UM_TO_WU * vis_boost;
}

// ── Organelle copy number by cycle phase ────────────────────────────────
// Mitochondria fission doubles the network through S/G2 (Mishra & Chan
// 2014 Nat Rev MCB). Golgi fragments into vesicles for inheritance at
// M-phase (Warren & Wickner 1996). These counts feed the render pipeline
// so organelle instances grow before division rather than appearing
// fully formed on both daughters from the first frame.
constexpr int MITO_COUNT_G1 = 3;
constexpr int MITO_COUNT_G2 = 6;   // fission: 3 → 6 during S/G2
constexpr int GOLGI_COUNT_G1 = 1;
constexpr int GOLGI_COUNT_G2 = 2;  // Golgi fragmentation for mitosis

// ── Scene ────────────────────────────────────────────────────────────────
constexpr float FLOOR_Y      = -7.0f;
constexpr float SCENE_BOUND  = 55.0f;
constexpr int   MAX_CELLS    = 2000;  // hard cap for memory
// Plate carrying capacity: π×55²/(π×2.15²) × 0.9 hex packing ≈ 590 monolayer cells
// Growth should plateau around this number via contact inhibition
constexpr int   MONOLAYER_CAPACITY = 600;
constexpr int   INIT_CELLS   = 5;
constexpr int   SINGLE_CELL_COUNT = 1;

// ── Timescales ───────────────────────────────────────────────────────────
// SINGLE SOURCE OF TRUTH for bio-time → wall-time mapping.
// At `timeScale = 1.0`, 1 wall-second = 180 bio-seconds = 3 bio-minutes.
// A full HeLa 24-hour cycle therefore takes 86400/180 = 480 wall-seconds
// (8 wall-minutes) at default speed. The user-facing speed slider scales
// this further: timeScale=60 → 1 cycle = 8 wall-seconds.
//
// Every per-bio-second rate constant in this file is multiplied by
// `bioDt(real_dt, timeScale)` at its call site. There is no second clock
// — `SLOW_DT_SCALE`, `MEDIUM_DT_SCALE`, `CDK_DT_SCALE`, `BIO_MIN_PER_SEC`
// are kept as plain numbers so the legacy code paths still compile while
// they are migrated, but new code MUST use `bioDt()`.
constexpr float BIO_TIME_SCALE  = 180.0f;  // bio-seconds per real-second
constexpr float BIO_MIN_PER_SEC = 3.0f;    // legacy alias; = BIO_TIME_SCALE/60

inline float bioDt(float real_dt, float timeScale) {
    return real_dt * timeScale * BIO_TIME_SCALE;
}

// Legacy ratios — kept so existing call sites compile during migration.
// New code must call `bioDt()` and use `*_PER_BIOSEC` constants.
//
// Anchor: at `timeScale = 1.0` we want the full HeLa cycle (24 bio-h) to
// take 480 wall-seconds (8 wall-min). The legacy `cycleTimer` accumulates
// in "slow units" with `sum(CYCLE_*_DUR) = 6.5`. So one slow-unit
// represents 480/6.5 ≈ 73.8 wall-seconds at timeScale=1, which means
// `SLOW_DT_SCALE = 6.5 / 480 = 0.01354 sdt per real-second`. Lowering
// from the original 0.06 stretches the cycle 4.4× — matching the bio-time
// target. CDK and metabolism rates run on `sdt` and `mdt`, so they
// automatically slow down in lockstep, preserving every previously tuned
// rate ratio between metabolism and cycle.
constexpr float FAST_DT_SCALE   = 1.0f;
constexpr float MEDIUM_DT_SCALE = 0.0497f;   // was 0.22, scaled to match
constexpr float SLOW_DT_SCALE   = 0.01354f;  // was 0.06; gives 24 bio-h cycle
constexpr float CDK_DT_SCALE    = 0.8f;      // unchanged: ratio kept

// ── Lag phase ────────────────────────────────────────────────────────────
// Newborn-cell adaptation interval before resuming full cycling. Real
// freshly-plated HeLa lag is 1–6 h (AAT Bioquest); we use 30 bio-min
// for newly-divided daughters so they don't immediately re-enter S.
constexpr float LAG_DURATION    = 0.3f;            // legacy slow-time units
constexpr float LAG_BIOSEC      = 30.0f * 60.0f;   // 30 bio-min

// ── Diffusion (Fick's 2nd law) ───────────────────────────────────────────
// Legacy normalized coefficients kept for the not-yet-replaced
// NutrientField path. New code uses MediumField D_um2_per_s table.
constexpr int   NUTRIENT_GRID_N   = 64;
constexpr float DIFF_O2_COEFF     = 1.8f;
constexpr float DIFF_GLC_COEFF    = 0.9f;
constexpr float DIFF_CO2_COEFF    = 1.5f;
constexpr float DIFF_PH_COEFF     = 0.6f;
constexpr float O2_CONSUME_BASE   = 0.0075f;
constexpr float GLC_CONSUME_BASE  = 0.0040f;
constexpr float CO2_PRODUCE_BASE  = 0.0060f;
constexpr float LACTATE_PRODUCE   = 0.0030f;

// ── Culture-medium composition (DMEM high-glucose, mM) ───────────────────
// Bake-once initial concentrations for the closed-system MediumField.
// Sources: Sigma-Aldrich D5796 datasheet, AAT Bioquest DMEM-HG recipe,
// R&D Systems M17250. FBS 10 % v/v supplements the growth-factor term.
// HeLa is the canonical reference; same recipe drives every cell line
// in the simulator.
namespace MediumComposition {
    // Energy substrates
    constexpr float DMEM_GLUCOSE_MM      = 25.0f;
    constexpr float DMEM_GLUTAMINE_MM    = 4.0f;
    constexpr float DMEM_PYRUVATE_MM     = 1.0f;
    // Pooled essential amino acids (Leu+Lys+Arg+Ile+Val+Thr+Phe+Met+His+
    // Trp+Tyr+Cys+Gly+Ser, summed to a single transport / synthesis pool)
    constexpr float DMEM_AA_POOL_MM      = 5.0f;
    // Dissolved gases at incubator equilibrium (37 °C, 18 % O₂, 5 % CO₂)
    constexpr float DMEM_O2_MM           = 0.20f;
    constexpr float DMEM_CO2_MM          = 1.20f;
    // Waste (starts at zero)
    constexpr float DMEM_LACTATE_MM      = 0.0f;
    // Bulk ions (Na+ + K+ + Cl- lumped — set/forget for osmotic balance)
    constexpr float DMEM_IONS_MM         = 155.0f;
    constexpr float DMEM_CALCIUM_MM      = 1.8f;
    // Acid-base
    constexpr float DMEM_PH              = 7.40f;
    // Growth factors (10 % FBS contribution: ~11 ng/mL IGF + ~3.7 ng/mL
    // FGF + EGF/TGF traces lumped together)
    constexpr float DMEM_GROWTH_FACTOR_NG_PER_ML = 50.0f;
    // Drug (always starts clean)
    constexpr float DMEM_DRUG_UM         = 0.0f;
    // Water — DMEM is mostly water; pure water is ~55.5 M, dissolved
    // species reduce this slightly. 55 000 mM ≈ 55 M (BNID 109506).
    constexpr float DMEM_WATER_MM        = 55000.0f;
}

// Diffusion coefficients in DILUTE WATER at 37 °C (µm² / s).
// Source list in plan file. Cytoplasmic values are 26 % of these for
// small molecules (Luby-Phelps 2000) — applied at the cell-particle
// site, not at the medium grid (which stays bulk-water).
namespace MediumDiffusion {
    constexpr float D_O2_UM2_PER_S        = 2100.0f;   // Goldstick 1976
    constexpr float D_CO2_UM2_PER_S       = 1900.0f;
    constexpr float D_GLUCOSE_UM2_PER_S   =  670.0f;   // Longsworth 1953
    constexpr float D_GLUTAMINE_UM2_PER_S =  760.0f;
    constexpr float D_PYRUVATE_UM2_PER_S  = 1100.0f;
    constexpr float D_LACTATE_UM2_PER_S   = 1030.0f;   // Stein 1986
    constexpr float D_AA_UM2_PER_S        =  800.0f;   // amino-acid avg
    constexpr float D_GROWTH_F_UM2_PER_S  =   80.0f;   // ~25 kDa Stokes-Einstein
    constexpr float D_IONS_UM2_PER_S      = 1400.0f;   // Na/K/Cl avg
    constexpr float D_CALCIUM_UM2_PER_S   =  790.0f;
    constexpr float D_H_UM2_PER_S         = 9300.0f;   // Grotthuss; effectively pH
    constexpr float D_WATER_UM2_PER_S     = 2200.0f;   // self-diffusion in water 25 °C
    constexpr float CYTOPLASM_SLOWDOWN    = 0.26f;     // Luby-Phelps 2000
}

// ── Osmosis & turgor (Verkman 2008, Stewart 2011, Sten-Knudsen 2002) ────
// Aquaporin-mediated water flux follows van't Hoff via osmotic pressure
// gradient. Lp expressed per-aquaporin per-bio-second so the per-cell rate
// is `Lp * aquaporinCount * Δosmolarity_mM`.
namespace Osmosis {
    // Per-aquaporin permeability — calibrated so 10⁵ aquaporins move
    // ~10 mM water per bio-second under a 10 mM gradient (~1 % volume
    // shift in 10 bio-min, roughly matching osmotic shock kinetics).
    constexpr float AQP_LP_PER_BIOSEC      = 1.0e-8f;   // mM_water / (count × mM_osmo × bio-s)
    constexpr float OSMO_ISOTONIC_MM       = 290.0f;    // cytosolic baseline
    // Bulk modulus mapping water inflow to turgor pressure.
    // 5 % volume change → 5 kPa pressure change.
    constexpr float BULK_MOD_PA            = 1.0e5f;
    constexpr float TURGOR_LYSE_PA         = 5.0e4f;    // safety cap
    // Default aquaporin count when not yet enumerated from membrane mosaic.
    constexpr int   AQUAPORIN_COUNT_DEFAULT = 200;
}

// ── Medium chemical particle binding / transport timing (visual only) ───
// The state-machine animation that takes a free medium particle, binds it
// to a membrane transporter, then streams it into the cell interior.
namespace MediumBinding {
    constexpr float UPTAKE_RADIUS_WU       = 0.6f;      // distance to start binding
    constexpr float BINDING_BIOSEC         = 60.0f;     // pulse duration (~1 bio-min)
    constexpr float TRANSPORT_BIOSEC       = 120.0f;    // membrane → interior animation
}

// ── Substrate adhesion (integrin / focal adhesion timeline) ─────────────
// Maturation timeline drawn from Cavalcanti-Adam 2007, Geiger 2001,
// Reinhart-King 2005, Gallant 2005. State machine durations in bio-sec.
namespace Adhesion {
    constexpr float FLOATING_BIOSEC          = 5.0f;             // 0–5 s after seed
    constexpr float CONTACT_BIOSEC           = 5.0f * 60.0f;     // 5 s – 5 min
    constexpr float SPREADING_BIOSEC         = 60.0f * 60.0f;    // 5 min – 1 h
    // Parameters reached at each end of state machine.
    constexpr float SPREAD_FACTOR_FLOATING   = 1.00f;
    constexpr float SPREAD_FACTOR_MATURE     = 1.30f;
    constexpr float STRENGTH_FLOATING        = 0.00f;
    constexpr float STRENGTH_MATURE          = 0.95f;
    // Focal adhesion cluster size grows from 5 → 50 over maturation.
    constexpr int   FA_COUNT_FLOATING        = 0;
    constexpr int   FA_COUNT_MATURE          = 50;
    // Round-up (during mitosis or after detachment) decay constant.
    constexpr float ROUND_UP_BIOSEC          = 5.0f * 60.0f;     // 5 bio-min
}

// ── Metabolism ───────────────────────────────────────────────────────────
// BIOMASS_SYNTH_K bumped 0.010 → 0.014 so biomass actually reaches ~2.0 by
// G2 end — matches HeLa cell volume doubling (4000 µm³ → 8000 µm³ over one
// cycle; Alberts Ch 17). With the comprehensive M-entry checkpoint now
// gating on biomass ≥ 1.70, the old rate left cells stuck at ~1.5 and
// failing the gate forever.
// Biomass synth rate tuned so G1 biomass accumulation (1.0 → 1.30) takes
// ~8 bio-h at full nutrients, giving a 24-bio-h total cycle that matches
// HeLa doubling time. Prior value 0.014 produced a ~60 bio-h G1.
constexpr float BIOMASS_SYNTH_K       = 0.055f;
constexpr float BIOMASS_DEGRADE_K     = 0.0015f;
constexpr float BIOMASS_DIVIDE_THRESH = 1.85f;  // soft gate (only for quality check)
constexpr float ROS_MUTATION_SCALE    = 0.03f;

// ── Hypoxia (Gillies 2012 Nat Rev Cancer) ────────────────────────────────
constexpr float HYPOXIA_MODERATE     = 0.15f;
constexpr float HYPOXIA_SEVERE       = 0.08f;
constexpr float HYPOXIA_NECROTIC_O2  = 0.035f;
constexpr float HYPOXIA_NECROTIC_TIME = 120.0f;

// ── Warburg (Vander Heiden 2009 Science 324:1029) ────────────────────────
constexpr float MITO_HEALTH_SWITCH   = 0.4f;
constexpr float MITO_HEALTH_RESTORE  = 0.7f;
constexpr float WARBURG_GLY_BOOST    = 1.6f;
constexpr float WARBURG_OXP_PENALTY  = 0.5f;

// ── Contact inhibition (Delarue 2018 Dev Cell) ──────────────────────────
constexpr float MECH_P21_COUPLING    = 0.018f;
constexpr float CONTACT_INHIBIT_NBRS = 6.0f;

// ── Cell cycle TARGETS (emergent, not enforced) ─────────────────────────
// Phase durations are no longer enforced by a wall-clock timer. They
// EMERGE from biomass synthesis kinetics, CDK ODE, replication-fork
// progression, and mitosis sub-stage advancement. Values below are the
// HeLa literature targets used to calibrate the underlying rate constants.
//
//   G1 = 11 bio-h   46 % of 24-h cycle  (BNID 108483; Casciari 1992)
//   S  =  8 bio-h   33 %                (BNID 109393; Chao 2019 eLife)
//   G2 =  4 bio-h   17 %                (BNID 109225)
//   M  =  1 bio-h    4 %                (BNID 109226; Held 2010)
constexpr float CYCLE_G1_TARGET_BIOSEC = 11.0f * 3600.0f;
constexpr float CYCLE_S_TARGET_BIOSEC  =  8.0f * 3600.0f;
constexpr float CYCLE_G2_TARGET_BIOSEC =  4.0f * 3600.0f;
constexpr float CYCLE_M_TARGET_BIOSEC  =  1.0f * 3600.0f;
constexpr float CYCLE_TOTAL_TARGET_BIOSEC = 24.0f * 3600.0f;
// Legacy "slow-time-unit" durations kept for the not-yet-migrated path
// in updateCellCycle (deleted once cycleTimer is fully retired).
constexpr float CYCLE_G1_DUR = 3.0f;
constexpr float CYCLE_S_DUR  = 2.2f;
constexpr float CYCLE_G2_DUR = 1.1f;
constexpr float CYCLE_M_DUR  = 0.2f;

// ── Per-bio-second biology rates ─────────────────────────────────────────
// Calibrated against literature where available. Each one is the dial
// the validator (scripts/validate_population.py) uses to retune emergent
// phase dwell times against the HeLa target window.
constexpr float BIOMASS_SYNTH_PER_BIOSEC          = 4.0e-5f;   // → +1.0 in ~7 h at full nutrients (tuned)
constexpr float BIOMASS_DECAY_PER_BIOSEC          = 4.0e-7f;   // ~0.001 / hour stress decay
constexpr float REPLICATION_FORK_BP_PER_BIOSEC    = 17.0f;     // 1 kb/min — Conti 2007 / scEdU-seq 2024
constexpr float TRANSLATION_AA_PER_BIOSEC         = 6.0f;      // per ribosome — Ingolia 2011
constexpr float MITO_FUSION_PER_BIOSEC            = 5.0f / 3600.0f;  // Liu 2009 — 5 cycles/h
constexpr float MITO_FISSION_PER_BIOSEC           = 5.0f / 3600.0f;
constexpr float MITO_BIOGENESIS_PER_BIOSEC        = 0.7f / (24.0f * 3600.0f); // ~half/day
constexpr float MITO_PHAGY_PER_BIOSEC             = 1.0f / (24.0f * 3600.0f); // 1/day basal

// ── DNA / Telomere (Blackburn 2001 Cell 106:661) ─────────────────────────
// TELO_INIT_LENGTH: middle of 5-15 kb human range (Harley 1990).
// TELO_LOSS_PER_DIV: 50-200 bp per division in somatic cells (Harley 1990 Nature).
//   Previously 20 — too low; updated to physiological.
// TELO_CRITICAL: Hayflick limit ~2-4 kb → senescence.
// DNA_BASE_ERROR_RATE: post-proofreading + MMR = 10⁻⁸ – 10⁻¹⁰ per bp
//   (Kunkel 2004). Pre-repair Pol δ rate is ~10⁻⁴.
// GENOME_MUTATION_RATE: dimensionless mutation probability applied to
//   heritable traits (glycolysisBias etc.) at division.
constexpr float TELO_INIT_LENGTH    = 10000.0f; // bp
constexpr float TELO_LOSS_PER_DIV   = 100.0f;   // bp per division (was 20 — Harley 1990: 50-200)
constexpr float TELO_CRITICAL       = 3000.0f;  // bp
constexpr float DNA_BASE_ERROR_RATE = 1.0e-8f;  // per bp (was 2.5e-4 — Kunkel 2004 post-repair)
constexpr float GENOME_MUTATION_RATE = 0.04f;   // dimensionless, for trait drift

// ── Physics ──────────────────────────────────────────────────────────────
// Cell radius anchored to HeLa: 10 µm median radius ÷ MICROMETERS_PER_WORLD_UNIT
// (5 µm/wu) = 2.0 wu. Variation is ±1 µm = ±0.2 wu (Harper & Brooks 2005 flow
// cytometry; data/reference/hela_reference.csv range 15-30 µm diameter).
constexpr float CELL_RADIUS_BASE = 2.00f;   // was 2.15 — now exactly 10 µm
constexpr float CELL_RADIUS_VAR  = 0.20f;   // was 0.35 — ±1 µm biological variation
constexpr float HERTZ_STIFFNESS  = 2.0f;
constexpr float CELL_DAMPING     = 0.92f;
constexpr float MOTILITY_SPEED   = 0.3f;

// ── Carrying capacity ────────────────────────────────────────────────────
constexpr int   CARRY_CAP = 900;

// ── Apoptosis (Green & Kroemer 2004 Science 305:626) ─────────────────────
constexpr float ATP_DANGER_THRESHOLD = 6.0f;
constexpr float ATP_DANGER_DURATION  = 90.0f;
constexpr float NECROTIC_ROS_RELEASE = 0.06f;

// ── Natural turnover at confluence ───────────────────────────────────────
// Real tissue: ~1-3% daily turnover even at homeostasis
constexpr float TURNOVER_AGE_THRESHOLD = 400.0f;  // medium-dt age units
constexpr float TURNOVER_PROB_PER_DT   = 0.0001f; // probability per medium-dt step

// ── Drug mechanism of action codes ───────────────────────────────────────
constexpr int MOA_ANTI_PROLIF   = 0;  // freezes cell cycle (CDK inhibition)
constexpr int MOA_PRO_APOPTOSIS = 1;  // triggers programmed cell death
constexpr int MOA_DNA_DAMAGE    = 2;  // crosslinks/intercalation → p53 → apoptosis
constexpr int MOA_MITO_TOXIN    = 3;  // collapses mitochondrial membrane potential
