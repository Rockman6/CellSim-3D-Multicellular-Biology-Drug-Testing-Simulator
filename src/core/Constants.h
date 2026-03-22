#pragma once

// ══════════════════════════════════════════════════════════════════════════
//  Constants.h — Biology constants calibrated against published data
//  Target: HeLa-like cells, ~20h doubling, G1:46% S:33% G2:17% M:4%
//  References: Novak & Tyson (2008), Casciari (1992), Green & Kroemer (2004)
// ══════════════════════════════════════════════════════════════════════════

// ── Scene ────────────────────────────────────────────────────────────────
constexpr float FLOOR_Y      = -7.0f;
constexpr float SCENE_BOUND  = 55.0f;
constexpr int   MAX_CELLS    = 2000;  // hard cap for memory
// Plate carrying capacity: π×55²/(π×2.15²) × 0.9 hex packing ≈ 590 monolayer cells
// Growth should plateau around this number via contact inhibition
constexpr int   MONOLAYER_CAPACITY = 600;
constexpr int   INIT_CELLS   = 5;

// ── Timescales ───────────────────────────────────────────────────────────
constexpr float FAST_DT_SCALE   = 1.0f;
constexpr float MEDIUM_DT_SCALE = 0.22f;
constexpr float SLOW_DT_SCALE   = 0.06f;
constexpr float BIO_MIN_PER_SEC = 3.0f;
constexpr float CDK_DT_SCALE    = 0.8f;    // Novak-Tyson CDK dt scaling factor

// ── Lag phase ────────────────────────────────────────────────────────────
// Ref: AAT Bioquest — lag phase 1-6h in mammalian culture
constexpr float LAG_DURATION    = 0.3f;    // slow-time units ≈ 3 bio-hours

// ── Diffusion (Fick's 2nd law) ───────────────────────────────────────────
constexpr int   NUTRIENT_GRID_N   = 64;
constexpr float DIFF_O2_COEFF     = 1.8f;
constexpr float DIFF_GLC_COEFF    = 0.9f;
constexpr float DIFF_CO2_COEFF    = 1.5f;
constexpr float DIFF_PH_COEFF     = 0.6f;
// Consumption rates — scaled for meaningful depletion at scale
constexpr float O2_CONSUME_BASE   = 0.0075f;   // was 0.0015
constexpr float GLC_CONSUME_BASE  = 0.0040f;   // was 0.0008
constexpr float CO2_PRODUCE_BASE  = 0.0060f;   // was 0.0012
constexpr float LACTATE_PRODUCE   = 0.0030f;   // was 0.0006

// ── Metabolism ───────────────────────────────────────────────────────────
constexpr float BIOMASS_SYNTH_K       = 0.010f;
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

// ── Cell cycle durations ─────────────────────────────────────────────────
// HeLa: G1≈11h, S≈8h, G2≈4h, M≈1h (Campisi 2007)
constexpr float CYCLE_G1_DUR = 3.0f;
constexpr float CYCLE_S_DUR  = 2.2f;
constexpr float CYCLE_G2_DUR = 1.1f;
constexpr float CYCLE_M_DUR  = 0.2f;

// ── DNA / Telomere (Blackburn 2001 Cell 106:661) ─────────────────────────
constexpr float TELO_INIT_LENGTH    = 10000.0f;
constexpr float TELO_LOSS_PER_DIV   = 20.0f;
constexpr float TELO_CRITICAL       = 3000.0f;
constexpr float DNA_BASE_ERROR_RATE = 2.5e-4f;
constexpr float GENOME_MUTATION_RATE = 0.04f;

// ── Physics ──────────────────────────────────────────────────────────────
constexpr float CELL_RADIUS_BASE = 2.15f;
constexpr float CELL_RADIUS_VAR  = 0.35f;
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
