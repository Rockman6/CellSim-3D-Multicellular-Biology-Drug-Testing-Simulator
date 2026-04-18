#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  StemCells.h — Stem cell biology and tissue homeostasis
//  Based on: Alberts Ch 22 "Stem Cells in Tissue Homeostasis and Regeneration"
//
//  Implements:
//    - Stem cell properties (self-renewal + differentiation)
//    - Stem cell niche (Wnt, Notch, BMP signals from support cells)
//    - Asymmetric vs symmetric division choice
//    - Transit-amplifying progenitor cells
//    - Tissue homeostasis (birth-death balance)
//    - Lgr5+ intestinal stem cells (crypt model)
//    - Hematopoietic stem cell hierarchy
//    - Satellite cells (quiescent muscle stem cells)
//    - Neural stem cells (hippocampus, SVZ)
//    - iPSC reprogramming (Oct4, Sox2, Klf4, Myc → pluripotency)
//    - Stem cell aging (telomere shortening, epigenetic drift)
//    - Regeneration (blastema formation in salamanders)
//
//  Units: relative levels (0-1), cell counts
// ══════════════════════════════════════════════════════════════════════════

enum StemCellType {
    SC_TOTIPOTENT = 0,     // Can make ALL cell types + placenta (zygote)
    SC_PLURIPOTENT,        // Can make all body cell types (ES/iPS)
    SC_MULTIPOTENT,        // Multiple lineages (HSC, neural SC)
    SC_OLIGOPOTENT,        // Few lineages (myeloid progenitor)
    SC_UNIPOTENT,          // One lineage (satellite cell → muscle)
    SC_PROGENITOR,         // Transit-amplifying (limited divisions)
    SC_DIFFERENTIATED,     // Terminally differentiated
};

enum TissueType {
    TISSUE_INTESTINE = 0,  // Lgr5+ crypt stem cells
    TISSUE_SKIN,           // Basal layer stem cells
    TISSUE_BLOOD,          // HSCs in bone marrow
    TISSUE_MUSCLE,         // Satellite cells
    TISSUE_NEURAL,         // Neural stem cells (SVZ, hippocampus)
    TISSUE_LIVER,          // Hepatocyte self-renewal
    TISSUE_GENERIC,
};

// ── Stem Cell Niche ─────────────────────────────────────────────────────
struct StemCellNiche {
    TissueType tissue = TISSUE_GENERIC;

    // Niche signals (from support cells)
    float Wnt_signal = 0.5f;       // Promotes self-renewal (intestine, hair follicle)
    float Notch_signal = 0.5f;     // Maintains stemness (intestine: Paneth cells)
    float BMP_signal = 0.0f;       // Promotes differentiation (opposes Wnt)
    float EGF_signal = 0.3f;       // Promotes proliferation
    float Noggin = 0.3f;           // BMP antagonist (maintains stemness)
    float R_spondin = 0.5f;        // Lgr5 ligand, amplifies Wnt

    // Physical niche properties
    float niche_capacity = 15.0f;  // Max stem cells the niche supports
    float stiffness = 1.0f;        // ECM stiffness (affects differentiation)
    float oxygen = 0.5f;           // O2 level (HSC niche is hypoxic)

    // Net stemness signal
    float stemness() const {
        float pro_stem = Wnt_signal * (1.0f + R_spondin) + Notch_signal * 0.5f + Noggin * 0.3f;
        float pro_diff = BMP_signal * 0.5f;
        return pro_stem / (pro_stem + pro_diff + 0.1f);
    }
};

// ── Stem Cell State ─────────────────────────────────────────────────────
struct StemCellState {
    StemCellType type = SC_MULTIPOTENT;

    // Self-renewal capacity
    float self_renewal_prob = 0.5f;    // Probability of symmetric self-renewal
    float remaining_divisions = 50.0f;  // Hayflick limit (telomere-dependent)

    // Pluripotency factors
    float Oct4 = 0.0f;                 // Pluripotency TF
    float Sox2 = 0.0f;                 // Pluripotency TF
    float Nanog = 0.0f;                // Pluripotency TF
    float Klf4 = 0.0f;                 // Reprogramming factor
    float Myc = 0.0f;                  // Reprogramming factor (also oncogene!)

    // Tissue-specific markers
    float Lgr5 = 0.0f;                 // Intestinal stem cell marker
    float CD34 = 0.0f;                 // HSC marker
    float Pax7 = 0.0f;                 // Satellite cell marker
    float Nestin = 0.0f;               // Neural stem cell marker

    // Quiescence
    bool quiescent = false;             // G0 state
    float quiescence_depth = 0.0f;      // How deeply quiescent (0-1)

    // Aging
    float epigenetic_age = 0.0f;        // Epigenetic clock (Horvath)
    float accumulated_damage = 0.0f;     // DNA damage burden

    // Division history
    int total_divisions = 0;
    int symmetric_self_renewals = 0;
    int asymmetric_divisions = 0;
    int differentiation_events = 0;
};

// ── Tissue Homeostasis ──────────────────────────────────────────────────
struct TissueHomeostasis {
    float stem_cell_count = 15.0f;      // Active stem cells
    float progenitor_count = 50.0f;     // Transit-amplifying cells
    float differentiated_count = 500.0f; // Mature cells
    float apoptotic_rate = 0.01f;       // Cell death rate (1/s)
    float birth_rate = 0.01f;           // New cell production rate
    float turnover_time = 86400.0f * 5; // ~5 days for intestine

    float homeostatic_balance() const {
        return birth_rate / (apoptotic_rate + 0.001f);
    }
};

// ══════════════════════════════════════════════════════════════════════════

class StemCellEngine {
public:
    StemCellState cell;
    StemCellNiche niche;
    TissueHomeostasis homeostasis;

    void init(TissueType tissue) {
        cell = StemCellState();
        niche = StemCellNiche();
        niche.tissue = tissue;
        homeostasis = TissueHomeostasis();

        switch (tissue) {
        case TISSUE_INTESTINE:
            cell.Lgr5 = 1.0f;
            niche.Wnt_signal = 0.8f; niche.Notch_signal = 0.6f;
            niche.R_spondin = 0.7f;
            homeostasis.turnover_time = 86400.0f * 5; // 5 days
            break;
        case TISSUE_BLOOD:
            cell.CD34 = 1.0f;
            cell.quiescent = true; // HSCs are mostly quiescent
            cell.quiescence_depth = 0.8f;
            niche.oxygen = 0.05f; // Hypoxic niche
            homeostasis.turnover_time = 86400.0f * 120; // ~4 months (RBC lifespan)
            break;
        case TISSUE_MUSCLE:
            cell.Pax7 = 1.0f;
            cell.quiescent = true;
            cell.quiescence_depth = 0.95f; // Very deeply quiescent
            cell.type = SC_UNIPOTENT;
            break;
        case TISSUE_NEURAL:
            cell.Nestin = 1.0f;
            niche.Notch_signal = 0.7f;
            homeostasis.turnover_time = 86400.0f * 30; // ~1 month (hippocampus)
            break;
        default:
            break;
        }
    }

    // ── Division decision ───────────────────────────────────────────
    enum DivisionOutcome { DIV_SYMMETRIC_RENEWAL, DIV_ASYMMETRIC, DIV_SYMMETRIC_DIFF, DIV_QUIESCENT };

    DivisionOutcome decideDivision() {
        if (cell.quiescent && cell.quiescence_depth > 0.5f) return DIV_QUIESCENT;
        if (cell.remaining_divisions <= 0) return DIV_QUIESCENT; // Senescent

        float stemness = niche.stemness();

        // Stochastic fate choice (biased by niche signals)
        // Ref: Alberts Ch 22 "Symmetric Divisions and Stochastic Fate Choice"
        float p_self_renew = stemness * 0.6f;
        float p_asymmetric = 0.3f;
        float p_diff = 1.0f - p_self_renew - p_asymmetric;

        float r = (float)rand() / RAND_MAX;
        if (r < p_self_renew) return DIV_SYMMETRIC_RENEWAL;
        if (r < p_self_renew + p_asymmetric) return DIV_ASYMMETRIC;
        return DIV_SYMMETRIC_DIFF;
    }

    void step(float dt) {
        dt = fminf(dt, 1.0f);

        // ── Niche signal integration ────────────────────────────────
        float stemness = niche.stemness();

        // Stem cell markers respond to niche
        if (niche.tissue == TISSUE_INTESTINE) {
            cell.Lgr5 = stemness;
        }

        // ── Quiescence entry/exit ───────────────────────────────────
        if (homeostasis.stem_cell_count > niche.niche_capacity) {
            cell.quiescence_depth += 0.01f * dt; // Too many → go quiescent
        }
        if (homeostasis.differentiated_count < homeostasis.differentiated_count * 0.8f) {
            cell.quiescence_depth -= 0.05f * dt; // Tissue loss → activate
        }
        cell.quiescence_depth = std::clamp(cell.quiescence_depth, 0.0f, 1.0f);
        cell.quiescent = cell.quiescence_depth > 0.5f;

        // ── Tissue homeostasis ──────────────────────────────────────
        // Differentiated cells die (apoptosis, shedding)
        float death = homeostasis.apoptotic_rate * homeostasis.differentiated_count * dt;
        homeostasis.differentiated_count -= death;

        // Progenitors differentiate
        float prog_diff = 0.1f * homeostasis.progenitor_count * dt;
        homeostasis.progenitor_count -= prog_diff;
        homeostasis.differentiated_count += prog_diff;

        // Stem cells produce progenitors
        float sc_division = 0.005f * homeostasis.stem_cell_count * (1.0f - cell.quiescence_depth) * dt;
        homeostasis.progenitor_count += sc_division * 2.0f; // Transit amplification

        // Progenitors amplify
        homeostasis.progenitor_count += 0.05f * homeostasis.progenitor_count * dt;

        homeostasis.birth_rate = sc_division + prog_diff;

        // Clamp
        homeostasis.stem_cell_count = fmaxf(1, homeostasis.stem_cell_count);
        homeostasis.progenitor_count = fmaxf(0, homeostasis.progenitor_count);
        homeostasis.differentiated_count = fmaxf(0, homeostasis.differentiated_count);

        // ── Aging ───────────────────────────────────────────────────
        cell.epigenetic_age += dt / 86400.0f / 365.0f; // Years
        cell.accumulated_damage += 0.0001f * dt;
    }

    // ── iPSC Reprogramming ──────────────────────────────────────────
    // Ref: Alberts Ch 22 "Induced Pluripotent Stem (iPS) Cells"
    // Yamanaka factors: Oct4, Sox2, Klf4, Myc
    void applyReprogrammingFactors(float oct4, float sox2, float klf4, float myc) {
        cell.Oct4 = oct4; cell.Sox2 = sox2; cell.Klf4 = klf4; cell.Myc = myc;
    }

    float reprogrammingProgress() const {
        // All four factors needed for reprogramming
        float factor_product = cell.Oct4 * cell.Sox2 * cell.Klf4 * cell.Myc;
        // Very inefficient (~0.1% success rate in real experiments)
        return fminf(1.0f, factor_product * 0.01f);
    }

    bool isPluripotent() const {
        return cell.Oct4 > 0.5f && cell.Sox2 > 0.5f && cell.Nanog > 0.3f;
    }
};
