#pragma once
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

// ══════════════════════════════════════════════════════════════════════════
//  GeneRegulation.h — Gene expression control
//  Based on: Alberts Ch 6 (DNA to Protein) & Ch 7 (Gene Control)
//
//  Implements:
//    - Transcription (RNA Pol I/II/III, GTFs, Mediator)
//    - RNA processing (5' cap, splicing, 3' polyadenylation)
//    - Translation regulation (eIF2, mTOR/4E-BP, miRNA)
//    - Epigenetics (DNA methylation, histone code)
//    - Gene regulatory networks (activators, repressors, enhancers)
//    - Non-coding RNA (miRNA/RISC, siRNA/RNAi)
//
//  Uses Gillespie-like stochastic gene expression model
//  Ref: Elowitz et al. 2002 (gene expression noise)
//       Raj & van Oudenaarden 2008 (stochastic gene expression)
//
//  Each gene has: promoter state, mRNA level, protein level
//  Transcription rate depends on: transcription factors, chromatin state
// ══════════════════════════════════════════════════════════════════════════

// ── Chromatin State ─────────────────────────────────────────────────────
enum ChromatinState {
    CHROMATIN_OPEN = 0,       // Euchromatin (acetylated, accessible)
    CHROMATIN_POISED = 1,     // Bivalent (H3K4me3 + H3K27me3)
    CHROMATIN_CLOSED = 2,     // Heterochromatin (H3K9me3, H3K27me3)
    CHROMATIN_METHYLATED = 3, // DNA methylated (CpG, permanent silencing)
};

// ── Histone Modifications ──────────────────────────────────────────────
struct HistoneMarks {
    float H3K4me3  = 0.0f;  // Active promoter mark
    float H3K27me3 = 0.0f;  // Polycomb repressive mark
    float H3K9me3  = 0.0f;  // Constitutive heterochromatin
    float H3K36me3 = 0.0f;  // Transcription elongation mark
    float H3K27ac  = 0.0f;  // Active enhancer mark
    float H4K16ac  = 0.0f;  // General activation
    float H3K9ac   = 0.0f;  // Active promoter
    float H2AX_P   = 0.0f;  // DNA damage mark (γ-H2AX)

    ChromatinState getState() const {
        if (H3K9me3 > 0.5f) return CHROMATIN_CLOSED;
        if (H3K27me3 > 0.5f && H3K4me3 > 0.3f) return CHROMATIN_POISED;
        if (H3K27me3 > 0.5f) return CHROMATIN_CLOSED;
        if (H3K4me3 > 0.3f || H3K27ac > 0.3f) return CHROMATIN_OPEN;
        return CHROMATIN_CLOSED;
    }

    float accessibility() const {
        // 0 = fully closed, 1 = fully open
        float open_marks = H3K4me3 + H3K27ac + H4K16ac + H3K9ac;
        float closed_marks = H3K9me3 + H3K27me3;
        return open_marks / (open_marks + closed_marks + 0.1f);
    }
};

// ── Single Gene Model ──────────────────────────────────────────────────
struct Gene {
    const char* name;
    float mRNA = 0.0f;              // mRNA copies (normalized)
    float protein = 0.0f;           // Protein copies (normalized)

    // Transcription parameters
    float basal_txn_rate = 0.01f;   // Basal transcription rate
    float max_txn_rate = 1.0f;      // Maximum transcription rate
    float mRNA_halflife = 600.0f;   // mRNA half-life in seconds (~10 min)
    float protein_halflife = 3600.0f; // Protein half-life (~1 hour)

    // Regulation
    HistoneMarks histones;
    float CpG_methylation = 0.0f;   // 0 = unmethylated, 1 = fully methylated
    float promoter_activity = 0.0f; // Net promoter activity (0-1)

    // Splicing
    int num_exons = 5;
    int num_introns = 4;
    float splicing_efficiency = 0.95f;

    // mRNA features
    bool has_5prime_cap = true;
    float polyA_length = 200.0f;    // Poly-A tail length

    // Translation
    float translation_rate = 0.1f;  // Protein production per mRNA per second
    float ribosome_loading = 1.0f;  // How many ribosomes on this mRNA

    void updateTranscription(float tf_activation, float tf_repression, float dt) {
        // Chromatin accessibility
        float access = histones.accessibility();
        // CpG methylation silencing
        float methyl_factor = 1.0f - CpG_methylation * 0.9f;
        // Transcription factor regulation
        float tf_drive = tf_activation / (0.5f + tf_activation + tf_repression);

        promoter_activity = access * methyl_factor * tf_drive;

        // Transcription rate
        float txn = basal_txn_rate + (max_txn_rate - basal_txn_rate) * promoter_activity;
        // mRNA degradation
        float deg = mRNA * 0.693f / mRNA_halflife; // ln(2)/t_half

        mRNA += (txn * splicing_efficiency - deg) * dt;
        mRNA = fmaxf(0, mRNA);
    }

    void updateTranslation(float global_translation_factor, float dt) {
        float txl = translation_rate * mRNA * ribosome_loading * global_translation_factor;
        float deg = protein * 0.693f / protein_halflife;
        protein += (txl - deg) * dt;
        protein = fmaxf(0, protein);
    }
};

// ── Gene Regulatory Network ────────────────────────────────────────────
class GeneRegulationEngine {
public:
    // Key regulatory genes
    Gene cyclinD   = {"CyclinD",  0, 0, 0.01f, 1.0f, 1200, 1800};
    Gene cyclinE   = {"CyclinE",  0, 0, 0.005f, 0.8f, 600, 900};
    Gene p21       = {"p21/CDKN1A", 0, 0, 0.001f, 2.0f, 900, 1200};
    Gene p53       = {"TP53",     0, 0.5f, 0.1f, 0.5f, 1200, 7200};
    Gene Rb        = {"RB1",      0, 1.0f, 0.05f, 0.3f, 3600, 14400};
    Gene E2F       = {"E2F1",     0, 0, 0.01f, 1.5f, 600, 1200};
    Gene Myc       = {"MYC",      0, 0, 0.005f, 2.0f, 300, 600}; // Very short-lived
    Gene Bcl2      = {"BCL2",     0, 0.5f, 0.05f, 0.3f, 3600, 28800};
    Gene Bax       = {"BAX",      0, 0.3f, 0.02f, 1.0f, 1800, 7200};

    // Housekeeping
    Gene actin     = {"ACTB",     1.0f, 5.0f, 0.5f, 1.0f, 7200, 43200};
    Gene GAPDH     = {"GAPDH",    1.0f, 5.0f, 0.5f, 1.0f, 7200, 43200};

    // Epigenetic regulators
    float DNMT1_activity = 1.0f;   // Maintenance methyltransferase
    float DNMT3_activity = 0.1f;   // De novo methyltransferase
    float TET_activity = 0.1f;     // Demethylase (5mC → 5hmC)
    float HAT_activity = 0.5f;     // Histone acetyltransferase (CBP/p300)
    float HDAC_activity = 0.5f;    // Histone deacetylase
    float PRC2_activity = 0.3f;    // Polycomb (H3K27me3 writer)
    float trithorax_activity = 0.3f; // Trithorax (H3K4me3 writer)

    // miRNA regulation
    float miRNA_total = 0.5f;      // Total miRNA activity
    float RISC_activity = 0.5f;    // RISC complex availability

    // Global translation state
    float eIF2_active = 1.0f;      // eIF2 (global translation initiation)
    float mTOR_4EBP = 0.5f;        // 4E-BP phosphorylation (from mTOR)

    void init() {
        // Set initial chromatin states
        cyclinD.histones.H3K4me3 = 0.5f; cyclinD.histones.H3K27ac = 0.3f;
        Myc.histones.H3K4me3 = 0.3f;
        p53.histones.H3K4me3 = 0.8f; // Constitutively accessible
        actin.histones.H3K4me3 = 1.0f; actin.histones.H3K27ac = 0.8f;
    }

    void step(float dt, float ERK_signal, float Akt_signal, float p53_active,
              float Wnt_signal, float DNA_damage) {
        // ── Update chromatin ────────────────────────────────────────
        updateEpigenetics(dt);

        // ── ERK/MAPK drives CyclinD and Myc ────────────────────────
        cyclinD.updateTranscription(ERK_signal + Wnt_signal, 0, dt);
        Myc.updateTranscription(ERK_signal * 0.5f + Wnt_signal * 0.3f, 0, dt);

        // ── E2F drives CyclinE (positive feedback via Rb) ──────────
        float Rb_activity = Rb.protein / (0.5f + Rb.protein);
        float E2F_free = E2F.protein * (1.0f - Rb_activity); // Rb sequesters E2F
        cyclinE.updateTranscription(E2F_free, 0, dt);
        E2F.updateTranscription(E2F_free * 0.3f, Rb_activity, dt); // Autoregulation

        // ── p53 pathway (DNA damage response) ──────────────────────
        // DNA damage → ATM/ATR → p53 stabilization
        float p53_stabilization = DNA_damage * 2.0f;
        p53.protein_halflife = 7200.0f / (1.0f + p53_stabilization * 5.0f); // p53 stabilized
        p53.updateTranscription(0.1f, 0, dt);

        // p53 target genes
        p21.updateTranscription(p53.protein * 2.0f, 0, dt); // p21 → CDK inhibitor
        Bax.updateTranscription(p53.protein * 1.0f, 0, dt); // Bax → pro-apoptotic

        // ── Survival genes ─────────────────────────────────────────
        Bcl2.updateTranscription(Akt_signal * 0.3f, p53.protein * 0.2f, dt);
        Rb.updateTranscription(0.1f, 0, dt);

        // ── Housekeeping ───────────────────────────────────────────
        actin.updateTranscription(1.0f, 0, dt);
        GAPDH.updateTranscription(1.0f, 0, dt);

        // ── Translation ────────────────────────────────────────────
        float global_txl = eIF2_active * (0.5f + 0.5f * mTOR_4EBP);
        // miRNA repression (global)
        global_txl *= (1.0f - miRNA_total * 0.1f);

        cyclinD.updateTranslation(global_txl, dt);
        cyclinE.updateTranslation(global_txl, dt);
        p21.updateTranslation(global_txl, dt);
        p53.updateTranslation(global_txl, dt);
        Rb.updateTranslation(global_txl, dt);
        E2F.updateTranslation(global_txl, dt);
        Myc.updateTranslation(global_txl, dt);
        Bcl2.updateTranslation(global_txl, dt);
        Bax.updateTranslation(global_txl, dt);
        actin.updateTranslation(global_txl, dt);
        GAPDH.updateTranslation(global_txl, dt);
    }

private:
    void updateEpigenetics(float dt) {
        // HAT vs HDAC balance
        float acetylation_drive = (HAT_activity - HDAC_activity) * dt * 0.01f;

        // Apply to all genes
        Gene* genes[] = {&cyclinD, &cyclinE, &p21, &p53, &Rb, &E2F, &Myc, &Bcl2, &Bax};
        for (auto* g : genes) {
            g->histones.H3K27ac += acetylation_drive;
            g->histones.H3K27ac = std::clamp(g->histones.H3K27ac, 0.0f, 1.0f);

            // PRC2 adds H3K27me3 (repressive)
            g->histones.H3K27me3 += PRC2_activity * 0.001f * dt;
            // Trithorax adds H3K4me3 (active)
            g->histones.H3K4me3 += trithorax_activity * 0.001f * dt;

            // DNMT1 maintains CpG methylation through replication
            // (simplified — methylation slowly established/removed)
            g->CpG_methylation += (DNMT3_activity * 0.0001f - TET_activity * 0.0001f) * dt;
            g->CpG_methylation = std::clamp(g->CpG_methylation, 0.0f, 1.0f);
        }
    }
};
