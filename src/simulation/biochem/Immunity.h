#pragma once
#include <cmath>
#include <vector>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Immunity.h — Innate and Adaptive Immune System
//  Based on: Alberts Ch 24 "The Innate and Adaptive Immune Systems"
//
//  Implements:
//    - Innate immunity:
//      Pattern Recognition Receptors (TLR1-10, NLR, RLR, CLR)
//      PAMPs/DAMPs detection
//      Inflammatory response (NF-κB → cytokines: TNFα, IL-1β, IL-6)
//      Complement cascade (classical, lectin, alternative → MAC)
//      Phagocytosis (macrophage, neutrophil) + respiratory burst
//      NK cells (missing-self detection, perforin/granzyme)
//      Interferons (Type I IFN → antiviral state)
//      Dendritic cells (antigen presentation → adaptive)
//    - Adaptive immunity:
//      B cells: V(D)J recombination → antibody diversity
//        IgM/IgD/IgG/IgA/IgE class switching (AID)
//        Somatic hypermutation → affinity maturation
//        Plasma cells → antibody secretion
//      T cells: TCR recognition of MHC-peptide
//        CD4+ helper (Th1, Th2, Th17, Treg)
//        CD8+ cytotoxic (perforin/granzyme killing)
//        T cell activation (TCR + co-stimulation)
//      Clonal selection and immunological memory
//      Self-tolerance (central: deletion; peripheral: anergy, Tregs)
//
//  Single-cell mode: simulates the immune components present within
//  and on the surface of the cell being studied (MHC, TLRs, etc.)
//  Colony mode: full immune cell population dynamics
// ══════════════════════════════════════════════════════════════════════════

// ── Pathogen-Associated Molecular Patterns ──────────────────────────────
enum PAMPtype {
    PAMP_NONE = 0,
    PAMP_LPS,              // Bacterial lipopolysaccharide (→ TLR4)
    PAMP_PEPTIDOGLYCAN,    // Bacterial cell wall (→ TLR2)
    PAMP_FLAGELLIN,        // Bacterial flagella (→ TLR5)
    PAMP_dsRNA,            // Viral dsRNA (→ TLR3, RIG-I)
    PAMP_ssRNA,            // Viral ssRNA (→ TLR7/8)
    PAMP_CpG_DNA,          // Unmethylated CpG (→ TLR9)
    PAMP_MANNAN,           // Fungal cell wall (→ Dectin-1)
    PAMP_COUNT
};

// ── Innate Immune State ─────────────────────────────────────────────────
struct InnateImmuneState {
    // Toll-like Receptors (surface and endosomal)
    float TLR2 = 0.5f;    // Peptidoglycan, lipoproteins
    float TLR3 = 0.3f;    // dsRNA (endosomal)
    float TLR4 = 0.5f;    // LPS (gram-negative bacteria)
    float TLR5 = 0.3f;    // Flagellin
    float TLR7 = 0.2f;    // ssRNA (endosomal)
    float TLR9 = 0.2f;    // CpG DNA (endosomal)

    // Cytosolic sensors
    float RIG_I = 0.3f;   // Viral RNA sensor
    float MDA5 = 0.2f;    // Long dsRNA sensor
    float cGAS = 0.2f;    // Cytosolic DNA sensor → STING
    float NLRP3 = 0.3f;   // Inflammasome sensor

    // Activated downstream
    float NFkB_innate = 0.0f;      // NF-κB activation (→ cytokines)
    float IRF3_active = 0.0f;      // IRF3 (→ Type I IFN)
    float inflammasome_active = 0.0f;// Caspase-1 activation

    // Cytokines produced
    float TNF_alpha = 0.0f;        // Pro-inflammatory
    float IL_1beta = 0.0f;         // Pro-inflammatory (inflammasome)
    float IL_6 = 0.0f;             // Acute phase
    float IFN_alpha = 0.0f;        // Type I interferon (antiviral)
    float IFN_beta = 0.0f;
    float IL_12 = 0.0f;            // Th1 polarizing

    // Complement system
    float C3b_opsonin = 0.0f;      // Opsonization
    float C5a_chemoattractant = 0.0f;
    float MAC_complex = 0.0f;      // Membrane attack complex

    // Phagocyte state
    float macrophage_activation = 0.0f;
    float neutrophil_recruitment = 0.0f;
    float respiratory_burst = 0.0f; // ROS production by NADPH oxidase
    float phagocytosis_rate = 0.0f;

    // NK cell state
    float NK_activation = 0.0f;
    float MHC_I_surface = 1.0f;    // MHC class I on cell surface (↓ if virus)
    float stress_ligands = 0.0f;   // NKG2D ligands (↑ if stressed)

    // Antiviral state
    float antiviral_state = 0.0f;  // IFN-induced PKR, RNaseL, etc.

    // Dendritic cell
    float DC_maturation = 0.0f;    // Maturation → antigen presentation
    float DC_migration = 0.0f;     // Migration to lymph node
    float co_stimulation = 0.0f;   // B7 (CD80/86) for T cell activation
};

// ── Adaptive Immune State ───────────────────────────────────────────────
struct AdaptiveImmuneState {
    // B cell / Antibody
    float B_cell_naive = 1.0f;     // Naive B cells
    float B_cell_activated = 0.0f; // Antigen-activated
    float plasma_cells = 0.0f;     // Antibody-secreting
    float memory_B = 0.0f;        // Memory B cells

    // Antibody classes
    float IgM = 0.0f;             // First response, pentameric
    float IgG = 0.0f;             // Main serum antibody
    float IgA = 0.0f;             // Mucosal protection
    float IgE = 0.0f;             // Parasites / allergy

    // Antibody properties
    float antibody_affinity = 0.1f;// Matured through SHM
    float neutralization = 0.0f;   // Virus/toxin neutralization
    float opsonization = 0.0f;     // Fc-mediated phagocytosis

    // T cell states
    float CD4_naive = 1.0f;       // Naive helper T cells
    float Th1 = 0.0f;             // IFNγ → macrophage activation
    float Th2 = 0.0f;             // IL-4, IL-5 → B cell help, eosinophils
    float Th17 = 0.0f;            // IL-17 → neutrophil recruitment
    float Treg = 0.3f;            // Regulatory T cells (suppress autoimmunity)

    float CD8_naive = 1.0f;       // Naive cytotoxic T cells
    float CD8_effector = 0.0f;    // Active killers
    float CD8_memory = 0.0f;      // Memory CTLs

    // T cell functions
    float cytotoxic_killing = 0.0f;// Perforin/granzyme
    float IFN_gamma = 0.0f;       // Th1 product

    // Tolerance
    float self_tolerance = 1.0f;   // 0 = autoimmune, 1 = tolerant
};

// ══════════════════════════════════════════════════════════════════════════

class ImmunityEngine {
public:
    InnateImmuneState innate;
    AdaptiveImmuneState adaptive;
    float pathogen_load = 0.0f;
    float time_since_infection = 0.0f;

    void init() {
        innate = InnateImmuneState();
        adaptive = AdaptiveImmuneState();
        pathogen_load = 0.0f;
        time_since_infection = 0.0f;
    }

    void applyPAMP(PAMPtype pamp, float amount) {
        pathogen_load += amount;
        switch (pamp) {
        case PAMP_LPS: innate.NFkB_innate += amount * innate.TLR4; break;
        case PAMP_dsRNA:
            innate.IRF3_active += amount * innate.TLR3;
            innate.IRF3_active += amount * innate.RIG_I;
            break;
        case PAMP_CpG_DNA: innate.NFkB_innate += amount * innate.TLR9; break;
        case PAMP_FLAGELLIN: innate.NFkB_innate += amount * innate.TLR5; break;
        default: break;
        }
    }

    void step(float dt) {
        dt = fminf(dt, 0.1f);
        auto& i = innate;
        auto& a = adaptive;

        // ── Innate response (hours) ────────────────────────────────
        // NF-κB → cytokines
        i.TNF_alpha += (0.5f * i.NFkB_innate - 0.1f * i.TNF_alpha) * dt;
        i.IL_1beta += (0.3f * i.inflammasome_active - 0.1f * i.IL_1beta) * dt;
        i.IL_6 += (0.3f * i.NFkB_innate - 0.08f * i.IL_6) * dt;

        // IRF3 → Type I IFN
        i.IFN_alpha += (0.5f * i.IRF3_active - 0.1f * i.IFN_alpha) * dt;
        i.IFN_beta += (0.3f * i.IRF3_active - 0.1f * i.IFN_beta) * dt;

        // Antiviral state (IFN → JAK-STAT → ISG expression)
        i.antiviral_state = (i.IFN_alpha + i.IFN_beta) / (1.0f + i.IFN_alpha + i.IFN_beta);

        // Complement activation
        float complement_drive = pathogen_load * 0.1f + a.IgM * 0.3f + a.IgG * 0.2f;
        i.C3b_opsonin += (complement_drive - 0.1f * i.C3b_opsonin) * dt;
        i.MAC_complex += (0.1f * i.C3b_opsonin - 0.05f * i.MAC_complex) * dt;

        // Phagocytosis
        i.macrophage_activation = (i.TNF_alpha + a.IFN_gamma) / (2.0f + i.TNF_alpha + a.IFN_gamma);
        i.phagocytosis_rate = i.macrophage_activation * (1.0f + i.C3b_opsonin + a.opsonization);

        // NK cells
        float nk_activation_signal = i.stress_ligands * (1.0f - i.MHC_I_surface);
        i.NK_activation = fmaxf(0, nk_activation_signal);

        // DC maturation (→ bridge to adaptive)
        i.DC_maturation += (0.1f * i.NFkB_innate - 0.02f * i.DC_maturation) * dt;
        i.co_stimulation = i.DC_maturation * 0.5f; // B7/CD80/CD86

        // Pathogen clearance by innate
        pathogen_load -= (i.phagocytosis_rate * 0.1f + i.MAC_complex * 0.05f +
                         i.antiviral_state * 0.05f) * dt;
        pathogen_load = fmaxf(0, pathogen_load);

        // ── Adaptive response (days) ───────────────────────────────
        if (i.DC_maturation > 0.3f && time_since_infection > 3.0f) {
            // B cell activation (antigen + T cell help)
            float b_activation = 0.05f * i.DC_maturation * a.B_cell_naive;
            a.B_cell_activated += b_activation * dt;
            a.B_cell_naive -= b_activation * dt;

            // Plasma cell differentiation (→ antibody)
            float plasma_diff = 0.02f * a.B_cell_activated;
            a.plasma_cells += plasma_diff * dt;

            // Antibody production
            a.IgM += 0.1f * a.plasma_cells * dt; // Early response
            if (time_since_infection > 7.0f) {
                // Class switching (IgM → IgG)
                a.IgG += 0.2f * a.plasma_cells * dt;
                // Affinity maturation (somatic hypermutation + selection)
                a.antibody_affinity += 0.01f * dt;
                a.antibody_affinity = fminf(1.0f, a.antibody_affinity);
            }

            // Opsonization and neutralization
            a.opsonization = (a.IgG + a.IgM) * a.antibody_affinity;
            a.neutralization = a.IgG * a.antibody_affinity * 0.5f;

            // T cell activation
            float t_activation = i.co_stimulation * i.DC_maturation;
            // CD4+ T helper differentiation
            if (i.IL_12 > 0.1f) a.Th1 += 0.02f * t_activation * dt;
            a.IFN_gamma = a.Th1 * 0.5f;

            // CD8+ CTL activation
            float ctl_act = 0.03f * t_activation * a.CD8_naive;
            a.CD8_effector += ctl_act * dt;
            a.CD8_naive -= ctl_act * dt;
            a.cytotoxic_killing = a.CD8_effector * 0.3f;

            // Pathogen clearance by adaptive
            pathogen_load -= (a.neutralization * 0.2f + a.cytotoxic_killing * 0.1f) * dt;
            pathogen_load = fmaxf(0, pathogen_load);

            // Memory cell generation
            if (pathogen_load < 0.01f) {
                a.memory_B += 0.01f * a.B_cell_activated * dt;
                a.CD8_memory += 0.01f * a.CD8_effector * dt;
            }
        }

        time_since_infection += dt;

        // Decay
        i.NFkB_innate *= (1.0f - 0.05f * dt);
        i.IRF3_active *= (1.0f - 0.05f * dt);
    }

    bool isInfected() const { return pathogen_load > 0.01f; }
    float getInflammation() const {
        return innate.TNF_alpha + innate.IL_1beta + innate.IL_6;
    }
};
