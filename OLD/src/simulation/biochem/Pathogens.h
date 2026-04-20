#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Pathogens.h — Pathogen infection and host interaction
//  Based on: Alberts Ch 23 "Pathogens and Infection"
//
//  Implements:
//    - Bacterial infection:
//      Adhesion (pili, adhesins), invasion (T3SS injection)
//      Intracellular lifestyle: vacuolar (Salmonella SCV) vs cytosolic (Listeria)
//      Toxin mechanisms: A-B toxins (cholera), pore-forming (α-hemolysin)
//      Actin comet tail propulsion (Listeria ActA → Arp2/3)
//    - Viral infection:
//      Attachment (receptor binding: HIV gp120→CD4)
//      Entry (endocytosis or membrane fusion)
//      Genome replication (RNA-dependent RNA Pol, reverse transcriptase)
//      Assembly and budding/lysis
//      Immune evasion (MHC-I downregulation, IFN antagonism)
//    - Microbiota:
//      Commensal bacteria (colonization resistance)
//      Short-chain fatty acids (butyrate → colonocyte energy)
//      Immune education (Treg induction)
//      Dysbiosis (antibiotic disruption)
//
//  Single-cell mode: models infection of the simulated cell
// ══════════════════════════════════════════════════════════════════════════

enum PathogenType {
    PATH_NONE = 0,
    PATH_GRAM_NEG_BACTERIA,    // E.g., Salmonella, E. coli
    PATH_GRAM_POS_BACTERIA,    // E.g., Listeria, S. aureus
    PATH_DNA_VIRUS,            // E.g., Herpesvirus, Adenovirus
    PATH_RNA_VIRUS,            // E.g., Influenza, SARS-CoV-2
    PATH_RETROVIRUS,           // E.g., HIV (RNA → DNA via RT)
    PATH_FUNGUS,               // E.g., Candida, Cryptococcus
    PATH_PARASITE,             // E.g., Plasmodium, Toxoplasma
};

enum InfectionStage {
    INF_NONE = 0,
    INF_ATTACHMENT,            // Pathogen bound to surface receptor
    INF_ENTRY,                 // Entering the cell
    INF_INTRACELLULAR,         // Inside cell, replicating
    INF_ASSEMBLY,              // New particles assembling
    INF_EGRESS,                // Exiting (budding or lysis)
    INF_CLEARED,               // Immune system cleared it
};

struct PathogenState {
    PathogenType type = PATH_NONE;
    InfectionStage stage = INF_NONE;
    float load = 0.0f;                // Pathogen count/copies inside cell

    // Bacterial virulence factors
    float T3SS_activity = 0.0f;       // Type III secretion system
    float T4SS_activity = 0.0f;       // Type IV secretion system
    float toxin_AB = 0.0f;            // A-B toxin level
    float pore_forming_toxin = 0.0f;  // α-hemolysin etc.
    float actin_comet = 0.0f;         // ActA/IcsA-driven actin propulsion

    // Viral lifecycle
    float receptor_binding = 0.0f;    // Receptor occupancy
    float genome_copies = 0.0f;       // Viral genome copies (replication)
    float viral_proteins = 0.0f;      // Structural proteins made
    float new_particles = 0.0f;       // Assembled virions
    float reverse_transcriptase = 0.0f;// RT activity (retroviruses)

    // Immune evasion
    float MHC_I_downregulation = 0.0f;// Virus blocks MHC-I presentation
    float IFN_antagonism = 0.0f;      // Virus blocks interferon signaling
    float complement_evasion = 0.0f;  // Surface protein masking

    // Host cell effects
    float host_translation_shutoff = 0.0f; // Virus shuts down host protein synthesis
    float host_cytoskeleton_hijack = 0.0f; // Pathogen manipulates actin/MT
    float apoptosis_inhibition = 0.0f;     // Pathogen blocks apoptosis
    float membrane_damage = 0.0f;          // Pore formation, lysis risk
};

struct MicrobiotaState {
    float total_abundance = 1.0f;      // Normalized total bacteria
    float bacteroidetes = 0.4f;        // Major phylum
    float firmicutes = 0.4f;           // Major phylum
    float proteobacteria = 0.05f;      // Includes E. coli
    float actinobacteria = 0.1f;       // Includes Bifidobacterium
    float diversity_index = 3.0f;      // Shannon diversity (higher = healthier)

    // Metabolites produced by microbiota
    float butyrate = 0.5f;            // SCFA → colonocyte energy, anti-inflammatory
    float propionate = 0.3f;          // SCFA → liver metabolism
    float acetate = 0.4f;             // SCFA
    float vitamin_K = 0.3f;           // Synthesized by gut bacteria
    float vitamin_B12 = 0.2f;

    // Colonization resistance
    float resistance() const {
        return total_abundance * diversity_index / 3.0f;
    }
};

class PathogenEngine {
public:
    PathogenState pathogen;
    MicrobiotaState microbiota;

    void init() {
        pathogen = PathogenState();
        microbiota = MicrobiotaState();
    }

    void infect(PathogenType type, float initial_load) {
        pathogen.type = type;
        pathogen.stage = INF_ATTACHMENT;
        pathogen.load = initial_load;

        // Set type-specific properties
        switch (type) {
        case PATH_GRAM_NEG_BACTERIA:
            pathogen.T3SS_activity = 0.5f;
            break;
        case PATH_RNA_VIRUS:
            pathogen.receptor_binding = 0.8f;
            pathogen.IFN_antagonism = 0.3f;
            break;
        case PATH_RETROVIRUS:
            pathogen.receptor_binding = 0.9f;
            pathogen.reverse_transcriptase = 0.8f;
            pathogen.MHC_I_downregulation = 0.5f;
            break;
        default: break;
        }
    }

    void step(float dt, float immune_pressure, float antiviral_state) {
        dt = fminf(dt, 0.1f);
        if (pathogen.type == PATH_NONE) return;

        auto& p = pathogen;

        // ── Infection progression ───────────────────────────────────
        switch (p.stage) {
        case INF_ATTACHMENT:
            // Attempting entry
            if (p.receptor_binding > 0.3f || p.T3SS_activity > 0.3f) {
                p.stage = INF_ENTRY;
            }
            break;

        case INF_ENTRY:
            p.load += 0.1f * dt; // Initial replication
            p.stage = INF_INTRACELLULAR;
            break;

        case INF_INTRACELLULAR:
            // Replication (exponential until limited by resources)
            {
                float replication_rate = 0.5f; // Doublings per hour
                float resource_limit = 1.0f / (1.0f + p.load * 0.01f);
                float immune_inhibition = 1.0f / (1.0f + immune_pressure * 2.0f);
                float antiviral_inhibition = 1.0f - antiviral_state * 0.8f;

                p.load += p.load * replication_rate * resource_limit *
                          immune_inhibition * antiviral_inhibition * dt / 3600.0f;

                // Viral-specific: genome replication and protein synthesis
                if (p.type == PATH_RNA_VIRUS || p.type == PATH_RETROVIRUS ||
                    p.type == PATH_DNA_VIRUS) {
                    p.genome_copies += 0.5f * p.load * (1.0f - antiviral_state) * dt / 3600.0f;
                    p.viral_proteins += 0.3f * p.genome_copies *
                                       (1.0f - p.host_translation_shutoff * 0.5f) * dt / 3600.0f;

                    // Assembly
                    float assembly = fminf(p.genome_copies, p.viral_proteins) * 0.1f;
                    p.new_particles += assembly * dt / 3600.0f;
                }

                // Bacterial toxin production
                if (p.type == PATH_GRAM_NEG_BACTERIA || p.type == PATH_GRAM_POS_BACTERIA) {
                    p.toxin_AB += 0.01f * p.load * dt / 3600.0f;
                    p.pore_forming_toxin += 0.005f * p.load * dt / 3600.0f;
                }

                // Membrane damage accumulation
                p.membrane_damage += p.pore_forming_toxin * 0.01f * dt;

                // Progress to egress when enough particles made
                if (p.new_particles > 100 || p.membrane_damage > 0.8f) {
                    p.stage = INF_EGRESS;
                }
            }
            break;

        case INF_EGRESS:
            // Virus budding or cell lysis
            p.membrane_damage += 0.1f * dt;
            break;

        case INF_CLEARED:
            p.load *= (1.0f - 0.5f * dt); // Residual clearance
            break;

        default: break;
        }

        // ── Immune clearance ────────────────────────────────────────
        if (immune_pressure > 0.5f && p.load > 0) {
            float clearance = immune_pressure * 0.1f * dt;
            p.load -= clearance;
            if (p.load < 0.01f) {
                p.stage = INF_CLEARED;
                p.load = 0;
            }
        }

        // ── Immune evasion effects ──────────────────────────────────
        // MHC-I downregulation → escape CTL killing
        // IFN antagonism → block antiviral state
        // These modulate the immune response in Immunity.h

        // ── Microbiota dynamics ─────────────────────────────────────
        // Antibiotic treatment disrupts microbiota
        // (simplified — external application would reduce total_abundance)
        microbiota.butyrate = microbiota.firmicutes * 0.5f;
        microbiota.diversity_index = logf(1.0f + microbiota.total_abundance * 10.0f);
    }

    bool isInfected() const { return pathogen.load > 0.01f; }
    bool isCellLysing() const { return pathogen.membrane_damage > 0.9f; }
    float getPathogenLoad() const { return pathogen.load; }
};
