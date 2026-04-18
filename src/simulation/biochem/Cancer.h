#pragma once
#include <cmath>
#include <algorithm>
#include <vector>

// ══════════════════════════════════════════════════════════════════════════
//  Cancer.h — Cancer biology simulation (Hallmarks of Cancer)
//  Based on: Alberts Ch 20 "Cancer"
//
//  Models the hallmarks as emergent properties of molecular mutations:
//    1. Sustaining proliferative signaling (Ras, Myc, EGFR amplification)
//    2. Evading growth suppressors (Rb, p53, APC loss)
//    3. Resisting cell death (Bcl-2 overexpression, p53 loss)
//    4. Enabling replicative immortality (telomerase reactivation)
//    5. Inducing angiogenesis (VEGF overexpression)
//    6. Activating invasion/metastasis (EMT, E-cadherin loss)
//    7. Deregulating metabolism (Warburg effect)
//    8. Avoiding immune destruction
//
//  Each mutation modifies specific pathway parameters in the other engines
//  Clonal evolution: fitness-based selection of mutant subclones
// ══════════════════════════════════════════════════════════════════════════

enum MutationType {
    MUT_NONE = 0,
    // Oncogenes (gain of function)
    MUT_RAS_G12V,          // Constitutively active Ras
    MUT_BRAF_V600E,        // Constitutively active B-Raf
    MUT_MYC_AMPLIFICATION, // Myc overexpression
    MUT_EGFR_AMPLIFICATION,// EGFR overexpression
    MUT_HER2_AMPLIFICATION,// HER2 overexpression
    MUT_PIK3CA_H1047R,     // Constitutively active PI3K
    MUT_BCR_ABL,           // Philadelphia chromosome fusion kinase
    MUT_CYCLIN_D_AMPLIFICATION, // Cyclin D1 overexpression
    MUT_BCL2_OVEREXPRESSION,    // Anti-apoptotic
    MUT_TERT_REACTIVATION,      // Telomerase reactivation

    // Tumor suppressors (loss of function)
    MUT_TP53_LOSS,         // p53 loss → no checkpoint, no apoptosis
    MUT_RB1_LOSS,          // Rb loss → unrestricted E2F
    MUT_APC_LOSS,          // APC loss → constitutive Wnt/β-catenin
    MUT_PTEN_LOSS,         // PTEN loss → constitutive PI3K/Akt
    MUT_BRCA1_LOSS,        // BRCA1 loss → HR repair deficiency
    MUT_BRCA2_LOSS,        // BRCA2 loss → HR repair deficiency
    MUT_P16_LOSS,          // p16/INK4a loss → CDK4/6 unrestricted
    MUT_ECADHERIN_LOSS,    // E-cadherin loss → EMT, invasion

    // Caretaker gene mutations
    MUT_MSH2_LOSS,         // MMR deficiency → hypermutator
    MUT_MLH1_LOSS,         // MMR deficiency (Lynch syndrome)

    MUT_COUNT
};

struct Mutation {
    MutationType type;
    bool heterozygous;     // One allele or both? (tumor suppressors need both = two-hit)
    float penetrance;      // How strongly it affects phenotype (0-1)
    int generation;        // When it arose (cell division count)
};

struct CancerHallmarks {
    float proliferative_signaling = 0.0f;  // Hallmark 1
    float growth_suppression_evasion = 0.0f;// Hallmark 2
    float death_resistance = 0.0f;          // Hallmark 3
    float replicative_immortality = 0.0f;   // Hallmark 4
    float angiogenesis_induction = 0.0f;    // Hallmark 5
    float invasion_metastasis = 0.0f;       // Hallmark 6
    float metabolic_reprogramming = 0.0f;   // Hallmark 7 (Warburg)
    float immune_evasion = 0.0f;            // Hallmark 8

    float overallMalignancy() const {
        return (proliferative_signaling + growth_suppression_evasion +
                death_resistance + replicative_immortality +
                angiogenesis_induction + invasion_metastasis +
                metabolic_reprogramming + immune_evasion) / 8.0f;
    }
};

class CancerEngine {
public:
    std::vector<Mutation> mutations;
    CancerHallmarks hallmarks;
    float mutator_phenotype = 1.0f;  // Mutation rate multiplier (>1 if MMR deficient)
    float fitness = 1.0f;            // Selective advantage
    int division_count = 0;
    bool is_transformed = false;     // Fully transformed cancer cell?

    void init() {
        mutations.clear();
        hallmarks = CancerHallmarks();
        mutator_phenotype = 1.0f;
        fitness = 1.0f;
        division_count = 0;
        is_transformed = false;
    }

    void addMutation(MutationType type, bool hetero = false) {
        Mutation m;
        m.type = type;
        m.heterozygous = hetero;
        m.penetrance = hetero ? 0.5f : 1.0f;
        m.generation = division_count;
        mutations.push_back(m);
        updateHallmarks();
    }

    // Called each cell division — chance of spontaneous mutation
    void onDivision(float base_mutation_rate) {
        division_count++;
        float rate = base_mutation_rate * mutator_phenotype;

        // Random driver mutation (very rare)
        if (randF() < rate * 0.001f) {
            MutationType type = (MutationType)(1 + (rand() % (MUT_COUNT - 1)));
            addMutation(type, true); // Start heterozygous
        }

        // Loss of heterozygosity (second hit for tumor suppressors)
        for (auto& m : mutations) {
            if (m.heterozygous && randF() < rate * 0.01f) {
                m.heterozygous = false;
                m.penetrance = 1.0f;
            }
        }

        updateHallmarks();
    }

    // Apply mutation effects to other simulation engines
    struct MutationEffects {
        float ras_constitutive = 0.0f;     // → Signaling.Ras_GTP override
        float pi3k_constitutive = 0.0f;    // → Signaling.PI3K_active override
        float p53_functional = 1.0f;       // → reduced checkpoint
        float rb_functional = 1.0f;        // → unrestricted E2F
        float pten_functional = 1.0f;      // → PIP3 accumulation
        float bcl2_level = 1.0f;           // → death resistance
        float telomerase_active = false;   // → unlimited replication
        float ecadherin_level = 1.0f;      // → invasion
        float mmr_functional = 1.0f;       // → mutator phenotype
        float wnt_constitutive = 0.0f;     // → β-catenin accumulation
        float warburg_glycolysis = 0.0f;   // → aerobic glycolysis
    };

    MutationEffects getEffects() const {
        MutationEffects e;
        for (auto& m : mutations) {
            float p = m.penetrance;
            switch (m.type) {
            case MUT_RAS_G12V: e.ras_constitutive += 0.8f * p; break;
            case MUT_BRAF_V600E: e.ras_constitutive += 0.5f * p; break;
            case MUT_PIK3CA_H1047R: e.pi3k_constitutive += 0.7f * p; break;
            case MUT_TP53_LOSS: e.p53_functional -= 0.9f * p; break;
            case MUT_RB1_LOSS: e.rb_functional -= 0.9f * p; break;
            case MUT_PTEN_LOSS: e.pten_functional -= 0.9f * p; break;
            case MUT_BCL2_OVEREXPRESSION: e.bcl2_level += 3.0f * p; break;
            case MUT_TERT_REACTIVATION: e.telomerase_active = true; break;
            case MUT_ECADHERIN_LOSS: e.ecadherin_level -= 0.8f * p; break;
            case MUT_APC_LOSS: e.wnt_constitutive += 0.8f * p; break;
            case MUT_MSH2_LOSS: case MUT_MLH1_LOSS:
                e.mmr_functional -= 0.9f * p; break;
            case MUT_MYC_AMPLIFICATION: e.warburg_glycolysis += 0.3f * p; break;
            default: break;
            }
        }
        e.p53_functional = fmaxf(0, e.p53_functional);
        e.rb_functional = fmaxf(0, e.rb_functional);
        e.pten_functional = fmaxf(0, e.pten_functional);
        e.ecadherin_level = fmaxf(0, e.ecadherin_level);
        e.mmr_functional = fmaxf(0, e.mmr_functional);
        return e;
    }

private:
    void updateHallmarks() {
        auto e = getEffects();
        hallmarks.proliferative_signaling = fminf(1, e.ras_constitutive + e.pi3k_constitutive + e.wnt_constitutive);
        hallmarks.growth_suppression_evasion = 1.0f - fminf(e.p53_functional, e.rb_functional);
        hallmarks.death_resistance = fminf(1, (e.bcl2_level - 1.0f) * 0.3f + (1.0f - e.p53_functional) * 0.5f);
        hallmarks.replicative_immortality = e.telomerase_active ? 1.0f : 0.0f;
        hallmarks.invasion_metastasis = 1.0f - e.ecadherin_level;
        hallmarks.metabolic_reprogramming = fminf(1, e.warburg_glycolysis);
        mutator_phenotype = 1.0f / fmaxf(0.1f, e.mmr_functional);
        fitness = 1.0f + hallmarks.overallMalignancy();
        is_transformed = hallmarks.overallMalignancy() > 0.5f;
    }

    float randF() { return (float)rand() / RAND_MAX; }
};
