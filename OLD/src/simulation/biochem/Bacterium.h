#pragma once
#include <simd/simd.h>
#include <string>
#include <vector>
#include <cstdint>

// ══════════════════════════════════════════════════════════════════════════
//  Bacterium.h — Agent-sized prokaryotic pathogen.
//
//  A Bacterium is a small cell in its own right. Unlike virions it can
//  replicate without hijacking a host (in the medium), adhere to or invade
//  host cells, and secrete chemicals into the fluid (toxins, QS molecules,
//  DNA via lysis / conjugation).
//
//  Structural fields drive BacteriumBuilder procedural rendering:
//  cocci/rod/vibrio/spiral outer shape, gram-pos thick wall OR gram-neg
//  inner+outer membrane with LPS fringe, flagella count, pilus count.
//
//  Infection discipline: bacteria bind host receptors via the same
//  BindingMatcher shape+chemistry score as virions. Secretion of toxins
//  into MediumField reuses the same chemical-release pathway as apoptotic
//  cells. Nothing is hardcoded by strain name.
// ══════════════════════════════════════════════════════════════════════════

enum class BacteriumShape : uint8_t {
    COCCUS      = 0,   // sphere
    ROD         = 1,   // bacillus
    VIBRIO      = 2,   // curved rod
    SPIRILLUM   = 3,   // short spiral
    SPIROCHETE  = 4,   // long spiral
};

enum class GramType : uint8_t {
    GRAM_POSITIVE = 0, // thick peptidoglycan, no outer membrane
    GRAM_NEGATIVE = 1, // inner + outer membrane, periplasm, LPS fringe
    ACID_FAST     = 2, // mycolic acid (TB)
};

enum class BacteriumStage : uint8_t {
    FREE        = 0,   // swimming in medium
    ADHERING    = 1,   // attached via adhesin to host
    INVADED     = 2,   // inside host cell (cytosol or vacuole)
    REPLICATING = 3,   // dividing either extracellular or intracellular
    LYSED       = 4,   // burst, contents released
};

struct BacteriumSpec {
    std::string id;                  // "ecoli_K12", "listeria_EGDe"
    std::string displayName;
    BacteriumShape shape = BacteriumShape::ROD;
    GramType       gram  = GramType::GRAM_NEGATIVE;
    float length_um      = 2.0f;
    float width_um       = 0.5f;
    float peptidoglycan_nm = 30.0f;
    bool has_outer_membrane = true;
    bool has_lps_fringe     = true;
    bool has_capsule        = false;
    int  flagella_count     = 4;    // peritrichous by default
    int  pili_count         = 50;

    // Adhesin descriptor — same Lipinski-like 5-vector as virion spikes.
    float adhesin_logP = 1.5f;
    float adhesin_mw   = 35000.0f;
    int   adhesin_hbd  = 8;
    int   adhesin_hba  = 12;
    int   adhesin_aromatic = 2;
    std::vector<std::string> preferredReceptors;

    // Replication
    float doublingTime_bioSec = 1200.0f;   // 20 bio-min (fast gram-neg)
    float maxExtracellularDensity = 50.0f; // bacteria / voxel before quorum

    // Secretion (toxin release into MediumField when adhering or invaded).
    float toxinRate_per_s = 0.0f;          // µM per bio-s per bacterium
    float toxinCytotoxicity = 0.02f;       // damageLevel increment per µM·s

    // Lysis-released chemistry (fraction of biomass as AA, glucose, DNA).
    float biomass_bm = 0.01f;
};

struct Bacterium {
    int specIdx = -1;
    simd_float3 position = {0,0,0};
    simd_float3 velocity = {0,0,0};
    BacteriumStage stage = BacteriumStage::FREE;
    int hostCellIdx = -1;
    float stageTimer_s = 0.0f;
    float divisionClock_s = 0.0f;
    float biomass_bm = 0.01f;

    // Run-and-tumble (Berg 1972) — tumble every ~1 bio-s on average,
    // reorient random, run straight otherwise.
    float tumbleTimer_s = 0.0f;
    simd_float3 runDir = {1,0,0};

    uint32_t uid = 0;
};
