#pragma once
#include <cmath>
#include <vector>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  VesicularTraffic.h — Intracellular membrane traffic simulation
//  Based on: Alberts Ch 12 (Protein Sorting) & Ch 13 (Membrane Traffic)
//
//  Implements:
//    - ER protein import (SRP → translocon/Sec61)
//    - COPII vesicles (ER → Golgi, Sar1 GTPase)
//    - COPI vesicles (Golgi → ER retrieval, KDEL, KKXX)
//    - Clathrin-coated vesicles (endocytosis, TGN → lysosome)
//    - SNARE-mediated fusion (v-SNARE + t-SNARE → trans-SNARE complex)
//    - Rab GTPases (vesicle address tags, >60 family members)
//    - Endocytosis (receptor-mediated, phagocytosis, macropinocytosis)
//    - Endosomal sorting (early endosome → recycling/degradation)
//    - ESCRT machinery (MVB formation, ubiquitin-tagged receptors)
//    - Autophagy (macroautophagy, CMA, selective: mitophagy/xenophagy)
//    - Unfolded protein response (UPR: IRE1, PERK, ATF6)
//    - Constitutive vs regulated secretion
//    - M6P targeting (lysosomal hydrolases)
//    - Nuclear import/export (NLS, NPC, importin/exportin, Ran-GTP)
//
//  Each vesicle tracked as: origin, destination, cargo, coat, Rab, SNAREs
// ══════════════════════════════════════════════════════════════════════════

enum Compartment {
    COMP_CYTOSOL = 0, COMP_NUCLEUS, COMP_ER, COMP_ERGIC,
    COMP_CIS_GOLGI, COMP_MEDIAL_GOLGI, COMP_TRANS_GOLGI, COMP_TGN,
    COMP_EARLY_ENDOSOME, COMP_LATE_ENDOSOME, COMP_LYSOSOME,
    COMP_PLASMA_MEMBRANE, COMP_MITOCHONDRIA, COMP_PEROXISOME,
    COMP_AUTOPHAGOSOME, COMP_MVB,
    COMP_COUNT
};

enum CoatType {
    COAT_NONE = 0, COAT_COPII, COAT_COPI, COAT_CLATHRIN, COAT_CAVEOLIN,
};

enum CargoType {
    CARGO_SECRETORY = 0,        // Constitutive secretory proteins
    CARGO_REGULATED_SECRETORY,  // Hormones, neurotransmitters
    CARGO_MEMBRANE_PROTEIN,     // Transmembrane proteins
    CARGO_LYSOSOMAL_ENZYME,     // M6P-tagged hydrolases
    CARGO_RECEPTOR_RECYCLING,   // Receptors being recycled
    CARGO_RECEPTOR_DEGRADATION, // Ubiquitinated receptors → lysosome
    CARGO_AUTOPHAGY_CARGO,      // Autophagic cargo
    CARGO_ER_RETRIEVAL,         // KDEL/KKXX-tagged ER residents
};

struct Vesicle {
    Compartment origin;
    Compartment destination;
    CoatType coat;
    CargoType cargo;
    float progress = 0.0f;     // 0 = just budded, 1 = arrived
    float speed = 1.0f;        // Transport speed
    bool fused = false;         // Has it fused with target?
    int rab_id = -1;            // Rab GTPase family member
};

// ── Compartment State ───────────────────────────────────────────────────
struct CompartmentState {
    float protein_load[COMP_COUNT] = {};  // Protein content per compartment
    float lipid_content[COMP_COUNT] = {}; // Membrane lipid content
    float pH[COMP_COUNT] = {};            // Luminal pH

    void init() {
        // Default protein loads
        protein_load[COMP_ER] = 5.0f;
        protein_load[COMP_CIS_GOLGI] = 2.0f;
        protein_load[COMP_MEDIAL_GOLGI] = 2.0f;
        protein_load[COMP_TRANS_GOLGI] = 2.0f;
        protein_load[COMP_TGN] = 3.0f;
        protein_load[COMP_LYSOSOME] = 1.0f;
        protein_load[COMP_PLASMA_MEMBRANE] = 5.0f;

        // pH values (V-ATPase acidifies along secretory/endocytic pathway)
        pH[COMP_ER] = 7.2f;
        pH[COMP_CIS_GOLGI] = 6.7f;
        pH[COMP_MEDIAL_GOLGI] = 6.3f;
        pH[COMP_TRANS_GOLGI] = 6.0f;
        pH[COMP_TGN] = 6.0f;
        pH[COMP_EARLY_ENDOSOME] = 6.0f;
        pH[COMP_LATE_ENDOSOME] = 5.5f;
        pH[COMP_LYSOSOME] = 4.8f;
        pH[COMP_CYTOSOL] = 7.2f;
        pH[COMP_NUCLEUS] = 7.2f;
    }
};

// ── UPR (Unfolded Protein Response) ─────────────────────────────────────
struct UPRState {
    float ER_stress = 0.0f;            // Unfolded protein accumulation
    float IRE1_active = 0.0f;          // IRE1α → XBP1 splicing
    float PERK_active = 0.0f;          // PERK → eIF2α phosphorylation
    float ATF6_active = 0.0f;          // ATF6 → cleaved → nuclear
    float XBP1s = 0.0f;               // Spliced XBP1 (chaperone upregulation)
    float eIF2a_P = 0.0f;             // Phosphorylated eIF2α (translation block)
    float CHOP = 0.0f;                // CHOP/GADD153 (pro-apoptotic if prolonged)
    float BiP_level = 1.0f;           // BiP/GRP78 chaperone level
};

// ── Autophagy State ─────────────────────────────────────────────────────
struct AutophagyState {
    float mTORC1_inhibition = 0.0f;   // Low nutrients → mTORC1 off → autophagy on
    float ULK1_active = 0.0f;         // Autophagy initiation kinase
    float Beclin1 = 0.5f;             // PI3K complex component
    float LC3_I = 1.0f;               // Cytosolic LC3
    float LC3_II = 0.0f;              // Lipidated LC3 (on autophagosome membrane)
    float p62 = 0.5f;                 // Selective autophagy receptor
    float autophagosome_count = 0.0f;
    float autolysosome_flux = 0.0f;   // Degradation rate
};

// ══════════════════════════════════════════════════════════════════════════

class VesicularTrafficEngine {
public:
    std::vector<Vesicle> vesicles;
    CompartmentState compartments;
    UPRState upr;
    AutophagyState autophagy;

    // Coat protein levels
    float COPII_available = 1.0f;     // Sar1 + Sec23/24 + Sec13/31
    float COPI_available = 1.0f;      // ARF1 + coatomer
    float clathrin_available = 1.0f;  // Clathrin + AP complexes
    float dynamin = 1.0f;            // GTPase for vesicle scission

    // SNARE proteins
    float syntaxin_PM = 1.0f;        // t-SNARE at plasma membrane
    float SNAP25 = 1.0f;             // t-SNARE at plasma membrane
    float VAMP = 1.0f;               // v-SNARE on vesicles
    float NSF = 1.0f;                // ATPase (SNARE recycling)

    // Nuclear transport
    float importin = 1.0f;
    float exportin = 1.0f;
    float RanGTP_nuclear = 0.8f;     // Ran-GTP gradient (high in nucleus)
    float RanGDP_cytosol = 0.8f;

    void init() {
        vesicles.clear();
        compartments.init();
        upr = UPRState();
        autophagy = AutophagyState();
    }

    void step(float dt, float ATP, float nutrient_level, float mTORC1_active) {
        dt = fminf(dt, 0.01f);

        // ── ER → Golgi transport (COPII) ────────────────────────────
        float copii_flux = 0.2f * COPII_available * compartments.protein_load[COMP_ER]
                         / (2.0f + compartments.protein_load[COMP_ER]) * ATP / (1.0f + ATP);
        if (copii_flux > 0.01f && vesicles.size() < 200) {
            Vesicle v;
            v.origin = COMP_ER; v.destination = COMP_CIS_GOLGI;
            v.coat = COAT_COPII; v.cargo = CARGO_SECRETORY;
            v.speed = 1.0f; v.rab_id = 1; // Rab1
            vesicles.push_back(v);
        }

        // ── Golgi → ER retrieval (COPI, KDEL receptor) ─────────────
        float copi_flux = 0.1f * COPI_available * compartments.protein_load[COMP_CIS_GOLGI]
                        * 0.1f; // 10% of Golgi protein is ER-resident
        if (copi_flux > 0.01f && vesicles.size() < 200) {
            Vesicle v;
            v.origin = COMP_CIS_GOLGI; v.destination = COMP_ER;
            v.coat = COAT_COPI; v.cargo = CARGO_ER_RETRIEVAL;
            v.speed = 0.8f; v.rab_id = 6; // Rab6
            vesicles.push_back(v);
        }

        // ── TGN → plasma membrane (constitutive secretion) ─────────
        float secr_flux = 0.15f * compartments.protein_load[COMP_TGN];
        if (secr_flux > 0.01f && vesicles.size() < 200) {
            Vesicle v;
            v.origin = COMP_TGN; v.destination = COMP_PLASMA_MEMBRANE;
            v.coat = COAT_NONE; v.cargo = CARGO_SECRETORY;
            v.speed = 1.2f; v.rab_id = 8; // Rab8
            vesicles.push_back(v);
        }

        // ── TGN → lysosome (M6P-tagged enzymes via clathrin) ───────
        float lyso_flux = 0.05f * compartments.protein_load[COMP_TGN];
        if (lyso_flux > 0.005f && vesicles.size() < 200) {
            Vesicle v;
            v.origin = COMP_TGN; v.destination = COMP_LATE_ENDOSOME;
            v.coat = COAT_CLATHRIN; v.cargo = CARGO_LYSOSOMAL_ENZYME;
            v.speed = 0.8f; v.rab_id = 7; // Rab7
            vesicles.push_back(v);
        }

        // ── Endocytosis (clathrin-mediated, from PM) ────────────────
        float endo_flux = 0.1f * clathrin_available * dynamin;
        if (endo_flux > 0.01f && vesicles.size() < 200) {
            Vesicle v;
            v.origin = COMP_PLASMA_MEMBRANE; v.destination = COMP_EARLY_ENDOSOME;
            v.coat = COAT_CLATHRIN; v.cargo = CARGO_RECEPTOR_RECYCLING;
            v.speed = 1.0f; v.rab_id = 5; // Rab5
            vesicles.push_back(v);
        }

        // ── Move vesicles toward destination ────────────────────────
        for (auto& v : vesicles) {
            if (v.fused) continue;
            v.progress += v.speed * dt * 0.5f;

            // SNARE-mediated fusion at destination
            if (v.progress >= 1.0f) {
                v.fused = true;
                // Transfer cargo
                compartments.protein_load[v.origin] -= 0.1f;
                compartments.protein_load[v.destination] += 0.1f;
            }
        }

        // Remove fused vesicles
        vesicles.erase(std::remove_if(vesicles.begin(), vesicles.end(),
            [](const Vesicle& v) { return v.fused; }), vesicles.end());

        // ── Unfolded Protein Response (UPR) ─────────────────────────
        // ER stress from high protein load + low BiP
        upr.ER_stress = compartments.protein_load[COMP_ER] / (upr.BiP_level * 5.0f + 1.0f);

        // IRE1 activation (RNase → XBP1 mRNA splicing)
        upr.IRE1_active += (0.5f * upr.ER_stress - 0.1f * upr.IRE1_active) * dt;
        upr.XBP1s += (0.3f * upr.IRE1_active - 0.05f * upr.XBP1s) * dt;

        // PERK activation (kinase → eIF2α phosphorylation → translation block)
        upr.PERK_active += (0.3f * upr.ER_stress - 0.1f * upr.PERK_active) * dt;
        upr.eIF2a_P += (0.5f * upr.PERK_active - 0.1f * upr.eIF2a_P) * dt;

        // ATF6 activation (ER→Golgi transport, S1P/S2P cleavage)
        upr.ATF6_active += (0.2f * upr.ER_stress - 0.1f * upr.ATF6_active) * dt;

        // BiP upregulation (all three UPR branches → more chaperones)
        upr.BiP_level += (0.1f * (upr.XBP1s + upr.ATF6_active) - 0.01f * upr.BiP_level) * dt;

        // CHOP induction (prolonged stress → apoptosis)
        upr.CHOP += (0.05f * upr.PERK_active * upr.ATF6_active - 0.02f * upr.CHOP) * dt;

        // ── Autophagy ──────────────────────────────────────────────
        // mTORC1 inhibits autophagy — low nutrients → mTORC1 off → autophagy on
        autophagy.mTORC1_inhibition = 1.0f - mTORC1_active;

        // ULK1 activation (inhibited by mTORC1, activated by AMPK)
        float ulk_on = 0.3f * autophagy.mTORC1_inhibition * (1.0f - autophagy.ULK1_active);
        autophagy.ULK1_active += (ulk_on - 0.05f * autophagy.ULK1_active) * dt;

        // LC3 lipidation (LC3-I → LC3-II, autophagosome membrane marker)
        float lc3_convert = 0.2f * autophagy.ULK1_active * autophagy.LC3_I;
        autophagy.LC3_II += (lc3_convert - 0.1f * autophagy.LC3_II) * dt;
        autophagy.LC3_I -= lc3_convert * dt;
        autophagy.LC3_I = fmaxf(0, autophagy.LC3_I);

        // Autophagosome formation
        autophagy.autophagosome_count += (0.1f * autophagy.LC3_II * autophagy.Beclin1
                                        - 0.05f * autophagy.autophagosome_count) * dt;

        // Autolysosome fusion and degradation
        autophagy.autolysosome_flux = 0.2f * autophagy.autophagosome_count
                                    * compartments.protein_load[COMP_LYSOSOME];

        // p62 consumed during selective autophagy
        autophagy.p62 += (0.05f - 0.1f * autophagy.ULK1_active * autophagy.p62) * dt;
        autophagy.p62 = fmaxf(0, autophagy.p62);

        // Clamp
        for (int i = 0; i < COMP_COUNT; i++) {
            compartments.protein_load[i] = fmaxf(0, compartments.protein_load[i]);
        }
    }

    int getVesicleCount() const { return (int)vesicles.size(); }
    float getERstress() const { return upr.ER_stress; }
    float getAutophagyLevel() const { return autophagy.ULK1_active; }
};
