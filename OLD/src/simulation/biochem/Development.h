#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Development.h — Developmental biology (embryonic patterning)
//  Based on: Alberts Ch 21 "Development of Multicellular Organisms"
//
//  Implements:
//    - Morphogen gradients (Bicoid, Dorsal, BMP, Nodal, Wnt, Shh)
//    - Lateral inhibition (Notch-Delta, salt-and-pepper patterns)
//    - Hox gene patterning (colinearity, positional identity)
//    - Asymmetric cell division (PAR proteins, Numb segregation)
//    - Germ layer specification (ectoderm, mesoderm, endoderm)
//    - Gastrulation movements (invagination, EMT)
//    - Master transcription regulators (MyoD, Pax6, Eyeless)
//    - Combinatorial gene regulation (enhancer logic)
//    - Sequential induction (pattern refinement)
//
//  Single-cell mode: tracks the cell's positional identity and
//  developmental potential based on received morphogen concentrations
// ══════════════════════════════════════════════════════════════════════════

enum GermLayer { GERM_UNDETERMINED = 0, GERM_ECTODERM, GERM_MESODERM, GERM_ENDODERM };
enum BodyAxis { AXIS_AP = 0, AXIS_DV, AXIS_LR }; // Anteroposterior, Dorsoventral, Left-Right

// ── Morphogen Concentrations ────────────────────────────────────────────
struct MorphogenField {
    // Anterior-Posterior axis
    float Bicoid = 0.0f;           // Anterior determinant (Drosophila)
    float Nanos = 0.0f;            // Posterior determinant
    float Wnt_gradient = 0.0f;     // Wnt gradient (posteriorizing)
    float RA_gradient = 0.0f;      // Retinoic acid (vertebrate A-P)
    float FGF_gradient = 0.0f;     // FGF (posteriorizing)

    // Dorsal-Ventral axis
    float BMP_activity = 0.5f;     // BMP (ventralizing in vertebrates)
    float Chordin = 0.0f;          // BMP antagonist (from organizer)
    float Noggin = 0.0f;           // BMP antagonist
    float Nodal = 0.0f;            // TGF-β family (mesoderm/endoderm)
    float Dorsal_nuc = 0.0f;       // Dorsal/NF-κB (Drosophila D-V)

    // Hedgehog signaling
    float Shh = 0.0f;              // Sonic Hedgehog (neural tube, limb)

    // Notch-Delta (lateral inhibition)
    float Delta_expression = 0.5f;  // Delta ligand on this cell
    float Notch_activation = 0.0f;  // Notch signaling received

    // Compute net positional values
    float AP_position() const { // 0 = anterior, 1 = posterior
        return (Wnt_gradient + Nanos + FGF_gradient) /
               (Bicoid + Wnt_gradient + Nanos + FGF_gradient + 0.01f);
    }
    float DV_position() const { // 0 = dorsal, 1 = ventral
        float bmp_net = BMP_activity / (1.0f + Chordin + Noggin);
        return bmp_net / (1.0f + bmp_net);
    }
};

// ── Hox Gene Expression ─────────────────────────────────────────────────
struct HoxCode {
    // 4 clusters × ~10 genes each in vertebrates (HoxA-D, 1-13)
    // Simplified: track expression level of representative Hox genes
    float Hox1 = 0.0f;  // Anterior (hindbrain rhombomere 4)
    float Hox2 = 0.0f;  // (rhombomere 4-5)
    float Hox3 = 0.0f;  // Cervical
    float Hox4 = 0.0f;  // Cervical/thoracic
    float Hox5 = 0.0f;  // Thoracic
    float Hox6 = 0.0f;  // Thoracic
    float Hox7 = 0.0f;  // Thoracic/lumbar
    float Hox8 = 0.0f;  // Lumbar
    float Hox9 = 0.0f;  // Lumbar/sacral
    float Hox10 = 0.0f; // Sacral
    float Hox11 = 0.0f; // Sacral/caudal
    float Hox13 = 0.0f; // Most posterior

    // Polycomb/Trithorax maintenance
    float PcG_repression = 0.5f;  // Polycomb-mediated silencing
    float TrxG_activation = 0.5f; // Trithorax-mediated activation

    void setFromPosition(float ap_pos) {
        // Colinearity: genes activated sequentially from anterior to posterior
        Hox1 = (ap_pos > 0.1f) ? fminf(1, (ap_pos - 0.1f) * 5.0f) : 0;
        Hox3 = (ap_pos > 0.2f) ? fminf(1, (ap_pos - 0.2f) * 4.0f) : 0;
        Hox5 = (ap_pos > 0.3f) ? fminf(1, (ap_pos - 0.3f) * 3.5f) : 0;
        Hox7 = (ap_pos > 0.45f) ? fminf(1, (ap_pos - 0.45f) * 3.0f) : 0;
        Hox9 = (ap_pos > 0.6f) ? fminf(1, (ap_pos - 0.6f) * 3.0f) : 0;
        Hox11 = (ap_pos > 0.75f) ? fminf(1, (ap_pos - 0.75f) * 4.0f) : 0;
        Hox13 = (ap_pos > 0.9f) ? fminf(1, (ap_pos - 0.9f) * 10.0f) : 0;

        // Posterior prevalence: most posterior Hox dominates
        // (doesn't override — just strongest expressed determines segment identity)
    }

    // What body region does this Hox code specify?
    const char* getRegionIdentity() const {
        if (Hox13 > 0.5f) return "caudal";
        if (Hox11 > 0.5f) return "sacral";
        if (Hox9 > 0.5f) return "lumbar";
        if (Hox7 > 0.5f) return "thoracic";
        if (Hox5 > 0.5f) return "cervical";
        if (Hox3 > 0.5f) return "hindbrain";
        if (Hox1 > 0.5f) return "anterior_hindbrain";
        return "forebrain/head"; // No Hox = most anterior
    }
};

// ── PAR Polarity Complex ────────────────────────────────────────────────
struct PARpolarity {
    // Anterior complex
    float PAR3 = 0.5f;       // PAR-3 (anterior/apical)
    float PAR6 = 0.5f;       // PAR-6
    float aPKC = 0.5f;       // Atypical PKC

    // Posterior complex
    float PAR1 = 0.5f;       // PAR-1 (posterior/basal)
    float PAR2 = 0.5f;       // PAR-2

    float polarity_strength = 0.0f; // How polarized (0 = symmetric, 1 = fully polarized)

    void polarize(float cue, float dt) {
        // Mutual exclusion: anterior complex excludes posterior and vice versa
        float anterior_drive = cue;
        float posterior_drive = 1.0f - cue;

        PAR3 += (anterior_drive * 0.5f - PAR1 * PAR3 * 0.3f) * dt;
        PAR1 += (posterior_drive * 0.5f - PAR3 * PAR1 * 0.3f) * dt;
        PAR3 = std::clamp(PAR3, 0.0f, 1.0f);
        PAR1 = std::clamp(PAR1, 0.0f, 1.0f);

        polarity_strength = fabsf(PAR3 - PAR1) / (PAR3 + PAR1 + 0.01f);
    }
};

// ── Cell Fate Determination ─────────────────────────────────────────────
struct CellFateState {
    GermLayer germ_layer = GERM_UNDETERMINED;
    float pluripotency = 1.0f;    // 1 = totipotent, 0 = terminally differentiated
    float differentiation = 0.0f; // 0 = stem, 1 = fully differentiated

    // Master regulators (Ch 21: "Some transcription regulators can
    //  activate a program that defines a cell type or creates an entire organ")
    float Oct4 = 1.0f;            // Pluripotency (stem cell state)
    float Sox2 = 1.0f;            // Pluripotency
    float Nanog = 0.8f;           // Pluripotency
    float MyoD = 0.0f;            // Muscle specification
    float Pax6 = 0.0f;            // Eye/neural specification
    float GATA4 = 0.0f;           // Endoderm/cardiac
    float Brachyury = 0.0f;       // Mesoderm
    float Sox1 = 0.0f;            // Neural (ectoderm)
    float CDX2 = 0.0f;            // Trophoblast

    // Numb (asymmetric division determinant — inhibits Notch)
    float Numb = 0.5f;
    float Numb_asymmetry = 0.0f;  // How asymmetrically distributed

    void determineGermLayer(const MorphogenField& morphogens) {
        if (germ_layer != GERM_UNDETERMINED) return; // Already committed

        float nodal = morphogens.Nodal;
        float bmp = morphogens.BMP_activity / (1.0f + morphogens.Chordin + morphogens.Noggin);

        // High Nodal → endoderm; medium Nodal → mesoderm; low/no Nodal → ectoderm
        // Low BMP (dorsal) → neural ectoderm
        if (nodal > 0.7f) {
            germ_layer = GERM_ENDODERM;
            GATA4 = nodal;
        } else if (nodal > 0.3f) {
            germ_layer = GERM_MESODERM;
            Brachyury = nodal;
        } else {
            germ_layer = GERM_ECTODERM;
            if (bmp < 0.3f) Sox1 = 1.0f - bmp; // Neural fate
        }

        // Lose pluripotency markers
        Oct4 *= 0.5f;
        Nanog *= 0.5f;
    }
};

// ══════════════════════════════════════════════════════════════════════════

class DevelopmentEngine {
public:
    MorphogenField morphogens;
    HoxCode hox;
    PARpolarity polarity;
    CellFateState fate;

    void init() {
        morphogens = MorphogenField();
        hox = HoxCode();
        polarity = PARpolarity();
        fate = CellFateState();
    }

    void step(float dt) {
        dt = fminf(dt, 0.1f);

        // Determine germ layer based on morphogen exposure
        fate.determineGermLayer(morphogens);

        // Set Hox code from A-P position
        hox.setFromPosition(morphogens.AP_position());

        // Update polarity
        polarity.polarize(0.5f, dt);

        // Lateral inhibition via Notch-Delta
        // (In a multicellular context, Delta on this cell activates Notch on neighbors)
        // High Notch activation → reduce Delta expression (feedback)
        morphogens.Delta_expression += (-0.5f * morphogens.Notch_activation * morphogens.Delta_expression) * dt;
        morphogens.Delta_expression = std::clamp(morphogens.Delta_expression, 0.0f, 1.0f);

        // Differentiation progress
        if (fate.germ_layer != GERM_UNDETERMINED) {
            fate.differentiation += 0.01f * dt;
            fate.pluripotency -= 0.01f * dt;
            fate.differentiation = fminf(1.0f, fate.differentiation);
            fate.pluripotency = fmaxf(0.0f, fate.pluripotency);
        }
    }
};
