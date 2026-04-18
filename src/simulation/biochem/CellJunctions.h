#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  CellJunctions.h — Cell-cell and cell-matrix adhesion
//  Based on: Alberts Ch 19 "Cell Junctions and the Extracellular Matrix"
//
//  Implements:
//    - Cadherins: E-cadherin, N-cadherin (homophilic, Ca2+-dependent)
//      + catenins (α, β, p120) linking to actin cytoskeleton
//    - Desmosomes: desmoglein/desmocollin + desmoplakin → IF
//    - Tight junctions: claudins, occludin, ZO-1 (barrier + fence)
//    - Gap junctions: connexins → connexons → ion/metabolite coupling
//    - Integrins: αβ heterodimers (inside-out, outside-in signaling)
//      + focal adhesions (FAK, paxillin, vinculin, talin)
//    - ECM components: collagen, fibronectin, laminin, proteoglycans
//    - Mechanotransduction: force → adhesion strengthening
//    - Selectins: transient adhesion (leukocyte rolling)
//
//  Units: relative levels (0-1 normalized), forces in pN
// ══════════════════════════════════════════════════════════════════════════

struct CadherinJunction {
    float E_cadherin = 1.0f;         // Epithelial cadherin
    float N_cadherin = 0.0f;         // Neural cadherin
    float P_cadherin = 0.0f;         // Placental cadherin
    float beta_catenin_junc = 0.5f;  // β-catenin at junction (not Wnt pool!)
    float alpha_catenin = 0.5f;      // α-catenin (actin linker)
    float p120_catenin = 0.5f;       // p120-catenin (cadherin stability)
    float vinculin_junc = 0.0f;      // Vinculin (recruited under tension)

    float adhesion_strength = 0.0f;  // Total adhesive force (pN)
    float contact_area = 0.0f;       // Cell-cell contact area (µm²)

    void update(float Ca_extracellular, float tension, float dt) {
        // Cadherin requires Ca2+ for rigid conformation
        float ca_factor = Ca_extracellular / (1.0f + Ca_extracellular);

        // Adhesion strength = cadherins × catenins × Ca2+
        adhesion_strength = (E_cadherin + N_cadherin) * ca_factor
                          * beta_catenin_junc * alpha_catenin * 100.0f;

        // Mechanosensing: tension exposes vinculin-binding site on α-catenin
        if (tension > 5.0f) {
            vinculin_junc += 0.1f * (tension / 50.0f) * dt;
            vinculin_junc = fminf(1.0f, vinculin_junc);
            // Vinculin strengthens junction
            adhesion_strength *= (1.0f + vinculin_junc * 0.5f);
        }
    }
};

struct Desmosome {
    float desmoglein = 0.5f;         // Desmosomal cadherin
    float desmocollin = 0.5f;        // Desmosomal cadherin
    float plakoglobin = 0.5f;        // γ-catenin
    float desmoplakin = 0.5f;        // Links to intermediate filaments
    float plakophilin = 0.3f;        // Plaque protein

    float mechanical_strength = 0.0f;

    void update() {
        mechanical_strength = desmoglein * desmocollin * desmoplakin * 200.0f;
    }
};

struct TightJunction {
    float claudin_1 = 0.5f;         // Barrier claudin
    float claudin_2 = 0.1f;         // Pore-forming claudin (cation pore)
    float occludin = 0.5f;          // Barrier protein
    float ZO_1 = 0.5f;             // Scaffold (links to actin)
    float ZO_2 = 0.3f;
    float JAM_A = 0.3f;            // Junctional adhesion molecule

    float barrier_integrity = 0.0f;  // 0 = leaky, 1 = tight seal
    float paracellular_permeability = 0.0f;

    void update() {
        barrier_integrity = (claudin_1 + occludin) * ZO_1 / 2.0f;
        barrier_integrity = std::clamp(barrier_integrity, 0.0f, 1.0f);
        // Claudin-2 creates cation-selective pores (increases permeability)
        paracellular_permeability = claudin_2 * 0.5f * (1.0f - barrier_integrity * 0.5f);
    }
};

struct GapJunction {
    float connexin43 = 0.5f;        // Cx43 (most common, heart, etc.)
    float connexin26 = 0.1f;        // Cx26 (inner ear)
    float connexin32 = 0.1f;        // Cx32 (liver)
    float connexon_count = 50.0f;   // Number of connexon channels

    float ion_coupling = 0.0f;      // Electrical coupling coefficient
    float metabolite_transfer = 0.0f;// Small molecule transfer rate

    // Gating: pH, Ca2+, voltage dependent
    float gate_open_fraction = 0.8f;

    void update(float pH, float Ca_cytosol, float voltage_diff) {
        // Low pH closes gap junctions
        float ph_gating = pH > 6.5f ? 1.0f : (pH - 5.5f) / 1.0f;
        // High Ca2+ closes gap junctions
        float ca_gating = 1.0f / (1.0f + Ca_cytosol / 0.5f);
        // Voltage gating
        float v_gating = 1.0f / (1.0f + fabsf(voltage_diff) / 50.0f);

        gate_open_fraction = ph_gating * ca_gating * v_gating;
        ion_coupling = connexon_count * gate_open_fraction * 0.01f;
        metabolite_transfer = connexon_count * gate_open_fraction * 0.005f;
    }
};

struct IntegrinAdhesion {
    // Integrin heterodimers
    float alpha5_beta1 = 0.5f;       // Fibronectin receptor
    float alpha6_beta4 = 0.3f;       // Laminin receptor (hemidesmosomes)
    float alphaV_beta3 = 0.2f;       // Vitronectin/fibronectin
    float alpha2_beta1 = 0.2f;       // Collagen receptor

    // Activation state (bent=inactive, extended=active)
    float integrin_active = 0.0f;    // Fraction in active conformation

    // Focal adhesion components
    float talin = 0.5f;             // Inside-out activation, actin linker
    float vinculin_FA = 0.0f;       // Recruited under tension
    float paxillin = 0.3f;          // Scaffold
    float FAK = 0.3f;               // Focal Adhesion Kinase
    float FAK_active = 0.0f;        // Autophosphorylated FAK
    float Src_active = 0.0f;        // Src kinase (recruited by FAK)
    float zyxin = 0.0f;             // Tension-sensitive

    // ECM binding
    float fibronectin_bound = 0.0f;
    float collagen_bound = 0.0f;
    float laminin_bound = 0.0f;

    // Outputs
    float adhesion_force = 0.0f;    // pN
    float survival_signal = 0.0f;   // Anoikis prevention

    void update(float tension, float talin_signal, float dt) {
        // Inside-out activation: talin binding to β-tail
        float activation = talin_signal * talin;
        integrin_active += (activation * (1.0f - integrin_active) - 0.1f * integrin_active) * dt;

        // ECM binding (proportional to active integrins × ECM)
        fibronectin_bound = integrin_active * alpha5_beta1;
        collagen_bound = integrin_active * alpha2_beta1;
        laminin_bound = integrin_active * alpha6_beta4;

        float total_bound = fibronectin_bound + collagen_bound + laminin_bound;

        // Focal adhesion assembly (force-dependent maturation)
        if (tension > 2.0f) {
            vinculin_FA += 0.1f * tension * 0.01f * dt;
            zyxin += 0.05f * tension * 0.01f * dt;
            vinculin_FA = fminf(1.0f, vinculin_FA);
            zyxin = fminf(1.0f, zyxin);
        }

        // FAK autophosphorylation (upon integrin clustering)
        FAK_active = FAK * total_bound / (0.3f + total_bound);
        Src_active = 0.5f * FAK_active;

        // Adhesion force
        adhesion_force = total_bound * (1.0f + vinculin_FA + zyxin) * 50.0f;

        // Survival signal (prevents anoikis)
        survival_signal = FAK_active * 0.3f + Src_active * 0.2f;
    }
};

// ══════════════════════════════════════════════════════════════════════════

class CellJunctionsEngine {
public:
    CadherinJunction adherens;
    Desmosome desmosome;
    TightJunction tightJunc;
    GapJunction gapJunc;
    IntegrinAdhesion integrin;

    // Selectins (for immune cell rolling)
    float E_selectin = 0.0f;         // On endothelium (inflammation-induced)
    float P_selectin = 0.0f;         // On platelets/endothelium
    float L_selectin = 0.0f;         // On leukocytes

    void init() {
        adherens = CadherinJunction();
        desmosome = Desmosome();
        tightJunc = TightJunction();
        gapJunc = GapJunction();
        integrin = IntegrinAdhesion();
    }

    void step(float dt, float Ca_ext, float tension, float pH, float Ca_cyt) {
        adherens.update(Ca_ext, tension, dt);
        desmosome.update();
        tightJunc.update();
        gapJunc.update(pH, Ca_cyt, 0);
        integrin.update(tension, 0.5f, dt);
    }

    float totalAdhesion() const {
        return adherens.adhesion_strength + desmosome.mechanical_strength +
               integrin.adhesion_force;
    }
};
