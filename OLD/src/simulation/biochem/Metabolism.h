#pragma once
#include <cmath>
#include <cstdio>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Metabolism.h — Complete cellular metabolism simulation
//  Based on: Alberts "Molecular Biology of the Cell" 7th Ed, Ch 2 & 14
//
//  Implements:
//    - Glycolysis (10 enzymatic steps, Embden-Meyerhof pathway)
//    - Pyruvate dehydrogenase complex
//    - Citric acid cycle (Krebs cycle, 8 steps)
//    - Oxidative phosphorylation (ETC Complexes I-IV + ATP synthase)
//    - Fermentation (lactic acid, ethanol pathways)
//    - β-oxidation of fatty acids
//    - Pentose phosphate pathway (oxidative branch)
//
//  All enzymes modeled with Michaelis-Menten kinetics:
//    v = Vmax * [S] / (Km + [S])
//  With allosteric regulation where applicable.
//
//  Kinetic parameters from BRENDA enzyme database (brenda-enzymes.org)
//  and published metabolomics literature.
//
//  Units: concentrations in mM, time in seconds, energy in kJ/mol
// ══════════════════════════════════════════════════════════════════════════

// ── Physical/Chemical Constants ─────────────────────────────────────────
static constexpr float R_GAS   = 8.314e-3f; // kJ/(mol·K)
static constexpr float TEMP_K  = 310.15f;   // 37°C in Kelvin
static constexpr float RT      = R_GAS * TEMP_K; // ~2.577 kJ/mol
static constexpr float FARADAY = 96.485f;   // kJ/(mol·V)
static constexpr float LN10    = 2.302585f;

// ── Standard Free Energies of Hydrolysis (kJ/mol) ──────────────────────
// Under standard cellular conditions (not standard state)
static constexpr float DG_ATP_HYDROLYSIS  = -54.0f;  // ATP → ADP + Pi
static constexpr float DG_GTP_HYDROLYSIS  = -54.0f;  // GTP → GDP + Pi
static constexpr float DG_NADH_OXIDATION  = -220.0f; // NADH → NAD+ (to O2)
static constexpr float DG_FADH2_OXIDATION = -181.0f; // FADH2 → FAD

// ── Metabolite Pool ─────────────────────────────────────────────────────
// All concentrations in mM (millimolar), typical mammalian cell values
// References: Park et al. 2016 Nature Chem Biol; Lehninger Biochemistry

struct MetabolitePool {
    // ── Glycolysis intermediates ────────────────────────────────────
    float glucose          = 5.0f;   // Extracellular: ~5mM (blood glucose)
    float glucose_6P       = 0.08f;  // Glucose-6-phosphate
    float fructose_6P      = 0.014f; // Fructose-6-phosphate
    float fructose_16BP    = 0.03f;  // Fructose-1,6-bisphosphate
    float DHAP             = 0.16f;  // Dihydroxyacetone phosphate
    float G3P              = 0.019f; // Glyceraldehyde-3-phosphate
    float BPG_13           = 0.001f; // 1,3-Bisphosphoglycerate
    float PG_3             = 0.12f;  // 3-Phosphoglycerate
    float PG_2             = 0.03f;  // 2-Phosphoglycerate
    float PEP              = 0.023f; // Phosphoenolpyruvate
    float pyruvate         = 0.08f;  // Pyruvate

    // ── Citric Acid Cycle intermediates ─────────────────────────────
    float acetyl_CoA       = 0.01f;  // Acetyl-CoA (mito matrix)
    float oxaloacetate     = 0.011f; // Oxaloacetate
    float citrate          = 0.20f;  // Citrate
    float isocitrate       = 0.02f;  // Isocitrate
    float alpha_KG         = 0.10f;  // α-Ketoglutarate
    float succinyl_CoA     = 0.04f;  // Succinyl-CoA
    float succinate        = 0.30f;  // Succinate
    float fumarate         = 0.05f;  // Fumarate
    float malate           = 0.15f;  // Malate
    float CoA_free         = 0.05f;  // Free Coenzyme A

    // ── Energy carriers ────────────────────────────────────────────
    float ATP              = 3.0f;   // Adenosine triphosphate
    float ADP              = 0.8f;   // Adenosine diphosphate
    float AMP              = 0.1f;   // Adenosine monophosphate
    float GTP              = 0.5f;   // Guanosine triphosphate
    float GDP              = 0.1f;   // Guanosine diphosphate
    float Pi               = 5.0f;   // Inorganic phosphate

    // ── Electron carriers ──────────────────────────────────────────
    float NAD_plus         = 1.5f;   // NAD+ (oxidized)
    float NADH             = 0.1f;   // NADH (reduced)
    float NADP_plus        = 0.002f; // NADP+ (oxidized)
    float NADPH            = 0.05f;  // NADPH (reduced)
    float FAD              = 0.1f;   // FAD (oxidized)
    float FADH2            = 0.01f;  // FADH2 (reduced)

    // ── Proton gradient (mitochondria) ─────────────────────────────
    float H_plus_matrix    = 1e-5f;  // [H+] in matrix (~pH 7.8 → ~1.6e-5 mM)
    float H_plus_IMS       = 6.3e-5f;// [H+] in intermembrane space (~pH 7.2)
    float membrane_potential = 0.18f; // ΔΨ in volts (~180 mV, matrix negative)

    // ── Other metabolites ──────────────────────────────────────────
    float O2               = 0.03f;  // Dissolved O2 (~30 µM intracellular)
    float CO2              = 1.2f;   // Dissolved CO2
    float H2O              = 55000.0f; // Water (constant, not tracked dynamically)
    float lactate          = 1.0f;   // Lactate
    float ethanol          = 0.0f;   // Ethanol (not normally in human cells)

    // ── Fatty acid metabolism ──────────────────────────────────────
    float palmitoyl_CoA    = 0.01f;  // Palmitoyl-CoA (C16 fatty acid)
    float fatty_acid_free  = 0.05f;  // Free fatty acids

    // ── Amino acid pool (simplified) ───────────────────────────────
    float amino_acid_pool  = 2.0f;   // Generic amino acid pool

    // ── Pentose phosphate pathway ──────────────────────────────────
    float glucose_6P_DH_active = 1.0f; // G6PDH activity level
    float ribose_5P        = 0.05f;  // Ribose-5-phosphate (for nucleotide synthesis)

    // ── Derived quantities ─────────────────────────────────────────
    float energyCharge() const {
        // Atkinson energy charge: (ATP + 0.5*ADP) / (ATP + ADP + AMP)
        float total = ATP + ADP + AMP;
        if (total < 1e-6f) return 0.5f;
        return (ATP + 0.5f * ADP) / total;
    }

    float NAD_ratio() const {
        // [NADH]/[NAD+] ratio — indicates redox state
        if (NAD_plus < 1e-6f) return 100.0f;
        return NADH / NAD_plus;
    }

    float pH_cytosol() const {
        // Approximate cytosolic pH from proton concentration
        return 7.2f; // Normally tightly regulated
    }

    float pH_matrix() const {
        if (H_plus_matrix < 1e-9f) return 8.0f;
        return -log10f(H_plus_matrix * 1e-3f); // Convert mM to M for pH
    }

    float protonMotiveForce() const {
        // PMF = ΔΨ - (2.303 * RT/F) * ΔpH
        float dpH = pH_matrix() - (-log10f(H_plus_IMS * 1e-3f));
        return membrane_potential - (2.303f * RT / FARADAY) * dpH;
    }

    float totalAdenylate() const { return ATP + ADP + AMP; }
};

// ── Enzyme Kinetics Helper ──────────────────────────────────────────────

static float michaelisMenten(float substrate, float Vmax, float Km) {
    return Vmax * substrate / (Km + substrate);
}

// Michaelis-Menten with competitive inhibitor
static float mmCompetitiveInhibition(float S, float I, float Vmax, float Km, float Ki) {
    return Vmax * S / (Km * (1.0f + I / Ki) + S);
}

// Hill equation for cooperative binding
static float hillEquation(float S, float Vmax, float K05, float n) {
    float Sn = powf(S, n);
    return Vmax * Sn / (powf(K05, n) + Sn);
}

// Reversible Michaelis-Menten (Haldane equation)
static float reversibleMM(float S, float P, float Vmax_f, float Km_S,
                           float Vmax_r, float Km_P) {
    float num = Vmax_f * S / Km_S - Vmax_r * P / Km_P;
    float denom = 1.0f + S / Km_S + P / Km_P;
    return num / denom;
}

// Clamp a value to prevent negative concentrations
static float clampPositive(float x) { return x > 0 ? x : 0; }

// ══════════════════════════════════════════════════════════════════════════
//  GLYCOLYSIS — 10 enzymatic steps
//  Glucose → 2 Pyruvate + 2 ATP(net) + 2 NADH
//  Location: Cytosol
//
//  Ref: Alberts Ch 2 "Stage 1 (Cytosol): Glycolysis"
//       Berg, Tymoczko & Stryer "Biochemistry" 9th Ed
// ══════════════════════════════════════════════════════════════════════════

struct GlycolysisEnzymes {
    // ── Step 1: Hexokinase (EC 2.7.1.1) ────────────────────────────
    // Glucose + ATP → Glucose-6-P + ADP
    // Irreversible, committed step
    // Inhibited by product (glucose-6-P)
    static constexpr float HK_Vmax  = 1.0f;    // mM/s (normalized)
    static constexpr float HK_Km_glucose = 0.1f; // mM
    static constexpr float HK_Km_ATP = 0.4f;     // mM
    static constexpr float HK_Ki_G6P = 0.5f;     // mM (product inhibition)

    static float hexokinase(const MetabolitePool& m) {
        // Bi-substrate with product inhibition
        float v = HK_Vmax * m.glucose / (HK_Km_glucose + m.glucose)
                          * m.ATP / (HK_Km_ATP + m.ATP);
        // Product inhibition by glucose-6-phosphate
        v *= HK_Ki_G6P / (HK_Ki_G6P + m.glucose_6P);
        return v;
    }

    // ── Step 2: Phosphoglucose Isomerase (EC 5.3.1.9) ──────────────
    // Glucose-6-P ⇌ Fructose-6-P
    // Near-equilibrium, reversible
    static constexpr float PGI_Vmax_f = 10.0f;  // Fast, near equilibrium
    static constexpr float PGI_Km_G6P = 0.5f;
    static constexpr float PGI_Vmax_r = 8.0f;
    static constexpr float PGI_Km_F6P = 0.15f;

    static float phosphoglucoseIsomerase(const MetabolitePool& m) {
        return reversibleMM(m.glucose_6P, m.fructose_6P,
                            PGI_Vmax_f, PGI_Km_G6P,
                            PGI_Vmax_r, PGI_Km_F6P);
    }

    // ── Step 3: Phosphofructokinase-1 (PFK-1) (EC 2.7.1.11) ────────
    // Fructose-6-P + ATP → Fructose-1,6-BP + ADP
    // THE major regulatory step of glycolysis
    // Allosteric activation by: AMP, ADP, fructose-2,6-BP
    // Allosteric inhibition by: ATP, citrate
    static constexpr float PFK_Vmax = 0.8f;
    static constexpr float PFK_Km_F6P = 0.1f;
    static constexpr float PFK_Km_ATP = 0.05f;
    static constexpr float PFK_Ki_ATP = 1.0f;    // High ATP inhibits
    static constexpr float PFK_Ka_AMP = 0.05f;   // AMP activates
    static constexpr float PFK_Ki_citrate = 0.5f; // Citrate inhibits

    static float phosphofructokinase1(const MetabolitePool& m) {
        // Base rate
        float v = PFK_Vmax * m.fructose_6P / (PFK_Km_F6P + m.fructose_6P)
                           * m.ATP / (PFK_Km_ATP + m.ATP);
        // ATP inhibition (high energy charge slows glycolysis)
        float atp_inhibition = PFK_Ki_ATP / (PFK_Ki_ATP + m.ATP);
        // AMP activation (low energy charge accelerates glycolysis)
        float amp_activation = 1.0f + 4.0f * m.AMP / (PFK_Ka_AMP + m.AMP);
        // Citrate inhibition (TCA cycle intermediates inhibit glycolysis)
        float citrate_inhibition = PFK_Ki_citrate / (PFK_Ki_citrate + m.citrate);

        return v * atp_inhibition * amp_activation * citrate_inhibition;
    }

    // ── Step 4: Aldolase (EC 4.1.2.13) ─────────────────────────────
    // Fructose-1,6-BP → DHAP + G3P
    // Reversible
    static constexpr float ALD_Vmax_f = 3.0f;
    static constexpr float ALD_Km = 0.01f;

    static float aldolase(const MetabolitePool& m) {
        return michaelisMenten(m.fructose_16BP, ALD_Vmax_f, ALD_Km);
    }

    // ── Step 5: Triose Phosphate Isomerase (TPI) (EC 5.3.1.1) ──────
    // DHAP ⇌ G3P
    // "Perfect enzyme" — diffusion-limited, near equilibrium
    // Equilibrium strongly favors DHAP, but G3P is consumed rapidly
    static constexpr float TPI_Vmax = 100.0f;   // Very fast
    static constexpr float TPI_Km_DHAP = 0.5f;
    static constexpr float TPI_Keq = 0.0476f;   // [G3P]/[DHAP] at equilibrium

    static float triosephosphateIsomerase(const MetabolitePool& m) {
        // Net flux toward G3P
        return TPI_Vmax * (m.DHAP - m.G3P / TPI_Keq) / (TPI_Km_DHAP + m.DHAP + m.G3P);
    }

    // ── Step 6: Glyceraldehyde-3-P Dehydrogenase (GAPDH) (EC 1.2.1.12)
    // G3P + NAD+ + Pi → 1,3-BPG + NADH + H+
    // First energy-harvesting step (NADH production)
    // Requires NAD+ — blocked if NAD+ is depleted
    static constexpr float GAPDH_Vmax = 5.0f;
    static constexpr float GAPDH_Km_G3P = 0.005f;
    static constexpr float GAPDH_Km_NAD = 0.09f;
    static constexpr float GAPDH_Km_Pi = 3.0f;

    static float G3P_dehydrogenase(const MetabolitePool& m) {
        return GAPDH_Vmax * m.G3P / (GAPDH_Km_G3P + m.G3P)
                          * m.NAD_plus / (GAPDH_Km_NAD + m.NAD_plus)
                          * m.Pi / (GAPDH_Km_Pi + m.Pi);
    }

    // ── Step 7: Phosphoglycerate Kinase (PGK) (EC 2.7.2.3) ─────────
    // 1,3-BPG + ADP → 3-PG + ATP
    // First substrate-level phosphorylation (ATP production!)
    static constexpr float PGK_Vmax = 10.0f;
    static constexpr float PGK_Km_BPG = 0.002f;
    static constexpr float PGK_Km_ADP = 0.1f;

    static float phosphoglycerateKinase(const MetabolitePool& m) {
        return PGK_Vmax * m.BPG_13 / (PGK_Km_BPG + m.BPG_13)
                        * m.ADP / (PGK_Km_ADP + m.ADP);
    }

    // ── Step 8: Phosphoglycerate Mutase (PGM) (EC 5.4.2.11) ────────
    // 3-PG ⇌ 2-PG
    static constexpr float PGM_Vmax = 8.0f;
    static constexpr float PGM_Km = 0.3f;

    static float phosphoglycerateMutase(const MetabolitePool& m) {
        return reversibleMM(m.PG_3, m.PG_2, PGM_Vmax, PGM_Km, PGM_Vmax * 0.8f, PGM_Km * 0.5f);
    }

    // ── Step 9: Enolase (EC 4.2.1.11) ──────────────────────────────
    // 2-PG ⇌ PEP + H2O
    static constexpr float ENO_Vmax = 6.0f;
    static constexpr float ENO_Km = 0.1f;

    static float enolase(const MetabolitePool& m) {
        return michaelisMenten(m.PG_2, ENO_Vmax, ENO_Km);
    }

    // ── Step 10: Pyruvate Kinase (PK) (EC 2.7.1.40) ────────────────
    // PEP + ADP → Pyruvate + ATP
    // Second substrate-level phosphorylation
    // Allosteric: activated by fructose-1,6-BP, inhibited by ATP
    static constexpr float PK_Vmax = 2.0f;
    static constexpr float PK_Km_PEP = 0.2f;
    static constexpr float PK_Km_ADP = 0.3f;
    static constexpr float PK_Ka_F16BP = 0.01f; // Feedforward activation
    static constexpr float PK_Ki_ATP = 3.0f;

    static float pyruvateKinase(const MetabolitePool& m) {
        float v = PK_Vmax * m.PEP / (PK_Km_PEP + m.PEP)
                          * m.ADP / (PK_Km_ADP + m.ADP);
        // Feedforward activation by fructose-1,6-BP (product of PFK)
        float f16bp_activation = 1.0f + 3.0f * m.fructose_16BP / (PK_Ka_F16BP + m.fructose_16BP);
        // ATP inhibition
        float atp_inh = PK_Ki_ATP / (PK_Ki_ATP + m.ATP);
        return v * f16bp_activation * atp_inh;
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  PYRUVATE DEHYDROGENASE COMPLEX
//  Pyruvate + CoA + NAD+ → Acetyl-CoA + CO2 + NADH
//  Location: Mitochondrial matrix
//  Links glycolysis to TCA cycle
//  Ref: Alberts Ch 2 "Stage 2: Pyruvate Oxidation"
// ══════════════════════════════════════════════════════════════════════════

struct PyruvateDehydrogenaseComplex {
    static constexpr float PDH_Vmax = 0.5f;
    static constexpr float PDH_Km_pyruvate = 0.05f;
    static constexpr float PDH_Km_CoA = 0.01f;
    static constexpr float PDH_Km_NAD = 0.05f;
    static constexpr float PDH_Ki_acetylCoA = 0.03f; // Product inhibition
    static constexpr float PDH_Ki_NADH = 0.04f;      // Product inhibition

    static float rate(const MetabolitePool& m) {
        float v = PDH_Vmax * m.pyruvate / (PDH_Km_pyruvate + m.pyruvate)
                           * m.CoA_free / (PDH_Km_CoA + m.CoA_free)
                           * m.NAD_plus / (PDH_Km_NAD + m.NAD_plus);
        // Product inhibition by acetyl-CoA and NADH
        v *= PDH_Ki_acetylCoA / (PDH_Ki_acetylCoA + m.acetyl_CoA);
        v *= PDH_Ki_NADH / (PDH_Ki_NADH + m.NADH);
        return v;
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  CITRIC ACID CYCLE (KREBS CYCLE / TCA CYCLE)
//  Acetyl-CoA → 2 CO2 + 3 NADH + 1 FADH2 + 1 GTP
//  Location: Mitochondrial matrix
//  8 enzymatic steps
//  Ref: Alberts Ch 2 "Stage 2: The Citric Acid Cycle"
// ══════════════════════════════════════════════════════════════════════════

struct CitricAcidCycleEnzymes {
    // ── Step 1: Citrate Synthase (EC 2.3.3.1) ──────────────────────
    // Acetyl-CoA + Oxaloacetate + H2O → Citrate + CoA
    // Irreversible, regulated
    static constexpr float CS_Vmax = 0.4f;
    static constexpr float CS_Km_acetylCoA = 0.005f;
    static constexpr float CS_Km_OAA = 0.005f;
    static constexpr float CS_Ki_ATP = 5.0f;
    static constexpr float CS_Ki_NADH = 0.05f;
    static constexpr float CS_Ki_citrate = 0.5f;

    static float citrateSynthase(const MetabolitePool& m) {
        float v = CS_Vmax * m.acetyl_CoA / (CS_Km_acetylCoA + m.acetyl_CoA)
                          * m.oxaloacetate / (CS_Km_OAA + m.oxaloacetate);
        v *= CS_Ki_ATP / (CS_Ki_ATP + m.ATP);
        v *= CS_Ki_NADH / (CS_Ki_NADH + m.NADH);
        v *= CS_Ki_citrate / (CS_Ki_citrate + m.citrate);
        return v;
    }

    // ── Step 2: Aconitase (EC 4.2.1.3) ─────────────────────────────
    // Citrate ⇌ Isocitrate (via cis-aconitate)
    static constexpr float ACO_Vmax = 5.0f;
    static constexpr float ACO_Km = 0.5f;

    static float aconitase(const MetabolitePool& m) {
        return reversibleMM(m.citrate, m.isocitrate, ACO_Vmax, ACO_Km, ACO_Vmax * 0.7f, ACO_Km * 0.3f);
    }

    // ── Step 3: Isocitrate Dehydrogenase (IDH) (EC 1.1.1.41) ───────
    // Isocitrate + NAD+ → α-Ketoglutarate + CO2 + NADH
    // KEY regulatory step — allosteric
    // Activated by ADP, Ca2+; inhibited by ATP, NADH
    static constexpr float IDH_Vmax = 0.3f;
    static constexpr float IDH_Km_isocitrate = 0.02f;
    static constexpr float IDH_Km_NAD = 0.08f;
    static constexpr float IDH_Ka_ADP = 0.05f;
    static constexpr float IDH_Ki_ATP = 2.0f;
    static constexpr float IDH_Ki_NADH = 0.03f;

    static float isocitrateDehydrogenase(const MetabolitePool& m) {
        float v = IDH_Vmax * m.isocitrate / (IDH_Km_isocitrate + m.isocitrate)
                           * m.NAD_plus / (IDH_Km_NAD + m.NAD_plus);
        float adp_activation = 1.0f + 3.0f * m.ADP / (IDH_Ka_ADP + m.ADP);
        float atp_inh = IDH_Ki_ATP / (IDH_Ki_ATP + m.ATP);
        float nadh_inh = IDH_Ki_NADH / (IDH_Ki_NADH + m.NADH);
        return v * adp_activation * atp_inh * nadh_inh;
    }

    // ── Step 4: α-Ketoglutarate Dehydrogenase (EC 1.2.4.2) ────────
    // α-KG + CoA + NAD+ → Succinyl-CoA + CO2 + NADH
    // Similar regulation to PDH
    static constexpr float KGDH_Vmax = 0.25f;
    static constexpr float KGDH_Km_aKG = 0.1f;
    static constexpr float KGDH_Km_CoA = 0.01f;
    static constexpr float KGDH_Km_NAD = 0.05f;
    static constexpr float KGDH_Ki_succinylCoA = 0.02f;
    static constexpr float KGDH_Ki_NADH = 0.02f;

    static float alphaKGDehydrogenase(const MetabolitePool& m) {
        float v = KGDH_Vmax * m.alpha_KG / (KGDH_Km_aKG + m.alpha_KG)
                            * m.CoA_free / (KGDH_Km_CoA + m.CoA_free)
                            * m.NAD_plus / (KGDH_Km_NAD + m.NAD_plus);
        v *= KGDH_Ki_succinylCoA / (KGDH_Ki_succinylCoA + m.succinyl_CoA);
        v *= KGDH_Ki_NADH / (KGDH_Ki_NADH + m.NADH);
        return v;
    }

    // ── Step 5: Succinyl-CoA Synthetase (EC 6.2.1.4) ──────────────
    // Succinyl-CoA + GDP + Pi → Succinate + CoA + GTP
    // Substrate-level phosphorylation (GTP!)
    static constexpr float SCS_Vmax = 1.0f;
    static constexpr float SCS_Km = 0.02f;

    static float succinylCoASynthetase(const MetabolitePool& m) {
        return michaelisMenten(m.succinyl_CoA, SCS_Vmax, SCS_Km)
             * m.GDP / (0.05f + m.GDP);
    }

    // ── Step 6: Succinate Dehydrogenase (Complex II) (EC 1.3.5.1) ──
    // Succinate + FAD → Fumarate + FADH2
    // Only TCA enzyme in inner mitochondrial membrane
    // Also Complex II of electron transport chain
    static constexpr float SDH_Vmax = 0.6f;
    static constexpr float SDH_Km = 0.3f;

    static float succinateDehydrogenase(const MetabolitePool& m) {
        return michaelisMenten(m.succinate, SDH_Vmax, SDH_Km)
             * m.FAD / (0.05f + m.FAD);
    }

    // ── Step 7: Fumarase (EC 4.2.1.2) ─────────────────────────────
    // Fumarate + H2O ⇌ Malate
    static constexpr float FUM_Vmax = 8.0f;
    static constexpr float FUM_Km = 0.05f;

    static float fumarase(const MetabolitePool& m) {
        return reversibleMM(m.fumarate, m.malate, FUM_Vmax, FUM_Km, FUM_Vmax * 0.5f, FUM_Km * 2.0f);
    }

    // ── Step 8: Malate Dehydrogenase (MDH) (EC 1.1.1.37) ──────────
    // Malate + NAD+ ⇌ Oxaloacetate + NADH + H+
    // Equilibrium strongly favors malate — driven by removal of OAA by citrate synthase
    static constexpr float MDH_Vmax = 5.0f;
    static constexpr float MDH_Km_malate = 0.15f;
    static constexpr float MDH_Km_NAD = 0.1f;

    static float malateDehydrogenase(const MetabolitePool& m) {
        return MDH_Vmax * m.malate / (MDH_Km_malate + m.malate)
                        * m.NAD_plus / (MDH_Km_NAD + m.NAD_plus)
             - MDH_Vmax * 0.8f * m.oxaloacetate / (0.01f + m.oxaloacetate)
                                * m.NADH / (0.05f + m.NADH);
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  OXIDATIVE PHOSPHORYLATION
//  NADH/FADH2 → O2 via ETC → proton gradient → ATP synthase
//  Location: Inner mitochondrial membrane
//  Ref: Alberts Ch 14 "The Proton Pumps of the Electron-Transport Chain"
//       Alberts Ch 14 "ATP Production in Mitochondria"
// ══════════════════════════════════════════════════════════════════════════

struct OxidativePhosphorylation {
    // ── Complex I: NADH Dehydrogenase (EC 7.1.1.2) ─────────────────
    // NADH + H+ + Q → NAD+ + QH2
    // Pumps 4 H+ from matrix to IMS
    static constexpr float CI_Vmax = 0.8f;
    static constexpr float CI_Km_NADH = 0.02f;
    static constexpr int   CI_H_pumped = 4;

    static float complexI(const MetabolitePool& m) {
        return CI_Vmax * m.NADH / (CI_Km_NADH + m.NADH)
                       * m.O2 / (0.003f + m.O2); // O2 dependence
    }

    // ── Complex II: Succinate Dehydrogenase ─────────────────────────
    // (Handled in TCA cycle — SDH feeds electrons to Q pool)
    // FADH2 + Q → FAD + QH2
    // Does NOT pump protons
    static constexpr int CI2_H_pumped = 0;

    static float complexII(const MetabolitePool& m) {
        return michaelisMenten(m.FADH2, 0.5f, 0.01f);
    }

    // ── Complex III: Cytochrome bc1 (EC 7.1.1.8) ───────────────────
    // QH2 + 2 Cyt c(ox) → Q + 2 Cyt c(red) + 2 H+(pumped)
    // Q cycle: effectively pumps 4 H+ per 2 electrons
    static constexpr float CIII_Vmax = 1.5f;
    static constexpr int   CIII_H_pumped = 4; // per 2 electrons

    static float complexIII(float QH2_flux) {
        return QH2_flux; // Rate follows upstream supply
    }

    // ── Complex IV: Cytochrome c Oxidase (EC 7.1.1.9) ──────────────
    // 4 Cyt c(red) + O2 + 8 H+(matrix) → 4 Cyt c(ox) + 2 H2O + 4 H+(IMS)
    // Pumps 2 H+ per electron (total 4 per O2)
    static constexpr float CIV_Vmax = 1.0f;
    static constexpr float CIV_Km_O2 = 0.001f; // Very high affinity for O2
    static constexpr int   CIV_H_pumped = 2;    // per electron

    static float complexIV(float cytC_flux, const MetabolitePool& m) {
        return cytC_flux * m.O2 / (CIV_Km_O2 + m.O2);
    }

    // ── ATP Synthase (Complex V) (EC 7.1.2.2) ──────────────────────
    // ADP + Pi + ~3.3 H+(IMS→matrix) → ATP + H2O
    // Rotary motor: F0 (proton channel) drives F1 (catalytic) rotation
    // 360° rotation = 3 ATP produced
    // ~10 c-subunits in human F0 → 10/3 ≈ 3.3 H+ per ATP
    static constexpr float ATPsyn_Vmax = 2.0f;
    static constexpr float ATPsyn_Km_ADP = 0.05f;
    static constexpr float ATPsyn_Km_Pi = 1.0f;
    static constexpr float H_per_ATP = 3.33f;  // 10 c-subunits / 3 catalytic sites

    static float atpSynthase(const MetabolitePool& m) {
        // Rate depends on ADP, Pi availability AND proton gradient
        float substrate = m.ADP / (ATPsyn_Km_ADP + m.ADP)
                        * m.Pi / (ATPsyn_Km_Pi + m.Pi);
        // Thermodynamic driving force from PMF
        float pmf = m.protonMotiveForce();
        float pmf_drive = fmaxf(0, pmf - 0.12f) / 0.08f; // Threshold ~120mV
        pmf_drive = fminf(pmf_drive, 1.0f);
        return ATPsyn_Vmax * substrate * pmf_drive;
    }

    // ── Total oxidative phosphorylation ─────────────────────────────
    // Per NADH: ~2.5 ATP (10 H+ pumped / ~4 H+ per ATP)
    // Per FADH2: ~1.5 ATP (6 H+ pumped / ~4 H+ per ATP)
    static constexpr float ATP_per_NADH  = 2.5f;
    static constexpr float ATP_per_FADH2 = 1.5f;
};

// ══════════════════════════════════════════════════════════════════════════
//  FERMENTATION
//  Regenerates NAD+ when O2 is unavailable
//  Ref: Alberts Ch 2 "Metabolism in the Absence of Oxygen: Fermentation"
// ══════════════════════════════════════════════════════════════════════════

struct Fermentation {
    // ── Lactate Dehydrogenase (LDH) (EC 1.1.1.27) ──────────────────
    // Pyruvate + NADH + H+ → Lactate + NAD+
    // Predominant in mammalian muscle under anaerobic conditions
    static constexpr float LDH_Vmax = 3.0f;
    static constexpr float LDH_Km_pyruvate = 0.1f;
    static constexpr float LDH_Km_NADH = 0.02f;

    static float lacticAcidFermentation(const MetabolitePool& m) {
        // Only significant when O2 is low
        float o2_suppression = 0.003f / (0.003f + m.O2); // Active when O2 < 3µM
        return LDH_Vmax * m.pyruvate / (LDH_Km_pyruvate + m.pyruvate)
                        * m.NADH / (LDH_Km_NADH + m.NADH)
                        * o2_suppression;
    }

    // ── Alcohol Dehydrogenase (ADH) (EC 1.1.1.1) ───────────────────
    // Pyruvate → Acetaldehyde + CO2 → Ethanol + NAD+
    // Yeast pathway (not normally active in human cells)
    static float alcoholicFermentation(const MetabolitePool& m) {
        // Not modeled for human cells — return 0
        return 0.0f;
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  β-OXIDATION OF FATTY ACIDS
//  Palmitoyl-CoA → 8 Acetyl-CoA + 7 NADH + 7 FADH2
//  Location: Mitochondrial matrix
//  Ref: Alberts Ch 2 "Fuel Choice"
// ══════════════════════════════════════════════════════════════════════════

struct FattyAcidOxidation {
    static constexpr float BETA_OX_Vmax = 0.1f;
    static constexpr float BETA_OX_Km = 0.005f;

    // Each cycle: removes 2 carbons as acetyl-CoA, produces 1 NADH + 1 FADH2
    // Palmitate (C16): 7 cycles → 8 acetyl-CoA + 7 NADH + 7 FADH2
    static float rate(const MetabolitePool& m) {
        return michaelisMenten(m.palmitoyl_CoA, BETA_OX_Vmax, BETA_OX_Km)
             * m.NAD_plus / (0.1f + m.NAD_plus)
             * m.FAD / (0.05f + m.FAD)
             * m.CoA_free / (0.01f + m.CoA_free);
    }

    static constexpr float NADH_per_cycle = 1.0f;
    static constexpr float FADH2_per_cycle = 1.0f;
    static constexpr float acetylCoA_per_cycle = 1.0f;
    static constexpr int   cycles_palmitate = 7;
};

// ══════════════════════════════════════════════════════════════════════════
//  PENTOSE PHOSPHATE PATHWAY (Oxidative Branch)
//  Glucose-6-P → Ribose-5-P + 2 NADPH + CO2
//  Provides NADPH for biosynthesis and ribose for nucleotides
//  Ref: Alberts Ch 2 "Activated Carrier Molecules" (NADPH section)
// ══════════════════════════════════════════════════════════════════════════

struct PentosePhosphatePathway {
    // Glucose-6-Phosphate Dehydrogenase (G6PDH) (EC 1.1.1.49)
    // Rate-limiting step, regulated by NADP+/NADPH ratio
    static constexpr float G6PDH_Vmax = 0.2f;
    static constexpr float G6PDH_Km_G6P = 0.05f;
    static constexpr float G6PDH_Km_NADP = 0.003f;

    static float rate(const MetabolitePool& m) {
        // Strongly regulated by NADP+ availability
        // When NADPH is consumed (by biosynthesis), NADP+ rises, driving this pathway
        return G6PDH_Vmax * m.glucose_6P / (G6PDH_Km_G6P + m.glucose_6P)
                          * m.NADP_plus / (G6PDH_Km_NADP + m.NADP_plus);
    }

    static constexpr float NADPH_per_G6P = 2.0f;
    static constexpr float ribose5P_per_G6P = 1.0f;
};

// ══════════════════════════════════════════════════════════════════════════
//  ADENYLATE KINASE (Myokinase)
//  2 ADP ⇌ ATP + AMP
//  Maintains adenylate equilibrium
// ══════════════════════════════════════════════════════════════════════════

struct AdenylateKinase {
    static constexpr float AK_Vmax = 50.0f;  // Very fast, near equilibrium
    static constexpr float AK_Keq = 0.44f;   // [ATP][AMP] / [ADP]^2

    static float rate(const MetabolitePool& m) {
        // Net direction depends on mass action ratio vs Keq
        float Q = (m.ATP * m.AMP) / fmaxf(m.ADP * m.ADP, 1e-9f);
        return AK_Vmax * (1.0f - Q / AK_Keq) * m.ADP / (0.1f + m.ADP);
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  MAIN METABOLISM ENGINE
//  Integrates all pathways and updates metabolite pools
// ══════════════════════════════════════════════════════════════════════════

class MetabolismEngine {
public:
    MetabolitePool pool;
    bool initialized = false;

    // Flux tracking (mM/s for each reaction)
    struct MetabolicFlux {
        // Glycolysis
        float hexokinase = 0, PGI = 0, PFK1 = 0, aldolase = 0, TPI = 0;
        float GAPDH = 0, PGK = 0, PGM = 0, enolase = 0, PK = 0;
        // Pyruvate handling
        float PDH = 0, LDH = 0;
        // TCA cycle
        float citrateSynthase = 0, aconitase = 0, IDH = 0, alphaKGDH = 0;
        float SCS = 0, SDH = 0, fumarase = 0, MDH = 0;
        // OxPhos
        float complexI = 0, complexII = 0, complexIII = 0, complexIV = 0;
        float atpSynthase = 0;
        // Other
        float betaOxidation = 0, PPP = 0, adenylateKinase = 0;
        // Summary
        float totalATPproduction = 0, totalATPconsumption = 0;
        float totalO2consumption = 0;
    };
    MetabolicFlux flux;

    // ATP consumption rate (from cellular processes: protein synthesis, ion pumps, etc.)
    float atpConsumptionRate = 0.5f; // mM/s — baseline cellular ATP demand
    float external_O2 = 0.03f;       // mM — set by caller from nutrient field
    void setExternalO2(float o2_mM) { external_O2 = o2_mM; }

    void init() {
        pool = MetabolitePool(); // Reset to defaults
        initialized = true;
    }

    // ── Step the metabolism by dt seconds ────────────────────────────
    void step(float dt) {
        if (!initialized) init();
        dt = fminf(dt, 0.01f); // Max 10ms step for stability

        // ── Compute all fluxes ──────────────────────────────────────
        // Glycolysis
        flux.hexokinase = GlycolysisEnzymes::hexokinase(pool);
        flux.PGI = GlycolysisEnzymes::phosphoglucoseIsomerase(pool);
        flux.PFK1 = GlycolysisEnzymes::phosphofructokinase1(pool);
        flux.aldolase = GlycolysisEnzymes::aldolase(pool);
        flux.TPI = GlycolysisEnzymes::triosephosphateIsomerase(pool);
        flux.GAPDH = GlycolysisEnzymes::G3P_dehydrogenase(pool);
        flux.PGK = GlycolysisEnzymes::phosphoglycerateKinase(pool);
        flux.PGM = GlycolysisEnzymes::phosphoglycerateMutase(pool);
        flux.enolase = GlycolysisEnzymes::enolase(pool);
        flux.PK = GlycolysisEnzymes::pyruvateKinase(pool);

        // Pyruvate fate
        flux.PDH = PyruvateDehydrogenaseComplex::rate(pool);
        flux.LDH = Fermentation::lacticAcidFermentation(pool);

        // TCA cycle
        flux.citrateSynthase = CitricAcidCycleEnzymes::citrateSynthase(pool);
        flux.aconitase = CitricAcidCycleEnzymes::aconitase(pool);
        flux.IDH = CitricAcidCycleEnzymes::isocitrateDehydrogenase(pool);
        flux.alphaKGDH = CitricAcidCycleEnzymes::alphaKGDehydrogenase(pool);
        flux.SCS = CitricAcidCycleEnzymes::succinylCoASynthetase(pool);
        flux.SDH = CitricAcidCycleEnzymes::succinateDehydrogenase(pool);
        flux.fumarase = CitricAcidCycleEnzymes::fumarase(pool);
        flux.MDH = CitricAcidCycleEnzymes::malateDehydrogenase(pool);

        // Oxidative phosphorylation
        flux.complexI = OxidativePhosphorylation::complexI(pool);
        flux.complexII = OxidativePhosphorylation::complexII(pool);
        float totalQH2 = flux.complexI + flux.complexII;
        flux.complexIII = OxidativePhosphorylation::complexIII(totalQH2);
        flux.complexIV = OxidativePhosphorylation::complexIV(flux.complexIII, pool);
        flux.atpSynthase = OxidativePhosphorylation::atpSynthase(pool);

        // Other pathways
        flux.betaOxidation = FattyAcidOxidation::rate(pool);
        flux.PPP = PentosePhosphatePathway::rate(pool);
        flux.adenylateKinase = AdenylateKinase::rate(pool);

        // ── Update metabolite pools ─────────────────────────────────
        // Glycolysis intermediates
        pool.glucose     -= flux.hexokinase * dt;
        pool.glucose_6P  += (flux.hexokinase - flux.PGI - flux.PPP) * dt;
        pool.fructose_6P += (flux.PGI - flux.PFK1) * dt;
        pool.fructose_16BP += (flux.PFK1 - flux.aldolase) * dt;
        pool.DHAP        += (flux.aldolase - flux.TPI) * dt;
        pool.G3P         += (flux.aldolase + flux.TPI - flux.GAPDH) * dt;
        pool.BPG_13      += (flux.GAPDH - flux.PGK) * dt;
        pool.PG_3        += (flux.PGK - flux.PGM) * dt;
        pool.PG_2        += (flux.PGM - flux.enolase) * dt;
        pool.PEP         += (flux.enolase - flux.PK) * dt;
        pool.pyruvate    += (flux.PK - flux.PDH - flux.LDH) * dt;

        // TCA intermediates
        pool.acetyl_CoA  += (flux.PDH - flux.citrateSynthase + flux.betaOxidation) * dt;
        pool.citrate     += (flux.citrateSynthase - flux.aconitase) * dt;
        pool.isocitrate  += (flux.aconitase - flux.IDH) * dt;
        pool.alpha_KG    += (flux.IDH - flux.alphaKGDH) * dt;
        pool.succinyl_CoA += (flux.alphaKGDH - flux.SCS) * dt;
        pool.succinate   += (flux.SCS - flux.SDH) * dt;
        pool.fumarate    += (flux.SDH - flux.fumarase) * dt;
        pool.malate      += (flux.fumarase - flux.MDH) * dt;
        pool.oxaloacetate += flux.MDH * dt; // Consumed by citrate synthase

        // CoA recycling
        pool.CoA_free    += (flux.citrateSynthase + flux.SCS - flux.PDH - flux.alphaKGDH) * dt;

        // ── Energy carriers ────────────────────────────────────────
        // ATP bookkeeping. Rate laws above referenced ATP/ADP as substrates
        // for their kinetics (Michaelis-Menten terms) but did NOT modify the
        // pool. Pool update is therefore applied only ONCE here from the
        // stoichiometry of each reaction.
        //
        // Glycolysis: HK consumes 1 ATP (glucose level), PFK consumes 1 ATP.
        // PGK and PK each produce 1 ATP, but operate on the triose pool
        // (G3P/PEP), which is already fed by BOTH trioses from aldolase.
        // No ×2 multiplier — double-counting was a bug.
        float atp_produced = flux.PGK + flux.PK + flux.SCS + flux.atpSynthase;
        float atp_consumed = flux.hexokinase + flux.PFK1 + atpConsumptionRate;

        flux.totalATPproduction  = atp_produced;
        flux.totalATPconsumption = atp_consumed;

        float d_atp = (atp_produced - atp_consumed) * dt;
        if (d_atp >= 0.0f) {
            // ATP synthesis is ADP-limited. If ADP hits zero, clamping later
            // would otherwise let ATP keep increasing and violate conservation.
            float synth = fminf(d_atp, pool.ADP);
            pool.ATP += synth;
            pool.ADP -= synth;
        } else {
            // ATP hydrolysis is limited by the ATP that is actually present.
            float hydro = fminf(-d_atp, pool.ATP);
            pool.ATP -= hydro;
            pool.ADP += hydro;
        }

        // NADH/NAD+ (cytosolic + mitochondrial combined for simplicity)
        float nadh_produced = flux.GAPDH * 2.0f + flux.PDH + flux.IDH + flux.alphaKGDH + flux.MDH;
        float nadh_consumed = flux.complexI + flux.LDH;
        pool.NADH += (nadh_produced - nadh_consumed) * dt;
        pool.NAD_plus += (nadh_consumed - nadh_produced) * dt;

        // FADH2/FAD
        float fadh2_produced = flux.SDH + flux.betaOxidation;
        float fadh2_consumed = flux.complexII;
        pool.FADH2 += (fadh2_produced - fadh2_consumed) * dt;
        pool.FAD += (fadh2_consumed - fadh2_produced) * dt;

        // NADPH (pentose phosphate pathway)
        pool.NADPH += flux.PPP * PentosePhosphatePathway::NADPH_per_G6P * dt;
        pool.NADP_plus -= flux.PPP * PentosePhosphatePathway::NADPH_per_G6P * dt;
        // NADPH consumption (biosynthesis — simplified)
        float nadph_consumption = 0.01f; // Baseline biosynthetic demand
        pool.NADPH -= nadph_consumption * dt;
        pool.NADP_plus += nadph_consumption * dt;

        // GTP/GDP — produced by SCS, consumed by background cellular demand
        // (translation, Ras/Rho/Rab GTP hydrolysis). Without consumption GTP
        // would rise without bound.
        float gtp_consumption = 0.05f * pool.GTP / (0.3f + pool.GTP);
        pool.GTP += (flux.SCS - gtp_consumption) * dt;
        pool.GDP += (gtp_consumption - flux.SCS) * dt;

        // Ribose-5-phosphate
        pool.ribose_5P += flux.PPP * PentosePhosphatePathway::ribose5P_per_G6P * dt;

        // ── Proton gradient & membrane potential ─────────────────
        // Treat ΔΨ as fast-equilibrium with net proton flux.
        // Protons pumped by ETC complexes:
        float h_pumped = flux.complexI  * OxidativePhosphorylation::CI_H_pumped
                       + flux.complexIII * OxidativePhosphorylation::CIII_H_pumped
                       + flux.complexIV  * OxidativePhosphorylation::CIV_H_pumped;
        // Protons consumed by ATP synthase:
        float h_consumed = flux.atpSynthase * OxidativePhosphorylation::H_per_ATP;
        // Proton leak (background dissipation across inner membrane):
        float h_leak = 0.2f * fmaxf(0.0f, pool.membrane_potential - 0.14f);
        float j_net = h_pumped - h_consumed - h_leak; // mM/s

        // [H+] in matrix and IMS — keep concentrations realistic by clamping
        // the pH swing rather than letting raw j_net drive them (they're
        // strongly buffered by phosphate/bicarbonate in reality).
        float h_scale = 1e-6f; // couples flux to Δ[H+] in mM
        pool.H_plus_IMS    = clampPositive(pool.H_plus_IMS    + j_net * h_scale * dt);
        pool.H_plus_matrix = clampPositive(pool.H_plus_matrix - j_net * h_scale * dt);

        // ΔΨ as a first-order relaxation driven by net current.
        // Time constant ~1 s, steady-state ΔΨ ≈ 180 mV when j_net > 0.
        float dpsi_target = 0.180f + 0.02f * tanhf(j_net * 2.0f);
        pool.membrane_potential += (dpsi_target - pool.membrane_potential) * dt;
        pool.membrane_potential = std::clamp(pool.membrane_potential, 0.10f, 0.22f);

        // ── O2 consumption ─────────────────────────────────────────
        flux.totalO2consumption = flux.complexIV * 0.25f; // 1 O2 per 4 electrons
        pool.O2 -= flux.totalO2consumption * dt;

        // ── CO2 production ─────────────────────────────────────────
        pool.CO2 += (flux.PDH + flux.IDH + flux.alphaKGDH) * dt;

        // ── Lactate ────────────────────────────────────────────────
        pool.lactate += flux.LDH * dt;

        // ── Adenylate kinase equilibrium ────────────────────────────
        float ak = flux.adenylateKinase * dt;
        if (ak >= 0.0f) {
            ak = fminf(ak, pool.ADP * 0.5f);
        } else {
            ak = -fminf(-ak, fminf(pool.ATP, pool.AMP));
        }
        pool.ADP -= 2.0f * ak;
        pool.ATP += ak;
        pool.AMP += ak;

        // ── Glucose import (from extracellular) ────────────────────
        // GLUT transporter — facilitated diffusion
        float glucose_import = 0.5f * (5.0f - pool.glucose) / (1.0f + pool.glucose);
        pool.glucose += glucose_import * dt;

        // ── O2 resupply — coupled to external field ────────────────
        // P·A for a 20 µm cell is on the order of 0.01 /s; use that as a
        // first-order relaxation toward externally provided O2.
        // external_O2 is set by the caller via setExternalO2(); defaults
        // to 0.03 mM if not set.
        float o2_resupply = 0.1f * (external_O2 - pool.O2);
        pool.O2 += o2_resupply * dt;

        // ── Clamp all concentrations to prevent negative values ─────
        pool.glucose      = clampPositive(pool.glucose);
        pool.glucose_6P   = clampPositive(pool.glucose_6P);
        pool.fructose_6P  = clampPositive(pool.fructose_6P);
        pool.fructose_16BP = clampPositive(pool.fructose_16BP);
        pool.DHAP         = clampPositive(pool.DHAP);
        pool.G3P          = clampPositive(pool.G3P);
        pool.BPG_13       = clampPositive(pool.BPG_13);
        pool.PG_3         = clampPositive(pool.PG_3);
        pool.PG_2         = clampPositive(pool.PG_2);
        pool.PEP          = clampPositive(pool.PEP);
        pool.pyruvate     = clampPositive(pool.pyruvate);
        pool.acetyl_CoA   = clampPositive(pool.acetyl_CoA);
        pool.oxaloacetate = clampPositive(pool.oxaloacetate);
        pool.citrate      = clampPositive(pool.citrate);
        pool.isocitrate   = clampPositive(pool.isocitrate);
        pool.alpha_KG     = clampPositive(pool.alpha_KG);
        pool.succinyl_CoA = clampPositive(pool.succinyl_CoA);
        pool.succinate    = clampPositive(pool.succinate);
        pool.fumarate     = clampPositive(pool.fumarate);
        pool.malate       = clampPositive(pool.malate);
        pool.CoA_free     = clampPositive(pool.CoA_free);
        pool.ATP          = clampPositive(pool.ATP);
        pool.ADP          = clampPositive(pool.ADP);
        pool.AMP          = clampPositive(pool.AMP);
        pool.GTP          = clampPositive(pool.GTP);
        pool.GDP          = clampPositive(pool.GDP);
        pool.NAD_plus     = clampPositive(pool.NAD_plus);
        pool.NADH         = clampPositive(pool.NADH);
        pool.NADP_plus    = clampPositive(pool.NADP_plus);
        pool.NADPH        = clampPositive(pool.NADPH);
        pool.FAD          = clampPositive(pool.FAD);
        pool.FADH2        = clampPositive(pool.FADH2);
        pool.Pi           = clampPositive(pool.Pi);
        pool.O2           = clampPositive(pool.O2);
        pool.CO2          = clampPositive(pool.CO2);
        pool.lactate      = clampPositive(pool.lactate);
        pool.ribose_5P    = clampPositive(pool.ribose_5P);
        pool.H_plus_matrix = clampPositive(pool.H_plus_matrix);
        pool.H_plus_IMS   = clampPositive(pool.H_plus_IMS);

        // Mass-balance diagnostics only — do NOT rescale.
        // If these drift, it signals a bug in the flux bookkeeping and
        // should be fixed at the source, not masked here.
        float total_aden = pool.ATP + pool.ADP + pool.AMP;
        float total_nad  = pool.NAD_plus + pool.NADH;
        flux_conservation_adenylate = total_aden;
        flux_conservation_nad = total_nad;
    }

    // Conservation diagnostics (should stay near 3.9 and 1.6 if balanced)
    float flux_conservation_adenylate = 0;
    float flux_conservation_nad = 0;

    // ── Get ATP level as percentage (for compatibility with existing Simulation.h) ──
    float getATPpercent() const {
        return (pool.ATP / 3.9f) * 100.0f; // 3.9 mM is ~100%
    }

    // ── Print diagnostic ────────────────────────────────────────────
    void printStatus() const {
        printf("[Metabolism] ATP=%.2f ADP=%.2f AMP=%.2f (EC=%.3f)\n",
               pool.ATP, pool.ADP, pool.AMP, pool.energyCharge());
        printf("  NADH/NAD+=%.3f  FADH2=%.3f  O2=%.4f  Lactate=%.2f\n",
               pool.NAD_ratio(), pool.FADH2, pool.O2, pool.lactate);
        printf("  PMF=%.1fmV  ATP synthase=%.4f mM/s\n",
               pool.protonMotiveForce() * 1000, flux.atpSynthase);
        printf("  ATP prod=%.4f cons=%.4f mM/s  O2 cons=%.5f\n",
               flux.totalATPproduction, flux.totalATPconsumption, flux.totalO2consumption);
    }
};
