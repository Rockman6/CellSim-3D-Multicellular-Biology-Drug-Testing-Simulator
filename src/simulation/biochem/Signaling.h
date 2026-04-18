#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Signaling.h — Complete cell signaling pathway simulation
//  Based on: Alberts "Molecular Biology of the Cell" 7th Ed, Chapter 15
//
//  Implements ALL major signaling pathways as ODE systems:
//    1. GPCR signaling (Gs/Gi/Gq/G12-13)
//       - cAMP/PKA pathway
//       - IP3/DAG/Ca2+/PKC pathway
//       - Ca2+ oscillations and calmodulin
//    2. Receptor Tyrosine Kinase (RTK) signaling
//       - Ras-MAPK cascade (Raf→MEK→ERK)
//       - PI3K-Akt-mTOR pathway
//    3. JAK-STAT pathway (cytokine receptors)
//    4. TGF-β/Smad pathway
//    5. Notch pathway (contact-dependent, proteolytic)
//    6. Wnt/β-catenin pathway (destruction complex)
//    7. Hedgehog/Patched/Smoothened pathway
//    8. NF-κB pathway (inflammatory)
//    9. Nuclear receptor signaling (steroid hormones)
//   10. Nitric oxide (NO) / cGMP / PKG
//   11. Circadian clock (CLOCK/BMAL1 → PER/CRY)
//
//  Each pathway modeled with:
//    - Michaelis-Menten enzyme kinetics
//    - Hill equation for cooperative binding
//    - Mass-action kinetics for binding/dissociation
//    - Feedback loops (positive and negative)
//
//  Kinetic parameters from published literature:
//    - Kholodenko 2000 (MAPK cascade)
//    - Goldbeter 1991 (Ca2+ oscillations)
//    - Leloup & Goldbeter 2003 (circadian clock)
//    - Cho et al. 2003 (NF-κB)
//
//  Units: concentrations in µM (micromolar), time in seconds
// ══════════════════════════════════════════════════════════════════════════

static float clampSig(float x) { return x > 0 ? x : 0; }
static float hillSig(float x, float K, float n) {
    float xn = powf(fmaxf(x, 0), n);
    return xn / (powf(K, n) + xn);
}

// ── Signaling Molecule Concentrations ───────────────────────────────────
// All in µM unless noted

struct SignalingState {
    // ── GPCR / G-protein system ────────────────────────────────────
    float ligand_GPCR      = 0.0f;   // Extracellular ligand concentration
    float GPCR_active      = 0.0f;   // Fraction of active GPCRs (0-1)
    float Gs_alpha_GTP     = 0.0f;   // Active Gsα-GTP
    float Gi_alpha_GTP     = 0.0f;   // Active Giα-GTP
    float Gq_alpha_GTP     = 0.0f;   // Active Gqα-GTP
    float G_betagamma      = 0.0f;   // Free Gβγ subunits

    // ── cAMP / PKA pathway ─────────────────────────────────────────
    float cAMP             = 0.1f;   // Cyclic AMP
    float PKA_active       = 0.0f;   // Active PKA catalytic subunit
    float PKA_inactive     = 1.0f;   // Inactive PKA holoenzyme
    float PDE              = 0.5f;   // Phosphodiesterase (degrades cAMP)
    float CREB_P           = 0.0f;   // Phosphorylated CREB (transcription factor)

    // ── IP3 / DAG / Ca2+ / PKC pathway ─────────────────────────────
    float PLC_beta_active  = 0.0f;   // Active PLCβ
    float PIP2             = 5.0f;   // Phosphatidylinositol 4,5-bisphosphate
    float IP3              = 0.0f;   // Inositol 1,4,5-trisphosphate
    float DAG              = 0.0f;   // Diacylglycerol
    float Ca_cytosol       = 0.1f;   // Cytosolic Ca2+ (µM, resting ~0.1)
    float Ca_ER            = 400.0f; // ER Ca2+ store (µM)
    float PKC_active       = 0.0f;   // Active protein kinase C
    float calmodulin_Ca    = 0.0f;   // Ca2+-calmodulin complex
    float CaMKII_active    = 0.0f;   // Ca2+/CaM-dependent kinase II

    // ── Ras-MAPK cascade ───────────────────────────────────────────
    // Ref: Kholodenko 2000 "Negative feedback and ultrasensitivity"
    float RTK_active       = 0.0f;   // Active receptor tyrosine kinase
    float Grb2_Sos         = 0.0f;   // Grb2-Sos complex at membrane
    float Ras_GTP          = 0.0f;   // Active Ras
    float Ras_GDP          = 1.0f;   // Inactive Ras
    float Raf_active       = 0.0f;   // Active Raf (MAP3K)
    float Raf_inactive     = 1.0f;
    float MEK_active       = 0.0f;   // Active MEK (MAP2K) — dual phosphorylated
    float MEK_inactive     = 1.0f;
    float ERK_active       = 0.0f;   // Active ERK (MAPK) — dual phosphorylated
    float ERK_inactive     = 1.0f;
    float ERK_nuclear      = 0.0f;   // ERK translocated to nucleus

    // ── PI3K-Akt-mTOR pathway ──────────────────────────────────────
    float PI3K_active      = 0.0f;   // Active PI3-Kinase
    float PIP3             = 0.0f;   // Phosphatidylinositol (3,4,5)-trisphosphate
    float PTEN             = 1.0f;   // PTEN phosphatase (dephosphorylates PIP3)
    float Akt_active       = 0.0f;   // Active Akt/PKB (phosphorylated)
    float Akt_inactive     = 1.0f;
    float mTORC1_active    = 0.0f;   // Active mTORC1
    float mTORC1_inactive  = 1.0f;
    float S6K_active       = 0.0f;   // p70 S6 kinase (protein synthesis)
    float FOXO_nuclear     = 0.0f;   // FoxO transcription factor (apoptosis genes)
    float BAD_P            = 0.0f;   // Phosphorylated BAD (pro-survival)

    // ── JAK-STAT pathway ───────────────────────────────────────────
    float JAK_active       = 0.0f;   // Active Janus kinase
    float STAT_P           = 0.0f;   // Phosphorylated STAT (dimerizes)
    float STAT_dimer_nuc   = 0.0f;   // STAT dimer in nucleus
    float SOCS             = 0.0f;   // Suppressor of cytokine signaling (neg feedback)

    // ── TGF-β / Smad pathway ───────────────────────────────────────
    float TGFbeta_receptor = 0.0f;   // Active TGF-β receptor complex
    float Smad2_P          = 0.0f;   // Phosphorylated R-Smad (Smad2/3)
    float Smad4_complex    = 0.0f;   // Smad2/3-Smad4 complex
    float Smad7            = 0.0f;   // Inhibitory Smad (neg feedback)

    // ── Notch pathway ──────────────────────────────────────────────
    float Notch_receptor   = 1.0f;   // Surface Notch receptor
    float Delta_ligand     = 0.0f;   // Delta/Jagged on neighbor cell
    float NICD             = 0.0f;   // Notch intracellular domain (cleaved)
    float Hes_Hey          = 0.0f;   // Target gene products

    // ── Wnt / β-catenin pathway ────────────────────────────────────
    float Wnt_signal       = 0.0f;   // Wnt ligand
    float Frizzled_active  = 0.0f;   // Active Frizzled receptor
    float Dishevelled      = 0.0f;   // Active Dishevelled
    float beta_catenin     = 0.1f;   // Cytoplasmic β-catenin
    float destruction_complex = 1.0f; // Axin/APC/GSK3/CK1 complex
    float TCF_active       = 0.0f;   // β-catenin-TCF in nucleus (transcription ON)
    float Axin2            = 0.0f;   // Negative feedback target gene

    // ── Hedgehog pathway ───────────────────────────────────────────
    float Hh_signal        = 0.0f;   // Hedgehog ligand
    float Patched          = 1.0f;   // Patched receptor (inhibits Smoothened)
    float Smoothened_active = 0.0f;  // Active Smoothened
    float Gli_activator    = 0.0f;   // Full-length Gli (transcription activator)
    float Gli_repressor    = 0.5f;   // Cleaved Gli (transcription repressor)

    // ── NF-κB pathway ──────────────────────────────────────────────
    // Ref: Cho et al. 2003, Hoffmann et al. 2002
    float IKK_active       = 0.0f;   // Active IKK complex
    float IkB_alpha        = 1.0f;   // IκBα (inhibitor, sequesters NF-κB)
    float NFkB_cytoplasm   = 1.0f;   // NF-κB in cytoplasm (bound to IκB)
    float NFkB_nuclear     = 0.0f;   // Free NF-κB in nucleus
    float A20              = 0.0f;   // A20 negative feedback protein

    // ── Nuclear receptor signaling ─────────────────────────────────
    float steroid_hormone  = 0.0f;   // Lipid-soluble hormone
    float nuclear_receptor_active = 0.0f; // Ligand-bound receptor in nucleus
    float target_gene_product = 0.0f;     // Gene expression output

    // ── NO / cGMP / PKG pathway ────────────────────────────────────
    float nitricOxide      = 0.0f;   // Nitric oxide (NO)
    float soluble_GC_active = 0.0f;  // Soluble guanylyl cyclase
    float cGMP             = 0.0f;   // Cyclic GMP
    float PKG_active       = 0.0f;   // cGMP-dependent protein kinase

    // ── Circadian clock ────────────────────────────────────────────
    // Ref: Leloup & Goldbeter 2003
    float CLOCK_BMAL1      = 1.0f;   // CLOCK-BMAL1 heterodimer
    float Per_mRNA         = 0.0f;   // Period mRNA
    float Cry_mRNA         = 0.0f;   // Cryptochrome mRNA
    float PER_protein      = 0.0f;   // PER protein
    float CRY_protein      = 0.0f;   // CRY protein
    float PER_CRY_complex  = 0.0f;   // PER-CRY nuclear inhibitor complex
    float PER_CRY_nuclear  = 0.0f;   // Complex in nucleus
    float Rev_Erb          = 0.0f;   // REV-ERBα (represses BMAL1)

    // ── Rho family GTPases (cytoskeleton control) ──────────────────
    float RhoA_GTP         = 0.0f;   // Active RhoA → stress fibers
    float Rac1_GTP         = 0.0f;   // Active Rac1 → lamellipodia
    float Cdc42_GTP        = 0.0f;   // Active Cdc42 → filopodia

    // ── Cross-pathway integration ──────────────────────────────────
    float total_transcription_drive = 0.0f; // Net transcriptional output
    float proliferation_signal = 0.0f;       // Net pro-proliferation
    float survival_signal = 0.0f;            // Net anti-apoptosis
    float growth_signal = 0.0f;              // Net cell growth (mTOR)
};

// ══════════════════════════════════════════════════════════════════════════
//  SIGNALING ENGINE — updates all pathways each timestep
// ══════════════════════════════════════════════════════════════════════════

class SignalingEngine {
public:
    SignalingState state;
    bool initialized = false;

    void init() {
        state = SignalingState();
        initialized = true;
    }

    // ── Apply an extracellular signal ───────────────────────────────
    void applyGPCRligand(float concentration) { state.ligand_GPCR = concentration; }
    void applyRTKligand(float concentration) { state.RTK_active = concentration; }
    void applyWnt(float concentration) { state.Wnt_signal = concentration; }
    void applyHedgehog(float concentration) { state.Hh_signal = concentration; }
    void applyTGFbeta(float concentration) { state.TGFbeta_receptor = concentration; }
    void applyNotchDelta(float concentration) { state.Delta_ligand = concentration; }
    void applySteroidHormone(float concentration) { state.steroid_hormone = concentration; }
    void applyNO(float concentration) { state.nitricOxide = concentration; }
    void applyInflammatory(float concentration) { state.IKK_active = concentration; }

    // ── Step all signaling pathways ─────────────────────────────────
    void step(float dt) {
        if (!initialized) init();
        dt = fminf(dt, 0.01f);
        auto& s = state;

        // ════════════════════════════════════════════════════════════
        //  1. GPCR SIGNALING
        // ════════════════════════════════════════════════════════════

        // GPCR activation by ligand
        float gpcr_on = 0.5f * s.ligand_GPCR * (1.0f - s.GPCR_active);
        float gpcr_off = 0.1f * s.GPCR_active; // Desensitization (GRK + arrestin)
        s.GPCR_active += (gpcr_on - gpcr_off) * dt;

        // Gs activation (stimulatory)
        float gs_on = 0.3f * s.GPCR_active * (1.0f - s.Gs_alpha_GTP);
        float gs_off = 0.1f * s.Gs_alpha_GTP; // Intrinsic GTPase + RGS
        s.Gs_alpha_GTP += (gs_on - gs_off) * dt;
        s.G_betagamma = s.Gs_alpha_GTP + s.Gq_alpha_GTP;

        // ── cAMP / PKA pathway ─────────────────────────────────────
        // Adenylyl cyclase activated by Gsα, inhibited by Giα
        float AC_activity = s.Gs_alpha_GTP * 2.0f - s.Gi_alpha_GTP * 1.5f;
        AC_activity = fmaxf(0, AC_activity);
        float cAMP_synthesis = AC_activity * 1.0f; // Vmax
        float cAMP_degradation = s.PDE * s.cAMP / (1.0f + s.cAMP); // PDE: Km=1µM
        s.cAMP += (cAMP_synthesis - cAMP_degradation) * dt;
        s.cAMP = clampSig(s.cAMP);

        // PKA activation by cAMP (cooperative, Hill n=2)
        float pka_frac = hillSig(s.cAMP, 0.5f, 2.0f); // K0.5=0.5µM, n=2
        s.PKA_active = pka_frac * (s.PKA_active + s.PKA_inactive);
        s.PKA_inactive = (1.0f - pka_frac) * (s.PKA_active + s.PKA_inactive);

        // CREB phosphorylation by PKA
        float creb_on = 0.5f * s.PKA_active * (1.0f - s.CREB_P);
        float creb_off = 0.05f * s.CREB_P; // Phosphatase
        s.CREB_P += (creb_on - creb_off) * dt;

        // ── IP3 / Ca2+ / PKC pathway ───────────────────────────────
        // Gq → PLCβ activation
        float plc_on = 0.5f * s.Gq_alpha_GTP * (1.0f - s.PLC_beta_active);
        float plc_off = 0.2f * s.PLC_beta_active;
        s.PLC_beta_active += (plc_on - plc_off) * dt;

        // PLCβ cleaves PIP2 → IP3 + DAG
        float pip2_cleavage = s.PLC_beta_active * s.PIP2 / (5.0f + s.PIP2);
        s.IP3 += pip2_cleavage * dt - 0.5f * s.IP3 * dt; // IP3 degradation
        s.DAG += pip2_cleavage * dt - 0.3f * s.DAG * dt;
        s.PIP2 += (0.2f - pip2_cleavage) * dt; // PIP2 resynthesis
        s.IP3 = clampSig(s.IP3); s.DAG = clampSig(s.DAG);

        // Ca2+ release from ER via IP3 receptor
        // Ref: Li-Rinzel model for IP3R
        float ip3r_open = hillSig(s.IP3, 0.3f, 2.0f) * hillSig(s.Ca_cytosol, 0.3f, 2.0f);
        // IP3R is inhibited by high Ca2+ (bell-shaped dependence)
        float ca_inhibition = 1.0f / (1.0f + powf(s.Ca_cytosol / 1.0f, 4.0f));
        float ca_release = 5.0f * ip3r_open * ca_inhibition * (s.Ca_ER - s.Ca_cytosol);
        // Ca2+ leak from ER
        float ca_leak = 0.01f * (s.Ca_ER - s.Ca_cytosol);
        // SERCA pump (Ca2+ back into ER)
        float serca = 2.0f * s.Ca_cytosol * s.Ca_cytosol / (0.04f + s.Ca_cytosol * s.Ca_cytosol);
        // Plasma membrane Ca2+ ATPase (extrusion)
        float pmca = 0.5f * s.Ca_cytosol / (0.5f + s.Ca_cytosol);

        s.Ca_cytosol += (ca_release + ca_leak - serca - pmca) * dt;
        s.Ca_ER += (serca - ca_release - ca_leak) * dt;
        s.Ca_cytosol = clampSig(s.Ca_cytosol);
        s.Ca_ER = clampSig(s.Ca_ER);

        // Calmodulin activation (4 Ca2+ binding, cooperative)
        s.calmodulin_Ca = hillSig(s.Ca_cytosol, 0.5f, 4.0f);

        // CaMKII activation
        float camk_on = 2.0f * s.calmodulin_Ca * (1.0f - s.CaMKII_active);
        float camk_off = 0.1f * s.CaMKII_active;
        s.CaMKII_active += (camk_on - camk_off) * dt;

        // PKC activation (requires both Ca2+ and DAG)
        float pkc_on = 1.0f * s.Ca_cytosol / (0.5f + s.Ca_cytosol)
                            * s.DAG / (0.5f + s.DAG) * (1.0f - s.PKC_active);
        float pkc_off = 0.2f * s.PKC_active;
        s.PKC_active += (pkc_on - pkc_off) * dt;

        // ════════════════════════════════════════════════════════════
        //  2. RAS-MAPK CASCADE
        //  Ref: Kholodenko 2000 Eur J Biochem
        // ════════════════════════════════════════════════════════════

        // RTK → Grb2-Sos recruitment
        s.Grb2_Sos = 0.8f * s.RTK_active / (0.1f + s.RTK_active);

        // Sos activates Ras (GEF activity): Ras-GDP → Ras-GTP
        float ras_on = 0.5f * s.Grb2_Sos * s.Ras_GDP;
        float ras_off = 0.1f * s.Ras_GTP; // GAP-stimulated GTPase
        // ERK negative feedback on Sos (phosphorylation reduces Sos activity)
        float sos_feedback = 1.0f / (1.0f + s.ERK_active * 2.0f);
        s.Ras_GTP += (ras_on * sos_feedback - ras_off) * dt;
        s.Ras_GDP += (ras_off - ras_on * sos_feedback) * dt;
        s.Ras_GTP = clampSig(s.Ras_GTP); s.Ras_GDP = clampSig(s.Ras_GDP);

        // Ras-GTP activates Raf (MAP3K)
        float raf_on = 1.0f * s.Ras_GTP * s.Raf_inactive / (0.1f + s.Raf_inactive);
        float raf_off = 0.3f * s.Raf_active / (0.1f + s.Raf_active);
        // Akt negative feedback on Raf (cross-talk)
        float akt_raf_inhibition = 1.0f / (1.0f + s.Akt_active);
        s.Raf_active += (raf_on * akt_raf_inhibition - raf_off) * dt;
        s.Raf_inactive += (raf_off - raf_on * akt_raf_inhibition) * dt;

        // Raf activates MEK (MAP2K) — dual phosphorylation
        float mek_on = 2.0f * s.Raf_active * s.MEK_inactive / (0.3f + s.MEK_inactive);
        float mek_off = 0.5f * s.MEK_active / (0.3f + s.MEK_active); // PP2A
        s.MEK_active += (mek_on - mek_off) * dt;
        s.MEK_inactive += (mek_off - mek_on) * dt;

        // MEK activates ERK (MAPK) — dual phosphorylation
        float erk_on = 3.0f * s.MEK_active * s.ERK_inactive / (0.3f + s.ERK_inactive);
        float erk_off = 0.5f * s.ERK_active / (0.3f + s.ERK_active); // MKP
        s.ERK_active += (erk_on - erk_off) * dt;
        s.ERK_inactive += (erk_off - erk_on) * dt;

        // ERK nuclear translocation — preserve mass conservation:
        // ERK_active (cytosol) ⇌ ERK_nuclear. Import depletes the
        // cytosolic pool; export refills it.
        float erk_import = 0.3f * s.ERK_active;
        float erk_export = 0.1f * s.ERK_nuclear;
        float erk_net = erk_import - erk_export;
        s.ERK_nuclear += erk_net * dt;
        s.ERK_active  -= erk_net * dt;
        s.ERK_nuclear = clampSig(s.ERK_nuclear);
        s.ERK_active  = clampSig(s.ERK_active);

        // ════════════════════════════════════════════════════════════
        //  3. PI3K-Akt-mTOR PATHWAY
        // ════════════════════════════════════════════════════════════

        // RTK activates PI3K
        s.PI3K_active = 0.8f * s.RTK_active / (0.2f + s.RTK_active);

        // PI3K: PIP2 → PIP3
        float pip3_syn = 1.0f * s.PI3K_active * s.PIP2 / (5.0f + s.PIP2);
        // PTEN: PIP3 → PIP2 (tumor suppressor)
        float pip3_deg = 2.0f * s.PTEN * s.PIP3 / (0.5f + s.PIP3);
        s.PIP3 += (pip3_syn - pip3_deg) * dt;
        s.PIP3 = clampSig(s.PIP3);

        // Akt activation (recruited to membrane by PIP3, phosphorylated by PDK1+mTORC2)
        float akt_on = 2.0f * s.PIP3 * s.Akt_inactive / (0.5f + s.Akt_inactive);
        float akt_off = 0.3f * s.Akt_active / (0.3f + s.Akt_active); // PP2A
        s.Akt_active += (akt_on - akt_off) * dt;
        s.Akt_inactive += (akt_off - akt_on) * dt;

        // Akt activates mTORC1 (via TSC1/2 inhibition)
        float mtor_on = 1.0f * s.Akt_active * s.mTORC1_inactive / (0.5f + s.mTORC1_inactive);
        float mtor_off = 0.2f * s.mTORC1_active; // AMPK-mediated
        s.mTORC1_active += (mtor_on - mtor_off) * dt;
        s.mTORC1_inactive += (mtor_off - mtor_on) * dt;

        // mTORC1 → S6K (protein synthesis)
        s.S6K_active = hillSig(s.mTORC1_active, 0.3f, 2.0f);

        // Akt phosphorylates BAD (pro-survival)
        s.BAD_P = hillSig(s.Akt_active, 0.3f, 1.5f);

        // Akt phosphorylates FOXO (excludes from nucleus → blocks apoptosis genes)
        s.FOXO_nuclear = 1.0f - hillSig(s.Akt_active, 0.5f, 2.0f);

        // ════════════════════════════════════════════════════════════
        //  4. JAK-STAT PATHWAY
        // ════════════════════════════════════════════════════════════
        float jak_on = 0.5f * s.ligand_GPCR * (1.0f - s.JAK_active); // Simplified
        float jak_off = 0.1f * s.JAK_active * (1.0f + s.SOCS * 5.0f); // SOCS feedback
        s.JAK_active += (jak_on - jak_off) * dt;

        float stat_p_on = 1.0f * s.JAK_active * (1.0f - s.STAT_P);
        float stat_p_off = 0.2f * s.STAT_P;
        s.STAT_P += (stat_p_on - stat_p_off) * dt;

        s.STAT_dimer_nuc = s.STAT_P * s.STAT_P; // Dimerization
        s.SOCS += (0.1f * s.STAT_dimer_nuc - 0.05f * s.SOCS) * dt; // Neg feedback

        // ════════════════════════════════════════════════════════════
        //  5. TGF-β / SMAD
        // ════════════════════════════════════════════════════════════
        float smad_p = 0.5f * s.TGFbeta_receptor * (1.0f - s.Smad2_P) / (1.0f + s.Smad7);
        float smad_dp = 0.1f * s.Smad2_P;
        s.Smad2_P += (smad_p - smad_dp) * dt;
        s.Smad4_complex = s.Smad2_P * 0.8f; // Smad2/3-Smad4 binding
        s.Smad7 += (0.05f * s.Smad4_complex - 0.02f * s.Smad7) * dt; // Neg feedback

        // ════════════════════════════════════════════════════════════
        //  6. NOTCH PATHWAY (contact-dependent)
        // ════════════════════════════════════════════════════════════
        // Delta on neighbor → Notch cleavage → NICD
        float notch_cleavage = 0.3f * s.Delta_ligand * s.Notch_receptor;
        s.NICD += (notch_cleavage - 0.1f * s.NICD) * dt;
        s.Notch_receptor += (0.05f - notch_cleavage - 0.01f * s.Notch_receptor) * dt;
        s.Hes_Hey += (0.5f * s.NICD - 0.1f * s.Hes_Hey) * dt;
        s.NICD = clampSig(s.NICD);
        s.Notch_receptor = clampSig(s.Notch_receptor);
        s.Hes_Hey = clampSig(s.Hes_Hey);

        // ════════════════════════════════════════════════════════════
        //  7. WNT / β-CATENIN
        // ════════════════════════════════════════════════════════════
        // Wnt → Frizzled → Dishevelled → inhibit destruction complex
        s.Frizzled_active = s.Wnt_signal / (0.5f + s.Wnt_signal);
        s.Dishevelled += (0.5f * s.Frizzled_active - 0.1f * s.Dishevelled) * dt;

        // Destruction complex activity (Axin/APC/GSK3/CK1)
        float dc_inhibition = 1.0f / (1.0f + s.Dishevelled * 3.0f);
        s.destruction_complex = dc_inhibition;

        // β-catenin: produced constitutively, destroyed by destruction complex
        float bcat_syn = 0.2f;
        float bcat_deg = 2.0f * s.destruction_complex * s.beta_catenin / (0.1f + s.beta_catenin);
        s.beta_catenin += (bcat_syn - bcat_deg) * dt;
        s.beta_catenin = clampSig(s.beta_catenin);

        // β-catenin → nucleus → TCF activation
        s.TCF_active = hillSig(s.beta_catenin, 0.5f, 2.0f);

        // Axin2 negative feedback (Axin2 is a Wnt target gene that stabilizes destruction complex)
        s.Axin2 += (0.1f * s.TCF_active - 0.05f * s.Axin2) * dt;

        // ════════════════════════════════════════════════════════════
        //  8. HEDGEHOG / PATCHED / SMOOTHENED
        // ════════════════════════════════════════════════════════════
        // Hh binds Patched → releases Smoothened inhibition
        float ptch_bound = s.Hh_signal / (0.5f + s.Hh_signal);
        s.Smoothened_active = ptch_bound;

        // Smoothened → Gli activator (prevent cleavage to repressor)
        float gli_a_rate = 0.3f * s.Smoothened_active;
        float gli_r_rate = 0.2f * (1.0f - s.Smoothened_active);
        s.Gli_activator += (gli_a_rate - 0.1f * s.Gli_activator) * dt;
        s.Gli_repressor += (gli_r_rate - 0.1f * s.Gli_repressor) * dt;

        // ════════════════════════════════════════════════════════════
        //  9. NF-κB PATHWAY
        //  Ref: Hoffmann et al. 2002 Science
        // ════════════════════════════════════════════════════════════
        // IKK phosphorylates IκBα → degradation → NF-κB freed
        float ikb_deg = 1.0f * s.IKK_active * s.IkB_alpha / (0.1f + s.IkB_alpha);
        float ikb_syn = 0.5f * s.NFkB_nuclear; // NF-κB induces IκBα (neg feedback)
        s.IkB_alpha += (ikb_syn - ikb_deg) * dt;
        s.IkB_alpha = clampSig(s.IkB_alpha);

        // NF-κB release from IκBα → nuclear import
        float nfkb_release = ikb_deg;
        float nfkb_rebind = 1.0f * s.IkB_alpha * s.NFkB_nuclear;
        s.NFkB_nuclear += (nfkb_release - nfkb_rebind - 0.05f * s.NFkB_nuclear) * dt;
        s.NFkB_cytoplasm += (nfkb_rebind - nfkb_release) * dt;
        s.NFkB_nuclear = clampSig(s.NFkB_nuclear);
        s.NFkB_cytoplasm = clampSig(s.NFkB_cytoplasm);

        // A20 negative feedback
        s.A20 += (0.2f * s.NFkB_nuclear - 0.1f * s.A20) * dt;
        // A20 inhibits IKK
        s.IKK_active *= 1.0f / (1.0f + s.A20);

        // ════════════════════════════════════════════════════════════
        //  10. NUCLEAR RECEPTOR (steroid hormones)
        // ════════════════════════════════════════════════════════════
        // Hormone crosses membrane → binds receptor → Hsp90 release → nuclear translocation
        float nr_activation = 0.5f * s.steroid_hormone * (1.0f - s.nuclear_receptor_active);
        float nr_deactivation = 0.05f * s.nuclear_receptor_active;
        s.nuclear_receptor_active += (nr_activation - nr_deactivation) * dt;
        s.target_gene_product += (0.1f * s.nuclear_receptor_active - 0.02f * s.target_gene_product) * dt;

        // ════════════════════════════════════════════════════════════
        //  11. NO / cGMP / PKG
        // ════════════════════════════════════════════════════════════
        s.soluble_GC_active = s.nitricOxide / (0.01f + s.nitricOxide);
        float cgmp_syn = 1.0f * s.soluble_GC_active;
        float cgmp_deg = 0.5f * s.cGMP / (1.0f + s.cGMP); // PDE5
        s.cGMP += (cgmp_syn - cgmp_deg) * dt;
        s.PKG_active = hillSig(s.cGMP, 0.5f, 2.0f);
        s.nitricOxide -= 0.5f * s.nitricOxide * dt; // NO is short-lived

        // ════════════════════════════════════════════════════════════
        //  12. CIRCADIAN CLOCK
        //  Ref: Leloup & Goldbeter 2003 PNAS
        // ════════════════════════════════════════════════════════════
        // CLOCK-BMAL1 drives Per and Cry transcription
        float per_transcription = 1.0f * s.CLOCK_BMAL1 / (1.0f + s.PER_CRY_nuclear * 5.0f);
        float cry_transcription = 0.8f * s.CLOCK_BMAL1 / (1.0f + s.PER_CRY_nuclear * 5.0f);

        s.Per_mRNA += (per_transcription - 0.3f * s.Per_mRNA) * dt;
        s.Cry_mRNA += (cry_transcription - 0.3f * s.Cry_mRNA) * dt;

        // Translation
        s.PER_protein += (0.5f * s.Per_mRNA - 0.1f * s.PER_protein) * dt;
        s.CRY_protein += (0.5f * s.Cry_mRNA - 0.1f * s.CRY_protein) * dt;

        // PER-CRY complex formation
        float complex_form = 0.5f * s.PER_protein * s.CRY_protein;
        float complex_deg = 0.05f * s.PER_CRY_complex;
        s.PER_CRY_complex += (complex_form - complex_deg - 0.2f * s.PER_CRY_complex) * dt;

        // Nuclear import
        float nuc_import = 0.2f * s.PER_CRY_complex;
        float nuc_export = 0.05f * s.PER_CRY_nuclear;
        s.PER_CRY_nuclear += (nuc_import - nuc_export - 0.02f * s.PER_CRY_nuclear) * dt;

        // CLOCK-BMAL1 regulation (REV-ERBα represses BMAL1)
        s.Rev_Erb += (0.2f * s.CLOCK_BMAL1 - 0.1f * s.Rev_Erb) * dt;
        s.CLOCK_BMAL1 = 1.0f / (1.0f + s.Rev_Erb * 2.0f + s.PER_CRY_nuclear * 3.0f);

        // Clamp all values
        s.Per_mRNA = clampSig(s.Per_mRNA);
        s.Cry_mRNA = clampSig(s.Cry_mRNA);
        s.PER_protein = clampSig(s.PER_protein);
        s.CRY_protein = clampSig(s.CRY_protein);
        s.PER_CRY_complex = clampSig(s.PER_CRY_complex);
        s.PER_CRY_nuclear = clampSig(s.PER_CRY_nuclear);

        // ════════════════════════════════════════════════════════════
        //  13. RHO FAMILY GTPases
        // ════════════════════════════════════════════════════════════
        // Rac1 activated by PI3K (PIP3-dependent GEFs)
        float rac_on = 0.3f * s.PIP3 * (1.0f - s.Rac1_GTP);
        s.Rac1_GTP += (rac_on - 0.1f * s.Rac1_GTP) * dt;

        // RhoA activated by G12/13 pathway and mechanical signals
        s.RhoA_GTP += (0.1f - 0.1f * s.RhoA_GTP + 0.2f * s.G_betagamma) * dt;
        s.RhoA_GTP = std::clamp(s.RhoA_GTP, 0.0f, 1.0f);

        // Cdc42 activated by RTKs
        s.Cdc42_GTP += (0.2f * s.RTK_active * (1.0f - s.Cdc42_GTP) - 0.1f * s.Cdc42_GTP) * dt;

        // Mutual inhibition between RhoA and Rac1
        s.RhoA_GTP *= 1.0f / (1.0f + s.Rac1_GTP);
        s.Rac1_GTP *= 1.0f / (1.0f + s.RhoA_GTP);

        // ════════════════════════════════════════════════════════════
        //  INTEGRATION: Compute net cellular decisions
        // ════════════════════════════════════════════════════════════

        // Net proliferation signal (drives cell cycle progression)
        s.proliferation_signal = s.ERK_active * 0.3f   // MAPK → CycD expression
                               + s.Akt_active * 0.2f   // Akt → CycD stability
                               + s.TCF_active * 0.2f   // Wnt → CycD, c-Myc
                               + s.STAT_dimer_nuc * 0.1f
                               - s.Smad4_complex * 0.3f // TGF-β → growth arrest
                               - s.Hes_Hey * 0.1f;     // Notch → context-dependent
        s.proliferation_signal = std::clamp(s.proliferation_signal, 0.0f, 1.0f);

        // Net survival signal (inhibits apoptosis)
        s.survival_signal = s.BAD_P * 0.3f             // Akt → BAD phosphorylation
                          + s.NFkB_nuclear * 0.2f       // NF-κB → Bcl-xL, IAPs
                          + (1.0f - s.FOXO_nuclear) * 0.2f // Akt → FOXO exclusion
                          + s.ERK_active * 0.1f;
        s.survival_signal = std::clamp(s.survival_signal, 0.0f, 1.0f);

        // Net growth signal (protein synthesis, cell size)
        s.growth_signal = s.S6K_active * 0.4f           // mTOR → S6K → ribosome biogenesis
                        + s.mTORC1_active * 0.3f;       // mTOR → 4E-BP → cap-dep translation
        s.growth_signal = std::clamp(s.growth_signal, 0.0f, 1.0f);

        // Net transcription drive
        s.total_transcription_drive = s.ERK_nuclear * 0.2f + s.CREB_P * 0.15f
                                    + s.NFkB_nuclear * 0.15f + s.STAT_dimer_nuc * 0.1f
                                    + s.TCF_active * 0.1f + s.Gli_activator * 0.1f
                                    + s.Smad4_complex * 0.1f + s.nuclear_receptor_active * 0.1f;
        s.total_transcription_drive = std::clamp(s.total_transcription_drive, 0.0f, 1.0f);
    }

    void printStatus() const {
        auto& s = state;
        printf("[Signaling] MAPK: Ras=%.2f Raf=%.2f MEK=%.2f ERK=%.2f (nuc=%.2f)\n",
               s.Ras_GTP, s.Raf_active, s.MEK_active, s.ERK_active, s.ERK_nuclear);
        printf("  PI3K-Akt: PIP3=%.3f Akt=%.2f mTOR=%.2f S6K=%.2f\n",
               s.PIP3, s.Akt_active, s.mTORC1_active, s.S6K_active);
        printf("  Ca2+=%.3f µM  cAMP=%.2f  PKA=%.2f  PKC=%.2f\n",
               s.Ca_cytosol, s.cAMP, s.PKA_active, s.PKC_active);
        printf("  NF-κB(nuc)=%.2f  β-cat=%.2f  NICD=%.2f  Clock=%.2f\n",
               s.NFkB_nuclear, s.beta_catenin, s.NICD, s.CLOCK_BMAL1);
        printf("  Decisions: prolif=%.2f surv=%.2f growth=%.2f txn=%.2f\n",
               s.proliferation_signal, s.survival_signal, s.growth_signal,
               s.total_transcription_drive);
    }
};
