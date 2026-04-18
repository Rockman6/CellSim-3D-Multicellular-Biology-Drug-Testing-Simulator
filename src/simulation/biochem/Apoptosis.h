#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Apoptosis.h — Programmed cell death simulation
//  Based on: Alberts Ch 18 "Cell Death"
//
//  Implements:
//    - Intrinsic (mitochondrial) pathway:
//      Bcl-2 family (Bax, Bak, Bcl-2, Bcl-xL, BH3-only: Bid, Bad, Bim, Puma, Noxa)
//      MOMP → cytochrome c release → apoptosome (Apaf-1) → caspase-9
//    - Extrinsic (death receptor) pathway:
//      Fas/TNF-R → FADD → DISC → caspase-8
//      Cross-talk: caspase-8 → tBid → MOMP
//    - Executioner caspases (caspase-3, -6, -7)
//    - IAPs (XIAP) and antagonists (Smac/DIABLO)
//    - Survival signals (PI3K-Akt → BAD phosphorylation)
//    - Phosphatidylserine exposure ("eat me" signal)
//
//  Ref: Albeck et al. 2008 PLoS Biol (EARM model)
//       Spencer & Bhatt 2020 (caspase dynamics)
//
//  Units: concentrations in µM, time in seconds
// ══════════════════════════════════════════════════════════════════════════

struct ApoptosisState {
    // ── Bcl-2 family (mitochondrial gatekeepers) ───────────────────
    // Pro-apoptotic effectors
    float Bax_inactive     = 1.0f;   // Cytosolic, monomeric Bax
    float Bax_active       = 0.0f;   // Activated Bax (oligomeric, in OMM)
    float Bak_inactive     = 0.5f;   // OMM-bound, inactive Bak
    float Bak_active       = 0.0f;   // Activated Bak (oligomeric)

    // Anti-apoptotic (pro-survival)
    float Bcl2             = 1.0f;   // Bcl-2 (OMM)
    float BclXL            = 0.5f;   // Bcl-xL
    float Mcl1             = 0.3f;   // Mcl-1 (short-lived, proteasome target)

    // BH3-only proteins (sensors of death signals)
    float Bid              = 0.5f;   // Full-length Bid
    float tBid             = 0.0f;   // Truncated Bid (activated by caspase-8)
    float Bad              = 0.3f;   // BAD (free, active form)
    float Bad_P            = 0.0f;   // Phosphorylated BAD (inactive, sequestered by 14-3-3)
    float Bim              = 0.2f;   // Bim (sequestered on microtubules normally)
    float Puma             = 0.0f;   // p53-upregulated modulator of apoptosis
    float Noxa             = 0.0f;   // p53-target, specifically inhibits Mcl-1

    // ── Mitochondrial events ───────────────────────────────────────
    float MOMP_pores       = 0.0f;   // Bax/Bak pore formation (0-1)
    float cytochrome_c_mito = 1.0f;  // Cytochrome c in mitochondria
    float cytochrome_c_cyto = 0.0f;  // Cytochrome c in cytosol
    float Smac_mito        = 0.5f;   // Smac/DIABLO in mitochondria
    float Smac_cyto        = 0.0f;   // Smac/DIABLO in cytosol

    // ── Apoptosome and caspase cascade ─────────────────────────────
    float Apaf1            = 0.5f;   // Apaf-1 monomer
    float apoptosome       = 0.0f;   // Assembled apoptosome (Apaf-1 heptamer + cyt c)

    // Initiator caspases
    float caspase9_inactive = 1.0f;  // Procaspase-9
    float caspase9_active   = 0.0f;  // Active caspase-9
    float caspase8_inactive = 1.0f;  // Procaspase-8
    float caspase8_active   = 0.0f;  // Active caspase-8

    // Executioner caspases
    float caspase3_inactive = 1.0f;  // Procaspase-3
    float caspase3_active   = 0.0f;  // Active caspase-3 (THE executioner)
    float caspase6_active   = 0.0f;  // Active caspase-6
    float caspase7_active   = 0.0f;  // Active caspase-7

    // ── Inhibitors of Apoptosis (IAPs) ─────────────────────────────
    float XIAP             = 1.0f;   // X-linked IAP (inhibits caspase-3, -7, -9)
    float XIAP_free        = 1.0f;   // XIAP not bound to Smac

    // ── Death receptor pathway ─────────────────────────────────────
    float FasL             = 0.0f;   // Fas ligand (from cytotoxic T cell)
    float Fas_active       = 0.0f;   // Activated Fas receptor (trimerized)
    float DISC             = 0.0f;   // Death-Inducing Signaling Complex
    float FADD             = 0.5f;   // Fas-Associated Death Domain protein

    // ── Execution phase markers ────────────────────────────────────
    float DNA_fragmentation = 0.0f;  // CAD activity (DNA ladder)
    float lamin_cleavage    = 0.0f;  // Nuclear envelope breakdown
    float PS_exposure       = 0.0f;  // Phosphatidylserine on outer leaflet
    float cell_shrinkage    = 0.0f;  // Volume loss
    float membrane_blebbing = 0.0f;  // Actin/fodrin cleavage

    // ── p53 input ──────────────────────────────────────────────────
    float p53_active       = 0.0f;   // Active p53 (from DNA damage response)

    // ── Survival inputs ────────────────────────────────────────────
    float Akt_survival     = 0.0f;   // From PI3K-Akt pathway
    float NFkB_survival    = 0.0f;   // From NF-κB pathway

    // ── Overall state ──────────────────────────────────────────────
    bool committed         = false;  // Past point of no return?
    float apoptosis_progress = 0.0f; // 0 = alive, 1 = fully apoptotic
};

class ApoptosisEngine {
public:
    ApoptosisState state;

    void init() { state = ApoptosisState(); }

    // ── Apply death signals ─────────────────────────────────────────
    void applyFasLigand(float conc) { state.FasL = conc; }
    void applyDNAdamage(float p53) { state.p53_active = p53; }
    void applySurvival(float akt, float nfkb) {
        state.Akt_survival = akt;
        state.NFkB_survival = nfkb;
    }

    void step(float dt) {
        dt = fminf(dt, 0.01f);
        auto& s = state;

        // ── p53 induces BH3-only proteins ──────────────────────────
        s.Puma += (0.5f * s.p53_active - 0.05f * s.Puma) * dt;
        s.Noxa += (0.3f * s.p53_active - 0.1f * s.Noxa) * dt;
        s.Bim += (0.1f * s.p53_active - 0.02f * s.Bim) * dt;

        // ── Survival: Akt phosphorylates BAD ───────────────────────
        float bad_p_rate = 0.5f * s.Akt_survival * s.Bad;
        float bad_dp_rate = 0.05f * s.Bad_P;
        s.Bad_P += (bad_p_rate - bad_dp_rate) * dt;
        s.Bad += (bad_dp_rate - bad_p_rate) * dt;

        // ── NF-κB upregulates Bcl-xL and XIAP ─────────────────────
        s.BclXL += (0.1f * s.NFkB_survival - 0.02f * s.BclXL) * dt;
        s.XIAP += (0.05f * s.NFkB_survival - 0.01f * s.XIAP) * dt;

        // ── BH3-only neutralize anti-apoptotic proteins ────────────
        // Bad binds Bcl-2/Bcl-xL; Noxa binds Mcl-1; tBid/Bim activate Bax/Bak
        float bcl2_neutralized = s.Bad * 0.3f + s.Puma * 0.2f + s.tBid * 0.5f;
        float effective_Bcl2 = fmaxf(0, s.Bcl2 - bcl2_neutralized);
        float effective_BclXL = fmaxf(0, s.BclXL - s.Bad * 0.2f);
        float effective_Mcl1 = fmaxf(0, s.Mcl1 - s.Noxa * 0.5f);
        float total_antiapoptotic = effective_Bcl2 + effective_BclXL + effective_Mcl1;

        // ── Bax/Bak activation ─────────────────────────────────────
        // Direct activators: tBid, Bim
        float bax_activation = (s.tBid * 1.0f + s.Bim * 0.5f + s.Puma * 0.3f)
                             / (1.0f + total_antiapoptotic * 3.0f); // Anti-apoptotic block
        float bax_on = bax_activation * s.Bax_inactive;
        float bak_on = bax_activation * 0.5f * s.Bak_inactive;
        s.Bax_active += bax_on * dt;
        s.Bax_inactive -= bax_on * dt;
        s.Bak_active += bak_on * dt;
        s.Bak_inactive -= bak_on * dt;

        // ── MOMP (pore formation) ──────────────────────────────────
        // Bax/Bak oligomerize in the OMM — threshold behavior
        float pore_drive = (s.Bax_active + s.Bak_active);
        float pore_form = 2.0f * pore_drive * pore_drive / (0.1f + pore_drive * pore_drive);
        s.MOMP_pores += pore_form * dt;
        s.MOMP_pores = std::clamp(s.MOMP_pores, 0.0f, 1.0f);

        // ── Cytochrome c and Smac release ──────────────────────────
        float release_rate = 5.0f * s.MOMP_pores;
        float cytc_release = release_rate * s.cytochrome_c_mito;
        float smac_release = release_rate * s.Smac_mito;
        s.cytochrome_c_cyto += cytc_release * dt;
        s.cytochrome_c_mito -= cytc_release * dt;
        s.Smac_cyto += smac_release * dt;
        s.Smac_mito -= smac_release * dt;

        // ── Apoptosome assembly ────────────────────────────────────
        // 7 Apaf-1 + 7 cytochrome c → apoptosome (cooperative)
        float apo_form = 1.0f * s.cytochrome_c_cyto * s.Apaf1
                       * s.cytochrome_c_cyto / (0.1f + s.cytochrome_c_cyto);
        s.apoptosome += apo_form * dt;

        // ── Caspase-9 activation (by apoptosome) ───────────────────
        float c9_on = 2.0f * s.apoptosome * s.caspase9_inactive;
        s.caspase9_active += c9_on * dt;
        s.caspase9_inactive -= c9_on * dt;
        // XIAP inhibits caspase-9
        s.caspase9_active *= s.XIAP_free < 0.1f ? 1.0f : (1.0f / (1.0f + s.XIAP_free));

        // ── Extrinsic pathway ──────────────────────────────────────
        // FasL → Fas trimerization → DISC formation
        s.Fas_active = s.FasL / (0.1f + s.FasL);
        s.DISC = s.Fas_active * s.FADD;
        float c8_on = 1.0f * s.DISC * s.caspase8_inactive;
        s.caspase8_active += c8_on * dt;
        s.caspase8_inactive -= c8_on * dt;

        // Cross-talk: caspase-8 cleaves Bid → tBid
        float bid_cleavage = 2.0f * s.caspase8_active * s.Bid;
        s.tBid += bid_cleavage * dt;
        s.Bid -= bid_cleavage * dt;

        // ── Smac antagonizes XIAP ──────────────────────────────────
        s.XIAP_free = fmaxf(0, s.XIAP - s.Smac_cyto);

        // ── Executioner caspase activation ─────────────────────────
        // Caspase-9 and caspase-8 both activate caspase-3
        float c3_on = (2.0f * s.caspase9_active + 1.0f * s.caspase8_active) * s.caspase3_inactive;
        // XIAP inhibits caspase-3
        float xiap_inhibition = s.XIAP_free / (0.1f + s.XIAP_free);
        s.caspase3_active += c3_on * (1.0f - xiap_inhibition) * dt;
        s.caspase3_inactive -= c3_on * dt;

        // Caspase-3 activates caspase-6 and -7
        s.caspase6_active += 0.5f * s.caspase3_active * dt;
        s.caspase7_active += 0.5f * s.caspase3_active * dt;

        // Positive feedback: caspase-3 further activates caspase-9
        s.caspase9_active += 0.3f * s.caspase3_active * s.caspase9_inactive * dt;

        // ── Execution phase ────────────────────────────────────────
        float exec = s.caspase3_active;
        s.DNA_fragmentation += 0.5f * exec * (1.0f - s.DNA_fragmentation) * dt;
        s.lamin_cleavage += 0.3f * exec * (1.0f - s.lamin_cleavage) * dt;
        s.PS_exposure += 0.8f * exec * (1.0f - s.PS_exposure) * dt;
        s.cell_shrinkage += 0.2f * exec * (1.0f - s.cell_shrinkage) * dt;
        s.membrane_blebbing += 0.4f * exec * (1.0f - s.membrane_blebbing) * dt;

        // ── Point of no return ─────────────────────────────────────
        if (s.caspase3_active > 0.5f) s.committed = true;

        // Overall progress
        s.apoptosis_progress = (s.DNA_fragmentation + s.lamin_cleavage + s.PS_exposure
                              + s.cell_shrinkage + s.membrane_blebbing) / 5.0f;

        // Clamp
        s.Bax_inactive = fmaxf(0, s.Bax_inactive);
        s.Bak_inactive = fmaxf(0, s.Bak_inactive);
        s.Bid = fmaxf(0, s.Bid);
        s.caspase9_inactive = fmaxf(0, s.caspase9_inactive);
        s.caspase8_inactive = fmaxf(0, s.caspase8_inactive);
        s.caspase3_inactive = fmaxf(0, s.caspase3_inactive);
        s.cytochrome_c_mito = fmaxf(0, s.cytochrome_c_mito);
        s.Smac_mito = fmaxf(0, s.Smac_mito);
    }

    bool isApoptotic() const { return state.committed; }
    float progress() const { return state.apoptosis_progress; }
};
