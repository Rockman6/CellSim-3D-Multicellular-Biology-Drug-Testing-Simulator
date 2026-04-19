#pragma once
#include "ChemicalEntity.h"
#include <cmath>

// ══════════════════════════════════════════════════════════════════════════
//  BindingMatcher.h — Tier-0 descriptor-based drug × target affinity.
//
//  Emergent MOA: the sim doesn't know a drug is "anti-proliferative" or
//  "pro-apoptotic". BindingMatcher computes a compatibility score
//  between the drug's physicochemical profile and the target's preferred
//  zone; mass-action kinetics + modulator cascades do the rest.
//
//  Three terms (all 0..1):
//    (1) Descriptor compatibility — Lipinski-style fit against
//        target's prefer_{logP,mw,hbd,hba,aromatic_rings} ranges.
//    (2) Fingerprint overlap — Tanimoto (stub; real impl in Phase-C).
//    (3) 3D pharmacophore RMSD — skipped at Tier 0.
//
//  Kd in mM = 10^(-6 * score); score = 1 → 1 µM, score = 0 → 1 M.
//  Returns a BindingAffinity populated with Kd + on/off rates.
// ══════════════════════════════════════════════════════════════════════════

// A minimal struct describing a target's preferred drug profile.
// Loaded from data/bioagents/targets.csv by BioagentRegistry.
struct TargetProfile {
    std::string id;
    float pref_logP_min = -3.0f, pref_logP_max = 5.0f;
    float pref_mw_min   = 100.0f, pref_mw_max   = 900.0f;
    int   pref_hbd      = 2;
    int   pref_hba      = 5;
    int   pref_aromatic = 1;
};

namespace BindingMatcher {

inline float hillFit(float value, float lo, float hi) {
    // 1.0 if value ∈ [lo, hi]; falls off outside.
    if (value >= lo && value <= hi) return 1.0f;
    float center = (lo + hi) * 0.5f;
    float width  = (hi - lo) * 0.5f + 1e-3f;
    float d = fabsf(value - center) / width;   // 0 at center, 1 at edge, >1 outside
    return expf(-(d - 1.0f) * (d - 1.0f));     // soft Gaussian falloff
}

inline float integerFit(int value, int preferred) {
    int d = value - preferred;
    return expf(-0.25f * (float)(d * d));      // sharp fit on small integer counts
}

// Main entry point. Consumes Tier-0 descriptor bundle from the drug
// and the TargetProfile from targets.csv.
inline BindingAffinity score(const ChemicalEntity& drug,
                              const TargetProfile& target) {
    float logP_fit = hillFit(drug.logP, target.pref_logP_min, target.pref_logP_max);
    float mw_fit   = hillFit(drug.mw,   target.pref_mw_min,   target.pref_mw_max);
    float hbd_fit  = integerFit(drug.hbd, target.pref_hbd);
    float hba_fit  = integerFit(drug.hba, target.pref_hba);
    float arom_fit = integerFit(drug.aromatic_rings, target.pref_aromatic);

    // Geometric mean of descriptor fits — any zero kills the score.
    float descriptor_score = powf(
        fmaxf(1e-6f, logP_fit * mw_fit * hbd_fit * hba_fit * arom_fit),
        1.0f / 5.0f);

    // Fingerprint term (placeholder — set to 0.5 when we don't have
    // target fingerprint yet; ignored in Tier-0).
    float fp_score = 0.5f;

    // Combined score weighted toward descriptor fit at Tier 0.
    float score = 0.7f * descriptor_score + 0.3f * fp_score;
    score = fmaxf(0.0f, fminf(1.0f, score));

    BindingAffinity aff;
    // Kd in mM — calibrated so realistic descriptor-fit scores give
    // pharmacologically-meaningful affinities:
    //   score 0.00 → 100 mM (no binding)
    //   score 0.25 →   1 mM (weak)
    //   score 0.50 →  10 µM (moderate, typical drug range)
    //   score 0.75 → 100 nM (strong)
    //   score 1.00 →   1 nM (exceptional, rarely reached by Tier 0 alone)
    // Tier 3 docking / Tier 4 MD will refine these into real values.
    aff.Kd_mM   = powf(10.0f, 2.0f - 8.0f * score);   // 10^(2 - 8·score)
    // Diffusion-limited on-rate capped at ~1e8 M^-1 s^-1 (Berg-Purcell)
    aff.k_on_per_uM = 1.0e-2f;      // per µM·s
    aff.k_off_per_s = aff.k_on_per_uM * aff.Kd_mM * 1e3f; // µM·Kd in µM
    aff.score   = score;
    aff.tier    = 0;
    return aff;
}

} // namespace BindingMatcher
