#pragma once
#include "Virion.h"
#include "Bacterium.h"
#include <cmath>
#include <string>

// ══════════════════════════════════════════════════════════════════════════
//  PathogenKinetics.h — Pure scoring functions for pathogen-host binding.
//
//  Same physics-first discipline as BindingMatcher.h: a virion spike
//  binds a host receptor because the 5-D Lipinski-like descriptor vector
//  of the spike is compatible with a reference vector for that receptor
//  class. No hardcoded tropism.
//
//  Output: a binding "score" in 0..1. Higher = tighter binding. A score
//  of 0 means Kd ≈ 100 mM (no engagement); 1 means Kd ≈ 1 nM. This score
//  is consumed by the virion / bacterium state machine inside Simulation.
// ══════════════════════════════════════════════════════════════════════════

namespace PathogenKinetics {

// Built-in receptor-class descriptor profiles. A host cell carries a
// *mosaic* — a density per class, multiplied by the per-class profile
// compatibility. Values are curated to approximate known receptor-
// pocket preferences and to differentiate receptor classes.
//
// id               logP    mw(Da)   hbd  hba  arom
// PT_RECEPTOR_GPCR  2.5    45000    5    12   3
// PT_RECEPTOR_RTK   2.0    70000    6    14   3
// PT_RECEPTOR_DEATH 1.2    50000    7    10   1
// PT_RECEPTOR_TLR   1.8    95000    9    18   4
// PT_RECEPTOR_CYTO  1.6    55000    6    12   2
// PT_RECEPTOR_INTEG 2.2    110000   10   20   3
// (ICAM / sialic-acid — a generic sugar-binding site)
// PT_GLYCOPROTEIN   0.5    20000    12   18   0
struct ReceptorProfile {
    const char* id;
    float logP, mw;
    int   hbd, hba, aromatic;
};

inline const ReceptorProfile* receptorProfiles(int& count) {
    static const ReceptorProfile P[] = {
        {"PT_RECEPTOR_GPCR",     2.5f, 45000.0f,  5, 12, 3},
        {"PT_RECEPTOR_RTK",      2.0f, 70000.0f,  6, 14, 3},
        {"PT_RECEPTOR_DEATH",    1.2f, 50000.0f,  7, 10, 1},
        {"PT_RECEPTOR_TLR",      1.8f, 95000.0f,  9, 18, 4},
        {"PT_RECEPTOR_CYTOKINE", 1.6f, 55000.0f,  6, 12, 2},
        {"PT_RECEPTOR_INTEGRIN", 2.2f,110000.0f, 10, 20, 3},
        {"PT_GLYCOPROTEIN",      0.5f, 20000.0f, 12, 18, 0},
    };
    count = (int)(sizeof(P) / sizeof(P[0]));
    return P;
}

// Gaussian fit (unit-less); wider = more tolerant. Reused from Part 1.
inline float hillFitLP(float value, float ref, float width) {
    float d = (value - ref) / (width + 1e-3f);
    return expf(-0.5f * d * d);
}
inline float intFit(int value, int ref) {
    int d = value - ref;
    return expf(-0.20f * (float)(d * d));
}

// Score a 5-D descriptor against one receptor profile.
inline float compatScore5D(float logP, float mw, int hbd, int hba, int arom,
                            const ReceptorProfile& r) {
    float a = hillFitLP(logP, r.logP, 1.2f);
    float b = hillFitLP(mw, r.mw, r.mw * 0.40f + 1.0f);
    float c = intFit(hbd, r.hbd);
    float d = intFit(hba, r.hba);
    float e = intFit(arom, r.aromatic);
    return powf(fmaxf(1e-6f, a * b * c * d * e), 1.0f / 5.0f);
}

// Best compatibility across a list of preferred receptors. Returns the
// max score and writes the winning receptor id. Density weighting is
// applied by the caller (multiplies score by receptor density).
inline float bestSpikeScore(float logP, float mw, int hbd, int hba, int arom,
                             const std::vector<std::string>& preferred,
                             std::string& outReceptorId) {
    int nProfiles;
    const ReceptorProfile* P = receptorProfiles(nProfiles);
    float bestS = 0.0f;
    outReceptorId.clear();
    // If no preference list, score against ALL profiles.
    if (preferred.empty()) {
        for (int i = 0; i < nProfiles; i++) {
            float s = compatScore5D(logP, mw, hbd, hba, arom, P[i]);
            if (s > bestS) { bestS = s; outReceptorId = P[i].id; }
        }
        return bestS;
    }
    for (const std::string& want : preferred) {
        for (int i = 0; i < nProfiles; i++) {
            if (want != P[i].id) continue;
            float s = compatScore5D(logP, mw, hbd, hba, arom, P[i]);
            if (s > bestS) { bestS = s; outReceptorId = P[i].id; }
        }
    }
    return bestS;
}

// Convenience wrappers for spec types.
inline float virionBindScore(const VirionSpec& v, std::string& rcvOut) {
    return bestSpikeScore(v.spike_logP, v.spike_mw, v.spike_hbd, v.spike_hba,
                           v.spike_aromatic, v.preferredReceptors, rcvOut);
}
inline float bacteriumBindScore(const BacteriumSpec& b, std::string& rcvOut) {
    return bestSpikeScore(b.adhesin_logP, b.adhesin_mw, b.adhesin_hbd, b.adhesin_hba,
                           b.adhesin_aromatic, b.preferredReceptors, rcvOut);
}

// Translate score (0..1) to a per-bio-second effective on-rate (1/s)
// multiplier on the spike density. Calibrated so score 0.7 → k_eff ~ 1/s
// meaning a typical contact event progresses binding within ~1 bio-s.
inline float scoreToOnRate(float score) {
    // k = 10^(3·(score - 0.7))  → 0.1/s at 0.3, 1/s at 0.7, 10/s at 1.0
    return powf(10.0f, 3.0f * (score - 0.7f));
}

} // namespace PathogenKinetics
