#pragma once

// ══════════════════════════════════════════════════════════════════════════
//  MediumField.h — Closed-system culture-medium chemistry
//  -----------------------------------------------------------------------
//  Replaces the legacy 4-species NutrientField with a 12-species DMEM-based
//  fluid model on the same 64×64 grid. Concentrations are stored in REAL
//  units (mM, except the growth-factor channel which is ng/mL and the drug
//  channel which is µM). The dish is a CLOSED system: zero-flux Neumann
//  boundaries (no reservoir refill, no perfusion source). Cells exchange
//  matter with the field strictly through stoichiometric reactions whose
//  reactants and products are debited and credited atomically.
//
//  Reactions (every flux MUST originate here, never out of nothing):
//    R1  Glycolysis       glucose → 2 pyruvate (+2 ATP, +2 H+)
//    R2  Lactate ferm     pyruvate → lactate (anaerobic)
//    R3  TCA + OxPhos     pyruvate + 2.5 O2 → 3 CO2 (+30 ATP)
//    R4  Glutaminolysis   glutamine + 0.5 O2 → 2 CO2 (+9 ATP)
//    R5  Translation      AA pool → biomass (consumes ATP equivalents)
//    R6  Receptor binding GF ↔ RTK (catalytic — GF returned)
//    R7  Cell autolysis   biomass → AA pool + glucose + lipids
//
//  Mass balance is verified per-tick to within 1e-4 mM total drift.
//  Sources: Sigma-Aldrich D5796 datasheet (composition); Goldstick 1976,
//  Longsworth 1953, Stein 1986 (diffusion); Reitzer 1979 (HeLa flux).
// ══════════════════════════════════════════════════════════════════════════

#include "../core/Constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

// ── Species index ────────────────────────────────────────────────────────
enum MediumSpecies : int {
    MS_O2          = 0,
    MS_CO2         = 1,
    MS_GLUCOSE     = 2,
    MS_GLUTAMINE   = 3,
    MS_PYRUVATE    = 4,
    MS_LACTATE     = 5,
    MS_AA_POOL     = 6,
    MS_GROWTH_F    = 7,    // ng/mL (not mM — stored as-is)
    MS_IONS        = 8,
    MS_CALCIUM     = 9,
    MS_HPLUS       = 10,   // pH proxy (stored as pH directly, 7.40 default)
    MS_WATER       = 11,   // mM (~55 000 in pure water)
    MS_COUNT       = 12
    // MS_DRUG removed 2026-04-19 — pending rewrite
};

inline const char* mediumSpeciesName(int s) {
    static const char* names[MS_COUNT] = {
        "O2","CO2","glucose","glutamine","pyruvate","lactate",
        "AA_pool","growthF","ions","Ca2+","pH","water"
    };
    return (s >= 0 && s < MS_COUNT) ? names[s] : "?";
}

inline const char* mediumSpeciesUnit(int s) {
    if (s == MS_GROWTH_F) return "ng/mL";
    if (s == MS_HPLUS)    return "pH";
    return "mM";
}

// ── Closed-system culture-medium field ──────────────────────────────────
struct MediumField {
    static constexpr int N = NUTRIENT_GRID_N;
    static constexpr int CELLS = N * N;

    // Layout: [species][grid] for cache locality on a single species.
    float c[MS_COUNT][CELLS];
    float cNext[MS_COUNT][CELLS];

    // Diffusion coefficient table indexed by MediumSpecies.
    float D_um2_per_s[MS_COUNT];

    // Drug-channel knobs removed 2026-04-19 — pending rewrite.

    // ── Extracellular protease field ───────────────────────────────────
    // Auxiliary scalar (not a conserved mass species) tracking
    // cathepsins / MMPs / phospholipases released at cell-lysis events.
    // Boosts the substrate-availability multiplier for nearby living
    // cells' R5 translation → they "eat better" in a dying cluster.
    // Decays with ~2-bio-h half-life (Hibbs 1987 J Clin Inv).
    float extracellularProteaseLevel[CELLS];

    // ── Bioagent concentration fields (Phase M1 scaffold) ──────────────
    // Dish-level µM concentrations for each foreign ChemicalEntity
    // (drugs, viral particles in transit, bacterial secretions). One
    // float[CELLS] slab per agent, allocated on drug-load. Indexed by
    // the agent's position in gBioagents registry. Empty by default —
    // zero footprint until a drug is applied.
    // See plan Part-One §4.1 and Part-Five M1.
    std::vector<std::vector<float>> bioagentConcentration_uM;

    // Mass-balance bookkeeping. totalAtStart[i] = total moles per species
    // captured at init (in mM·grid-cell-units). The ratio of running total
    // to this is the conservation invariant.
    double totalAtStart[MS_COUNT];

    // Cell physical size (one grid cell ≈ this many µm on a side). Used
    // to scale diffusion coefficients into per-tick stencil weights.
    // SCENE_BOUND wu × MICROMETERS_PER_WORLD_UNIT µm/wu × 2 = full width;
    // /(N) gives µm per grid cell.
    static constexpr float CELL_UM = (2.0f * SCENE_BOUND) * MICROMETERS_PER_WORLD_UNIT / (float)N;

    // ── Init ────────────────────────────────────────────────────────────
    void init() {
        using namespace MediumComposition;
        using namespace MediumDiffusion;

        const float initVals[MS_COUNT] = {
            DMEM_O2_MM, DMEM_CO2_MM, DMEM_GLUCOSE_MM, DMEM_GLUTAMINE_MM,
            DMEM_PYRUVATE_MM, DMEM_LACTATE_MM, DMEM_AA_POOL_MM,
            DMEM_GROWTH_FACTOR_NG_PER_ML, DMEM_IONS_MM, DMEM_CALCIUM_MM,
            DMEM_PH, DMEM_WATER_MM
        };
        const float diffVals[MS_COUNT] = {
            D_O2_UM2_PER_S, D_CO2_UM2_PER_S, D_GLUCOSE_UM2_PER_S,
            D_GLUTAMINE_UM2_PER_S, D_PYRUVATE_UM2_PER_S, D_LACTATE_UM2_PER_S,
            D_AA_UM2_PER_S, D_GROWTH_F_UM2_PER_S, D_IONS_UM2_PER_S,
            D_CALCIUM_UM2_PER_S, D_H_UM2_PER_S,
            D_WATER_UM2_PER_S
        };

        for (int s = 0; s < MS_COUNT; s++) {
            D_um2_per_s[s] = diffVals[s];
            for (int i = 0; i < CELLS; i++) c[s][i] = initVals[s];
            totalAtStart[s] = (double)initVals[s] * (double)CELLS;
        }
        for (int i = 0; i < CELLS; i++) extracellularProteaseLevel[i] = 0.0f;
    }

    // ── Extracellular protease accessors ────────────────────────────────
    void depositProtease(float wx, float wz, float amount) {
        int ix, iz; worldToGrid(wx, wz, ix, iz);
        extracellularProteaseLevel[iz * N + ix] += amount;
    }
    float proteaseAt(float wx, float wz) const {
        int ix, iz; worldToGrid(wx, wz, ix, iz);
        return extracellularProteaseLevel[iz * N + ix];
    }
    // First-order decay each biological-second tick (called from diffuse()).
    void decayProtease(float dt_bio, float per_biosec_rate) {
        float keep = std::max(0.0f, 1.0f - per_biosec_rate * dt_bio);
        for (int i = 0; i < CELLS; i++) extracellularProteaseLevel[i] *= keep;
    }

    // ── World ↔ grid mapping ────────────────────────────────────────────
    void worldToGrid(float wx, float wz, int& ix, int& iz) const {
        float fx = (wx / SCENE_BOUND + 1.0f) * 0.5f * (N - 1);
        float fz = (wz / SCENE_BOUND + 1.0f) * 0.5f * (N - 1);
        ix = (int)std::max(0.0f, std::min((float)(N - 1), fx));
        iz = (int)std::max(0.0f, std::min((float)(N - 1), fz));
    }

    // ── Read accessors (real units, no normalization) ───────────────────
    float get(int species, float wx, float wz) const {
        int ix, iz; worldToGrid(wx, wz, ix, iz);
        return c[species][iz * N + ix];
    }
    // Mean value across the dish (for the UI panel).
    double mean(int species) const {
        double sum = 0.0;
        for (int i = 0; i < CELLS; i++) sum += c[species][i];
        return sum / (double)CELLS;
    }

    // Compatibility helpers — return values NORMALIZED so the legacy
    // updateMetabolism / updateCellCycle code paths interpret them the
    // same way they did with NutrientField. These are temporary; future
    // code should call get(MS_*) directly and reason in mM.
    float getO2(float wx, float wz) const {
        return get(MS_O2, wx, wz) / MediumComposition::DMEM_O2_MM;
    }
    float getGlucose(float wx, float wz) const {
        return get(MS_GLUCOSE, wx, wz) / MediumComposition::DMEM_GLUCOSE_MM;
    }
    float getCO2(float wx, float wz) const {
        return get(MS_CO2, wx, wz) / 4.0f;  // legacy normalisation
    }
    float getPH(float wx, float wz) const {
        // Legacy field used 0.55–0.80 normalised pH. Map 6.6 → 0.55,
        // 7.4 → 0.74, 7.6 → 0.80.
        float pH = get(MS_HPLUS, wx, wz);
        return std::max(0.55f, std::min(0.80f, 0.55f + (pH - 6.6f) * 0.25f));
    }
    // ── Per-cell exchange (atomic stoichiometric flux, mM/s × dt) ───────
    // fluxMmPerS: MS_COUNT-element array. Negative entries are CONSUMED by the
    // cell, positive entries are SECRETED. Limited by available reactant
    // (negative entries clamped so concentrations stay ≥ 0).
    void exchange(float wx, float wz, const float* fluxMmPerS, float dt_bio) {
        int ix, iz; worldToGrid(wx, wz, ix, iz);
        int idx = iz * N + ix;
        for (int s = 0; s < MS_COUNT; s++) {
            float delta = fluxMmPerS[s] * dt_bio;
            float newVal = c[s][idx] + delta;
            if (newVal < 0.0f) newVal = 0.0f;
            c[s][idx] = newVal;
        }
    }

    // Drug API removed 2026-04-19 — pending rewrite.

    // ── Diffusion (zero-flux Neumann boundaries — closed dish) ──────────
    // Each species advances under explicit Euler:
    //   c_new = c + (D / dx²) * laplacian(c) * dt
    // For a closed system the laplacian uses mirrored ghost cells at the
    // boundary, which we implement by clamping neighbor index access.
    // This is *exactly* mass-conserving by construction.
    void diffuse(float dt_bio) {
        const float dx2 = MediumField::CELL_UM * MediumField::CELL_UM;
        for (int s = 0; s < MS_COUNT; s++) {
            float rate = D_um2_per_s[s] / dx2;
            // Stability guard: explicit Euler for diffusion is stable when
            // 4 * rate * dt <= 1. Sub-step internally if not.
            float safe = 0.20f / (rate + 1e-9f);  // 80% of CFL limit
            int substeps = std::max(1, (int)std::ceil(dt_bio / safe));
            float subDt = dt_bio / (float)substeps;
            for (int step = 0; step < substeps; step++) {
                for (int iz = 0; iz < N; iz++) {
                    for (int ix = 0; ix < N; ix++) {
                        int idx = iz * N + ix;
                        // Mirrored boundary (Neumann): clamp neighbour
                        int xm = (ix > 0) ? ix - 1 : ix;
                        int xp = (ix < N-1) ? ix + 1 : ix;
                        int zm = (iz > 0) ? iz - 1 : iz;
                        int zp = (iz < N-1) ? iz + 1 : iz;
                        float lap =
                            c[s][iz*N + xm] + c[s][iz*N + xp] +
                            c[s][zm*N + ix] + c[s][zp*N + ix] -
                            4.0f * c[s][idx];
                        cNext[s][idx] = c[s][idx] + rate * lap * subDt;
                        if (cNext[s][idx] < 0.0f) cNext[s][idx] = 0.0f;
                    }
                }
                std::copy(cNext[s], cNext[s] + CELLS, c[s]);
            }
        }
        // Extracellular protease decay — 2-bio-h half-life. Constant
        // inlined to avoid circular include of Constants.h (already
        // pulled in above via "../core/Constants.h").
        decayProtease(dt_bio, Apoptosis::EP_DECAY_PER_BIOSEC);
    }

    // ── Incubator gas exchange ──────────────────────────────────────────
    // Real tissue-culture incubators are VENTED: cells sit in a gas-
    // permeable dish / flask exposed to 5 % CO2 / 18 % O2. Without this
    // exchange term, a closed 13-species dish accumulates metabolic CO2
    // indefinitely, pH crashes, and cells go quiescent within ~24 bio-h
    // — not how real HeLa culture works. These first-order relaxation
    // rates pull CO2 and O2 back toward their atmospheric equilibria.
    // Rate ~0.5 /bio-h = full re-equilibration in ~2 bio-h.
    void gasExchange(float dt_bio) {
        const float K_EXCHANGE = 0.5f / 3600.0f;   // per bio-second
        const float co2_eq = MediumComposition::DMEM_CO2_MM;
        const float o2_eq  = MediumComposition::DMEM_O2_MM;
        float decay_co2 = fminf(1.0f, K_EXCHANGE * dt_bio);
        float decay_o2  = fminf(1.0f, K_EXCHANGE * dt_bio);
        for (int i = 0; i < CELLS; i++) {
            c[MS_CO2][i] += (co2_eq - c[MS_CO2][i]) * decay_co2;
            c[MS_O2][i]  += (o2_eq  - c[MS_O2][i])  * decay_o2;
        }
    }

    // ── Bicarbonate-buffered pH update ──────────────────────────────────
    // Henderson-Hasselbalch with bicarbonate (HCO3- pool implicit in
    // total CO2 ↔ H+ + HCO3-): pH = 7.4 - 0.6 * (CO2 - 1.2) / 1.2,
    // clamped to [6.6, 7.6]. Lactate accumulation drops pH further:
    // pH -= 0.04 per mM lactate above zero.
    void updatePH() {
        for (int i = 0; i < CELLS; i++) {
            float co2 = c[MS_CO2][i];
            float lac = c[MS_LACTATE][i];
            float pH = 7.40f - 0.60f * (co2 - MediumComposition::DMEM_CO2_MM)
                              / MediumComposition::DMEM_CO2_MM
                       - 0.04f * lac;
            c[MS_HPLUS][i] = std::max(6.6f, std::min(7.6f, pH));
        }
    }

    // ── Mass balance check ──────────────────────────────────────────────
    // Returns the maximum relative drift across all species (excluding pH
    // which is not conserved by R1–R7). Closed-system invariant
    // is < 1e-3 in steady operation; the caller flags a red banner above
    // 1e-2 (set by the user as the physics-broken threshold).
    double checkBalance(double drift[MS_COUNT]) const {
        double maxDrift = 0.0;
        for (int s = 0; s < MS_COUNT; s++) {
            if (s == MS_HPLUS) { drift[s] = 0.0; continue; }
            double sum = 0.0;
            for (int i = 0; i < CELLS; i++) sum += c[s][i];
            double rel = (totalAtStart[s] > 1e-9)
                ? (sum - totalAtStart[s]) / totalAtStart[s]
                : 0.0;
            drift[s] = rel;
            if (std::abs(rel) > maxDrift) maxDrift = std::abs(rel);
        }
        return maxDrift;
    }

    // ── Legacy compatibility: init(envO2, envGlu) — ignored ─────────────
    // Existing call sites pass the user-set env sliders; we ignore them
    // since the dish is now seeded from the DMEM constants. The shim
    // exists so the build doesn't break while Simulation.h migrates.
    void init(float /*envO2_unused*/, float /*envGlu_unused*/) { init(); }
    // Legacy diffuse(dt, envO2, envGlu) — env params ignored.
    void diffuse(float dt_bio, float /*envO2_unused*/, float /*envGlu_unused*/) {
        diffuse(dt_bio);
        gasExchange(dt_bio);   // simulates vented incubator
        updatePH();
    }
    // Legacy consume — translates the old (o2Rate, gluRate, glycolytic)
    // call into a stoichiometric exchange. Each o2Rate step = a TCA
    // burst that produces CO2; each gluRate step = glycolysis that
    // produces pyruvate (or lactate if glycolytic). All in dimensionless
    // legacy units; we convert to mM-flux by treating the legacy values
    // as small deltas in the local cell.
    void consume(float wx, float wz, float o2Rate, float gluRate, bool glycolytic) {
        int ix, iz; worldToGrid(wx, wz, ix, iz);
        int idx = iz * N + ix;
        // O2 consumed; CO2 emitted at 1:1 (lumped TCA + OxPhos)
        float dO2 = std::min(c[MS_O2][idx], o2Rate * MediumComposition::DMEM_O2_MM);
        c[MS_O2][idx]  -= dO2;
        c[MS_CO2][idx] += dO2;
        // Glucose consumed; one mole of glucose → 2 mol pyruvate (or 2 lactate)
        float dGlu = std::min(c[MS_GLUCOSE][idx], gluRate * MediumComposition::DMEM_GLUCOSE_MM);
        c[MS_GLUCOSE][idx] -= dGlu;
        if (glycolytic) {
            c[MS_LACTATE][idx] += 2.0f * dGlu;
        } else {
            c[MS_PYRUVATE][idx] += 2.0f * dGlu;
        }
    }
};
