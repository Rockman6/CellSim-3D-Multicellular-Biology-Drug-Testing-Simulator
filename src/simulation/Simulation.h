#pragma once

#include "../core/Constants.h"
#include "CellCycleProgram.h"
#include "MediumField.h"
#include <simd/simd.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <cstdio>

// ══════════════════════════════════════════════════════════════════════════
//  Simulation.h — Cell biology simulation + pharmacology platform
//
//  Biology: CDK/Cyclin ODE, metabolism, Warburg, hypoxia, contact inhibition
//  Pharmacology: Drug diffusion, uptake, PD damage model, IC50/Hill response
//
//  Literature:
//    Novak & Tyson (2008) Nat Rev MCB — CDK/Cyclin oscillator
//    Casciari et al (1992) Biotechnol Bioeng — Fick diffusion
//    Green & Kroemer (2004) Science — ATP-dependent apoptosis
//    Delarue et al (2018) Dev Cell — mechanical contact inhibition
//    Blackburn (2001) Cell — telomere erosion / senescence
//
//  Drug PD model adapted from PhysiPKPD (BSD-3 license)
//    Bergman et al., GigaByte 2023. github.com/drbergman/PhysiPKPD
//  Cell cycle drug arrest: PhysiCell (BSD-3)
//    Ghaffarizadeh et al., PLoS Comput Biol 2018
//  IC50 values: PMC2751448, PMC1501422, PMC5342590
//  Drug parameters: BioModels Database (CC0), ebi.ac.uk/biomodels
// ══════════════════════════════════════════════════════════════════════════

static float randf() { return (float)rand() / (float)RAND_MAX; }
// clampf is defined in CentralDogma.h (included above via CellCycleProgram.h)
static float biomassFromCellSizeFactor(float sizeFactor) {
    // In updateMetabolism(): size = 0.6 + biomass * 0.4.
    // Keep daughters internally consistent with the rendered daughter size.
    return clampf((sizeFactor - 0.6f) / 0.4f, 0.4f, 2.3f);
}

// Resource split ratio from physical particle distribution
struct SplitStats {
    float cytoplasmicRatioA = 0.5f; // Mass-weighted non-nuclear particle share for daughter A
    int dnaA = 0, dnaB = 0;         // DNA node counts (verification only)
};

// ── Drug definition ─────────────────────────────────────────────────────
// PK/PD model: PhysiPKPD framework (Bergman et al. 2023)
// dC_int/dt = uptake × C_ext − efflux × C_int
// dDamage/dt = Hill(C_int, EC50, n) × maxEffect − repairRate × Damage
// Effect depends on MOA: cycle arrest, apoptosis, DNA damage, mito toxicity
struct Drug {
    const char* name;
    // Pharmacokinetics
    float diffusionCoeff;   // Fick diffusion in tissue (sim units)
    float decayRate;        // spontaneous degradation per dt
    float uptakeRate;       // cell internalization rate
    float effluxRate;       // MDR pump export rate
    // Pharmacodynamics (Hill equation)
    float EC50;             // half-maximal effective conc (µM)
    float hillCoeff;        // Hill coefficient (dose-response steepness)
    float maxEffect;        // max fractional effect (0-1)
    int   mechanism;        // MOA code (MOA_ANTI_PROLIF etc.)
    int   mechanism2;       // secondary MOA (-1 = none)
    // Repair
    float damageRepairRate; // cellular repair of drug damage
    float resistanceMutRate;// mutation prob per division
};

// Pre-built drug library with published parameters
// Ref: PMC2751448 (cisplatin/paclitaxel), PMC1501422 (doxorubicin)
static const Drug DRUG_LIBRARY[] = {
    // name              diffuse decay  uptake efflux EC50   hill  maxEff MOA1              MOA2              repair resistMut
    {"None (Control)",   0,      0,     0,     0,     999,   1.0f, 0,     -1,               -1,               0,     0},
    {"Cisplatin",        1.2f,   0.002f,0.04f, 0.005f,2.0f,  1.5f, 0.9f,  MOA_DNA_DAMAGE,   -1,               0.008f,0.002f},
    {"Doxorubicin",      0.8f,   0.003f,0.06f, 0.01f, 0.5f,  2.0f, 0.95f, MOA_DNA_DAMAGE,   MOA_MITO_TOXIN,   0.005f,0.003f},
    {"Paclitaxel",       0.6f,   0.001f,0.08f, 0.008f,0.01f, 2.5f, 0.85f, MOA_ANTI_PROLIF,  -1,               0.01f, 0.001f},
    {"5-Fluorouracil",   1.5f,   0.004f,0.05f, 0.006f,5.0f,  1.2f, 0.8f,  MOA_ANTI_PROLIF,  MOA_DNA_DAMAGE,   0.012f,0.002f},
};
static const int DRUG_COUNT = sizeof(DRUG_LIBRARY) / sizeof(Drug);

// ── Legacy NutrientField (REPLACED by MediumField — kept dead so any
// stray reference produces a compile error rather than silent drift)
#if 0
struct NutrientField {
    static constexpr int N = NUTRIENT_GRID_N;
    float o2[N*N], glucose[N*N], co2[N*N], pH_field[N*N];
    float o2Next[N*N], gluNext[N*N], co2Next[N*N], pHNext[N*N];
    // Drug concentration field (µM)
    float drug[N*N], drugNext[N*N];
    float drugDiffCoeff = 0;  // set when drug is applied
    float drugDecayRate = 0;

    void init(float envO2, float envGlu) {
        for (int iz = 0; iz < N; iz++)
            for (int ix = 0; ix < N; ix++) {
                float dx = (float)(ix-N/2)/(N/2), dz = (float)(iz-N/2)/(N/2);
                float r = sqrtf(dx*dx+dz*dz);
                float edge = 0.55f + 0.45f * r;
                int idx = iz*N+ix;
                o2[idx] = envO2 * edge;
                glucose[idx] = envGlu * edge;
                co2[idx] = (1-edge)*0.15f + 0.02f;
                pH_field[idx] = 0.74f + (1-edge)*(-0.04f);
                drug[idx] = 0; drugNext[idx] = 0;
            }
    }

    void worldToGrid(float wx, float wz, int& ix, int& iz) const {
        ix = (int)clampf((wx/SCENE_BOUND+1)*0.5f*(N-1), 0, N-1);
        iz = (int)clampf((wz/SCENE_BOUND+1)*0.5f*(N-1), 0, N-1);
    }

    float getO2(float wx, float wz) const { int ix,iz; worldToGrid(wx,wz,ix,iz); return o2[iz*N+ix]; }
    float getGlucose(float wx, float wz) const { int ix,iz; worldToGrid(wx,wz,ix,iz); return glucose[iz*N+ix]; }
    float getCO2(float wx, float wz) const { int ix,iz; worldToGrid(wx,wz,ix,iz); return co2[iz*N+ix]; }
    float getPH(float wx, float wz) const { int ix,iz; worldToGrid(wx,wz,ix,iz); return pH_field[iz*N+ix]; }
    float getDrug(float wx, float wz) const { int ix,iz; worldToGrid(wx,wz,ix,iz); return drug[iz*N+ix]; }

    void consumeDrug(float wx, float wz, float amount) {
        int ix,iz; worldToGrid(wx,wz,ix,iz);
        drug[iz*N+ix] = fmaxf(0, drug[iz*N+ix] - amount);
    }

    // Apply drug uniformly across entire plate
    void applyDrugUniform(float concentration) {
        for (int i = 0; i < N*N; i++) drug[i] = concentration;
    }

    // Inject drug at a world position (Gaussian blob)
    void injectDrug(float wx, float wz, float concentration, float radius) {
        for (int iz = 0; iz < N; iz++)
            for (int ix = 0; ix < N; ix++) {
                float gx = (float)ix/(N-1)*2*SCENE_BOUND - SCENE_BOUND;
                float gz = (float)iz/(N-1)*2*SCENE_BOUND - SCENE_BOUND;
                float dx = gx-wx, dz = gz-wz;
                float dist2 = dx*dx+dz*dz;
                float sigma2 = radius*radius;
                drug[iz*N+ix] += concentration * expf(-dist2/(2*sigma2));
            }
    }

    // Wash out all drug
    void washOut() {
        for (int i = 0; i < N*N; i++) drug[i] = 0;
    }

    void consume(float wx, float wz, float o2Rate, float gluRate, bool glycolytic) {
        int ix,iz; worldToGrid(wx,wz,ix,iz);
        int idx = iz*N+ix;
        o2[idx] = fmaxf(0, o2[idx] - o2Rate);
        glucose[idx] = fmaxf(0, glucose[idx] - gluRate);
        float aerobicFrac = glycolytic ? 0.2f : 1.0f;
        co2[idx] = fminf(1.0f, co2[idx] + CO2_PRODUCE_BASE * aerobicFrac);
        float lactateProd = glycolytic ? LACTATE_PRODUCE*2.0f : LACTATE_PRODUCE*0.3f;
        // Henderson-Hasselbalch with bicarbonate buffer (beta ≈ 25 mM/pH-unit in blood).
        // Acid load δ[H+] is proportional to lactateProd; buffering factor 0.4
        // dampens the raw pH drop. Normalised pH here stays in [0.55, 0.80].
        const float BUFFER_STRENGTH = 0.40f;
        pH_field[idx] = fmaxf(0.55f, pH_field[idx] - lactateProd * BUFFER_STRENGTH);
    }

    void diffuse(float dt, float envO2, float envGlu) {
        float edgeCO2 = 0.03f, edgePH = 0.74f;
        for (int iz = 0; iz < N; iz++)
            for (int ix = 0; ix < N; ix++) {
                int idx = iz*N+ix;
                bool onEdge = ix==0||ix==N-1||iz==0||iz==N-1;
                if (onEdge) {
                    o2Next[idx]  = envO2*0.55f + o2[idx]*0.45f;
                    gluNext[idx] = envGlu*0.55f + glucose[idx]*0.45f;
                    co2Next[idx] = edgeCO2 + (co2[idx]-edgeCO2)*0.55f;
                    pHNext[idx]  = edgePH + (pH_field[idx]-edgePH)*0.55f;
                    drugNext[idx] = drug[idx] * 0.90f; // drug washes out at edges
                    continue;
                }
                auto lap = [&](const float* f, int i) {
                    return f[i-1]+f[i+1]+f[i-N]+f[i+N]-4*f[i];
                };
                o2Next[idx]  = clampf((o2[idx]+DIFF_O2_COEFF*lap(o2,idx)*dt)*0.9998f, 0, 1);
                gluNext[idx] = clampf((glucose[idx]+DIFF_GLC_COEFF*lap(glucose,idx)*dt)*0.9998f, 0, 1);
                co2Next[idx] = clampf((co2[idx]+DIFF_CO2_COEFF*lap(co2,idx)*dt)*0.9998f, 0, 1);
                pHNext[idx]  = clampf(pH_field[idx]+DIFF_PH_COEFF*lap(pH_field,idx)*dt, 0.55f, 0.80f);
                // Drug diffusion + decay
                if (drugDiffCoeff > 0) {
                    float lapD = drug[idx-1]+drug[idx+1]+drug[idx-N]+drug[idx+N]-4*drug[idx];
                    drugNext[idx] = fmaxf(0, (drug[idx]+drugDiffCoeff*lapD*dt) * (1.0f-drugDecayRate));
                } else {
                    drugNext[idx] = drug[idx];
                }
            }
        std::copy(o2Next, o2Next+N*N, o2);
        std::copy(gluNext, gluNext+N*N, glucose);
        std::copy(co2Next, co2Next+N*N, co2);
        std::copy(pHNext, pHNext+N*N, pH_field);
        if (drugDiffCoeff > 0) std::copy(drugNext, drugNext+N*N, drug);
    }
};
#endif

// ── CDK/Cyclin ODE (Novak-Tyson 7-variable model) ──────────────────────
struct CDKState {
    float CycD=0.05f, Rb=0.90f, E2F=0.02f;
    float CycE=0.01f, CycA=0.01f, CycB=0.01f, p21=0.05f;

    // Unsynchronized-population init: an asynchronously-seeded HeLa
    // culture has ~46% G1, ~33% S, ~17% G2, ~4% M (Casciari 1992).
    // Without this, every cell in a fresh colony starts at CycE≈0.01
    // and takes 30+ bio-h to build up to the G1/S threshold together —
    // producing a big synchronized first-division burst at ~50 bio-h
    // that doesn't match real HeLa growth curves.
    void randomize() {
        float u = randf();
        if (u < 0.46f) {
            // G1 — anywhere along G1's CycD/Rb trajectory.
            float g1u = randf();   // 0 = fresh G1, 1 = about to enter S
            CycD = 0.10f + g1u * 0.50f;
            Rb   = 0.95f - g1u * 0.40f;      // Rb drops 0.95 → 0.55 through G1
            E2F  = 0.02f + g1u * 0.15f;
            CycE = 0.02f + g1u * 0.18f;
            CycA = 0.01f + g1u * 0.05f;
            CycB = 0.01f + randf() * 0.02f;
        } else if (u < 0.79f) {
            // S — CycE high, CycA ramping.
            float su = randf();
            CycD = 0.25f + randf() * 0.20f;
            Rb   = 0.35f - su * 0.15f;
            E2F  = 0.40f + randf() * 0.30f;
            CycE = 0.20f + randf() * 0.25f;
            CycA = 0.15f + su * 0.25f;
            CycB = 0.02f + randf() * 0.05f;
        } else if (u < 0.96f) {
            // G2 — CycA high, CycB rising.
            float g2u = randf();
            CycD = 0.20f + randf() * 0.20f;
            Rb   = 0.20f;
            E2F  = 0.55f;
            CycE = 0.15f + randf() * 0.10f;
            CycA = 0.40f + randf() * 0.15f;
            CycB = 0.08f + g2u * 0.15f;
        } else {
            // M — CycB above threshold, mitosis imminent.
            CycD = 0.30f;
            Rb   = 0.15f;
            E2F  = 0.50f;
            CycE = 0.10f;
            CycA = 0.55f;
            CycB = 0.28f + randf() * 0.10f;
        }
        p21 = 0.05f + randf() * 0.05f;
    }

    void step(float dt_bio, float growthSignal) {
        float gs = growthSignal;
        float p21e = fminf(1.0f, p21*1.5f);

        // CycD synthesis boosted for faster G1 transit (target: G1=46% of cycle)
        float dCycD = (gs*0.45f - CycD*(0.20f+(1-gs)*0.10f)) * dt_bio;
        float CDK4act = CycD*(1-p21e*0.5f);
        float CDK2Eact = CycE*(1-p21e);
        float dRb = (0.08f*(1-Rb) - (CDK4act*0.60f+CDK2Eact*0.40f)*Rb) * dt_bio;
        float RbP = 1-Rb;
        float dE2F = (RbP*0.50f*(1+E2F*1.2f) - (Rb*0.40f+0.10f)*E2F) * dt_bio;
        float dCycE = (E2F*0.45f*(1-CycA) - CycE*(0.15f+CycA*0.55f)) * dt_bio * (1-p21e*0.8f);
        float APC_Cdh1 = fmaxf(0, 1-(CDK2Eact+CycA)*1.2f);
        float APC_Cdc20 = fmaxf(0, CycB-0.25f)*2.0f; // lower APC threshold
        float dCycA = (E2F*0.40f*fmaxf(0,CycE-0.12f) - CycA*(0.05f+APC_Cdh1*0.35f+APC_Cdc20)) * dt_bio * (1-p21e*0.6f);
        // Cdc25 + CycB synthesis tuned so CycB rises past 0.25 in ~4 bio-h
        // of G2 (matches HeLa G2 duration). Rates were 4× too slow vs the
        // sdt scale change so CycB never crossed the M-entry threshold.
        float Cdc25 = fmaxf(0, CycA-0.20f)*2.0f;
        float dCycB = (0.55f*Cdc25*(1+CycB*0.8f) - CycB*(0.04f+APC_Cdc20*2.5f)) * dt_bio * (1-p21e*0.3f);
        float dp21 = -p21*0.04f*dt_bio;

        CycD=clampf(CycD+dCycD,0,1.5f); Rb=clampf(Rb+dRb,0,1);
        E2F=clampf(E2F+dE2F,0,1); CycE=clampf(CycE+dCycE,0,1.5f);
        CycA=clampf(CycA+dCycA,0,1.5f); CycB=clampf(CycB+dCycB,0,1.5f);
        p21=clampf(p21+dp21,0,1);
    }

    // Phase thresholds calibrated for G1:46% S:33% G2:17% M:4%
    int getPhase() const {
        if (CycB>0.25f) return 3; // M  (was 0.30 — never reached)
        if (CycA>0.30f) return 2; // G2 (was 0.40)
        if (CycA>0.08f||CycE>0.18f) return 1; // S (was 0.10/0.25)
        return 0; // G1
    }
    bool readyForS() const { return CycE>0.18f&&Rb<0.50f&&p21<0.50f; }
    bool readyForM() const { return CycB>0.25f&&p21<0.35f&&CycA>0.30f; }
    void resetForNewCycle(float gs) {
        CycD=0.08f*gs; Rb=0.90f; E2F=0.03f;
        CycE=0.01f; CycA=0.01f; CycB=0.01f;
        p21=fmaxf(0.02f, p21*0.2f);
    }
};

// ── Cell state ──────────────────────────────────────────────────────────
static constexpr int SIM_FATE_UNDETERMINED = 0;
static constexpr int SIM_FATE_PROLIF       = 1;
static constexpr int SIM_FATE_QUIESCENT    = 2;
static constexpr int SIM_FATE_APOPTOTIC    = 3;

struct SimCell {
    simd_float3 position, velocity;
    float radius, size;
    bool alive;

    float ATP, stress, ROS, biomass, damageLevel, age;
    float mitoPotential, mitoHealth;
    bool glycolytic; float warburgTimer;
    float hypoxiaTimer, hypoxiaIntensity; bool necrotic;

    CDKState cdk; int phase;
    float cycleTimer, cycleProgress;
    bool checkpointG1Passed, checkpointG2Passed;
    float g1WaitTimer, g2WaitTimer;

    int fate; float fateScores[3]; float fateTimer; bool fateLocked;
    float atpDangerTimer;

    float glycolysisBias, prolifBias, rosTolerance, repairRate;
    int generation, cloneId; float telomere; bool senescent;
    int cellUid = 0; // Unique per cell instance — never shared, never reused

    bool divisionPending; float localPressure;
    float motileAngle, motileSpeed;
    int apoptosisPhase; float apoTimer;
    float adaptationTimer;
    float divisionCooldown = 0; // Prevents immediate re-entry into M-phase after division
    float postDivisionRecovery = 0; // Seconds of post-mitotic settling / render blend

    // Drug response (PhysiPKPD model)
    float drugInternal;     // internalized drug concentration
    float drugDamage;       // accumulated drug-induced damage (0-1)
    bool  drugResistant;    // MDR resistance mutation

    // ── Substrate adhesion (integrin / focal adhesion maturation) ────
    // Drives Z-position drift, motility damping, and a "spread" deform.
    // Newborn / rounded-up cells start at adhesionStrength = 0 and walk
    // through the Adhesion::* timeline as they sit on the dish floor.
    // Cells re-entering mitosis round up: spreadFactor decays toward 1.0
    // and adhesionStrength decays toward 0 over Adhesion::ROUND_UP_BIOSEC.
    float adhesionStrength = 0.0f;       // 0 (floating) → 0.95 (mature FA cluster)
    float adhesionTimer    = 0.0f;       // bio-seconds since attachment began
    int   focalAdhesionCount = 0;        // 0 → ~50 dots rendered under the cell
    float spreadFactor     = 1.0f;       // 1.0 spheroid → 1.30 spread disk

    // ── Osmosis / water flux / turgor pressure ───────────────────────
    // Water moves across the membrane down its osmolarity gradient via
    // aquaporins (Verkman 2008). Net flux drives elastic swelling /
    // shrinking that the user sees as a slight cell-radius change. The
    // real biology is per-bio-second; visual effect ranges ±20 %.
    float osmoCytoMM   = Osmosis::OSMO_ISOTONIC_MM; // cytosolic osmolarity (mM)
    float waterMM      = MediumComposition::DMEM_WATER_MM; // intracellular [H2O]
    float turgorPa     = 0.0f;                      // hydrostatic pressure above ext
    int   aquaporinCount = Osmosis::AQUAPORIN_COUNT_DEFAULT; // membrane mosaic count

    // ── Mitochondrial network dynamics ─────────────────────────────
    // HeLa cells carry 383-882 mitochondria per cell (Posakony 1977
    // JCB — serial-section EM morphometry on synchronized cells;
    // count approximately doubles from G1 → G2 and partitions at
    // division). Mammalian somatic range is 80-2000 (Alberts Ch 14);
    // hepatocytes 1000-2000, fibroblasts 100-500, muscle/neurons >1000.
    // Dynamics:
    //   * DRP1 fission — recruited at ER-mito contact sites by MFF /
    //     MiD49/51; rises in stress, G2/M, low ATP (Friedman 2011 Sci).
    //   * MFN1/2 + OPA1 fusion — requires GTP; dominant in G1/S with
    //     high ATP and low stress (Chen 2003 JCB).
    //   * PINK1/Parkin mitophagy — activated by ΔΨm loss; removes
    //     damaged mitos (Pickrell & Youle 2015 Neuron).
    //   * PGC-1α biogenesis — AMPK → SIRT1 → PGC-1α → NRF1/2 + TFAM
    //     (Scarpulla 2011 BBA); driven by ATP deficit.
    int   mitoCount       = 500;   // HeLa median count (Posakony 1977)
    float mitoHealthFrac  = 1.0f;  // fraction passing ΔΨm threshold (0-1)
    float mitoNetworkedFrac = 0.5f;// fraction in the fused tubular network
    float drp1Active      = 0.20f; // fission regulator
    float mfnActive       = 0.60f; // fusion regulator (MFN1/2 + OPA1 lumped)
    float pink1Active     = 0.05f; // mitophagy signal
    float pgc1aActive     = 0.30f; // biogenesis driver
    static constexpr int MITO_TARGET_DEFAULT = 500;  // HeLa median
    static constexpr int MITO_MIN = 200;             // fibroblast low end
    static constexpr int MITO_MAX = 4000;            // high-demand tissues

    // ── Biologically REAL molecule and macromolecule counts ─────────
    // HeLa cell volume ≈ 4 pL (4×10⁻¹² L). At 1 mM the cell holds
    // ~2.4×10⁹ molecules (1e-3 × 4e-15 × Nₐ). These counts drive
    // biology (protein synthesis rate scales with ribosome count,
    // ATP pool with molecule count, etc.) while particle rendering
    // shows a representative sample, not one sphere per molecule.
    //
    // Source anchors (Milo 2013 Bioessays; book.bionumbers.org;
    // hela_reference.csv for intracellular mM):
    //   Ribosomes     4–10 × 10⁶  (Wolf & Schlessinger 1977; nucleolus output)
    //   ATP           7.2 × 10⁹   (3 mM × 4 pL × Nₐ)
    //   Glucose       1.2 × 10¹⁰  (5 mM)
    //   NADH          2.4 × 10⁸   (0.1 mM)
    //   Amino acids   2.4 × 10¹⁰  (10 mM total free pool)
    //   Water         1.3 × 10¹⁴  (55 M — cell 70% water by mass)
    //   tRNA          1.0 × 10⁷   (Dittmar 2006 RNA)
    //   mRNA          3.6 × 10⁵   (Schwanhäusser 2011 Nature)
    //   Ca²⁺ free     2.4 × 10⁵   (100 nM cytosolic)
    //
    // We use double so we can hold counts up to ~10¹⁵ without loss.
    //
    // ── Genome ──────────────────────────────────────────────────────
    double genomeBp          = 6.4e9;   // human diploid, 3.2 Gb × 2 (post-S)
    double nucleosomes       = 3.0e7;   // 30 million nucleosomes/cell
                                        // (Luger 1997; ~1 per 200 bp)
    // ── Central-dogma machinery ────────────────────────────────────
    double ribosomeCount     = 6.0e6;   // Wolf & Schlessinger 1977
    double rnaPolII          = 4.0e4;   // Jackson 1998 / Osheim 2002
    double rnaPolI           = 2.0e2;   // nucleolar rRNA (~200 per cell)
    double rnaPolIII         = 1.0e3;   // tRNA/5S transcription
    double spliceosomes      = 1.0e5;   // Shepard 2009
    double nuclearPores      = 4.0e3;   // Ribbeck & Görlich 2001
    double chaperones        = 1.0e6;   // HSP70 + HSP90 + BiP ≈ 1% protein
    double proteasomes       = 1.0e6;   // Eytan 1989; ~1% protein pool
    double copiiVesicles     = 1.0e3;   // Bonfanti 1998 Cell
    double secretoryVesicles = 3.0e3;   // varies; HeLa moderate secretor
    double lysosomes         = 3.0e2;   // Saftig & Klumperman 2009
    double peroxisomes       = 3.0e2;   // Smith & Aitchison 2013
    double lipidDroplets     = 1.0e2;   // HeLa basal ~100 (Fujimoto 2014)
    // ── Cytoskeleton ────────────────────────────────────────────────
    double microtubulesTotal = 3.0e2;   // ~300 interphase MTs per HeLa
    double actinFilaments    = 5.0e4;   // Pollard 2000 estimates
    // ── Replication state ───────────────────────────────────────────
    double replicationForks  = 0.0;     // only nonzero during S-phase
    // ── Metabolites ────────────────────────────────────────────────
    double atpMolecules      = 7.2e9;   // 3 mM × 4 pL × Nₐ
    double glucoseMolecules  = 1.2e10;  // 5 mM
    double nadhMolecules     = 2.4e8;   // 0.1 mM
    double aaMolecules       = 2.4e10;  // total free pool, 10 mM
    double waterMolecules    = 1.3e14;  // 55 M
    double tRNACount         = 1.0e7;   // Dittmar 2006
    double mRNACount         = 3.6e5;   // Schwanhäusser 2011 Nature
    double calciumFree       = 2.4e5;   // 100 nM cytosolic
    double cAMPMolecules     = 2.4e4;   // 10 nM basal (Beavo & Brunton 2002)
    // ── Total protein pool (drives biomass) ────────────────────────
    double totalProteins     = 4.0e9;   // Milo 2013 Bioessays (2-4 × 10⁹)

    // Per-cell biological machinery: replication + mitosis state.
    // Before this existed, only cells[0] had a real replication program
    // (single global gCDogma) and a real mitosis state machine (single
    // global gMitosis). Every cell now owns its own.
    CellCycleProgram program;

    void init(simd_float3 pos, int idx) {
        position=pos; velocity={0,0,0};
        radius=CELL_RADIUS_BASE+randf()*CELL_RADIUS_VAR;
        size=0.85f+randf()*0.15f; alive=true;
        ATP=75+randf()*15; stress=5; ROS=0;
        biomass=1.0f+randf()*0.1f; damageLevel=randf()*0.05f; age=randf()*30;
        mitoPotential=170+randf()*10; mitoHealth=1.0f;
        glycolytic=false; warburgTimer=0;
        hypoxiaTimer=0; hypoxiaIntensity=0; necrotic=false;
        cdk.randomize();
        // Align cycleTimer + checkpoints with the CDK phase so a cell
        // initialized (via randomize()) in S / G2 / M actually progresses
        // through its cycle. Without this, a cell with CycB > 0.25 reads
        // phase=M from the CDK but its cycleTimer is 0 and its G1/G2
        // checkpoints are false — so it sits in M forever, never dividing.
        phase = cdk.getPhase();
        const float g1end = CYCLE_G1_DUR;
        const float send  = g1end + CYCLE_S_DUR;
        const float g2end = send  + CYCLE_G2_DUR;
        if (phase == 0) {
            cycleTimer           = randf() * g1end;
            checkpointG1Passed   = false;
            checkpointG2Passed   = false;
        } else if (phase == 1) {
            cycleTimer           = g1end + randf() * CYCLE_S_DUR;
            checkpointG1Passed   = true;
            checkpointG2Passed   = false;
        } else if (phase == 2) {
            cycleTimer           = send + randf() * CYCLE_G2_DUR;
            checkpointG1Passed   = true;
            checkpointG2Passed   = false;
        } else {  // M
            cycleTimer           = g2end;
            checkpointG1Passed   = true;
            checkpointG2Passed   = true;
        }
        cycleProgress = 0;
        g1WaitTimer = 0; g2WaitTimer = 0;
        // For cells seeded in G2 / M, their DNA is already replicated so
        // the mitosis-entry gate can pass. Without this, M-phase cells
        // stall at startMitosisProgram() because dnaCheckpointPassed
        // never flips true (replicationProgress still 0).
        if (phase >= 2) {
            program.ensureCDogmaInitialized();
            program.cdogma.replicationProgress = 1.0f;
            program.cdogma.replicationQuality  = 1.0f;
            program.cdogma.chk1Signal          = 0.0f;
            program.cdogma.replicationActive   = false;
        }
        fate=SIM_FATE_UNDETERMINED; fateScores[0]=fateScores[1]=fateScores[2]=0;
        fateTimer=0; fateLocked=false; atpDangerTimer=0;
        glycolysisBias=0.8f+randf()*0.4f; prolifBias=0.8f+randf()*0.4f;
        rosTolerance=0.8f+randf()*0.4f; repairRate=0.8f+randf()*0.4f;
        generation=0; cloneId=idx; telomere=TELO_INIT_LENGTH; senescent=false;
        divisionPending=false; localPressure=0;
        motileAngle=randf()*M_PI*2; motileSpeed=MOTILITY_SPEED*(0.5f+randf());
        apoptosisPhase=0; apoTimer=0;
        adaptationTimer=0; divisionCooldown=0; postDivisionRecovery=0;
        drugInternal=0; drugDamage=0; drugResistant=false;
        // Newborn cells start floating (adhesion = 0) and progress through
        // the Adhesion timeline naturally via updateAdhesion().
        adhesionStrength = 0.0f;
        adhesionTimer    = 0.0f;
        focalAdhesionCount = 0;
        spreadFactor     = 1.0f;
        // Osmosis baseline — isotonic with DMEM at start (290 mOsm; cell
        // is in equilibrium with the dish until uptake / secretion shifts
        // either side of the membrane).
        osmoCytoMM   = Osmosis::OSMO_ISOTONIC_MM;
        waterMM      = MediumComposition::DMEM_WATER_MM;
        turgorPa     = 0.0f;
        aquaporinCount = Osmosis::AQUAPORIN_COUNT_DEFAULT;
        // Mitochondrial count at seeding — HeLa 500 ±150 matches
        // Posakony 1977 JCB range (383-882) and accounts for cell-cycle
        // phase: G1 low end, G2 high end.
        mitoCount = MITO_TARGET_DEFAULT + (int)((randf() - 0.5f) * 300.0f);
        if (mitoCount < MITO_MIN) mitoCount = MITO_MIN;
        mitoHealthFrac = 0.95f + randf() * 0.05f;
        mitoNetworkedFrac = 0.5f;
        drp1Active = 0.20f; mfnActive = 0.60f;
        pink1Active = 0.05f; pgc1aActive = 0.30f;

        // Biologically real per-cell counts — values are the literature
        // MEDIAN with symmetric jitter around the median so an average
        // cell initializes at the target count, not below it.
        // Sources: Milo 2013 Bioessays; BioNumbers; Alberts MBoC 7e.
        auto jitter = [](float pct) -> float { return 1.0f - pct + 2.0f * pct * randf(); };
        // Genome (diploid 2N; doubles to 4N in S-phase)
        genomeBp          = 6.4e9;                        // fixed human diploid
        nucleosomes       = 3.0e7  * jitter(0.10f);
        // Central-dogma machinery
        ribosomeCount     = 6.0e6  * jitter(0.30f);       // median 6M (Wolf 1977)
        rnaPolII          = 4.0e4  * jitter(0.30f);
        rnaPolI           = 2.0e2  * jitter(0.40f);
        rnaPolIII         = 1.0e3  * jitter(0.40f);
        spliceosomes      = 1.0e5  * jitter(0.25f);
        nuclearPores      = 4.0e3  * jitter(0.15f);
        chaperones        = 1.0e6  * jitter(0.25f);
        proteasomes       = 1.0e6  * jitter(0.25f);
        copiiVesicles     = 1.0e3  * jitter(0.40f);
        secretoryVesicles = 3.0e3  * jitter(0.40f);
        lysosomes         = 3.0e2  * jitter(0.40f);
        peroxisomes       = 3.0e2  * jitter(0.40f);
        lipidDroplets     = 1.0e2  * jitter(0.70f);       // highly variable
        // Cytoskeleton
        microtubulesTotal = 3.0e2  * jitter(0.30f);
        actinFilaments    = 5.0e4  * jitter(0.25f);
        replicationForks  = 0.0;                          // 0 outside S
        // Metabolites
        atpMolecules     = 7.2e9  * jitter(0.25f);
        glucoseMolecules = 1.2e10 * jitter(0.30f);
        nadhMolecules    = 2.4e8  * jitter(0.25f);
        aaMolecules      = 2.4e10 * jitter(0.25f);
        waterMolecules   = 1.3e14 * jitter(0.05f);
        tRNACount        = 1.0e7  * jitter(0.20f);
        mRNACount        = 3.6e5  * jitter(0.30f);
        calciumFree      = 2.4e5  * jitter(0.50f);
        cAMPMolecules    = 2.4e4  * jitter(0.70f);
        totalProteins    = 4.0e9  * jitter(0.30f);
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  Checkpoint predicates
//  -----------------------------------------------------------------------
//  Each predicate returns (passed, reason). On failure, `reason` is a
//  literal describing the FIRST predicate that failed, in the order
//  evaluated. The audit log writes the verbatim string so "why didn't
//  this cell divide" always has a citable answer.
//  -----------------------------------------------------------------------
struct CheckpointResult { bool passed; const char* reason; };

// G1 / S restriction-point gate. Per Bartek 2004, Geng 2003, Sherr 2002.
static CheckpointResult evalG1S(const SimCell& c,
                                float local_glucose_mM,
                                float local_glutamine_mM,
                                bool contactArrested) {
    if (c.cdk.Rb >= 0.40f)        return {false, "Rb still active"};
    if (c.cdk.CycE <= 0.18f)      return {false, "CycE below restriction-pt"};
    if (c.cdk.p21 >= 0.35f)       return {false, "p21 high (stress / damage)"};
    if (c.biomass <= 1.30f)       return {false, "insufficient biomass"};
    if (c.ATP <= 25.0f)           return {false, "ATP < 25"};
    if (local_glucose_mM <= 1.0f) return {false, "glucose-starved"};
    if (local_glutamine_mM <= 0.05f) return {false, "glutamine-starved"};
    if (c.damageLevel >= 0.25f)   return {false, "DNA damage above threshold"};
    if (contactArrested)          return {false, "contact-inhibited (Hippo)"};
    return {true, "G1S_OK"};
}

// G2 / M entry gate. Per Bartek 2004, Lindqvist 2009, Hartwell 1989.
// `cdogma` is the cell's per-cell CentralDogmaState (replication program).
static CheckpointResult evalG2M(const SimCell& c, const CentralDogmaState& cdogma) {
    if (cdogma.replicationProgress < 0.999f)
        return {false, "replication incomplete"};
    if (cdogma.countActiveReplicationForks() > 0)
        return {false, "forks still active"};
    if (cdogma.chk1Signal >= 0.35f)
        return {false, "Chk1 high (replication stress)"};
    if (cdogma.escapedErrors >= 5)
        return {false, "too many uncorrected errors"};
    if (c.cdk.CycB <= 0.25f)      return {false, "CycB below activation"};
    if (c.cdk.p21 >= 0.35f)       return {false, "p21 high (p53 / DDR)"};
    if (c.damageLevel >= 0.40f)   return {false, "DSBs above threshold"};
    if (c.ATP <= 30.0f)           return {false, "ATP insufficient for mitosis"};
    if (c.biomass <= 1.70f)       return {false, "biomass < 1.70"};
    return {true, "G2M_OK"};
}

// Spindle-Assembly Checkpoint release. Anaphase only when EVERY chromosome
// is amphitelically attached, under tension, and MCC has dissolved.
// Per Musacchio 2015, Cimini 2001.
static CheckpointResult evalSAC(const MitosisState& mitosis) {
    if (mitosis.mccLevel >= 0.10f)
        return {false, "MCC still inhibiting APC"};
    for (int i = 0; i < MitosisState::NUM_CHROMO; i++) {
        const ChromosomeState& ch = mitosis.chromosomes[i];
        if (ch.attachmentA < 0.85f || ch.attachmentB < 0.85f)
            return {false, "unattached kinetochore"};
        if (ch.attachmentType != ATTACH_AMPHITELIC)
            return {false, "merotelic / syntelic"};
        if (ch.tension < 0.30f)
            return {false, "low tension"};
    }
    return {true, "SAC_OK"};
}

// ══════════════════════════════════════════════════════════════════════════
class Simulation {
public:
    std::vector<SimCell> cells;
    // The dish chemistry. Closed system seeded from MediumComposition::DMEM_*.
    // Member name kept as `nutrients` so existing call sites compile through
    // the migration; future PRs can rename to `medium`.
    MediumField nutrients;
    float bioTime=0;
    // Legacy env sliders — kept so the (about-to-be-deleted) UI keeps
    // compiling. `MediumField::init(envO2, envGlu)` ignores them.
    float envO2=0.80f, envGlucose=0.70f;
    float timeScale=1.0f;
    float lastExecutedScaledDt=0.0f;
    float pendingScaledDt=0.0f;
    bool paused=false;
    SimMode mode=MODE_SINGLE_CELL;
    float primaryReplicationProgress = 1.0f;
    float primaryReplicationCheckpoint = 0.0f;
    float primaryReplicationQuality = 1.0f;
    bool primaryReplicationReady = true;

    int statAlive=0, statProlif=0, statQuiescent=0, statApoptotic=0, statNecrotic=0;
    float statAvgATP=0;
    int statDivisions=0, statDeaths=0;
    int statPhases[4]={};
    int statGlycolytic=0;
    int nextCloneId=0;
    int nextCellUid=0;
    int allocateCellUid() { return nextCellUid++; }

    static constexpr float MAX_SUBSTEP_DT = 0.05f;
    static constexpr float MAX_SCALED_DT_PER_FRAME_SINGLE = 2.00f;
    static constexpr float MAX_SCALED_DT_PER_FRAME_SINGLE_LINEAGE = 0.90f;
    static constexpr float MAX_SCALED_DT_PER_FRAME_COLONY = 0.20f;

    // Division events for cross-boundary interior splitting
    struct DivisionEvent {
        int parentCellUid;
        int daughterACellUid, daughterBCellUid;
        simd_float3 posA, posB;
    };
    std::vector<DivisionEvent> pendingDivisions;

    bool visualMitosisActive = false; // Set by main.mm each frame

    // Derive a daughter cell from an IMMUTABLE parent snapshot.
    // Both daughters should be derived from the same snapshot to prevent double resource loss.
    SimCell deriveDaughter(const SimCell& original, simd_float3 pos, simd_float3 vel,
                           float sizeFactor, const SplitStats& stats, bool isA) {
        SimCell d = original; // Full copy from immutable snapshot
        d.position = pos;
        d.velocity = vel;
        d.cellUid = allocateCellUid();
        bool singleMode = (mode == MODE_SINGLE_CELL);

        // RESET (cell cycle restarts):
        d.phase = 0;
        d.cdk.resetForNewCycle(singleMode ? 1.12f : 1.0f);
        d.cycleTimer = 0; d.cycleProgress = 0;
        d.checkpointG1Passed = false; d.checkpointG2Passed = false;
        d.g1WaitTimer = 0; d.g2WaitTimer = 0;
        d.divisionPending = false;
        d.divisionCooldown = singleMode ? 1.5f : 8.0f;
        // Reset the daughter's replication + mitosis program — without this
        // the daughter inherits the parent's post-S-phase state (replication
        // progress ≈ 1.0, all origins fired, forks drained). That would let
        // it walk straight from G1 into G2/M on the next cycle with no real
        // replication, breaking the central invariant "no cell enters M
        // until its OWN genome has been re-replicated."
        d.program.resetForNewCycle();
        // Leave cdogmaInitialized as-is (inherited from parent via the
        // copy); the reset clears the replication fields but keeps the
        // codingMRNA/transcription scaffolding that init() built.
        // This timer is interpreted in real wall-clock seconds, not scaled sim
        // time, so the daughters remain visually continuous even at 20x/70x.
        d.postDivisionRecovery = singleMode ? 6.0f : 4.0f;
        d.fateLocked = false;
        d.fate = SIM_FATE_PROLIF;
        d.fateScores[0] = 8.0f;
        d.fateScores[1] = 0.0f;
        d.fateScores[2] = 0.0f;
        d.fateTimer = 0.0f;
        d.localPressure = 0.0f;

        // LINEAGE-UPDATED:
        d.generation = original.generation + 1;
        d.telomere = original.telomere - TELO_LOSS_PER_DIV;
        if (d.telomere < TELO_CRITICAL) d.senescent = true;

        // PARTITIONED by weighted cytoplasmic ratio:
        float ratio = isA ? stats.cytoplasmicRatioA : (1.0f - stats.cytoplasmicRatioA);
        if (ratio <= 0.0f || ratio >= 1.0f) ratio = 0.5f;

        // ATP and similar fast pools equilibrate much faster than structural cargo.
        // In single-cell teaching mode, heavily compress asymmetry so daughters
        // don't emerge metabolically crippled and stall for an unnaturally long time.
        float fastPoolRatio = ratio;
        if (singleMode) {
            fastPoolRatio = 0.5f + (ratio - 0.5f) * 0.20f;
            fastPoolRatio = clampf(fastPoolRatio, 0.40f, 0.60f);
        }

        d.ATP = fmaxf(original.ATP * fastPoolRatio, 55.0f); // floor to avoid instant ATP-starvation
        d.ROS = original.ROS * ratio;
        d.stress = original.stress * (singleMode
            ? (0.5f + (ratio - 0.5f) * 0.40f)
            : ratio);
        d.drugInternal = original.drugInternal * ratio;
        d.drugDamage = original.drugDamage * ratio;
        // Mitochondria partition between daughters by the same cytoplasmic
        // ratio (they ride the cytoplasm at cytokinesis). Biogenesis will
        // restore each daughter toward the target over the next cycle.
        d.mitoCount = (int)fmaxf(
            (float)SimCell::MITO_MIN,
            (float)original.mitoCount * ratio);
        d.mitoHealthFrac = original.mitoHealthFrac;
        d.mitoNetworkedFrac = 0.3f; // newly divided cells start more fragmented
        d.drp1Active = 0.30f; // slightly elevated just after division
        d.mfnActive = 0.45f;
        d.pink1Active = original.pink1Active * 0.5f;
        d.pgc1aActive = 0.40f; // biogenesis elevated post-division

        // Real per-cell molecule/macromolecule counts are partitioned by
        // the same cytoplasmic ratio — this is what ACTUALLY happens in a
        // cell: whichever ribosomes/mRNAs/ATP were on the daughter-A side
        // of the cleavage furrow go to daughter A. No program enforcement
        // or re-spawn. Biogenesis/consumption from updateMetabolism then
        // regrows the count over the next cycle.
        // Genome: each daughter gets full 2N diploid (sister chromatids split).
        d.genomeBp          = 6.4e9;
        d.nucleosomes       = original.nucleosomes       * ratio;
        // Central-dogma machinery partitions with cytoplasm / nucleus split.
        d.ribosomeCount     = original.ribosomeCount     * ratio;
        d.rnaPolII          = original.rnaPolII          * ratio;
        d.rnaPolI           = original.rnaPolI           * ratio;
        d.rnaPolIII         = original.rnaPolIII         * ratio;
        d.spliceosomes      = original.spliceosomes      * ratio;
        d.nuclearPores      = original.nuclearPores      * ratio;
        d.chaperones        = original.chaperones        * ratio;
        d.proteasomes       = original.proteasomes       * ratio;
        d.copiiVesicles     = original.copiiVesicles     * ratio;
        d.secretoryVesicles = original.secretoryVesicles * ratio;
        d.lysosomes         = original.lysosomes         * ratio;
        d.peroxisomes       = original.peroxisomes       * ratio;
        d.lipidDroplets     = original.lipidDroplets     * ratio;
        d.microtubulesTotal = original.microtubulesTotal * ratio;
        d.actinFilaments    = original.actinFilaments    * ratio;
        d.replicationForks  = 0.0;                       // not in S yet
        d.atpMolecules      = original.atpMolecules      * ratio;
        d.glucoseMolecules  = original.glucoseMolecules  * ratio;
        d.nadhMolecules     = original.nadhMolecules     * ratio;
        d.aaMolecules       = original.aaMolecules       * ratio;
        d.waterMolecules    = original.waterMolecules    * ratio;
        d.tRNACount         = original.tRNACount         * ratio;
        d.mRNACount         = original.mRNACount         * ratio;
        d.calciumFree       = original.calciumFree       * ratio;
        d.cAMPMolecules     = original.cAMPMolecules     * ratio;
        d.totalProteins     = original.totalProteins     * ratio;

        // SIZE/BIOMASS: size is source of truth, biomass derived
        d.size = original.size * sizeFactor;
        d.biomass = biomassFromCellSizeFactor(d.size);

        // INHERITED-with-cleanup: damage is halved at division (daughters
        // get half each) rather than both inheriting full damage, which
        // would accumulate without bound and cause random cell death.
        d.damageLevel = original.damageLevel * 0.5f;
        // Reset transient danger/stress timers so a daughter is not
        // immediately flagged apoptotic from the parent's timer state.
        d.atpDangerTimer = 0;
        d.hypoxiaTimer = 0;
        d.apoptosisPhase = 0;
        d.apoTimer = 0;
        d.necrotic = false;

        // INHERITED as-is (mutations accumulate across generations):
        // glycolysisBias, prolifBias, rosTolerance, repairRate, drugResistant
        return d;
    }

    // Drug system
    int   activeDrugIdx = 0;       // index into DRUG_LIBRARY (0 = none)
    float drugConcentration = 0;   // applied concentration (µM)
    float statViability = 100.0f;  // % cells alive vs pre-drug count
    int   preDrugCount = 0;        // cell count when drug was applied
    float statAvgDrugDamage = 0;

    // Mitosis visualization flag (set by renderer when visual mitosis completes)
    bool mitosisVisualizationComplete = false;

    void applyDrugUniform(int drugIdx, float conc) {
        activeDrugIdx = drugIdx;
        drugConcentration = conc;
        preDrugCount = statAlive > 0 ? statAlive : 1;
        const Drug& d = DRUG_LIBRARY[drugIdx];
        nutrients.drugDiffCoeff = d.diffusionCoeff;
        nutrients.drugDecayRate = d.decayRate;
        nutrients.applyDrugUniform(conc);
    }

    void injectDrug(int drugIdx, float conc, float wx, float wz) {
        activeDrugIdx = drugIdx;
        drugConcentration = conc;
        if (preDrugCount == 0) preDrugCount = statAlive > 0 ? statAlive : 1;
        const Drug& d = DRUG_LIBRARY[drugIdx];
        nutrients.drugDiffCoeff = d.diffusionCoeff;
        nutrients.drugDecayRate = d.decayRate;
        nutrients.injectDrug(wx, wz, conc, 8.0f); // 8-unit radius Gaussian
    }

    void washOutDrug() {
        nutrients.washOut();
        nutrients.drugDiffCoeff = 0;
        nutrients.drugDecayRate = 0;
        for (auto& c : cells) { c.drugInternal = 0; c.drugDamage *= 0.5f; }
    }

    void init() { init(mode); }

    void init(SimMode m) {
        mode = m;
        cells.clear(); nutrients.init(envO2, envGlucose);
        bioTime = 0;
        lastExecutedScaledDt = 0;
        pendingScaledDt = 0;
        statDivisions = 0; statDeaths = 0;
        activeDrugIdx = 0; drugConcentration = 0;
        statViability = 100.0f; preDrugCount = 0; statAvgDrugDamage = 0;

        int count = (mode == MODE_SINGLE_CELL) ? SINGLE_CELL_COUNT : INIT_CELLS;
        nextCloneId = count;
        nextCellUid = 0;
        pendingDivisions.clear();
        for (int i = 0; i < count; i++) {
            float a = (float)i / fmaxf(1.0f, (float)count) * 2 * M_PI;
            float r = (mode == MODE_SINGLE_CELL) ? 0.0f : 5.0f + randf() * 3.0f;
            simd_float3 pos = {r * cosf(a), FLOOR_Y + 2.2f, r * sinf(a)};
            SimCell c; c.init(pos, i);
            c.cellUid = allocateCellUid();
            if (mode == MODE_SINGLE_CELL) {
                // Teaching / inspection mode should start from one healthy,
                // division-capable founder cell rather than a randomly aged
                // cell that can silently drift into apoptosis before any
                // lineage behavior is visible.
                c.age = 0.0f;
                c.ATP = fmaxf(c.ATP, 82.0f);
                c.stress = fminf(c.stress, 3.0f);
                c.ROS = fminf(c.ROS, 2.0f);
                c.damageLevel = fminf(c.damageLevel, 0.01f);
                c.fate = SIM_FATE_PROLIF;
                c.fateScores[0] = 8.0f;
                c.fateScores[1] = 0.0f;
                c.fateScores[2] = 0.0f;
                c.fateTimer = 0.0f;
                c.fateLocked = false;
                c.apoptosisPhase = 0;
                c.apoTimer = 0.0f;
                c.atpDangerTimer = 0.0f;
            }
            cells.push_back(c);
        }
    }

    void update(float dt) {
        lastExecutedScaledDt = 0.0f;
        if (paused||dt<=0) return;

        for (auto& c : cells) {
            c.postDivisionRecovery = fmaxf(0.0f, c.postDivisionRecovery - dt);
        }

        float requestedDt = dt * timeScale;
        pendingScaledDt += requestedDt;
        float responsiveCap = MAX_SCALED_DT_PER_FRAME_COLONY;
        if (mode == MODE_SINGLE_CELL) {
            responsiveCap = (cells.size() <= 2)
                ? MAX_SCALED_DT_PER_FRAME_SINGLE
                : MAX_SCALED_DT_PER_FRAME_SINGLE_LINEAGE;
        }
        float totalDt = fminf(pendingScaledDt, responsiveCap);
        if (totalDt <= 0.0f) return;
        pendingScaledDt -= totalDt;

        float subDt = fminf(totalDt, MAX_SUBSTEP_DT);
        int steps = (int)ceilf(totalDt / fmaxf(subDt, 1e-6f));
        subDt = totalDt / fmaxf((float)steps, 1.0f);

        lastExecutedScaledDt = totalDt;
        bioTime += totalDt*BIO_MIN_PER_SEC*60;
        for (int s=0; s<steps; s++) {
            updatePhysics(subDt);
            for (auto& c:cells) { if(!c.alive) continue;
                updateMetabolism(c,subDt); updateDrugResponse(c,subDt);
                updateMitochondria(c,subDt);
                updateAdhesion(c,subDt);
                updateCellCycle(c,subDt); updateMitosisProgram(c,subDt);
                // Per-cell replication / transcription tick. WITHOUT this
                // every cell is stuck at replicationProgress = 0 and the
                // G2/M checkpoint never opens. Dt clamped to a safe range
                // so replication is stable across timeScales.
                c.program.ensureCDogmaInitialized();
                float dogmaDt = fminf(fmaxf(subDt, 1.0f/240.0f), 1.0f/30.0f);
                c.program.cdogma.update(dogmaDt, c.phase);
                updateFate(c,subDt); updateApoptosis(c,subDt);
            }
            nutrients.diffuse(subDt, envO2, envGlucose);
        }
        processDivisions();
        cells.erase(std::remove_if(cells.begin(),cells.end(),
            [](const SimCell& c){return !c.alive;}), cells.end());
        updateStats();
    }

private:
    bool singleCellPostDivisionDisplay() const {
        return false; // All cells run full cycle — no freezing after division
    }

    bool isPrimaryCell(const SimCell& c) const {
        return !cells.empty() && &c == &cells[0];
    }

    bool usesPrimaryVisualMitosis(const SimCell& c) const {
        // Only the focused primary cell (cells[0]) uses the detailed visual
        // mitosis machinery in main.mm (gMitosis: condensed chromatin,
        // spindle, chromosome alignment, furrow). Every other cell — even
        // in single-cell teaching mode — runs its own per-cell mitosis
        // program and divides through processDivisions. Previously this
        // function returned true for ALL cells in single mode, which meant
        // background cells reached M-phase (DNA=1524) and then sat there
        // forever because neither the visual path nor the per-cell path
        // would finalize their division.
        return mode == MODE_SINGLE_CELL && isPrimaryCell(c);
    }

    void syncMitosisCheckpointInputs(SimCell& c) {
        c.program.ensureCDogmaInitialized();
        c.program.mitosis.dnaDuplicationProgress =
            clampf(c.program.cdogma.replicationProgress, 0.0f, 1.0f);
        c.program.mitosis.dnaCheckpointPassed =
            c.program.cdogma.replicationReadyForM();
        c.program.mitosis.replicationQuality =
            c.program.cdogma.replicationQuality;
        // Background cells do not yet own a full per-cell particle interior
        // inside Simulation itself, so their mitosis program uses the
        // replicated genome state as the biological duplication signal.
        c.program.mitosis.particlesDuplicated =
            c.program.mitosis.dnaDuplicationProgress >= 0.995f;
    }

    void startMitosisProgram(SimCell& c) {
        float cellR = c.radius * c.size;
        c.program.mitosis.start(c.position, cellR);
        syncMitosisCheckpointInputs(c);
        c.program.logEvent(bioTime, 3, "MITOSIS_START");
    }

    void updateMitosisProgram(SimCell& c, float dt) {
        if (!c.alive || c.apoptosisPhase > 0 || c.senescent) return;
        if (c.postDivisionRecovery > 0.0f) return;
        if (usesPrimaryVisualMitosis(c)) return;
        if (c.phase != 3 && !c.program.mitosis.active) return;

        syncMitosisCheckpointInputs(c);
        if (!c.program.mitosis.active) {
            if (!c.checkpointG2Passed || !c.program.mitosis.dnaCheckpointPassed) return;
            startMitosisProgram(c);
        }

        c.phase = 3;
        c.cycleProgress = 1.0f;

        if (!c.program.mitosis.postDivisionComplete()) {
            float cellR = c.radius * c.size;
            c.program.mitosis.update(dt, c.position, cellR);
            c.cycleTimer = fmaxf(c.cycleTimer,
                                 CYCLE_G1_DUR + CYCLE_S_DUR + CYCLE_G2_DUR +
                                 fminf(c.program.mitosis.totalProgress, CYCLE_M_DUR));
            if (c.program.mitosis.postDivisionComplete()) {
                c.divisionPending = true;
                c.program.logEvent(bioTime, 3, "CYTOKINESIS_OK");
            }
        } else {
            c.divisionPending = true;
            c.program.mitosis.postDivisionTimer += dt;
        }
    }

    void updatePhysics(float dt) {
        int n=(int)cells.size();
        for (auto& c:cells) c.localPressure=0;
        for (int i=0;i<n;i++) {
            auto& a=cells[i]; if(!a.alive) continue;
            float recoveryT = clampf(a.postDivisionRecovery / 6.0f, 0.0f, 1.0f);
            float motilityBlend = 1.0f - recoveryT * 0.82f;
            a.motileAngle += (randf()-0.5f)*0.3f*dt*(0.35f + 0.65f*motilityBlend);
            float spd=a.motileSpeed*(a.fate==SIM_FATE_PROLIF?1.5f:0.5f)*motilityBlend;
            if (a.necrotic) spd*=0.1f;
            a.velocity.x += cosf(a.motileAngle)*spd*dt*0.5f;
            a.velocity.z += sinf(a.motileAngle)*spd*dt*0.5f;
            for (int j=i+1;j<n;j++) {
                auto& b=cells[j]; if(!b.alive) continue;
                float dx=a.position.x-b.position.x, dz=a.position.z-b.position.z;
                float dist=sqrtf(dx*dx+dz*dz);
                // Use the true rendered radii for cell-cell exclusion so
                // post-mitosis daughters don't visually overlap and hide one
                // another's nucleus/DNA.
                float minD=(a.radius*a.size + b.radius*b.size)*0.98f;
                if (dist<minD&&dist>0.01f) {
                    float overlap=minD-dist;
                    float force=HERTZ_STIFFNESS*powf(overlap,1.5f);
                    float nx=dx/dist, nz=dz/dist;
                    a.velocity.x+=nx*force*dt; a.velocity.z+=nz*force*dt;
                    b.velocity.x-=nx*force*dt; b.velocity.z-=nz*force*dt;
                    a.localPressure+=overlap; b.localPressure+=overlap;
                }
            }
            a.velocity.x*=powf(CELL_DAMPING,dt*60);
            a.velocity.z*=powf(CELL_DAMPING,dt*60);
            a.position.x+=a.velocity.x*dt; a.position.z+=a.velocity.z*dt;
            float dfc=sqrtf(a.position.x*a.position.x+a.position.z*a.position.z);
            if (dfc>SCENE_BOUND-a.radius) {
                float nx=a.position.x/dfc, nz=a.position.z/dfc;
                a.position.x=nx*(SCENE_BOUND-a.radius);
                a.position.z=nz*(SCENE_BOUND-a.radius);
                float vn=a.velocity.x*nx+a.velocity.z*nz;
                if(vn>0){a.velocity.x-=nx*vn*1.5f; a.velocity.z-=nz*vn*1.5f;}
            }
            // Substrate anchor: spread cells sit slightly lower (flatter)
            // and rounded cells sit higher. baseY shifts with spreadFactor
            // so the visual matches the adhesion state.
            float spread = a.spreadFactor;
            float baseY = FLOOR_Y + a.radius * a.size * 0.85f / spread;
            if (a.postDivisionRecovery > 0.0f) {
                float bobY = baseY + sinf(bioTime*0.0001f+(float)i*0.7f)*0.12f;
                float settleBlend = 1.0f - recoveryT;
                a.velocity.y += (-0.6f - 1.2f * recoveryT) * dt;
                a.velocity.y *= powf(0.92f, dt * 60.0f);
                a.position.y += a.velocity.y * dt;
                if (a.position.y < baseY) {
                    a.position.y = baseY;
                    if (a.velocity.y < 0.0f) a.velocity.y *= -0.18f;
                }
                a.position.y += (bobY - a.position.y) * settleBlend * 0.12f;
            } else {
                // Strong elastic anchor toward baseY scaled by adhesion
                // strength (mature focal adhesions are stiff; floating
                // cells just flop around). Per Geiger 2001 / Gallant 2005.
                float anchorY = baseY + sinf(bioTime*0.0001f+(float)i*0.7f)*0.12f;
                float k = 0.30f + 4.0f * a.adhesionStrength;     // anchor stiffness
                a.position.y += (anchorY - a.position.y) * fminf(1.0f, k * dt);
                // Adhered cells also damp lateral velocity (FA grip).
                float damp = 1.0f - a.adhesionStrength * 0.40f;
                a.velocity.x *= damp;
                a.velocity.z *= damp;
            }
        }
    }

    // ── Mitochondrial network dynamics ──────────────────────────────
    // Models steady-state mitochondrial count via a balance of fission,
    // fusion, mitophagy, and biogenesis. Not an exact reproduction — the
    // regulators (DRP1, MFN1/2, OPA1, PINK1, PGC-1α) are lumped rather
    // than individually integrated — but captures the biology:
    //
    //  * stress + G2/M phase → DRP1 active → fission (count ↑)
    //  * healthy G1 with high ATP → fusion (count ↓, more networked)
    //  * membrane-potential loss + damage → PINK1 → mitophagy (count ↓)
    //  * low ATP / AMPK activation → PGC-1α → biogenesis (count ↑)
    //
    // Sources: Mishra & Chan 2014 Nat Rev MCB (dynamics overview);
    //          Friedman 2011 Science (ER-mito contact → DRP1);
    //          Chen 2003 JCB (MFN1/2); Pickrell & Youle 2015 Neuron;
    //          Scarpulla 2011 BBA (PGC-1α biogenesis).
    // ── Substrate adhesion (integrin-driven focal-adhesion maturation) ──
    // Drives the cell along Cavalcanti-Adam 2007 / Geiger 2001 timeline:
    // floating → contact → spreading → mature, with FA cluster growing
    // 0 → ~50 punctate dots and spreadFactor 1.00 → 1.30. The maturation
    // *decays* during mitosis when the cell rounds up (Maddox 2003) so
    // daughters re-spread on the dish naturally.
    void updateAdhesion(SimCell& c, float dt) {
        if (!c.alive) return;
        float bio_dt = bioDt(dt, timeScale);

        // Round-up: cells in late G2 and M lose adhesion (mitotic round-up,
        // Maddox 2003). Apoptotic / necrotic cells also detach.
        bool roundedUp = (c.phase == 3) || c.necrotic || c.apoptosisPhase > 0;
        if (roundedUp) {
            float decay = bio_dt / Adhesion::ROUND_UP_BIOSEC;
            c.adhesionStrength = fmaxf(0.0f, c.adhesionStrength - decay);
            c.spreadFactor     = fmaxf(Adhesion::SPREAD_FACTOR_FLOATING,
                                       c.spreadFactor - decay * 0.30f);
            c.adhesionTimer    = 0.0f;
            c.focalAdhesionCount = (int)(Adhesion::FA_COUNT_MATURE * c.adhesionStrength);
            return;
        }

        // Maturation timeline. spreadFactor and adhesionStrength interpolate
        // along the four-stage Cavalcanti-Adam / Geiger schedule.
        c.adhesionTimer += bio_dt;
        float t = c.adhesionTimer;
        if (t < Adhesion::FLOATING_BIOSEC) {
            c.spreadFactor     = Adhesion::SPREAD_FACTOR_FLOATING;
            c.adhesionStrength = 0.0f;
        } else if (t < Adhesion::CONTACT_BIOSEC) {
            float u = (t - Adhesion::FLOATING_BIOSEC)
                      / (Adhesion::CONTACT_BIOSEC - Adhesion::FLOATING_BIOSEC);
            c.spreadFactor     = 1.00f + 0.05f * u;
            c.adhesionStrength = 0.15f * u;
        } else if (t < Adhesion::SPREADING_BIOSEC) {
            float u = (t - Adhesion::CONTACT_BIOSEC)
                      / (Adhesion::SPREADING_BIOSEC - Adhesion::CONTACT_BIOSEC);
            c.spreadFactor     = 1.05f + (Adhesion::SPREAD_FACTOR_MATURE - 1.05f) * u;
            c.adhesionStrength = 0.15f + (0.65f - 0.15f) * u;
        } else {
            c.spreadFactor     = Adhesion::SPREAD_FACTOR_MATURE;
            c.adhesionStrength = fminf(Adhesion::STRENGTH_MATURE,
                                       0.65f + (t - Adhesion::SPREADING_BIOSEC) / 3600.0f * 0.30f);
        }
        c.focalAdhesionCount = (int)(Adhesion::FA_COUNT_MATURE * c.adhesionStrength);
    }

    void updateMitochondria(SimCell& c, float dt) {
        if (!c.alive) return;
        float mdt = dt * MEDIUM_DT_SCALE;

        // ── Signalling inputs (0-1 each) ─────────────────────────────
        auto sigmoid = [](float x) { return 1.0f / (1.0f + expf(-x)); };

        // DRP1 rises with stress + G2/M phase + low ATP (high ADP).
        float phaseG2M = (c.phase >= 2) ? 0.45f : 0.0f;
        float atpDeficit = clampf((70.0f - c.ATP) / 40.0f, 0.0f, 1.0f);
        c.drp1Active = sigmoid(c.stress / 60.0f + phaseG2M + atpDeficit * 0.6f - 1.2f);

        // MFN / OPA1 rise in G1 / S with low stress and sufficient ATP.
        float phaseG1S = (c.phase <= 1) ? 0.35f : 0.0f;
        c.mfnActive = sigmoid(-c.stress / 80.0f + phaseG1S + (c.ATP - 40.0f) / 50.0f);

        // PINK1 stabilized when mito health drops (proxy: mitoPotential
        // below threshold and/or damage accumulating).
        float healthDeficit = 1.0f - clampf(c.mitoHealthFrac, 0.0f, 1.0f);
        c.pink1Active = sigmoid(c.damageLevel * 2.0f + healthDeficit * 1.5f + c.ROS / 8.0f - 1.5f);

        // PGC-1α biogenesis signal — activated by AMPK (ATP low) and
        // in response to mitochondrial depletion.
        float countDeficit = clampf((float)(SimCell::MITO_TARGET_DEFAULT - c.mitoCount)
                                    / (float)SimCell::MITO_TARGET_DEFAULT, 0.0f, 1.0f);
        c.pgc1aActive = sigmoid(atpDeficit * 1.8f + countDeficit * 1.2f - 0.9f);

        // ── Event rates per second ───────────────────────────────────
        // Rates scale with current count (fission/fusion/mitophagy are
        // per-mitochondrion events). Biogenesis is a rate toward target.
        float k_fission    = 0.030f * c.drp1Active;
        float k_fusion     = 0.025f * c.mfnActive * c.mfnActive; // fusion needs 2
        float k_mitophagy  = 0.010f * c.pink1Active * healthDeficit;
        float k_biogenesis = 0.015f * c.pgc1aActive;

        float N = (float)c.mitoCount;
        // dN/dt = (+fission - fusion - mitophagy) * N + biogenesis * target
        float dN = (k_fission - k_fusion - k_mitophagy) * N
                 + k_biogenesis * (float)SimCell::MITO_TARGET_DEFAULT;
        N += dN * mdt;

        // Clamp to physiological bounds.
        c.mitoCount = (int)clampf(N, (float)SimCell::MITO_MIN, (float)SimCell::MITO_MAX);

        // Networked fraction tracks MFN/DRP1 balance (fusion wins → more
        // networked, fission wins → more fragmented).
        float netTarget = clampf(c.mfnActive / fmaxf(c.mfnActive + c.drp1Active, 0.001f), 0.0f, 1.0f);
        c.mitoNetworkedFrac += (netTarget - c.mitoNetworkedFrac) * mdt * 0.2f;

        // Health fraction: damage + ROS erode it; biogenesis + fusion
        // restore it (fusion re-homogenizes damaged components).
        float healthDecay = (c.ROS * 0.003f + c.damageLevel * 0.002f) * mdt;
        float healthRepair = (c.mfnActive * 0.002f + c.pgc1aActive * 0.001f) * mdt;
        c.mitoHealthFrac = clampf(c.mitoHealthFrac - healthDecay + healthRepair, 0.0f, 1.0f);
    }

    void updateMetabolism(SimCell& c, float dt) {
        float mdt=dt*MEDIUM_DT_SCALE;
        bool postDivisionDisplay = (mode == MODE_SINGLE_CELL && c.postDivisionRecovery > 0.0f);
        float lockedSize = c.size;
        float localO2=nutrients.getO2(c.position.x,c.position.z);
        float localGlu=nutrients.getGlucose(c.position.x,c.position.z);
        float localPH=nutrients.getPH(c.position.x,c.position.z);

        // ── CONTACT INHIBITION: G0 quiescence at confluence ─────────
        // Ref: Wikipedia Contact Inhibition; Delarue 2018 Dev Cell
        // At high density, cells enter G0: reduce metabolism, stop dividing, stay alive
        // This is the KEY mechanism for stationary phase (fill plate, stop growing)
        //
        // Quiescent cells have dramatically reduced metabolic demand:
        //   - O2 consumption: 30% of proliferating (Ref: Guppy 2002)
        //   - Glucose consumption: 25% of proliferating
        //   - ROS production: 20% of proliferating
        //   - This allows the colony to be self-sustaining at confluence
        bool isQuiescent = (c.fate == SIM_FATE_QUIESCENT);
        float metabolicActivity = isQuiescent ? 0.3f : (c.fate == SIM_FATE_PROLIF ? 1.0f : 0.7f);

        // ── Hypoxia tiered response ─────────────────────────────────
        float o2Def=fmaxf(0,HYPOXIA_MODERATE-localO2)/HYPOXIA_MODERATE;
        c.hypoxiaIntensity=c.hypoxiaIntensity*0.95f+o2Def*0.05f;
        if(localO2<HYPOXIA_SEVERE) {
            c.hypoxiaTimer+=mdt;
            if(c.hypoxiaTimer>HYPOXIA_NECROTIC_TIME&&localO2<HYPOXIA_NECROTIC_O2) c.necrotic=true;
        } else {
            c.hypoxiaTimer=fmaxf(0,c.hypoxiaTimer-mdt*0.5f);
            if(c.necrotic&&localO2>HYPOXIA_MODERATE) c.necrotic=false;
        }

        if(c.necrotic) {
            c.ATP=fmaxf(0,c.ATP-1.5f*mdt);
            c.stress=fminf(100,c.stress+0.8f*mdt);
            c.size=fminf(2.0f,c.size+0.001f*mdt);
        }

        // ── Warburg switch ──────────────────────────────────────────
        bool forcedGly=localO2<HYPOXIA_SEVERE;
        if(!forcedGly) {
            if(c.mitoHealth<MITO_HEALTH_SWITCH){c.warburgTimer+=mdt;if(c.warburgTimer>2.0f&&!c.glycolytic) c.glycolytic=true;}
            else{c.warburgTimer=fmaxf(0,c.warburgTimer-mdt*0.3f);if(c.mitoHealth>MITO_HEALTH_RESTORE&&c.glycolytic&&c.warburgTimer<0.3f) c.glycolytic=false;}
        } else if(!c.glycolytic) c.glycolytic=true;

        // ── Dual-pathway ATP ────────────────────────────────────────
        float glyB=c.glycolytic?WARBURG_GLY_BOOST:1.0f;
        float oxpM=c.glycolytic?WARBURG_OXP_PENALTY*c.mitoHealth:1.0f;
        if(forcedGly){glyB=2.2f;oxpM=0.1f;}
        float gly=localGlu*1.8f*c.glycolysisBias*glyB;
        // OxPhos output scales with mitochondrial capacity: ATP is
        // produced by the mito network, so more mitos (× health) =
        // more ATP. Normalized so a HeLa cell at the 500 target and
        // full health has OxPhos factor 1.0 — matches the previous
        // tuning. Empty/damaged network collapses OxPhos.
        // Ref: Brand & Nicholls 2011 Biochem J (coupled ATP output
        // proportional to active mito mass).
        float mitoCapacity = (float)c.mitoCount / (float)SimCell::MITO_TARGET_DEFAULT;
        mitoCapacity *= clampf(c.mitoHealthFrac, 0.0f, 1.0f);
        mitoCapacity = clampf(mitoCapacity, 0.15f, 2.5f); // guard rails
        float oxp=localGlu*localO2*6.5f*oxpM*mitoCapacity;

        // Quiescent cells have much lower metabolic costs (G0 state)
        float basalC = 0.3f * mdt * metabolicActivity;
        float cycleC = (c.phase==1?0.15f:c.phase==3?0.25f:0.05f) * mdt * metabolicActivity;
        float motC = (c.fate==SIM_FATE_PROLIF?0.12f:0.02f) * mdt;
        float repC = c.damageLevel>0 ? c.damageLevel*0.08f*mdt : 0;
        c.ATP=clampf(c.ATP+(gly+oxp)*mdt*metabolicActivity-basalC-cycleC-motC-repC, 0, 100);

        // ── Biologically real molecule counts (Milo 2013 Bioessays) ──
        // Driven by the same metabolic pathways; counts grow 1× → 2×
        // from G1→G2 and halve at division. Rates tuned so that over
        // a sim-scale cell cycle the integrated change ≈ 2×.
        // Phase multiplier: cycle progresses from 0 (early G1) to 1 (M).
        float cycleGrowth = clampf(c.cycleProgress, 0.0f, 1.0f);
        // First-order relaxation to the biological target count:
        //   dN/dt = k × (target − N)
        // So counts regrow toward the literature median after division
        // halved them, and stop at the target (don't blow up forever).
        // k chosen per-molecule based on turnover time: fast pools
        // (mRNA) converge in seconds, slow pools (ribosomes) in minutes.
        auto relax = [mdt](double& N, double target, float k) {
            N += (double)k * mdt * (target - N);
        };
        float mitoCap = (float)c.mitoCount / (float)SimCell::MITO_TARGET_DEFAULT
                        * clampf(c.mitoHealthFrac, 0.1f, 1.5f);
        float rateScale = metabolicActivity * fmaxf(mitoCap, 0.3f);
        // Proliferating cells have 2× synthesis target; quiescent stays
        // at 1× maintenance. Scale targets based on cycle progress so
        // ribosomes, proteins, and organelles grow through G1→G2.
        float cycleMul = 1.0f + clampf(c.cycleProgress, 0.0f, 1.0f) *
                         (c.fate == SIM_FATE_PROLIF ? 1.0f : 0.0f);

        relax(c.ribosomeCount, 6.0e6 * cycleMul,   0.040f * rateScale);
        relax(c.rnaPolII,      4.0e4 * cycleMul,   0.050f * rateScale);
        relax(c.rnaPolI,       2.0e2 * cycleMul,   0.050f * rateScale);
        relax(c.rnaPolIII,     1.0e3 * cycleMul,   0.050f * rateScale);
        relax(c.spliceosomes,  1.0e5 * cycleMul,   0.050f * rateScale);
        relax(c.nuclearPores,  4.0e3 * cycleMul,   0.030f * rateScale);
        relax(c.chaperones,    1.0e6 * cycleMul,   0.040f * rateScale);
        relax(c.proteasomes,   1.0e6 * cycleMul,   0.040f * rateScale);
        relax(c.copiiVesicles, 1.0e3 * cycleMul,   0.080f * rateScale);
        relax(c.secretoryVesicles, 3.0e3 * cycleMul, 0.080f * rateScale);
        relax(c.lysosomes,     3.0e2 * cycleMul,   0.020f * rateScale);
        relax(c.peroxisomes,   3.0e2 * cycleMul,   0.020f * rateScale);
        relax(c.lipidDroplets, 1.0e2,              0.015f * rateScale);
        relax(c.microtubulesTotal, 3.0e2 * cycleMul, 0.060f * rateScale);
        relax(c.actinFilaments,   5.0e4 * cycleMul, 0.060f * rateScale);
        relax(c.totalProteins, 4.0e9 * cycleMul,   0.040f * rateScale);
        // mRNA: faster turnover — converges in ~10 s sim time.
        relax(c.mRNACount,     3.6e5 * cycleMul,   0.100f * metabolicActivity);
        // tRNA: stable pool, moderate relaxation.
        relax(c.tRNACount,     1.0e7 * cycleMul,   0.050f * rateScale);
        // ATP: consumed at ~10⁹/s (turnover = whole pool every ~2 s);
        // balanced with production below. We track a count derived from
        // c.ATP (mM) and cell volume so the two stay consistent.
        // Cell volume ~4 pL; 1 mM ≈ 2.4×10⁹ molecules.
        // c.ATP is in "charge units" 0–100 mapped to 0–4 mM so each unit
        // ≈ 2.4×10⁷ ATP molecules.
        c.atpMolecules = (double)c.ATP * 7.5e7;  // 100 units → 7.5×10⁹
        // Glucose consumption: ~3×10⁸ molecules per second per cell.
        c.glucoseMolecules = clampf(
            c.glucoseMolecules * (1.0f - 0.002f * metabolicActivity * mdt)
                              + localGlu * 1.0e8 * mdt,
            1.0e9, 5.0e10);
        // NADH shuttles in TCA; turnover fast; count tracks production.
        c.nadhMolecules = clampf(
            c.nadhMolecules * (1.0f - 0.015f * mdt) + oxp * 5.0e6 * mdt,
            1.0e7, 1.0e9);
        // Amino acid pool — consumed by translation proportional to ribosome count.
        float aaConsume = (float)c.ribosomeCount * 1e-7f * metabolicActivity * mdt;
        c.aaMolecules = clampf(
            c.aaMolecules * (1.0f - aaConsume * 0.05f)
                         + localGlu * 1.5e8 * metabolicActivity * mdt,
            1.0e9, 1.0e11);
        // Water exchange is enormous (~10¹⁴ flux) but pool stays near
        // steady state. Fix small variance.
        c.waterMolecules = 1.3e14 * (1.0f + 0.02f * sinf(bioTime * 0.05f));
        // Free cytosolic Ca²⁺ is dynamic (oscillates with signaling).
        c.calciumFree = 2.4e5 * (0.8f + 0.4f * sinf(bioTime * 0.7f + c.cellUid * 0.31f));
        // cAMP: rises on GPCR activation (here we just oscillate gently).
        c.cAMPMolecules = 2.4e4 * (0.7f + 0.6f * sinf(bioTime * 0.12f + c.cellUid * 0.19f));

        // ── Genome doubling through S-phase ─────────────────────────
        // genomeBp goes 2N → 4N during S, halves at division.
        if (c.phase == 1) {
            float rep = clampf(c.program.cdogma.replicationProgress, 0.0f, 1.0f);
            c.genomeBp = 6.4e9 * (1.0 + rep);  // 1× → 2×
            c.nucleosomes = 3.0e7 * (1.0 + rep);
            // Replication forks peak mid-S at ~30,000 active forks.
            float forkShape = 4.0f * rep * (1.0f - rep);  // parabola, peaks at rep=0.5
            c.replicationForks = 3.0e4 * forkShape;
        } else {
            c.replicationForks = 0.0;
        }

        // ── Central-dogma machinery growth (proliferating cells) ─────
        if (c.fate == SIM_FATE_PROLIF) {
            float mitoCap = (float)c.mitoCount / (float)SimCell::MITO_TARGET_DEFAULT;
            mitoCap *= clampf(c.mitoHealthFrac, 0.1f, 1.5f);
            float growthRate = 0.012f * mitoCap * (0.3f + 0.7f * cycleGrowth);
            // RNA polymerases: made from proteins, grow with ribosome biogenesis.
            c.rnaPolII   = clampf(c.rnaPolII   * (1.0f + growthRate * mdt), 1.0e4, 2.0e5);
            c.rnaPolI    = clampf(c.rnaPolI    * (1.0f + growthRate * mdt), 5.0e1, 1.0e3);
            c.rnaPolIII  = clampf(c.rnaPolIII  * (1.0f + growthRate * mdt), 2.0e2, 5.0e3);
            c.spliceosomes = clampf(c.spliceosomes * (1.0f + growthRate * mdt), 5.0e4, 5.0e5);
            // Nuclear pores assemble at ~2000/hour in proliferating cells.
            c.nuclearPores = clampf(c.nuclearPores * (1.0f + growthRate * 0.5f * mdt), 1.5e3, 1.0e4);
            // Chaperones / proteasomes scale with total protein synth.
            c.chaperones = clampf(c.chaperones * (1.0f + growthRate * mdt), 3.0e5, 3.0e6);
            c.proteasomes = clampf(c.proteasomes * (1.0f + growthRate * mdt), 3.0e5, 3.0e6);
            // Total cellular protein doubles per cycle (~4e9 → 8e9).
            c.totalProteins = clampf(c.totalProteins * (1.0f + growthRate * mdt), 1.0e9, 1.0e10);
            // Cytoskeleton grows proportionally with cell.
            c.microtubulesTotal = clampf(c.microtubulesTotal * (1.0f + growthRate * 0.7f * mdt), 1.5e2, 8.0e2);
            c.actinFilaments    = clampf(c.actinFilaments    * (1.0f + growthRate * 0.7f * mdt), 2.0e4, 1.5e5);
            // Membrane traffic scales with secretory demand.
            c.copiiVesicles     = clampf(c.copiiVesicles     * (1.0f + growthRate * 0.5f * mdt), 5.0e2, 5.0e3);
            c.secretoryVesicles = clampf(c.secretoryVesicles * (1.0f + growthRate * 0.5f * mdt), 1.0e3, 1.0e4);
            // Digestive / metabolic organelles — slow growth.
            c.lysosomes  = clampf(c.lysosomes  * (1.0f + growthRate * 0.3f * mdt), 1.5e2, 8.0e2);
            c.peroxisomes = clampf(c.peroxisomes * (1.0f + growthRate * 0.3f * mdt), 1.5e2, 8.0e2);
            // Lipid droplets depend on excess acetyl-CoA (simplify: slow drift)
            c.lipidDroplets = clampf(c.lipidDroplets * (1.0f + growthRate * 0.2f * mdt), 3.0e1, 5.0e2);
        }

        // ── Nutrient exchange (closed-system stoichiometric reactions) ─
        // Build a 12-species flux array (mM / bio-second) from the cell's
        // current biology, then atomically debit/credit the local grid
        // cell. Every entry comes from a balanced reaction (R1–R7); we
        // never inject mass out of nothing.
        //
        // Rate constants are calibrated so a confluent dish (~600 cells)
        // exhausts its glucose pool over ~24 bio-h, matching DMEM
        // refresh practice (Reitzer 1979 / BNID 105005).
        const float MAX_FLUX_MM_PER_BIOSEC = 5.0e-4f;          // saturating per-cell rate
        float bioDtLocal = bioDt(dt, timeScale);
        float consumeFactor = metabolicActivity * (c.fate==SIM_FATE_PROLIF ? 1.8f : 1.0f);
        float exMitoCap = clampf((float)c.mitoCount / (float)SimCell::MITO_TARGET_DEFAULT, 0.0f, 2.0f)
                          * clampf(c.mitoHealthFrac, 0.1f, 1.5f);

        // Substrate availability (Michaelis-Menten saturation).
        float gluLocal = nutrients.get(MS_GLUCOSE,   c.position.x, c.position.z);
        float gln      = nutrients.get(MS_GLUTAMINE, c.position.x, c.position.z);
        float o2L      = nutrients.get(MS_O2,        c.position.x, c.position.z);
        float aaL      = nutrients.get(MS_AA_POOL,   c.position.x, c.position.z);
        float pyrL     = nutrients.get(MS_PYRUVATE,  c.position.x, c.position.z);

        auto mm = [](float S, float Km) { return S / (S + Km + 1e-9f); };

        float r_glycolysis  = MAX_FLUX_MM_PER_BIOSEC * consumeFactor * mm(gluLocal, 5.0f);
        float r_oxphos      = MAX_FLUX_MM_PER_BIOSEC * consumeFactor * exMitoCap
                              * mm(o2L, 0.05f) * mm(pyrL, 0.5f);
        float r_glutamin    = 0.4f * MAX_FLUX_MM_PER_BIOSEC * consumeFactor * mm(gln, 0.5f);
        float r_translation = 0.5f * MAX_FLUX_MM_PER_BIOSEC * consumeFactor * mm(aaL, 0.5f)
                              * clampf((float)c.ribosomeCount / 1.0e6f, 0.1f, 2.0f);
        float r_warburg     = c.glycolytic ? 0.8f : 0.1f;     // fraction of pyruvate → lactate

        float flux[MS_COUNT] = {0};
        // R1 Glycolysis: glucose → 2 pyruvate
        flux[MS_GLUCOSE]  -= r_glycolysis;
        flux[MS_PYRUVATE] += 2.0f * r_glycolysis;
        // R2 Lactate fermentation: pyruvate → lactate (Warburg cells dump more)
        float r_lactate = r_warburg * (2.0f * r_glycolysis + r_oxphos);
        flux[MS_PYRUVATE] -= r_lactate;
        flux[MS_LACTATE]  += r_lactate;
        // R3 TCA + OxPhos: pyruvate + 2.5 O2 → 3 CO2
        flux[MS_PYRUVATE] -= r_oxphos;
        flux[MS_O2]       -= 2.5f * r_oxphos;
        flux[MS_CO2]      += 3.0f * r_oxphos;
        // R4 Glutaminolysis: glutamine + 0.5 O2 → 2 CO2
        flux[MS_GLUTAMINE] -= r_glutamin;
        flux[MS_O2]        -= 0.5f * r_glutamin;
        flux[MS_CO2]       += 2.0f * r_glutamin;
        // R5 Translation: AA pool → biomass (consumed only)
        flux[MS_AA_POOL]   -= r_translation;
        // R6 Receptor binding — catalytic, no net flux
        // R7 Cell autolysis — handled at apoptosis / necrosis site
        // R8 Aquaporin osmosis — water moves down osmolarity gradient.
        //    Sum extracellular non-water non-pH non-drug species, compare
        //    to cytosolic baseline; positive Δ pulls water INTO the cell
        //    (cell swells); negative Δ pushes water OUT (cell shrinks).
        float osmoExt = 0.0f;
        for (int s = 0; s < MS_COUNT; s++) {
            if (s == MS_HPLUS || s == MS_DRUG || s == MS_WATER || s == MS_GROWTH_F) continue;
            osmoExt += nutrients.get(s, c.position.x, c.position.z);
        }
        float dOsmo = osmoExt - c.osmoCytoMM;
        float Jw_mM = Osmosis::AQP_LP_PER_BIOSEC * (float)c.aquaporinCount * dOsmo;
        flux[MS_WATER] -= Jw_mM;       // water leaves medium (rises if dOsmo>0)
        c.waterMM      += Jw_mM * bioDtLocal;
        // Turgor: deviation from baseline 55 000 mM water → elastic
        // pressure. Capped to ±5 kPa to prevent runaway swelling.
        float volRatio = c.waterMM / MediumComposition::DMEM_WATER_MM;
        c.turgorPa = clampf(Osmosis::BULK_MOD_PA * (volRatio - 1.0f),
                            -Osmosis::TURGOR_LYSE_PA, Osmosis::TURGOR_LYSE_PA);

        nutrients.exchange(c.position.x, c.position.z, flux, bioDtLocal);
        nutrients.updatePH();

        // ── Stress homeostasis ──────────────────────────────────────
        // Quiescent cells have lower stress because they consume less
        float pHS=fmaxf(0,(0.68f-localPH)/0.06f)*4.0f*mdt;
        float atpS=fmaxf(0,(15-c.ATP)*0.03f)*mdt;
        // Recovery scales with metabolic rest — quiescent cells recover faster
        float recovery = (isQuiescent ? 3.0f : 1.8f) * mdt;
        c.stress=clampf(c.stress + atpS + pHS + (1-localO2)*0.3f*mdt - recovery, 0, 100);

        // ── ROS ─────────────────────────────────────────────────────
        // Quiescent cells produce much less ROS (reduced metabolic activity)
        float rosP = (c.fate==SIM_FATE_PROLIF ? 0.08f*c.prolifBias : 0.015f) * mdt * metabolicActivity;
        float rosC = 0.06f*c.rosTolerance*(c.ATP/100)*mdt;
        c.ROS=clampf(c.ROS+rosP-rosC, 0, 100);

        // ── Damage ──────────────────────────────────────────────────
        float dI = c.ROS*0.0008f*mdt + (c.stress>80 ? (c.stress-80)*0.0001f*mdt : 0);
        float dO = c.repairRate*(c.ATP/100)*0.012f*mdt;
        c.damageLevel=clampf(c.damageLevel+dI-dO, 0, 2);

        // ── Mito health ─────────────────────────────────────────────
        c.mitoHealth=clampf(c.mitoHealth+0.002f*(c.ATP/100)*mdt-c.ROS*0.00006f*mdt, 0, 1);
        c.mitoPotential=clampf(c.mitoPotential+(180-c.mitoPotential)*0.003f*mdt-c.ROS*0.04f*mdt, 40, 220);

        // ── ATP danger → apoptosis (only after sustained collapse) ──
        if(c.ATP<ATP_DANGER_THRESHOLD&&c.apoptosisPhase==0){
            c.atpDangerTimer+=mdt;
            if(c.atpDangerTimer>=ATP_DANGER_DURATION&&!c.fateLocked)
                {c.fate=SIM_FATE_APOPTOTIC;c.apoptosisPhase=1;c.apoTimer=0;c.fateLocked=true;}
        } else c.atpDangerTimer=fmaxf(0,c.atpDangerTimer-mdt*0.5f);

        // ── Emergency apoptosis triggers (genuine damage only) ──────
        if(c.apoptosisPhase==0) {
            if(c.mitoPotential<45&&c.ROS>95){c.fate=SIM_FATE_APOPTOTIC;c.apoptosisPhase=1;c.apoTimer=0;c.fateLocked=true;}
            if(c.damageLevel>1.2f&&!c.fateLocked){c.fate=SIM_FATE_APOPTOTIC;c.apoptosisPhase=1;c.apoTimer=0;c.fateLocked=true;}
        }

        // ── Biomass growth ──────────────────────────────────────────
        if(!c.necrotic&&c.apoptosisPhase==0) {
            float presP=fmaxf(0,1.0f-c.localPressure*0.12f);
            float syn=BIOMASS_SYNTH_K*(c.ATP/100)*localGlu;
            float deg=BIOMASS_DEGRADE_K*(c.stress/100);
            float growthMult = (c.fate==SIM_FATE_PROLIF||c.fate==SIM_FATE_UNDETERMINED) ? 1.5f : 0.3f;
            c.biomass=clampf(c.biomass+(syn-deg)*mdt*presP*growthMult, 0.4f, 2.3f);
        }

        // ── Mechanical quorum → p27/p21 induction ───────────────────
        // Ref: Delarue 2018 — pressure induces CDK inhibitors → G1 arrest
        // This is the PRIMARY mechanism for contact inhibition at confluence
        if(c.localPressure>0.8f) {
            float p21Induction = c.localPressure * MECH_P21_COUPLING * mdt;
            c.cdk.p21 = fminf(1.0f, c.cdk.p21 + p21Induction);
        }

        // ── Natural turnover at confluence ──────────────────────────
        // Real tissue: ~1-3% daily cell death even in homeostasis
        if (mode != MODE_SINGLE_CELL && c.age > TURNOVER_AGE_THRESHOLD && c.apoptosisPhase == 0) {
            float deathProb = TURNOVER_PROB_PER_DT * (c.age - TURNOVER_AGE_THRESHOLD) * 0.001f;
            if (randf() < deathProb * mdt) {
                c.fate = SIM_FATE_APOPTOTIC; c.apoptosisPhase = 1; c.apoTimer = 0;
                c.fateLocked = true;
            }
        }

        // Base size from biomass + turgor-driven elastic swelling.
        // 5 kPa turgor → 5 % radius change (Stewart 2011 Nat Rev MCB).
        float turgorFrac = clampf(c.turgorPa / Osmosis::BULK_MOD_PA, -0.20f, 0.20f);
        c.size=clampf((0.6f+c.biomass*0.4f) * (1.0f + turgorFrac), 0.5f, 1.8f);
        if (postDivisionDisplay && c.apoptosisPhase == 0 && !c.necrotic) {
            // In single-cell teaching mode, after mitosis we want two healthy,
            // persistent daughters rather than a second round of growth/death.
            c.size = clampf(lockedSize, 0.5f, 1.8f);
            c.biomass = biomassFromCellSizeFactor(c.size);
            c.ATP = fmaxf(c.ATP, 35.0f);
            c.stress = fminf(c.stress, 20.0f);
            c.damageLevel = fminf(c.damageLevel, 0.10f);
            c.mitoPotential = fmaxf(c.mitoPotential, 150.0f);
        }
        c.age+=mdt;
    }

    // ── Drug pharmacodynamics (PhysiPKPD model) ─────────────────────
    // Ref: Bergman et al., GigaByte 2023 (BSD-3)
    // 1. Uptake: dC_int/dt = uptake × C_ext − efflux × C_int
    // 2. Damage: Hill(C_int, EC50, n) × maxEffect − repair × Damage
    // 3. Effect: depends on MOA (mechanism of action)
    void updateDrugResponse(SimCell& c, float dt) {
        if (activeDrugIdx <= 0) return; // no drug active
        const Drug& d = DRUG_LIBRARY[activeDrugIdx];
        float mdt = dt * MEDIUM_DT_SCALE;

        // Drug uptake from local environment
        float localDrug = nutrients.getDrug(c.position.x, c.position.z);
        float resistFactor = c.drugResistant ? 0.2f : 1.0f; // MDR resistance
        float uptake = d.uptakeRate * localDrug * resistFactor * mdt;
        float efflux = d.effluxRate * c.drugInternal * mdt;
        c.drugInternal = fmaxf(0, c.drugInternal + uptake - efflux);
        nutrients.consumeDrug(c.position.x, c.position.z, uptake);

        // Damage accumulation (Hill equation dose-response)
        // Hill(C, EC50, n) = C^n / (EC50^n + C^n)
        float Cn = powf(c.drugInternal, d.hillCoeff);
        float EC50n = powf(d.EC50, d.hillCoeff);
        float hillResponse = Cn / (EC50n + Cn + 1e-12f);
        float damageIn = hillResponse * d.maxEffect * mdt;
        float damageOut = d.damageRepairRate * c.drugDamage * mdt;
        c.drugDamage = clampf(c.drugDamage + damageIn - damageOut, 0, 1);

        // Apply MOA effects
        auto applyMOA = [&](int moa) {
            if (moa < 0) return;
            switch (moa) {
                case MOA_ANTI_PROLIF:
                    // CDK inhibition: freeze cycle, induce p21
                    if (c.drugDamage > 0.3f) {
                        c.cdk.p21 = fminf(1.0f, c.cdk.p21 + c.drugDamage * 0.1f * mdt);
                    }
                    break;
                case MOA_PRO_APOPTOSIS:
                    // Direct apoptosis trigger at high damage
                    if (c.drugDamage > 0.7f && c.apoptosisPhase == 0) {
                        c.fate = SIM_FATE_APOPTOTIC;
                        c.apoptosisPhase = 1; c.apoTimer = 0;
                        c.fateLocked = true;
                    }
                    break;
                case MOA_DNA_DAMAGE:
                    // Adds to existing DNA damage model → p53 → apoptosis
                    c.damageLevel += c.drugDamage * 0.05f * mdt;
                    break;
                case MOA_MITO_TOXIN:
                    // Collapse mitochondrial membrane potential
                    c.mitoPotential -= c.drugDamage * 2.0f * mdt;
                    c.mitoPotential = fmaxf(40, c.mitoPotential);
                    break;
            }
        };
        applyMOA(d.mechanism);
        applyMOA(d.mechanism2);
    }

    void updateCellCycle(SimCell& c, float dt) {
        if(c.apoptosisPhase>0||c.senescent) return;

        // Fresh daughters should settle in early G1 before immediately
        // resuming full growth/cycling. This keeps the lineage continuous and
        // avoids the fake-looking "small daughters -> big M-phase cells" snap.
        if (mode == MODE_SINGLE_CELL && c.postDivisionRecovery > 0.0f) {
            c.phase = 0;
            c.cycleProgress = 0;
            c.cycleTimer = 0;
            c.divisionPending = false;
            c.checkpointG1Passed = false;
            c.checkpointG2Passed = false;
            c.g1WaitTimer = 0;
            c.g2WaitTimer = 0;
            // Do NOT raise divisionCooldown here. The previous code forced
            // it to at least 18 slow-time units, which (with SLOW_DT_SCALE
            // =0.06) meant ~5 real minutes before the cell could enter M
            // again. The daughter's cooldown is already set appropriately
            // in deriveDaughter() (1.5 single / 8 colony), and
            // postDivisionRecovery itself already holds phase=0 for 6 s.
            return;
        }

        // Once a cell has entered mitosis, the cell-cycle timer no longer
        // owns progression. The per-cell mitosis program does.
        if (c.program.mitosis.active) {
            c.phase = 3;
            c.cycleProgress = 1.0f;
            c.divisionPending = c.program.mitosis.postDivisionComplete();
            c.checkpointG1Passed = false;
            c.checkpointG2Passed = true;
            return;
        }

        float sdt=dt*SLOW_DT_SCALE;

        // ── Lag phase: cells adapt before first division ────────────
        if(c.adaptationTimer < LAG_DURATION) {
            c.adaptationTimer += sdt;
            c.cycleProgress = 0;
            return; // no cycle progression during lag
        }

        float localO2=nutrients.getO2(c.position.x,c.position.z);
        float localGlu=nutrients.getGlucose(c.position.x,c.position.z);
        float localPH=nutrients.getPH(c.position.x,c.position.z);
        bool pHOk=localPH>0.64f;

        // ── CONTACT INHIBITION: Hippo-YAP G1 arrest ────────────────
        // Ref: Nature Comms 2018 — YAP/TAZ → cytoplasm at high density
        // When mechanical pressure is high, cell cycle HALTS in G1.
        // This is the primary mechanism that stops growth at confluence.
        //
        // Count immediate neighbors for this cell
        int neighborCount = 0;
        for (auto& o : cells) {
            if (&o == &c || !o.alive) continue;
            float dx = c.position.x - o.position.x;
            float dz = c.position.z - o.position.z;
            if (sqrtf(dx*dx + dz*dz) < c.radius * 2.5f) neighborCount++;
        }

        // Confluence arrest: if surrounded by >= 5 neighbors, cell cycle STOPS
        // This mimics the Hippo pathway: mechanical compression → YAP nuclear exclusion → G1 arrest
        bool contactArrested = (neighborCount >= 5) || (c.localPressure > 1.0f);

        // Growth signal (Michaelis-Menten)
        float o2mm=localO2/(localO2+0.10f);
        float glumm=localGlu/(localGlu+0.15f);
        float atpmm=(c.ATP/100)/((c.ATP/100)+0.25f);
        float gs=fminf(1.0f, o2mm*0.40f+glumm*0.30f+atpmm*0.20f+0.05f);

        // CDK ODE: only advances if NOT contact-arrested
        if (!contactArrested) {
            c.cdk.step(sdt * CDK_DT_SCALE, gs);
        } else {
            // Contact arrest → p21/p27 induced, CycD suppressed (Hippo pathway)
            c.cdk.p21 = fminf(1.0f, c.cdk.p21 + 0.05f * sdt);
            c.cdk.CycD = fmaxf(0, c.cdk.CycD - 0.02f * sdt);
        }
        if(c.damageLevel>0.2f) c.cdk.p21=fminf(1.0f,c.cdk.p21+c.damageLevel*0.35f*0.05f);

        c.phase=c.cdk.getPhase();
        float g1end=CYCLE_G1_DUR, se=g1end+CYCLE_S_DUR, g2end=se+CYCLE_G2_DUR;

        // Timer only advances if NOT contact-arrested
        if (!contactArrested) {
            c.cycleTimer += sdt;
        }
        // If arrested, force cell to stay in G1 and convert to quiescent
        if (contactArrested && c.cycleTimer < g1end) {
            if (c.fate != SIM_FATE_QUIESCENT && c.fate != SIM_FATE_APOPTOTIC) {
                c.fate = SIM_FATE_QUIESCENT;
                c.fateLocked = true;
            }
        }

        // ── Sharp predicate-based G1/S checkpoint ──────────────────
        // Replaces the old stochastic Hill product. Cell either passes the
        // restriction point (all gates met) or holds at end-of-G1, with
        // the *exact* failing predicate written to the audit log so the
        // user can answer "why didn't this cell divide?" with a citation.
        // No randf() — the gate is fully deterministic given the state.
        if(c.cycleTimer<g1end) {
            c.cycleProgress=c.cycleTimer/g1end;
            if(!c.checkpointG1Passed && !contactArrested) {
                float gluLocalMM = nutrients.get(MS_GLUCOSE,   c.position.x, c.position.z);
                float glnLocalMM = nutrients.get(MS_GLUTAMINE, c.position.x, c.position.z);
                CheckpointResult g1s = evalG1S(c, gluLocalMM, glnLocalMM, contactArrested);
                if (g1s.passed) {
                    c.checkpointG1Passed = true;
                    c.program.logEvent(bioTime, 0, "G1S_PASS", g1s.reason);
                } else {
                    c.g1WaitTimer += sdt;
                    c.cycleTimer = g1end * 0.99f;
                    // Log only intermittently to avoid filling the ring
                    // every frame with the same reason.
                    if (c.g1WaitTimer > CYCLE_G1_DUR * 0.10f &&
                        ((int)(c.g1WaitTimer / sdt) % 256) == 0) {
                        c.program.logEvent(bioTime, 0, "G1S_HOLD", g1s.reason);
                    }
                }
            }
        } else if(c.cycleTimer<se) {
            c.cycleProgress=(c.cycleTimer-g1end)/CYCLE_S_DUR;
            // Per-cell S-phase gate: every cell has its own replication
            // program (c.program.cdogma). S cannot complete visually faster
            // than replication actually progresses. This replaces the
            // previous "primary-only" special case (which used the global
            // primaryReplicationProgress) and extends the same biology to
            // every lineage in the dish.
            //
            // Reference: Chao et al. 2019 eLife single-cell tracking —
            // S-phase duration in HeLa is ~8 h and scales with replication
            // fork speed + origin firing, not a wall clock.
            if (c.program.cdogmaInitialized) {
                float replProgress = clampf(c.program.cdogma.replicationProgress, 0.0f, 0.999f);
                c.cycleProgress = fminf(c.cycleProgress, replProgress);
                bool replReady = c.program.cdogma.replicationReadyForM();
                if (!replReady && c.cycleProgress > 0.985f) {
                    // Hold at end of S until the real replication program
                    // finishes all forks and clears chk1Signal.
                    c.cycleTimer = g1end + CYCLE_S_DUR * 0.985f;
                    c.cycleProgress = 0.985f;
                    c.checkpointG2Passed = false;
                }
            }
        } else if(c.cycleTimer<g2end) {
            c.cycleProgress=(c.cycleTimer-se)/CYCLE_G2_DUR;
            if(!c.checkpointG2Passed) {
                // Sharp predicate-based G2/M gate (replaces stochastic Hill product).
                // Holds the cycle clock at G2 boundary until every required
                // condition is met — no probabilistic noise, no safety valve.
                c.program.ensureCDogmaInitialized();
                CheckpointResult g2m = evalG2M(c, c.program.cdogma);
                if (g2m.passed) {
                    c.checkpointG2Passed = true;
                    c.program.logEvent(bioTime, 2, "G2M_PASS", g2m.reason);
                } else {
                    c.g2WaitTimer += sdt;
                    c.cycleTimer = g2end * 0.99f;
                    c.cycleProgress = fminf(c.cycleProgress, 0.995f);
                    if (c.g2WaitTimer > CYCLE_G2_DUR * 0.10f &&
                        ((int)(c.g2WaitTimer / sdt) % 256) == 0) {
                        c.program.logEvent(bioTime, 2, "G2M_HOLD", g2m.reason);
                    }
                }
            }
        } else {
            // ── Division cooldown: prevent immediate re-entry ───────
            if (c.divisionCooldown > 0) {
                c.divisionCooldown -= sdt;
                c.cycleTimer = g2end * 0.99f; // Hold at G2/M boundary
                c.cycleProgress = 0.95f;
                return; // Don't enter M-phase yet
            }

            c.cycleProgress = 1.0f;
            c.phase = 3; // M phase

            if (c.program.mitosis.active) {
                c.cycleTimer = fmaxf(c.cycleTimer, g2end + fminf(c.program.mitosis.totalProgress, CYCLE_M_DUR));
                return;
            }

            // ── M-entry licensing via the SHARP `evalG2M` predicate ──
            // All biology gates collapse into one call returning
            // `{passed, reason}`. The first failing predicate is what gets
            // written to the audit log verbatim — same string the user
            // would read in the in-app Audit panel — so "why didn't this
            // cell divide" always has a citable answer.
            //
            // Refused licenses with reasons that are NOT covered by
            // evalG2M (space, dish full, fate, senescence, tetraploid)
            // are logged inline.
            bool hasSpace = neighborCount < 4 && c.localPressure < 1.5f;
            bool canDiv = (c.fate==SIM_FATE_PROLIF || c.fate==SIM_FATE_UNDETERMINED);
            bool notTetraploid = !c.program.tetraploid;
            c.program.ensureCDogmaInitialized();
            CheckpointResult m_entry = evalG2M(c, c.program.cdogma);

            if (canDiv && !c.senescent && hasSpace && notTetraploid &&
                m_entry.passed && (int)cells.size() < MAX_CELLS) {
                c.divisionPending = false;
                if(c.fate==SIM_FATE_UNDETERMINED) {
                    c.fate = SIM_FATE_PROLIF; c.fateLocked = true;
                }
                c.cycleTimer = fmaxf(c.cycleTimer, g2end);
                c.program.logEvent(bioTime, 3, "M_LICENSED", m_entry.reason);
                return;
            }

            // Log the specific reason that refused the license.
            const char* refuseReason = m_entry.reason;
            if (!canDiv)                     refuseReason = "fate_nonprolif";
            else if (c.senescent)            refuseReason = "senescent";
            else if (!notTetraploid)         refuseReason = "tetraploid_from_slippage";
            else if (!hasSpace)              refuseReason = "no_space";
            else if ((int)cells.size() >= MAX_CELLS) refuseReason = "dish_full";
            c.program.logEvent(bioTime, 3, "DIV_REFUSED", refuseReason);

            // Reset for next cycle
            c.cycleTimer = 0; c.cycleProgress = 0;
            c.checkpointG1Passed = false; c.checkpointG2Passed = false;
            c.g1WaitTimer = 0; c.g2WaitTimer = 0;
            c.cdk.resetForNewCycle(gs);
        }
    }

    void updateFate(SimCell& c, float dt) {
        if(c.apoptosisPhase>0||c.senescent) return;
        if (mode == MODE_SINGLE_CELL && c.postDivisionRecovery > 0.0f) return;

        // Fate timer in MEDIUM timescale
        float fdt = dt * MEDIUM_DT_SCALE;
        c.fateTimer += fdt;

        float localO2=nutrients.getO2(c.position.x,c.position.z);
        float localGlu=nutrients.getGlucose(c.position.x,c.position.z);
        float nutrition=localO2*localGlu;

        // ── REVERSIBLE G0: quiescent cells re-enter cycle if space opens ──
        // Ref: PMC6496145 — "Contact inhibition is a reversible form of cell cycle arrest"
        // If a quiescent cell has few neighbors + good nutrients → re-enter proliferation
        if(c.fate == SIM_FATE_QUIESCENT && c.fateLocked) {
            int nbrs = 0;
            for(auto& o:cells){if(&o==&c||!o.alive) continue;
                float dx2=c.position.x-o.position.x,dz2=c.position.z-o.position.z;
                if(sqrtf(dx2*dx2+dz2*dz2)<c.radius*2.8f) nbrs++;}
            bool spaceAvailable = nbrs < 3 && c.localPressure < 0.3f;
            if(spaceAvailable && nutrition > 0.20f && c.ATP > 40) {
                c.fate = SIM_FATE_PROLIF;
                c.fateScores[0] = 8; c.fateScores[1] = 0;
                c.cdk.p21 = fmaxf(0, c.cdk.p21 - 0.3f); // release p27/p21 block
            }
        }

        int neighbors=0;
        for(auto& o:cells){if(&o==&c||!o.alive) continue;
            float dx=c.position.x-o.position.x,dz=c.position.z-o.position.z;
            if(sqrtf(dx*dx+dz*dz)<c.radius*2.8f) neighbors++;}
        bool crowded = neighbors >= (int)CONTACT_INHIBIT_NBRS;

        // ── Contact inhibition: density-dependent G0 arrest ─────────
        // Ref: Wikipedia Contact Inhibition; Nature Comms (YAP/TAZ axis)
        // At confluence, cells enter reversible G0 via p27/Hippo pathway
        // This is NOT stress-induced death — it's an active survival program

        // Proliferation: needs nutrients, space, energy
        bool canProlif = nutrition > 0.10f && !crowded && c.ATP > 18 && c.stress < 72;
        if(canProlif) c.fateScores[0] += fdt * 0.5f * c.prolifBias;
        else c.fateScores[0] *= powf(0.990f, fdt*60);

        // Quiescence: crowding is the DOMINANT signal (contact inhibition)
        // Also triggered by moderate stress where cell can survive by resting
        bool quiesceSignal = crowded || (c.localPressure > 0.5f) || (c.stress > 40 && c.ATP > 15);
        float quiesceRate = crowded ? 1.5f : (c.localPressure > 0.5f ? 0.8f : 0.3f);
        if(quiesceSignal) c.fateScores[1] += fdt * quiesceRate;
        else c.fateScores[1] *= powf(0.990f, fdt*60);

        // Death: only from genuine severe damage, NOT from crowding stress
        // Apoptosis is a last resort, not the primary response to confluence
        if(c.ATP<4||c.damageLevel>0.8f) c.fateScores[2]+=fdt*0.4f;
        else c.fateScores[2] *= powf(0.990f, fdt*60);

        for(int i=0;i<3;i++) c.fateScores[i]=fminf(25,c.fateScores[i]);

        // Commit fate — timer in medium-dt units (~12-30 medium-time units)
        float commitTime = 12.0f + randf()*18.0f;
        if(!c.fateLocked && c.fateTimer > commitTime) {
            float maxS=fmaxf(c.fateScores[0],fmaxf(c.fateScores[1],c.fateScores[2]));
            if(maxS > 4) { // lower threshold for faster commitment
                float T=4.0f;
                float e0=expf((c.fateScores[0]-maxS)/T);
                float e1=expf((c.fateScores[1]-maxS)/T);
                float e2=expf((c.fateScores[2]-maxS)/T);
                float tot=e0+e1+e2, r=randf()*tot;
                if(r<e0) c.fate=SIM_FATE_PROLIF;
                else if(r<e0+e1) c.fate=SIM_FATE_QUIESCENT;
                else {c.fate=SIM_FATE_APOPTOTIC;c.apoptosisPhase=1;c.apoTimer=0;}
                c.fateLocked=true;
            }
        }
    }

    void updateApoptosis(SimCell& c, float dt) {
        if(c.apoptosisPhase==0) {
            if(c.necrotic&&c.ATP<2){c.alive=false;statDeaths++;} return;
        }
        c.apoTimer+=dt;
        c.size=fmaxf(0.1f,c.size-0.003f*dt);
        c.ATP=fmaxf(0,c.ATP-0.5f*dt);
        if(c.apoTimer>8.0f){c.alive=false;statDeaths++;}
    }

    void processDivisions() {
        // cells[0] in single-cell mode is finalized by the detailed visual
        // mitosis in main.mm (gMitosis). Here we only route its
        // divisionPending flag back into the "hold in M-phase for visual"
        // state that main.mm watches. EVERY OTHER cell — including
        // background cells in single mode — finalizes through the normal
        // deriveDaughter() loop below. Previously this function early-
        // returned in single mode, which left daughter B (and every
        // subsequent background cell) permanently stuck at DNA=1524 /
        // phase=3 with no one to spawn the daughters.
        if (mode == MODE_SINGLE_CELL && !cells.empty()) {
            auto& m = cells[0];
            if (m.divisionPending) {
                m.divisionPending = false;
                m.phase = 3; // Stay in M-phase for visual mitosis
            }
            if (mitosisVisualizationComplete) {
                mitosisVisualizationComplete = false;
            }
        }

        // Colony mode and non-primary single-cell lineages use deriveDaughter()
        // once their own mitosis program reaches cytokinesis completion.
        // In single-cell mode, start at index 1 (cells[0] handled above).
        int startIdx = (mode == MODE_SINGLE_CELL) ? 1 : 0;
        std::vector<SimCell> daughters;
        daughters.reserve(cells.size());

        for (int i = startIdx; i < (int)cells.size(); i++) {
            auto& m = cells[i];
            bool mitosisReady = m.program.mitosis.postDivisionComplete();
            if ((!m.divisionPending && !mitosisReady) || !m.alive) continue;

            float ang = randf() * M_PI * 2.0f;
            float off = m.radius * 1.85f * m.size;
            simd_float3 dp = {
                m.position.x + cosf(ang) * off,
                m.position.y,
                m.position.z + sinf(ang) * off
            };
            // If the candidate daughter would land outside the dish, try a few
            // alternate angles before giving up so the parent doesn't silently
            // get its divisionPending flag cleared with no offspring produced.
            int tries = 0;
            while (sqrtf(dp.x*dp.x + dp.z*dp.z) > SCENE_BOUND - 3.0f && tries < 8) {
                ang = randf() * M_PI * 2.0f;
                dp.x = m.position.x + cosf(ang) * off;
                dp.z = m.position.z + sinf(ang) * off;
                tries++;
            }
            if (sqrtf(dp.x*dp.x + dp.z*dp.z) > SCENE_BOUND - 3.0f) {
                // Can't find a valid spot — keep the parent's pending flag so
                // it retries next frame rather than losing its division.
                m.program.logEvent(bioTime, 3, "DIV_REFUSED", "dish_edge");
                continue;
            }
            m.divisionPending = false;

            const SimCell parentSnapshot = m;
            int parentUid = m.cellUid;
            SplitStats stats;
            stats.cytoplasmicRatioA = 0.5f;
            float sf = 0.794f;
            simd_float3 pushVel = {cosf(ang) * 1.2f, 0, sinf(ang) * 1.2f};

            SimCell d = deriveDaughter(parentSnapshot, dp, pushVel, sf, stats, false);
            if (mode != MODE_SINGLE_CELL) {
                // Add small mutation to daughter's traits (colony evolution)
                d.glycolysisBias = clampf(d.glycolysisBias + (randf() - 0.5f) * GENOME_MUTATION_RATE * 2, 0.3f, 2.2f);
                d.prolifBias = clampf(d.prolifBias + (randf() - 0.5f) * GENOME_MUTATION_RATE * 2, 0.3f, 2.2f);
                d.rosTolerance = clampf(d.rosTolerance + (randf() - 0.5f) * GENOME_MUTATION_RATE * 2, 0.3f, 2.2f);
                d.repairRate = clampf(d.repairRate + (randf() - 0.5f) * GENOME_MUTATION_RATE * 2, 0.2f, 2.0f);
                d.cloneId = (fabsf(d.glycolysisBias - parentSnapshot.glycolysisBias) > 0.08f)
                    ? nextCloneId++
                    : parentSnapshot.cloneId;
            }

            m = deriveDaughter(parentSnapshot, parentSnapshot.position,
                               {-cosf(ang) * 0.5f, 0, -sinf(ang) * 0.5f},
                               sf, stats, true);
            m.cellUid = parentUid; // parent keeps UID

            pendingDivisions.push_back({parentUid, m.cellUid, d.cellUid, m.position, d.position});
            daughters.push_back(d);
            statDivisions++;
            m.program.logEvent(bioTime, 0, "DIV_FINALIZE");
        }

        for (auto& d : daughters) cells.push_back(d);
    }

    void updateStats() {
        statAlive=0;statProlif=0;statQuiescent=0;statApoptotic=0;statNecrotic=0;statGlycolytic=0;
        float sumATP=0; statPhases[0]=statPhases[1]=statPhases[2]=statPhases[3]=0;
        for(auto& c:cells){if(!c.alive) continue;
            statAlive++; sumATP+=c.ATP;
            if(c.fate==SIM_FATE_PROLIF) statProlif++;
            else if(c.fate==SIM_FATE_QUIESCENT) statQuiescent++;
            else if(c.fate==SIM_FATE_APOPTOTIC) statApoptotic++;
            if(c.necrotic) statNecrotic++;
            if(c.glycolytic) statGlycolytic++;
            if(c.phase>=0&&c.phase<=3) statPhases[c.phase]++;
        }
        statAvgATP=statAlive>0?sumATP/statAlive:0;
        // Drug stats
        if (preDrugCount > 0) statViability = (float)statAlive / (float)preDrugCount * 100.0f;
        float sumDrugDmg = 0;
        for (auto& c : cells) if (c.alive) sumDrugDmg += c.drugDamage;
        statAvgDrugDamage = statAlive > 0 ? sumDrugDmg / statAlive : 0;
    }
};
