#pragma once

#include "../core/Constants.h"
#include "CellCycleProgram.h"
#include "MediumField.h"
#include "biochem/Apoptosis.h"
#include "biochem/Bioagent.h"
#include "biochem/TargetLibrary.h"
#include "biochem/BindingMatcher.h"
#include "biochem/Virion.h"
#include "biochem/Bacterium.h"
#include "biochem/PathogenRegistry.h"
#include "biochem/PathogenKinetics.h"
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

// Drug PK/PD model + library removed 2026-04-19 (pending rewrite).

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

        // CDK rates tuned for a 20 bio-h CTC-HeLa cycle (Fluo-N2DL-HeLa
        // 125→674 in 46 h = ~20 h doubling). Rate constants bumped 1.5×
        // over the 24 h tuning so G1/S/G2 transits complete on the
        // faster schedule.
        float dCycD = (gs*1.95f - CycD*(0.90f+(1-gs)*0.45f)) * dt_bio;
        float CDK4act = CycD*(1-p21e*0.5f);
        float CDK2Eact = CycE*(1-p21e);
        float dRb = (0.38f*(1-Rb) - (CDK4act*2.70f+CDK2Eact*1.80f)*Rb) * dt_bio;
        float RbP = 1-Rb;
        float dE2F = (RbP*2.25f*(1+E2F*1.2f) - (Rb*1.80f+0.45f)*E2F) * dt_bio;
        float dCycE = (E2F*2.00f*(1-CycA) - CycE*(0.65f+CycA*2.45f)) * dt_bio * (1-p21e*0.8f);
        float APC_Cdh1 = fmaxf(0, 1-(CDK2Eact+CycA)*1.2f);
        float APC_Cdc20 = fmaxf(0, CycB-0.25f)*2.0f;
        float dCycA = (E2F*1.80f*fmaxf(0,CycE-0.12f) - CycA*(0.22f+APC_Cdh1*1.60f+APC_Cdc20)) * dt_bio * (1-p21e*0.6f);
        // Cdc25 + CycB synthesis bumped 1.5× to cross the M-entry
        // threshold in ~3 bio-h of G2 on the 20 h CTC schedule.
        float Cdc25 = fmaxf(0, CycA-0.20f)*2.0f;
        float dCycB = (2.50f*Cdc25*(1+CycB*0.8f) - CycB*(0.18f+APC_Cdc20*2.5f)) * dt_bio * (1-p21e*0.3f);
        float dp21 = -p21*0.18f*dt_bio;

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

    // ── Phase G2: TP53 ↔ MDM2 stress-response loop ─────────────────────
    // Units: normalised protein concentration (arbitrary; baseline ≈ 0.05
    // represents "healthy cell" MDM2-dominated steady state). Calibrated
    // so that sustained high DNA damage produces ~5 h oscillations
    // matching Purvis 2012 ATM-driven p53 pulsing, with protein counts
    // staying in [0, ~2.0] range.
    //
    // Rate constants:
    //   k_p53_basal    — basal p53 synthesis (protein/s)
    //   k_mdm2_from_p53— p53-driven MDM2 transcription (per p53·s)
    //   k_mdm2_deg     — MDM2 basal degradation (1/s, ≈ 30 min t½
    //                     Momand 1992)
    //   k_p53_by_mdm2  — MDM2-mediated p53 ubiquitination + proteasome
    //                     degradation rate (per MDM2·s). K_d = 200 nM
    //                     for p53↔MDM2 binding (Picksley 1994, Bottger
    //                     1997); under normal conditions this keeps p53
    //                     t½ ≈ 30 min (Haupt 1997, Kubbutat 1997).
    //   atm_stabilize_k— ATM activates at damageLevel > ATM_DAMAGE_K;
    //                     ATM phosphorylates p53-Ser15 which blocks
    //                     MDM2 binding → p53 accumulates.
    //
    // Downstream: p53 transcriptionally induces p21 (CDKN1A) which
    // flows into the existing cdk.p21 pool, pushing cells into G1
    // arrest — the central mechanism of the DNA-damage checkpoint.
    //
    // Sources:
    //   Purvis TJ et al. 2012 Science 336:1440 — p53 pulsing dynamics
    //   Haupt Y et al. 1997 Nature 387:296   — MDM2 targets p53 for deg.
    //   Kubbutat MH et al. 1997 Nature 387:299 — same
    //   Picksley SM et al. 1994 Oncogene 9:2523 — p53↔MDM2 Kd
    //   Bottger V et al. 1997 J Mol Biol 269:744 — refined Kd
    //   Momand J et al. 1992 Cell 69:1237 — MDM2 half-life
    //   Lev Bar-Or R et al. 2000 PNAS 97:11250 — early oscillator model
    float p53_protein = 0.089f;
    float MDM2_protein = 0.21f;
    // Phase G6.2: Nutlin-3 / MDM2-antagonist orthogonal perturbation.
    // Vassilev 2004 Science 303:844 showed Nutlin-3 binds the p53-
    // binding pocket of MDM2, blocking ubiquitination of p53 WITHOUT
    // DNA damage. This flag simulates saturating Nutlin binding:
    // effective MDM2-mediated p53 degradation is set to zero (i.e.
    // p53 is no longer cleared even though damageLevel=0). It is
    // the strongest orthogonal test that the p53 axis is driven by
    // the real MDM2 machinery — if p53 rises under Nutlin without
    // ATM activation, the code must be going through the protein
    // loop, not a damage shortcut.
    bool mdm2_inhibited = false;
    // Phase G3: MDM2 mRNA as explicit intermediate between p53
    // transactivation and MDM2 protein synthesis. Required for
    // Purvis-2012-style oscillation; Geva-Zatorsky 2006 showed the
    // 3-variable (p53, mRNA, protein) delay model captures the ~5 h
    // pulsing period observed in single HeLa cells under sustained
    // γ-IR damage. Without this, a 2-variable (p53, MDM2) ODE settles
    // to a stable elevated p53 — which is what G2 delivered.
    float MDM2_mRNA = 0.21f;
    float ATM_active = 0.0f; // 0..1, rises when damageLevel > threshold
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
    // Legacy integer phase/timer. Kept on the struct (preserves field
    // layout for code paths that still reference them) but the new
    // `apo` engine is the source of truth. `apoptosisPhase` 0 → alive,
    // any >0 → dying/dead, per the nine existing call-sites that gate
    // on `c.apoptosisPhase > 0`. Rewritten logic sets this to a value
    // that mirrors `apoPhase` so those checks still work during the
    // transition period.
    int apoptosisPhase; float apoTimer;
    float adaptationTimer;

    // ── Multi-threshold apoptosis (Apoptosis.h engine) ─────────────────
    // Each cell owns its own engine. `apoPhase` is the coarse visual
    // state (ALIVE..CLEARED); `apo` holds the full Bcl-2/Bax/caspase ODE
    // state. See [src/simulation/biochem/Apoptosis.h](src/simulation/biochem/Apoptosis.h)
    // and the `Apoptosis::` block in Constants.h for constants.
    ApoptosisEngine apo;
    Apoptosis::Phase apoPhase = Apoptosis::ALIVE;
    // Mass ledgers — the physics-of-matter principle: every pool that
    // exists in a live cell must go *somewhere* when that cell dies.
    // Units: biomass-equivalent (1 biomass ≈ 200 pg dry mass).
    float membraneMass_bm      = 0.0f;  // set at init to 0.2 × biomass
    float receptorMass_bm      = 0.0f;  // set at init to 0.07 × biomass
    float initialBiomassAtDeath = 0.0f; // snapshot for release math
    float initialMembraneAtDeath = 0.0f;
    float initialReceptorAtDeath = 0.0f;
    // Trigger-input timers (bio-seconds; drive crowding / replicative curves).
    float chronicPressureBioSec = 0.0f; // elapsed under localPressure ≥ 3
    float senescenceBioSec      = 0.0f; // elapsed while senescent == true
    float fasLExposure          = 0.0f; // ng/mL accumulated from DAMPs / drug
    // Lysosomal pools — accept mass from engulfed apoptotic bodies
    // (efferocytosis) and digest it back into this cell's biomass.
    float lysosomalLoad_cyto    = 0.0f;
    float lysosomalLoad_mem     = 0.0f;
    float lysosomalLoad_rec     = 0.0f;
    // DAMP-reception window: seconds since last local lysis event within
    // DAMP_NEIGHBOR_RADIUS_WU. Drives transient stress/ROS bumps.
    float damageSensedBioSec    = 0.0f;
    // Body-spawn latch: main.mm sets this true once it has fragmented
    // this cell into apoptotic bodies (with the cell's remaining mass
    // moved into body ledgers). Without a renderer (headless), it stays
    // false → updateApoptosis takes the mass-release fallback path on
    // the BODIES transition so closed-system conservation still holds.
    bool bodiesSpawned          = false;
    float divisionCooldown = 0; // Prevents immediate re-entry into M-phase after division
    float postDivisionRecovery = 0; // Seconds of post-mitotic settling / render blend

    // ── Bioagent inventory (Phase M2 scaffold) ────────────────────────
    // Intracellular concentration of each foreign molecule that has
    // crossed the membrane, keyed by its position in gBioagents registry.
    // Sparse — only populated on drug-load. Empty vectors cost ~24 bytes
    // per cell when no drugs are present.
    //   drugIntra_uM[bioagentIdx]              — µM of drug inside cytoplasm
    //   targetBound_uM[targetIdx][drugIdx]     — bound drug on this target
    //   targetCount[targetIdx]                 — total copies of target
    //   targetOccupancy[targetIdx]             — 0..1 bound fraction
    // See plan Part-One §4.2 and Part-Five M3.
    std::vector<float> drugIntra_uM;
    std::vector<std::vector<float>> targetBound_uM;
    std::vector<int>   targetCount;
    std::vector<float> targetOccupancy;

    // ── Pathogen inventory (Phase 7 §44) ────────────────────────────
    // Intracellular load of each virion species + bacterium species.
    //   intraVirions[specIdx]    — viral genome copies inside cytosol
    //   intraBacteria[specIdx]   — bacterial load inside host vacuole / cytosol
    //   assembledVirions[idx]    — new virions ready to bud / burst
    //   pathogensReleased        — latch so apoptosis spill runs once
    // Vectors stay empty (≈24 B) if no pathogen ever enters the cell.
    std::vector<float> intraVirions;
    std::vector<float> intraBacteria;
    std::vector<float> assembledVirions;
    bool pathogensReleased = false;

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
        // Phase G2: stress-response pools initialised to deterministic
        // fixed baseline (MDM2-dominated steady state). Per-cell
        // heterogeneity in p53 pulse amplitude (Lahav 2004, Geva-Zatorsky
        // 2006) emerges from the per-cell damageLevel variance set above
        // — not from RNG draws at init time, so the RNG stream and the
        // 5.1% seq02 baseline are preserved bit-for-bit. When G3 wires
        // MDM2 and p53 into real transcription via CentralDogma, intrinsic
        // stochastic gene-expression noise will supply the variance
        // professionally (Elowitz 2002, Raj 2008).
        // Phase G2+G3: baseline set to analytical steady-state of the
        // 3-variable (p53, MDM2_mRNA, MDM2_protein) ODE with
        // damageLevel≈0. All three sit at 0.21 (except p53 at 0.089)
        // under the tuned rate constants: since k_translate == k_prot_deg,
        // MDM2_protein_ss = MDM2_mRNA_ss at steady state. Initialising
        // exactly at the fixed point avoids a burn-in transient and
        // preserves the 5.1% seq02 calibration bit-for-bit.
        p53_protein  = 0.089f;
        MDM2_mRNA    = 0.21f;
        MDM2_protein = 0.21f;
        ATM_active   = 0.0f;
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
        // Strongly bias initial cycleTimer toward LATE-IN-PHASE
        // (cube-root distribution) so the seeded population looks like a
        // steady-state mid-log-phase culture. Real HeLa at imaging start
        // had been growing for days; many cells are close to their next
        // division. Uniform random under-represents late-phase cells.
        float lateBias = powf(randf(), 1.0f/3.0f);
        if (phase == 0) {
            cycleTimer           = lateBias * g1end;
            checkpointG1Passed   = false;
            checkpointG2Passed   = false;
        } else if (phase == 1) {
            cycleTimer           = g1end + lateBias * CYCLE_S_DUR;
            checkpointG1Passed   = true;
            checkpointG2Passed   = false;
        } else if (phase == 2) {
            cycleTimer           = send + lateBias * CYCLE_G2_DUR;
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
        // Apoptosis engine reset, plus membrane + receptor mass pools
        // allocated proportionally to biomass. See Apoptosis:: block.
        apo.init();
        apoPhase = Apoptosis::ALIVE;
        membraneMass_bm = Apoptosis::MEMBRANE_MASS_PER_BIOMASS * biomass;
        receptorMass_bm = Apoptosis::RECEPTOR_MASS_PER_BIOMASS * biomass;
        initialBiomassAtDeath = 0.0f;
        initialMembraneAtDeath = 0.0f;
        initialReceptorAtDeath = 0.0f;
        chronicPressureBioSec = 0.0f;
        senescenceBioSec = 0.0f;
        fasLExposure = 0.0f;
        lysosomalLoad_cyto = 0.0f;
        lysosomalLoad_mem = 0.0f;
        lysosomalLoad_rec = 0.0f;
        damageSensedBioSec = 0.0f;
        bodiesSpawned = false;
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
    // p21 threshold raised 0.35 → 0.70. HeLa cells have non-zero basal
    // p21 that fluctuates; only strong p53-induced p21 (~0.7+) genuinely
    // stalls mitosis entry. The old 0.35 threshold blocked routine
    // cycling once p21 accumulated even slightly above baseline.
    if (c.cdk.p21 >= 0.70f)       return {false, "p21 high (p53 / DDR)"};
    if (c.damageLevel >= 0.40f)   return {false, "DSBs above threshold"};
    if (c.ATP <= 20.0f)           return {false, "ATP insufficient for mitosis"};
    // biomass gate relaxed from 1.70 → 1.20. Real HeLa doubles protein
    // mass over the cycle but daughter biomass at mitosis entry need
    // only be ~60 % of parent M-entry biomass for clean division; the
    // strict 1.70 threshold prevented daughters from ever re-entering M.
    if (c.biomass <= 1.20f)       return {false, "biomass < 1.20"};
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

    // ── Free-medium pathogen pools (Phase 7 §44.7) ──────────────────
    // Virions and bacteria drifting in the fluid. They bind nearby live
    // cells via shape-compatibility scoring; upon host lysis fresh
    // pathogens are appended here. Inoculation enters through
    // `seedVirion()` / `seedBacterium()`.
    std::vector<Virion>    gFreeVirions;
    std::vector<Bacterium> gFreeBacteria;
    uint32_t nextPathogenUid = 1;

    // Telemetry counters populated each tick.
    int statVirionsFree    = 0;
    int statVirionsBound   = 0;
    int statVirionsIntra   = 0;
    int statBacteriaFree   = 0;
    int statBacteriaIntra  = 0;
    int statInfectedCells  = 0;

    // Find the cell nearest to (cx, cz) that we'll inoculate around.
    // Returns nullptr if no cells exist.
    const SimCell* cellNearest(float cx, float cz) const {
        const SimCell* best = nullptr;
        float bestD2 = 1e30f;
        for (const auto& c : cells) {
            if (!c.alive) continue;
            float dx = c.position.x - cx;
            float dz = c.position.z - cz;
            float d2 = dx*dx + dz*dz;
            if (d2 < bestD2) { bestD2 = d2; best = &c; }
        }
        return best;
    }

    // Inoculate N virions in an annulus OUTSIDE the nearest cell's visual
    // radius so they actually drift in from the fluid rather than popping
    // in through the membrane. Y matches the cell's Y so they sit in the
    // same slab as the host. spreadRadius ignored (back-compat); the
    // annulus thickness is a fixed multiple of cell radius.
    void seedVirion(int specIdx, float cx, float cz, int count, float /*spreadRadius*/) {
        if (specIdx < 0 || specIdx >= (int)gPathogens().virionSpecs.size()) return;
        const SimCell* host = cellNearest(cx, cz);
        float hy = host ? host->position.y : (FLOOR_Y + 2.0f);
        float hr = host ? fmaxf(1.0f, host->radius * host->size) : 2.0f;
        // Annulus: [hr + 0.5, hr + 3.0] wu outside the membrane.
        float rMin = hr + 0.5f;
        float rMax = hr + 3.0f;
        for (int i = 0; i < count; i++) {
            Virion v;
            v.specIdx = specIdx;
            float ang = 6.2831853f * ((float)rand() / (float)RAND_MAX);
            float u   = (float)rand() / (float)RAND_MAX;
            float r   = sqrtf(rMin*rMin + u * (rMax*rMax - rMin*rMin));
            v.position = simd_make_float3(cx + r*cosf(ang), hy,
                                          cz + r*sinf(ang));
            v.velocity = simd_make_float3(0.0f, 0.0f, 0.0f);
            v.stage = VirionStage::FREE;
            v.uid = nextPathogenUid++;
            gFreeVirions.push_back(v);
        }
    }
    void seedBacterium(int specIdx, float cx, float cz, int count, float /*spreadRadius*/) {
        if (specIdx < 0 || specIdx >= (int)gPathogens().bacteriumSpecs.size()) return;
        const BacteriumSpec& bs = gPathogens().bacteriumSpecs[specIdx];
        const SimCell* host = cellNearest(cx, cz);
        float hy = host ? host->position.y : (FLOOR_Y + 2.0f);
        float hr = host ? fmaxf(1.0f, host->radius * host->size) : 2.0f;
        float rMin = hr + 1.0f;
        float rMax = hr + 4.0f;
        for (int i = 0; i < count; i++) {
            Bacterium b;
            b.specIdx = specIdx;
            float ang = 6.2831853f * ((float)rand() / (float)RAND_MAX);
            float u   = (float)rand() / (float)RAND_MAX;
            float r   = sqrtf(rMin*rMin + u * (rMax*rMax - rMin*rMin));
            b.position = simd_make_float3(cx + r*cosf(ang), hy,
                                          cz + r*sinf(ang));
            b.biomass_bm = bs.biomass_bm;
            float a2 = 6.2831853f * ((float)rand() / (float)RAND_MAX);
            b.runDir = simd_make_float3(cosf(a2), 0.0f, sinf(a2));
            b.uid = nextPathogenUid++;
            gFreeBacteria.push_back(b);
        }
    }

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
        // Daughter starts partway through G1 (20 % of G1 elapsed).
        // Represents the fact that immediately after cytokinesis, HeLa
        // cells already have considerable CycD / protein synthesis
        // machinery inherited from the parent. Zero-init forces daughters
        // to "cold start" their cycle, which artificially slows growth.
        // Daughter starts 20 % through G1 — HeLa daughters inherit
        // considerable CycD / protein machinery from the parent, so a
        // zero-init forces an artificial "cold start" lag.
        d.cycleTimer = 0.20f * CYCLE_G1_DUR; d.cycleProgress = 0;
        d.checkpointG1Passed = false; d.checkpointG2Passed = false;
        d.g1WaitTimer = 0; d.g2WaitTimer = 0;
        d.divisionPending = false;
        // Cooldown (sdt-units) before a fresh daughter can re-attempt
        // mitosis. Colony cooldown 8.0 → 0.5 because daughters otherwise
        // pile up in G2 with cycleTimer pinned at g2end × 0.99 for
        // ~14 bio-h, blocking the observed CTC doubling rates.
        d.divisionCooldown = singleMode ? 1.5f : 0.5f;
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
        // Multi-threshold apoptosis — daughter starts fresh. Membrane
        // and receptor mass scale with the daughter's new biomass.
        d.apo.init();
        d.apoPhase = Apoptosis::ALIVE;
        d.membraneMass_bm = Apoptosis::MEMBRANE_MASS_PER_BIOMASS * d.biomass;
        d.receptorMass_bm = Apoptosis::RECEPTOR_MASS_PER_BIOMASS * d.biomass;
        d.initialBiomassAtDeath = 0.0f;
        d.initialMembraneAtDeath = 0.0f;
        d.initialReceptorAtDeath = 0.0f;
        d.chronicPressureBioSec = 0.0f;
        d.senescenceBioSec = 0.0f;
        d.fasLExposure = 0.0f;
        d.lysosomalLoad_cyto = 0.0f;
        d.lysosomalLoad_mem = 0.0f;
        d.lysosomalLoad_rec = 0.0f;
        d.damageSensedBioSec = 0.0f;
        d.bodiesSpawned = false;

        // INHERITED as-is (mutations accumulate across generations):
        // glycolysisBias, prolifBias, rosTolerance, repairRate
        return d;
    }

    // ── Chemistry-first drug system (MOA-free) ────────────────────────
    // List of drugs currently present in the dish. Each entry holds a
    // registry index, dish-level µM concentration, and per-cell-target
    // Kd precomputed by BindingMatcher. Empty until applyDrug() fires.
    struct AppliedDrug {
        int    entityIdx = -1;     // into gBioagents
        float  dishConc_uM = 0.0f; // applied concentration
        std::vector<BindingAffinity> affinityPerTarget; // sized = gTargets.count()
    };
    std::vector<AppliedDrug> appliedDrugs;

    // Resize per-cell bioagent vectors after a drug is added, so indices
    // stay stable.
    void syncCellBioagentInventory() {
        int nDrugs = (int)appliedDrugs.size();
        int nTarg  = gTargets.count();
        for (auto& c : cells) {
            c.drugIntra_uM.resize(nDrugs, 0.0f);
            c.targetBound_uM.assign(nTarg, std::vector<float>(nDrugs, 0.0f));
            c.targetCount.assign(nTarg, 100000);     // default target abundance
            c.targetOccupancy.assign(nTarg, 0.0f);
        }
    }

    // Apply a drug uniformly across the dish. drugId must be a registry id.
    void applyDrug(const std::string& drugId, float conc_uM) {
        const ChemicalEntity* drug = gBioagents.get(drugId);
        if (!drug) {
            printf("[Drug] unknown drug '%s'\n", drugId.c_str());
            return;
        }
        // Find the index in gBioagents.
        int idx = -1;
        for (int i = 0; i < (int)gBioagents.all().size(); i++) {
            if (gBioagents.all()[i].id == drugId) { idx = i; break; }
        }
        if (idx < 0) return;

        AppliedDrug ad;
        ad.entityIdx = idx;
        ad.dishConc_uM = conc_uM;
        // Precompute drug × target affinity once.
        ad.affinityPerTarget.resize(gTargets.count());
        for (int t = 0; t < gTargets.count(); t++) {
            ad.affinityPerTarget[t] = BindingMatcher::score(*drug, gTargets.at(t).profile);
        }
        appliedDrugs.push_back(std::move(ad));
        syncCellBioagentInventory();

        printf("[Drug] Applied %s at %.3f µM. Top affinities:\n",
               drugId.c_str(), conc_uM);
        for (int t = 0; t < gTargets.count(); t++) {
            const auto& aff = appliedDrugs.back().affinityPerTarget[t];
            if (aff.score > 0.30f)
                printf("  %-24s  Kd=%.3f mM  score=%.2f\n",
                       gTargets.idAt(t).c_str(), aff.Kd_mM, aff.score);
        }
    }

    // One-tick update of drug pharmacokinetics + target binding +
    // function modulators. Called from the main update loop.
    void updateDrugPK(float dt) {
        if (appliedDrugs.empty()) return;
        float dt_biosec = dt * SLOW_DT_SCALE * 3600.0f;
        for (auto& c : cells) {
            if (!c.alive) continue;
            for (int di = 0; di < (int)appliedDrugs.size(); di++) {
                AppliedDrug& ad = appliedDrugs[di];
                const ChemicalEntity& d = gBioagents.all()[ad.entityIdx];
                // Lipinski-based permeability: logP > 0 favours entry,
                // large MW + high TPSA hurt. Base rate ≈ 0.02 /bio-s.
                float perm = 0.02f * expf(0.3f * (d.logP - 1.0f))
                                   * expf(-d.mw / 700.0f)
                                   * expf(-d.tpsa / 150.0f);
                perm = fmaxf(1e-5f, fminf(0.2f, perm));
                // Flux from dish to cytoplasm (simple linear, µM/s).
                float dC = perm * (ad.dishConc_uM - c.drugIntra_uM[di]);
                c.drugIntra_uM[di] = fmaxf(0.0f,
                    c.drugIntra_uM[di] + dC * dt_biosec);

                // Bind each target: dB/dt = k_on * free_drug * free_target - k_off * bound
                //
                // Convert target-copy count to µM in a 4 pL HeLa cell:
                //   µM = copies / (NA × V × 1e-6)
                //      = copies / (6.02e23 × 4e-12 × 1e-6)
                //      = copies / 2.41e6
                for (int ti = 0; ti < (int)ad.affinityPerTarget.size(); ti++) {
                    const BindingAffinity& aff = ad.affinityPerTarget[ti];
                    float drugFree = c.drugIntra_uM[di];
                    float bound    = c.targetBound_uM[ti][di];
                    float total    = (float)c.targetCount[ti];
                    float totalUM  = total / 2.41e6f;  // Avogadro-corrected
                    float free_t   = fmaxf(0.0f, totalUM - bound);
                    float onFlux   = aff.k_on_per_uM  * drugFree * free_t;
                    float offFlux  = aff.k_off_per_s * bound;
                    float delta    = (onFlux - offFlux) * dt_biosec;
                    c.targetBound_uM[ti][di] = fmaxf(0.0f, bound + delta);
                    if (totalUM > 1e-9f) {
                        c.targetOccupancy[ti] = fminf(1.0f,
                            c.targetBound_uM[ti][di] / totalUM);
                    }
                }
            }
            // Fire every target's modulator with its combined occupancy
            // (sum across all drugs bound to this target).
            for (int ti = 0; ti < gTargets.count(); ti++) {
                float occ = c.targetOccupancy[ti];
                if (occ > 0.01f && gTargets.at(ti).modulator) {
                    gTargets.at(ti).modulator(c, occ, dt_biosec);
                }
            }
        }
    }

    // ── Pathogen one-tick update (Phase 7 §44) ──────────────────────
    //
    // Discipline:
    //   * No tropism / MOA flag is consulted anywhere. A virion binds
    //     because its 5-D spike descriptor matches a host receptor
    //     profile (PathogenKinetics::virionBindScore).
    //   * Stage machine is driven entirely by physical conditions
    //     (bound-spike count, elapsed dwell, intracellular load).
    //   * On host apoptosis we do NOT drive it here — the existing
    //     apoptosis engine transitions to BODIES; our spill hook
    //     (see releasePathogens on BODIES entry) runs once.
    //
    // Called once per inner sub-step after updateDrugPK.
    void updatePathogens(float dt) {
        if (gFreeVirions.empty() && gFreeBacteria.empty()) {
            // Still iterate cells to drive intracellular replication
            // when no free pool exists.
            bool anyIntra = false;
            for (auto& c : cells) {
                if (!c.intraVirions.empty() || !c.intraBacteria.empty()) { anyIntra = true; break; }
            }
            if (!anyIntra) return;
        }
        const float dt_biosec = dt * SLOW_DT_SCALE * 3600.0f;

        // -- 0. Advance stage timers for all virions (needed for BOUND→ENTERING etc.) --
        for (auto& v : gFreeVirions) v.stageTimer_s += dt_biosec;
        for (auto& b : gFreeBacteria) b.stageTimer_s += dt_biosec;
        // -- 1. Free-medium transport (Brownian + simple advection) --
        // Brownian: σ = sqrt(2D·dt) per axis. For a ~100 nm virion in
        // water at 37 °C, Stokes-Einstein D ≈ 4 µm²/s = 0.16 wu²/s
        // (MICROMETERS_PER_WORLD_UNIT = 5). We use σ per bio-s in wu.
        const float BROWN_SIGMA = 0.40f;   // wu · bio-s^(-1/2)
        // Contact sampling: ~2 × median cell radius (HeLa ≈ 10 µm ≈ 2 wu)
        // so a virion drifting in the cell's vicinity sees it.
        const float CONTACT_RADIUS = 6.0f;
        for (auto& v : gFreeVirions) {
            if (v.stage != VirionStage::FREE) continue;
            // Proper random walk: displacement per step ~ σ · sqrt(dt).
            float u1 = (float)rand() / (float)RAND_MAX;
            float u2 = (float)rand() / (float)RAND_MAX;
            // Box-Muller → standard normal.
            float r0 = sqrtf(-2.0f * logf(fmaxf(1e-12f, u1)));
            float g1 = r0 * cosf(6.2831853f * u2);
            float g2 = r0 * sinf(6.2831853f * u2);
            float step = BROWN_SIGMA * sqrtf(dt_biosec);
            v.position.x += g1 * step;
            v.position.z += g2 * step;
        }
        for (auto& b : gFreeBacteria) {
            if (b.stage != BacteriumStage::FREE) continue;
            // Run-and-tumble. Tumble every ~1 bio-s on avg.
            b.tumbleTimer_s += dt_biosec;
            if (b.tumbleTimer_s > 1.0f) {
                float a = 6.2831853f * ((float)rand() / (float)RAND_MAX);
                b.runDir = simd_make_float3(cosf(a), 0.0f, sinf(a));
                b.tumbleTimer_s = 0.0f;
            }
            float speed = 0.6f; // wu/bio-s
            b.position += b.runDir * speed * dt_biosec;
            b.divisionClock_s += dt_biosec;
        }

        // -- 2. Binding step: for each free virion / bacterium, find
        //       closest live cell within CONTACT_RADIUS, compute the
        //       shape-compatibility score, accumulate bound spikes.
        for (auto& v : gFreeVirions) {
            if (v.stage != VirionStage::FREE) continue;
            if (v.specIdx < 0) continue;
            const VirionSpec& vs = gPathogens().virionSpecs[v.specIdx];
            int bestIdx = -1;
            float bestD2 = CONTACT_RADIUS * CONTACT_RADIUS;
            for (int i = 0; i < (int)cells.size(); i++) {
                SimCell& c = cells[i];
                // Any cell that hasn't yet fragmented is a viable host:
                // real viruses/bacteria happily infect PRIMED and MOMP
                // cells, and sometimes accelerate their demise. Only
                // cells already in BODIES/CLEARED are skipped.
                if (!c.alive || c.apoPhase >= Apoptosis::FRAGMENTATION) continue;
                float dx = c.position.x - v.position.x;
                float dz = c.position.z - v.position.z;
                float d2 = dx*dx + dz*dz;
                if (d2 < bestD2) { bestD2 = d2; bestIdx = i; }
            }
            if (bestIdx < 0) continue;
            std::string rcvId;
            float score = PathogenKinetics::virionBindScore(vs, rcvId);
            // Hard gate: a virion whose spike doesn't match any host
            // receptor profile simply bounces — no productive contact.
            if (score < 0.35f) continue;
            // A binding event progresses mass-action style:
            //   d[bound]/dt = k · spikes_total · (1 − bound_frac)
            float k_eff = PathogenKinetics::scoreToOnRate(score);
            float capacity = (float)vs.spikesPerVirion;
            float growth = k_eff * (capacity - v.boundSpikes) * dt_biosec;
            v.boundSpikes = fminf(capacity, v.boundSpikes + growth);
            // Binding threshold → entry. At least 25 % of spikes engaged
            // AND at least 20 absolute contacts. Real flu HA needs ~10-20
            // engaged spikes for fusion (Matlin 1982, Imai 2009).
            float threshold = fmaxf(20.0f, capacity * 0.25f);
            if (v.boundSpikes >= threshold) {
                v.stage = VirionStage::BOUND;
                v.hostCellIdx = bestIdx;
                v.stageTimer_s = 0.0f;
            }
        }
        for (auto& b : gFreeBacteria) {
            if (b.stage != BacteriumStage::FREE) continue;
            if (b.specIdx < 0) continue;
            const BacteriumSpec& bs = gPathogens().bacteriumSpecs[b.specIdx];
            int bestIdx = -1;
            float bestD2 = CONTACT_RADIUS * CONTACT_RADIUS;
            for (int i = 0; i < (int)cells.size(); i++) {
                SimCell& c = cells[i];
                // Any cell that hasn't yet fragmented is a viable host:
                // real viruses/bacteria happily infect PRIMED and MOMP
                // cells, and sometimes accelerate their demise. Only
                // cells already in BODIES/CLEARED are skipped.
                if (!c.alive || c.apoPhase >= Apoptosis::FRAGMENTATION) continue;
                float dx = c.position.x - b.position.x;
                float dz = c.position.z - b.position.z;
                float d2 = dx*dx + dz*dz;
                if (d2 < bestD2) { bestD2 = d2; bestIdx = i; }
            }
            if (bestIdx < 0) {
                // Extracellular division capped at a hard global ceiling
                // so the dish can't runaway-fill with bacteria. Each
                // generation doubles — 2 → 4 → 8 → … — so we cap at
                // MAX_FREE_BACTERIA and just skip division once reached.
                constexpr int MAX_FREE_BACTERIA = 2000;
                if ((int)gFreeBacteria.size() >= MAX_FREE_BACTERIA) continue;
                if (b.divisionClock_s > bs.doublingTime_bioSec) {
                    Bacterium nb = b;
                    nb.uid = nextPathogenUid++;
                    nb.divisionClock_s = 0.0f;
                    b.divisionClock_s = 0.0f;
                    float a = 6.2831853f * ((float)rand() / (float)RAND_MAX);
                    nb.position.x += 0.5f * cosf(a);
                    nb.position.z += 0.5f * sinf(a);
                    gFreeBacteria.push_back(nb); // safe: appended at end
                }
                continue;
            }
            std::string rcvId;
            float score = PathogenKinetics::bacteriumBindScore(bs, rcvId);
            // Same strict gate as virions — no score, no productive contact.
            if (score < 0.35f) continue;
            float k_eff = PathogenKinetics::scoreToOnRate(score);
            (void)k_eff;
            b.stage = BacteriumStage::ADHERING;
            b.hostCellIdx = bestIdx;
            b.stageTimer_s = 0.0f;
        }

        // -- 3. Entry → replication transitions --
        for (auto& v : gFreeVirions) {
            if (v.stage == VirionStage::BOUND && v.stageTimer_s > 60.0f) {
                v.stage = VirionStage::ENTERING;
                v.stageTimer_s = 0.0f;
            }
            if (v.stage == VirionStage::ENTERING && v.stageTimer_s > 300.0f) {
                // Deposit one genome into host intracellular pool.
                if (v.hostCellIdx >= 0 && v.hostCellIdx < (int)cells.size()) {
                    SimCell& h = cells[v.hostCellIdx];
                    int nSpec = (int)gPathogens().virionSpecs.size();
                    if ((int)h.intraVirions.size() < nSpec) h.intraVirions.resize(nSpec, 0.0f);
                    if ((int)h.assembledVirions.size() < nSpec) h.assembledVirions.resize(nSpec, 0.0f);
                    h.intraVirions[v.specIdx] += 1.0f;
                }
                v.stage = VirionStage::UNCOATING;  // effectively retired from free list
            }
        }
        for (auto& b : gFreeBacteria) {
            if (b.stage == BacteriumStage::ADHERING && b.stageTimer_s > 120.0f) {
                b.stage = BacteriumStage::INVADED;
                b.stageTimer_s = 0.0f;
                if (b.hostCellIdx >= 0 && b.hostCellIdx < (int)cells.size()) {
                    SimCell& h = cells[b.hostCellIdx];
                    int nSpec = (int)gPathogens().bacteriumSpecs.size();
                    if ((int)h.intraBacteria.size() < nSpec) h.intraBacteria.resize(nSpec, 0.0f);
                    h.intraBacteria[b.specIdx] += 1.0f;
                }
            }
        }

        // -- 4. Intracellular replication + host stress --
        int nVSpec = (int)gPathogens().virionSpecs.size();
        int nBSpec = (int)gPathogens().bacteriumSpecs.size();
        for (auto& c : cells) {
            if (!c.alive) continue;
            if ((int)c.intraVirions.size()   < nVSpec) c.intraVirions.resize(nVSpec, 0.0f);
            if ((int)c.assembledVirions.size()< nVSpec) c.assembledVirions.resize(nVSpec, 0.0f);
            if ((int)c.intraBacteria.size()  < nBSpec) c.intraBacteria.resize(nBSpec, 0.0f);
            bool infected = false;
            for (int s = 0; s < nVSpec; s++) {
                if (c.intraVirions[s] <= 0.0f) continue;
                infected = true;
                const VirionSpec& vs = gPathogens().virionSpecs[s];
                // Autocatalytic replication (clamped by available machinery).
                float dG = vs.replicationRate_per_s * c.intraVirions[s] * dt_biosec;
                // Competition for ribosomes: slowed by translationRateMul.
                c.intraVirions[s] += dG;
                // Assembly once threshold crossed: 10 % of genome pool
                // converts to completed virions per bio-s.
                if (c.intraVirions[s] > vs.assemblyThreshold) {
                    float dA = 0.002f * c.intraVirions[s] * dt_biosec;
                    c.intraVirions[s]    -= dA;
                    c.assembledVirions[s] += dA;
                }
                // Pro-apoptotic pressure via existing drug_pro_apop hook.
                c.apo.state.p53_active = fminf(1.0f,
                    c.apo.state.p53_active
                      + vs.cytotoxicity_per_copy * c.intraVirions[s] * dt_biosec / 3600.0f);
                c.stress = fminf(100.0f,
                    c.stress + vs.cytotoxicity_per_copy * c.intraVirions[s] * dt_biosec * 0.5f);
                // Continuous budding: export a fraction of assembled pool.
                constexpr int MAX_FREE_VIRIONS = 3000;
                if (c.assembledVirions[s] > 1.0f
                    && (int)gFreeVirions.size() < MAX_FREE_VIRIONS) {
                    float bud = vs.budRate_per_s * c.assembledVirions[s] * dt_biosec;
                    bud = fminf(bud, c.assembledVirions[s]);
                    c.assembledVirions[s] -= bud;
                    // Spawn a handful of free virions near this cell.
                    int nBud = (int)bud;
                    for (int k = 0; k < nBud && k < 10; k++) {
                        Virion nv;
                        nv.specIdx = s;
                        float a = 6.2831853f * ((float)rand() / (float)RAND_MAX);
                        float rr = c.radius * 1.05f + 0.2f * ((float)rand() / (float)RAND_MAX);
                        nv.position = simd_make_float3(
                            c.position.x + rr * cosf(a), 0.0f,
                            c.position.z + rr * sinf(a));
                        nv.stage = VirionStage::FREE;
                        nv.uid = nextPathogenUid++;
                        gFreeVirions.push_back(nv);
                    }
                }
            }
            for (int s = 0; s < nBSpec; s++) {
                if (c.intraBacteria[s] <= 0.0f) continue;
                infected = true;
                const BacteriumSpec& bs = gPathogens().bacteriumSpecs[s];
                // Intracellular bacterial replication.
                float ratePerS = 0.693f / bs.doublingTime_bioSec;
                c.intraBacteria[s] += ratePerS * c.intraBacteria[s] * dt_biosec;
                // Toxin secretion into host cytosol → damageLevel.
                c.damageLevel = fminf(2.0f,
                    c.damageLevel + bs.toxinRate_per_s * bs.toxinCytotoxicity
                                    * c.intraBacteria[s] * dt_biosec);
                c.stress = fminf(100.0f,
                    c.stress + 0.01f * bs.toxinRate_per_s * c.intraBacteria[s] * dt_biosec);
            }
            if (infected) statInfectedCells++;
        }
    }

    // ── Pathogen spill hook: called by updateApoptosis when the host
    //    transitions into BODIES. Dumps intracellular virions + bacteria
    //    into the free-medium pools at the host's position.
    void releasePathogens(SimCell& c) {
        if (c.pathogensReleased) return;
        c.pathogensReleased = true;
        int nVSpec = (int)gPathogens().virionSpecs.size();
        int nBSpec = (int)gPathogens().bacteriumSpecs.size();
        for (int s = 0; s < nVSpec && s < (int)c.assembledVirions.size(); s++) {
            float load = c.assembledVirions[s];
            if (load < 1e-3f) continue;
            int nRelease = (int)fminf(load, 2000.0f);
            for (int k = 0; k < nRelease; k++) {
                Virion nv;
                nv.specIdx = s;
                float a = 6.2831853f * ((float)rand() / (float)RAND_MAX);
                float rr = c.radius * (1.0f + 1.5f * ((float)rand() / (float)RAND_MAX));
                nv.position = simd_make_float3(
                    c.position.x + rr * cosf(a), 0.0f,
                    c.position.z + rr * sinf(a));
                nv.stage = VirionStage::FREE;
                nv.uid = nextPathogenUid++;
                gFreeVirions.push_back(nv);
            }
            c.assembledVirions[s] = 0.0f;
            c.intraVirions[s]     = 0.0f;
        }
        for (int s = 0; s < nBSpec && s < (int)c.intraBacteria.size(); s++) {
            float load = c.intraBacteria[s];
            if (load < 0.5f) continue;
            int nRelease = (int)fminf(load, 500.0f);
            for (int k = 0; k < nRelease; k++) {
                Bacterium nb;
                nb.specIdx = s;
                float a = 6.2831853f * ((float)rand() / (float)RAND_MAX);
                float rr = c.radius * (1.0f + 1.8f * ((float)rand() / (float)RAND_MAX));
                nb.position = simd_make_float3(
                    c.position.x + rr * cosf(a), 0.0f,
                    c.position.z + rr * sinf(a));
                nb.biomass_bm = gPathogens().bacteriumSpecs[s].biomass_bm;
                float a2 = 6.2831853f * ((float)rand() / (float)RAND_MAX);
                nb.runDir = simd_make_float3(cosf(a2), 0.0f, sinf(a2));
                nb.uid = nextPathogenUid++;
                gFreeBacteria.push_back(nb);
            }
            c.intraBacteria[s] = 0.0f;
        }
    }

    // Mitosis visualization flag (set by renderer when visual mitosis completes)
    bool mitosisVisualizationComplete = false;

    void init() { init(mode); }

    void init(SimMode m) {
        mode = m;
        cells.clear(); nutrients.init(envO2, envGlucose);
        // One-shot: populate the bioagent registry from disk and
        // register built-in targets. Idempotent (returns 0 on repeat).
        static bool s_registryLoaded = false;
        if (!s_registryLoaded) {
            gBioagents.loadFromDisk("data");
            gPathogens().loadVirionsYaml("data/pathogens/virions.yaml");
            gPathogens().loadBacteriaYaml("data/pathogens/bacteria.yaml");
            s_registryLoaded = true;
        }
        gFreeVirions.clear();
        gFreeBacteria.clear();
        appliedDrugs.clear();
        bioTime = 0;
        lastExecutedScaledDt = 0;
        pendingScaledDt = 0;
        statDivisions = 0; statDeaths = 0;

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
                updateMetabolism(c,subDt);
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
            // Drug PK + target binding + modulators — runs after per-cell
            // biology so modulator adjustments take effect next tick.
            updateDrugPK(subDt);
            // Virion + bacterium physics: Brownian drift, shape-based
            // binding, entry, intracellular replication, continuous
            // budding. Apoptotic spill is handled inside updateApoptosis
            // via releasePathogens() on the BODIES transition.
            statInfectedCells = 0;
            updatePathogens(subDt);
            nutrients.diffuse(subDt, envO2, envGlucose);
        }
        processDivisions();
        cells.erase(std::remove_if(cells.begin(),cells.end(),
            [](const SimCell& c){return !c.alive;}), cells.end());
        // Retire virions that have entered a host (stage beyond FREE/BOUND).
        // Keep FREE and BOUND in the live vector so they render and can
        // still bind. Anything beyond ENTERING is consumed by its host.
        gFreeVirions.erase(std::remove_if(gFreeVirions.begin(), gFreeVirions.end(),
            [](const Virion& v){
                return v.stage == VirionStage::UNCOATING
                    || v.stage == VirionStage::CLEARED;
            }), gFreeVirions.end());
        // Bacteria: retire once invaded.
        gFreeBacteria.erase(std::remove_if(gFreeBacteria.begin(), gFreeBacteria.end(),
            [](const Bacterium& b){
                return b.stage == BacteriumStage::INVADED
                    || b.stage == BacteriumStage::LYSED;
            }), gFreeBacteria.end());
        // Recompute pathogen telemetry.
        statVirionsFree = 0; statVirionsBound = 0;
        for (const auto& v : gFreeVirions) {
            if (v.stage == VirionStage::FREE)  statVirionsFree++;
            if (v.stage == VirionStage::BOUND) statVirionsBound++;
        }
        statBacteriaFree = 0;
        for (const auto& b : gFreeBacteria) {
            if (b.stage == BacteriumStage::FREE) statBacteriaFree++;
        }
        statVirionsIntra = 0; statBacteriaIntra = 0;
        for (const auto& c : cells) {
            for (float x : c.intraVirions)    if (x > 0.0f) statVirionsIntra++;
            for (float x : c.intraBacteria)   if (x > 0.0f) statBacteriaIntra++;
        }
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
            // Substrate anchor: cell rests on the floor. spreadFactor is
            // a visual spreading cue only — we do NOT divide baseY by it
            // (that earlier bug sunk cells INTO the floor when spread
            // rose, so the bottom half disappeared beneath the substrate).
            float baseY = FLOOR_Y + a.radius * a.size * 0.85f;
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
            if (s == MS_HPLUS || s == MS_WATER || s == MS_GROWTH_F) continue;
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
        // pH stress only fires below pH 6.7 (normalized 0.60) — HeLa
        // tolerates 6.8–7.8 without stress response. Rate reduced 4× so
        // brief fluctuations don't pin stress at 100.
        float pHS=fmaxf(0,(0.60f-localPH)/0.06f)*1.0f*mdt;
        float atpS=fmaxf(0,(15-c.ATP)*0.03f)*mdt;
        // Recovery rate increased so cells exit stress when conditions
        // improve. Bumped from 1.8→4.0 for PROLIF cells.
        float recovery = (isQuiescent ? 4.0f : 4.0f) * mdt;
        c.stress=clampf(c.stress + atpS + pHS + (1-localO2)*0.10f*mdt - recovery, 0, 100);

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

        // ── ATP danger timer (feeds the multi-threshold cascade) ────
        // No instant fate-lock here — the timer itself is one of eleven
        // inputs to ApoTriggers. updateApoptosis() builds the full
        // trigger set and calls apo.step() with it.
        if(c.ATP<ATP_DANGER_THRESHOLD){
            c.atpDangerTimer+=mdt;
        } else c.atpDangerTimer=fmaxf(0,c.atpDangerTimer-mdt*0.5f);

        // Chronic-crowding timer: accumulates while the cell is
        // squashed (localPressure ≥ 3). Input to the apoptosis cascade.
        if (c.localPressure >= 3.0f) c.chronicPressureBioSec += mdt;
        else c.chronicPressureBioSec = fmaxf(0, c.chronicPressureBioSec - mdt * 0.25f);

        // Replicative-exhaustion timer: accumulates only after the cell
        // has actually entered senescence via telomere shortening.
        if (c.senescent) c.senescenceBioSec += mdt;
        // DAMP exposure (receiver): decay back toward 0 between hits.
        c.damageSensedBioSec = fmaxf(0.0f, c.damageSensedBioSec - mdt);
        c.fasLExposure = fmaxf(0.0f, c.fasLExposure - mdt * 0.01f);

        // ── Biomass growth ──────────────────────────────────────────
        if(!c.necrotic&&c.apoptosisPhase==0) {
            float presP=fmaxf(0,1.0f-c.localPressure*0.12f);
            float syn=BIOMASS_SYNTH_K*(c.ATP/100)*localGlu;
            float deg=BIOMASS_DEGRADE_K*(c.stress/100);
            float growthMult = (c.fate==SIM_FATE_PROLIF||c.fate==SIM_FATE_UNDETERMINED) ? 1.5f : 0.3f;
            c.biomass=clampf(c.biomass+(syn-deg)*mdt*presP*growthMult, 0.4f, 2.3f);
        }

        // ── Mechanical quorum → p27/p21 induction ───────────────────
        // Ref: Delarue 2018 — pressure induces CDK inhibitors at true
        // confluence. Threshold raised 1.5 → 3.0 because Hertz cell-cell
        // contact in a ~300-cell dish gives baseline pressure 1.5–2.5
        // even mid-log-phase. Real HeLa contact arrest kicks in only
        // when cells genuinely can't spread (monolayer confluent).
        if(c.localPressure>3.0f) {
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

        // ── Background tracking-loss / FOV-emigration / baseline apop ──
        // Density-MODULATED: emigration out of a microscope FOV is
        // dominated by motile isolated cells. A densely packed cell
        // cannot migrate, its neighbours block it, so loss rate is
        // effectively ZERO at high `localPressure`. Model:
        //   p = base × exp(-localPressure)
        // At pressure 0 (isolated): full base rate.
        // At pressure 3 (crowded):  ~5 % of base rate.
        // This is what lets seq01 (sparse) and seq02 (dense) both
        // converge below 7.5 % mean |rel_err| with a single parameter.
        if (TRACKING_LOSS_PROB_PER_BIOSEC > 0.0f
            && mode != MODE_SINGLE_CELL
            && c.alive && c.apoptosisPhase == 0) {
            float dt_biosec = dt * SLOW_DT_SCALE * 3600.0f;
            float densityDamp = expf(-fmaxf(0.0f, c.localPressure));
            float p = TRACKING_LOSS_PROB_PER_BIOSEC * densityDamp * dt_biosec;
            if (randf() < p) {
                releaseAllMass(c);
                c.alive = false;
                statDeaths++;
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

    // updateDrugResponse removed 2026-04-19 — pending rewrite.

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
            if (sqrtf(dx*dx + dz*dz) < c.radius * 1.6f) neighborCount++;
        }

        // Confluence arrest: contact inhibition fires only at true confluence
        // (~8 neighbors in hex-packed monolayer), not at 30 % density. Real
        // HeLa Hippo arrest kicks in above ~80 % monolayer coverage, so
        // require 7+ close neighbors AND elevated pressure to trigger.
        bool contactArrested = (neighborCount >= 7 && c.localPressure > 0.4f)
                            || (c.localPressure > 1.5f);

        // Growth signal (Michaelis-Menten)
        float o2mm=localO2/(localO2+0.10f);
        float glumm=localGlu/(localGlu+0.15f);
        float atpmm=(c.ATP/100)/((c.ATP/100)+0.25f);
        float gs=fminf(1.0f, o2mm*0.40f+glumm*0.30f+atpmm*0.20f+0.05f);

        // ── Phase G2: TP53 ↔ MDM2 stress-response loop ─────────────────
        // Runs as a *parallel observable track* alongside the existing
        // damage→p21 shortcut (below). Rates are calibrated such that
        // under sustained damageLevel ≈ 0.3 the p53/MDM2 pair oscillates
        // with period ≈ 4-6 h matching Purvis 2012 ATM-driven HeLa
        // p53 pulsing. At default steady-state (damageLevel ≲ 0.08),
        // p53 and MDM2 sit at their basal level (≈ 0.05-0.15) with
        // vanishing impact on cell-cycle state — hence zero drift on
        // the 5.1% seq02 baseline. A follow-up PR (G3) replaces the
        // direct damage→p21 shortcut on line below with p53-mediated
        // p21 induction once the dynamics are independently validated.
        {
            // ATM activation: sigmoid on damageLevel; threshold 0.10
            // matches γ-H2AX onset at ~1-2 DSB per nucleus (Rothkamm 2003).
            // Time base: 1 sdt-unit ≈ 1 bio-hour given
            //   SLOW_DT_SCALE=0.055 and BIO_TIME_SCALE=180.
            // Rate constants below are therefore "per bio-hour".
            float atm_target =
                1.0f / (1.0f + expf(-20.0f * (c.damageLevel - 0.10f)));
            // First-order relaxation, τ=0.02 sdt (ATM reaches ~0.9
            // within ~15 bio-min; real ATM auto-phos is ~5 min but
            // using 0.02 keeps the activation on the same time-scale
            // as the downstream oscillator period (Purvis ~5.5 h)
            // rather than saturating in one tick).
            const float atm_tau_sdt = 0.02f;
            c.ATM_active += (atm_target - c.ATM_active) * (sdt / atm_tau_sdt);
            c.ATM_active = clampf(c.ATM_active, 0.0f, 1.0f);

            // p53 production (basal) + MDM2-mediated Michaelis degradation.
            // K_d = 0.20 normalised ≙ 200 nM p53↔MDM2 (Picksley 1994,
            // Bottger 1997). atm_shield = 0.70 — Ser15-P p53 retains
            // ~30 % MDM2 vulnerability, which empirically is what lets
            // sustained-damage oscillations emerge rather than a locked-
            // on plateau (Loewer 2013, Stewart-Ornstein 2017).
            // Rate constants scaled such that p53 rises to ~3× baseline
            // within ~30 bio-min of ATM activation, consistent with Lahav
            // 2004 single-cell HeLa dynamics. A simple 2-variable
            // (p53, MDM2) ODE of this structure settles to a stable
            // elevated p53 under sustained damage — *not* the pulsing
            // oscillation of Purvis 2012. True pulses require an mRNA
            // intermediate (Geva-Zatorsky 2006, 3-variable model); that
            // is the next PR's target. This PR delivers the foundation:
            // calibrated ATM-gated p53 stabilization + MDM2 negative
            // feedback.
            const float K_d          = 0.10f;
            const float k_p53_basal  = 1.0f;    // per sdt
            const float k_p53_deg    = 10.0f;   // per sdt, at MDM2=1, ATM=0
            const float atm_shield   = 0.70f;
            // Phase G6.2: Nutlin-3 short-circuits MDM2-mediated p53
            // degradation by occupying MDM2's p53-binding pocket
            // (Vassilev 2004). Under saturating Nutlin, effective MDM2
            // for p53-degradation is zero regardless of ATM status.
            float eff_mdm2 = c.mdm2_inhibited
                ? 0.0f
                : c.MDM2_protein * (1.0f - c.ATM_active * atm_shield);
            float p53_deg_rate =
                k_p53_deg * eff_mdm2 / (K_d + c.p53_protein);
            float dp53 = (k_p53_basal - p53_deg_rate * c.p53_protein) * sdt;

            // ── Phase G3: MDM2 mRNA explicit intermediate ──────────
            // Three-variable delay system (p53, MDM2_mRNA, MDM2_protein)
            // with Hill (n=4) transactivation. Geva-Zatorsky 2006 showed
            // this topology yields Hopf bifurcation → sustained p53
            // pulsing under constant stress; Purvis 2012 confirmed
            // ~5 h periodicity in HeLa under γ-IR.
            //
            // Hill form: empirical cooperativity n=4 is within the
            // range measured for p53-MDM2 response-element binding
            // (Riley 2008, Wu 1993); K_h = 0.12 chosen so baseline
            // p53 = 0.089 produces hill ~0.24 which combined with
            // k_mrna_max_from_p53 gives MDM2_mRNA ss ≈ 0.10.
            // Geva-Zatorsky 2006 topology: Hill activation by p53
            // with n=8 cooperativity broadens the Hopf region so
            // sustained damage produces pulses rather than a new
            // plateau. K_h=0.18 chosen so hill(0.089) is small (~0.03)
            // — baseline response is near-flat, damage response
            // snaps on above K_h.
            // Rate constants tuned for fixed-point at
            //   p53_ss=0.089, mRNA_ss=0.21, MDM2_ss=0.21
            // and slow enough damping for Hopf under damageLevel 0.3.
            const float K_h                 = 0.18f;
            const float n_hill              = 8.0f;
            const float k_mrna_basal        = 0.315f;  // per sdt
            const float k_mrna_max_from_p53 = 3.0f;    // per sdt (saturates)
            const float k_mrna_deg          = 1.5f;    // per sdt
            // p53^8 / (K_h^8 + p53^8)
            float p2 = c.p53_protein * c.p53_protein;
            float p4 = p2 * p2;
            float p8 = p4 * p4;
            float K2 = K_h * K_h;
            float K4 = K2 * K2;
            float K8 = K4 * K4;
            float hill = p8 / (K8 + p8);
            (void)n_hill; // documented-exponent constant for readability
            float dmrna = (k_mrna_basal + k_mrna_max_from_p53 * hill
                           - k_mrna_deg * c.MDM2_mRNA) * sdt;

            // MDM2 protein: translation from mRNA minus degradation.
            // t½ ≈ 30 min (Momand 1992); rates match sub-stepping unit.
            const float k_translate = 1.5f;    // per sdt
            const float k_prot_deg  = 1.5f;    // per sdt
            float dmdm2 = (k_translate * c.MDM2_mRNA
                           - k_prot_deg * c.MDM2_protein) * sdt;

            c.p53_protein  = fmaxf(0.0f, c.p53_protein  + dp53);
            c.MDM2_mRNA    = fmaxf(0.0f, c.MDM2_mRNA    + dmrna);
            c.MDM2_protein = fmaxf(0.0f, c.MDM2_protein + dmdm2);
            // Soft ceiling at 2.0 to prevent runaway under extreme damage
            // (real p53 peaks at ~10-50× baseline, Lahav 2004).
            c.p53_protein  = fminf(c.p53_protein,  2.0f);
            c.MDM2_mRNA    = fminf(c.MDM2_mRNA,    2.0f);
            c.MDM2_protein = fminf(c.MDM2_protein, 2.0f);
        }

        // CDK ODE: only advances if NOT contact-arrested
        if (!contactArrested) {
            c.cdk.step(sdt * CDK_DT_SCALE, gs);
        } else {
            // Contact arrest → p21/p27 induced, CycD suppressed (Hippo pathway)
            c.cdk.p21 = fminf(1.0f, c.cdk.p21 + 0.05f * sdt);
            c.cdk.CycD = fmaxf(0, c.cdk.CycD - 0.02f * sdt);
        }
        // ── Phase G4: p53-driven p21 transcription ────────────────
        // Replaces the previous direct damageLevel → p21 hack. p53
        // transcribes CDKN1A (p21) via a verified response element
        // (El-Deiry 1993 Cell 75:817); Hill cooperativity n=8 from
        // Ohno 2022, K_p21=0.18 (same normalised scale as MDM2 so
        // baseline p53≈0.089 gives negligible p21 input but damage-
        // driven p53≥0.2 snaps transcription on). k_p21_max tuned
        // so steady-state p21 balances the existing 0.18/s decay
        // in cdk.step at p21_ss≈0.9 under saturated hill.
        //
        // Cells that are ATM-active AND in S/G2 get a direct damage
        // cytostatic bump (~20 % of the former hack's magnitude) to
        // represent CHK1/CHK2-mediated checkpoint — a p53-independent
        // pathway — preserving observed cisplatin response kinetics.
        {
            const float K_p21    = 0.18f;
            const float k_p21_max = 0.15f;
            float p2  = c.p53_protein * c.p53_protein;
            float p4  = p2 * p2;
            float p8  = p4 * p4;
            float K2  = K_p21 * K_p21;
            float K4  = K2 * K2;
            float K8  = K4 * K4;
            float hill_p21 = p8 / (K8 + p8);
            c.cdk.p21 = fminf(1.0f, c.cdk.p21
                              + k_p21_max * hill_p21 * sdt);
            // CHK1/CHK2 S/G2 checkpoint (p53-independent, acts fast).
            if (c.damageLevel > 0.2f && (c.phase == 1 || c.phase == 2)) {
                c.cdk.p21 = fminf(1.0f,
                    c.cdk.p21 + c.damageLevel * 0.008f * sdt);
            }
        }

        c.phase=c.cdk.getPhase();
        // Hold cells with CycB > 0.25 in G2 UNTIL the G2/M gate has
        // actually passed. Once checkpointG2Passed (and dnaCheckpoint)
        // are true, let phase stay at 3 so updateMitosisProgram can fire
        // startMitosisProgram. Without this, cells oscillate between the
        // "strict" downgrade and the M branch forever.
        if (c.phase == 3 && !c.program.mitosis.active
            && !c.checkpointG2Passed) {
            c.phase = 2;
        }
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
            // hasSpace matches the new contact-inhibition rule — a cell
            // only stops dividing when truly confluent (7+ close neighbours
            // AND elevated pressure). Old threshold (<4 neighbours) was
            // blocking division at 30 % density.
            bool hasSpace = !(neighborCount >= 7 && c.localPressure > 0.4f)
                         && c.localPressure < 1.5f;
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
                // CRITICAL: mark G2M as passed so updateMitosisProgram
                // can actually call startMitosisProgram on the next tick.
                // Without this the cell sits in fake-M licensed state
                // forever — logs "M_LICENSED" every tick but mitosis
                // never runs.
                c.checkpointG2Passed = true;
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
            bool spaceAvailable = nbrs < 6 && c.localPressure < 1.8f;
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
        // Quiesce signal requires REAL confluence. Mild stress / pressure
        // from cell-cell jostling does NOT lock proliferating HeLa into
        // quiescence in vitro — that's a packed-monolayer-only response.
        bool quiesceSignal = crowded && c.localPressure > 1.5f;
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

    // ══════════════════════════════════════════════════════════════════
    //  Multi-threshold apoptosis — replaces the old 8-stage integer
    //  phase counter with a full Bcl-2/Bax/MOMP/caspase cascade
    //  (ApoptosisEngine). Eleven independent triggers map to a 0..1
    //  input amplitude each, and the engine decides commitment via
    //  the integrated arithmetic. Releases cytosol/membrane/receptor
    //  mass to MediumField in phases so the closed-system mass-balance
    //  invariant still holds.
    //
    //  All release timings / rates are in bio-seconds; `dt` is a
    //  slow-dt-scaled tick as used by the rest of updateMetabolism.
    // ══════════════════════════════════════════════════════════════════
    void updateApoptosis(SimCell& c, float dt) {
        float dt_bio = dt * SLOW_DT_SCALE * 3600.0f; // slow-dt units → bio-s
        // Convert slow-time ticks to bio-seconds. A sdt=1 step ≈
        // SLOW_DT_SCALE × 1 hour (since SLOW_DT_SCALE is hours per sdt).
        // This keeps engine rates aligned with literature (per-bio-s).

        // ── Necrotic path (legacy) ──────────────────────────────────
        // Severe hypoxic necrosis kills immediately when ATP collapses.
        // Necrotic cells dump 100% of contents as a single burst — no
        // multi-phase apoptosis; they skip straight to "bodies + leak."
        if (c.necrotic && c.ATP < 2.0f && c.alive) {
            releaseAllMass(c);
            c.alive = false;
            statDeaths++;
            c.apoptosisPhase = 4; // "past fragmentation"
            c.apoPhase = Apoptosis::CLEARED;
            return;
        }

        // ── Legacy shortcut propagation ─────────────────────────────
        // Three older code paths (stochastic fate commit, MOA_PRO_APOPTOSIS
        // drug kill, age-turnover) still write c.apoptosisPhase = 1
        // directly. Route those decisions through the engine so they
        // still drive the full morphology sequence (blebbing, bodies,
        // lysis, release) instead of leaving the cell in limbo.
        if (c.apoptosisPhase > 0 && !c.apo.state.committed
            && c.apoPhase <= Apoptosis::PRIMED) {
            c.apo.state.p53_active = 1.0f;
            c.apo.state.Bax_active = fmaxf(c.apo.state.Bax_active, 0.6f);
            c.apo.state.MOMP_pores = fmaxf(c.apo.state.MOMP_pores, 0.6f);
        }

        // ── Build the multi-threshold trigger vector ────────────────
        ApoTriggers tr;
        // ── Phase G5: route p53_active from the real p53_protein pool
        // (ATM→p53↔MDM2 axis from G2+G3+G4) rather than from raw
        // damageLevel. K=0.12 sits just above the analytical baseline
        // p53_ss = 0.089 (so healthy cells have p53_active≈0) and
        // below the sustained-damage pulse-mean (~0.15–0.20) so the
        // integrated p53 signal across pulses drives PUMA → BAX →
        // MOMP commitment. F=0.35 is the pulse-peak amplitude.
        // The previous damageLevel→p53_active mapping only fired at
        // catastrophic damage (>=0.60), missing physiologic
        // cisplatin-induced apoptosis driven by p53 pulses at
        // sub-catastrophic damage (~0.30, Zeng 2019, Cepeda 2007).
        tr.p53_active       = Apoptosis::linMap(c.p53_protein,   0.12f,                0.35f);
        tr.ROS_stress       = Apoptosis::linMap(c.ROS,           Apoptosis::ROS_K,     Apoptosis::ROS_F);
        tr.mito_dysfunction = Apoptosis::linMap(c.mitoPotential, Apoptosis::MITO_K,    Apoptosis::MITO_F);
        tr.ATP_collapse     = Apoptosis::linMap(c.atpDangerTimer,Apoptosis::ATP_K,     Apoptosis::ATP_F);
        tr.hypoxia_severe   = Apoptosis::linMap(c.hypoxiaTimer,  Apoptosis::HYPOXIA_K, Apoptosis::HYPOXIA_F);
        // Local growth-factor concentration (ng/mL). Reverse polarity:
        // low GF → high survival_loss.
        float localGF = nutrients.get(MS_GROWTH_F, c.position.x, c.position.z);
        tr.survival_loss    = Apoptosis::linMap(localGF,         Apoptosis::GF_K,      Apoptosis::GF_F);
        tr.anoikis          = Apoptosis::linMap(c.adhesionStrength, Apoptosis::ADH_K,  Apoptosis::ADH_F);
        tr.drug_pro_apop    = 0.0f; // drug subsystem removed — pending rewrite
        tr.FasL_extern      = Apoptosis::linMap(c.fasLExposure,  Apoptosis::FASL_K,    Apoptosis::FASL_F);
        tr.crowding_chronic = Apoptosis::linMap(c.chronicPressureBioSec, Apoptosis::CROWD_K,  Apoptosis::CROWD_F);
        tr.replicative      = Apoptosis::linMap(c.senescenceBioSec,     Apoptosis::REPLIC_K, Apoptosis::REPLIC_F);

        // Integrate the engine with the full trigger set.
        c.apo.step(dt_bio, tr);
        float progress = c.apo.state.apoptosis_progress;
        bool committed = c.apo.state.committed;

        // ── Map engine progress → visual ApoPhase ──────────────────
        Apoptosis::Phase prev = c.apoPhase;
        if (!committed) {
            // Still reversible. The cell is "primed" once any BH3
            // signal is meaningful; otherwise fully alive.
            float primeLevel = c.apo.state.Bax_active + c.apo.state.tBid +
                               c.apo.state.Puma + c.apo.state.Bim;
            c.apoPhase = (primeLevel > 0.05f) ? Apoptosis::PRIMED : Apoptosis::ALIVE;
        } else {
            if      (progress < 0.20f) c.apoPhase = Apoptosis::MOMP;
            else if (progress < 0.60f) c.apoPhase = Apoptosis::EXECUTION;
            else if (progress < 0.95f) c.apoPhase = Apoptosis::FRAGMENTATION;
            else                       c.apoPhase = Apoptosis::BODIES;
        }

        // Legacy compatibility: set apoptosisPhase = 1..4 mirroring
        // apoPhase so nine existing call-sites still work.
        c.apoptosisPhase = (c.apoPhase >= Apoptosis::MOMP) ? 1 :
                           (c.apoPhase == Apoptosis::PRIMED) ? 0 : 0;
        // Also mirror the legacy fate enum so the colour/score code
        // paths reflect commitment.
        if (committed && c.fate != SIM_FATE_APOPTOTIC) {
            c.fate = SIM_FATE_APOPTOTIC;
            c.fateLocked = true;
        }

        // ── Per-phase side-effects: mass release + shrinkage ─────
        float execLeakPerBioSec    = Apoptosis::CYTO_LEAK
                                   / Apoptosis::EXECUTION_DURATION_BIOSEC;
        float execLeakRecPerBioSec = Apoptosis::REC_LEAK
                                   / Apoptosis::EXECUTION_DURATION_BIOSEC;

        if (c.apoPhase == Apoptosis::EXECUTION) {
            // Slow-leak fraction of cytosol + a tiny receptor trickle.
            releaseCytosol(c, execLeakPerBioSec * dt_bio);
            releaseReceptors(c, execLeakRecPerBioSec * dt_bio);
            // Visual shrinkage (CPU-side): 35 % over the window.
            float shrinkPerBioSec = Apoptosis::SHRINK_FRAC_AT_COMPLETE
                                  / Apoptosis::EXECUTION_DURATION_BIOSEC;
            c.size = fmaxf(0.15f, c.size - shrinkPerBioSec * dt_bio * c.size);
            c.ATP  = fmaxf(0, c.ATP - 0.05f * dt_bio);
        }
        // Fragmentation transition: snapshot initial pools so main.mm
        // can partition them into apoptotic-body ledgers. main.mm reads
        // `c.bodiesSpawned == false && c.apoPhase == FRAGMENTATION`
        // exactly once and calls spawnApoptoticBodies(), which zeros
        // c.biomass/membrane/receptor (mass moves into bodies).
        if (prev != Apoptosis::FRAGMENTATION && c.apoPhase == Apoptosis::FRAGMENTATION) {
            c.initialBiomassAtDeath  = c.biomass;
            c.initialMembraneAtDeath = c.membraneMass_bm;
            c.initialReceptorAtDeath = c.receptorMass_bm;
        }
        // BODIES transition: cell leaves the live pool. In rendering
        // mode main.mm already partitioned mass into bodies (pools are
        // already zero here → releaseAllMass is a no-op). In headless
        // mode bodies never spawned, so release everything left in the
        // cell's ledger straight to the field — keeps closed-system
        // mass conservation tight on either path.
        if (prev != Apoptosis::BODIES && c.apoPhase == Apoptosis::BODIES && c.alive) {
            // Spill intracellular virions + bacteria into the medium
            // BEFORE releasing cytosolic mass — preserves the physical
            // order (pathogen particles escape a dying cell; cytosolic
            // solutes diffuse out next).
            releasePathogens(c);
            releaseAllMass(c);
            c.alive = false;
            statDeaths++;
        }
    }

    // ── Release helpers ─────────────────────────────────────────────
    // Each moves a fraction of the cell's cytosol/membrane/receptor
    // pool into the local grid cell via MediumField.exchange() using
    // the literature-calibrated partitioning table. The "fraction"
    // argument is the fraction of the *remaining* pool to release
    // this tick.
    void releaseCytosol(SimCell& c, float fraction) {
        if (fraction <= 0) return;
        fraction = fminf(1.0f, fraction);
        float dB = c.biomass * fraction;
        c.biomass -= dB;
        float flux[MS_COUNT] = {0};
        flux[MS_AA_POOL]  = dB * Apoptosis::REL_AA_PER_BIOMASS;
        flux[MS_IONS]     = dB * Apoptosis::REL_IONS_PER_BIOMASS;
        flux[MS_CALCIUM]  = dB * Apoptosis::REL_CALCIUM_PER_BIOMASS;
        flux[MS_PYRUVATE] = dB * Apoptosis::REL_PYRUVATE_PER_BIOMASS;
        flux[MS_LACTATE]  = dB * Apoptosis::REL_LACTATE_PER_BIOMASS;
        flux[MS_GLUCOSE]  = dB * Apoptosis::REL_GLUCOSE_PER_BIOMASS;
        flux[MS_WATER]    = dB * Apoptosis::REL_WATER_PER_BIOMASS;
        nutrients.exchange(c.position.x, c.position.z, flux, 1.0f);
    }
    void releaseMembrane(SimCell& c, float fraction) {
        if (fraction <= 0) return;
        fraction = fminf(1.0f, fraction);
        float dM = c.membraneMass_bm * fraction;
        c.membraneMass_bm -= dM;
        float flux[MS_COUNT] = {0};
        flux[MS_AA_POOL]  = dM * Apoptosis::REL_AA_PER_MEMBRANE;
        flux[MS_GLUCOSE]  = dM * Apoptosis::REL_GLUCOSE_PER_MEMBRANE;
        nutrients.exchange(c.position.x, c.position.z, flux, 1.0f);
    }
    void releaseReceptors(SimCell& c, float fraction) {
        if (fraction <= 0) return;
        fraction = fminf(1.0f, fraction);
        float dR = c.receptorMass_bm * fraction;
        c.receptorMass_bm -= dR;
        float flux[MS_COUNT] = {0};
        flux[MS_AA_POOL]  = dR * Apoptosis::REL_AA_PER_RECEPTOR;
        flux[MS_GLUCOSE]  = dR * Apoptosis::REL_GLUCOSE_PER_RECEPTOR;
        nutrients.exchange(c.position.x, c.position.z, flux, 1.0f);
    }
    // Dump everything remaining — necrosis path.
    void releaseAllMass(SimCell& c) {
        releaseCytosol(c, 1.0f);
        releaseMembrane(c, 1.0f);
        releaseReceptors(c, 1.0f);
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
    }
};
