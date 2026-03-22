#pragma once

#include "../core/Constants.h"
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
static float clampf(float v, float lo, float hi) { return fmaxf(lo, fminf(hi, v)); }

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

// ── Nutrient Diffusion Field (Fick 2nd law, 4-species) ──────────────────
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
        pH_field[idx] = fmaxf(0.60f, pH_field[idx] - lactateProd);
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

// ── CDK/Cyclin ODE (Novak-Tyson 7-variable model) ──────────────────────
struct CDKState {
    float CycD=0.05f, Rb=0.90f, E2F=0.02f;
    float CycE=0.01f, CycA=0.01f, CycB=0.01f, p21=0.05f;

    void randomize() {
        CycD=0.05f+randf()*0.10f; Rb=0.90f+randf()*0.08f; E2F=0.02f+randf()*0.06f;
        CycE=0.01f+randf()*0.04f; CycA=0.01f+randf()*0.03f; CycB=0.01f+randf()*0.02f;
        p21=0.05f+randf()*0.05f;
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
        // Cdc25 lower threshold + CycB faster synthesis → M phase actually reached
        float Cdc25 = fmaxf(0, CycA-0.20f)*2.0f;
        float dCycB = (0.15f*Cdc25*(1+CycB*0.8f) - CycB*(0.04f+APC_Cdc20*2.5f)) * dt_bio * (1-p21e*0.3f);
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

    bool divisionPending; float localPressure;
    float motileAngle, motileSpeed;
    int apoptosisPhase; float apoTimer;
    float adaptationTimer;

    // Drug response (PhysiPKPD model)
    float drugInternal;     // internalized drug concentration
    float drugDamage;       // accumulated drug-induced damage (0-1)
    bool  drugResistant;    // MDR resistance mutation

    void init(simd_float3 pos, int idx) {
        position=pos; velocity={0,0,0};
        radius=CELL_RADIUS_BASE+randf()*CELL_RADIUS_VAR;
        size=0.85f+randf()*0.15f; alive=true;
        ATP=75+randf()*15; stress=5; ROS=0;
        biomass=1.0f+randf()*0.1f; damageLevel=randf()*0.05f; age=randf()*30;
        mitoPotential=170+randf()*10; mitoHealth=1.0f;
        glycolytic=false; warburgTimer=0;
        hypoxiaTimer=0; hypoxiaIntensity=0; necrotic=false;
        cdk.randomize(); phase=0;
        cycleTimer=randf()*2.0f; cycleProgress=0;
        checkpointG1Passed=false; checkpointG2Passed=false;
        g1WaitTimer=0; g2WaitTimer=0;
        fate=SIM_FATE_UNDETERMINED; fateScores[0]=fateScores[1]=fateScores[2]=0;
        fateTimer=0; fateLocked=false; atpDangerTimer=0;
        glycolysisBias=0.8f+randf()*0.4f; prolifBias=0.8f+randf()*0.4f;
        rosTolerance=0.8f+randf()*0.4f; repairRate=0.8f+randf()*0.4f;
        generation=0; cloneId=idx; telomere=TELO_INIT_LENGTH; senescent=false;
        divisionPending=false; localPressure=0;
        motileAngle=randf()*M_PI*2; motileSpeed=MOTILITY_SPEED*(0.5f+randf());
        apoptosisPhase=0; apoTimer=0;
        adaptationTimer=0;
        drugInternal=0; drugDamage=0; drugResistant=false;
    }
};

// ══════════════════════════════════════════════════════════════════════════
class Simulation {
public:
    std::vector<SimCell> cells;
    NutrientField nutrients;
    float bioTime=0;
    float envO2=0.80f, envGlucose=0.70f;
    float timeScale=1.0f;
    bool paused=false;

    int statAlive=0, statProlif=0, statQuiescent=0, statApoptotic=0, statNecrotic=0;
    float statAvgATP=0;
    int statDivisions=0, statDeaths=0;
    int statPhases[4]={};
    int statGlycolytic=0;
    int nextCloneId=0;

    // Drug system
    int   activeDrugIdx = 0;       // index into DRUG_LIBRARY (0 = none)
    float drugConcentration = 0;   // applied concentration (µM)
    float statViability = 100.0f;  // % cells alive vs pre-drug count
    int   preDrugCount = 0;        // cell count when drug was applied
    float statAvgDrugDamage = 0;

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

    void init() {
        cells.clear(); nutrients.init(envO2, envGlucose);
        nextCloneId=INIT_CELLS;
        for (int i=0; i<INIT_CELLS; i++) {
            float a=(float)i/INIT_CELLS*2*M_PI;
            float r=5.0f+randf()*3.0f;
            simd_float3 pos={r*cosf(a), FLOOR_Y+2.2f, r*sinf(a)};
            SimCell c; c.init(pos, i); cells.push_back(c);
        }
    }

    void update(float dt) {
        if (paused||dt<=0) return;
        float totalDt=dt*timeScale;
        float subDt=fminf(totalDt, 0.05f);
        int steps=(int)ceilf(totalDt/subDt);
        subDt=totalDt/steps;
        bioTime += totalDt*BIO_MIN_PER_SEC*60;
        for (int s=0; s<steps; s++) {
            updatePhysics(subDt);
            for (auto& c:cells) { if(!c.alive) continue;
                updateMetabolism(c,subDt); updateDrugResponse(c,subDt);
                updateCellCycle(c,subDt); updateFate(c,subDt); updateApoptosis(c,subDt);
            }
            nutrients.diffuse(subDt, envO2, envGlucose);
        }
        processDivisions();
        cells.erase(std::remove_if(cells.begin(),cells.end(),
            [](const SimCell& c){return !c.alive;}), cells.end());
        updateStats();
    }

private:
    void updatePhysics(float dt) {
        int n=(int)cells.size();
        for (auto& c:cells) c.localPressure=0;
        for (int i=0;i<n;i++) {
            auto& a=cells[i]; if(!a.alive) continue;
            a.motileAngle += (randf()-0.5f)*0.3f*dt;
            float spd=a.motileSpeed*(a.fate==SIM_FATE_PROLIF?1.5f:0.5f);
            if (a.necrotic) spd*=0.1f;
            a.velocity.x += cosf(a.motileAngle)*spd*dt*0.5f;
            a.velocity.z += sinf(a.motileAngle)*spd*dt*0.5f;
            for (int j=i+1;j<n;j++) {
                auto& b=cells[j]; if(!b.alive) continue;
                float dx=a.position.x-b.position.x, dz=a.position.z-b.position.z;
                float dist=sqrtf(dx*dx+dz*dz);
                float minD=(a.radius+b.radius)*fminf(a.size,b.size)*0.9f;
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
            a.position.y = FLOOR_Y+a.radius*a.size*0.85f+sinf(bioTime*0.0001f+(float)i*0.7f)*0.12f;
        }
    }

    void updateMetabolism(SimCell& c, float dt) {
        float mdt=dt*MEDIUM_DT_SCALE;
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
        float oxp=localGlu*localO2*6.5f*oxpM;

        // Quiescent cells have much lower metabolic costs (G0 state)
        float basalC = 0.3f * mdt * metabolicActivity;
        float cycleC = (c.phase==1?0.15f:c.phase==3?0.25f:0.05f) * mdt * metabolicActivity;
        float motC = (c.fate==SIM_FATE_PROLIF?0.12f:0.02f) * mdt;
        float repC = c.damageLevel>0 ? c.damageLevel*0.08f*mdt : 0;
        c.ATP=clampf(c.ATP+(gly+oxp)*mdt*metabolicActivity-basalC-cycleC-motC-repC, 0, 100);

        // ── Nutrient consumption ────────────────────────────────────
        // Each cell consumes based on its metabolic activity
        // Quiescent cells consume much less → colony sustainable at confluence
        // No artificial density multiplier — real per-cell consumption creates natural depletion
        float consumeFactor = metabolicActivity * (c.fate==SIM_FATE_PROLIF ? 1.8f : 1.0f);
        nutrients.consume(c.position.x, c.position.z,
            O2_CONSUME_BASE * consumeFactor * mdt,
            GLC_CONSUME_BASE * consumeFactor * mdt, c.glycolytic);

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
        if (c.age > TURNOVER_AGE_THRESHOLD && c.apoptosisPhase == 0) {
            float deathProb = TURNOVER_PROB_PER_DT * (c.age - TURNOVER_AGE_THRESHOLD) * 0.001f;
            if (randf() < deathProb * mdt) {
                c.fate = SIM_FATE_APOPTOTIC; c.apoptosisPhase = 1; c.apoTimer = 0;
                c.fateLocked = true;
            }
        }

        c.size=clampf(0.6f+c.biomass*0.4f, 0.5f, 1.8f);
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

        if(c.cycleTimer<g1end) {
            c.cycleProgress=c.cycleTimer/g1end;
            if(!c.checkpointG1Passed && !contactArrested) {
                bool pass=c.cdk.readyForS()||(c.cdk.CycE>0.12f&&c.ATP>15&&pHOk&&c.ROS<70&&c.stress<90);
                if(pass) c.checkpointG1Passed=true;
                else { c.g1WaitTimer+=sdt;
                       // Safety valve: but NOT if contact-arrested
                       if(c.g1WaitTimer>CYCLE_G1_DUR*2.0f && !contactArrested) c.checkpointG1Passed=true;
                       c.cycleTimer=g1end*0.99f; }
            }
        } else if(c.cycleTimer<se) {
            c.cycleProgress=(c.cycleTimer-g1end)/CYCLE_S_DUR;
        } else if(c.cycleTimer<g2end) {
            c.cycleProgress=(c.cycleTimer-se)/CYCLE_G2_DUR;
            if(!c.checkpointG2Passed) {
                bool pass=c.cdk.readyForM()||(c.cdk.CycB>0.10f&&c.ATP>12&&pHOk&&c.stress<85);
                if(pass) c.checkpointG2Passed=true;
                else { c.g2WaitTimer+=sdt;
                       if(c.g2WaitTimer>CYCLE_G2_DUR*3.5f) c.checkpointG2Passed=true;
                       c.cycleTimer=g2end*0.99f; }
            }
        } else {
            c.cycleProgress = 1.0f;
            c.phase = 3; // M phase

            // ── SPACE CHECK before division ─────────────────────────
            // Cell must have physical room for a daughter cell
            // If surrounded, skip division → reset cycle (like a failed mitosis)
            bool hasSpace = neighborCount < 4 && c.localPressure < 1.5f;
            bool canDiv = (c.fate==SIM_FATE_PROLIF || c.fate==SIM_FATE_UNDETERMINED);
            if(canDiv && !c.senescent && hasSpace && (int)cells.size()<MAX_CELLS) {
                c.divisionPending = true;
                if(c.fate==SIM_FATE_UNDETERMINED) {
                    c.fate = SIM_FATE_PROLIF; c.fateLocked = true;
                }
            }

            // Reset for next cycle
            c.cycleTimer = 0; c.cycleProgress = 0;
            c.checkpointG1Passed = false; c.checkpointG2Passed = false;
            c.g1WaitTimer = 0; c.g2WaitTimer = 0;
            c.cdk.resetForNewCycle(gs);
        }
    }

    void updateFate(SimCell& c, float dt) {
        if(c.apoptosisPhase>0||c.senescent) return;

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
        std::vector<SimCell> daughters;
        for(auto& m:cells) {
            if(!m.divisionPending||!m.alive) continue;
            m.divisionPending=false;
            float ang=randf()*M_PI*2, off=m.radius*1.85f*m.size;
            simd_float3 dp={m.position.x+cosf(ang)*off, m.position.y, m.position.z+sinf(ang)*off};
            if(sqrtf(dp.x*dp.x+dp.z*dp.z)>SCENE_BOUND-3) continue;
            SimCell d; d.init(dp, nextCloneId++);
            float split=0.48f+randf()*0.04f;
            d.ATP=fmaxf(35,m.ATP*split*1.3f); d.stress=fminf(15,m.stress*0.25f);
            d.ROS=m.ROS*split; d.damageLevel=m.damageLevel*split;
            d.biomass=1.0f; d.mitoHealth=m.mitoHealth*(0.45f+randf()*0.1f);
            d.glycolysisBias=clampf(m.glycolysisBias+(randf()-0.5f)*GENOME_MUTATION_RATE*2, 0.3f, 2.2f);
            d.prolifBias=clampf(m.prolifBias+(randf()-0.5f)*GENOME_MUTATION_RATE*2, 0.3f, 2.2f);
            d.rosTolerance=clampf(m.rosTolerance+(randf()-0.5f)*GENOME_MUTATION_RATE*2, 0.3f, 2.2f);
            d.repairRate=clampf(m.repairRate+(randf()-0.5f)*GENOME_MUTATION_RATE*2, 0.2f, 2.0f);
            d.generation=m.generation+1;
            d.cloneId=(fabsf(d.glycolysisBias-m.glycolysisBias)>0.08f)?nextCloneId++:m.cloneId;
            d.telomere=m.telomere-TELO_LOSS_PER_DIV*(0.8f+randf()*0.4f);
            if(d.telomere<TELO_CRITICAL) d.senescent=true;
            m.biomass=fmaxf(0.6f,m.biomass*(1-split)); m.ATP=fmaxf(30,m.ATP*(1-split));
            m.telomere-=TELO_LOSS_PER_DIV*(0.8f+randf()*0.4f);
            if(m.telomere<TELO_CRITICAL) m.senescent=true;
            d.velocity.x=cosf(ang)*1.2f; d.velocity.z=sinf(ang)*1.2f;
            m.velocity.x-=cosf(ang)*0.5f; m.velocity.z-=sinf(ang)*0.5f;
            daughters.push_back(d); statDivisions++;
        }
        for(auto& d:daughters) cells.push_back(d);
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
