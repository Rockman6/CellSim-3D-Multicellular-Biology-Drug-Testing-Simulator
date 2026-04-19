#import <Cocoa/Cocoa.h>
#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#import <UniformTypeIdentifiers/UniformTypeIdentifiers.h>

#define GLFW_INCLUDE_NONE
#define GLFW_EXPOSE_NATIVE_COCOA
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_metal.h"
#include "implot.h"

#include "gpu/MetalContext.h"
#include "gpu/ShaderTypes.h"
#include "render/Camera.h"
#include "render/MeshLibrary.h"
#include "render/GLBLoader.h"
#include "simulation/Simulation.h"
#include "simulation/TelemetryLog.h"

static TelemetryLog gTelemetryLog;
static bool gTelemetryOpened = false;
#include "simulation/CentralDogma.h"
#include "simulation/IntracellularPhysics.h"
#include "simulation/biochem/CellBiologyEngine.h"
#include "molgen/MoleculeCache.h"
#include "molgen/MolGenRunner.h"
#include "molgen/ProteinMeshCache.h"

#include <cstdio>
#include <sys/stat.h>
#include <errno.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <map>
#include <set>

// ══════════════════════════════════════════════════════════════════════════
//  CellSim — Full Simulation + ImGui Research UI
// ══════════════════════════════════════════════════════════════════════════

static const int INIT_WIDTH  = 1600;
static const int INIT_HEIGHT = 1000;

static Camera       gCamera;
static MetalContext  gCtx;
static MeshLibrary  gMeshes;
static Simulation   gSim;
static std::string  gModelDir;
static SimMode      gSimMode = MODE_SINGLE_CELL;

// ── Startup Setup overlay ──────────────────────────────────────────────
// Shown immediately after app launch. User picks a preset template or
// manually tunes the initial conditions (cell count, medium composition,
// timescale target), then clicks Start to begin the sim. The old
// "Mode [Single Cell / Colony]" switcher is gone — there's just ONE
// configurable cell mode now. Colony-specific internals are kept in the
// Simulation class so headless validation still works with the same code.
struct SimSetup {
    bool   showOverlay      = true;   // visible until user clicks Start
    int    initCells        = 5;      // seed count (1 = classic single-cell)
    float  glucoseMM        = 25.0f;
    float  glutamineMM      = 6.0f;
    float  aaPoolMM         = 7.0f;
    float  o2MM             = 0.20f;
    float  co2MM            = 2.40f;
    float  pH               = 7.40f;
    float  growthFactorNgML = 50.0f;
    float  initialTimeScale = 1.0f;   // 1× = 3 bio-min / wall-sec
    int    templateIdx      = 0;      // which preset is highlighted
};
static SimSetup gSetup;

// Molecule rendering
static MoleculeCache gMolCache;
static MolGenRunner  gMolGen;
static std::string   gMoleculeDir;
static std::string   gProjectRoot;

// GPU buffers for molecule rendering
struct MolRenderData {
    id<MTLBuffer> atomBuffer = nil;
    id<MTLBuffer> bondBuffer = nil;
    int atomCount = 0;
    int bondCount = 0;
    std::string currentMolId;

    void clear() {
        atomBuffer = nil; bondBuffer = nil;
        atomCount = 0; bondCount = 0;
        currentMolId.clear();
    }
};
static MolRenderData gMolRender;
static int gSelectedMolecule = 0;

// Protein structure cache
static ProteinMeshCache gProtCache;
static std::string gProteinDir;

// Central dogma + intracellular physics.
// gCDogma is a "view" that mirrors gSim.cells[focusIdx].program.cdogma so
// the existing teaching/animation code can keep reading a single global.
// The canonical per-cell replication state lives on each SimCell.
// gBoundDogmaCellUid tracks which cell gCDogma is currently mirroring so
// bindDogmaStateToPrimaryCell can write the mirror back when focus changes.
static CentralDogmaState gCDogma;
static MitosisState gMitosis;
static int gBoundDogmaCellUid = -1;
static int gBoundMitosisCellUid = -1;

// Per-cell intracellular physics — every cell gets its own particle system
struct CellInterior {
    IntracellularPhysics phys;
    bool initialized = false;
    float postMitoticCondense = 0.0f;
};
static std::map<int, CellInterior> gCellInteriors;

struct OrganelleMotionState {
    bool initialized = false;
    float lastUpdateTime = -1.0f;
    simd_float3 golgiPos = {0, 0, -0.3f};
    simd_float3 golgiVel = {0, 0, 0};
    simd_float3 mitoPos[3];
    simd_float3 mitoVel[3];
};
static std::map<int, OrganelleMotionState> gOrganelleMotionStates;

// Primary cell (cell[0]) uses gIntraPhys for backwards compat
static IntracellularPhysics gIntraPhys;
static bool gIntraInitialized = false;

// Master cell biology engine (16 biochemistry modules)
static CellBiologyEngine gCellBio;

// Daughter cell physics state (post-division)
static simd_float3 sDaughterPosA = {0,0,0}, sDaughterPosB = {0,0,0};
static simd_float3 sDaughterVelA = {0,0,0}, sDaughterVelB = {0,0,0};
static float sDaughterR = 0;
static bool sDaughtersInitialized = false;
static int sMitosisDNAStartCount = 0;
static int sMitosisExpectedDNACount = 0;
static std::string gRunSessionTag;
static std::string gRunLaunchTimeText;
static std::string gMitosisLogPath;
static std::string gDiagLogPath;
static std::string gCentrosomeLogPath;

// ── Mitosis log ─────────────────────────────────────────────────────────
static FILE* sMitosisLog = nullptr;
static int sMitosisLogPhase = -1;
static std::string makeRunSessionTag() {
    time_t now = time(nullptr);
    struct tm tmNow = {};
    localtime_r(&now, &tmNow);
    char buf[64];
    strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", &tmNow);
    return std::string(buf);
}

static std::string makeRunLaunchText() {
    time_t now = time(nullptr);
    struct tm tmNow = {};
    localtime_r(&now, &tmNow);
    char buf[96];
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tmNow);
    return std::string(buf);
}

static void initializeSessionLogging(const std::string& projectRoot) {
    if (gRunSessionTag.empty()) gRunSessionTag = makeRunSessionTag();
    if (gRunLaunchTimeText.empty()) gRunLaunchTimeText = makeRunLaunchText();

    std::string logsDir = projectRoot.empty() ? "/Users/henry/CellSim/logs" : (projectRoot + "/logs");
    gMitosisLogPath = logsDir + "/mitosis_" + gRunSessionTag + ".log";
    gDiagLogPath = logsDir + "/diag_" + gRunSessionTag + ".log";
    gCentrosomeLogPath = logsDir + "/centrosome_" + gRunSessionTag + ".log";

    std::string latestPath = logsDir + "/latest_session.txt";
    FILE* f = fopen(latestPath.c_str(), "w");
    if (f) {
        fprintf(f, "session=%s\n", gRunSessionTag.c_str());
        fprintf(f, "launch_time=%s\n", gRunLaunchTimeText.c_str());
        fprintf(f, "mitosis_log=%s\n", gMitosisLogPath.c_str());
        fprintf(f, "diag_log=%s\n", gDiagLogPath.c_str());
        fprintf(f, "centrosome_log=%s\n", gCentrosomeLogPath.c_str());
        fclose(f);
    }
}

static void mitosisLog(const char* fmt, ...) {
    if (!sMitosisLog) {
        const char* path = gMitosisLogPath.empty()
            ? "/Users/henry/CellSim/logs/mitosis_log.txt"
            : gMitosisLogPath.c_str();
        sMitosisLog = fopen(path, "w");
        if (!sMitosisLog) return;
        fprintf(sMitosisLog, "session=%s launch=%s\n",
            gRunSessionTag.empty() ? "unknown" : gRunSessionTag.c_str(),
            gRunLaunchTimeText.empty() ? "unknown" : gRunLaunchTimeText.c_str());
    }
    va_list ap; va_start(ap, fmt);
    vfprintf(sMitosisLog, fmt, ap);
    va_end(ap);
    fflush(sMitosisLog);
}

static int particleMitosisHalf(const IntraParticle& p) {
    int half = (p.mitosisHalf >= 0) ? p.mitosisHalf : ((p.position.x >= 0) ? 0 : 1);
    if (half < 0 || half > 1) half = 0;
    return half;
}

static bool isNuclearParticleType(ParticleType type) {
    return type == PT_DNA_NODE || type == PT_DNA_POLYMERASE || type == PT_RNA_POLYMERASE ||
           type == PT_SPLICEOSOME || type == PT_NUCLEAR_PORE;
}

static bool canCrossCleavageFurrow(const IntraParticle& p) {
    return p.type == PT_MOLECULE && (p.tag == TAG_WATER || p.tag == TAG_LIPID);
}

static bool isCleavageLockedParticle(const IntraParticle& p) {
    if (!p.active) return false;
    if (isNuclearParticleType(p.type)) return false;
    return !canCrossCleavageFurrow(p);
}

static bool mitosisCleavageCompartmentsActive() {
    return gMitosis.active &&
           gMitosis.phase >= MITO_METAPHASE &&
           gMitosis.phase <= MITO_CYTOKINESIS;
}

static float mitosisCytoplasmPartitionProgress() {
    switch (gMitosis.phase) {
        case MITO_METAPHASE:
            return 0.14f;
        case MITO_ANAPHASE:
            return 0.34f + gMitosis.chromatidSeparation * 0.22f;
        case MITO_TELOPHASE:
            return 0.62f + gMitosis.nucleusAlpha * 0.18f;
        case MITO_CYTOKINESIS:
            return 0.86f + gMitosis.furrowDepth * 0.14f;
        default:
            return 0.0f;
    }
}

static simd_float3 clampPointToSphere(simd_float3 point, simd_float3 center, float radius) {
    float dx = point.x - center.x;
    float dy = point.y - center.y;
    float dz = point.z - center.z;
    float dist = sqrtf(dx * dx + dy * dy + dz * dz);
    if (dist <= radius || dist < 0.0001f) return point;
    float invDist = 1.0f / dist;
    return {
        center.x + dx * invDist * radius,
        center.y + dy * invDist * radius,
        center.z + dz * invDist * radius
    };
}

static simd_float3 mitosisCytoplasmCompartmentCenter(int half, float cellR) {
    float sign = (half == 0) ? 1.0f : -1.0f;
    float progress = mitosisCytoplasmPartitionProgress();
    float targetX = cellR * (0.08f + 0.28f * progress);
    return {sign * targetX, 0, 0};
}

static float mitosisCytoplasmCompartmentRadius(float cellR) {
    float progress = mitosisCytoplasmPartitionProgress();
    float radius = cellR * (0.96f - 0.20f * progress);
    return fmaxf(radius, cellR * 0.58f);
}

static simd_float3 partitionPointIntoDaughterHalf(simd_float3 point, int half, float cellR,
                                                  float keepLocal = 0.72f) {
    simd_float3 center = mitosisCytoplasmCompartmentCenter(half, cellR);
    simd_float3 local = {
        point.x * keepLocal * 0.42f,
        point.y * keepLocal,
        point.z * keepLocal
    };
    float sign = (half == 0) ? 1.0f : -1.0f;
    local.x += sign * fabsf(point.x) * 0.12f;
    simd_float3 shifted = {
        center.x + local.x,
        center.y + local.y,
        center.z + local.z
    };
    return clampPointToSphere(shifted, center, mitosisCytoplasmCompartmentRadius(cellR) * 0.82f);
}

static void enforceCleavageCompartmentBarrier(IntracellularPhysics& phys, float cellR) {
    if (!mitosisCleavageCompartmentsActive()) return;

    float progress = mitosisCytoplasmPartitionProgress();
    float minAbsX = cellR * (0.02f + 0.22f * progress);
    for (auto& p : phys.particles) {
        if (!isCleavageLockedParticle(p)) continue;
        int half = particleMitosisHalf(p);
        float sign = (half == 0) ? 1.0f : -1.0f;
        simd_float3 center = mitosisCytoplasmCompartmentCenter(half, cellR);
        float radius = p.confineRadius;
        if (radius <= 0.0f) radius = mitosisCytoplasmCompartmentRadius(cellR);

        intraProjectInsideConfinement(p);
        if (sign * p.position.x < minAbsX) {
            p.position.x = sign * minAbsX;
        }
        if (sign * p.velocity.x < 0.0f) {
            p.velocity.x *= 0.15f;
        }
        p.position = clampPointToSphere(p.position, center, fmaxf(radius - p.radius * 0.5f, 0.0f));
    }
}

static float primaryPostMitoticCondenseBlend() {
    if (gSim.cells.empty()) return 0.0f;
    return fmaxf(0.0f, fminf(1.0f, gSim.cells[0].postDivisionRecovery / 6.0f));
}

// Central-dogma state now lives on each SimCell (cell.program.cdogma).
// gCDogma is a mirror of the currently-focused cell's program.cdogma so the
// existing teaching/animation code in main.mm keeps working without touching
// every gCDogma.xxx call site. Binding saves gCDogma back into the previously
// focused cell and loads the new primary cell's state into gCDogma.

static SimCell* findCellByUid(int cellUid) {
    if (cellUid < 0) return nullptr;
    for (auto& c : gSim.cells) {
        if (c.cellUid == cellUid) return &c;
    }
    return nullptr;
}

static void persistBoundDogmaState() {
    if (gBoundDogmaCellUid < 0) return;
    SimCell* bound = findCellByUid(gBoundDogmaCellUid);
    if (bound) bound->program.cdogma = gCDogma;
}

static void persistBoundMitosisState() {
    if (gBoundMitosisCellUid < 0) return;
    SimCell* bound = findCellByUid(gBoundMitosisCellUid);
    if (bound) bound->program.mitosis = gMitosis;
}

static void bindDogmaStateToPrimaryCell(bool forceReload = false) {
    if (gSim.cells.empty()) {
        persistBoundDogmaState();
        gBoundDogmaCellUid = -1;
        gCDogma.init();
        return;
    }

    int primaryUid = gSim.cells[0].cellUid;
    gSim.cells[0].program.ensureCDogmaInitialized();
    if (!forceReload && gBoundDogmaCellUid == primaryUid) return;

    persistBoundDogmaState();
    gCDogma = gSim.cells[0].program.cdogma;
    gBoundDogmaCellUid = primaryUid;
}

static void bindMitosisStateToPrimaryCell(bool forceReload = false) {
    if (gSim.cells.empty()) {
        persistBoundMitosisState();
        gBoundMitosisCellUid = -1;
        gMitosis = MitosisState{};
        return;
    }

    int primaryUid = gSim.cells[0].cellUid;
    if (!forceReload && gBoundMitosisCellUid == primaryUid) return;

    persistBoundMitosisState();
    gMitosis = gSim.cells[0].program.mitosis;
    gBoundMitosisCellUid = primaryUid;
}

static void initializeDogmaStatesFromSimulation() {
    gBoundDogmaCellUid = -1;
    gBoundMitosisCellUid = -1;
    gCDogma.init();
    gMitosis = MitosisState{};
    for (auto& c : gSim.cells) {
        if (!c.alive) continue;
        c.program.ensureCDogmaInitialized();
    }
    if (!gSim.cells.empty()) {
        gSim.cells[0].program.cdogma = gCDogma;
        gBoundDogmaCellUid = gSim.cells[0].cellUid;
        gMitosis = gSim.cells[0].program.mitosis;
        gBoundMitosisCellUid = gSim.cells[0].cellUid;
    }
}

// Advance every non-primary cell's replication/transcription state.
// The primary cell's state is stepped via gCDogma in the main simulation
// loop (see gCDogma.update(...) downstream) and then mirrored back to
// cells[0].program.cdogma via persistBoundDogmaState().
static void stepBackgroundDogmaStates(float dt) {
    float dogmaDt = fminf(fmaxf(dt, 1.0f / 240.0f), 1.0f / 30.0f);
    int primaryUid = gSim.cells.empty() ? -1 : gSim.cells[0].cellUid;
    for (auto& c : gSim.cells) {
        if (!c.alive) continue;
        c.program.ensureCDogmaInitialized();
        if (c.cellUid == primaryUid) continue;
        c.program.cdogma.update(dogmaDt, c.phase);
    }
}

static bool interiorGenomeAlreadyReplicated(const IntracellularPhysics& phys) {
    int dnaOriginal = 0, dnaSister = 0;
    for (const auto& p : phys.particles) {
        if (!p.active || p.type != PT_DNA_NODE) continue;
        if (p.mitosisHalf == 1) dnaSister++;
        else dnaOriginal++;
    }
    return dnaOriginal > 0 && dnaSister >= dnaOriginal;
}

static float mitosisCenterWeight(const IntraParticle& p) {
    switch (p.type) {
        case PT_DNA_NODE:
            return 6.0f;
        case PT_RNA_POLYMERASE:
        case PT_SPLICEOSOME:
        case PT_NUCLEAR_PORE:
            return 3.0f;
        case PT_RIBOSOME_SMALL:
        case PT_RIBOSOME_LARGE:
        case PT_CHAPERONE:
        case PT_VESICLE_COPII:
        case PT_VESICLE_SECRETORY:
            return 1.2f;
        case PT_MOLECULE:
            return 0.03f;
        default:
            return 1.8f;
    }
}

struct DNADistribution {
    int total = 0;
    int half[2] = {0, 0};
};

static DNADistribution countActiveDNADistribution() {
    DNADistribution dna = {};
    for (const auto& p : gIntraPhys.particles) {
        if (!p.active || p.type != PT_DNA_NODE) continue;
        dna.total++;
        dna.half[particleMitosisHalf(p)]++;
    }
    return dna;
}

static void countActiveDNAVisuals(const IntracellularPhysics& phys, int& totalDNA, int& sisterDNA) {
    totalDNA = 0;
    sisterDNA = 0;
    for (const auto& p : phys.particles) {
        if (!p.active || p.type != PT_DNA_NODE) continue;
        totalDNA++;
        if (p.mitosisHalf == 1) sisterDNA++;
    }
}

static bool hasPreparedSPhaseReplicatedDNAVisuals(const IntracellularPhysics& phys) {
    for (const auto& p : phys.particles) {
        if (p.type == PT_DNA_NODE && p.mitosisHalf == 1) return true;
    }
    return false;
}

static bool hasDuplicatedNuclearAuxForDivision(const IntracellularPhysics& phys) {
    for (const auto& p : phys.particles) {
        if (!p.active || p.mitosisHalf != 1) continue;
        if (isNuclearParticleType(p.type) && p.type != PT_DNA_NODE) return true;
    }
    return false;
}

static bool interiorPreparedForPhysicalDivision(const IntracellularPhysics& phys) {
    return interiorGenomeAlreadyReplicated(phys) && hasDuplicatedNuclearAuxForDivision(phys);
}

static void prepareGenomeForSPhaseReplicationVisuals(IntracellularPhysics& phys, float cellR) {
    if (hasPreparedSPhaseReplicatedDNAVisuals(phys)) return;

    int originalCount = (int)phys.particles.size();
    for (int i = 0; i < originalCount; i++) {
        auto& p = phys.particles[i];
        if (!p.active || p.type != PT_DNA_NODE) continue;
        p.mitosisHalf = 0;
    }

    std::vector<IntraParticle> source = phys.particles;
    std::vector<int> duplicateIndexMap(originalCount, -1);
    int duplicateCount = 0;
    for (int i = 0; i < originalCount; i++) {
        const auto& src = source[i];
        if (!src.active || src.type != PT_DNA_NODE) continue;
        duplicateCount++;
    }
    phys.particles.reserve(originalCount + duplicateCount);

    for (int i = 0; i < originalCount; i++) {
        const auto& src = source[i];
        if (!src.active || src.type != PT_DNA_NODE) continue;

        IntraParticle dup = src;
        dup.active = false;
        dup.mitosisHalf = 1;
        dup.velocity = {0, 0, 0};
        dup.position.x -= cellR * 0.006f;
        dup.home.x -= cellR * 0.006f;
        dup.position.z += cellR * 0.028f;
        dup.home.z += cellR * 0.028f;
        dup.spawnPos = dup.home;
        dup.glowIntensity = 1.0f;
        duplicateIndexMap[i] = (int)phys.particles.size();
        phys.particles.push_back(dup);
    }

    for (int i = 0; i < originalCount; i++) {
        int dupIdx = duplicateIndexMap[i];
        if (dupIdx < 0) continue;
        auto& dup = phys.particles[dupIdx];
        int oldTether = source[i].tetherId;
        dup.tetherId =
            (oldTether >= 0 && oldTether < originalCount) ? duplicateIndexMap[oldTether] : -1;
        int oldJobTarget = source[i].jobTargetId;
        dup.jobTargetId =
            (oldJobTarget >= 0 && oldJobTarget < originalCount) ? duplicateIndexMap[oldJobTarget] : -1;
        if (dup.job == JOB_ATTACH && dup.jobTargetId < 0) dup.job = JOB_IDLE;
    }
}

static void syncGenomeReplicationVisuals(IntracellularPhysics& phys,
                                         const CentralDogmaState& dogma,
                                         float cellR,
                                         bool freezeDuringMitosis = false) {
    if (!dogma.sPhaseProgramStarted || freezeDuringMitosis) return;

    prepareGenomeForSPhaseReplicationVisuals(phys, cellR);
    for (auto& p : phys.particles) {
        if (p.type != PT_DNA_NODE || p.mitosisHalf != 1) continue;
        bool replicated = dogma.isBaseReplicated(p.stateIndex);
        p.active = replicated;
        if (replicated) {
            p.glowIntensity = 0.95f + 0.35f * dogma.replicationProgress;
        }
    }
}

struct MitosisHalfCenters {
    simd_float3 local[2];
    bool valid[2];
};

static MitosisHalfCenters computeHalfCentersForInterior(const IntracellularPhysics& phys,
                                                        float fallbackX,
                                                        bool nuclearOnly = false) {
    MitosisHalfCenters centers = {};
    centers.local[0] = { 0, 0, 0 };
    centers.local[1] = { 0, 0, 0 };
    centers.valid[0] = false;
    centers.valid[1] = false;

    float weight[2] = {0, 0};
    for (const auto& p : phys.particles) {
        if (!p.active) continue;
        if (nuclearOnly && !isNuclearParticleType(p.type)) continue;
        int half = particleMitosisHalf(p);
        float w = mitosisCenterWeight(p);
        centers.local[half].x += p.position.x * w;
        centers.local[half].y += p.position.y * w;
        centers.local[half].z += p.position.z * w;
        weight[half] += w;
    }

    for (int half = 0; half < 2; half++) {
        if (weight[half] > 0.0f) {
            float invW = 1.0f / weight[half];
            centers.local[half].x *= invW;
            centers.local[half].y *= invW;
            centers.local[half].z *= invW;
            centers.valid[half] = true;
        } else {
            centers.local[half] = {(half == 0) ? fallbackX : -fallbackX, 0, 0};
        }
    }

    return centers;
}

static MitosisHalfCenters computeMitosisHalfCenters(float fallbackX, bool nuclearOnly = false) {
    return computeHalfCentersForInterior(gIntraPhys, fallbackX, nuclearOnly);
}

static float mitosisNuclearTargetX(float cellR) {
    float phaseDrive = 0.0f;
    switch (gMitosis.phase) {
        case MITO_PROMETAPHASE: phaseDrive = 0.10f + gMitosis.spindleAssembly * 0.06f; break;
        case MITO_METAPHASE:    phaseDrive = 0.16f; break;
        case MITO_ANAPHASE:     phaseDrive = 0.34f + gMitosis.chromatidSeparation * 0.30f; break;
        case MITO_TELOPHASE:    phaseDrive = 0.78f; break;
        case MITO_CYTOKINESIS:  phaseDrive = 1.00f; break;
        default: break;
    }
    float poleAbsX = fmaxf(fabsf(gMitosis.poleA.x), fabsf(gMitosis.poleB.x));
    float daughterDriven = cellR * (0.18f + 0.40f * phaseDrive);
    float spindleDriven = poleAbsX + cellR * (0.05f + 0.08f * phaseDrive);
    float targetX = fmaxf(daughterDriven, spindleDriven);
    return fminf(targetX, cellR * 0.68f);
}

static float mitosisDaughterCenterBlend() {
    switch (gMitosis.phase) {
        case MITO_PROMETAPHASE: return 0.0f;
        case MITO_METAPHASE:    return 0.0f;
        case MITO_ANAPHASE:     return 0.16f + gMitosis.chromatidSeparation * 0.28f;
        case MITO_TELOPHASE:    return 0.58f + gMitosis.nucleusAlpha * 0.24f;
        case MITO_CYTOKINESIS:  return 0.86f + gMitosis.furrowDepth * 0.14f;
        default:                return 0.0f;
    }
}

// Persistent chromatin-condensation blend. This variable decays smoothly
// from whatever level mitosis set it to, toward 0 (fully decondensed),
// regardless of whether gMitosis.active is true. Rendering uses this value
// for DNA-particle tinting and chromosome-rod opacity so chromosomes fade
// continuously into the interphase nucleus instead of vanishing the frame
// cytokinesis ends.
static float gChromatinDecondenseState = 0.0f;

static float mitosisCondenseBlend() {
    float live = gMitosis.active
        ? fmaxf(0.0f, fminf(1.0f, gMitosis.chromatinCondensation))
        : 0.0f;
    // Track the maximum of live condensation and the persistent state,
    // then decay the persistent state toward live. This lets the value
    // rise instantly when mitosis starts but decay over several frames
    // after cytokinesis.
    if (live > gChromatinDecondenseState) gChromatinDecondenseState = live;
    else gChromatinDecondenseState += (live - gChromatinDecondenseState) * 0.02f;
    return gChromatinDecondenseState;
}

static int mitosisChromosomeSlot(const IntraParticle& p) {
    if (p.stateIndex <= 0) return 0;
    int slot = (p.stateIndex * MitosisState::NUM_CHROMO) / fmax(HBB_LENGTH, 1);
    if (slot < 0) slot = 0;
    if (slot >= MitosisState::NUM_CHROMO) slot = MitosisState::NUM_CHROMO - 1;
    return slot;
}

static simd_float3 mitosisDaughterCenter(int half, float cellR) {
    float sign = (half == 0) ? 1.0f : -1.0f;
    float targetX = mitosisNuclearTargetX(cellR);
    float blend = mitosisDaughterCenterBlend();
    return {sign * targetX * blend, 0, 0};
}

static simd_float3 mitosisChromosomeAnchor(const IntraParticle& p, float cellR) {
    int half = particleMitosisHalf(p);
    int slot = mitosisChromosomeSlot(p);
    const auto& chr = gMitosis.chromosomes[slot];
    simd_float3 chromo = (half == 0) ? chr.position : chr.sisterPosition;
    simd_float3 daughter = mitosisDaughterCenter(half, cellR);

    float daughterBlend = mitosisDaughterCenterBlend();
    simd_float3 anchor = {
        chromo.x * (1.0f - daughterBlend) + daughter.x * daughterBlend,
        chromo.y * (1.0f - daughterBlend) + daughter.y * daughterBlend,
        chromo.z * (1.0f - daughterBlend) + daughter.z * daughterBlend
    };

    float slotAngle = ((float)slot / (float)MitosisState::NUM_CHROMO) * 2.0f * M_PI;
    anchor.y += sinf(slotAngle) * cellR * 0.035f;
    anchor.z += cosf(slotAngle) * cellR * 0.045f;
    return anchor;
}

static simd_float3 mitosisDNAHomeTarget(const IntraParticle& p, float cellR) {
    simd_float3 anchor = mitosisChromosomeAnchor(p, cellR);
    int slot = mitosisChromosomeSlot(p);
    const auto& chr = gMitosis.chromosomes[slot];
    auto normalizeLocal = [](simd_float3 v, simd_float3 fallback) -> simd_float3 {
        float d = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
        if (d < 1.0e-5f) return fallback;
        float invD = 1.0f / d;
        return {v.x * invD, v.y * invD, v.z * invD};
    };
    auto crossLocal = [](simd_float3 a, simd_float3 b) -> simd_float3 {
        return {
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        };
    };
    simd_float3 fallbackY = {0, 1, 0};
    simd_float3 fallbackX = {1, 0, 0};
    simd_float3 fallbackZ = {0, 0, 1};
    simd_float3 axis = normalizeLocal(chr.axis, fallbackY);
    simd_float3 radial = normalizeLocal(crossLocal(axis, fallbackZ), fallbackX);
    simd_float3 binormal = normalizeLocal(crossLocal(axis, radial), fallbackZ);

    float slotSpan = fmaxf((float)HBB_LENGTH / (float)MitosisState::NUM_CHROMO, 1.0f);
    float slotOffset = fmodf((float)fmax(p.stateIndex, 0), slotSpan);
    float slotLocal = (slotOffset + 0.5f) / slotSpan - 0.5f;
    float condense = mitosisCondenseBlend();
    float armExtent = cellR * (0.22f - 0.11f * condense);
    float coilRadius = cellR * (0.065f * (1.0f - condense) + 0.014f);
    float coilHeight = cellR * (0.048f * (1.0f - condense) + 0.012f);
    float centromerePinch = 0.30f + 0.70f * fminf(fabsf(slotLocal) * 2.0f, 1.0f);
    float coil = slotLocal * 10.0f * M_PI + (float)(p.stateIndex % 17) * 0.31f;

    simd_float3 offset = {
        axis.x * (slotLocal * armExtent) + radial.x * cosf(coil) * coilRadius * centromerePinch + binormal.x * sinf(coil) * coilHeight * centromerePinch,
        axis.y * (slotLocal * armExtent) + radial.y * cosf(coil) * coilRadius * centromerePinch + binormal.y * sinf(coil) * coilHeight * centromerePinch,
        axis.z * (slotLocal * armExtent) + radial.z * cosf(coil) * coilRadius * centromerePinch + binormal.z * sinf(coil) * coilHeight * centromerePinch
    };
    return {
        anchor.x + offset.x,
        anchor.y + offset.y,
        anchor.z + offset.z
    };
}

static simd_float3 mitosisAuxNuclearHomeTarget(const IntraParticle& p, float cellR) {
    simd_float3 anchor = mitosisChromosomeAnchor(p, cellR);
    int slot = mitosisChromosomeSlot(p);
    const auto& chr = gMitosis.chromosomes[slot];
    auto normalizeLocal = [](simd_float3 v, simd_float3 fallback) -> simd_float3 {
        float d = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
        if (d < 1.0e-5f) return fallback;
        float invD = 1.0f / d;
        return {v.x * invD, v.y * invD, v.z * invD};
    };
    auto crossLocal = [](simd_float3 a, simd_float3 b) -> simd_float3 {
        return {
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        };
    };
    simd_float3 fallbackY = {0, 1, 0};
    simd_float3 fallbackX = {1, 0, 0};
    simd_float3 fallbackZ = {0, 0, 1};
    simd_float3 axis = normalizeLocal(chr.axis, fallbackY);
    simd_float3 radial = normalizeLocal(crossLocal(axis, fallbackY), fallbackX);
    simd_float3 binormal = normalizeLocal(crossLocal(axis, radial), fallbackZ);
    float relax = 1.0f - mitosisCondenseBlend() * 0.70f;
    simd_float3 local = {
        radial.x * p.spawnPos.x * relax * 0.22f + binormal.x * p.spawnPos.z * relax * 0.22f + axis.x * p.spawnPos.y * relax * 0.10f,
        radial.y * p.spawnPos.x * relax * 0.22f + binormal.y * p.spawnPos.z * relax * 0.22f + axis.y * p.spawnPos.y * relax * 0.10f,
        radial.z * p.spawnPos.x * relax * 0.22f + binormal.z * p.spawnPos.z * relax * 0.22f + axis.z * p.spawnPos.y * relax * 0.10f
    };
    return {anchor.x + local.x, anchor.y + local.y, anchor.z + local.z};
}

static void restoreInterphaseNuclearLayout(IntracellularPhysics& phys, float cellR,
                                          float condenseBlend = 0.0f) {
    condenseBlend = fmaxf(0.0f, fminf(1.0f, condenseBlend));
    float decondenseBlend = 1.0f - condenseBlend;
    float nucR = cellR * 0.45f;
    float helixR = nucR * 0.30f;
    float helixH = nucR * 1.60f;
    int turns = 8;
    std::map<int, int> dnaPairVisit;
    int poreIdx = 0;
    int spliceIdx = 0;
    int polIdx = 0;
    int dnaPolIdx = 0;
    auto clampToRadius = [](simd_float3 p, float maxR) -> simd_float3 {
        float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
        if (d <= maxR || d < 1.0e-5f) return p;
        float s = maxR / d;
        return {p.x * s, p.y * s, p.z * s};
    };

    for (auto& p : phys.particles) {
        if (!p.active) continue;

        p.confineCenter = {0, 0, 0};
        if (isNuclearParticleType(p.type)) p.confineRadius = nucR;

        if (p.type == PT_DNA_NODE) {
            int visit = dnaPairVisit[p.stateIndex]++;
            float t = (float)p.stateIndex / fmaxf((float)HBB_LENGTH, 1.0f);
            float angle = t * turns * 2.0f * M_PI + ((visit % 2) ? M_PI : 0.0f);
            simd_float3 h = {cosf(angle) * helixR,
                             (t - 0.5f) * helixH,
                             sinf(angle) * helixR};
            simd_float3 targetHome = {
                h.x + (p.home.x - h.x) * condenseBlend,
                h.y + (p.home.y - h.y) * condenseBlend,
                h.z + (p.home.z - h.z) * condenseBlend
            };
            p.home = targetHome;
            p.spawnPos = targetHome;
            p.homeK = 1.10f + condenseBlend * 0.90f;
            float posLerp = 0.08f + decondenseBlend * 0.57f;
            p.position = {
                p.position.x + (targetHome.x - p.position.x) * posLerp,
                p.position.y + (targetHome.y - p.position.y) * posLerp,
                p.position.z + (targetHome.z - p.position.z) * posLerp
            };
            p.position = clampToRadius(p.position, nucR * 0.88f);
            float velDamp = 0.25f + condenseBlend * 0.15f;
            p.velocity = {p.velocity.x * velDamp, p.velocity.y * velDamp, p.velocity.z * velDamp};
            p.mitosisHalf = -1;
        } else if (p.type == PT_DNA_POLYMERASE) {
            float a = (float)dnaPolIdx * 1.57f + 0.35f;
            simd_float3 h = {cosf(a) * nucR * 0.16f,
                             ((float)dnaPolIdx - 1.5f) * nucR * 0.05f,
                             sinf(a) * nucR * 0.16f};
            simd_float3 targetHome = {
                h.x + (p.home.x - h.x) * condenseBlend,
                h.y + (p.home.y - h.y) * condenseBlend,
                h.z + (p.home.z - h.z) * condenseBlend
            };
            p.home = targetHome;
            p.spawnPos = targetHome;
            p.homeK = 0.34f + condenseBlend * 0.16f;
            float posLerp = 0.14f + decondenseBlend * 0.62f;
            p.position = {
                p.position.x + (targetHome.x - p.position.x) * posLerp,
                p.position.y + (targetHome.y - p.position.y) * posLerp,
                p.position.z + (targetHome.z - p.position.z) * posLerp
            };
            p.velocity = {p.velocity.x * 0.18f, p.velocity.y * 0.18f, p.velocity.z * 0.18f};
            p.mitosisHalf = -1;
            dnaPolIdx++;
        } else if (p.type == PT_RNA_POLYMERASE) {
            float a = (float)polIdx * 1.9f;
            simd_float3 h = {sinf(a) * nucR * 0.22f,
                             ((float)polIdx - 1.0f) * nucR * 0.10f,
                             cosf(a) * nucR * 0.22f};
            simd_float3 targetHome = {
                h.x + (p.home.x - h.x) * condenseBlend,
                h.y + (p.home.y - h.y) * condenseBlend,
                h.z + (p.home.z - h.z) * condenseBlend
            };
            p.home = targetHome;
            p.spawnPos = targetHome;
            p.homeK = 0.30f + condenseBlend * 0.18f;
            float posLerp = 0.16f + decondenseBlend * 0.64f;
            p.position = {
                p.position.x + (targetHome.x - p.position.x) * posLerp,
                p.position.y + (targetHome.y - p.position.y) * posLerp,
                p.position.z + (targetHome.z - p.position.z) * posLerp
            };
            p.velocity = {p.velocity.x * 0.20f, p.velocity.y * 0.20f, p.velocity.z * 0.20f};
            p.mitosisHalf = -1;
            polIdx++;
        } else if (p.type == PT_SPLICEOSOME) {
            simd_float3 h = {(spliceIdx ? 0.16f : -0.16f) * cellR,
                             0.10f * cellR,
                             (spliceIdx ? 0.10f : -0.10f) * cellR};
            simd_float3 targetHome = {
                h.x + (p.home.x - h.x) * condenseBlend,
                h.y + (p.home.y - h.y) * condenseBlend,
                h.z + (p.home.z - h.z) * condenseBlend
            };
            p.home = targetHome;
            p.spawnPos = targetHome;
            p.homeK = 0.40f + condenseBlend * 0.10f;
            float posLerp = 0.18f + decondenseBlend * 0.60f;
            p.position = {
                p.position.x + (targetHome.x - p.position.x) * posLerp,
                p.position.y + (targetHome.y - p.position.y) * posLerp,
                p.position.z + (targetHome.z - p.position.z) * posLerp
            };
            p.velocity = {p.velocity.x * 0.18f, p.velocity.y * 0.18f, p.velocity.z * 0.18f};
            p.mitosisHalf = -1;
            spliceIdx++;
        } else if (p.type == PT_NUCLEAR_PORE) {
            float a = (float)poreIdx / 6.0f * 2.0f * M_PI;
            simd_float3 h = {cosf(a) * nucR * 0.95f,
                             sinf(a * 0.5f) * nucR * 0.3f,
                             sinf(a) * nucR * 0.95f};
            simd_float3 targetHome = {
                h.x + (p.home.x - h.x) * condenseBlend,
                h.y + (p.home.y - h.y) * condenseBlend,
                h.z + (p.home.z - h.z) * condenseBlend
            };
            p.home = targetHome;
            p.spawnPos = targetHome;
            p.homeK = 0.20f + condenseBlend * 0.08f;
            float posLerp = 0.18f + decondenseBlend * 0.70f;
            p.position = {
                p.position.x + (targetHome.x - p.position.x) * posLerp,
                p.position.y + (targetHome.y - p.position.y) * posLerp,
                p.position.z + (targetHome.z - p.position.z) * posLerp
            };
            p.velocity = {p.velocity.x * 0.12f, p.velocity.y * 0.12f, p.velocity.z * 0.12f};
            p.mitosisHalf = -1;
            poreIdx++;
        }
    }
}

static void retargetNuclearHomesForMitosis(IntracellularPhysics& phys, float cellR) {
    bool activeMitosis = gMitosis.active && gMitosis.particlesDuplicated &&
        gMitosis.phase >= MITO_PROMETAPHASE && gMitosis.phase <= MITO_CYTOKINESIS;

    if (!activeMitosis) {
        float condenseBlend = primaryPostMitoticCondenseBlend();
        if (condenseBlend > 0.02f) {
            restoreInterphaseNuclearLayout(phys, cellR, condenseBlend);
        } else {
            for (auto& p : phys.particles) {
                if (!p.active || !isNuclearParticleType(p.type)) continue;
                p.home = p.spawnPos;
                if (p.type == PT_DNA_NODE) p.homeK = 1.10f;
            }
        }
        return;
    }

    for (auto& p : phys.particles) {
        if (!p.active || !isNuclearParticleType(p.type)) continue;

        p.home = p.spawnPos;
        if (p.type == PT_DNA_NODE) p.homeK = 2.0f;
    }

    float targetX = mitosisNuclearTargetX(cellR);
    for (auto& p : phys.particles) {
        if (!p.active || !isNuclearParticleType(p.type)) continue;
        int half = particleMitosisHalf(p);

        if (p.type == PT_DNA_NODE) {
            p.home = mitosisDNAHomeTarget(p, cellR);
            float phaseBoost = 0.80f + mitosisCondenseBlend() * 0.65f;
            p.homeK = phaseBoost;
        } else {
            p.home = mitosisAuxNuclearHomeTarget(p, cellR);
            p.homeK = 0.22f + mitosisCondenseBlend() * 0.22f;
        }

        // Daughter-side nuclear confinement should move outward gradually
        // across the spindle phases instead of pinning everything to one line.
        simd_float3 daughterCenter = mitosisDaughterCenter(half, cellR);
        p.confineCenter = daughterCenter;
        float confineRadius = cellR * (0.28f - mitosisCondenseBlend() * 0.10f);
        if (gMitosis.phase >= MITO_TELOPHASE) {
            confineRadius = cellR * (0.20f + 0.18f * gMitosis.nucleusAlpha);
        }
        p.confineRadius = confineRadius;
        if (p.type == PT_DNA_NODE) {
            p.confineRadius *= (gMitosis.phase >= MITO_TELOPHASE) ? 0.88f : 0.72f;
        }
    }
}

static CellInterior splitMitosisDaughterInterior(
    const IntracellularPhysics& source,
    int half,
    simd_float3 localCenter,
    float daughterR,
    float postMitoticCondense = 0.0f
) {
    CellInterior daughter;
    daughter.phys.init({0, 0, 0}, daughterR);
    daughter.phys.params = source.params;
    daughter.phys.params.cellCenter = {0, 0, 0};
    daughter.phys.params.cellRadius = daughterR;
    daughter.phys.params.cellVelocity = {0, 0, 0};
    daughter.phys.skipCollisions = false;
    daughter.phys.attractionFields.clear();
    daughter.phys.particles.clear();
    daughter.phys.spatialHash.cellSize = fmaxf(daughterR * 0.08f, 0.1f);

    float sourceR = fmaxf(source.params.cellRadius, 0.1f);
    float daughterScale = daughterR / sourceR;
    float daughterNucR = daughterR * 0.45f;
    std::vector<int> indexMap(source.particles.size(), -1);
    std::vector<int> oldIndices;
    oldIndices.reserve(source.particles.size());

    for (int i = 0; i < (int)source.particles.size(); i++) {
        const auto& p = source.particles[i];
        if (!p.active) continue;
        // Each particle was tagged with mitosisHalf during duplication/tagging.
        // Nuclear particles were duplicated, so each half has its own full set.
        // Non-nuclear particles were tagged by X position.
        if (particleMitosisHalf(p) != half) continue;

        IntraParticle q = p;
        q.position = {p.position.x - localCenter.x,
                      p.position.y - localCenter.y,
                      p.position.z - localCenter.z};
        q.home = {p.home.x - localCenter.x,
                  p.home.y - localCenter.y,
                  p.home.z - localCenter.z};
        q.spawnPos = {p.spawnPos.x - localCenter.x,
                      p.spawnPos.y - localCenter.y,
                      p.spawnPos.z - localCenter.z};
        q.jobTarget = {p.jobTarget.x - localCenter.x,
                       p.jobTarget.y - localCenter.y,
                       p.jobTarget.z - localCenter.z};
        q.confineCenter = {0, 0, 0};
        q.confineRadius = isNuclearParticleType(q.type)
            ? daughterNucR
            : daughterR;
        q.radius *= daughterScale;
        q.baseRadius *= daughterScale;
        q.mass = fmaxf(q.radius * q.radius * 10.0f, 0.01f);
        q.tetherLen *= daughterScale;
        q.velocity = {p.velocity.x * 0.35f + ipRandGauss() * 0.015f,
                      p.velocity.y * 0.35f + ipRandGauss() * 0.015f,
                      p.velocity.z * 0.35f + ipRandGauss() * 0.015f};
        // Preserve the particle's job state so transport vesicles, attached
        // motors, and consumed particles keep their roles after division.
        // Only the relative coordinates (jobTarget, spawnPos, home) shift
        // with localCenter, which already happened above.
        q.job = p.job;
        q.jobProgress = p.jobProgress;
        q.jobTimer = p.jobTimer;
        q.jobSpeed = p.jobSpeed;
        // IMPORTANT: keep spawnPos as the organelle the particle originated at
        // (already shifted above). Do NOT overwrite it with the current home —
        // that would make consumed vesicles respawn at their last position
        // instead of at the Golgi/ER etc.
        q.mitosisHalf = -1;

        indexMap[i] = (int)daughter.phys.particles.size();
        daughter.phys.particles.push_back(q);
        oldIndices.push_back(i);
    }

    for (int i = 0; i < (int)daughter.phys.particles.size(); i++) {
        int oldIdx = oldIndices[i];
        auto& q = daughter.phys.particles[i];
        int oldTether = source.particles[oldIdx].tetherId;
        int oldJobTarget = source.particles[oldIdx].jobTargetId;
        q.tetherId = (oldTether >= 0 && oldTether < (int)indexMap.size()) ? indexMap[oldTether] : -1;
        q.jobTargetId = (oldJobTarget >= 0 && oldJobTarget < (int)indexMap.size()) ? indexMap[oldJobTarget] : -1;
        if (q.job == JOB_ATTACH && q.jobTargetId < 0) q.job = JOB_IDLE;
    }

    for (const auto& field : source.attractionFields) {
        AttractionField shifted = field;
        shifted.center = {field.center.x - localCenter.x,
                          field.center.y - localCenter.y,
                          field.center.z - localCenter.z};
        shifted.radius *= daughterScale;
        daughter.phys.attractionFields.push_back(shifted);
    }

    daughter.postMitoticCondense = fmaxf(0.0f, fminf(1.0f, postMitoticCondense));
    restoreInterphaseNuclearLayout(daughter.phys, daughterR, daughter.postMitoticCondense);
    daughter.initialized = true;
    return daughter;
}

static void prepareInteriorForPhysicalDivision(IntracellularPhysics& phys, float cellR) {
    // Idempotent guard: if THIS interior's nuclear aux has already been
    // duplicated (mitosisHalf==1 copies present), don't do it again — repeated
    // calls would duplicate nuclear particles exponentially. A per-interior
    // check is required because this function runs on the primary cell's
    // gIntraPhys AND on every background cell's interior; a global flag would
    // let one cell's call block another cell's legitimate prep.
    if (hasDuplicatedNuclearAuxForDivision(phys)) return;

    bool sPhaseDNAPrepared = hasPreparedSPhaseReplicatedDNAVisuals(phys);
    int originalCount = (int)phys.particles.size();

    // Tag existing particles first so cytoplasmic contents split by position.
    for (int i = 0; i < originalCount; i++) {
        auto& p = phys.particles[i];
        if (!p.active) continue;
        if (p.type == PT_DNA_NODE && p.mitosisHalf == 1) continue;
        bool isNuclear = isNuclearParticleType(p.type);
        if (p.type == PT_DNA_NODE) p.mitosisHalf = 0;
        else p.mitosisHalf = isNuclear ? 0 : ((p.position.x >= 0) ? 0 : 1);
    }

    std::vector<IntraParticle> source = phys.particles;
    std::vector<int> duplicateIndexMap(originalCount, -1);
    int duplicateCount = 0;
    for (int i = 0; i < originalCount; i++) {
        const auto& src = source[i];
        if (!src.active || !isNuclearParticleType(src.type)) continue;
        if (src.mitosisHalf == 1) continue;
        if (src.type == PT_DNA_NODE && sPhaseDNAPrepared) continue;
        duplicateCount++;
    }
    phys.particles.reserve(originalCount + duplicateCount);

    // Duplicate nuclear particles so each daughter inherits a full genome set.
    for (int i = 0; i < originalCount; i++) {
        const auto& src = source[i];
        if (!src.active || !isNuclearParticleType(src.type)) continue;
        if (src.mitosisHalf == 1) continue;
        if (src.type == PT_DNA_NODE && sPhaseDNAPrepared) continue;

        IntraParticle dup = src;
        dup.mitosisHalf = 1;
        dup.velocity = {0, 0, 0};

        if (dup.type == PT_DNA_NODE) {
            dup.position.z += cellR * 0.02f;
            dup.home.z += cellR * 0.02f;
        }

        dup.spawnPos = dup.home;
        duplicateIndexMap[i] = (int)phys.particles.size();
        phys.particles.push_back(dup);
    }

    for (int i = 0; i < originalCount; i++) {
        int dupIdx = duplicateIndexMap[i];
        if (dupIdx < 0) continue;
        auto& dup = phys.particles[dupIdx];
        int oldTether = source[i].tetherId;
        dup.tetherId =
            (oldTether >= 0 && oldTether < originalCount) ? duplicateIndexMap[oldTether] : -1;
        int oldJobTarget = source[i].jobTargetId;
        dup.jobTargetId =
            (oldJobTarget >= 0 && oldJobTarget < originalCount) ? duplicateIndexMap[oldJobTarget] : -1;
        if (dup.job == JOB_ATTACH && dup.jobTargetId < 0) dup.job = JOB_IDLE;
    }

    // If DNA sister copies were staged during S-phase, make sure they are all
    // visible before mitosis proper begins.
    if (sPhaseDNAPrepared) {
        for (auto& p : phys.particles) {
            if (p.type == PT_DNA_NODE && p.mitosisHalf == 1) {
                p.active = true;
            }
        }
    }
}

static void mitosisLogSnapshot(const char* event, float cellR) {
    mitosisLog("[%s] phase=%d furrow=%.2f poleSep=%.2f timer=%.3f\n",
        event, gMitosis.phase, gMitosis.furrowDepth, gMitosis.poleSeparation, gMitosis.phaseTimer);
    mitosisLog("  poleA=(%.3f,%.3f,%.3f) poleB=(%.3f,%.3f,%.3f)\n",
        gMitosis.poleA.x, gMitosis.poleA.y, gMitosis.poleA.z,
        gMitosis.poleB.x, gMitosis.poleB.y, gMitosis.poleB.z);
    mitosisLog("  cellR=%.3f daughterR=%.3f\n", cellR, cellR * 0.794f);
    mitosisLog("  dnaProgress=%.2f checkpoint=%d replProgress=%.3f forks=%d unresolved=%d chk1=%.2f ready=%d\n",
        gMitosis.dnaDuplicationProgress, gMitosis.dnaCheckpointPassed ? 1 : 0,
        gCDogma.replicationProgress, gCDogma.countActiveReplicationForks(),
        gCDogma.unresolvedReplicationErrors, gCDogma.chk1Signal,
        gCDogma.replicationReadyForM() ? 1 : 0);
    mitosisLog("  envBreak=%.2f spindleAsm=%.2f spindleDis=%.2f attach=%.2f align=%.2f checkpoint=%.2f chromatidSep=%.2f ring=%.2f\n",
        gMitosis.nuclearEnvelopeBreakdown, gMitosis.spindleAssembly,
        gMitosis.spindleDisassembly, gMitosis.kinetochoreAttachment,
        gMitosis.metaphaseAlignment, gMitosis.spindleCheckpoint,
        gMitosis.chromatidSeparation, gMitosis.contractileRingAssembly);
    // Particle stats per half
    int countH[2] = {}, activeCount = 0;
    float minX[2] = {1e9,1e9}, maxX[2] = {-1e9,-1e9};
    float avgDist[2] = {};
    for (auto& p : gIntraPhys.particles) {
        if (!p.active) continue;
        activeCount++;
        int h = particleMitosisHalf(p);
        countH[h]++;
        if (p.position.x < minX[h]) minX[h] = p.position.x;
        if (p.position.x > maxX[h]) maxX[h] = p.position.x;
        float d = sqrtf(p.position.x*p.position.x + p.position.y*p.position.y + p.position.z*p.position.z);
        avgDist[h] += d;
    }
    mitosisLog("  active=%d  half0: n=%d X=[%.2f,%.2f] avgD=%.2f  half1: n=%d X=[%.2f,%.2f] avgD=%.2f\n",
        activeCount,
        countH[0], minX[0], maxX[0], countH[0]>0 ? avgDist[0]/countH[0] : 0,
        countH[1], minX[1], maxX[1], countH[1]>0 ? avgDist[1]/countH[1] : 0);
    if (sDaughtersInitialized) {
        mitosisLog("  sDaughterA=(%.2f,%.2f,%.2f) sDaughterB=(%.2f,%.2f,%.2f) dR=%.3f\n",
            sDaughterPosA.x, sDaughterPosA.y, sDaughterPosA.z,
            sDaughterPosB.x, sDaughterPosB.y, sDaughterPosB.z, sDaughterR);
    }
    mitosisLog("  cells=%d\n", (int)gSim.cells.size());
}

// Rendering
static std::vector<CellInstance> gCellInstances;
static id<MTLBuffer> gCellBuffer = nil;
static std::map<std::string, LoadedMesh> gOrganelleMeshes;
static int gSelectedCell = 0;
static int gSelectedCellUid = -1;
static bool gFollowCell = true; // Track active single-cell division by default
static float gFollowZoom = 1.0f; // Zoom multiplier for follow mode (mouse-wheel driven)
static bool gSoloFocusCell = true;
static std::vector<int> gRenderCellSourceIndices;
static float gPostMitosisPairCameraTimer = 0.0f;
static int gPostMitosisPairUidA = -1;
static int gPostMitosisPairUidB = -1;
static simd_float3 gPostMitosisPairA = {0, 0, 0};
static simd_float3 gPostMitosisPairB = {0, 0, 0};
static float gPostMitosisPairRadius = 0.0f;

static float clamp01f(float v) {
    return fmaxf(0.0f, fminf(1.0f, v));
}

static float smooth01f(float v) {
    v = clamp01f(v);
    return v * v * (3.0f - 2.0f * v);
}

static simd_float3 normalizeVec3(simd_float3 v, simd_float3 fallback) {
    float d = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    if (d < 1.0e-5f) return fallback;
    float invD = 1.0f / d;
    return {v.x * invD, v.y * invD, v.z * invD};
}

static simd_float3 crossVec3(simd_float3 a, simd_float3 b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static constexpr float BIO_TOTAL_ADENYLATE_MILLIMOLAR = 3.9f;
static constexpr float BIO_ATP_MILLIMOLAR_AT_100_PERCENT = BIO_TOTAL_ADENYLATE_MILLIMOLAR;

static float simATPPercentToBioATP(float atpPercent) {
    float clamped = fmaxf(0.0f, fminf(100.0f, atpPercent));
    return BIO_ATP_MILLIMOLAR_AT_100_PERCENT * (clamped / 100.0f);
}

static float bioATPToSimPercent(float atpMillimolar) {
    float pct = (atpMillimolar / BIO_ATP_MILLIMOLAR_AT_100_PERCENT) * 100.0f;
    return fmaxf(0.0f, fminf(100.0f, pct));
}

static void syncPrimarySimATPIntoBioEngine() {
    if (!gCellBio.initialized || gSim.cells.empty()) return;

    auto& pool = gCellBio.metabolism.pool;
    float targetATP = simATPPercentToBioATP(gSim.cells[0].ATP);
    float remaining = fmaxf(0.0f, BIO_TOTAL_ADENYLATE_MILLIMOLAR - targetATP);
    float oldOther = pool.ADP + pool.AMP;

    pool.ATP = targetATP;
    if (oldOther > 1.0e-6f) {
        float scale = remaining / oldOther;
        pool.ADP *= scale;
        pool.AMP *= scale;
    } else {
        // Default mammalian pool split: ~0.8 mM ADP and ~0.1 mM AMP.
        pool.ADP = remaining * (0.8f / 0.9f);
        pool.AMP = remaining * (0.1f / 0.9f);
    }
}

static simd_float3 chromosomeSlotColor(int slot) {
    float hue = (float)slot / (float)fmax(MitosisState::NUM_CHROMO, 1);
    return {
        0.5f + 0.5f * sinf(hue * 6.28f),
        0.5f + 0.5f * sinf(hue * 6.28f + 2.09f),
        0.5f + 0.5f * sinf(hue * 6.28f + 4.19f)
    };
}

static float focusCellScore(const SimCell& c, int idx) {
    if (!c.alive) return -1.0e9f;

    float score = 0.0f;
    if (idx == gSelectedCell) score += 15.0f;
    if (idx == 0 && gMitosis.active) score += 1000.0f;
    if (c.program.mitosis.active) score += 550.0f;
    if (c.program.mitosis.postDivisionComplete()) score += 220.0f;
    if (c.divisionPending) score += 400.0f;
    score += c.phase * 120.0f;
    score += c.cycleProgress * 35.0f;
    score += c.ATP * 0.25f;
    score -= c.stress * 0.35f;
    score -= c.damageLevel * 90.0f;
    if (c.fate == SIM_FATE_PROLIF) score += 60.0f;
    if (c.fate == SIM_FATE_QUIESCENT) score -= 30.0f;
    if (c.fate == SIM_FATE_APOPTOTIC || c.apoptosisPhase > 0) score -= 1200.0f;
    if (c.necrotic) score -= 1500.0f;
    return score;
}

static int bestFocusCellIndex() {
    if (gSim.cells.empty()) return 0;

    int bestIdx = 0;
    float bestScore = focusCellScore(gSim.cells[0], 0);
    for (int i = 1; i < (int)gSim.cells.size(); i++) {
        float score = focusCellScore(gSim.cells[i], i);
        if (score > bestScore) {
            bestScore = score;
            bestIdx = i;
        }
    }
    return bestIdx;
}

static bool focusCellIsUsable(const SimCell& c) {
    if (!c.alive) return false;
    if (c.necrotic) return false;
    if (c.apoptosisPhase > 0) return false;
    if (c.fate == SIM_FATE_APOPTOTIC) return false;
    return true;
}

static void selectCellIndex(int idx) {
    if (gSim.cells.empty()) {
        gSelectedCell = 0;
        gSelectedCellUid = -1;
        return;
    }
    if (idx < 0) idx = 0;
    if (idx >= (int)gSim.cells.size()) idx = (int)gSim.cells.size() - 1;
    gSelectedCell = idx;
    gSelectedCellUid = gSim.cells[idx].cellUid;
}

// ── Directional cell switching ──────────────────────────────────────────
// When in follow mode the user can press W/A/S/D to jump to the nearest
// cell in that screen-space direction. We pick the cell whose displacement
// from the current tracked cell has the best-aligned projection onto the
// desired screen-space axis (camera-right for left/right, camera-forward
// for up/down).
//
// Scoring: for a direction vector d and candidate cell at world position
// p relative to current cell, score = max(0, (p·d) / |p|) - 0.2·(|p|/R).
// The angle term prefers cells aligned with the direction, and the
// distance term gently prefers nearby cells.
static int nearestCellInDirection(int currentIdx, simd_float3 dir_world) {
    if (gSim.cells.size() <= 1 || currentIdx < 0) return currentIdx;
    const auto& cur = gSim.cells[currentIdx];
    float dlen = sqrtf(dir_world.x*dir_world.x + dir_world.y*dir_world.y + dir_world.z*dir_world.z);
    if (dlen < 1e-4f) return currentIdx;
    simd_float3 d = {dir_world.x/dlen, dir_world.y/dlen, dir_world.z/dlen};

    int best = -1;
    float bestScore = -1e9f;
    for (int i = 0; i < (int)gSim.cells.size(); i++) {
        if (i == currentIdx) continue;
        const auto& c = gSim.cells[i];
        if (!c.alive) continue;
        simd_float3 v = {c.position.x - cur.position.x,
                         c.position.y - cur.position.y,
                         c.position.z - cur.position.z};
        float vlen = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
        if (vlen < 1e-4f) continue;
        float cos_angle = (v.x*d.x + v.y*d.y + v.z*d.z) / vlen;
        if (cos_angle < 0.1f) continue; // skip cells behind or perpendicular
        float score = cos_angle - 0.02f * vlen;
        if (score > bestScore) { bestScore = score; best = i; }
    }
    return best >= 0 ? best : currentIdx;
}

static int findCellIndexByUid(int uid) {
    if (uid < 0) return -1;
    for (int i = 0; i < (int)gSim.cells.size(); i++) {
        if (gSim.cells[i].cellUid == uid) return i;
    }
    return -1;
}

static bool trackedPostMitosisPairExists() {
    return gPostMitosisPairUidA >= 0 && gPostMitosisPairUidB >= 0;
}

static bool trackedPostMitosisPairActive() {
    if (!trackedPostMitosisPairExists()) return false;
    int idxA = findCellIndexByUid(gPostMitosisPairUidA);
    int idxB = findCellIndexByUid(gPostMitosisPairUidB);
    if (idxA < 0 || idxB < 0) return false;
    return gPostMitosisPairCameraTimer > 0.0f ||
           gSim.cells[idxA].postDivisionRecovery > 0.0f ||
           gSim.cells[idxB].postDivisionRecovery > 0.0f;
}

static bool postMitosisPairHoldActive() {
    return gPostMitosisPairCameraTimer > 0.0f && trackedPostMitosisPairExists();
}

static bool sourceIdxMatchesTrackedPostMitosisPair(int sourceIdx) {
    if (!trackedPostMitosisPairActive()) return false;
    if (sourceIdx < 0 || sourceIdx >= (int)gSim.cells.size()) return false;
    int uid = gSim.cells[sourceIdx].cellUid;
    return uid == gPostMitosisPairUidA || uid == gPostMitosisPairUidB;
}

static void refreshTrackedPostMitosisPairFromSimulation() {
    if (!trackedPostMitosisPairExists()) return;

    int idxA = findCellIndexByUid(gPostMitosisPairUidA);
    int idxB = findCellIndexByUid(gPostMitosisPairUidB);
    if (idxA < 0 || idxB < 0) {
        gPostMitosisPairCameraTimer = 0.0f;
        gPostMitosisPairUidA = -1;
        gPostMitosisPairUidB = -1;
        return;
    }

    const auto& a = gSim.cells[idxA];
    const auto& b = gSim.cells[idxB];
    gPostMitosisPairA = a.position;
    gPostMitosisPairB = b.position;
    gPostMitosisPairRadius = fmaxf(a.radius * a.size, b.radius * b.size);

    if (gPostMitosisPairCameraTimer <= 0.0f &&
        a.postDivisionRecovery <= 0.0f &&
        b.postDivisionRecovery <= 0.0f) {
        gPostMitosisPairUidA = -1;
        gPostMitosisPairUidB = -1;
    }
}

static int bestQueuedDivisionCellIndex() {
    int bestIdx = -1;
    float bestScore = -1.0e9f;
    for (int i = 0; i < (int)gSim.cells.size(); i++) {
        const auto& c = gSim.cells[i];
        bool readyForVisualDivision =
            c.divisionPending ||
            c.phase == 3 ||
            c.program.mitosis.active ||
            c.program.mitosis.postDivisionComplete();
        if (!c.alive || !readyForVisualDivision) continue;
        float score = focusCellScore(c, i);
        if (score > bestScore) {
            bestScore = score;
            bestIdx = i;
        }
    }
    return bestIdx;
}

static void promoteCellToPrimary(int idx) {
    if (idx <= 0 || idx >= (int)gSim.cells.size()) return;
    if (gMitosis.active) return;

    persistBoundDogmaState();
    persistBoundMitosisState();

    // Persist the current primary interior before demoting that cell.
    if (gIntraInitialized && !gSim.cells.empty()) {
        CellInterior storedPrimary;
        storedPrimary.initialized = true;
        storedPrimary.phys = gIntraPhys;
        gCellInteriors[gSim.cells[0].cellUid] = std::move(storedPrimary);
    }

    std::swap(gSim.cells[0], gSim.cells[idx]);

    // Load the promoted cell's real interior if we already have one. If not,
    // defer initialization to uploadCellInterior() this frame.
    auto it = gCellInteriors.find(gSim.cells[0].cellUid);
    if (it != gCellInteriors.end() && it->second.initialized) {
        gIntraPhys = it->second.phys;
        gIntraInitialized = true;
        gCellInteriors.erase(it);
    } else {
        gIntraInitialized = false;
    }

    bindDogmaStateToPrimaryCell(true);
    bindMitosisStateToPrimaryCell(true);

    selectCellIndex(0);
    if (!trackedPostMitosisPairActive()) {
        gPostMitosisPairCameraTimer = 0.0f;
        gPostMitosisPairUidA = -1;
        gPostMitosisPairUidB = -1;
    }
    if (gCellBio.initialized) {
        syncPrimarySimATPIntoBioEngine();
    }
}

static void ensurePrimaryCellBinding() {
    if (gMitosis.active || gSim.cells.empty()) return;

    int desiredIdx = findCellIndexByUid(gSelectedCellUid);
    if (desiredIdx < 0 ||
        desiredIdx >= (int)gSim.cells.size() ||
        !focusCellIsUsable(gSim.cells[desiredIdx])) {
        desiredIdx = bestFocusCellIndex();
    }
    if (desiredIdx < 0 || desiredIdx >= (int)gSim.cells.size()) return;

    // In follow mode, don't stay glued to an old selected UID if another
    // lineage has reached the replicated M-phase boundary and is waiting to
    // enter the primary visual mitosis pipeline. Without this, second-gen
    // cells can sit indefinitely in phase 3 with DNA=1524 while the camera
    // and primary slot remain attached to a less interesting cell.
    if (gSimMode == MODE_SINGLE_CELL && gFollowCell) {
        int readyIdx = bestQueuedDivisionCellIndex();
        if (readyIdx >= 0 && readyIdx < (int)gSim.cells.size()) {
            desiredIdx = readyIdx;
        } else {
            int bestIdx = bestFocusCellIndex();
            if (bestIdx >= 0 && bestIdx < (int)gSim.cells.size()) {
                float currentScore = focusCellScore(gSim.cells[desiredIdx], desiredIdx);
                float bestScore = focusCellScore(gSim.cells[bestIdx], bestIdx);
                if (bestScore > currentScore + 180.0f) {
                    desiredIdx = bestIdx;
                }
            }
        }
    }

    if (desiredIdx > 0) {
        promoteCellToPrimary(desiredIdx);
    } else {
        selectCellIndex(0);
        bindDogmaStateToPrimaryCell();
        bindMitosisStateToPrimaryCell();
        if (gCellBio.initialized) {
            syncPrimarySimATPIntoBioEngine();
        }
    }
}

static int activeFocusCellIndex() {
    if (gSim.cells.empty()) return 0;

    // Track selection by stable cell UID so dead-cell erases and daughter
    // insertions don't silently jump the focus to a different lineage.
    int uidIdx = findCellIndexByUid(gSelectedCellUid);
    if (uidIdx >= 0) {
        gSelectedCell = uidIdx;
    } else {
        selectCellIndex(bestFocusCellIndex());
    }

    if (gSimMode == MODE_SINGLE_CELL && gFollowCell) {
        // While visual mitosis is active, keep the focus on the visual lineage.
        // Daughter A preserves cell[0]'s UID across finalization.
        if (gMitosis.active && !gSim.cells.empty() && focusCellIsUsable(gSim.cells[0])) {
            selectCellIndex(0);
            return gSelectedCell;
        }

        if (!focusCellIsUsable(gSim.cells[gSelectedCell])) {
            selectCellIndex(bestFocusCellIndex());
        } else {
            gSelectedCellUid = gSim.cells[gSelectedCell].cellUid;
        }
    } else {
        gSelectedCellUid = gSim.cells[gSelectedCell].cellUid;
    }
    return gSelectedCell;
}

static bool soloFocusEnabled() {
    return gSimMode == MODE_SINGLE_CELL && gSoloFocusCell;
}

static bool shouldRenderSourceCellIndex(int sourceIdx) {
    if (!soloFocusEnabled()) return true;
    if (gSimMode == MODE_SINGLE_CELL && sourceIdx >= 0 && (int)gSim.cells.size() <= 2) return true;
    int focusIdx = activeFocusCellIndex();
    if (sourceIdx == focusIdx) return true;
    if (sourceIdx >= 0 && sourceIdxMatchesTrackedPostMitosisPair(sourceIdx)) return true;
    if (sourceIdx == -2 && focusIdx == 0 &&
        gMitosis.active && gMitosis.phase == MITO_COMPLETE) return true;
    return false;
}

static bool shouldRenderCellShellSourceIndex(int sourceIdx) {
    // Solo-focus mode is meant to reduce detailed interior clutter and keep the
    // camera locked to one lineage, not to make real live cells disappear from
    // the colony shell renderer. Keep all simulation-backed cells visible.
    if (sourceIdx >= 0) return true;
    return shouldRenderSourceCellIndex(sourceIdx);
}

static int renderSourceIndex(int renderIdx) {
    if (renderIdx >= 0 && renderIdx < (int)gRenderCellSourceIndices.size()) {
        return gRenderCellSourceIndices[renderIdx];
    }
    return renderIdx;
}
static bool gShowMitosisDebugOverlay = true;
static bool gDebugOverlayShowLabels = true;
static bool gDebugOverlayShowNucleusMarkers = true;
static bool gDebugOverlayShowCollisionRings = true;
static bool gDebugOverlayShowWarnings = true;

// Time-series recording (growing vectors, no circular buffer issues)
struct TimeSeriesData {
    std::vector<float> time;        // bio-hours
    std::vector<float> population;
    std::vector<float> proliferating;
    std::vector<float> quiescent;
    std::vector<float> apoptotic;
    std::vector<float> necrotic;
    std::vector<float> avgATP;
    std::vector<float> avgStress;
    std::vector<float> glycolyticPct;
    std::vector<float> phaseG1, phaseS, phaseG2, phaseM;
    std::vector<float> divisions;   // cumulative
    std::vector<float> deaths;      // cumulative
    float sampleTimer = 0;

    void sample(const Simulation& sim) {
        float bioH = sim.bioTime / 3600.0f;
        time.push_back(bioH);
        population.push_back((float)sim.statAlive);
        proliferating.push_back((float)sim.statProlif);
        quiescent.push_back((float)sim.statQuiescent);
        apoptotic.push_back((float)sim.statApoptotic);
        necrotic.push_back((float)sim.statNecrotic);
        avgATP.push_back(sim.statAvgATP);
        // Compute avg stress
        float sumStress = 0;
        for (auto& c : sim.cells) if (c.alive) sumStress += c.stress;
        avgStress.push_back(sim.statAlive > 0 ? sumStress / sim.statAlive : 0);
        glycolyticPct.push_back(sim.statAlive > 0 ? (float)sim.statGlycolytic / sim.statAlive * 100 : 0);
        float total = fmaxf(1, (float)sim.statAlive);
        phaseG1.push_back(sim.statPhases[0] / total * 100);
        phaseS.push_back(sim.statPhases[1] / total * 100);
        phaseG2.push_back(sim.statPhases[2] / total * 100);
        phaseM.push_back(sim.statPhases[3] / total * 100);
        divisions.push_back((float)sim.statDivisions);
        deaths.push_back((float)sim.statDeaths);
    }

    int count() const { return (int)time.size(); }

    void exportCSV(const std::string& path) const {
        FILE* f = fopen(path.c_str(), "w");
        if (!f) { printf("[Export] Failed to write %s\n", path.c_str()); return; }
        fprintf(f, "bio_time_h,population,proliferating,quiescent,apoptotic,necrotic,"
                   "avg_ATP,avg_stress,glycolytic_pct,"
                   "phase_G1_pct,phase_S_pct,phase_G2_pct,phase_M_pct,"
                   "cumulative_divisions,cumulative_deaths\n");
        for (int i = 0; i < count(); i++) {
            fprintf(f, "%.4f,%g,%g,%g,%g,%g,%.2f,%.2f,%.1f,%.1f,%.1f,%.1f,%.1f,%g,%g\n",
                    time[i], population[i], proliferating[i], quiescent[i],
                    apoptotic[i], necrotic[i], avgATP[i], avgStress[i], glycolyticPct[i],
                    phaseG1[i], phaseS[i], phaseG2[i], phaseM[i],
                    divisions[i], deaths[i]);
        }
        fclose(f);
        printf("[Export] Saved %d rows to %s\n", count(), path.c_str());
    }


    // Per-cell snapshot export
    void exportCellSnapshot(const std::string& path, const Simulation& sim) const {
        FILE* f = fopen(path.c_str(), "w");
        if (!f) { printf("[Export] Failed to write %s\n", path.c_str()); return; }
        fprintf(f, "cell_id,pos_x,pos_z,phase,fate,ATP,stress,ROS,biomass,damage,"
                   "telomere,generation,clone_id,glycolytic,mito_health,mito_potential,"
                   "CycD,Rb,E2F,CycE,CycA,CycB,p21,pressure,necrotic,senescent\n");
        for (int i = 0; i < (int)sim.cells.size(); i++) {
            const auto& c = sim.cells[i];
            if (!c.alive) continue;
            const char* phaseN[] = {"G1","S","G2","M"};
            const char* fateN[] = {"undetermined","proliferating","quiescent","apoptotic"};
            fprintf(f, "%d,%.2f,%.2f,%s,%s,%.1f,%.1f,%.1f,%.3f,%.4f,"
                       "%.0f,%d,%d,%d,%.3f,%.1f,"
                       "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%d,%d\n",
                    i, c.position.x, c.position.z,
                    phaseN[c.phase], fateN[c.fate],
                    c.ATP, c.stress, c.ROS, c.biomass, c.damageLevel,
                    c.telomere, c.generation, c.cloneId,
                    c.glycolytic ? 1 : 0, c.mitoHealth, c.mitoPotential,
                    c.cdk.CycD, c.cdk.Rb, c.cdk.E2F, c.cdk.CycE, c.cdk.CycA, c.cdk.CycB, c.cdk.p21,
                    c.localPressure, c.necrotic ? 1 : 0, c.senescent ? 1 : 0);
        }
        fclose(f);
        printf("[Export] Cell snapshot: %d cells to %s\n", sim.statAlive, path.c_str());
    }
};
static TimeSeriesData gTS;
static std::string gExportDir;

// Per-organelle uniforms
struct OrgUniforms {
    simd_float4x4 viewProjection;
    simd_float4x4 model;
    simd_float3   cameraPos;
    float          time;
    simd_float3   lightDir;
    float          pad0;
    simd_float3   baseColor;
    float          emissiveIntensity;
    simd_float3   emissiveColor;
    float          pad1;
    simd_float3   cellCenter;
    float          cellRadius;
    float          furrowDepth;    // matches CellRender furrow for organelle clamping
    // Selected-cell glow pass: when > 0, shader emits a strong sinusoidal
    // pulse using baseColor that pierces through other geometry (drawn
    // with depth test always-pass).  Vanilla draw uses 0.
    float          glowBoost;
    float          pad3;
    float          pad4;
};

// ── Matrix helpers ──────────────────────────────────────────────────────
static simd_float4x4 mat4_translation(simd_float3 t) {
    return (simd_float4x4){{
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {t.x,t.y,t.z,1}
    }};
}
static simd_float4x4 mat4_scale(float s) {
    return (simd_float4x4){{
        {s,0,0,0}, {0,s,0,0}, {0,0,s,0}, {0,0,0,1}
    }};
}
static simd_float4x4 mat4_rotateY(float a) {
    float c=cosf(a),s=sinf(a);
    return (simd_float4x4){{{c,0,s,0},{0,1,0,0},{-s,0,c,0},{0,0,0,1}}};
}
static simd_float4x4 mat4_rotateX(float a) {
    float c=cosf(a),s=sinf(a);
    return (simd_float4x4){{{1,0,0,0},{0,c,-s,0},{0,s,c,0},{0,0,0,1}}};
}
static simd_float4x4 mat4_rotateZ(float a) {
    float c=cosf(a),s=sinf(a);
    return (simd_float4x4){{{c,-s,0,0},{s,c,0,0},{0,0,1,0},{0,0,0,1}}};
}

static int countDNAInPhysics(const IntracellularPhysics& phys) {
    int dna = 0;
    for (const auto& p : phys.particles) {
        if (p.active && p.type == PT_DNA_NODE) dna++;
    }
    return dna;
}

static bool computeNucleusLocalCenter(const IntracellularPhysics& phys, simd_float3& out) {
    simd_float3 sum = {0, 0, 0};
    float weight = 0.0f;
    for (const auto& p : phys.particles) {
        if (!p.active || !isNuclearParticleType(p.type)) continue;
        float w = (p.type == PT_DNA_NODE) ? 4.0f : 1.0f;
        sum.x += p.position.x * w;
        sum.y += p.position.y * w;
        sum.z += p.position.z * w;
        weight += w;
    }
    if (weight <= 0.0f) return false;
    out = {sum.x / weight, sum.y / weight, sum.z / weight};
    return true;
}

static bool projectWorldToScreen(simd_float3 world, const simd_float4x4& viewProjection, ImVec2& out) {
    simd_float4 clip = simd_mul(viewProjection, simd_make_float4(world.x, world.y, world.z, 1.0f));
    if (clip.w <= 0.001f) return false;
    float invW = 1.0f / clip.w;
    float ndcX = clip.x * invW;
    float ndcY = clip.y * invW;
    if (ndcX < -1.3f || ndcX > 1.3f || ndcY < -1.3f || ndcY > 1.3f) return false;
    ImVec2 display = ImGui::GetIO().DisplaySize;
    out.x = (ndcX * 0.5f + 0.5f) * display.x;
    out.y = (1.0f - (ndcY * 0.5f + 0.5f)) * display.y;
    return true;
}

static float projectScreenRadius(simd_float3 center, float radius, const simd_float4x4& viewProjection) {
    ImVec2 a, b;
    if (!projectWorldToScreen(center, viewProjection, a)) return -1.0f;
    simd_float3 edge = {center.x + radius, center.y, center.z};
    if (!projectWorldToScreen(edge, viewProjection, b)) return -1.0f;
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    return sqrtf(dx * dx + dy * dy);
}

static void drawOverlayLabel(ImDrawList* dl, ImVec2 anchor, ImU32 accent, const char* text) {
    ImVec2 size = ImGui::CalcTextSize(text);
    ImVec2 pad = {6.0f, 4.0f};
    ImVec2 min = {anchor.x - size.x * 0.5f - pad.x, anchor.y - size.y - pad.y * 2.0f};
    ImVec2 max = {anchor.x + size.x * 0.5f + pad.x, anchor.y};
    dl->AddRectFilled(min, max, IM_COL32(3, 9, 20, 220), 6.0f);
    dl->AddRect(min, max, accent, 6.0f, 0, 1.5f);
    dl->AddText({min.x + pad.x, min.y + pad.y}, IM_COL32(230, 240, 255, 255), text);
}

static void drawCrossMarker(ImDrawList* dl, ImVec2 p, float r, ImU32 col) {
    dl->AddCircleFilled(p, 3.0f, col);
    dl->AddLine({p.x - r, p.y}, {p.x + r, p.y}, col, 1.5f);
    dl->AddLine({p.x, p.y - r}, {p.x, p.y + r}, col, 1.5f);
}

static void drawMitosisDebugOverlay(const simd_float4x4& viewProjection) {
    if (!gShowMitosisDebugOverlay || gSim.cells.empty()) return;

    ImDrawList* dl = ImGui::GetForegroundDrawList();
    auto drawCellDebug = [&](int idx, const SimCell& c, const IntracellularPhysics* phys) {
        if (!shouldRenderSourceCellIndex(idx)) return;
        ImVec2 cellScreen;
        if (!projectWorldToScreen(c.position, viewProjection, cellScreen)) return;

        float cellRadius = c.radius * c.size;
        float screenRadius = projectScreenRadius(c.position, cellRadius, viewProjection);
        if (gDebugOverlayShowCollisionRings && screenRadius > 3.0f) {
            ImU32 ringColor = (idx == 0) ? IM_COL32(90, 190, 255, 220) : IM_COL32(110, 255, 180, 220);
            dl->AddCircle(cellScreen, screenRadius, ringColor, 48, 1.8f);
        }

        bool hasInterior = (phys != nullptr);
        int dnaCount = hasInterior ? countDNAInPhysics(*phys) : -1;
        simd_float3 nucleusWorld = c.position;
        bool hasNucleus = false;
        if (hasInterior) {
            simd_float3 nucleusLocal = {0, 0, 0};
            if (computeNucleusLocalCenter(*phys, nucleusLocal)) {
                nucleusWorld = {c.position.x + nucleusLocal.x,
                                c.position.y + nucleusLocal.y,
                                c.position.z + nucleusLocal.z};
                hasNucleus = true;
            }
        }

        if (gDebugOverlayShowNucleusMarkers && hasNucleus) {
            ImVec2 nucleusScreen;
            if (projectWorldToScreen(nucleusWorld, viewProjection, nucleusScreen)) {
                ImU32 nucleusColor = (dnaCount > 0) ? IM_COL32(255, 210, 90, 255) : IM_COL32(255, 90, 90, 255);
                drawCrossMarker(dl, nucleusScreen, 7.0f, nucleusColor);
                dl->AddLine(cellScreen, nucleusScreen, IM_COL32(255, 255, 255, 90), 1.0f);
            }
        }

        if (gDebugOverlayShowLabels) {
            char text[256];
            const char* daughterTag = "";
            if (gSimMode == MODE_SINGLE_CELL) {
                if (idx == 0) daughterTag = "A ";
                else if (idx == 1) daughterTag = "B ";
            }
            // Biological counts formatted in scientific notation so the
            // full order-of-magnitude reads cleanly (ribosomes 10⁶,
            // ATP 10⁹, water 10¹⁴ — Milo 2013 / BioNumbers).
            auto fmtSci = [](double v, char* buf, size_t n) {
                if (v >= 1e12)      snprintf(buf, n, "%.1fT", v / 1e12);
                else if (v >= 1e9)  snprintf(buf, n, "%.1fG", v / 1e9);
                else if (v >= 1e6)  snprintf(buf, n, "%.1fM", v / 1e6);
                else if (v >= 1e3)  snprintf(buf, n, "%.1fk", v / 1e3);
                else                snprintf(buf, n, "%.0f",  v);
            };
            char sATP[16], sGlu[16], sRib[16], sTRNA[16], sMRNA[16], sAA[16];
            fmtSci(c.atpMolecules,     sATP,  sizeof(sATP));
            fmtSci(c.glucoseMolecules, sGlu,  sizeof(sGlu));
            fmtSci(c.ribosomeCount,    sRib,  sizeof(sRib));
            fmtSci(c.tRNACount,        sTRNA, sizeof(sTRNA));
            fmtSci(c.mRNACount,        sMRNA, sizeof(sMRNA));
            fmtSci(c.aaMolecules,      sAA,   sizeof(sAA));
            snprintf(text, sizeof(text),
                     "%sc%d u%d g%d p%d | DNA %s%d mito %d | "
                     "ATP %s gluc %s rib %s tRNA %s mRNA %s aa %s",
                     daughterTag, idx, c.cellUid, c.generation, c.phase,
                     hasInterior ? "" : "?", hasInterior ? dnaCount : 0,
                     c.mitoCount,
                     sATP, sGlu, sRib, sTRNA, sMRNA, sAA);
            ImU32 accent = IM_COL32(110, 210, 255, 230);
            if (!hasInterior) accent = IM_COL32(255, 90, 90, 230);
            else if (dnaCount <= 0) accent = IM_COL32(255, 150, 70, 230);
            drawOverlayLabel(dl, {cellScreen.x, cellScreen.y - fmaxf(screenRadius, 20.0f) - 10.0f}, accent, text);
        }

        if (gDebugOverlayShowWarnings && (!hasInterior || dnaCount <= 0)) {
            const char* warning = hasInterior ? "DNA MISSING" : "INTERIOR MISSING";
            drawOverlayLabel(dl, {cellScreen.x, cellScreen.y + fmaxf(screenRadius, 20.0f) + 26.0f},
                             IM_COL32(255, 80, 80, 240), warning);
        }
    };

    for (int i = 0; i < (int)gSim.cells.size(); i++) {
        const SimCell& c = gSim.cells[i];
        if (!c.alive) continue;
        const IntracellularPhysics* phys = nullptr;
        if (i == 0) {
            phys = gIntraInitialized ? &gIntraPhys : nullptr;
        } else {
            auto it = gCellInteriors.find(c.cellUid);
            if (it != gCellInteriors.end() && it->second.initialized) {
                phys = &it->second.phys;
            }
        }
        drawCellDebug(i, c, phys);
    }

    // Preview labels removed. MITO_COMPLETE is a single-frame state;
    // both daughters exist as real cells in gSim.cells and are labeled by
    // the normal drawCellDebug() call above.
}

static float organelleVecLength(simd_float3 v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

static float organelleDot(simd_float3 a, simd_float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static simd_float3 organelleNormalizeOr(simd_float3 v, simd_float3 fallback) {
    float len = organelleVecLength(v);
    if (len > 0.0001f) {
        float inv = 1.0f / len;
        return {v.x * inv, v.y * inv, v.z * inv};
    }
    return fallback;
}

static uint32_t organelleSeedBits(int key, int slot) {
    uint32_t x = (uint32_t)(key * 1103515245u) ^ (uint32_t)(slot * 2246822519u);
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

static float organelleSeed01(int key, int slot) {
    return (float)(organelleSeedBits(key, slot) & 0x00FFFFFFu) / 16777215.0f;
}

static float organelleSeedSigned(int key, int slot) {
    return organelleSeed01(key, slot) * 2.0f - 1.0f;
}

static simd_float3 organelleNoiseVec(int key, int slot, float t) {
    float p0 = organelleSeed01(key, slot * 3 + 0) * 2.0f * M_PI;
    float p1 = organelleSeed01(key, slot * 3 + 1) * 2.0f * M_PI;
    float p2 = organelleSeed01(key, slot * 3 + 2) * 2.0f * M_PI;
    return {
        sinf(t * (0.63f + 0.07f * organelleSeed01(key, slot + 11)) + p0) +
            0.45f * cosf(t * 1.31f + p1),
        cosf(t * (0.57f + 0.05f * organelleSeed01(key, slot + 17)) + p1) +
            0.35f * sinf(t * 1.17f + p2),
        sinf(t * (0.71f + 0.06f * organelleSeed01(key, slot + 23)) + p2) +
            0.40f * cosf(t * 1.09f + p0)
    };
}

static void confineOrganelleBody(simd_float3& pos, simd_float3& vel,
                                 float minR, float maxR) {
    simd_float3 fallback = {1, 0, 0};
    simd_float3 dir = organelleNormalizeOr(pos, fallback);
    float len = organelleVecLength(pos);

    if (len < minR) {
        pos = {dir.x * minR, dir.y * minR, dir.z * minR};
        float vn = organelleDot(vel, dir);
        if (vn < 0.0f) {
            vel = {vel.x - dir.x * vn * 1.5f,
                   vel.y - dir.y * vn * 1.5f,
                   vel.z - dir.z * vn * 1.5f};
        }
    } else if (len > maxR) {
        pos = {dir.x * maxR, dir.y * maxR, dir.z * maxR};
        float vn = organelleDot(vel, dir);
        if (vn > 0.0f) {
            vel = {vel.x - dir.x * vn * 1.6f,
                   vel.y - dir.y * vn * 1.6f,
                   vel.z - dir.z * vn * 1.6f};
        }
    }
}

static void separateOrganelleBodies(simd_float3& posA, simd_float3& velA,
                                    simd_float3& posB, simd_float3& velB,
                                    float minDist) {
    simd_float3 delta = {posA.x - posB.x, posA.y - posB.y, posA.z - posB.z};
    float len = organelleVecLength(delta);
    if (len < 0.0001f) {
        delta = {minDist, 0, 0};
        len = minDist;
    }
    if (len >= minDist) return;

    simd_float3 fallback = {1, 0, 0};
    simd_float3 dir = organelleNormalizeOr(delta, fallback);
    float push = (minDist - len) * 0.5f;
    posA = {posA.x + dir.x * push, posA.y + dir.y * push, posA.z + dir.z * push};
    posB = {posB.x - dir.x * push, posB.y - dir.y * push, posB.z - dir.z * push};
    velA = {velA.x + dir.x * push * 2.0f, velA.y + dir.y * push * 2.0f, velA.z + dir.z * push * 2.0f};
    velB = {velB.x - dir.x * push * 2.0f, velB.y - dir.y * push * 2.0f, velB.z - dir.z * push * 2.0f};
}

static void initializeOrganelleMotionState(int cellKey, OrganelleMotionState& state) {
    state = {};
    state.initialized = true;
    state.lastUpdateTime = -1.0f;

    OrganelleConfigs base = OrganelleConfigs::defaults();
    state.golgiPos = {
        base.golgi.position.x + organelleSeedSigned(cellKey, 1) * 0.10f,
        organelleSeedSigned(cellKey, 2) * 0.05f,
        base.golgi.position.z + organelleSeedSigned(cellKey, 3) * 0.08f
    };
    confineOrganelleBody(state.golgiPos, state.golgiVel, 0.28f, 0.52f);

    for (int m = 0; m < 3; m++) {
        simd_float3 seedVec = {
            organelleSeedSigned(cellKey, 20 + m * 7 + 0),
            organelleSeedSigned(cellKey, 20 + m * 7 + 1) * 0.55f,
            organelleSeedSigned(cellKey, 20 + m * 7 + 2)
        };
        simd_float3 fallback = {1, 0, 0};
        simd_float3 dir = organelleNormalizeOr(
            seedVec, fallback);
        float radial = 0.42f + organelleSeed01(cellKey, 20 + m * 7 + 3) * 0.14f;
        state.mitoPos[m] = {dir.x * radial, dir.y * radial, dir.z * radial};
        state.mitoVel[m] = {
            organelleSeedSigned(cellKey, 20 + m * 7 + 4) * 0.04f,
            organelleSeedSigned(cellKey, 20 + m * 7 + 5) * 0.03f,
            organelleSeedSigned(cellKey, 20 + m * 7 + 6) * 0.04f
        };
        confineOrganelleBody(state.mitoPos[m], state.mitoVel[m], 0.38f, 0.62f);
    }
}

static OrganelleMotionState& updateOrganelleMotionState(int cellKey, float time) {
    auto& state = gOrganelleMotionStates[cellKey];
    if (!state.initialized) {
        initializeOrganelleMotionState(cellKey, state);
    }

    float dt = 1.0f / 60.0f;
    if (state.lastUpdateTime >= 0.0f) {
        dt = time - state.lastUpdateTime;
        if (dt < 0.0f) dt = 0.0f;
        if (dt > 0.05f) dt = 0.05f;
    }
    state.lastUpdateTime = time;
    if (dt <= 0.00001f) return state;

    simd_float3 golgiBiasRaw = {
        0.15f + organelleSeedSigned(cellKey, 31) * 0.25f,
        organelleSeedSigned(cellKey, 32) * 0.12f,
        -0.35f + organelleSeedSigned(cellKey, 33) * 0.18f
    };
    simd_float3 golgiFallback = {0, 0, -1};
    simd_float3 golgiBias = organelleNormalizeOr(golgiBiasRaw, golgiFallback);
    float golgiRadius = 0.34f + 0.08f * (0.5f + 0.5f *
        sinf(time * 0.21f + organelleSeed01(cellKey, 34) * 2.0f * M_PI));
    simd_float3 golgiTarget = {
        golgiBias.x * golgiRadius,
        golgiBias.y * golgiRadius,
        golgiBias.z * golgiRadius
    };
    simd_float3 golgiNoise = organelleNoiseVec(cellKey, 4, time * 1.2f);
    state.golgiVel = {
        state.golgiVel.x + (golgiTarget.x - state.golgiPos.x) * (1.7f * dt) + golgiNoise.x * (0.11f * dt),
        state.golgiVel.y + (golgiTarget.y - state.golgiPos.y) * (1.7f * dt) + golgiNoise.y * (0.11f * dt),
        state.golgiVel.z + (golgiTarget.z - state.golgiPos.z) * (1.7f * dt) + golgiNoise.z * (0.11f * dt)
    };
    float golgiDamp = powf(0.86f, dt * 60.0f);
    state.golgiVel = {state.golgiVel.x * golgiDamp,
                      state.golgiVel.y * golgiDamp,
                      state.golgiVel.z * golgiDamp};
    state.golgiPos = {state.golgiPos.x + state.golgiVel.x * dt,
                      state.golgiPos.y + state.golgiVel.y * dt,
                      state.golgiPos.z + state.golgiVel.z * dt};
    confineOrganelleBody(state.golgiPos, state.golgiVel, 0.30f, 0.54f);

    for (int m = 0; m < 3; m++) {
        simd_float3 drift = organelleNoiseVec(cellKey, 10 + m, time * (1.4f + 0.15f * (float)m));
        simd_float3 radialRaw = {
            state.mitoPos[m].x + drift.x * 0.20f,
            state.mitoPos[m].y + drift.y * 0.20f,
            state.mitoPos[m].z + drift.z * 0.20f
        };
        simd_float3 radialFallback = {
            (m == 0) ? 1.0f : -0.5f,
            (m == 1) ? 0.2f : 0.0f,
            (m == 2) ? -0.8f : 0.6f
        };
        simd_float3 radialDir = organelleNormalizeOr(
            radialRaw, radialFallback);
        float radialTarget = 0.44f + 0.10f * organelleSeed01(cellKey, 50 + m) +
            0.04f * sinf(time * (0.35f + 0.05f * (float)m) + organelleSeed01(cellKey, 60 + m) * 2.0f * M_PI);
        simd_float3 target = {
            radialDir.x * radialTarget,
            radialDir.y * radialTarget,
            radialDir.z * radialTarget
        };
        state.mitoVel[m] = {
            state.mitoVel[m].x + (target.x - state.mitoPos[m].x) * (1.1f * dt) + drift.x * (0.16f * dt),
            state.mitoVel[m].y + (target.y - state.mitoPos[m].y) * (1.1f * dt) + drift.y * (0.16f * dt),
            state.mitoVel[m].z + (target.z - state.mitoPos[m].z) * (1.1f * dt) + drift.z * (0.16f * dt)
        };
        float mitoDamp = powf(0.90f, dt * 60.0f);
        state.mitoVel[m] = {state.mitoVel[m].x * mitoDamp,
                            state.mitoVel[m].y * mitoDamp,
                            state.mitoVel[m].z * mitoDamp};
        state.mitoPos[m] = {state.mitoPos[m].x + state.mitoVel[m].x * dt,
                            state.mitoPos[m].y + state.mitoVel[m].y * dt,
                            state.mitoPos[m].z + state.mitoVel[m].z * dt};
        confineOrganelleBody(state.mitoPos[m], state.mitoVel[m], 0.40f, 0.63f);
    }

    for (int a = 0; a < 3; a++) {
        separateOrganelleBodies(state.golgiPos, state.golgiVel,
                                state.mitoPos[a], state.mitoVel[a], 0.23f);
        for (int b = a + 1; b < 3; b++) {
            separateOrganelleBodies(state.mitoPos[a], state.mitoVel[a],
                                    state.mitoPos[b], state.mitoVel[b], 0.24f);
        }
    }

    confineOrganelleBody(state.golgiPos, state.golgiVel, 0.30f, 0.54f);
    for (int m = 0; m < 3; m++) {
        confineOrganelleBody(state.mitoPos[m], state.mitoVel[m], 0.40f, 0.63f);
    }

    return state;
}

static OrganelleMotionState makeSplitOrganelleMotionState(const OrganelleMotionState& parent,
                                                          float xBias) {
    OrganelleMotionState out = parent;
    out.initialized = true;
    out.lastUpdateTime = -1.0f;

    float posScale = 0.86f;
    out.golgiPos = {
        parent.golgiPos.x * posScale + xBias * 0.10f,
        parent.golgiPos.y * posScale,
        parent.golgiPos.z * posScale
    };
    out.golgiVel = {
        parent.golgiVel.x * 0.35f + xBias * 0.05f,
        parent.golgiVel.y * 0.35f,
        parent.golgiVel.z * 0.35f
    };
    confineOrganelleBody(out.golgiPos, out.golgiVel, 0.30f, 0.54f);

    for (int m = 0; m < 3; m++) {
        float mitoBias = xBias * (0.07f + 0.02f * (float)m);
        out.mitoPos[m] = {
            parent.mitoPos[m].x * posScale + mitoBias,
            parent.mitoPos[m].y * posScale,
            parent.mitoPos[m].z * posScale
        };
        out.mitoVel[m] = {
            parent.mitoVel[m].x * 0.35f + mitoBias * 0.4f,
            parent.mitoVel[m].y * 0.35f,
            parent.mitoVel[m].z * 0.35f
        };
        confineOrganelleBody(out.mitoPos[m], out.mitoVel[m], 0.40f, 0.63f);
    }

    return out;
}

static void splitOrganelleMotionState(int parentUid, int daughterAUid, int daughterBUid) {
    OrganelleMotionState parentState = {};
    auto it = gOrganelleMotionStates.find(parentUid);
    if (it != gOrganelleMotionStates.end() && it->second.initialized) {
        parentState = it->second;
    } else {
        initializeOrganelleMotionState(parentUid, parentState);
    }

    gOrganelleMotionStates[daughterAUid] = makeSplitOrganelleMotionState(parentState, 0.12f);
    gOrganelleMotionStates[daughterBUid] = makeSplitOrganelleMotionState(parentState, -0.12f);
}

static int organelleMotionKeyForRenderIndex(int i) {
    int sourceIdx = renderSourceIndex(i);
    if (sourceIdx >= 0 && sourceIdx < (int)gSim.cells.size()) return gSim.cells[sourceIdx].cellUid;
    if (sourceIdx == -2 && gMitosis.active && gMitosis.phase == MITO_COMPLETE) return -1002;
    if (gMitosis.active && gMitosis.phase == MITO_COMPLETE) return -1000 - i;
    return -2000 - i;
}

static OrganelleConfigs animatedOrganelleConfigsForCell(int cellKey, float time) {
    OrganelleConfigs cfg = OrganelleConfigs::defaults();
    OrganelleMotionState& state = updateOrganelleMotionState(cellKey, time);

    float p0 = organelleSeed01(cellKey, 80) * 2.0f * M_PI;
    float p1 = organelleSeed01(cellKey, 81) * 2.0f * M_PI;
    float p2 = organelleSeed01(cellKey, 82) * 2.0f * M_PI;

    // Nucleus and ER stay near the center, while Golgi and mitochondria roam.
    cfg.nucleus.position = {0.0f, 0.006f * sinf(time * 0.45f + p0), 0.0f};
    cfg.nucleus.rotation.y += 0.05f * sinf(time * 0.18f + p1);

    cfg.smoothER.position = {
        0.26f + 0.03f * sinf(time * 0.42f + p1),
        0.02f * cosf(time * 0.31f + p0),
        0.11f + 0.03f * cosf(time * 0.37f + p2)
    };
    cfg.smoothER.rotation.y += 0.12f * sinf(time * 0.24f + p2);

    cfg.roughER.position = {
        0.02f * sinf(time * 0.28f + p2),
        0.016f * cosf(time * 0.35f + p1),
        0.24f + 0.025f * sinf(time * 0.33f + p0)
    };
    cfg.roughER.rotation.y += 0.10f * sinf(time * 0.21f + p0);

    cfg.golgi.position = state.golgiPos;
    cfg.golgi.rotation.y += time * 0.34f + p1;
    cfg.golgi.rotation.z += 0.18f * sinf(time * 0.52f + p2);

    for (int m = 0; m < 3; m++) {
        cfg.mito[m].position = state.mitoPos[m];
        cfg.mito[m].rotation.x += time * (0.55f + 0.06f * (float)m) + organelleSeed01(cellKey, 90 + m) * 2.0f * M_PI;
        cfg.mito[m].rotation.y += time * (0.42f + 0.08f * (float)m) + organelleSeed01(cellKey, 95 + m) * 2.0f * M_PI;
        cfg.mito[m].rotation.z += 0.22f * sinf(time * (0.70f + 0.05f * (float)m) + organelleSeed01(cellKey, 100 + m) * 2.0f * M_PI);
    }

    return cfg;
}

// ── Organelle model matrix ──────────────────────────────────────────────
static simd_float4x4 organelleModelMatrix(
    simd_float3 cellPos, float cellSize,
    const OrganelleConfig& cfg, float timeOffset
) {
    float t = timeOffset;
    // Simple transform — GPU shader handles membrane containment.
    // No CPU clamping needed.
    simd_float3 jitter = {sinf(t*0.3f)*0.02f, cosf(t*0.5f)*0.015f, sinf(t*0.4f)*0.02f};
    simd_float3 localPos = cfg.position + jitter;

    simd_float4x4 m = mat4_translation(cellPos);
    m = simd_mul(m, mat4_scale(cellSize));
    m = simd_mul(m, mat4_translation(localPos));
    m = simd_mul(m, mat4_rotateX(cfg.rotation.x));
    m = simd_mul(m, mat4_rotateY(cfg.rotation.y));
    m = simd_mul(m, mat4_rotateZ(cfg.rotation.z));
    m = simd_mul(m, mat4_scale(cfg.scale));
    return m;
}

// ── Draw one GLB organelle ──────────────────────────────────────────────
static void drawGLBOrganelle(
    id<MTLRenderCommandEncoder> enc, const OrganelleConfig& cfg,
    simd_float3 cellPos, float cellSize, const Uniforms& baseUni, float toff
) {
    auto it = gOrganelleMeshes.find(cfg.filename);
    if (it == gOrganelleMeshes.end() || !it->second.valid) return;
    const LoadedMesh& mesh = it->second;

    OrgUniforms ou;
    ou.viewProjection = baseUni.viewProjection;
    ou.model = organelleModelMatrix(cellPos, cellSize, cfg, toff);
    ou.cameraPos = baseUni.cameraPos;
    ou.time = baseUni.time;
    ou.lightDir = baseUni.lightDir;
    ou.pad0 = 0;
    ou.baseColor = cfg.color;
    ou.emissiveIntensity = cfg.emissiveIntensity;
    ou.emissiveColor = cfg.emissive;
    ou.pad1 = 0;
    ou.cellCenter = cellPos;
    ou.cellRadius = cellSize;
    // Only apply furrow deformation to the primary cell currently in
    // visual mitosis. Background cells use the sphere clamp only.
    ou.furrowDepth = gMitosis.active ? gMitosis.furrowDepth : 0.0f;
    ou.glowBoost = 0.0f;     // vanilla pass — no pierce-through glow
    ou.pad3 = ou.pad4 = 0.0f;

    [enc setVertexBuffer:mesh.vertexBuffer offset:0 atIndex:0];
    [enc setVertexBytes:&ou length:sizeof(OrgUniforms) atIndex:1];
    [enc setFragmentBytes:&ou length:sizeof(OrgUniforms) atIndex:1];
    [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                    indexCount:mesh.indexCount indexType:MTLIndexTypeUInt32
                   indexBuffer:mesh.indexBuffer indexBufferOffset:0];
}

// Selected-cell pierce-through glow pass: re-render an organelle with
// glowBoost > 0 and depth-test always-pass. The shader emits a strong
// sinusoidal pulse using baseColor; renders on TOP of any other geometry
// so the user sees the organelle's outline through cell membrane,
// fluid, etc. Use only for the focused / selected cell.
static void drawOrganelleGlowPass(id<MTLRenderCommandEncoder> enc,
                                  const Uniforms& baseUni,
                                  simd_float3 cellPos, float cellSize,
                                  const OrganelleConfig& cfg,
                                  float toff,
                                  float glowBoost) {
    auto it = gOrganelleMeshes.find(cfg.filename);
    if (it == gOrganelleMeshes.end() || !it->second.valid) return;
    const LoadedMesh& mesh = it->second;

    OrgUniforms ou;
    ou.viewProjection = baseUni.viewProjection;
    ou.model = organelleModelMatrix(cellPos, cellSize, cfg, toff);
    ou.cameraPos = baseUni.cameraPos;
    ou.time = baseUni.time;
    ou.lightDir = baseUni.lightDir;
    ou.pad0 = 0;
    ou.baseColor = cfg.color;
    ou.emissiveIntensity = cfg.emissiveIntensity;
    ou.emissiveColor = cfg.emissive;
    ou.pad1 = 0;
    ou.cellCenter = cellPos;
    ou.cellRadius = cellSize;
    ou.furrowDepth = gMitosis.active ? gMitosis.furrowDepth : 0.0f;
    ou.glowBoost = glowBoost;
    ou.pad3 = ou.pad4 = 0.0f;

    [enc setVertexBuffer:mesh.vertexBuffer offset:0 atIndex:0];
    [enc setVertexBytes:&ou length:sizeof(OrgUniforms) atIndex:1];
    [enc setFragmentBytes:&ou length:sizeof(OrgUniforms) atIndex:1];
    [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                    indexCount:mesh.indexCount indexType:MTLIndexTypeUInt32
                   indexBuffer:mesh.indexBuffer indexBufferOffset:0];
}

// ── Sync simulation → rendering instances ───────────────────────────────
static void syncCellInstances() {
    int sourceCount = (int)gSim.cells.size();
    gCellInstances.clear();
    gRenderCellSourceIndices.clear();
    gCellInstances.reserve(sourceCount + 1);
    gRenderCellSourceIndices.reserve(sourceCount + 1);
    for (int i = 0; i < sourceCount; i++) {
        if (!shouldRenderCellShellSourceIndex(i)) continue;
        auto& c = gSim.cells[i];
        CellInstance inst{};
        inst.position = c.position;
        inst.radius = c.radius * c.size;
        inst.phase = (float)c.phase;
        inst.lodLevel = 2;
        // Only the primary detailed-mitosis cell should show a cleavage furrow.
        // Applying the global mitosis furrow to every colony cell makes the
        // whole scene appear to reload/deform in sync, which reads as fake.
        inst.furrowDepth = (gMitosis.active && i == 0) ? gMitosis.furrowDepth : 0.0f;

        // Phase + fate dependent colors
        float glow = 1.0f;
        switch (c.fate) {
            case SIM_FATE_PROLIF:
                inst.color = {0.0f, 0.13f, 0.20f, 0.22f}; glow = 2.5f; break;
            case SIM_FATE_QUIESCENT:
                inst.color = {0.0f, 0.06f, 0.13f, 0.12f}; glow = 0.7f; break;
            case SIM_FATE_APOPTOTIC:
                inst.color = {0.15f, 0.02f, 0.0f, 0.28f}; glow = 1.8f; break;
            default: // UNDETERMINED — color by phase
                switch (c.phase) {
                    case 0: inst.color = {0.0f,0.08f,0.15f,0.17f}; glow=1.0f; break;
                    case 1: inst.color = {0.0f,0.10f,0.18f,0.20f}; glow=1.8f; break;
                    case 2: inst.color = {0.06f,0.04f,0.0f,0.17f}; glow=1.5f; break;
                    case 3: inst.color = {0.15f,0.02f,0.0f,0.25f}; glow=2.0f; break;
                }
        }
        // Override for necrotic cells (orange-red glow, swollen)
        if (c.necrotic) {
            inst.color = {0.25f, 0.04f, 0.0f, 0.35f};
            glow = 1.2f;
        }
        // Multi-threshold apoptosis visual — progress drives
        // shrinkage, PS-flip warm tint, fresnel ghost fade. No shader
        // changes required: we bake the effect into the existing
        // CellInstance color / radius / glowIntensity fields.
        float apoP = c.apo.state.apoptosis_progress;
        if (apoP > 0.01f) {
            // Shrinkage — up to SHRINK_FRAC_AT_COMPLETE (35%) at progress=1.
            inst.radius *= (1.0f - Apoptosis::SHRINK_FRAC_AT_COMPLETE * apoP);
            // PS-flip warm tint pulled from the engine (0..1).
            float ps = c.apo.psExposure();
            // Blend toward warm red-brown; alpha drops (cell ghosts).
            inst.color.x = inst.color.x * (1.0f - 0.55f * apoP) + 0.55f * apoP;
            inst.color.y = inst.color.y * (1.0f - 0.55f * apoP) + 0.18f * apoP;
            inst.color.z = inst.color.z * (1.0f - 0.55f * apoP) + 0.10f * apoP;
            // Tint further warmed by PS exposure — MBoC "eat-me" signal.
            inst.color.x = fminf(1.0f, inst.color.x + 0.20f * ps);
            inst.color.y = fminf(1.0f, inst.color.y + 0.08f * ps);
            // Alpha ghost fade from default 0.22 → 0.10 at progress 1.0.
            inst.color.w = fmaxf(0.10f, inst.color.w * (1.0f - 0.55f * apoP));
            // Glow drops as mitos depolarise — no ATP left to fluoresce.
            glow *= (1.0f - 0.60f * apoP);
            // Caspase-driven blebbing feeds furrow depth so the existing
            // shader gives it a transient pinch where no mitosis is.
            // Tie this to an oscillating phase so bleb cycles at ~1 Hz.
            float blebOsc = 0.5f + 0.5f * sinf((float)glfwGetTime() * 6.28f
                                                + (float)c.cellUid * 0.37f);
            inst.furrowDepth = fmaxf(inst.furrowDepth,
                                     Apoptosis::BLEB_AMPLITUDE_FRAC
                                     * apoP * blebOsc);
        }
        // NB: previously the renderer blended a warm-brown tint onto cells
        // whose postDivisionRecovery was non-zero. That tint had no biological
        // referent — daughter cells don't turn brown for 6 seconds after
        // cytokinesis. The timer is still used by updatePhysics to reduce
        // motility for a short period (cortical reorganisation), but it is
        // now invisible to the renderer.
        inst.glowIntensity = glow;
        gCellInstances.push_back(inst);
        gRenderCellSourceIndices.push_back(i);
    }

    int n = (int)gCellInstances.size();

    // ── Cell division rendering ────────────────────────────────────
    // During mitosis (prophase→cytokinesis): SINGLE cell with furrow deformation
    // At MITO_COMPLETE: split into TWO daughters with REAL physics (Hertz repulsion)

    // ── Cell division rendering (Option 1: one source of truth) ──
    // Finalization now happens atomically at the instant cytokinesis
    // ends (CentralDogma.h postDivisionComplete()), which pushes daughter B
    // into gSim.cells and updates daughter A in place. Both cells are then
    // rendered by the normal syncCellInstances() path — no "preview" injection,
    // no separate physics loop, no visual swap.
    //
    // During MITO_PROPHASE..MITO_CYTOKINESIS the cell renders as a single
    // sphere with a furrow deformation (inst.furrowDepth). At the moment
    // cytokinesis completes, finalization replaces that single cell with
    // two real cells in gSim.cells. The next frame's syncCellInstances
    // draws both from the same pipeline with the same colors.
    if (!gMitosis.active || gMitosis.phase != MITO_COMPLETE) {
        sDaughtersInitialized = false;
    }

    // Resize buffer if needed
    size_t bufSize = n * sizeof(CellInstance);
    if (!gCellBuffer || gCellBuffer.length < bufSize) {
        gCellBuffer = [gCtx.device() newBufferWithLength:std::max(bufSize, (size_t)256)
                                                options:MTLResourceStorageModeShared];
    }
    if (n > 0) memcpy(gCellBuffer.contents, gCellInstances.data(), bufSize);
}

// ── GLFW Callbacks ──────────────────────────────────────────────────────
static void mouseButtonCB(GLFWwindow* w, int b, int a, int m) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    double x, y; glfwGetCursorPos(w, &x, &y);
    gCamera.onMouseButton(b, a, x, y);
}
static void cursorPosCB(GLFWwindow* w, double x, double y) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    gCamera.onMouseMove(x, y);
}
static void scrollCB(GLFWwindow* w, double xo, double yo) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    if (gFollowCell) {
        // In follow mode, scroll adjusts the follow-zoom multiplier rather than
        // the raw camera distance (which is overwritten every frame by followCell).
        float factor = powf(1.15f, -(float)yo);
        gFollowZoom = fmaxf(0.1f, fminf(20.0f, gFollowZoom * factor));
    } else {
        gCamera.onScroll(-yo * 3.0);
    }
}
static void framebufferCB(GLFWwindow* w, int width, int height) {
    gCamera.onResize(width, height);
    gCtx.metalLayer().drawableSize = CGSizeMake(width, height);
    gCtx.recreateDepthTexture(width, height);
}

// ── Load models ─────────────────────────────────────────────────────────
static void loadModels() {
    struct Spec { const char* f; uint32_t maxTris; };
    // rough ER.glb is 70MB with 500k+ triangles — decimating to 20k saves
    // hundreds of MB of GPU memory with no visible quality loss at cell scale.
    // Rough ER has 3M+ triangles originally. At 20k (step=153) the mesh
    // decimates into visible islands — you can see individual triangles
    // as tiny disconnected pieces. Keeping 120k retains enough triangle
    // density to read as a continuous network (Palade/Porter 1955 JCB)
    // while staying within GPU memory.
    Spec specs[] = {
        {"nucleus.glb",           20000},
        {"smooth ER.glb",         20000},
        {"rough ER.glb",         120000},
        {"golgi apparatus.glb",   20000},
        {"mitochondria.glb",      20000}
    };
    for (auto& s : specs) {
        std::string path = gModelDir + "/" + s.f;
        auto mesh = GLBLoader::load(gCtx.device(), path, s.maxTris);
        if (mesh.valid) gOrganelleMeshes[s.f] = mesh;
    }
}

// ── Upload molecule + protein synthesis data to GPU ─────────────────────
struct GPUAtomInstance {
    simd_float3 position;
    float radius;
    simd_float3 color;
    float pad;
};

struct GPUBondInstance {
    simd_float3 posA;
    float radius;
    simd_float3 posB;
    float pad;
    simd_float3 color;
    float pad2;
};

static float molHash(int seed, int idx) {
    int h = seed * 2654435761 + idx * 340573321;
    return (float)(h & 0xFFFF) / 65535.0f;
}

// 20 µm cell diameter in Å (1 µm = 10 000 Å). HeLa median diameter per
// BioNumbers BNID 100434. The old value (20 000 Å = 2 µm) was a 10× error
// that scaled every PDB-sourced molecule atom 10× too large relative to
// the cell, making ATP look almost as big as a mitochondrion.
static constexpr float kReferenceCellDiameterAngstrom = 200000.0f;

static simd_float3 rotateY(const simd_float3& v, float cosR, float sinR) {
    return {v.x * cosR - v.z * sinR, v.y, v.x * sinR + v.z * cosR};
}

static float chemistryScaleFromRealSize(const MoleculeData& mol, float cellR,
                                        float visibilityBoost,
                                        float minHalfExtentFrac = 0.0f) {
    if (!mol.valid()) return 0.0f;
    float scale = (2.0f * cellR * visibilityBoost) / kReferenceCellDiameterAngstrom;
    if (minHalfExtentFrac > 0.0f) {
        float minScale = (cellR * minHalfExtentFrac) / fmaxf(mol.maxExtent, 1.0f);
        scale = fmaxf(scale, minScale);
    }
    return scale;
}

// Add a colored sphere to the atom list
static void addSphere(std::vector<GPUAtomInstance>& atoms, simd_float3 pos, float r,
                       float cr, float cg, float cb) {
    GPUAtomInstance a;
    a.position = pos; a.radius = r;
    a.color = {cr, cg, cb}; a.pad = 0;
    atoms.push_back(a);
}

// Add a bond/stick between two points
static void addStick(std::vector<GPUBondInstance>& bonds, simd_float3 a, simd_float3 b,
                      float r, float cr, float cg, float cb) {
    GPUBondInstance bi;
    bi.posA = a; bi.posB = b; bi.radius = r;
    bi.color = {cr, cg, cb}; bi.pad = 0; bi.pad2 = 0;
    bonds.push_back(bi);
}

// Forward declaration — defined further below; needed here so the medium
// chemical render path can call it.
static void addMoleculeGeometry(std::vector<GPUAtomInstance>& atoms,
                                std::vector<GPUBondInstance>& bonds,
                                const MoleculeData& mol, simd_float3 pos,
                                float scale, float tintR, float tintG, float tintB,
                                float rotAngle, int atomStride,
                                float minAtomRadius, float tintMix,
                                float atomRadiusMul, float bondRadiusMul,
                                float unfoldAmount, float repulsionAmount,
                                int seed);

// ── Culture medium chemicals (volumetric floating particles) ────────────
// A persistent swarm of dim, real-molecule renderings that fill the dish
// volume above the cells. Each particle is the actual molecule (glucose
// ring, water tri-atom, calcium ion bead) at the same scale used by the
// cell-interior renderer — so when a particle binds a transporter and
// streams into the cell, its visual size never jumps.
//
// State machine: FREE particles drift with a slow curl-noise current
// field; when they reach a cell with a matching transporter they go
// ATTRACTED → BINDING (membrane pulse) → TRANSPORT (fades into the
// cell) → DESPAWN (teleport to dish edge, become FREE again). This is
// purely visual — the chemistry is already balanced by MediumField's
// stoichiometric exchange called from updateMetabolism.
//
// Density-driven layering (heavy sugars sink, light gases rise) gives
// the fluid its 3-D character; a curl-noise stream function drives the
// horizontal current so particles travel along curved paths.
//
// References:
//   Berg & Purcell 1977 Biophys J — diffusion-limited binding kinetics
//   Stickle 2006 J Biotechnol      — convection in static dish culture
//   Verkman 2008 Annu Rev Med      — aquaporin water flux

enum MediumChemState : uint8_t {
    MC_FREE       = 0,
    MC_ATTRACTED  = 1,
    MC_BINDING    = 2,
    MC_TRANSPORT  = 3,
    MC_DESPAWN    = 4
};

struct MediumChemical {
    simd_float3 home;            // baseline anchor for vertical density layer
    simd_float3 position;        // current world position (advected)
    int   species;               // MediumSpecies enum
    int   targetCellIdx;         // -1 if FREE; else index in gSim.cells
    MediumChemState state;
    float stateTimer;            // bio-seconds in current state
    float phase;                 // per-particle wobble phase
    simd_float3 bindPoint;       // captured membrane point during binding
};

static std::vector<MediumChemical> gMediumChemicals;

// Specification for one species worth of medium particles.
// `sdfId` looks up the cached MoleculeData (MoleculeCache::get).
// nullptr → ad-hoc geometry built at init (currently O2, CO2).
//
// SIZE POLICY: rendered radius derives from the same formula the cell
// interior uses: `radiusUm × UM_TO_WU × visBoost`, where `visBoost` is
// shared with the interior pipeline (`MOL_VIS_BOOST_SMALL` for metabolites,
// `MOL_VIS_BOOST_MACROMOL` for proteins). This guarantees a glucose
// molecule outside the cell renders at the SAME size as glucose inside —
// no awkward jump on uptake. Real molecular sizes (~1 nm in a 10 µm cell)
// produce sub-pixel atoms at typical camera distance; that is the cost
// of physical accuracy and the user has accepted it.
struct MediumChemSpec {
    int          species;
    const char*  sdfId;          // cache key, e.g. "glucose", "water"
    float        radiusUm;       // physical molecular radius (MoleculeRadiusUm)
    float        visBoost;       // SAME boost the cell interior uses
    float        tintR, tintG, tintB;
    float        densityRel;     // 0=light/top, 1=heavy/bottom
    int          count;
};

// Colours match cell-interior spawn data EXACTLY. Particle COUNTS are
// proportional to the real DMEM HG concentrations so the dish reflects
// the actual molarity ratios users see in published media datasheets:
//   glucose 25 mM  >> AA 5 mM  >> glutamine 4 mM  >> calcium 1.8 mM
//   >> CO2 1.2 mM  >> pyruvate 1 mM  >> lactate ~1 mM  >> O2 0.20 mM
// A small floor (~50) is applied to species below 1 mM so they remain
// visible in the swarm. Water is rendered separately at fixed count
// (it dominates real composition at 55 000 mM, but a 5000-particle
// water bath would drown out everything else visually).
// Sources: Sigma D5796 / AAT Bioquest DMEM-HG datasheet, Verkman 2008.
static const MediumChemSpec MEDIUM_SPECS[] = {
    // species       sdfId         radiusUm                                  visBoost                  R     G     B     dens count
    { MS_O2,         nullptr,      0.00015f,                                 MOL_VIS_BOOST_SMALL,      0.55f,0.85f,1.00f, 0.10f,   60 },  // 0.20 mM
    { MS_CO2,        nullptr,      0.00018f,                                 MOL_VIS_BOOST_SMALL,      0.55f,0.55f,0.55f, 0.30f,  120 },  // 1.20 mM
    { MS_GLUCOSE,    "glucose",    MoleculeRadiusUm::GLUCOSE,                MOL_VIS_BOOST_SMALL,      0.90f,0.70f,0.10f, 0.95f, 2500 },  // 25 mM (dominant)
    { MS_GLUTAMINE,  "glycine",    MoleculeRadiusUm::AMINO_ACID,             MOL_VIS_BOOST_SMALL,      0.85f,0.35f,0.35f, 0.55f,  400 },  // 4 mM
    { MS_PYRUVATE,   "glycine",    MoleculeRadiusUm::AMINO_ACID,             MOL_VIS_BOOST_SMALL,      0.85f,0.35f,0.35f, 0.45f,  100 },  // 1 mM
    { MS_LACTATE,    "glycine",    MoleculeRadiusUm::AMINO_ACID,             MOL_VIS_BOOST_SMALL,      0.85f,0.35f,0.35f, 0.50f,  100 },  // ~1 mM (rises with Warburg)
    { MS_AA_POOL,    "glycine",    MoleculeRadiusUm::AMINO_ACID,             MOL_VIS_BOOST_SMALL,      0.85f,0.35f,0.35f, 0.40f,  500 },  // 5 mM essential AA pool
    { MS_GROWTH_F,   "camp",       MoleculeRadiusUm::CHAPERONE,              MOL_VIS_BOOST_MACROMOL,   0.50f,0.90f,0.50f, 0.20f,   50 },  // ~50 ng/mL (FBS)
    { MS_CALCIUM,    "calcium_ion",MoleculeRadiusUm::CALCIUM_ION,            MOL_VIS_BOOST_SMALL,      0.30f,0.70f,0.70f, 0.50f,  180 },  // 1.8 mM
    { MS_WATER,      "water",      MoleculeRadiusUm::WATER,                  MOL_VIS_BOOST_SMALL,      0.80f,0.85f,0.90f, 0.50f, 1000 },  // bath (true ~55 000 mM)
};
static const int MEDIUM_SPEC_COUNT = sizeof(MEDIUM_SPECS) / sizeof(MEDIUM_SPECS[0]);

// Total fluid height above floor (wu). 10 wu ≈ 50 µm — about 10× the
// cell radius, so the dish reads as a shallow film of medium covering
// the cells, not a deep aquarium that hides them.
static constexpr float MEDIUM_FLUID_HEIGHT = 10.0f;

// ── Fluid volume (procedural cylinder mesh + dedicated shader) ───────
// One CONTINUOUS cylinder mesh that exactly matches the dish footprint.
// Top ring is densely tessellated as a disk so the FluidRender vertex
// shader can displace each top vertex with sum-of-sines waves — yielding
// an animated ocean-like surface. Side walls and bottom cap stay flat.
//
// Mesh layout (vertex.position encoding):
//    .xz  = unit-disk position (radius 1)
//    .y   = ring tag: 0.0 = bottom, 1.0 = top (used by vertex shader to
//           place the vertex between floorY and floorY + fluidHeight)
//    .uv.y = role tag: 0=side, 1=top-cap, 2=bottom-cap
//
// The shader applies the dish radius and fluid-height uniforms, plus
// time-driven wave displacement only on the top cap.
struct FluidVertex {
    simd_float3 position;
    simd_float3 normal;
    simd_float2 uv;
};

struct FluidUniforms {
    simd_float4x4 viewProjection;
    simd_float3   cameraPos;     float time;
    simd_float3   lightDir;      float floorY;
    simd_float3   waterColor;    float waterAlpha;
    float radius;
    float fluidHeight;
    float waveAmp;
    float waveSpeed;
};

static id<MTLBuffer> gFluidVertexBuffer = nil;
static id<MTLBuffer> gFluidIndexBuffer  = nil;
static int           gFluidIndexCount   = 0;

static void buildFluidMesh() {
    // Cylinder with N radial segments around the perimeter, plus a
    // tessellated top cap (TR rings × N segments) for high-resolution
    // waves. Bottom cap is a simple fan (waves don't need it).
    const int N  = 96;     // perimeter segments
    const int TR = 16;     // top-cap rings (axis → rim)

    std::vector<FluidVertex> verts;
    std::vector<uint32_t>    idx;

    auto pushVert = [&](float xUnit, float zUnit, float ringTag,
                        simd_float3 normal, float radial, float role) {
        FluidVertex v;
        v.position = {xUnit, ringTag, zUnit};
        v.normal   = normal;
        v.uv       = {radial, role};
        verts.push_back(v);
    };

    // ── Side walls: 2 rings × N verts ─────────────────────────────────
    int sideBottomBase = (int)verts.size();
    for (int i = 0; i < N; i++) {
        float a = 2.0f * (float)M_PI * (float)i / (float)N;
        float xu = cosf(a), zu = sinf(a);
        simd_float3 normal = {xu, 0.0f, zu};
        pushVert(xu, zu, 0.0f, normal, 1.0f, 0.0f);  // bottom side
    }
    int sideTopBase = (int)verts.size();
    for (int i = 0; i < N; i++) {
        float a = 2.0f * (float)M_PI * (float)i / (float)N;
        float xu = cosf(a), zu = sinf(a);
        simd_float3 normal = {xu, 0.0f, zu};
        pushVert(xu, zu, 1.0f, normal, 1.0f, 0.0f);  // top side
    }
    // Side quad indices
    for (int i = 0; i < N; i++) {
        int i2 = (i + 1) % N;
        uint32_t a = sideBottomBase + i;
        uint32_t b = sideBottomBase + i2;
        uint32_t c = sideTopBase + i2;
        uint32_t d = sideTopBase + i;
        idx.push_back(a); idx.push_back(b); idx.push_back(c);
        idx.push_back(a); idx.push_back(c); idx.push_back(d);
    }

    // ── Top cap: tessellated disk (rings × N segments) ────────────────
    // Centre vertex at axis.
    int topCenterIdx = (int)verts.size();
    pushVert(0.0f, 0.0f, 1.0f, simd_float3{0,1,0}, 0.0f, 1.0f);
    int topRingBase[TR + 1];
    topRingBase[0] = topCenterIdx;
    for (int r = 1; r <= TR; r++) {
        float radial = (float)r / (float)TR;
        topRingBase[r] = (int)verts.size();
        for (int i = 0; i < N; i++) {
            float a = 2.0f * (float)M_PI * (float)i / (float)N;
            float xu = cosf(a) * radial, zu = sinf(a) * radial;
            pushVert(xu, zu, 1.0f, simd_float3{0, 1, 0}, radial, 1.0f);
        }
    }
    // Inner ring fan: centre → ring 1
    for (int i = 0; i < N; i++) {
        int i2 = (i + 1) % N;
        idx.push_back(topCenterIdx);
        idx.push_back(topRingBase[1] + i);
        idx.push_back(topRingBase[1] + i2);
    }
    // Outer rings: quads between consecutive rings
    for (int r = 1; r < TR; r++) {
        for (int i = 0; i < N; i++) {
            int i2 = (i + 1) % N;
            uint32_t a = topRingBase[r] + i;
            uint32_t b = topRingBase[r] + i2;
            uint32_t c = topRingBase[r + 1] + i2;
            uint32_t d = topRingBase[r + 1] + i;
            idx.push_back(a); idx.push_back(b); idx.push_back(c);
            idx.push_back(a); idx.push_back(c); idx.push_back(d);
        }
    }

    // ── Bottom cap: simple fan ────────────────────────────────────────
    int bottomCenterIdx = (int)verts.size();
    pushVert(0.0f, 0.0f, 0.0f, simd_float3{0, -1, 0}, 0.0f, 2.0f);
    int bottomRingBase = (int)verts.size();
    for (int i = 0; i < N; i++) {
        float a = 2.0f * (float)M_PI * (float)i / (float)N;
        float xu = cosf(a), zu = sinf(a);
        pushVert(xu, zu, 0.0f, simd_float3{0, -1, 0}, 1.0f, 2.0f);
    }
    for (int i = 0; i < N; i++) {
        int i2 = (i + 1) % N;
        idx.push_back(bottomCenterIdx);
        idx.push_back(bottomRingBase + i2);
        idx.push_back(bottomRingBase + i);
    }

    gFluidVertexBuffer = [gCtx.device() newBufferWithBytes:verts.data()
                          length:verts.size() * sizeof(FluidVertex)
                          options:MTLResourceStorageModeShared];
    gFluidIndexBuffer = [gCtx.device() newBufferWithBytes:idx.data()
                          length:idx.size() * sizeof(uint32_t)
                          options:MTLResourceStorageModeShared];
    gFluidIndexCount = (int)idx.size();
}

static void initFluidHaze() {
    if (!gFluidVertexBuffer) buildFluidMesh();
}

// (Old per-frame haze function removed — fluid is now drawn through its
//  own dedicated pipeline below in the main render pass.)

// Ad-hoc MoleculeData for diatomics not in the SDF cache.
static MoleculeData gO2Mol;
static MoleculeData gCO2Mol;

static void buildAdHocMolecules() {
    // O2 — two oxygens 1.21 Å apart.
    gO2Mol = MoleculeData{};
    Atom o1; o1.position = {-0.605f, 0.0f, 0.0f}; o1.element = 8;
    o1.radius = 1.52f; o1.color[0] = 0.95f; o1.color[1] = 0.35f; o1.color[2] = 0.30f;
    Atom o2 = o1; o2.position = {0.605f, 0.0f, 0.0f};
    gO2Mol.atoms.push_back(o1); gO2Mol.atoms.push_back(o2);
    Bond b; b.atomA = 0; b.atomB = 1; b.order = 2;
    gO2Mol.bonds.push_back(b);
    gO2Mol.computeBounds();

    // CO2 — O=C=O linear, C-O bond 1.16 Å.
    gCO2Mol = MoleculeData{};
    Atom c; c.position = {0.0f, 0.0f, 0.0f}; c.element = 6;
    c.radius = 1.70f; c.color[0] = 0.30f; c.color[1] = 0.30f; c.color[2] = 0.30f;
    Atom oa = o1; oa.position = {-1.16f, 0.0f, 0.0f};
    Atom ob = o1; ob.position = { 1.16f, 0.0f, 0.0f};
    gCO2Mol.atoms.push_back(c);
    gCO2Mol.atoms.push_back(oa);
    gCO2Mol.atoms.push_back(ob);
    Bond ba; ba.atomA = 0; ba.atomB = 1; ba.order = 2;
    Bond bb; bb.atomA = 0; bb.atomB = 2; bb.order = 2;
    gCO2Mol.bonds.push_back(ba); gCO2Mol.bonds.push_back(bb);
    gCO2Mol.computeBounds();
}

// Look up the MoleculeData for a spec — cache for SDF, ad-hoc for others.
static const MoleculeData* mediumChemMol(const MediumChemSpec& spec) {
    if (spec.sdfId) return gMolCache.get(spec.sdfId);
    if (spec.species == MS_O2)  return &gO2Mol;
    if (spec.species == MS_CO2) return &gCO2Mol;
    return nullptr;
}

// Dynamic drug-particle specs — one MediumChemSpec per applied drug so
// drug molecules render through the same ball-and-stick pipeline that
// handles ATP/glucose. `species` values ≥ 1000 route here instead of
// the fixed MEDIUM_SPECS table.  Populated by applyDrugVisuals() when
// Simulation::applyDrug is called.
static std::vector<MediumChemSpec> gDrugSpecs;
static constexpr int DRUG_SPECIES_BASE = 1000;

// Find the MediumChemSpec for a given species id. Searches the fixed
// metabolite table first, then the dynamic drug table.
static const MediumChemSpec* findMediumSpec(int species) {
    if (species >= DRUG_SPECIES_BASE) {
        int idx = species - DRUG_SPECIES_BASE;
        if (idx >= 0 && idx < (int)gDrugSpecs.size()) return &gDrugSpecs[idx];
        return nullptr;
    }
    for (int i = 0; i < MEDIUM_SPEC_COUNT; i++) {
        if (MEDIUM_SPECS[i].species == species) return &MEDIUM_SPECS[i];
    }
    return nullptr;
}

// Apply a drug from the gBioagents registry AND spawn visible particles
// in the medium swarm so the user watches the drug drift through the
// dish, bind cells' membranes, and stream inside. Reuses the existing
// MediumChemical state machine — no new animation code.
static void applyDrugVisuals(const std::string& drugId, float conc_uM) {
    // Look up the drug in the registry (loaded from data/bioagents/drugs.csv).
    const ChemicalEntity* drug = gBioagents.get(drugId);
    if (!drug) {
        printf("[DrugLab] unknown drug '%s'\n", drugId.c_str());
        return;
    }
    // Call the backend to update binding + affinity matrix.
    gSim.applyDrug(drugId, conc_uM);

    // Register a MediumChemSpec for this drug so rendering can find it.
    // One spec per drug id; subsequent applies just update the count.
    int drugSpecIdx = -1;
    for (int i = 0; i < (int)gDrugSpecs.size(); i++) {
        if (gDrugSpecs[i].sdfId && drugId == gDrugSpecs[i].sdfId) {
            drugSpecIdx = i; break;
        }
    }
    if (drugSpecIdx < 0) {
        MediumChemSpec ds{};
        ds.species   = DRUG_SPECIES_BASE + (int)gDrugSpecs.size();
        ds.sdfId     = strdup(drugId.c_str());  // held for sim lifetime
        ds.radiusUm  = MoleculeRadiusUm::GLUCOSE;
        ds.visBoost  = MOL_VIS_BOOST_SMALL;
        ds.tintR     = 1.00f;  // drugs render bright yellow so they're
        ds.tintG     = 0.85f;  // visually distinct from metabolites
        ds.tintB     = 0.15f;
        ds.densityRel = 0.50f;
        ds.count      = 0;     // count filled below per-application
        gDrugSpecs.push_back(ds);
        drugSpecIdx = (int)gDrugSpecs.size() - 1;
    }

    // Spawn ~120 visible particles for the first 10 µM, scaled linearly.
    // Capped so extreme concentrations don't blow the particle budget.
    int nSpawn = (int)fminf(400.0f, fmaxf(30.0f, 12.0f * conc_uM));
    auto urand = []() { return (float)rand() / (float)RAND_MAX; };
    for (int n = 0; n < nSpawn; n++) {
        float r = sqrtf(urand()) * SCENE_BOUND * 0.95f;
        float theta = urand() * 2.0f * (float)M_PI;
        float yFrac = 0.35f + (urand() - 0.5f) * 0.50f;
        yFrac = fmaxf(0.05f, fminf(0.95f, yFrac));
        MediumChemical mc{};
        mc.home = {r * cosf(theta),
                   FLOOR_Y + yFrac * MEDIUM_FLUID_HEIGHT,
                   r * sinf(theta)};
        mc.position = mc.home;
        mc.species  = DRUG_SPECIES_BASE + drugSpecIdx;
        mc.targetCellIdx = -1;
        mc.state = MC_FREE;
        mc.stateTimer = 0.0f;
        mc.phase = urand() * 6.28318f;
        mc.bindPoint = {0, 0, 0};
        gMediumChemicals.push_back(mc);
    }
    printf("[DrugLab] Spawned %d visible particles of %s at %.2f µM\n",
           nSpawn, drugId.c_str(), conc_uM);
}

// Multi-scale curl-noise XZ-plane current field. Divergence-free, so
// particles drift along closed streamlines (incompressible fluid).
// Three octaves layered: a slow large-scale gyre + medium eddies +
// fine turbulence — gives the swirly motion of an unstirred dish with
// thermal-convection currents (Stickle 2006 J Biotechnol). Magnitude
// scales with timeScale so currents accelerate alongside biology when
// the user speeds up time.
static simd_float2 mediumFlow(float x, float z, float t) {
    auto stream = [&](float xq, float zq) {
        // Octave 1: slow large gyre (period ~80 wu)
        float p = sinf(0.035f * xq + 0.06f * t) * cosf(0.035f * zq + 0.05f * t);
        // Octave 2: medium-scale eddies
        p += 0.55f * sinf(0.085f * xq - 0.10f * t) * cosf(0.090f * zq + 0.08f * t);
        // Octave 3: fine turbulence
        p += 0.25f * sinf(0.20f  * xq + 0.18f * t) * sinf(0.18f  * zq - 0.15f * t);
        return p;
    };
    const float h = 0.40f;
    float dpsi_dx = (stream(x + h, z) - stream(x - h, z)) / (2.0f * h);
    float dpsi_dz = (stream(x, z + h) - stream(x, z - h)) / (2.0f * h);
    // Velocity = (∂ψ/∂z, -∂ψ/∂x) — divergence-free incompressible flow.
    return simd_float2{dpsi_dz * 0.80f, -dpsi_dx * 0.80f};
}

static void initMediumChemicals() {
    buildAdHocMolecules();
    gMediumChemicals.clear();
    int total = 0;
    for (int i = 0; i < MEDIUM_SPEC_COUNT; i++) total += MEDIUM_SPECS[i].count;
    gMediumChemicals.reserve(total);

    auto urand = []() { return (float)rand() / (float)RAND_MAX; };

    for (int i = 0; i < MEDIUM_SPEC_COUNT; i++) {
        const MediumChemSpec& spec = MEDIUM_SPECS[i];
        for (int n = 0; n < spec.count; n++) {
            // Uniform random in disk of radius SCENE_BOUND.
            float r = sqrtf(urand()) * SCENE_BOUND * 0.95f;
            float theta = urand() * 2.0f * (float)M_PI;
            // Density-driven Y: heavy near floor, light near top, with
            // ±20 % spread around the band so the layers blend.
            float yFrac = 1.0f - spec.densityRel;
            yFrac += (urand() - 0.5f) * 0.40f;
            yFrac = fmaxf(0.05f, fminf(0.95f, yFrac));
            MediumChemical mc;
            mc.home = {r * cosf(theta),
                       FLOOR_Y + yFrac * MEDIUM_FLUID_HEIGHT,
                       r * sinf(theta)};
            mc.position = mc.home;
            mc.species = spec.species;
            mc.targetCellIdx = -1;
            mc.state = MC_FREE;
            mc.stateTimer = 0.0f;
            mc.phase = urand() * 6.28318f;
            mc.bindPoint = {0, 0, 0};
            gMediumChemicals.push_back(mc);
        }
    }
}

// Per-frame state-machine tick. `dt` is wall seconds; we convert to
// bio-seconds via bioDt() so binding / transport durations stay calibrated
// to the user's timeScale.
static void updateMediumChemicals(float dt) {
    if (gMediumChemicals.empty()) return;
    auto urand = []() { return (float)rand() / (float)RAND_MAX; };
    float t = (float)glfwGetTime();
    float bio_dt = bioDt(dt, gSim.timeScale);

    for (MediumChemical& mc : gMediumChemicals) {
        const MediumChemSpec* spec = findMediumSpec(mc.species);
        if (!spec) { mc.state = MC_DESPAWN; }

        switch (mc.state) {
        case MC_FREE: {
            // Advect by the curl-noise current; current scales with timeScale.
            simd_float2 v = mediumFlow(mc.position.x, mc.position.z, t);
            float currentScale = dt * (1.0f + 0.2f * gSim.timeScale);
            mc.position.x += v.x * currentScale;
            mc.position.z += v.y * currentScale;
            // Vertical relaxation toward density layer with light wobble.
            mc.position.y += (mc.home.y - mc.position.y) * fminf(1.0f, 0.5f * dt);
            mc.position.y += sinf(t * 0.83f + mc.phase) * 0.05f * dt;
            // Wrap around dish bounds.
            float r2 = mc.position.x * mc.position.x + mc.position.z * mc.position.z;
            float bound = SCENE_BOUND * 0.95f;
            if (r2 > bound * bound) {
                float r = sqrtf(r2);
                mc.position.x *= bound / r;
                mc.position.z *= bound / r;
            }

            // Look for a nearby cell with matching transporter.
            int   bestIdx = -1;
            float bestDist2 = MediumBinding::UPTAKE_RADIUS_WU
                             * MediumBinding::UPTAKE_RADIUS_WU;
            for (int ci = 0; ci < (int)gSim.cells.size(); ci++) {
                const SimCell& cell = gSim.cells[ci];
                if (!cell.alive || cell.necrotic) continue;
                float cellR = cell.radius * cell.size;
                float dx = mc.position.x - cell.position.x;
                float dy = mc.position.y - cell.position.y;
                float dz = mc.position.z - cell.position.z;
                float d2 = dx*dx + dy*dy + dz*dz;
                // Distance to membrane = distance to centre - cellR.
                float dMem = sqrtf(d2) - cellR;
                if (dMem < MediumBinding::UPTAKE_RADIUS_WU && dMem * dMem < bestDist2) {
                    bestDist2 = dMem * dMem;
                    bestIdx = ci;
                }
            }
            if (bestIdx >= 0) {
                mc.targetCellIdx = bestIdx;
                mc.state = MC_ATTRACTED;
                mc.stateTimer = 0.0f;
            }
            break;
        }
        case MC_ATTRACTED: {
            if (mc.targetCellIdx < 0
                || mc.targetCellIdx >= (int)gSim.cells.size()
                || !gSim.cells[mc.targetCellIdx].alive) {
                mc.state = MC_FREE; mc.targetCellIdx = -1; break;
            }
            const SimCell& cell = gSim.cells[mc.targetCellIdx];
            float cellR = cell.radius * cell.size;
            float dx = cell.position.x - mc.position.x;
            float dy = cell.position.y - mc.position.y;
            float dz = cell.position.z - mc.position.z;
            float dist = sqrtf(dx*dx + dy*dy + dz*dz);
            if (dist < 0.001f) dist = 0.001f;
            // Aim point: just outside the membrane along the line from
            // particle to cell centre.
            simd_float3 mem = {
                cell.position.x - dx / dist * cellR,
                cell.position.y - dy / dist * cellR,
                cell.position.z - dz / dist * cellR
            };
            float lerpT = fminf(1.0f, 5.0f * dt);
            mc.position.x += (mem.x - mc.position.x) * lerpT;
            mc.position.y += (mem.y - mc.position.y) * lerpT;
            mc.position.z += (mem.z - mc.position.z) * lerpT;
            float dMem = sqrtf((mem.x - mc.position.x) * (mem.x - mc.position.x)
                             + (mem.y - mc.position.y) * (mem.y - mc.position.y)
                             + (mem.z - mc.position.z) * (mem.z - mc.position.z));
            if (dMem < 0.05f) {
                mc.bindPoint = mc.position;
                mc.state = MC_BINDING;
                mc.stateTimer = 0.0f;
            }
            break;
        }
        case MC_BINDING: {
            mc.stateTimer += bio_dt;
            // Particle locked to membrane during binding pulse.
            if (mc.stateTimer >= MediumBinding::BINDING_BIOSEC) {
                mc.state = MC_TRANSPORT;
                mc.stateTimer = 0.0f;
            }
            break;
        }
        case MC_TRANSPORT: {
            if (mc.targetCellIdx < 0
                || mc.targetCellIdx >= (int)gSim.cells.size()
                || !gSim.cells[mc.targetCellIdx].alive) {
                mc.state = MC_DESPAWN; break;
            }
            const SimCell& cell = gSim.cells[mc.targetCellIdx];
            mc.stateTimer += bio_dt;
            float u = mc.stateTimer / MediumBinding::TRANSPORT_BIOSEC;
            if (u >= 1.0f) { mc.state = MC_DESPAWN; break; }
            // Smoothstep ease-in toward the cell centre interior region.
            float s = u * u * (3.0f - 2.0f * u);
            simd_float3 inside = {
                cell.position.x + (mc.bindPoint.x - cell.position.x) * 0.30f,
                cell.position.y + (mc.bindPoint.y - cell.position.y) * 0.30f,
                cell.position.z + (mc.bindPoint.z - cell.position.z) * 0.30f
            };
            mc.position.x = mc.bindPoint.x + (inside.x - mc.bindPoint.x) * s;
            mc.position.y = mc.bindPoint.y + (inside.y - mc.bindPoint.y) * s;
            mc.position.z = mc.bindPoint.z + (inside.z - mc.bindPoint.z) * s;
            break;
        }
        case MC_DESPAWN: {
            // Respawn at the dish edge as a fresh FREE particle of the
            // same species — keeps the swarm size constant.
            float r = SCENE_BOUND * 0.95f;
            float theta = urand() * 2.0f * (float)M_PI;
            float yFrac = 1.0f - spec->densityRel + (urand() - 0.5f) * 0.40f;
            yFrac = fmaxf(0.05f, fminf(0.95f, yFrac));
            mc.home = {r * cosf(theta), FLOOR_Y + yFrac * MEDIUM_FLUID_HEIGHT, r * sinf(theta)};
            mc.position = mc.home;
            mc.targetCellIdx = -1;
            mc.state = MC_FREE;
            mc.stateTimer = 0.0f;
            break;
        }
        }
    }
}

// Render the medium chemical swarm using the *real* cached molecule
// geometry for each species, at the same scale the cell-interior render
// uses — so when a glucose particle binds and streams into the cell, its
// visual size never jumps. State controls brightness so FREE particles
// read as nearly-transparent fluid; bound / transporting ones briefly
// pop and glow.
static void renderMediumChemicals(std::vector<GPUAtomInstance>& atoms,
                                  std::vector<GPUBondInstance>& bonds,
                                  const Simulation& sim,
                                  float wallTime) {
    if (gMediumChemicals.empty()) return;
    // Per-species concentration ratio for a subtle live-modulation. Drops
    // colour brightness when species is depleted, raises it when produced.
    const float speciesInit[MS_COUNT] = {
        MediumComposition::DMEM_O2_MM,
        MediumComposition::DMEM_CO2_MM,
        MediumComposition::DMEM_GLUCOSE_MM,
        MediumComposition::DMEM_GLUTAMINE_MM,
        MediumComposition::DMEM_PYRUVATE_MM,
        1.0f,                                   // lactate baseline
        MediumComposition::DMEM_AA_POOL_MM,
        MediumComposition::DMEM_GROWTH_FACTOR_NG_PER_ML,
        MediumComposition::DMEM_IONS_MM,
        MediumComposition::DMEM_CALCIUM_MM,
        MediumComposition::DMEM_PH,
        MediumComposition::DMEM_WATER_MM
    };

    for (const MediumChemical& mc : gMediumChemicals) {
        if (mc.state == MC_DESPAWN) continue;
        const MediumChemSpec* spec = findMediumSpec(mc.species);
        if (!spec) continue;
        const MoleculeData* mol = mediumChemMol(*spec);
        if (!mol || !mol->valid()) continue;

        // IDENTICAL formula to the cell-interior render at main.mm:5343:
        //   minHalfExtentFrac = fminf(0.018, (renderRadius * 1.45) / cellR)
        //   molScale = chemistryScaleFromRealSize(mol, cellR, 18.0, minHalfExtentFrac)
        // The per-particle `renderRadius` matches interior spawn radius
        // (molRadiusWu) so the floor of `minHalfExtentFrac` keeps the
        // molecule at the same minimum visual size as it would be inside
        // the cell. Result: a glucose in the medium is byte-identical in
        // size to glucose in the cytoplasm.
        float cellR = CELL_RADIUS_BASE;
        float renderRadius = spec->radiusUm * UM_TO_WU * spec->visBoost;
        float minHalfExtentFrac = fminf(0.018f, (renderRadius * 1.45f)
                                                 / fmaxf(cellR, 0.001f));
        float scale = chemistryScaleFromRealSize(*mol, cellR, 18.0f,
                                                 minHalfExtentFrac);

        // State-driven brightness multiplier. With real molecular sizes
        // (~0.003 wu) we need full colour intensity for the tiny atoms to
        // register as visible specks against the pink fluid backdrop —
        // alpha is still capped at 0.22 by the molecule shader so they
        // remain transparent. Depletion still modulates local intensity.
        float bright = 1.0f;
        switch (mc.state) {
        case MC_FREE: {
            // Modulate by local concentration so depletion is visible.
            float local = (float)sim.nutrients.get(mc.species,
                              mc.position.x, mc.position.z);
            float init = speciesInit[mc.species] > 1e-6f
                         ? speciesInit[mc.species] : 1.0f;
            float ratio = (mc.species == MS_LACTATE)
                ? fminf(1.0f, local / init)
                : fmaxf(0.50f, fminf(1.30f, local / init));
            bright = 1.0f * ratio;
            break;
        }
        case MC_ATTRACTED:
            bright = 1.20f;
            break;
        case MC_BINDING: {
            float pulse = 0.95f + 0.30f * sinf(wallTime * 25.0f + mc.phase);
            bright = pulse;
            break;
        }
        case MC_TRANSPORT: {
            float u = mc.stateTimer / MediumBinding::TRANSPORT_BIOSEC;
            bright = 0.85f - 0.30f * u;          // fades into "interior" alpha
            break;
        }
        default: bright = 0.0f;
        }
        if (bright <= 0.001f) continue;

        float rotAngle = mc.phase + wallTime * 0.10f;
        // Strong sinusoidal yellow shine — each particle pulses at its
        // own phase like a firefly. Pulse 0.6..2.0 of base so the
        // brightest peak is well above 1.0 — molecules clearly glow
        // through the dark medium fluid backdrop.
        float shine = 0.6f + 1.4f * (0.5f + 0.5f * sinf(wallTime * 3.5f + mc.phase * 4.0f));
        float yR = 1.00f * shine;
        float yG = 0.95f * shine;
        float yB = 0.40f * shine;
        // Stack species tint × bright × shine + yellow shine so each
        // particle reads as a glowing dot but keeps its species colour.
        float r = spec->tintR * bright * shine + yR * 0.95f;
        float g = spec->tintG * bright * shine + yG * 0.95f;
        float b = spec->tintB * bright * shine + yB * 0.95f;
        addMoleculeGeometry(atoms, bonds, *mol, mc.position,
                            scale, r, g, b,
                            rotAngle, /*atomStride*/1, /*minAtomRadius*/0.0f,
                            /*tintMix*/0.38f, /*atomRadiusMul*/0.56f,
                            /*bondRadiusMul*/0.20f,
                            /*unfoldAmount*/0.0f, /*repulsionAmount*/0.0f,
                            /*seed*/(int)(mc.phase * 1000.0f));
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Apoptotic bodies — companion object to MediumChemical
//  -----------------------------------------------------------------------
//  A dying cell fragments into 5-15 membrane-bound bodies (real 1-5 µm,
//  CellSim 0.10-0.50 wu radius). Each body:
//    • drifts with the same curl-noise current as the chemical swarm,
//    • decomposes its cytosolic / membrane / receptor pools into the
//      local grid cell via MediumField.exchange(),
//    • shrinks as it loses mass,
//    • 6+ bio-h later develops caspase-cleaved plasma-membrane nanopores,
//      Na⁺ rushes in, water follows, volume rises 20-30 % → osmotic
//      lysis (secondary necrosis; Silva 2010, Chen 2016).
//  Bodies within EFFEROCYTOSIS_RADIUS_WU of a live cell can be engulfed
//  into that cell's lysosomal pool (handled on the SimCell side).
// ══════════════════════════════════════════════════════════════════════════

enum ApoBodyRenderPhase : uint8_t {
    ABODY_DRIFT = 0,  // normal decomposition
    ABODY_SWELL = 1,  // osmotic lysis ramp-up
    ABODY_BURST = 2   // lysis rupture — dump & despawn
};

struct ApoptoticBody {
    simd_float3 position;
    simd_float3 velocity;
    float  radius;                // wu
    float  radius0;               // spawn radius (for shrink ratio)
    float  remainingBiomass;      // cytosolic-equivalent (biomass units)
    float  remainingMembrane;     // biomass-equivalent units
    float  remainingReceptor;     // biomass-equivalent units
    float  initialBiomass;
    float  initialMembrane;
    float  initialReceptor;
    float  ageBioSec;
    float  osmoticSwell;          // 0..0.30
    float  membraneIntegrity;     // 1.0 → 0.0
    uint8_t kind;                 // 0=body, 1=nuclear frag, 2=mito frag
    uint8_t phase;                // ApoBodyRenderPhase
    int    originCellUid;
    float  spin;
    float  densityRel;
    // For DAMP signalling — host cell's position is carried so DAMP
    // events can flag a radius of neighbours even after the host cell
    // is gone from gSim.cells.
    simd_float3 lastBurstPos;     // filled on rupture
    bool   burstConsumed;
};
static std::vector<ApoptoticBody> gApoBodies;

// Pending-DAMP-ring: a transient event emitted at body-burst. Each tick
// we flag any neighbouring live cell inside DAMP_NEIGHBOR_RADIUS_WU with
// a short ROS/stress/fasL bump. Events consumed on the next SimCell
// update call (after dispatch), cleared here at end of each frame.
struct DampEvent {
    simd_float3 position;
    float strength;      // scales trigger bumps (0.5..1.5)
    float ageBioSec;
};
static std::vector<DampEvent> gDampEvents;

// Spawn 5-15 apoptotic bodies from a dying cell. Called by the sim
// tick when cell.apoPhase first crosses into BODIES.
static void spawnApoptoticBodies(SimCell& c) {
    auto urand = []() { return (float)rand() / (float)RAND_MAX; };
    int N = Apoptosis::BODIES_PER_CELL_MIN +
            (int)(urand() * (Apoptosis::BODIES_PER_CELL_MAX
                              - Apoptosis::BODIES_PER_CELL_MIN + 1));
    // Fragmentation dump: move FRAG fractions of each pool into bodies.
    float cytoDumpTotal = c.biomass        * Apoptosis::CYTO_FRAG;
    float memDumpTotal  = c.membraneMass_bm * Apoptosis::MEM_FRAG;
    float recDumpTotal  = c.receptorMass_bm * Apoptosis::REC_FRAG;
    c.biomass        -= cytoDumpTotal;
    c.membraneMass_bm -= memDumpTotal;
    c.receptorMass_bm -= recDumpTotal;
    float perBodyCyto = cytoDumpTotal / (float)N;
    float perBodyMem  = memDumpTotal  / (float)N;
    float perBodyRec  = recDumpTotal  / (float)N;
    // Also snap-release the cytosol that *doesn't* go to bodies (the
    // remaining portion after leak+frag). That's Apoptosis::CYTO_BODY
    // for the body ledger, plus whatever is still in the cell's own
    // biomass goes into the field over the body-decomposition phase.
    // To keep the closed-system invariant exact, we also transfer
    // everything that's left in c.biomass / membrane / receptor into
    // the body pool here (any residual not yet released).
    float residualCyto = c.biomass;
    float residualMem  = c.membraneMass_bm;
    float residualRec  = c.receptorMass_bm;
    c.biomass = 0.0f;
    c.membraneMass_bm = 0.0f;
    c.receptorMass_bm = 0.0f;
    float extraPerBody_cyto = residualCyto / (float)N;
    float extraPerBody_mem  = residualMem  / (float)N;
    float extraPerBody_rec  = residualRec  / (float)N;

    for (int i = 0; i < N; i++) {
        ApoptoticBody b{};
        // Jittered sphere just outside the cell boundary.
        float theta = urand() * (float)M_PI * 2.0f;
        float phi   = acosf(2.0f * urand() - 1.0f);
        float r0    = c.radius * c.size * 1.05f;
        simd_float3 dir = {
            sinf(phi) * cosf(theta),
            cosf(phi),
            sinf(phi) * sinf(theta)
        };
        b.position = {
            c.position.x + dir.x * r0,
            fmaxf(FLOOR_Y + 0.15f, c.position.y + dir.y * r0 * 0.4f),
            c.position.z + dir.z * r0
        };
        // Radial outward burst velocity 0.5-1.0 wu/s.
        float speed = 0.5f + urand() * 0.5f;
        b.velocity = {dir.x * speed, dir.y * speed * 0.2f, dir.z * speed};
        // Size: clamp each body radius into [BODY_RADIUS_MIN, MAX].
        float rfrac = Apoptosis::BODY_RADIUS_MIN_WU +
                      urand() * (Apoptosis::BODY_RADIUS_MAX_WU
                                 - Apoptosis::BODY_RADIUS_MIN_WU);
        b.radius   = rfrac;
        b.radius0  = rfrac;
        b.remainingBiomass  = perBodyCyto + extraPerBody_cyto;
        b.remainingMembrane = perBodyMem  + extraPerBody_mem;
        b.remainingReceptor = perBodyRec  + extraPerBody_rec;
        b.initialBiomass    = b.remainingBiomass;
        b.initialMembrane   = b.remainingMembrane;
        b.initialReceptor   = b.remainingReceptor;
        b.ageBioSec         = 0.0f;
        b.osmoticSwell      = 0.0f;
        b.membraneIntegrity = 1.0f;
        // Kind mixing: ~1 nuclear frag, 1-2 mito frags, rest generic.
        b.kind = (i == 0) ? 1 : ((i <= 2) ? 2 : 0);
        b.phase = ABODY_DRIFT;
        b.originCellUid = c.cellUid;
        b.spin = urand() * 6.28318f;
        b.densityRel = 0.65f + urand() * 0.10f;
        b.lastBurstPos = {0,0,0};
        b.burstConsumed = false;
        gApoBodies.push_back(b);
    }
    // Mass is now fully partitioned into bodies. When updateApoptosis
    // next transitions this cell into BODIES it will call
    // releaseAllMass(c), which is a no-op since c.biomass / membrane /
    // receptor are already zero.
}

// Per-frame update (called next to updateMediumChemicals()). Drift,
// decompose, possibly undergo osmotic lysis. All rates in bio-seconds.
static void updateApoptoticBodies(float dt) {
    if (gApoBodies.empty()) return;
    auto urand = []() { return (float)rand() / (float)RAND_MAX; };
    float t = (float)glfwGetTime();
    float bio_dt = bioDt(dt, gSim.timeScale);

    for (size_t i = 0; i < gApoBodies.size(); ) {
        ApoptoticBody& b = gApoBodies[i];
        b.ageBioSec += bio_dt;
        // Advect by curl-noise + self velocity + sink.
        simd_float2 v = mediumFlow(b.position.x, b.position.z, t);
        float currentScale = dt * (1.0f + 0.2f * gSim.timeScale);
        b.position.x += v.x * currentScale + b.velocity.x * dt;
        b.position.z += v.y * currentScale + b.velocity.z * dt;
        b.position.y += b.velocity.y * dt;
        b.velocity.x *= 0.92f;
        b.velocity.y *= 0.92f;
        b.velocity.z *= 0.92f;
        b.position.y -= Apoptosis::BODY_SINK_WU_PER_BIOSEC
                        * b.densityRel * bio_dt;
        if (b.position.y < FLOOR_Y + b.radius * 0.5f) {
            b.position.y = FLOOR_Y + b.radius * 0.5f;
        }
        // Wrap XZ inside dish.
        float r2 = b.position.x*b.position.x + b.position.z*b.position.z;
        float bound = SCENE_BOUND * 0.95f;
        if (r2 > bound * bound) {
            float rNorm = sqrtf(r2);
            b.position.x *= bound / rNorm;
            b.position.z *= bound / rNorm;
        }

        // Decomposition (phase 0): each pool decays at its own rate;
        // released mass goes to the MediumField via the same partitioning
        // table the cell's releaseCytosol/Membrane/Receptors helpers use,
        // plus a small extracellular-protease deposit so nearby cells
        // "eat better" around dying clusters.
        auto decomposeFraction = [&](float poolInitial,
                                     float &poolRemaining,
                                     float biosec) {
            if (poolRemaining <= 0.0f) return 0.0f;
            float dFrac = fminf(1.0f, bio_dt / biosec);
            float drelease = poolRemaining * dFrac;
            poolRemaining -= drelease;
            return drelease;
        };

        if (b.phase == ABODY_DRIFT) {
            float dcyto = decomposeFraction(b.initialBiomass,
                                            b.remainingBiomass,
                                            Apoptosis::DECOMPOSITION_BIOSEC);
            float dmem  = decomposeFraction(b.initialMembrane,
                                            b.remainingMembrane,
                                            Apoptosis::MEMBRANE_DECOMP_BIOSEC);
            float drec  = decomposeFraction(b.initialReceptor,
                                            b.remainingReceptor,
                                            Apoptosis::RECEPTOR_DECOMP_BIOSEC);
            if (dcyto > 0 || dmem > 0 || drec > 0) {
                float flux[MS_COUNT] = {0};
                flux[MS_AA_POOL]  = dcyto * Apoptosis::REL_AA_PER_BIOMASS
                                  + dmem  * Apoptosis::REL_AA_PER_MEMBRANE
                                  + drec  * Apoptosis::REL_AA_PER_RECEPTOR;
                flux[MS_IONS]     = dcyto * Apoptosis::REL_IONS_PER_BIOMASS;
                flux[MS_CALCIUM]  = dcyto * Apoptosis::REL_CALCIUM_PER_BIOMASS;
                flux[MS_PYRUVATE] = dcyto * Apoptosis::REL_PYRUVATE_PER_BIOMASS;
                flux[MS_LACTATE]  = dcyto * Apoptosis::REL_LACTATE_PER_BIOMASS;
                flux[MS_GLUCOSE]  = dcyto * Apoptosis::REL_GLUCOSE_PER_BIOMASS
                                  + dmem  * Apoptosis::REL_GLUCOSE_PER_MEMBRANE
                                  + drec  * Apoptosis::REL_GLUCOSE_PER_RECEPTOR;
                flux[MS_WATER]    = dcyto * Apoptosis::REL_WATER_PER_BIOMASS;
                gSim.nutrients.exchange(b.position.x, b.position.z, flux, 1.0f);
                gSim.nutrients.depositProtease(b.position.x, b.position.z,
                    (dcyto + dmem) * Apoptosis::EP_RELEASE_PER_LYSED_BIOMASS * 0.5f);
            }
            // Shrink proportional to cube-root of mass.
            float massFrac = (b.remainingBiomass + b.remainingMembrane
                              + b.remainingReceptor)
                            / fmaxf(1e-6f, b.initialBiomass
                                            + b.initialMembrane
                                            + b.initialReceptor);
            b.radius = b.radius0 * cbrtf(fmaxf(0.05f, massFrac));

            // Once old enough, begin osmotic-lysis nanopore formation.
            if (b.ageBioSec > Apoptosis::SECONDARY_NECROSIS_PORE_START_BIOSEC) {
                b.phase = ABODY_SWELL;
            }
        } else if (b.phase == ABODY_SWELL) {
            // Integrity decay cascade.
            float integrityDrop = bio_dt / Apoptosis::PORE_FORM_BIOSEC;
            b.membraneIntegrity = fmaxf(0.0f, b.membraneIntegrity - integrityDrop);
            // Osmotic water influx (real Chen 2016 mechanism, here
            // lumped into a geometric swell coefficient). Internal
            // osmolarity rises as ions leak — model by pushing swell
            // up whenever integrity < 0.5.
            if (b.membraneIntegrity < 0.5f) {
                float dOsmo = 1.0f;   // saturating internal:external mismatch
                b.osmoticSwell += Apoptosis::LP_PORE_FAILED
                                  * dOsmo
                                  * (1.0f - b.membraneIntegrity)
                                  * bio_dt
                                  * (1.0f / Apoptosis::SWELL_COEFFICIENT);
                // Radius: original × (1 + swell)^(1/3) to preserve
                // mass/volume relationship as best a simple body can.
                b.radius = b.radius0 * cbrtf(1.0f + b.osmoticSwell);
            }
            if (b.osmoticSwell >= Apoptosis::LYSIS_RUPTURE_SWELL) {
                b.phase = ABODY_BURST;
                b.lastBurstPos = b.position;
            }
        } else if (b.phase == ABODY_BURST) {
            if (!b.burstConsumed) {
                // One-tick 100 % dump of everything remaining.
                float dcyto = b.remainingBiomass;
                float dmem  = b.remainingMembrane;
                float drec  = b.remainingReceptor;
                b.remainingBiomass = 0;
                b.remainingMembrane = 0;
                b.remainingReceptor = 0;
                float flux[MS_COUNT] = {0};
                flux[MS_AA_POOL]  = dcyto * Apoptosis::REL_AA_PER_BIOMASS
                                  + dmem  * Apoptosis::REL_AA_PER_MEMBRANE
                                  + drec  * Apoptosis::REL_AA_PER_RECEPTOR;
                flux[MS_IONS]     = dcyto * Apoptosis::REL_IONS_PER_BIOMASS;
                flux[MS_CALCIUM]  = dcyto * Apoptosis::REL_CALCIUM_PER_BIOMASS
                                  + Apoptosis::DAMP_CA_BURST_MM;
                flux[MS_PYRUVATE] = dcyto * Apoptosis::REL_PYRUVATE_PER_BIOMASS;
                flux[MS_LACTATE]  = dcyto * Apoptosis::REL_LACTATE_PER_BIOMASS;
                flux[MS_GLUCOSE]  = dcyto * Apoptosis::REL_GLUCOSE_PER_BIOMASS
                                  + dmem  * Apoptosis::REL_GLUCOSE_PER_MEMBRANE
                                  + drec  * Apoptosis::REL_GLUCOSE_PER_RECEPTOR;
                flux[MS_WATER]    = dcyto * Apoptosis::REL_WATER_PER_BIOMASS;
                gSim.nutrients.exchange(b.position.x, b.position.z, flux, 1.0f);
                gSim.nutrients.depositProtease(b.position.x, b.position.z,
                    (dcyto + dmem) * Apoptosis::EP_RELEASE_PER_LYSED_BIOMASS);
                // Emit DAMP event.
                gDampEvents.push_back({b.position, 1.0f, 0.0f});
                b.burstConsumed = true;
                // Keep a brief visual fade by allowing the body to
                // stay in ABODY_BURST for LYSIS_DUMP_BIOSEC before
                // despawning — use ageBioSec offset.
            }
            // Hold for LYSIS_DUMP_BIOSEC then despawn.
            if (b.ageBioSec > Apoptosis::SECONDARY_NECROSIS_PORE_START_BIOSEC
                            + Apoptosis::OSMO_SWELL_BIOSEC
                            + Apoptosis::LYSIS_DUMP_BIOSEC) {
                gApoBodies.erase(gApoBodies.begin() + i);
                continue;
            }
        }

        // Normal despawn when fully decomposed.
        float massFracNow = (b.remainingBiomass + b.remainingMembrane
                              + b.remainingReceptor)
                            / fmaxf(1e-6f, b.initialBiomass
                                            + b.initialMembrane
                                            + b.initialReceptor);
        if (massFracNow < 0.02f && b.phase != ABODY_BURST) {
            gApoBodies.erase(gApoBodies.begin() + i);
            continue;
        }
        i++;
    }

    // Decay DAMP events so they don't accumulate forever.
    for (size_t i = 0; i < gDampEvents.size(); ) {
        gDampEvents[i].ageBioSec += bio_dt;
        if (gDampEvents[i].ageBioSec > 60.0f) {
            gDampEvents.erase(gDampEvents.begin() + i);
        } else i++;
    }
}

// Render apoptotic bodies using the same addStick/addMoleculeGeometry-free
// atom emitter — bodies are simple translucent spheres (represented as a
// single "atom" at their position) with warm tints and a size that
// matches their remainingBiomass. Colour varies by `kind`.
static void renderApoptoticBodies(std::vector<GPUAtomInstance>& atoms) {
    if (gApoBodies.empty()) return;
    float t = (float)glfwGetTime();
    for (const ApoptoticBody& b : gApoBodies) {
        float fade = 1.0f;
        if (b.phase == ABODY_BURST) {
            // Quick fade over the LYSIS_DUMP window.
            float u = (b.ageBioSec - Apoptosis::SECONDARY_NECROSIS_PORE_START_BIOSEC
                                    - Apoptosis::OSMO_SWELL_BIOSEC)
                      / Apoptosis::LYSIS_DUMP_BIOSEC;
            fade = fmaxf(0.0f, 1.0f - u);
        } else {
            fade = fmaxf(0.2f, fminf(1.0f,
                         (b.remainingBiomass + b.remainingMembrane
                          + b.remainingReceptor)
                         / fmaxf(1e-6f, b.initialBiomass
                                          + b.initialMembrane
                                          + b.initialReceptor)));
        }
        // Kind → colour (warm for generic; deep red for nuclear;
        // green-tinged for mito cytochrome-c residue). Alpha is baked
        // into the RGB via fade — the atom shader caps alpha itself.
        float rC, gC, bC;
        switch (b.kind) {
        case 1: // nuclear fragment — deep red
            rC = 0.55f; gC = 0.08f; bC = 0.14f; break;
        case 2: // mito fragment — green cyt-c residue
            rC = 0.25f; gC = 0.55f; bC = 0.22f; break;
        default: // generic body — warm brown
            rC = 0.40f; gC = 0.18f; bC = 0.12f; break;
        }
        // When swelling, add a bright flash via the rim — approximate
        // by boosting saturation proportional to osmoticSwell.
        if (b.phase == ABODY_SWELL) {
            float flash = b.osmoticSwell / Apoptosis::LYSIS_RUPTURE_SWELL;
            rC = fminf(1.0f, rC + 0.35f * flash);
            gC = fminf(1.0f, gC + 0.12f * flash);
        }
        // Bake fade into brightness (~alpha substitute).
        rC *= fade; gC *= fade; bC *= fade;
        // Gentle sine wobble so bodies don't look frozen.
        float wobble = 1.0f + 0.03f * sinf(t * 0.9f + b.spin);
        GPUAtomInstance ai;
        ai.position = b.position;
        ai.radius = b.radius * wobble;
        ai.color = {rC, gC, bC};
        ai.pad = 0.0f;
        atoms.push_back(ai);
    }
}

// Each tick, deliver DAMP bumps to live neighbours within the radius.
// Called from the main sim update after updateApoptoticBodies().
static void applyDampEventsToNeighbors(float dt) {
    if (gDampEvents.empty() || gSim.cells.empty()) return;
    float bio_dt = bioDt(dt, gSim.timeScale);
    float R = Apoptosis::DAMP_NEIGHBOR_RADIUS_WU;
    float R2 = R * R;
    for (DampEvent& e : gDampEvents) {
        for (SimCell& c : gSim.cells) {
            if (!c.alive || c.apoPhase != Apoptosis::ALIVE) continue;
            float dx = c.position.x - e.position.x;
            float dy = c.position.y - e.position.y;
            float dz = c.position.z - e.position.z;
            float d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > R2) continue;
            float falloff = 1.0f - sqrtf(d2) / R;
            float bump = e.strength * falloff * bio_dt / 30.0f; // over 30 bio-s
            c.ROS = fminf(100.0f, c.ROS
                + Apoptosis::DAMP_NEIGHBOR_ROS_BOOST * bump);
            c.stress = fminf(100.0f, c.stress
                + Apoptosis::DAMP_NEIGHBOR_STRESS_BOOST * bump);
            c.fasLExposure = fminf(50.0f, c.fasLExposure
                + Apoptosis::DAMP_NEIGHBOR_FASL_NG_PER_ML * bump);
            c.damageSensedBioSec = fmaxf(c.damageSensedBioSec, 60.0f);
        }
    }
}

// Efferocytosis: live cells engulf nearby apoptotic bodies.
static void updateEfferocytosis(float dt) {
    if (gApoBodies.empty() || gSim.cells.empty()) return;
    float bio_dt = bioDt(dt, gSim.timeScale);
    for (SimCell& c : gSim.cells) {
        if (!c.alive || c.apoPhase != Apoptosis::ALIVE) continue;
        float cellR = c.radius * c.size;
        for (ApoptoticBody& b : gApoBodies) {
            if (b.phase == ABODY_BURST) continue;
            float dx = c.position.x - b.position.x;
            float dy = c.position.y - b.position.y;
            float dz = c.position.z - b.position.z;
            float d = sqrtf(dx*dx + dy*dy + dz*dz) - cellR - b.radius;
            if (d > Apoptosis::EFFEROCYTOSIS_RADIUS_WU) continue;
            if (d < 0) d = 0;
            float rate = Apoptosis::EFFEROCYTOSIS_RATE_PER_BIOSEC
                       * c.adhesionStrength
                       * (1.0f - d / Apoptosis::EFFEROCYTOSIS_RADIUS_WU);
            float dFrac = fminf(1.0f, rate * bio_dt);
            float moveCyto = b.remainingBiomass  * dFrac;
            float moveMem  = b.remainingMembrane * dFrac;
            float moveRec  = b.remainingReceptor * dFrac;
            b.remainingBiomass  -= moveCyto;
            b.remainingMembrane -= moveMem;
            b.remainingReceptor -= moveRec;
            c.lysosomalLoad_cyto += moveCyto;
            c.lysosomalLoad_mem  += moveMem;
            c.lysosomalLoad_rec  += moveRec;
        }
    }
}

// Lysosomal digestion: each tick, drain cell lysosomal pools back into
// cell biomass (RECYCLE_EFFICIENCY fraction); the remainder lost to
// CO2/H2O returning to the field so mass balance still closes.
static void updateLysosomalDigestion(float dt) {
    if (gSim.cells.empty()) return;
    float bio_dt = bioDt(dt, gSim.timeScale);
    for (SimCell& c : gSim.cells) {
        if (!c.alive) continue;
        auto digest = [&](float &pool, float biosec) {
            if (pool <= 0) return 0.0f;
            float frac = fminf(1.0f, bio_dt / biosec);
            float d = pool * frac;
            pool -= d;
            return d;
        };
        float dCyto = digest(c.lysosomalLoad_cyto, Apoptosis::LYSO_DIGEST_CYTO_BIOSEC);
        float dMem  = digest(c.lysosomalLoad_mem,  Apoptosis::LYSO_DIGEST_MEM_BIOSEC);
        float dRec  = digest(c.lysosomalLoad_rec,  Apoptosis::LYSO_DIGEST_REC_BIOSEC);
        float total = dCyto + dMem + dRec;
        if (total <= 0) continue;
        float recycled = total * Apoptosis::RECYCLE_EFFICIENCY;
        float lost     = total - recycled;
        c.biomass += recycled;
        c.ATP     = fminf(100.0f, c.ATP + recycled * Apoptosis::ATP_PER_RECYCLED_BM);
        // Return oxidative loss to field as CO2 + H2O — closes mass balance.
        float flux[MS_COUNT] = {0};
        flux[MS_CO2]   = lost * Apoptosis::CO2_PER_RECYCLED_LOSS;
        flux[MS_WATER] = lost * Apoptosis::H2O_PER_RECYCLED_LOSS;
        gSim.nutrients.exchange(c.position.x, c.position.z, flux, 1.0f);
    }
}

static void addMoleculeGeometry(std::vector<GPUAtomInstance>& atoms,
                                std::vector<GPUBondInstance>& bonds,
                                const MoleculeData& mol, simd_float3 pos,
                                float scale, float tintR, float tintG, float tintB,
                                float rotAngle = 0.0f, int atomStride = 1,
                                float minAtomRadius = 0.0f, float tintMix = 0.45f,
                                float atomRadiusMul = 0.55f, float bondRadiusMul = 0.22f,
                                float unfoldAmount = 0.0f, float repulsionAmount = 0.0f,
                                int seed = 0) {
    if (!mol.valid()) return;
    if (atomStride < 1) atomStride = 1;

    float cosR = cosf(rotAngle), sinR = sinf(rotAngle);
    std::vector<int> atomRemap(mol.atoms.size(), -1);

    for (int i = 0; i < (int)mol.atoms.size(); i += atomStride) {
        const auto& a = mol.atoms[i];
        simd_float3 centered = {
            a.position.x - mol.center.x,
            a.position.y - mol.center.y,
            a.position.z - mol.center.z
        };

        float len = sqrtf(centered.x * centered.x + centered.y * centered.y + centered.z * centered.z);
        float hx = molHash(seed + 11, i * 3 + 0) * 2.0f - 1.0f;
        float hy = molHash(seed + 17, i * 3 + 1) * 2.0f - 1.0f;
        float hz = molHash(seed + 23, i * 3 + 2) * 2.0f - 1.0f;
        float hLen = sqrtf(hx * hx + hy * hy + hz * hz);
        simd_float3 noiseDir = (hLen > 0.0001f)
            ? simd_float3{hx / hLen, hy / hLen, hz / hLen}
            : simd_float3{0.0f, 1.0f, 0.0f};
        simd_float3 outward = (len > 0.0001f)
            ? simd_float3{centered.x / len, centered.y / len, centered.z / len}
            : noiseDir;

        float inflate = mol.maxExtent * unfoldAmount * (0.30f + 0.70f * molHash(seed + 41, i));
        float repel = mol.maxExtent * repulsionAmount * (0.35f + 0.65f * molHash(seed + 67, i));
        simd_float3 local = {
            centered.x + outward.x * inflate + noiseDir.x * repel,
            centered.y + outward.y * inflate + noiseDir.y * repel,
            centered.z + outward.z * inflate + noiseDir.z * repel
        };
        local = {local.x * scale, local.y * scale, local.z * scale};
        local = rotateY(local, cosR, sinR);

        GPUAtomInstance ai;
        ai.position = {pos.x + local.x, pos.y + local.y, pos.z + local.z};
        ai.radius = fmaxf(minAtomRadius, a.radius * scale * atomRadiusMul);
        ai.color = {
            a.color[0] * (1.0f - tintMix) + tintR * tintMix,
            a.color[1] * (1.0f - tintMix) + tintG * tintMix,
            a.color[2] * (1.0f - tintMix) + tintB * tintMix
        };
        ai.pad = 0.0f;
        atomRemap[i] = (int)atoms.size();
        atoms.push_back(ai);
    }

    for (const auto& bond : mol.bonds) {
        if (bond.atomA < 0 || bond.atomB < 0 ||
            bond.atomA >= (int)atomRemap.size() || bond.atomB >= (int)atomRemap.size()) {
            continue;
        }
        int idxA = atomRemap[bond.atomA];
        int idxB = atomRemap[bond.atomB];
        if (idxA < 0 || idxB < 0) continue;
        float radius = fmaxf(minAtomRadius * 0.40f,
                             scale * bondRadiusMul * (0.55f + 0.18f * (float)bond.order));
        addStick(bonds, atoms[idxA].position, atoms[idxB].position, radius,
                 tintR * 0.62f, tintG * 0.62f, tintB * 0.62f);
    }
}

// ── Initialize intracellular particle system ────────────────────────────
static simd_float3 gLastCellPos = {0, 0, 0}; // Track cell movement

// ── Render a real protein backbone at a position ────────────────────────
static void addProteinBackbone(std::vector<GPUAtomInstance>& atoms,
                                std::vector<GPUBondInstance>& bonds,
                                const MoleculeData& mol, simd_float3 pos,
                                float cellR, float tintR, float tintG, float tintB,
                                float rotAngle = 0, int stride = 1,
                                float foldProgress = 1.0f) {
    if (!mol.valid()) return;
    float clampedFold = fmaxf(0.0f, fminf(foldProgress, 1.0f));
    float scale = chemistryScaleFromRealSize(mol, cellR, 2.5f, 0.012f);
    float minAtomRadius = fmaxf(cellR * 0.0018f, scale * 0.8f);
    addMoleculeGeometry(atoms, bonds, mol, pos, scale, tintR, tintG, tintB,
                        rotAngle, stride, minAtomRadius, 0.72f, 0.62f, 0.24f,
                        (1.0f - clampedFold) * 0.78f,
                        (1.0f - clampedFold) * 0.24f,
                        (int)(rotAngle * 100.0f) + (int)mol.atoms.size());
}

// Map particle type to protein PDB filename
static const char* particleToProtein(ParticleType type) {
    switch (type) {
        case PT_TRNA:            return "tRNA-Phe.pdb";
        case PT_RIBOSOME_LARGE:  return "ribosome_80S.pdb";
        case PT_RNA_POLYMERASE:  return "rna_pol_II.pdb";
        case PT_SPLICEOSOME:     return "spliceosome.pdb";
        case PT_CHAPERONE:       return "hsp70_chaperone.pdb";
        case PT_DNA_POLYMERASE:  return "dna_polymerase.pdb";
        // Nup98 C-terminal autoproteolytic domain (PDB 2Q5X). Real NPC
        // contains ~30 different Nups (~500 copies total, 120 nm diameter);
        // this is one representative channel nucleoporin with FG-repeat
        // spatial extent — Hoelz 2011 Annu Rev Biochem.
        case PT_NUCLEAR_PORE:    return "nup98.pdb";
        // 1KX5: canonical nucleosome core particle (147 bp DNA wrapped
        // on histone octamer H2A/H2B/H3/H4). Rendered as Cα trace of the
        // octamer plus DNA phosphate backbone — Luger 1997 Nature.
        case PT_DNA_NODE:        return "nucleosome_1kx5.pdb";
        default: return nullptr;
    }
}

static const char* particleTagToMolecule(MolTag tag) {
    switch (tag) {
        case TAG_ATP:        return "atp";
        case TAG_GLUCOSE:    return "glucose";
        case TAG_NADH:       return "nadh";
        case TAG_AMINO_ACID: return "glycine";
        case TAG_WATER:      return "water";
        case TAG_CALCIUM:    return "calcium_ion";
        default: return nullptr;
    }
}

// ── Initialize particles for a non-primary cell's interior ──────────────
// Same layout as initIntracellularParticles but into a specific CellInterior
static void initCellInterior(CellInterior& ci, float cellSize) {
    float cs = cellSize;
    float cellR = cs;
    float nucR = cs * 0.45f;
    simd_float3 origin = {0, 0, 0};
    ci.phys.init(origin, cellR);

    simd_float3 nucleusPos = origin;
    float helixR = nucR * 0.3f;
    float helixH = nucR * 1.6f;
    int turns = 8;
    int dnaStride = 2;
    int dnaNodes = HBB_LENGTH / dnaStride;

    // DNA double helix
    for (int i = 0; i < dnaNodes; i++) {
        float t = (float)i / dnaNodes;
        float angle = t * turns * 2.0f * M_PI;
        float y = (t - 0.5f) * helixH;
        int bpIdx = i * dnaStride;
        bool ex = (bpIdx < HBB_LENGTH) ? isExon(bpIdx) : false;
        char base = (bpIdx < HBB_LENGTH) ? HBB_SEQUENCE[bpIdx] : 'A';
        BaseColor bc = getBaseColor(base, ex);
        BaseColor bc2 = getBaseColor(complement(base), ex);

        // Particle sphere radii: physical µm × UM_TO_WU × visibility boost.
        // DNA node = 10 nm nucleosome (Alberts Ch 4). See MoleculeRadiusUm.
        const float rDNA = molRadiusWu(MoleculeRadiusUm::DNA_NODE, MOL_VIS_BOOST_MACROMOL);
        simd_float3 pA = {cosf(angle) * helixR, y, sinf(angle) * helixR};
        int idA = ci.phys.addParticle(PT_DNA_NODE, pA, rDNA, bc.r, bc.g, bc.b);
        ci.phys.setHome(idA, pA, 2.0f);
        ci.phys.setConfinement(idA, nucleusPos, nucR);

        simd_float3 pB = {cosf(angle + M_PI) * helixR, y, sinf(angle + M_PI) * helixR};
        int idB = ci.phys.addParticle(PT_DNA_NODE, pB, rDNA, bc2.r, bc2.g, bc2.b);
        ci.phys.setHome(idB, pB, 2.0f);
        ci.phys.setConfinement(idB, nucleusPos, nucR);

        if (i > 0) {
            ci.phys.setTether(idA, idA - 2, cs * 0.025f, 5.0f);
            ci.phys.setTether(idB, idB - 2, cs * 0.025f, 5.0f);
        }
        ci.phys.particles[idA].stateIndex = bpIdx;
        ci.phys.particles[idB].stateIndex = bpIdx;
    }

    // RNA Polymerase
    for (int i = 0; i < 3; i++) {
        float t = (float)i * 0.3f;
        simd_float3 p = {sinf(t * 4) * nucR * 0.2f, (t - 0.5f) * helixH * 0.5f, cosf(t * 4) * nucR * 0.2f};
        int id = ci.phys.addParticle(PT_RNA_POLYMERASE, p,
            molRadiusWu(MoleculeRadiusUm::RNA_POL, MOL_VIS_BOOST_MACROMOL),
            0.6f, 0.2f, 0.8f);
        ci.phys.setConfinement(id, nucleusPos, nucR);
    }

    // DNA Polymerase (small replisome markers)
    for (int i = 0; i < 4; i++) {
        float angle = (float)i / 4.0f * 2.0f * M_PI;
        simd_float3 p = {
            cosf(angle) * nucR * 0.12f,
            ((i < 2) ? -0.04f : 0.04f) * cs,
            sinf(angle) * nucR * 0.12f
        };
        int id = ci.phys.addParticle(PT_DNA_POLYMERASE, p,
            molRadiusWu(MoleculeRadiusUm::DNA_POL, MOL_VIS_BOOST_MACROMOL),
            0.55f, 0.78f, 1.0f);
        ci.phys.setConfinement(id, nucleusPos, nucR);
        ci.phys.setHome(id, p, 0.7f);
    }

    // Spliceosome (2)
    for (int i = 0; i < 2; i++) {
        simd_float3 p = {(i ? 0.15f : -0.15f) * cs, 0.1f * cs, (i ? 0.1f : -0.1f) * cs};
        int id = ci.phys.addParticle(PT_SPLICEOSOME, p,
            molRadiusWu(MoleculeRadiusUm::SPLICEOSOME, MOL_VIS_BOOST_MACROMOL),
            0.8f, 0.6f, 0.2f);
        ci.phys.setConfinement(id, nucleusPos, nucR);
        ci.phys.setHome(id, p, 0.5f);
    }

    // Nuclear pores (6) — 120 nm NPC (Alber 2007 Nature)
    for (int i = 0; i < 6; i++) {
        float a = (float)i / 6.0f * 2.0f * M_PI;
        simd_float3 p = {cosf(a) * nucR * 0.95f, sinf(a * 0.5f) * nucR * 0.3f, sinf(a) * nucR * 0.95f};
        ci.phys.addParticle(PT_NUCLEAR_PORE, p,
            molRadiusWu(MoleculeRadiusUm::NUCLEAR_PORE, MOL_VIS_BOOST_MACROMOL),
            0.4f, 0.7f, 0.4f);
    }

    // Ribosomes — 1000 visible, real HeLa ~6 × 10⁶.
    simd_float3 roughERPos = {0, 0, cellR * 0.35f};
    for (int i = 0; i < 1000; i++) {
        float a = (float)i / 60.0f * 2.0f * M_PI;
        float rr = cellR * 0.15f + ipRandf() * cellR * 0.40f;
        simd_float3 rPos = {roughERPos.x + cosf(a) * rr,
                            (ipRandf() - 0.5f) * cellR * 0.35f,
                            roughERPos.z + sinf(a) * rr * 0.7f};
        int idS = ci.phys.addParticle(PT_RIBOSOME_SMALL, rPos,
            molRadiusWu(MoleculeRadiusUm::RIBOSOME_SMALL, MOL_VIS_BOOST_MACROMOL),
            0.5f, 0.4f, 0.75f);
        ci.phys.setHome(idS, rPos, 0.4f);
        simd_float3 o = {0,0,0}; ci.phys.setConfinement(idS, o, cellR);
        simd_float3 lPos = {rPos.x, rPos.y + cs * 0.015f, rPos.z};
        int idL = ci.phys.addParticle(PT_RIBOSOME_LARGE, lPos,
            molRadiusWu(MoleculeRadiusUm::RIBOSOME_LARGE, MOL_VIS_BOOST_MACROMOL),
            0.55f, 0.45f, 0.8f);
        ci.phys.setHome(idL, lPos, 0.4f);
        ci.phys.setConfinement(idL, o, cellR);
        ci.phys.setTether(idL, idS, cs * 0.015f, 2.0f);
    }

    // tRNA (1500 — visible sample, real ~10⁷)
    for (int i = 0; i < 1500; i++) {
        float u2 = ipRandf()*2-1, phi2 = ipRandf()*2*M_PI;
        float r2 = cbrtf(ipRandf()) * cellR * 0.85f;
        float st = sqrtf(1-u2*u2);
        simd_float3 p = {r2*st*cosf(phi2), r2*u2, r2*st*sinf(phi2)};
        int id = ci.phys.addParticle(PT_TRNA, p,
            molRadiusWu(MoleculeRadiusUm::TRNA, MOL_VIS_BOOST_MACROMOL),
            0.3f, 0.8f, 0.4f);
        simd_float3 o = {0,0,0}; ci.phys.setConfinement(id, o, cellR);
    }

    // Chaperones (1000 visible — real ~10⁶)
    for (int i = 0; i < 1000; i++) {
        simd_float3 p = {(ipRandf()-0.5f)*cellR*1.5f,
                         (ipRandf()-0.5f)*cellR*1.5f,
                         (ipRandf()-0.5f)*cellR*1.5f};
        int id = ci.phys.addParticle(PT_CHAPERONE, p,
            molRadiusWu(MoleculeRadiusUm::CHAPERONE, MOL_VIS_BOOST_MACROMOL),
            0.7f, 0.5f, 0.3f);
        ci.phys.setHome(id, p, 0.3f);
        simd_float3 o = {0,0,0}; ci.phys.setConfinement(id, o, cellR);
    }

    // COPII vesicles (300 visible — real ~10³)
    simd_float3 golgiPos = {0, 0, -cellR * 0.42f};
    for (int i = 0; i < 300; i++) {
        float u2 = ipRandf()*2-1, phi2 = ipRandf()*2*M_PI;
        float r2 = cbrtf(ipRandf()) * cellR * 0.1f;
        float st = sqrtf(1-u2*u2);
        simd_float3 p = {roughERPos.x + r2*st*cosf(phi2), r2*u2, roughERPos.z + r2*st*sinf(phi2)};
        int id = ci.phys.addParticle(PT_VESICLE_COPII, p,
            molRadiusWu(MoleculeRadiusUm::COPII_VESICLE, MOL_VIS_BOOST_MACROMOL),
            0.8f, 0.7f, 0.2f);
        simd_float3 o = {0,0,0}; ci.phys.setConfinement(id, o, cellR);
        ci.phys.particles[id].spawnPos = roughERPos;
        ci.phys.assignGoto(id, golgiPos, 0.012f);
    }

    // Secretory vesicles (600 visible — real ~10³-10⁴)
    for (int i = 0; i < 600; i++) {
        float u2 = ipRandf()*2-1, phi2 = ipRandf()*2*M_PI;
        float r2 = cbrtf(ipRandf()) * cellR * 0.08f;
        float st = sqrtf(1-u2*u2);
        simd_float3 p = {golgiPos.x + r2*st*cosf(phi2), r2*u2, golgiPos.z + r2*st*sinf(phi2)};
        int id = ci.phys.addParticle(PT_VESICLE_SECRETORY, p,
            molRadiusWu(MoleculeRadiusUm::SECRETORY_VESICLE, MOL_VIS_BOOST_MACROMOL),
            0.6f, 0.8f, 0.4f);
        simd_float3 o = {0,0,0}; ci.phys.setConfinement(id, o, cellR);
        ci.phys.particles[id].spawnPos = golgiPos;
        float angle = ipRandf() * 2.0f * M_PI;
        simd_float3 memTarget = {cosf(angle)*cellR*0.92f, 0, sinf(angle)*cellR*0.92f};
        ci.phys.assignGoto(id, memTarget, 0.01f);
    }

    // Free molecules — radii anchored to real molecular dimensions (µm)
    // via MoleculeRadiusUm + MOL_VIS_BOOST_SMALL (30×) so atoms remain
    // hoverable at normal camera distance. Physics uses the rendered radius
    // for spatial-hash spacing, so the boost also avoids sub-grid spurious
    // collisions. Sources: PDB SDF files in assets/molecules/.
    struct MolSpec { float r, g, b; int count; float radius; simd_float3 center; float spread; MolTag tag; };
    simd_float3 mitoPos = {cellR*0.3f, 0, 0};
    const float rATP     = molRadiusWu(MoleculeRadiusUm::ATP,        MOL_VIS_BOOST_SMALL);
    const float rGlucose = molRadiusWu(MoleculeRadiusUm::GLUCOSE,    MOL_VIS_BOOST_SMALL);
    const float rNADH    = molRadiusWu(MoleculeRadiusUm::NADH,       MOL_VIS_BOOST_SMALL);
    const float rAA      = molRadiusWu(MoleculeRadiusUm::AMINO_ACID, MOL_VIS_BOOST_SMALL);
    const float rCAMP    = molRadiusWu(MoleculeRadiusUm::CAMP,       MOL_VIS_BOOST_SMALL);
    const float rWater   = molRadiusWu(MoleculeRadiusUm::WATER,      MOL_VIS_BOOST_SMALL);
    const float rCalcium = molRadiusWu(MoleculeRadiusUm::CALCIUM_ION, MOL_VIS_BOOST_SMALL);
    // Background cell metabolite sample — smaller than primary (100-
    // 800 range) because far-away cells don't need the full visual
    // density. Still orders of magnitude above the previous 20-150.
    MolSpec specs[] = {
        {0.2f, 0.8f, 0.3f,  800, rATP,     mitoPos,    cellR*0.85f, TAG_ATP},
        {0.9f, 0.7f, 0.1f,  600, rGlucose, nucleusPos, cellR*0.85f, TAG_GLUCOSE},
        {0.1f, 0.3f, 0.9f,  300, rNADH,    mitoPos,    cellR*0.55f, TAG_NADH},
        {0.9f, 0.3f, 0.3f,  800, rAA,      roughERPos, cellR*0.80f, TAG_AMINO_ACID},
        {0.5f, 0.9f, 0.5f,  100, rCAMP,    nucleusPos, cellR*0.85f, TAG_CAMP},
        {0.8f, 0.85f, 0.9f,1200, rWater,   nucleusPos, cellR*0.90f, TAG_WATER},
        {0.3f, 0.7f, 0.7f,  200, rCalcium, roughERPos, cellR*0.70f, TAG_CALCIUM},
    };
    for (auto& spec : specs) {
        for (int i = 0; i < spec.count; i++) {
            float u2 = ipRandf() * 2.0f - 1.0f, phi2 = ipRandf() * 2.0f * M_PI;
            float r2 = cbrtf(ipRandf()) * spec.spread;
            float sinTh = sqrtf(1.0f - u2 * u2);
            simd_float3 p = {spec.center.x + r2 * sinTh * cosf(phi2),
                             spec.center.y + r2 * u2,
                             spec.center.z + r2 * sinTh * sinf(phi2)};
            int id = ci.phys.addParticle(PT_MOLECULE, p, spec.radius, spec.r, spec.g, spec.b);
            simd_float3 o = {0,0,0}; ci.phys.setConfinement(id, o, cellR);
            ci.phys.particles[id].tag = spec.tag;
            ci.phys.particles[id].spawnPos = spec.center;
        }
    }

    // Attraction fields — same as primary cell
    ci.phys.attractionFields.clear();
    ci.phys.attractionFields.push_back({mitoPos, cellR*0.4f, 0.015f, TAG_NADH, PT_COUNT});
    ci.phys.attractionFields.push_back({mitoPos, cellR*0.3f, 0.008f, TAG_GLUCOSE, PT_COUNT});
    ci.phys.attractionFields.push_back({roughERPos, cellR*0.35f, 0.012f, TAG_AMINO_ACID, PT_COUNT});
    ci.phys.attractionFields.push_back({roughERPos, cellR*0.4f, 0.010f, TAG_NONE, PT_TRNA});
    ci.phys.attractionFields.push_back({{0,0,0}, nucR*1.2f, 0.020f, TAG_NONE, PT_RNA_POLYMERASE});
    ci.phys.attractionFields.push_back({roughERPos, cellR*0.3f, 0.010f, TAG_CALCIUM, PT_COUNT});
    ci.phys.attractionFields.push_back({roughERPos, cellR*0.3f, 0.008f, TAG_NONE, PT_CHAPERONE});
    ci.phys.attractionFields.push_back({golgiPos, cellR*0.3f, 0.012f, TAG_NONE, PT_VESICLE_COPII});

    ci.initialized = true;
}

static simd_float3 randNearInit(simd_float3 c, float r) {
    float u = ipRandf()*2-1, phi = ipRandf()*2*M_PI;
    float rad = cbrtf(ipRandf()) * r;
    float st = sqrtf(1-u*u);
    return {c.x + rad*st*cosf(phi), c.y + rad*u, c.z + rad*st*sinf(phi)};
}

static void initIntracellularParticles(simd_float3 cellPos, float cellSize) {
    float cs = cellSize;
    float cellR = cs;
    // Nucleus confinement = 80% of visual nucleus radius
    // Nucleus GLB: bbox=8.93, scale=0.10 → visual radius = cs * 0.10 * 4.47 = cs * 0.447
    // DNA confinement = 80% of that = cs * 0.358
    float nucR  = cs * 0.45f; // ~cs * 0.358

    // ALL positions in LOCAL space: cell center = origin {0,0,0}
    // World = cellPos + local (added only at render time)
    simd_float3 origin = {0, 0, 0};
    gIntraPhys.init(origin, cellR);
    gLastCellPos = cellPos;
    gCDogma.init();

    simd_float3 nucleusPos = origin;
    float innerR = cellR * 0.5f;
    simd_float3 roughERPos = {0, 0, innerR * 0.55f};
    simd_float3 golgiPos   = {0, 0, -innerR * 0.6f};
    simd_float3 mitoPos[3] = {
        {innerR * 0.45f, innerR * 0.15f, innerR * 0.25f},
        {-innerR * 0.35f, -innerR * 0.15f, 0},
        {-innerR * 0.25f, 0, -innerR * 0.35f}
    };

    // ── DNA double helix (real HBB sequence) ────────────────────────
    // One bp every ~3.4Å, we use ~2bp per visual node for 763bp → ~380 nodes
    int dnaStride = 2;  // 2 bp per visual node
    int dnaNodes = HBB_LENGTH / dnaStride;
    float helixR = nucR * 0.3f;
    float helixH = nucR * 1.6f;
    int turns = 8;

    int firstDnaA = -1;
    for (int i = 0; i < dnaNodes; i++) {
        float t = (float)i / dnaNodes;
        float angle = t * turns * 2.0f * M_PI;
        float y = nucleusPos.y + (t - 0.5f) * helixH;

        int bpIdx = i * dnaStride;
        bool ex = (bpIdx < HBB_LENGTH) ? isExon(bpIdx) : false;
        char base = (bpIdx < HBB_LENGTH) ? HBB_SEQUENCE[bpIdx] : 'A';
        BaseColor bc = getBaseColor(base, ex);
        BaseColor bc2 = getBaseColor(complement(base), ex);

        // DNA node ≈ 10 nm nucleosome (Luger 1997 Nature) — see MoleculeRadiusUm.
        const float rDNA = molRadiusWu(MoleculeRadiusUm::DNA_NODE, MOL_VIS_BOOST_MACROMOL);
        // Strand A
        simd_float3 pA = {nucleusPos.x + cosf(angle) * helixR, y,
                          nucleusPos.z + sinf(angle) * helixR};
        int idA = gIntraPhys.addParticle(PT_DNA_NODE, pA, rDNA, bc.r, bc.g, bc.b);
        gIntraPhys.setHome(idA, pA, 2.0f);
        gIntraPhys.setConfinement(idA, nucleusPos, nucR);
        if (firstDnaA < 0) firstDnaA = idA;

        // Strand B (complementary)
        simd_float3 pB = {nucleusPos.x + cosf(angle + M_PI) * helixR, y,
                          nucleusPos.z + sinf(angle + M_PI) * helixR};
        int idB = gIntraPhys.addParticle(PT_DNA_NODE, pB, rDNA, bc2.r, bc2.g, bc2.b);
        gIntraPhys.setHome(idB, pB, 2.0f);
        gIntraPhys.setConfinement(idB, nucleusPos, nucR);

        // Tether A to previous A (backbone spring)
        if (i > 0) {
            gIntraPhys.setTether(idA, idA - 2, cs * 0.025f, 5.0f);
            gIntraPhys.setTether(idB, idB - 2, cs * 0.025f, 5.0f);
        }

        // Store bp index for transcription lookup
        gIntraPhys.particles[idA].stateIndex = bpIdx;
        gIntraPhys.particles[idB].stateIndex = bpIdx;
    }

    // ── RNA Polymerase II (3 instances) ─────────────────────────────
    for (int i = 0; i < 3; i++) {
        float t = (float)i * 0.3f;
        simd_float3 p = {nucleusPos.x + sinf(t * 4) * nucR * 0.2f,
                         nucleusPos.y + (t - 0.5f) * helixH * 0.5f,
                         nucleusPos.z + cosf(t * 4) * nucR * 0.2f};
        int id = gIntraPhys.addParticle(PT_RNA_POLYMERASE, p,
            molRadiusWu(MoleculeRadiusUm::RNA_POL, MOL_VIS_BOOST_MACROMOL),
            0.6f, 0.2f, 0.8f);
        gIntraPhys.setConfinement(id, nucleusPos, nucR);
        gIntraPhys.particles[id].glowIntensity = 0.8f;
    }

    // ── DNA Polymerase (4 small replisome markers in nucleus) ───────
    for (int i = 0; i < 4; i++) {
        float angle = (float)i / 4.0f * 2.0f * M_PI;
        simd_float3 p = {
            nucleusPos.x + cosf(angle) * nucR * 0.12f,
            nucleusPos.y + ((i < 2) ? -0.04f : 0.04f) * cs,
            nucleusPos.z + sinf(angle) * nucR * 0.12f
        };
        int id = gIntraPhys.addParticle(PT_DNA_POLYMERASE, p,
            molRadiusWu(MoleculeRadiusUm::DNA_POL, MOL_VIS_BOOST_MACROMOL),
            0.55f, 0.78f, 1.0f);
        gIntraPhys.setConfinement(id, nucleusPos, nucR);
        gIntraPhys.setHome(id, p, 0.7f);
        gIntraPhys.particles[id].glowIntensity = 0.55f;
    }

    // ── Spliceosome (2 instances) ───────────────────────────────────
    for (int i = 0; i < 2; i++) {
        simd_float3 p = {nucleusPos.x + (i ? 0.15f : -0.15f) * cs,
                         nucleusPos.y + 0.1f * cs,
                         nucleusPos.z + (i ? 0.1f : -0.1f) * cs};
        int id = gIntraPhys.addParticle(PT_SPLICEOSOME, p,
            molRadiusWu(MoleculeRadiusUm::SPLICEOSOME, MOL_VIS_BOOST_MACROMOL),
            0.8f, 0.6f, 0.2f);
        gIntraPhys.setConfinement(id, nucleusPos, nucR);
        gIntraPhys.setHome(id, p, 0.5f);
    }

    // ── Nuclear pores (6 on nucleus surface) ────────────────────────
    for (int i = 0; i < 6; i++) {
        float angle = (float)i / 6.0f * 2.0f * M_PI;
        simd_float3 p = {nucleusPos.x + cosf(angle) * nucR * 0.95f,
                         nucleusPos.y + sinf(angle * 0.5f) * nucR * 0.3f,
                         nucleusPos.z + sinf(angle) * nucR * 0.95f};
        int id = gIntraPhys.addParticle(PT_NUCLEAR_PORE, p,
            molRadiusWu(MoleculeRadiusUm::NUCLEAR_PORE, MOL_VIS_BOOST_MACROMOL),
            0.4f, 0.7f, 0.4f);
        gIntraPhys.particles[id].glowIntensity = 0.3f;
    }

    // ── Ribosomes — push visible count up to the max the particle
    //    system can carry per cell. Real HeLa has ~6 × 10⁶; we render
    //    a dense sample of 1500 here. Each represents ~4000 real
    //    ribosomes. Perf cap: 1500 × 1 cell = 1500 particles.
    for (int i = 0; i < 1500; i++) {
        float angle = (float)i / 1500.0f * 2.0f * M_PI;
        float rr = cellR * 0.15f + ipRandf() * cellR * 0.65f;
        simd_float3 rPos = {roughERPos.x + cosf(angle) * rr,
                            roughERPos.y + (ipRandf() - 0.5f) * cellR * 0.75f,
                            roughERPos.z + sinf(angle) * rr * 0.9f};
        int idS = gIntraPhys.addParticle(PT_RIBOSOME_SMALL, rPos,
            molRadiusWu(MoleculeRadiusUm::RIBOSOME_SMALL, MOL_VIS_BOOST_MACROMOL),
            0.5f, 0.4f, 0.75f);
        gIntraPhys.setHome(idS, rPos, 0.4f);
        gIntraPhys.setConfinement(idS, cellPos, cellR);
        simd_float3 lPos = {rPos.x, rPos.y + cs * 0.015f, rPos.z};
        int idL = gIntraPhys.addParticle(PT_RIBOSOME_LARGE, lPos,
            molRadiusWu(MoleculeRadiusUm::RIBOSOME_LARGE, MOL_VIS_BOOST_MACROMOL),
            0.55f, 0.45f, 0.8f);
        gIntraPhys.setHome(idL, lPos, 0.4f);
        gIntraPhys.setConfinement(idL, cellPos, cellR);
        gIntraPhys.setTether(idL, idS, cs * 0.015f, 2.0f);
    }

    // ── tRNA molecules — sample of 2500 out of real ~10⁷ per cell.
    //    Each visible bead represents ~4000 real tRNAs.
    for (int i = 0; i < 2500; i++) {
        float u = ipRandf() * 2.0f - 1.0f;
        float phi = ipRandf() * 2.0f * M_PI;
        float r = cbrtf(ipRandf()) * cellR * 0.85f; // spread full cell
        float sinTh = sqrtf(1.0f - u * u);
        simd_float3 p = {r * sinTh * cosf(phi),
                         r * u,
                         r * sinTh * sinf(phi)};
        int id = gIntraPhys.addParticle(PT_TRNA, p,
            molRadiusWu(MoleculeRadiusUm::TRNA, MOL_VIS_BOOST_MACROMOL),
            0.3f, 0.8f, 0.4f);
        simd_float3 o={0,0,0}; gIntraPhys.setConfinement(id, o, cellR);
        gIntraPhys.particles[id].glowIntensity = 0.4f;
    }

    // ── Chaperones — 1500 visible; real Hsp70+90+BiP ~10⁶.
    //    Each bead represents ~700 real chaperones.
    for (int i = 0; i < 1500; i++) {
        simd_float3 p = {(ipRandf() - 0.5f) * cellR * 1.5f,
                         (ipRandf() - 0.5f) * cellR * 1.5f,
                         (ipRandf() - 0.5f) * cellR * 1.5f};
        int id = gIntraPhys.addParticle(PT_CHAPERONE, p,
            molRadiusWu(MoleculeRadiusUm::CHAPERONE, MOL_VIS_BOOST_MACROMOL),
            0.7f, 0.5f, 0.3f);
        gIntraPhys.setHome(id, p, 0.3f);
        simd_float3 o={0,0,0}; gIntraPhys.setConfinement(id, o, cellR);
    }

    // ── COPII vesicles — 400 visible, closer to real 10³ per cell.
    for (int i = 0; i < 400; i++) {
        simd_float3 p = randNearInit(roughERPos, cellR * 0.1f);
        int id = gIntraPhys.addParticle(PT_VESICLE_COPII, p,
            molRadiusWu(MoleculeRadiusUm::COPII_VESICLE, MOL_VIS_BOOST_MACROMOL),
            0.8f, 0.7f, 0.2f);
        simd_float3 o={0,0,0}; gIntraPhys.setConfinement(id, o, cellR);
        gIntraPhys.particles[id].spawnPos = roughERPos;
        gIntraPhys.assignGoto(id, randNearInit(golgiPos, cellR*0.05f), 0.012f);
    }

    // ── Secretory vesicles — 800 visible (real pool ~10³-10⁴)
    for (int i = 0; i < 800; i++) {
        simd_float3 p = randNearInit(golgiPos, cellR * 0.08f);
        int id = gIntraPhys.addParticle(PT_VESICLE_SECRETORY, p,
            molRadiusWu(MoleculeRadiusUm::SECRETORY_VESICLE, MOL_VIS_BOOST_MACROMOL),
            0.6f, 0.8f, 0.4f);
        simd_float3 o={0,0,0}; gIntraPhys.setConfinement(id, o, cellR);
        gIntraPhys.particles[id].spawnPos = golgiPos;
        float angle = ipRandf() * 2.0f * M_PI;
        simd_float3 memTarget = {cosf(angle)*cellR*0.92f, 0, sinf(angle)*cellR*0.92f};
        gIntraPhys.assignGoto(id, memTarget, 0.01f);
    }

    // ── Free molecules (200+ tiny dots, with tags for job system) ──────
    // Radii anchored to PDB dimensions via MoleculeRadiusUm (µm) plus
    // MOL_VIS_BOOST_SMALL so 1 nm metabolites remain visible in a 10 µm
    // cell at normal camera distance.
    struct MolSpec {
        float r, g, b; int count; float radius;
        simd_float3 center; float spread;
        MolTag tag;
    };
    const float rATP     = molRadiusWu(MoleculeRadiusUm::ATP,        MOL_VIS_BOOST_SMALL);
    const float rGlucose = molRadiusWu(MoleculeRadiusUm::GLUCOSE,    MOL_VIS_BOOST_SMALL);
    const float rNADH    = molRadiusWu(MoleculeRadiusUm::NADH,       MOL_VIS_BOOST_SMALL);
    const float rAA      = molRadiusWu(MoleculeRadiusUm::AMINO_ACID, MOL_VIS_BOOST_SMALL);
    const float rCAMP    = molRadiusWu(MoleculeRadiusUm::CAMP,       MOL_VIS_BOOST_SMALL);
    const float rWater   = molRadiusWu(MoleculeRadiusUm::WATER,      MOL_VIS_BOOST_SMALL);
    const float rCalcium = molRadiusWu(MoleculeRadiusUm::CALCIUM_ION, MOL_VIS_BOOST_SMALL);
    // Metabolite visible sample — pushed much higher so the cytoplasm
    // is packed with particles. Real counts are 10⁹-10¹⁴ which is not
    // renderable, so each visible bead represents many real molecules.
    // Scale factors (visible : real): ATP 1:2M, glucose 1:6M,
    // NADH 1:300k, AA 1:2M, water 1:4×10¹⁰, Ca 1:500.
    MolSpec molSpecs[] = {
        {0.2f,  0.8f,  0.3f, 2500, rATP,     mitoPos[0], cellR*0.85f, TAG_ATP},
        {0.2f,  0.8f,  0.3f, 1200, rATP,     mitoPos[1], cellR*0.55f, TAG_ATP},
        {0.9f,  0.7f,  0.1f, 2000, rGlucose, nucleusPos, cellR*0.85f, TAG_GLUCOSE},
        {0.1f,  0.3f,  0.9f,  800, rNADH,    mitoPos[2], cellR*0.55f, TAG_NADH},
        {0.9f,  0.3f,  0.3f, 2500, rAA,      roughERPos, cellR*0.80f, TAG_AMINO_ACID},
        {0.5f,  0.9f,  0.5f,  300, rCAMP,    nucleusPos, cellR*0.85f, TAG_CAMP},
        {0.8f,  0.85f, 0.9f, 4000, rWater,   nucleusPos, cellR*0.90f, TAG_WATER},
        {0.3f,  0.7f,  0.7f,  500, rCalcium, roughERPos, cellR*0.70f, TAG_CALCIUM},
    };
    for (auto& spec : molSpecs) {
        for (int i = 0; i < spec.count; i++) {
            float u2 = ipRandf() * 2.0f - 1.0f;
            float phi2 = ipRandf() * 2.0f * M_PI;
            float r2 = cbrtf(ipRandf()) * spec.spread;
            float sinTh2 = sqrtf(1.0f - u2 * u2);
            simd_float3 p = {spec.center.x + r2 * sinTh2 * cosf(phi2),
                             spec.center.y + r2 * u2,
                             spec.center.z + r2 * sinTh2 * sinf(phi2)};
            int id = gIntraPhys.addParticle(PT_MOLECULE, p, spec.radius, spec.r, spec.g, spec.b);
            simd_float3 o={0,0,0}; gIntraPhys.setConfinement(id, o, cellR);
            gIntraPhys.particles[id].tag = spec.tag;
            gIntraPhys.particles[id].spawnPos = spec.center; // respawn near origin
        }
    }

    // ── Plasma membrane mosaic — fluid-mosaic model (Singer-Nicolson 1972) ─
    // Real membrane has ~5 million phospholipids per µm² of surface. Our
    // 10 µm cell has ~1300 µm² surface → 6.5 × 10⁹ lipids. We render a
    // dense visible sample: 2500 phospholipid head dots scattered on the
    // membrane surface + several hundred embedded receptor particles
    // (GPCRs, RTKs, cytokine, ion channels, pumps, adhesion). Each is
    // clamped to the membrane sphere (radius 0.98 × cellR). They diffuse
    // in 2D via existing Langevin with confinement at 0.98R.
    // Ref: Alberts MBoC 7e Ch 10 (membrane structure) and Ch 15 (signaling).
    {
        // Place membrane proteins SLIGHTLY OUTSIDE the membrane sphere
        // (at 1.02 × cellR) so they stick out through the translucent
        // blue cell membrane and are clearly visible as dots on the
        // outside fluid-mosaic surface. The cell render sphere is at
        // exactly cellR, so 1.02 puts them just on the outer face.
        float memR = cellR * 1.02f;
        auto sampleMembranePoint = [&]() -> simd_float3 {
            // Uniform on sphere surface (Marsaglia 1972)
            float u = ipRandf() * 2.0f - 1.0f;
            float phi = ipRandf() * 2.0f * M_PI;
            float st = sqrtf(fmaxf(0.0f, 1.0f - u * u));
            return {st * cosf(phi) * memR, u * memR, st * sinf(phi) * memR};
        };
        // Helper to spawn a membrane-embedded particle with distinctive
        // color, bigger radius, and strong emissive so they pop against
        // the translucent organelle haze.
        auto spawnMem = [&](ParticleType type, int count, float relRadius,
                             float r, float g, float b, float glow) {
            for (int i = 0; i < count; i++) {
                simd_float3 p = sampleMembranePoint();
                int id = gIntraPhys.addParticle(type, p,
                    cs * relRadius, r, g, b);
                simd_float3 zero = {0, 0, 0};
                gIntraPhys.setConfinement(id, zero, memR);
                gIntraPhys.setHome(id, p, 0.10f);
                gIntraPhys.particles[id].glowIntensity = glow;
            }
        };
        // Low-key aesthetic palette: deep, desaturated tones, low glow.
        // Reduced counts (~40%) so the mosaic reads as a sparse field
        // of muted dots — elegant rather than busy.
        spawnMem(PT_PHOSPHOLIPID,    2800, 0.005f, 0.28f, 0.38f, 0.48f, 0.10f);
        spawnMem(PT_RECEPTOR_GPCR,    250, 0.010f, 0.22f, 0.32f, 0.55f, 0.25f); // deep blue-slate
        spawnMem(PT_RECEPTOR_RTK,     180, 0.010f, 0.55f, 0.45f, 0.18f, 0.25f); // burnt gold
        spawnMem(PT_RECEPTOR_CYTOKINE,120, 0.010f, 0.50f, 0.28f, 0.42f, 0.22f); // dusky rose
        spawnMem(PT_RECEPTOR_STK,     100, 0.010f, 0.42f, 0.38f, 0.20f, 0.20f); // olive-dark
        spawnMem(PT_RECEPTOR_DEATH,    50, 0.010f, 0.55f, 0.22f, 0.22f, 0.28f); // rust
        spawnMem(PT_RECEPTOR_TLR,      80, 0.010f, 0.22f, 0.25f, 0.50f, 0.22f); // indigo
        spawnMem(PT_RECEPTOR_NOTCH,    80, 0.010f, 0.28f, 0.48f, 0.40f, 0.20f); // forest mint
        spawnMem(PT_RECEPTOR_FRIZZLED, 80, 0.010f, 0.50f, 0.35f, 0.50f, 0.22f); // dark mauve
        spawnMem(PT_RECEPTOR_PATCHED,  50, 0.010f, 0.55f, 0.35f, 0.25f, 0.22f); // terracotta
        spawnMem(PT_LDL_RECEPTOR,      80, 0.009f, 0.42f, 0.30f, 0.22f, 0.18f); // bronze
        spawnMem(PT_TRANSFERRIN_R,     80, 0.009f, 0.52f, 0.38f, 0.28f, 0.18f); // sienna
        spawnMem(PT_ION_CHANNEL,      320, 0.011f, 0.25f, 0.48f, 0.32f, 0.28f); // pine
        spawnMem(PT_AQUAPORIN,        200, 0.008f, 0.30f, 0.48f, 0.58f, 0.20f); // slate blue
        spawnMem(PT_ION_PUMP,         280, 0.011f, 0.55f, 0.50f, 0.22f, 0.28f); // antique brass
        spawnMem(PT_EXCHANGER,        140, 0.010f, 0.50f, 0.38f, 0.22f, 0.22f); // umber
        spawnMem(PT_TRANSPORTER,      220, 0.010f, 0.45f, 0.32f, 0.22f, 0.20f); // walnut
        spawnMem(PT_GAP_JUNCTION,      60, 0.012f, 0.42f, 0.45f, 0.48f, 0.20f); // graphite
        spawnMem(PT_ADHESION,         180, 0.010f, 0.42f, 0.30f, 0.55f, 0.22f); // dark lavender
        spawnMem(PT_MHC,              150, 0.009f, 0.55f, 0.50f, 0.42f, 0.18f); // parchment
    }

    // ── Set up organelle attraction fields ────────────────────────
    gIntraPhys.attractionFields.clear();
    // Mitochondria attract NADH, ATP (recycle), glucose
    for (int m = 0; m < 3; m++) {
        gIntraPhys.attractionFields.push_back({mitoPos[m], cellR*0.4f, 0.015f, TAG_NADH, PT_COUNT});
        gIntraPhys.attractionFields.push_back({mitoPos[m], cellR*0.3f, 0.008f, TAG_GLUCOSE, PT_COUNT});
    }
    // Ribosomes attract amino acids and tRNA
    gIntraPhys.attractionFields.push_back({roughERPos, cellR*0.35f, 0.012f, TAG_AMINO_ACID, PT_COUNT});
    gIntraPhys.attractionFields.push_back({roughERPos, cellR*0.4f, 0.010f, TAG_NONE, PT_TRNA});
    // Nucleus attracts RNA Pol
    gIntraPhys.attractionFields.push_back({nucleusPos, nucR*1.2f, 0.020f, TAG_NONE, PT_RNA_POLYMERASE});
    // ER attracts calcium and chaperones
    gIntraPhys.attractionFields.push_back({roughERPos, cellR*0.3f, 0.010f, TAG_CALCIUM, PT_COUNT});
    gIntraPhys.attractionFields.push_back({roughERPos, cellR*0.3f, 0.008f, TAG_NONE, PT_CHAPERONE});
    // Golgi attracts COPII vesicles
    gIntraPhys.attractionFields.push_back({golgiPos, cellR*0.3f, 0.012f, TAG_NONE, PT_VESICLE_COPII});

    gIntraInitialized = true;
    printf("[CellSim] Intracellular physics: %d particles, %d attraction fields\n",
           (int)gIntraPhys.particles.size(), (int)gIntraPhys.attractionFields.size());
}

static void driveInteriorActivity(IntracellularPhysics& phys, float cellR, float nucR,
                                  int cellKey, float time) {
    const float simDt = 1.0f / 60.0f;
    simd_float3 localOrigin = {0, 0, 0};
    simd_float3 nucPos = localOrigin;
    simd_float3 erPos = {0, 0, cellR * 0.35f};
    simd_float3 golgiP = {0, 0, -cellR * 0.42f};
    simd_float3 mitoP[3] = {
        {cellR*0.3f, cellR*0.12f, cellR*0.18f},
        {-cellR*0.24f, -cellR*0.12f, 0},
        {-cellR*0.18f, 0, -cellR*0.24f}
    };
    if (!gMitosis.active && cellKey >= 0) {
        OrganelleMotionState& motion = updateOrganelleMotionState(cellKey, time);
        golgiP = {motion.golgiPos.x * cellR, motion.golgiPos.y * cellR, motion.golgiPos.z * cellR};
        for (int m = 0; m < 3; m++) {
            mitoP[m] = {motion.mitoPos[m].x * cellR,
                        motion.mitoPos[m].y * cellR,
                        motion.mitoPos[m].z * cellR};
        }
    }

    auto randNear = [](simd_float3 c, float r) -> simd_float3 {
        float u = ipRandf()*2-1, phi = ipRandf()*2*M_PI;
        float rad = cbrtf(ipRandf()) * r;
        float st = sqrtf(1-u*u);
        return {c.x + rad*st*cosf(phi), c.y + rad*u, c.z + rad*st*sinf(phi)};
    };

    phys.attractionFields.clear();
    for (int m = 0; m < 3; m++) {
        phys.attractionFields.push_back({mitoP[m], cellR*0.4f, 0.015f, TAG_NADH, PT_COUNT});
        phys.attractionFields.push_back({mitoP[m], cellR*0.3f, 0.008f, TAG_GLUCOSE, PT_COUNT});
    }
    phys.attractionFields.push_back({erPos, cellR*0.35f, 0.012f, TAG_AMINO_ACID, PT_COUNT});
    phys.attractionFields.push_back({erPos, cellR*0.4f, 0.010f, TAG_NONE, PT_TRNA});
    phys.attractionFields.push_back({nucPos, nucR*1.2f, 0.020f, TAG_NONE, PT_RNA_POLYMERASE});
    phys.attractionFields.push_back({erPos, cellR*0.3f, 0.010f, TAG_CALCIUM, PT_COUNT});
    phys.attractionFields.push_back({erPos, cellR*0.3f, 0.008f, TAG_NONE, PT_CHAPERONE});
    phys.attractionFields.push_back({golgiP, cellR*0.3f, 0.012f, TAG_NONE, PT_VESICLE_COPII});

    auto findNearest = [&](simd_float3 from, ParticleType type, MolJob job) -> int {
        int best = -1;
        float bestD = 1e9f;
        for (int i = 0; i < (int)phys.particles.size(); i++) {
            auto& q = phys.particles[i];
            if (!q.active || q.type != type || q.job != job) continue;
            float dx = q.position.x - from.x;
            float dy = q.position.y - from.y;
            float dz = q.position.z - from.z;
            float d = dx*dx + dy*dy + dz*dz;
            if (d < bestD) { bestD = d; best = i; }
        }
        return best;
    };

    for (int i = 0; i < (int)phys.particles.size(); i++) {
        auto& p = phys.particles[i];
        if (!p.active) continue;

        if (p.type == PT_MOLECULE && p.tag == TAG_ATP) {
            if (p.job == JOB_IDLE) {
                simd_float3 dest;
                float r = ipRandf();
                if (r < 0.4f) {
                    int rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_IDLE);
                    dest = (rib >= 0) ? phys.particles[rib].position : randNear(erPos, cellR*0.3f);
                } else if (r < 0.7f) {
                    dest = randNear(erPos, cellR * 0.2f);
                } else {
                    dest = randNear(localOrigin, cellR * 0.6f);
                }
                phys.assignGoto(i, dest, 0.025f);
            } else if (p.job == JOB_GOTO_TARGET && phys.reachedTarget(i, 0.08f)) {
                p.glowIntensity = 2.0f;
                phys.assignConsume(i);
                p.spawnPos = mitoP[(int)(ipRandf() * 2.99f)];
            }
        }

        if (p.type == PT_MOLECULE && p.tag == TAG_GLUCOSE) {
            if (p.job == JOB_IDLE) {
                phys.assignGoto(i, randNear(localOrigin, cellR * 0.3f), 0.015f);
            } else if (p.job == JOB_GOTO_TARGET && phys.reachedTarget(i, 0.08f)) {
                phys.assignConsume(i);
                float a = ipRandf() * 2 * M_PI;
                p.spawnPos = {cosf(a)*cellR*0.88f, 0, sinf(a)*cellR*0.88f};
            }
        }

        if (p.type == PT_MOLECULE && p.tag == TAG_AMINO_ACID) {
            if (p.job == JOB_IDLE) {
                int rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_FOLLOW_PATH);
                if (rib >= 0) phys.assignGoto(i, phys.particles[rib].position, 0.02f);
                else phys.assignGoto(i, randNear(erPos, cellR * 0.25f), 0.008f);
            } else if (p.job == JOB_GOTO_TARGET && phys.reachedTarget(i, 0.06f)) {
                phys.assignConsume(i);
                p.spawnPos = randNear(erPos, cellR * 0.25f);
            }
        }

        if (p.type == PT_TRNA) {
            if (p.job == JOB_IDLE) {
                int rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_FOLLOW_PATH);
                if (rib < 0) rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_IDLE);
                if (rib >= 0) phys.assignGoto(i, phys.particles[rib].position, 0.03f);
                else phys.assignGoto(i, randNear(erPos, cellR * 0.3f), 0.01f);
            } else if (p.job == JOB_GOTO_TARGET && phys.reachedTarget(i, 0.07f)) {
                p.job = JOB_IDLE;
                p.jobTimer = 0;
            }
        }

        if (p.type == PT_RNA_POLYMERASE) {
            if (p.job == JOB_IDLE && p.jobTimer > 2.0f) {
                phys.assignFollowPath(i, 0.008f);
            }
            if (p.job == JOB_FOLLOW_PATH) {
                p.jobProgress += simDt * p.jobSpeed;
                int dnaCount = 0;
                for (auto& q : phys.particles) {
                    if (q.type == PT_DNA_NODE && q.active) dnaCount++;
                }
                int targetNode = (int)(p.jobProgress * dnaCount * 0.5f);
                int nodeIdx = 0;
                for (int j = 0; j < (int)phys.particles.size(); j++) {
                    auto& q = phys.particles[j];
                    if (q.type == PT_DNA_NODE && q.active) {
                        if (nodeIdx / 2 == targetNode) {
                            simd_float3 offset = {ipRandGauss() * 0.01f,
                                                  ipRandGauss() * 0.01f + 0.02f,
                                                  ipRandGauss() * 0.01f};
                            p.position = {q.position.x + offset.x,
                                          q.position.y + offset.y,
                                          q.position.z + offset.z};
                            break;
                        }
                        nodeIdx++;
                    }
                }
                if (p.jobProgress >= 1.0f) {
                    p.job = JOB_IDLE;
                    p.jobProgress = 0;
                    p.jobTimer = 0;
                }
            }
        }

        if (p.type == PT_VESICLE_COPII) {
            if (p.job == JOB_IDLE) {
                phys.assignGoto(i, randNear(golgiP, cellR * 0.05f), 0.012f);
            } else if (p.job == JOB_GOTO_TARGET && phys.reachedTarget(i, 0.06f)) {
                phys.assignConsume(i);
                p.spawnPos = randNear(erPos, cellR * 0.1f);
            }
        }

        if (p.type == PT_VESICLE_SECRETORY) {
            if (p.job == JOB_IDLE) {
                float a = ipRandf() * 2 * M_PI;
                simd_float3 memTarget = {cosf(a)*cellR*0.92f, 0, sinf(a)*cellR*0.92f};
                phys.assignGoto(i, memTarget, 0.01f);
            } else if (p.job == JOB_GOTO_TARGET) {
                float dx = p.position.x, dy = p.position.y, dz = p.position.z;
                float dist = sqrtf(dx*dx + dy*dy + dz*dz);
                if (dist > cellR * 0.88f) {
                    phys.assignConsume(i);
                    p.spawnPos = randNear(golgiP, cellR * 0.08f);
                }
            }
        }

        if (p.type == PT_CHAPERONE && p.job == JOB_IDLE) {
            int poly = findNearest(p.position, PT_POLYPEPTIDE, JOB_IDLE);
            if (poly >= 0) {
                phys.assignGoto(i, phys.particles[poly].position, 0.015f);
            }
        }

        if (p.type == PT_MOLECULE && (p.tag == TAG_NADH || p.tag == TAG_CALCIUM || p.tag == TAG_WATER)) {
            if (p.job == JOB_IDLE && p.jobTimer > 3.0f) {
                simd_float3 dest = randNear(p.tag == TAG_NADH ? mitoP[0] : localOrigin,
                                            p.tag == TAG_WATER ? cellR * 0.7f : cellR * 0.3f);
                phys.assignGoto(i, dest, 0.008f);
            } else if (p.job == JOB_GOTO_TARGET && phys.reachedTarget(i, 0.1f)) {
                p.job = JOB_IDLE;
                p.jobTimer = 0;
            }
        }
    }
}

// Build GPU data from physics particles
static void uploadCellInterior(simd_float3 cellPos, float cellSize, float time, float dt) {
    if (!gIntraInitialized) {
        initIntracellularParticles(cellPos, cellSize);
    }

    float simDt = fminf(fmaxf(dt, 1.0f / 240.0f), 1.0f / 30.0f);

    // ── LOCAL COORDINATE SYSTEM ─────────────────────────────────────
    // ALL particle positions are in CELL-LOCAL space (relative to cell center).
    // Cell center = origin (0,0,0) for particles.
    // NO position shifting, NO teleportation, NO 瞬移.
    // Membrane is a sphere of radius cellR centered at origin in local space.
    // World position = cellPos + localPosition (computed only for rendering).
    float cellR = cellSize;
    // Nucleus sizing anchored to 1 wu = 5 µm (MICROMETERS_PER_WORLD_UNIT).
    // HeLa nucleus is ~50% of cell diameter (Alberts MBoC Ch 4), so at
    // cellR = 2.0 wu (10 µm) the nucleus visual radius is 1.0 wu = 5 µm.
    // GLB scale in GLBLoader.h is picked to match this: nucleus.glb bbox
    // 9.40 × scale 0.106 / 2 ≈ 0.50 × cellR. DNA is confined 10% inside
    // the visual membrane so backbone sphere radii don't poke through.
    float nucVisualRadius = cellR * 0.50f;
    float nucR = cellR * 0.45f;

    // Physics params: membrane at origin, radius = cellR
    gIntraPhys.params.cellCenter = {0, 0, 0}; // LOCAL origin
    gIntraPhys.params.cellRadius = cellR;
    gIntraPhys.params.cellVelocity = {0, 0, 0}; // No cell velocity in local frame

    // Update confinement: all centered at local origin
    bool cleavageCompartmentsActive = mitosisCleavageCompartmentsActive();
    for (auto& p : gIntraPhys.particles) {
        if (!p.active) continue;
        p.confineCenter = {0, 0, 0}; // Local origin
        if (isNuclearParticleType(p.type)) {
            if (gMitosis.active && gMitosis.particlesDuplicated &&
                gMitosis.phase >= MITO_PROMETAPHASE && gMitosis.phase <= MITO_CYTOKINESIS) {
                float targetX = mitosisNuclearTargetX(cellR);
                float sign = (particleMitosisHalf(p) == 0) ? 1.0f : -1.0f;
                p.confineCenter = {sign * targetX, 0, 0};
                p.confineRadius = nucR * 0.86f;
            } else {
                p.confineRadius = nucR;
            }
        } else if (cleavageCompartmentsActive && isCleavageLockedParticle(p)) {
            if (p.mitosisHalf < 0) {
                p.mitosisHalf = (p.position.x >= 0.0f) ? 0 : 1;
            }
            p.confineCenter = mitosisCytoplasmCompartmentCenter(particleMitosisHalf(p), cellR);
            float compRadius = mitosisCytoplasmCompartmentRadius(cellR);
            if (p.type == PT_VESICLE_COPII || p.type == PT_VESICLE_SECRETORY) {
                compRadius *= 0.84f;
            } else if (p.type == PT_RIBOSOME_SMALL || p.type == PT_RIBOSOME_LARGE ||
                       p.type == PT_TRNA || p.type == PT_CHAPERONE) {
                compRadius *= 0.88f;
            }
            p.confineRadius = compRadius;
        } else {
            p.confineRadius = cellR;
        }
    }

    // Poles are computed inside gMitosis.update() — don't override them here

    gLastCellPos = cellPos;

    // ── JOB ASSIGNMENT LOGIC (every frame) ────────────────────────────
    // cellR already defined above
    int cellPhase = gSim.cells.empty() ? 0 : gSim.cells[0].phase;
    gCDogma.update(simDt, cellPhase);

    // ── Membrane-receptor signaling activity ─────────────────────
    // Each receptor pulses periodically to indicate ligand binding.
    // Occasionally a receptor "fires" (stronger flash) and spawns a
    // second-messenger particle that drifts toward the cell interior
    // — this is the visual of a real ligand → receptor → cascade event.
    //
    // Per-type fire rates reflect biological activity levels:
    //   GPCR        — frequent (β-adrenergic at rest ~0.5 Hz basal tone)
    //   RTK         — slower (growth factor pulses are rarer)
    //   Cytokine R  — episodic
    //   Ion channel — fast (gating kinetics ms timescale)
    //   Ion pump    — continuous (Na/K ATPase ~100 Hz)
    //   Adhesion    — static (structural, no firing)
    {
        float gTime = time;
        // per-type base pulse phase speed and fire probability per second
        struct PulseSpec { float baseGlow, pulseAmp, pulseHz, fireProb; };
        auto specFor = [](ParticleType t) -> PulseSpec {
            switch (t) {
                case PT_RECEPTOR_GPCR:     return {0.90f, 0.45f, 0.8f, 0.18f};
                case PT_RECEPTOR_RTK:      return {0.85f, 0.50f, 0.4f, 0.10f};
                case PT_RECEPTOR_CYTOKINE: return {0.80f, 0.50f, 0.5f, 0.07f};
                case PT_RECEPTOR_STK:      return {0.70f, 0.40f, 0.3f, 0.04f};
                case PT_RECEPTOR_DEATH:    return {0.80f, 0.35f, 0.6f, 0.01f};
                case PT_RECEPTOR_TLR:      return {0.70f, 0.45f, 0.5f, 0.03f};
                case PT_RECEPTOR_NOTCH:    return {0.70f, 0.30f, 0.3f, 0.02f};
                case PT_RECEPTOR_FRIZZLED: return {0.70f, 0.40f, 0.4f, 0.04f};
                case PT_RECEPTOR_PATCHED:  return {0.70f, 0.40f, 0.4f, 0.03f};
                case PT_LDL_RECEPTOR:      return {0.70f, 0.30f, 0.5f, 0.10f};
                case PT_TRANSFERRIN_R:     return {0.70f, 0.30f, 0.5f, 0.12f};
                case PT_ION_CHANNEL:       return {0.80f, 0.55f, 2.0f, 0.30f};
                case PT_ION_PUMP:          return {0.95f, 0.25f, 1.5f, 0.20f};
                case PT_EXCHANGER:         return {0.85f, 0.30f, 1.2f, 0.18f};
                case PT_TRANSPORTER:       return {0.80f, 0.30f, 0.8f, 0.15f};
                case PT_AQUAPORIN:         return {0.70f, 0.25f, 3.0f, 0.00f};
                case PT_GAP_JUNCTION:      return {0.80f, 0.20f, 0.2f, 0.05f};
                case PT_MHC:               return {0.70f, 0.15f, 0.2f, 0.00f};
                case PT_PHOSPHOLIPID:      return {0.30f, 0.10f, 0.3f, 0.00f};
                case PT_ADHESION:          return {0.75f, 0.25f, 0.2f, 0.02f};
                default: return {0.0f, 0.0f, 0.0f, 0.0f};
            }
        };

        auto isMembraneParticle = [](ParticleType t) {
            return t == PT_RECEPTOR_GPCR || t == PT_RECEPTOR_RTK ||
                   t == PT_RECEPTOR_CYTOKINE || t == PT_RECEPTOR_STK ||
                   t == PT_RECEPTOR_DEATH || t == PT_RECEPTOR_TLR ||
                   t == PT_RECEPTOR_NOTCH || t == PT_RECEPTOR_FRIZZLED ||
                   t == PT_RECEPTOR_PATCHED || t == PT_LDL_RECEPTOR ||
                   t == PT_TRANSFERRIN_R || t == PT_ION_CHANNEL ||
                   t == PT_ION_PUMP || t == PT_EXCHANGER ||
                   t == PT_TRANSPORTER || t == PT_AQUAPORIN ||
                   t == PT_GAP_JUNCTION || t == PT_MHC ||
                   t == PT_PHOSPHOLIPID || t == PT_ADHESION;
        };

        // Second-messenger spawn cap so firing events don't flood the sim.
        // Each fire event emits a short-lived cAMP/Ca²⁺/IP3 particle that
        // drifts from the membrane toward the cell interior.
        static int sMemMessengerCount = 0;
        constexpr int MAX_MEM_MESSENGERS = 800;

        for (int i = 0; i < (int)gIntraPhys.particles.size(); i++) {
            auto& p = gIntraPhys.particles[i];
            if (!p.active) continue;
            if (!isMembraneParticle(p.type)) continue;

            PulseSpec sp = specFor(p.type);
            // Background pulse — each particle has its own phase so the
            // membrane shimmers organically instead of blinking in sync.
            float phase = (float)i * 0.37f + p.position.x * 1.3f + p.position.z * 1.7f;
            float breathe = sinf(gTime * sp.pulseHz * 2.0f * M_PI + phase) * 0.5f + 0.5f;
            float glow = sp.baseGlow + sp.pulseAmp * breathe;

            // Stochastic fire event — receptor "binds" a ligand and flashes.
            // Decay of glow from fire event tracked in jobTimer.
            bool fired = false;
            if (sp.fireProb > 0.0f) {
                float firePerFrame = sp.fireProb * simDt;
                if (ipRandf() < firePerFrame) {
                    p.jobTimer = 1.0f; // peak flash
                    fired = true;
                }
            }
            if (p.jobTimer > 0.01f) {
                glow += p.jobTimer * 1.2f;          // bright flash
                p.jobTimer *= (1.0f - simDt * 3.0f); // decay over ~0.3 s
            }
            p.glowIntensity = glow;

            // Chemistry — on a fire event, the receptor/channel/pump
            // triggers a downstream effect: spawn a second-messenger bead
            // drifting inward. Messenger type depends on the membrane
            // protein's pathway:
            //   GPCR (Gs family) → cAMP
            //   GPCR (Gq) / ion channel / exchanger → Ca²⁺
            //   RTK → GF signal (sent to MAPK cascade; we emit Ca²⁺ proxy)
            //   Cytokine R → cAMP proxy (JAK activation mark)
            //   TLR/Death R → Ca²⁺ (inflammation / apoptosis spike)
            //   Ion pump → ATP consumption event (burst inward)
            //   Transporter → glucose / AA inward drift
            //   LDL/Transferrin R → endocytic uptake (no inward bead, just flash)
            // Hard-cap total particle count. Also disable messenger
            // spawning until we ship a reusable pool — the unbounded
            // vector growth from ~8 spawns/frame reaches 30k+ particles
            // within seconds and crashes at ~1.5M GPU atoms.
            constexpr int MAX_TOTAL_PARTICLES = 30000;
            constexpr bool ENABLE_MESSENGER_SPAWN = false;
            if (ENABLE_MESSENGER_SPAWN && fired &&
                sMemMessengerCount < MAX_MEM_MESSENGERS &&
                (int)gIntraPhys.particles.size() < MAX_TOTAL_PARTICLES) {
                MolTag tag = TAG_CAMP_SIGNAL;
                float cr = 0.2f, cg = 0.9f, cb = 1.0f;
                float rad = cellR * 0.003f;
                switch (p.type) {
                    case PT_RECEPTOR_GPCR:      tag = TAG_CAMP_SIGNAL;  cr=0.2f; cg=0.7f; cb=1.0f; break;
                    case PT_RECEPTOR_RTK:       tag = TAG_CALCIUM;      cr=1.0f; cg=0.7f; cb=0.2f; break;
                    case PT_RECEPTOR_CYTOKINE:  tag = TAG_CAMP_SIGNAL;  cr=1.0f; cg=0.4f; cb=0.8f; break;
                    case PT_RECEPTOR_TLR:       tag = TAG_CALCIUM;      cr=0.3f; cg=0.4f; cb=1.0f; break;
                    case PT_RECEPTOR_DEATH:     tag = TAG_CALCIUM;      cr=1.0f; cg=0.2f; cb=0.2f; break;
                    case PT_ION_CHANNEL:        tag = TAG_CALCIUM;      cr=0.3f; cg=1.0f; cb=0.6f; break;
                    case PT_EXCHANGER:          tag = TAG_CALCIUM;      cr=1.0f; cg=0.6f; cb=0.2f; break;
                    case PT_ION_PUMP:           tag = TAG_ATP;          cr=0.3f; cg=1.0f; cb=0.3f; break;
                    case PT_TRANSPORTER:        tag = TAG_GLUCOSE;      cr=0.9f; cg=0.8f; cb=0.2f; break;
                    default: break;
                }
                simd_float3 inward = { -p.position.x * 0.15f,
                                       -p.position.y * 0.15f,
                                       -p.position.z * 0.15f };
                simd_float3 spawnPos = { p.position.x * 0.92f,
                                         p.position.y * 0.92f,
                                         p.position.z * 0.92f };
                int mid = gIntraPhys.addParticle(PT_MOLECULE, spawnPos, rad, cr, cg, cb);
                if (mid >= 0 && mid < (int)gIntraPhys.particles.size()) {
                    auto& m = gIntraPhys.particles[mid];
                    m.tag = tag;
                    m.velocity = inward;
                    m.spawnPos = spawnPos;
                    m.jobTimer = 0.0f;
                    m.glowIntensity = 1.4f; // bright burst
                    simd_float3 zero = {0, 0, 0};
                    gIntraPhys.setConfinement(mid, zero, cellR);
                    sMemMessengerCount++;
                }
            }
        }
        // Soft decay of the messenger cap so over time new spawns are
        // allowed as old ones are consumed by metabolism/physics.
        sMemMessengerCount = (int)(sMemMessengerCount * 0.985f);
    }
    // Mirror the primary cell's live gCDogma state back to its SimCell so the
    // simulation's per-cell checkpoint code (updateCellCycle in Simulation.h)
    // can gate G1→S and G2→M decisions on the same replication progress
    // shown in the UI. Without this mirror the primary's program.cdogma
    // would lag by however long since the last focus change.
    if (!gSim.cells.empty()) {
        gSim.cells[0].program.cdogma = gCDogma;
        gSim.cells[0].program.cdogmaInitialized = true;
    }
    syncGenomeReplicationVisuals(gIntraPhys, gCDogma, cellR, gMitosis.active);

    // Key positions for job targets — all LOCAL space (origin = cell center)
    simd_float3 localOrigin = {0, 0, 0};
    simd_float3 nucPos = localOrigin;
    simd_float3 erPos = {0, 0, cellR * 0.35f};
    simd_float3 golgiP = {0, 0, -cellR * 0.42f};
    simd_float3 mitoP[3] = {
        {cellR*0.3f, cellR*0.12f, cellR*0.18f},
        {-cellR*0.24f, -cellR*0.12f, 0},
        {-cellR*0.18f, 0, -cellR*0.24f}
    };
    if (!gMitosis.active && !gSim.cells.empty()) {
        OrganelleMotionState& motion = updateOrganelleMotionState(gSim.cells[0].cellUid, time);
        golgiP = {motion.golgiPos.x * cellR, motion.golgiPos.y * cellR, motion.golgiPos.z * cellR};
        for (int m = 0; m < 3; m++) {
            mitoP[m] = {motion.mitoPos[m].x * cellR,
                        motion.mitoPos[m].y * cellR,
                        motion.mitoPos[m].z * cellR};
        }
    }

    // Helper: random point near a position
    auto randNear = [](simd_float3 c, float r) -> simd_float3 {
        float u = ipRandf()*2-1, phi = ipRandf()*2*M_PI;
        float rad = cbrtf(ipRandf()) * r;
        float st = sqrtf(1-u*u);
        return {c.x + rad*st*cosf(phi), c.y + rad*u, c.z + rad*st*sinf(phi)};
    };

    auto anchoredPointForParticle = [&](const IntraParticle& particle,
                                        simd_float3 anchor,
                                        float keepLocal = 0.72f) -> simd_float3 {
        if (!cleavageCompartmentsActive || !isCleavageLockedParticle(particle)) return anchor;
        return partitionPointIntoDaughterHalf(anchor, particleMitosisHalf(particle), cellR, keepLocal);
    };

    auto randNearParticleAnchor = [&](const IntraParticle& particle,
                                      simd_float3 anchor,
                                      float spread,
                                      float keepLocal = 0.72f) -> simd_float3 {
        simd_float3 base = anchoredPointForParticle(particle, anchor, keepLocal);
        float compartmentSpread = spread;
        if (cleavageCompartmentsActive && isCleavageLockedParticle(particle)) {
            compartmentSpread = fminf(spread, mitosisCytoplasmCompartmentRadius(cellR) * 0.42f);
        }
        return randNear(base, compartmentSpread);
    };

    // Keep particle traffic aligned with the visible large organelles.
    gIntraPhys.attractionFields.clear();
    if (cleavageCompartmentsActive) {
        for (int half = 0; half < 2; half++) {
            for (int m = 0; m < 3; m++) {
                simd_float3 mitoAnchor = partitionPointIntoDaughterHalf(mitoP[m], half, cellR, 0.68f);
                gIntraPhys.attractionFields.push_back({mitoAnchor, cellR*0.26f, 0.015f, TAG_NADH, PT_COUNT});
                gIntraPhys.attractionFields.push_back({mitoAnchor, cellR*0.22f, 0.008f, TAG_GLUCOSE, PT_COUNT});
            }
            simd_float3 erAnchor = partitionPointIntoDaughterHalf(erPos, half, cellR, 0.74f);
            simd_float3 golgiAnchor = partitionPointIntoDaughterHalf(golgiP, half, cellR, 0.70f);
            gIntraPhys.attractionFields.push_back({erAnchor, cellR*0.28f, 0.012f, TAG_AMINO_ACID, PT_COUNT});
            gIntraPhys.attractionFields.push_back({erAnchor, cellR*0.30f, 0.010f, TAG_NONE, PT_TRNA});
            gIntraPhys.attractionFields.push_back({erAnchor, cellR*0.24f, 0.010f, TAG_CALCIUM, PT_COUNT});
            gIntraPhys.attractionFields.push_back({erAnchor, cellR*0.24f, 0.008f, TAG_NONE, PT_CHAPERONE});
            gIntraPhys.attractionFields.push_back({golgiAnchor, cellR*0.24f, 0.012f, TAG_NONE, PT_VESICLE_COPII});
        }
        gIntraPhys.attractionFields.push_back({nucPos, nucR*1.2f, 0.020f, TAG_NONE, PT_RNA_POLYMERASE});
    } else {
        for (int m = 0; m < 3; m++) {
            gIntraPhys.attractionFields.push_back({mitoP[m], cellR*0.4f, 0.015f, TAG_NADH, PT_COUNT});
            gIntraPhys.attractionFields.push_back({mitoP[m], cellR*0.3f, 0.008f, TAG_GLUCOSE, PT_COUNT});
        }
        gIntraPhys.attractionFields.push_back({erPos, cellR*0.35f, 0.012f, TAG_AMINO_ACID, PT_COUNT});
        gIntraPhys.attractionFields.push_back({erPos, cellR*0.4f, 0.010f, TAG_NONE, PT_TRNA});
        gIntraPhys.attractionFields.push_back({nucPos, nucR*1.2f, 0.020f, TAG_NONE, PT_RNA_POLYMERASE});
        gIntraPhys.attractionFields.push_back({erPos, cellR*0.3f, 0.010f, TAG_CALCIUM, PT_COUNT});
        gIntraPhys.attractionFields.push_back({erPos, cellR*0.3f, 0.008f, TAG_NONE, PT_CHAPERONE});
        gIntraPhys.attractionFields.push_back({golgiP, cellR*0.3f, 0.012f, TAG_NONE, PT_VESICLE_COPII});
    }

    // Helper: find nearest particle of type with job
    auto findNearest = [&](simd_float3 from, ParticleType type, MolJob job) -> int {
        int best = -1; float bestD = 1e9f;
        for (int i = 0; i < (int)gIntraPhys.particles.size(); i++) {
            auto& q = gIntraPhys.particles[i];
            if (!q.active || q.type != type || q.job != job) continue;
            float dx = q.position.x-from.x, dy = q.position.y-from.y, dz = q.position.z-from.z;
            float d = dx*dx+dy*dy+dz*dz;
            if (d < bestD) { bestD = d; best = i; }
        }
        return best;
    };

    std::vector<int> primaryRNAPolParticles;
    std::vector<int> primaryDNAPolParticles;
    std::vector<int> primaryRibosomeParticles;
    std::vector<int> primaryTRNAParticles;
    primaryRNAPolParticles.reserve(3);
    primaryDNAPolParticles.reserve(4);
    primaryRibosomeParticles.reserve(6);
    primaryTRNAParticles.reserve(CDOGMA_MAX_TRNA_POOL);
    for (int i = 0; i < (int)gIntraPhys.particles.size(); i++) {
        const auto& p = gIntraPhys.particles[i];
        if (!p.active) continue;
        if (p.type == PT_RNA_POLYMERASE && (int)primaryRNAPolParticles.size() < 3) {
            primaryRNAPolParticles.push_back(i);
        } else if (p.type == PT_DNA_POLYMERASE && (int)primaryDNAPolParticles.size() < 4) {
            primaryDNAPolParticles.push_back(i);
        } else if (p.type == PT_RIBOSOME_LARGE && (int)primaryRibosomeParticles.size() < 6) {
            primaryRibosomeParticles.push_back(i);
        } else if (p.type == PT_TRNA && (int)primaryTRNAParticles.size() < CDOGMA_MAX_TRNA_POOL) {
            primaryTRNAParticles.push_back(i);
        }
    }

    // ── mRNA bead chain per polymerase ───────────────────────────
    // Every ~60 bp of transcription, spawn a PT_PRE_MRNA bead at the
    // polymerase position and tether it to the previously-spawned bead,
    // so the growing transcript looks like a chain trailing behind the
    // RNA Pol II. The last exon bead in an active transcript becomes
    // PT_MRNA_NODE (mature) when splicing flags it done.
    // Cap per-transcript bead count to keep particle counts sane.
    static std::map<int, int> sLastSpawnedBP;    // slot → last bp we spawned at
    static std::map<int, int> sPrevBeadIdx;      // slot → last bead particle idx
    static std::map<int, int> sBeadCount;        // slot → beads alive for this tx
    constexpr int MRNA_BP_PER_BEAD = 60;
    constexpr int MRNA_MAX_BEADS_PER_TX = 24;

    for (int slot = 0; slot < (int)primaryRNAPolParticles.size(); slot++) {
        auto& p = gIntraPhys.particles[primaryRNAPolParticles[slot]];
        const auto& ts = gCDogma.transcription[slot];
        if (ts.active) {
            p.job = JOB_FOLLOW_PATH;
            p.jobProgress = ts.rnaPolPosition;
            p.jobSpeed = 0.0f;
            p.jobTimer = 0.0f;
            p.glowIntensity = 1.3f;

            // Emit a bead when transcription has advanced far enough.
            int lastBP = sLastSpawnedBP.count(slot) ? sLastSpawnedBP[slot] : -MRNA_BP_PER_BEAD;
            int count = sBeadCount.count(slot) ? sBeadCount[slot] : 0;
            if (ts.currentBP - lastBP >= MRNA_BP_PER_BEAD && count < MRNA_MAX_BEADS_PER_TX) {
                // Color intron beads dim / desaturated, exon beads bright.
                bool inExon = isExon(ts.currentBP);
                float r = inExon ? 0.85f : 0.4f;
                float g = inExon ? 0.35f : 0.25f;
                float b = inExon ? 0.95f : 0.45f;
                simd_float3 beadPos = {
                    p.position.x + (ipRandf() - 0.5f) * 0.02f,
                    p.position.y + (ipRandf() - 0.5f) * 0.02f,
                    p.position.z + (ipRandf() - 0.5f) * 0.02f
                };
                // Bead radius: single RNA nucleotide ~0.7 nm → 0.00035 µm
                // radius. Boost like other macromolecule-adjacent beads.
                float beadR = molRadiusWu(MoleculeRadiusUm::TRNA, MOL_VIS_BOOST_MACROMOL) * 0.6f;
                int idx = gIntraPhys.addParticle(PT_PRE_MRNA, beadPos, beadR, r, g, b);
                simd_float3 origin = {0, 0, 0};
                gIntraPhys.setConfinement(idx, origin, cellR);
                // Tether to previous bead of this transcript for chain feel.
                int prev = sPrevBeadIdx.count(slot) ? sPrevBeadIdx[slot] : -1;
                if (prev >= 0 && prev < (int)gIntraPhys.particles.size()) {
                    gIntraPhys.setTether(idx, prev, 0.05f, 2.0f);
                }
                gIntraPhys.particles[idx].stateIndex = ts.currentBP;
                gIntraPhys.particles[idx].mitosisHalf = -1;
                sPrevBeadIdx[slot] = idx;
                sLastSpawnedBP[slot] = ts.currentBP;
                sBeadCount[slot] = count + 1;
            }
        } else {
            p.job = JOB_IDLE;
            p.jobProgress = 0.0f;
            p.glowIntensity = 0.6f;
            // Transcription ended — reset the per-slot spawn trackers so
            // a future transcription event starts a fresh bead chain.
            if (sLastSpawnedBP.count(slot)) {
                sLastSpawnedBP[slot] = -MRNA_BP_PER_BEAD;
                sPrevBeadIdx[slot] = -1;
                sBeadCount[slot] = 0;
            }
        }
    }

    // Promote mature mRNA beads to PT_MRNA_NODE (post-splicing) so their
    // color flips from pre-mRNA dim to mRNA bright. This is a cosmetic
    // mapping — the real splicing state lives in gCDogma.splicing[].
    for (auto& q : gIntraPhys.particles) {
        if (!q.active || q.type != PT_PRE_MRNA) continue;
        // Once the transcript this bead came from has any completed
        // splicing event, mark as mature. Simple proxy: check bead age.
        if (q.jobTimer > 3.0f && isExon(q.stateIndex)) {
            q.type = PT_MRNA_NODE;
            q.colorG = fminf(1.0f, q.colorG + 0.30f);
        }
    }

    auto findDNAWorldHome = [&](int bpIndex, int strandPick) -> simd_float3 {
        simd_float3 fallback = {0, 0, 0};
        int seen = 0;
        for (const auto& q : gIntraPhys.particles) {
            if (!q.active || q.type != PT_DNA_NODE || q.stateIndex != bpIndex) continue;
            if (seen == strandPick) return q.position;
            seen++;
            fallback = q.position;
        }
        return fallback;
    };

    if (!gMitosis.active) {
        for (int slot = 0; slot < (int)primaryDNAPolParticles.size(); slot++) {
            auto& p = gIntraPhys.particles[primaryDNAPolParticles[slot]];
            const auto* fork = gCDogma.getReplicationFork(slot);
            if (fork) {
                int targetBase = gCDogma.replicationForkBase(slot);
                int strandPick = (fork->direction > 0) ? 0 : 1;
                simd_float3 dnaPos = findDNAWorldHome(targetBase, strandPick);
                simd_float3 offset = {
                    0.0f,
                    cellR * 0.015f * ((slot & 1) ? 1.0f : -1.0f),
                    cellR * 0.012f * (fork->direction < 0 ? -1.0f : 1.0f)
                };
                simd_float3 target = {dnaPos.x + offset.x, dnaPos.y + offset.y, dnaPos.z + offset.z};
                p.job = JOB_FOLLOW_PATH;
                p.jobProgress = gCDogma.replicationProgress;
                p.jobSpeed = 0.0f;
                p.jobTimer = 0.0f;
                p.home = target;
                p.position = {
                    p.position.x + (target.x - p.position.x) * 0.30f,
                    p.position.y + (target.y - p.position.y) * 0.30f,
                    p.position.z + (target.z - p.position.z) * 0.30f
                };
                p.velocity = {p.velocity.x * 0.25f, p.velocity.y * 0.25f, p.velocity.z * 0.25f};
                if (fork->proofreading) {
                    p.glowIntensity = 1.65f;
                    p.position.y += cellR * 0.006f * sinf(time * 18.0f + slot);
                } else {
                    p.glowIntensity = 1.00f + 0.35f * fork->recruitment;
                }
            } else {
                p.job = JOB_IDLE;
                p.jobProgress = 0.0f;
                p.glowIntensity = 0.45f + 0.20f * gCDogma.polymeraseRecruitment;
            }
        }
    }

    for (int slot = 0; slot < (int)primaryRibosomeParticles.size(); slot++) {
        auto& p = gIntraPhys.particles[primaryRibosomeParticles[slot]];
        const auto& tr = gCDogma.translation[slot];
        if (tr.active) {
            p.job = JOB_FOLLOW_PATH;
            p.jobProgress = tr.ribosomePosition;
            p.jobTimer = (float)tr.polypeptideLength / 24.0f;
            p.glowIntensity = 0.9f + 0.45f * tr.codonProgress;
        } else {
            p.job = JOB_IDLE;
            p.jobProgress = 0.0f;
            p.glowIntensity = 0.6f;
        }
    }

    int drivenTRNASlot = 0;
    for (const auto& cargo : gCDogma.trnaPool) {
        if (!cargo.active || drivenTRNASlot >= (int)primaryTRNAParticles.size()) continue;
        int particleIdx = primaryTRNAParticles[drivenTRNASlot++];
        auto& p = gIntraPhys.particles[particleIdx];
        p.glowIntensity = 0.9f + 0.7f * cargo.shuttle;
        if (cargo.ribosomeIndex >= 0 && cargo.ribosomeIndex < (int)primaryRibosomeParticles.size()) {
            simd_float3 ribPos = gIntraPhys.particles[primaryRibosomeParticles[cargo.ribosomeIndex]].position;
            gIntraPhys.assignGoto(particleIdx, ribPos, 0.04f);
        }
    }

    // When particle count is elevated (S-phase/mitosis duplication), skip the
    // per-molecule job assignment loop. It calls findNearest() which is O(n)
    // per call × O(n) molecules = O(n²) total → freezes at 1800+ particles.
    // Existing jobs keep running; only new job assignments are deferred.
    bool skipJobAssignment = ((int)gIntraPhys.particles.size() > 1200);

    for (int i = 0; i < (int)gIntraPhys.particles.size(); i++) {
        auto& p = gIntraPhys.particles[i];
        if (!p.active) continue;
        if (skipJobAssignment && p.type == PT_MOLECULE) continue;

        // ── ATP: mito → random organelle → consumed → respawn at mito ──
        if (p.type == PT_MOLECULE && p.tag == TAG_ATP) {
            if (p.job == JOB_IDLE) {
                // Pick a random destination (ribosome, ER, or random cytoplasm point)
                simd_float3 dest;
                float r = ipRandf();
                if (r < 0.4f) {
                    int rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_IDLE);
                    dest = (rib >= 0)
                        ? gIntraPhys.particles[rib].position
                        : randNearParticleAnchor(p, erPos, cellR * 0.3f, 0.74f);
                } else if (r < 0.7f) {
                    dest = randNearParticleAnchor(p, erPos, cellR * 0.2f, 0.74f);
                } else {
                    dest = randNearParticleAnchor(p, localOrigin, cellR * 0.6f, 0.62f);
                }
                gIntraPhys.assignGoto(i, dest, 0.025f);
            } else if (p.job == JOB_GOTO_TARGET && gIntraPhys.reachedTarget(i, 0.08f)) {
                // Arrived → flash and consume
                p.glowIntensity = 2.0f; // bright flash
                gIntraPhys.assignConsume(i);
                p.spawnPos = anchoredPointForParticle(
                    p, mitoP[(int)(ipRandf() * 2.99f)], 0.68f); // respawn at random mito
            }
        }

        // ── Glucose: membrane → center → consumed → respawn at edge ────
        if (p.type == PT_MOLECULE && p.tag == TAG_GLUCOSE) {
            if (p.job == JOB_IDLE) {
                gIntraPhys.assignGoto(i,
                                      randNearParticleAnchor(p, localOrigin, cellR * 0.3f, 0.62f),
                                      0.015f);
            } else if (p.job == JOB_GOTO_TARGET && gIntraPhys.reachedTarget(i, 0.08f)) {
                gIntraPhys.assignConsume(i);
                // Respawn at cell membrane edge
                float a = ipRandf() * 2 * M_PI;
                simd_float3 membranePoint = {cosf(a)*cellR*0.88f, 0, sinf(a)*cellR*0.88f};
                p.spawnPos = anchoredPointForParticle(p, membranePoint, 0.96f);
            }
        }

        // ── Amino acids: float near ER → go to active ribosome ─────────
        if (p.type == PT_MOLECULE && p.tag == TAG_AMINO_ACID) {
            if (p.job == JOB_IDLE) {
                // Find an active ribosome (translating)
                int rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_FOLLOW_PATH);
                if (rib >= 0) {
                    gIntraPhys.assignGoto(i, gIntraPhys.particles[rib].position, 0.02f);
                } else {
                    // Just drift near ER
                    gIntraPhys.assignGoto(i,
                                          randNearParticleAnchor(p, erPos, cellR * 0.25f, 0.74f),
                                          0.008f);
                }
            } else if (p.job == JOB_GOTO_TARGET && gIntraPhys.reachedTarget(i, 0.06f)) {
                gIntraPhys.assignConsume(i);
                p.spawnPos = randNearParticleAnchor(p, erPos, cellR * 0.25f, 0.74f);
            }
        }

        // ── tRNA: pick up amino acid → deliver to ribosome → release ───
        if (p.type == PT_TRNA) {
            if (p.job == JOB_IDLE) {
                // Find nearest active ribosome
                int rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_FOLLOW_PATH);
                if (rib < 0) rib = findNearest(p.position, PT_RIBOSOME_LARGE, JOB_IDLE);
                if (rib >= 0) {
                    gIntraPhys.assignGoto(i, gIntraPhys.particles[rib].position, 0.03f);
                } else {
                    gIntraPhys.assignGoto(i,
                                          randNearParticleAnchor(p, erPos, cellR * 0.3f, 0.74f),
                                          0.01f);
                }
            } else if (p.job == JOB_GOTO_TARGET && gIntraPhys.reachedTarget(i, 0.07f)) {
                // Delivered, go idle briefly then find another ribosome
                p.job = JOB_IDLE;
                p.jobTimer = 0;
            }
        }

        // ── RNA Pol II: follow DNA backbone ────────────────────────────
        if (p.type == PT_RNA_POLYMERASE) {
            bool dogmaControlled = false;
            for (int slot = 0; slot < (int)primaryRNAPolParticles.size(); slot++) {
                if (primaryRNAPolParticles[slot] == i) {
                    dogmaControlled = true;
                    if (!gCDogma.transcription[slot].active) {
                        p.job = JOB_IDLE;
                        p.jobProgress = 0.0f;
                        p.jobTimer = 0.0f;
                        break;
                    }
                }
            }
            if (!dogmaControlled && p.job == JOB_IDLE && p.jobTimer > 2.0f) {
                gIntraPhys.assignFollowPath(i, 0.008f);
            }
            if (p.job == JOB_FOLLOW_PATH) {
                if (!dogmaControlled) p.jobProgress += simDt * p.jobSpeed;
                // Find corresponding DNA node
                int dnaCount = 0;
                for (auto& q : gIntraPhys.particles)
                    if (q.type == PT_DNA_NODE && q.active) dnaCount++;
                int targetNode = (int)(p.jobProgress * dnaCount * 0.5f); // strand A only
                int nodeIdx = 0;
                for (int j = 0; j < (int)gIntraPhys.particles.size(); j++) {
                    auto& q = gIntraPhys.particles[j];
                    if (q.type == PT_DNA_NODE && q.active) {
                        if (nodeIdx / 2 == targetNode) {
                            // Move RNA Pol to this DNA node
                            simd_float3 offset = {ipRandGauss() * 0.01f, ipRandGauss() * 0.01f + 0.02f, ipRandGauss() * 0.01f};
                            p.position = {q.position.x + offset.x, q.position.y + offset.y, q.position.z + offset.z};
                            break;
                        }
                        nodeIdx++;
                    }
                }
                if (p.jobProgress >= 1.0f) {
                    p.job = JOB_IDLE; p.jobProgress = 0; p.jobTimer = 0;
                }
            }
        }

        // ── COPII vesicle: ER → Golgi → recycle ───────────────────────
        if (p.type == PT_VESICLE_COPII) {
            if (p.job == JOB_IDLE) {
                gIntraPhys.assignGoto(i,
                                      randNearParticleAnchor(p, golgiP, cellR * 0.05f, 0.70f),
                                      0.012f);
            } else if (p.job == JOB_GOTO_TARGET && gIntraPhys.reachedTarget(i, 0.06f)) {
                // Arrived at Golgi → consume and respawn at ER
                gIntraPhys.assignConsume(i);
                p.spawnPos = randNearParticleAnchor(p, erPos, cellR * 0.1f, 0.74f);
            }
        }

        // ── Secretory vesicle: Golgi → membrane → exocytosis → recycle ─
        if (p.type == PT_VESICLE_SECRETORY) {
            if (p.job == JOB_IDLE) {
                float a = ipRandf() * 2 * M_PI;
                simd_float3 memTarget = {cosf(a)*cellR*0.92f, 0, sinf(a)*cellR*0.92f};
                memTarget = anchoredPointForParticle(p, memTarget, 0.96f);
                gIntraPhys.assignGoto(i, memTarget, 0.01f);
            } else if (p.job == JOB_GOTO_TARGET) {
                float dx = p.position.x, dy = p.position.y, dz = p.position.z; // local space, origin=0
                float dist = sqrtf(dx*dx+dy*dy+dz*dz);
                if (dist > cellR * 0.88f) {
                    // Near membrane → exocytosis: shrink then respawn at Golgi
                    gIntraPhys.assignConsume(i);
                    p.spawnPos = randNearParticleAnchor(p, golgiP, cellR * 0.08f, 0.70f);
                }
            }
        }

        // ── Chaperone: drift near ER, attracted to polypeptide ─────────
        if (p.type == PT_CHAPERONE && p.job == JOB_IDLE) {
            int poly = findNearest(p.position, PT_POLYPEPTIDE, JOB_IDLE);
            if (poly >= 0) {
                gIntraPhys.assignGoto(i, gIntraPhys.particles[poly].position, 0.015f);
            }
        }

        // ── NADH/Ca/water: gentle directed motion ──────────────────────
        if (p.type == PT_MOLECULE && (p.tag == TAG_NADH || p.tag == TAG_CALCIUM || p.tag == TAG_WATER)) {
            if (p.job == JOB_IDLE && p.jobTimer > 3.0f) {
                simd_float3 baseAnchor = (p.tag == TAG_NADH) ? mitoP[0] : localOrigin;
                simd_float3 dest = randNearParticleAnchor(
                    p, baseAnchor,
                    p.tag == TAG_WATER ? cellR * 0.7f : cellR * 0.3f,
                    p.tag == TAG_NADH ? 0.68f : 0.62f);
                gIntraPhys.assignGoto(i, dest, 0.008f);
            } else if (p.job == JOB_GOTO_TARGET && gIntraPhys.reachedTarget(i, 0.1f)) {
                p.job = JOB_IDLE; p.jobTimer = 0;
            }
        }
    }

    // Finish any remaining division prep during late S/G2 so visible mitosis
    // starts from a fully duplicated, readable genome state.
    // Guard: only run ONCE. Without this guard, if interiorPreparedForPhysicalDivision
    // returns false after the first call (e.g. aux particles not yet tagged), it
    // would re-duplicate nuclear particles every frame → exponential particle growth.
    {
        // Per-interior prep is already idempotent inside
        // prepareInteriorForPhysicalDivision(), so we do not need a global
        // sticky flag here. The old global guard broke promoted second-lineage
        // cells: a new primary cell could inherit "prep already done" from the
        // previous primary and then sit forever at DNA=1524 / phase=3 with no
        // spindle animation because its own interior never got nuclear-aux
        // duplication.
        if (!gMitosis.active && gCDogma.replicationProgress >= 0.995f) {
            prepareInteriorForPhysicalDivision(gIntraPhys, cellR);
        }
    }

    // ── MITOSIS STATE MANAGEMENT ────────────────────────────────────
    // Only start mitosis if cell is in M-phase, no active mitosis, AND no cooldown
    bool hasCooldown = (!gSim.cells.empty() && gSim.cells[0].divisionCooldown > 0);
    bool interiorDivisionReady = interiorPreparedForPhysicalDivision(gIntraPhys);
    if (cellPhase == 3 && !gMitosis.active && !hasCooldown && gIntraInitialized &&
        gCDogma.replicationReadyForM()) {
        if (!interiorDivisionReady) {
            prepareInteriorForPhysicalDivision(gIntraPhys, cellR);
            interiorDivisionReady = interiorPreparedForPhysicalDivision(gIntraPhys);
        }
        gMitosis.start(cellPos, cellR);
        bool genomePreReplicated = interiorGenomeAlreadyReplicated(gIntraPhys);
        gMitosis.particlesDuplicated = interiorDivisionReady;
        gMitosis.dnaCheckpointPassed = genomePreReplicated && gCDogma.replicationReadyForM();
        gMitosis.replicationQuality = gCDogma.replicationQuality;
        gMitosis.survivingDaughter = 0; // Cell[0] always tracks daughter A (primary)
        DNADistribution dnaStart = countActiveDNADistribution();
        if (genomePreReplicated) {
            sMitosisExpectedDNACount = dnaStart.total > 0 ? dnaStart.total : 1;
            sMitosisDNAStartCount = std::max(sMitosisExpectedDNACount / 2, 1);
        } else {
            sMitosisDNAStartCount = dnaStart.total > 0 ? dnaStart.total : 1;
            sMitosisExpectedDNACount = sMitosisDNAStartCount * 2;
        }
        // Start fresh log
        if (sMitosisLog) { fclose(sMitosisLog); sMitosisLog = nullptr; }
        sMitosisLogPhase = -1;
        mitosisLog("═══════════════════════════════════════════════\n");
        mitosisLog("MITOSIS STARTED — cells=%d cellR=%.3f\n", (int)gSim.cells.size(), cellR);
        mitosisLog("DNA checkpoint baseline=%d expected=%d\n",
            sMitosisDNAStartCount, sMitosisExpectedDNACount);
        mitosisLog("═══════════════════════════════════════════════\n");
        mitosisLogSnapshot("start", cellR);
    }
    // When nuclear particles are duplicated (~1029→~1800+), the
    // particle-particle collision detection becomes O(n²) on the
    // clustered nucleus and freezes the app.
    // Skip collisions whenever particle count exceeds the normal
    // interphase level (~1100). All other physics still runs.
    gIntraPhys.skipCollisions = ((int)gIntraPhys.particles.size() > 1200);

    if (gMitosis.active) {
        gMitosis.update(simDt, cellPos, cellR);

        // Log phase transitions
        if (gMitosis.phase != sMitosisLogPhase) {
            const char* phaseNames[] = {"NONE","PROPHASE","PROMETA","META","ANA","TELO","CYTO","COMPLETE"};
            mitosisLog("\n====== PHASE TRANSITION → %s ======\n",
                gMitosis.phase < 8 ? phaseNames[gMitosis.phase] : "?");
            mitosisLogSnapshot("transition", cellR);
            sMitosisLogPhase = gMitosis.phase;
        }

        // Fallback safety: if visible mitosis somehow starts before the
        // interior is fully prepared, finish the missing duplication work.
        // Guard with particlesDuplicated to prevent calling every frame.
        if (!gMitosis.particlesDuplicated && gCDogma.replicationProgress >= 0.995f) {
            prepareInteriorForPhysicalDivision(gIntraPhys, cellR);
            // Mark as done regardless of the check result to prevent infinite re-calls
            gMitosis.particlesDuplicated = true;
            gMitosis.dnaCheckpointPassed =
                interiorGenomeAlreadyReplicated(gIntraPhys) && gCDogma.replicationReadyForM();
            DNADistribution dnaRep = countActiveDNADistribution();
            sMitosisExpectedDNACount = dnaRep.total > 0 ? dnaRep.total : sMitosisExpectedDNACount;
            mitosisLog("S-phase replication completed before mitosis: DNA=%d\n", dnaRep.total);
        }

        // All mitosis positions are in LOCAL space (origin={0,0,0})

        if (!gMitosis.particlesDuplicated && gMitosis.phase == MITO_PROPHASE
            && gMitosis.chromatinCondensation > 0.15f) {
            int beforeCount = (int)gIntraPhys.particles.size();
            prepareInteriorForPhysicalDivision(gIntraPhys, cellR);
            gMitosis.particlesDuplicated = true; // Mark done unconditionally
            printf("[CellSim] Mitosis: prepared nuclear particles %d → %d\n",
                   beforeCount, (int)gIntraPhys.particles.size());
        }

        // During prophase→cytokinesis: gently separate the two DNA sets
        // Half=0 drifts toward +X, half=1 drifts toward -X
        // Non-nuclear particles are unaffected — they keep their normal movement
        if (gMitosis.particlesDuplicated &&
            gMitosis.phase >= MITO_PROMETAPHASE && gMitosis.phase <= MITO_CYTOKINESIS) {
            float springK = 0.050f;
            if (gMitosis.phase == MITO_METAPHASE) springK = 0.070f;
            if (gMitosis.phase == MITO_ANAPHASE) springK = 0.115f;
            if (gMitosis.phase == MITO_TELOPHASE) springK = 0.140f;
            if (gMitosis.phase == MITO_CYTOKINESIS) springK = 0.170f;
            for (auto& p : gIntraPhys.particles) {
                if (!p.active) continue;
                bool isNuclear = isNuclearParticleType(p.type);
                if (!isNuclear) continue;

                simd_float3 delta = {p.home.x - p.position.x,
                                     p.home.y - p.position.y,
                                     p.home.z - p.position.z};
                p.velocity.x += delta.x * springK;
                p.velocity.y += delta.y * (springK * 0.55f);
                p.velocity.z += delta.z * (springK * 0.55f);

                // If a genome lingers near the cleavage plane late in mitosis,
                // give it an extra shove toward its assigned daughter side.
                float signedX = ((p.mitosisHalf == 0) ? 1.0f : -1.0f) * p.position.x;
                if (gMitosis.phase >= MITO_ANAPHASE && signedX < fabsf(p.home.x) * 0.55f) {
                    p.velocity.x += ((p.mitosisHalf == 0) ? 1.0f : -1.0f) * cellR * 0.035f;
                }
            }
        }

        // ══════════════════════════════════════════════════════════════
        //  HANDS-OFF MITOSIS — particles stay natural, no position override
        //
        //  During mitosis, particles keep doing what they normally do:
        //  Brownian motion, jobs, confinement to cell membrane. The visual
        //  mitosis (furrow, spindle, chromosomes) is handled by the shader
        //  and MitosisState. Particles are NOT condensed, lerped, or killed.
        //
        //  At finalization (postDivisionComplete), splitMitosisDaughterInterior
        //  divides them by X position: left half → daughter B, right → daughter A.
        //  Each daughter gets exactly the particles that were on its side.
        // ══════════════════════════════════════════════════════════════
        DNADistribution dnaNow = countActiveDNADistribution();
        float dnaProgress = 0.0f;
        if (sMitosisExpectedDNACount > 0) {
            dnaProgress = fminf((float)dnaNow.total / (float)sMitosisExpectedDNACount, 1.0f);
        }
        gMitosis.dnaDuplicationProgress = dnaProgress;
        gMitosis.dnaCheckpointPassed =
            gMitosis.particlesDuplicated &&
            dnaNow.total >= sMitosisExpectedDNACount &&
            gCDogma.replicationReadyForM();
        gMitosis.replicationQuality = gCDogma.replicationQuality;

        float daughterR = cellR * 0.794f;
        float separationX = cellR * 0.35f;

        // ── POST-DIVISION: keep both daughter interiors visible for the
        // split-cell visualization, then finalize into 2 persistent cells. ──
        if (gMitosis.phase == MITO_COMPLETE) {
            gMitosis.postDivisionTimer += simDt;

            // Compute daughter positions / radius from the current cell geometry.
            // sDaughterPosA/B and sDaughterR are consumed immediately below by
            // the finalization block (Option 1: no preview, no delay).
            if (!gSim.cells.empty() && !sDaughtersInitialized) {
                auto& c0 = gSim.cells[0];
                float baseR = c0.radius * c0.size;
                float dR = baseR * 0.794f; // volume-conserving ∛2 split
                MitosisHalfCenters halfCenters =
                    computeMitosisHalfCenters(baseR * 0.35f, true);
                sDaughterR = dR;
                sDaughterPosA = {c0.position.x + halfCenters.local[0].x,
                                 c0.position.y + halfCenters.local[0].y,
                                 c0.position.z + halfCenters.local[0].z};
                sDaughterPosB = {c0.position.x + halfCenters.local[1].x,
                                 c0.position.y + halfCenters.local[1].y,
                                 c0.position.z + halfCenters.local[1].z};
                simd_float3 parentVel = c0.velocity;
                sDaughterVelA = {parentVel.x + 0.15f, parentVel.y, parentVel.z};
                sDaughterVelB = {parentVel.x - 0.15f, parentVel.y, parentVel.z};
                sDaughtersInitialized = true;
            }
        }

        // Atomic finalization — fires on first frame of MITO_COMPLETE.
        // The two daughters are immediately real cells in gSim.cells; the
        // normal render pipeline (syncCellInstances) draws them from the
        // same code path that draws every other cell. No preview, no swap.
        if (gMitosis.postDivisionComplete()) {
            mitosisLog("\n====== DIVISION FINALIZED ======\n");
            mitosisLogSnapshot("finalized", cellR);

            MitosisHalfCenters finalHalfCenters =
                computeMitosisHalfCenters(separationX, true);
            simd_float3 localCenterA = finalHalfCenters.local[0];
            simd_float3 localCenterB = finalHalfCenters.local[1];
            float daughterPostMitoticCondense = mitosisCondenseBlend();
            CellInterior daughterAInterior =
                splitMitosisDaughterInterior(gIntraPhys, 0, localCenterA, sDaughterR,
                                             daughterPostMitoticCondense);
            CellInterior daughterBInterior =
                splitMitosisDaughterInterior(gIntraPhys, 1, localCenterB, sDaughterR,
                                             daughterPostMitoticCondense);

            // Cell[0] becomes daughter A (the "primary" cell we track interior for)
            // cellSize in rendering = radius * size. We want cellSize = sDaughterR.
            // So keep radius = CELL_RADIUS_BASE, size = sDaughterR / CELL_RADIUS_BASE
            float daughterSizeFactor = sDaughterR / CELL_RADIUS_BASE;
            // Compute weighted split stats from the actual particle distribution
            SplitStats stats;
            {
                float massA = 0, massB = 0;
                int dnaA = 0, dnaB = 0;
                for (auto& p : gIntraPhys.particles) {
                    if (!p.active) continue;
                    bool isNuclear = (p.type == PT_DNA_NODE || p.type == PT_RNA_POLYMERASE ||
                                      p.type == PT_SPLICEOSOME || p.type == PT_NUCLEAR_PORE);
                    int half = particleMitosisHalf(p);
                    if (isNuclear) {
                        if (half == 0) dnaA++; else dnaB++;
                    } else {
                        if (half == 0) massA += p.mass; else massB += p.mass;
                    }
                }
                float totalMass = massA + massB;
                stats.cytoplasmicRatioA = (totalMass > 1e-6f) ? massA / totalMass : 0.5f;
                stats.dnaA = dnaA; stats.dnaB = dnaB;
            }

            if (!gSim.cells.empty()) {
                // IMMUTABLE snapshot — both daughters derive from this
                const SimCell parentSnapshot = gSim.cells[0];
                float sizeRatio = sDaughterR / (parentSnapshot.radius * parentSnapshot.size);

                SimCell daughterA = gSim.deriveDaughter(parentSnapshot, sDaughterPosA, {0,0,0}, sizeRatio, stats, true);
                SimCell daughterB = gSim.deriveDaughter(parentSnapshot, sDaughterPosB, {0,0,0}, sizeRatio, stats, false);

                // Propagate SAC outcome to daughters. Mitotic slippage
                // marks the daughters aneuploid (damage bump) but we
                // DON'T hard-gate future division on a single slippage
                // event — real cells with mild aneuploidy can still
                // divide (they're just more likely to trigger apoptosis
                // via p53). Leaving `tetraploid=true` on would permanently
                // lock the lineage and starve the simulation, so we only
                // set it when slippage is severe AND there are lagging
                // chromatids (clear evidence of major karyotype damage).
                bool severeSlippage = gMitosis.mitoticSlippage &&
                                      gMitosis.laggingChromatidCount >= 3;
                if (severeSlippage) {
                    daughterA.program.tetraploid = true;
                    daughterB.program.tetraploid = true;
                }
                if (gMitosis.mitoticSlippage) {
                    // Damage bump regardless; p53 gating in the next
                    // cycle will handle the response. Reasonable cells
                    // will repair over one cycle and continue.
                    daughterA.damageLevel = fminf(1.0f, daughterA.damageLevel + 0.20f);
                    daughterB.damageLevel = fminf(1.0f, daughterB.damageLevel + 0.20f);
                }
                if (gMitosis.misSegregationCount > 0) {
                    // Each mis-segregated chromosome adds a small damage
                    // tick. Micronuclei from lagging chromatids trigger
                    // cGAS-STING + ATR (Mackenzie 2017 Nature).
                    float dmg = 0.02f * (float)gMitosis.misSegregationCount;
                    daughterA.damageLevel = fminf(1.0f, daughterA.damageLevel + dmg * 0.5f);
                    daughterB.damageLevel = fminf(1.0f, daughterB.damageLevel + dmg * 0.5f);
                }

                // Parent (cell[0]) becomes daughter A — keep its UID
                daughterA.cellUid = parentSnapshot.cellUid;
                gSim.cells[0] = daughterA;
                gSim.cells.push_back(daughterB);

                mitosisLog("  Daughter A: uid=%d gen=%d telo=%.0f ATP=%.1f biomass=%.3f dna=%d\n",
                    daughterA.cellUid, daughterA.generation, daughterA.telomere, daughterA.ATP, daughterA.biomass, stats.dnaA);
                mitosisLog("  Daughter B: uid=%d gen=%d telo=%.0f ATP=%.1f biomass=%.3f dna=%d\n",
                    daughterB.cellUid, daughterB.generation, daughterB.telomere, daughterB.ATP, daughterB.biomass, stats.dnaB);
                mitosisLog("  Split ratio: %.2f (A) / %.2f (B)\n", stats.cytoplasmicRatioA, 1.0f - stats.cytoplasmicRatioA);
                mitosisLog("  SAC outcome: mcc=%.3f misSeg=%d lagging=%d slippage=%d auroraB=%.2f\n",
                    gMitosis.mccLevel, gMitosis.misSegregationCount, gMitosis.laggingChromatidCount,
                    gMitosis.mitoticSlippage ? 1 : 0, gMitosis.auroraBActivity);

                // Persist interiors keyed by cellUid
                gIntraPhys = daughterAInterior.phys;
                gIntraInitialized = true;
                gCellInteriors[daughterB.cellUid] = daughterBInterior;
                splitOrganelleMotionState(parentSnapshot.cellUid,
                                          daughterA.cellUid,
                                          daughterB.cellUid);

                // deriveDaughter() already reset each daughter's per-cell
                // replication + mitosis program without wiping inherited
                // lineage state. Just mark the programs as resident.
                gSim.cells[0].program.cdogmaInitialized = true;
                // cells.back() is the newly-pushed daughter B
                gSim.cells.back().program.cdogmaInitialized = true;
                gCDogma = gSim.cells[0].program.cdogma;
                gBoundDogmaCellUid = daughterA.cellUid;

                // Resync gCellBio so it doesn't overwrite daughter A's ATP
                syncPrimarySimATPIntoBioEngine();

                gPostMitosisPairUidA = daughterA.cellUid;
                gPostMitosisPairUidB = daughterB.cellUid;
            }

            gMitosis.active = false;
            gMitosis.furrowDepth = 0;
            gMitosis.poleSeparation = 0;
            gMitosis.phase = MITO_NONE;
            gBoundMitosisCellUid = gSim.cells.empty() ? -1 : gSim.cells[0].cellUid;
            gChromatinDecondenseState = 0.0f;
            gSim.mitosisVisualizationComplete = true;
            sDaughtersInitialized = false;
            sMitosisDNAStartCount = 0;
            sMitosisExpectedDNACount = 0;
            gLastCellPos = gSim.cells.empty() ? cellPos : gSim.cells[0].position;

            // Hold the camera on the daughter pair long enough for the split
            // to read as one continuous event rather than snapping back to the
            // wider colony before the daughters visually settle.
            gPostMitosisPairCameraTimer = 2.75f;
            gPostMitosisPairA = sDaughterPosA;
            gPostMitosisPairB = sDaughterPosB;
            gPostMitosisPairRadius = sDaughterR;
        }
    }
    if (cellPhase != 3 && gMitosis.active && gMitosis.phase < MITO_CYTOKINESIS) {
        // Cell left M-phase prematurely (shouldn't happen, but safety)
    }

    retargetNuclearHomesForMitosis(gIntraPhys, cellR);

    // During mitosis, also bias non-nuclear particles toward their assigned
    // daughter half so cytoplasmic cargo (vesicles, free molecules) follows
    // the furrow rather than clustering at the equator and getting trapped
    // in the bridge when cytokinesis pinches in.
    if (gMitosis.active && gMitosis.phase >= MITO_METAPHASE && gMitosis.phase <= MITO_CYTOKINESIS) {
        float furrowProgress = 0.0f;
        if (gMitosis.phase == MITO_METAPHASE)   furrowProgress = 0.10f;
        else if (gMitosis.phase == MITO_ANAPHASE)   furrowProgress = 0.30f;
        else if (gMitosis.phase == MITO_TELOPHASE)  furrowProgress = 0.55f;
        else if (gMitosis.phase == MITO_CYTOKINESIS) furrowProgress = 0.85f;

        float targetSeparation = cellR * 0.35f * furrowProgress;
        for (auto& p : gIntraPhys.particles) {
            if (!p.active || isNuclearParticleType(p.type)) continue;
            if (p.job == JOB_CONSUMED) continue; // let consumed particles respawn normally

            // Assign half once per particle by current X side, and remember it.
            if (p.mitosisHalf < 0) {
                p.mitosisHalf = (p.position.x >= 0.0f) ? 0 : 1;
            }
            int half = p.mitosisHalf;
            float destX = (half == 0) ? +targetSeparation : -targetSeparation;

            // Gentle spring pull — strong enough to clear the equator before
            // the furrow seals, weak enough that active transport jobs still
            // dominate particle motion.
            float pullStrength = 0.015f + 0.06f * furrowProgress;
            p.velocity.x += (destX - p.position.x) * pullStrength;

            // Only push non-job-driven cytoplasm (free molecules, water, ATP
            // pool). Particles with active transport jobs keep their targets.
            if (p.job != JOB_IDLE) continue;

            // Shift idle particles' home so they settle on the correct side.
            p.home.x = p.home.x * 0.92f + destX * 0.08f;
        }
    }

    // Step Langevin dynamics
    gIntraPhys.step(simDt);
    enforceCleavageCompartmentBarrier(gIntraPhys, cellR);

    // ── DIAGNOSTIC LOG (every 120 frames → ~2 seconds) ──────────────
    {
        static int diagFrame = 0;
        diagFrame++;
        bool tracePostDivision = postMitosisPairHoldActive();
        if (!tracePostDivision) {
            for (const auto& c : gSim.cells) {
                if (c.postDivisionRecovery > 0.0f) {
                    tracePostDivision = true;
                    break;
                }
            }
        }
        int diagInterval = tracePostDivision ? 15 : 120;
        if (diagFrame % diagInterval == 0) {
            const char* diagPath = gDiagLogPath.empty()
                ? "/tmp/cellsim_diag.log"
                : gDiagLogPath.c_str();
            FILE* logF = fopen(diagPath, "a");
            if (logF) {
                fprintf(logF, "\n=== CellSim Diagnostic Frame %d ===\n", diagFrame);
                fprintf(logF, "cellPos=(%.2f, %.2f, %.2f)  cellSize=%.3f  cellR=%.3f\n",
                        cellPos.x, cellPos.y, cellPos.z, cellSize, cellR);
                fprintf(logF, "cellVisualRadius=%.3f (cellSize itself)\n", cellSize);
                fprintf(logF, "physCellR=%.3f (confinement radius)\n", gIntraPhys.params.cellRadius);
                fprintf(logF, "mitosis: active=%d phase=%d furrow=%.2f poleSep=%.2f\n",
                        gMitosis.active, gMitosis.phase, gMitosis.furrowDepth, gMitosis.poleSeparation);
                int visibleDNA = 0, visibleSisterDNA = 0;
                countActiveDNAVisuals(gIntraPhys, visibleDNA, visibleSisterDNA);
                fprintf(logF,
                        "replication: progress=%.3f bases=%d/%d forks=%d unresolved=%d chk1=%.2f quality=%.2f ready=%d visibleDNA=%d sisterDNA=%d\n",
                        gCDogma.replicationProgress,
                        gCDogma.replicatedBaseCount(), HBB_LENGTH,
                        gCDogma.countActiveReplicationForks(),
                        gCDogma.unresolvedReplicationErrors,
                        gCDogma.chk1Signal,
                        gCDogma.replicationQuality,
                        gCDogma.replicationReadyForM() ? 1 : 0,
                        visibleDNA, visibleSisterDNA);
                if (gMitosis.active) {
                    fprintf(logF, "  poleA=(%.2f,%.2f,%.2f) poleB=(%.2f,%.2f,%.2f)\n",
                            gMitosis.poleA.x, gMitosis.poleA.y, gMitosis.poleA.z,
                            gMitosis.poleB.x, gMitosis.poleB.y, gMitosis.poleB.z);
                    fprintf(logF, "  orgDup=%.2f orgMig=%.2f orgOffsetY=%.3f\n",
                            gMitosis.organelleDuplication, gMitosis.organelleMigration,
                            gMitosis.organelleOffsetY());
                    fprintf(logF, "  envBreak=%.2f spindleAsm=%.2f spindleDis=%.2f attach=%.2f align=%.2f checkpoint=%.2f chromatidSep=%.2f ring=%.2f\n",
                            gMitosis.nuclearEnvelopeBreakdown, gMitosis.spindleAssembly,
                            gMitosis.spindleDisassembly, gMitosis.kinetochoreAttachment,
                            gMitosis.metaphaseAlignment, gMitosis.spindleCheckpoint,
                            gMitosis.chromatidSeparation, gMitosis.contractileRingAssembly);
                    fprintf(logF, "  daughterA=(%.2f,%.2f,%.2f) daughterB=(%.2f,%.2f,%.2f)\n",
                            gMitosis.daughterA.x, gMitosis.daughterA.y, gMitosis.daughterA.z,
                            gMitosis.daughterB.x, gMitosis.daughterB.y, gMitosis.daughterB.z);
                }

                // Count particles by type and check how many are outside cell
                int typeCounts[PT_COUNT] = {};
                int outsideCount = 0;
                float maxDist = 0;
                for (auto& p : gIntraPhys.particles) {
                    if (!p.active) continue;
                    typeCounts[p.type]++;
                    // Particles are in LOCAL space (origin = cell center)
                    float dx = p.position.x;
                    float dy = p.position.y;
                    float dz = p.position.z;
                    float dist = sqrtf(dx*dx + dy*dy + dz*dz);
                    if (dist > maxDist) maxDist = dist;
                    if (dist > cellSize) outsideCount++;
                }
                fprintf(logF, "\nParticles: %d total, %d OUTSIDE cell (cellSize=%.2f)\n",
                        (int)gIntraPhys.particles.size(), outsideCount, cellSize);
                fprintf(logF, "MaxDist from center: %.3f (limit=%.3f)\n", maxDist, cellSize);

                const char* typeNames[] = {
                    "ORGANELLE","DNA_NODE","DNA_POL","RNA_POL","SPLICEO","PRE_MRNA",
                    "MRNA","NUC_PORE","RIB_S","RIB_L","TRNA","POLYPEP",
                    "CHAPERONE","VES_COPII","VES_SECR","MOLECULE"
                };
                for (int t = 0; t < PT_COUNT; t++) {
                    if (typeCounts[t] > 0) {
                        // Find min/max Y for this type
                        float minY = 1e9, maxY = -1e9, sumDist = 0;
                        int cnt = 0, outCnt = 0;
                        for (auto& p : gIntraPhys.particles) {
                            if (!p.active || p.type != t) continue;
                            if (p.position.y < minY) minY = p.position.y;
                            if (p.position.y > maxY) maxY = p.position.y;
                            float dx = p.position.x - p.confineCenter.x;
                            float dy = p.position.y - p.confineCenter.y;
                            float dz = p.position.z - p.confineCenter.z;
                            float d = sqrtf(dx*dx+dy*dy+dz*dz);
                            sumDist += d;
                            cnt++;
                            if (d > p.confineRadius) outCnt++;
                        }
                        fprintf(logF, "  %-12s n=%3d  Y=[%.2f,%.2f]  avgDist=%.3f  outside=%d\n",
                                typeNames[t], typeCounts[t], minY, maxY,
                                cnt>0?sumDist/cnt:0, outCnt);
                    }
                }

                // Cell instances (for rendering)
                fprintf(logF, "\nRender focus: follow=%d solo=%d selectedIdx=%d selectedUid=%d activeIdx=%d\n",
                        gFollowCell ? 1 : 0,
                        gSoloFocusCell ? 1 : 0,
                        gSelectedCell, gSelectedCellUid, activeFocusCellIndex());
                fprintf(logF, "PostMitosisPair: timer=%.2f uidA=%d uidB=%d center=(%.2f,%.2f,%.2f) r=%.3f\n",
                        gPostMitosisPairCameraTimer,
                        gPostMitosisPairUidA,
                        gPostMitosisPairUidB,
                        (gPostMitosisPairA.x + gPostMitosisPairB.x) * 0.5f,
                        (gPostMitosisPairA.y + gPostMitosisPairB.y) * 0.5f,
                        (gPostMitosisPairA.z + gPostMitosisPairB.z) * 0.5f,
                        gPostMitosisPairRadius);
                fprintf(logF, "\nCellInstances: %d\n", (int)gCellInstances.size());
                for (int ci = 0; ci < (int)gCellInstances.size(); ci++) {
                    auto& inst = gCellInstances[ci];
                    fprintf(logF, "  cell[%d] pos=(%.2f,%.2f,%.2f) r=%.3f furrow=%.2f\n",
                            ci, inst.position.x, inst.position.y, inst.position.z,
                            inst.radius, inst.furrowDepth);
                }

                fprintf(logF, "\nSimCells: %d\n", (int)gSim.cells.size());
                for (int ci = 0; ci < (int)gSim.cells.size(); ci++) {
                    auto& c = gSim.cells[ci];
                    fprintf(logF,
                            "  sim[%d] alive=%d uid=%d clone=%d phase=%d fate=%d size=%.3f biomass=%.3f ATP=%.1f recovery=%.2f pos=(%.2f,%.2f,%.2f)\n",
                            ci, c.alive ? 1 : 0, c.cellUid, c.cloneId, c.phase, c.fate,
                            c.radius * c.size, c.biomass, c.ATP, c.postDivisionRecovery,
                            c.position.x, c.position.y, c.position.z);
                    if (ci > 0) {
                        auto it = gCellInteriors.find(c.cellUid);
                        if (it != gCellInteriors.end() && it->second.initialized) {
                            int dnaCount = 0;
                            for (const auto& p : it->second.phys.particles) {
                                if (p.active && p.type == PT_DNA_NODE) dnaCount++;
                            }
                            fprintf(logF, "    interior dna=%d particles=%d\n",
                                    dnaCount, (int)it->second.phys.particles.size());
                        } else {
                            fprintf(logF, "    interior missing\n");
                        }
                    }
                }
                if (gCellBio.initialized) {
                    auto& pool = gCellBio.metabolism.pool;
                    fprintf(logF,
                            "\nBioPrimary: ATP_mM=%.3f ADP_mM=%.3f AMP_mM=%.3f totalAden=%.3f ATP%%=%.1f\n",
                            pool.ATP, pool.ADP, pool.AMP, pool.totalAdenylate(),
                            bioATPToSimPercent(pool.ATP));
                }

                fclose(logF);
            }
        }
    }

    // ── CSV EXPORT (every 600 frames ≈ 10s wall-clock) ────────────────
    // Writes a single row per sample to logs/export/cellsim_timeseries.csv
    // Columns match data/reference/hela_reference.csv for easy comparison.
    {
        static int csvFrame = 0;
        csvFrame++;
        if (csvFrame % 600 == 0) {
            static bool csvHeaderWritten = false;
            std::string csvPath = "/Users/henry/CellSim/logs/export/cellsim_timeseries.csv";
            FILE* csvF = fopen(csvPath.c_str(), csvHeaderWritten ? "a" : "w");
            if (csvF) {
                if (!csvHeaderWritten) {
                    fprintf(csvF,
                        "wall_sec,bio_sec,bio_hours,"
                        "cell_count,phase_G1,phase_S,phase_G2,phase_M,"
                        "ATP_mM,ADP_mM,AMP_mM,energy_charge,NAD_plus_mM,NADH_mM,NAD_ratio,"
                        "glucose_mM,O2_mM,pyruvate_mM,lactate_mM,citrate_mM,"
                        "glycolysis_flux,TCA_flux,OxPhos_flux,O2_consumption,"
                        "membrane_potential_mV,"
                        "telomere_bp,generation,damage_level,"
                        "mean_ATP_sim,mean_biomass\n");
                    csvHeaderWritten = true;
                }
                float wallSec = csvFrame / 60.0f;
                float bioSec = gSim.bioTime;
                float bioHours = bioSec / 3600.0f;

                // Phase counts
                int phG1=0, phS=0, phG2=0, phM=0;
                float meanATP=0, meanBiomass=0;
                float meanTelo=0; int meanGen=0; float meanDmg=0;
                int nAlive = 0;
                for (auto& c : gSim.cells) {
                    if (!c.alive) continue;
                    nAlive++;
                    switch(c.phase) {
                        case 0: phG1++; break;
                        case 1: phS++; break;
                        case 2: phG2++; break;
                        case 3: phM++; break;
                    }
                    meanATP += c.ATP;
                    meanBiomass += c.biomass;
                    meanTelo += c.telomere;
                    meanGen += c.generation;
                    meanDmg += c.damageLevel;
                }
                if (nAlive > 0) {
                    meanATP /= nAlive;
                    meanBiomass /= nAlive;
                    meanTelo /= nAlive;
                    meanGen /= nAlive;
                    meanDmg /= nAlive;
                }

                // Biochemistry engine values (cell[0])
                auto& pool = gCellBio.metabolism.pool;
                auto& flux = gCellBio.metabolism.flux;
                float Vm = gCellBio.membrane.getMembranePotential();

                fprintf(csvF,
                    "%.1f,%.1f,%.3f,"
                    "%d,%d,%d,%d,%d,"
                    "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,"
                    "%.4f,%.5f,%.4f,%.4f,%.4f,"
                    "%.5f,%.5f,%.5f,%.6f,"
                    "%.1f,"
                    "%.0f,%d,%.4f,"
                    "%.1f,%.3f\n",
                    wallSec, bioSec, bioHours,
                    nAlive, phG1, phS, phG2, phM,
                    pool.ATP, pool.ADP, pool.AMP, pool.energyCharge(),
                    pool.NAD_plus, pool.NADH, pool.NAD_ratio(),
                    pool.glucose, pool.O2, pool.pyruvate, pool.lactate, pool.citrate,
                    flux.totalATPproduction, flux.citrateSynthase, flux.atpSynthase, flux.totalO2consumption,
                    Vm,
                    meanTelo, (int)meanGen, meanDmg,
                    meanATP, meanBiomass);
                fclose(csvF);
            }
        }
    }

    // ── Convert particles to GPU instances ──────────────────────────
    // Particle positions are always local. During MITO_COMPLETE, each half
    // is re-centered into its own daughter cell in world space.
    auto toWorldLocal = [&](simd_float3 local) -> simd_float3 {
        return {local.x + cellPos.x, local.y + cellPos.y, local.z + cellPos.z};
    };
    auto particleToWorld = [&](const IntraParticle& p) -> simd_float3 {
        if (gMitosis.active && gMitosis.phase == MITO_COMPLETE && sDaughtersInitialized) {
            int half = particleMitosisHalf(p);
            simd_float3 daughterCenter = (half == 0) ? sDaughterPosA : sDaughterPosB;
            simd_float3 localCenter = p.confineCenter;
            return {daughterCenter.x + (p.position.x - localCenter.x),
                    daughterCenter.y + (p.position.y - localCenter.y),
                    daughterCenter.z + (p.position.z - localCenter.z)};
        }
        return toWorldLocal(p.position);
    };

    std::vector<GPUAtomInstance> allAtoms;
    std::vector<GPUBondInstance> allBonds;
    allAtoms.reserve(gIntraPhys.particles.size() + 500 + (int)gMediumChemicals.size());
    allBonds.reserve(gIntraPhys.particles.size());

    // ── Medium chemicals (drawn FIRST so they survive truncation) ────
    // The fluid container itself is drawn via its own dedicated pipeline
    // (FluidRender.metal) in the main render pass — this loop only
    // populates the molecule pipeline with the chemical particles that
    // float inside the fluid.
    renderMediumChemicals(allAtoms, allBonds, gSim, (float)glfwGetTime());
    // Apoptotic bodies use the same atom-instance pipeline so they are
    // drawn in the same pass as medium molecules — no extra shader.
    renderApoptoticBodies(allAtoms);
    bool renderPrimaryInterior = !soloFocusEnabled() || activeFocusCellIndex() == 0;
    int focusCellIdx = activeFocusCellIndex();

    int prevDnaA = -1, prevDnaB = -1;
    float primaryPostMitoticCondense = primaryPostMitoticCondenseBlend();
    float primaryDNACondenseBlend = gMitosis.active ? mitosisCondenseBlend() : primaryPostMitoticCondense;

    if (renderPrimaryInterior) {
        for (int i = 0; i < (int)gIntraPhys.particles.size(); i++) {
            auto& p = gIntraPhys.particles[i];
            if (!p.active) continue;

            simd_float3 worldPos = particleToWorld(p);

            const char* protFile = particleToProtein(p.type);
            const MoleculeData* protMol = protFile ? gProtCache.get(protFile) : nullptr;
            const char* molId = (p.type == PT_MOLECULE) ? particleTagToMolecule(p.tag) : nullptr;
            const MoleculeData* smallMol = molId ? gMolCache.get(molId) : nullptr;
            // Always render atoms when a PDB is available. Idle proteins
            // use a higher stride so the perf cost scales. Step G of the
            // sizing refactor — every particle type with known atomic
            // structure emerges from real atoms, not a flat sphere.
            // PT_DNA_NODE is exempt: 760+ nodes × 1266 atoms/nucleosome
            // would be 960k atoms — we keep those as abstract beads and
            // render one representative nucleosome separately (below).
            //
            // During mitosis, the ribosomes / Pol II / spliceosomes / free
            // molecules should NOT drown the spindle in a colorful blob
            // at the equator. Biologically, transcription and translation
            // both halt in mitosis (Prescott & Bender 1962; Sivan 2007)
            // and the nucleolus disassembles. We suppress the atomic detail
            // of non-nuclear machinery during mitosis and fall back to
            // simple spheres. The spindle and chromosomes remain bright.
            bool inMitosis = gMitosis.active && gMitosis.phase >= MITO_PROMETAPHASE
                                             && gMitosis.phase <= MITO_TELOPHASE;
            bool showFullProtein = protMol && p.type != PT_DNA_NODE && !inMitosis;
            bool showChemicalMolecule = smallMol && p.type == PT_MOLECULE && !inMitosis;
            // During active mitosis use the live mitotic blend; after split,
            // use the same per-daughter recovery blend as the secondary cells.
            float condenseBlend = primaryDNACondenseBlend;
            bool mitoticDNA = p.type == PT_DNA_NODE && condenseBlend > 0.02f;
            bool replicatedSisterDNA = p.type == PT_DNA_NODE && !gMitosis.active && p.mitosisHalf == 1;
            float renderRadius = p.radius;
            float renderCr = p.colorR, renderCg = p.colorG, renderCb = p.colorB;
            if (mitoticDNA) {
                int slot = mitosisChromosomeSlot(p);
                simd_float3 tint = chromosomeSlotColor(slot);
                float condense = condenseBlend;
                float slotSpan = fmaxf((float)HBB_LENGTH / (float)MitosisState::NUM_CHROMO, 1.0f);
                float slotOffset = fmodf((float)fmax(p.stateIndex, 0), slotSpan);
                float slotLocal = (slotOffset + 0.5f) / slotSpan - 0.5f;
                float centromere = 1.0f - fminf(fabsf(slotLocal) * 2.2f, 1.0f);
                renderRadius *= 1.0f + condense * 1.35f + centromere * 0.20f;
                renderCr = p.colorR * (1.0f - condense * 0.35f) + tint.x * (0.45f + condense * 0.55f);
                renderCg = p.colorG * (1.0f - condense * 0.35f) + tint.y * (0.45f + condense * 0.55f);
                renderCb = p.colorB * (1.0f - condense * 0.35f) + tint.z * (0.45f + condense * 0.55f);
                renderCr = fminf(renderCr + centromere * 0.18f, 1.0f);
                renderCg = fminf(renderCg + centromere * 0.18f, 1.0f);
            } else if (replicatedSisterDNA) {
                float replTint = 0.35f + 0.55f * gCDogma.replicationProgress;
                renderRadius *= 1.08f;
                renderCr = p.colorR * (1.0f - 0.45f * replTint) + 0.95f * 0.45f * replTint;
                renderCg = p.colorG * (1.0f - 0.45f * replTint) + 0.85f * 0.45f * replTint;
                renderCb = p.colorB * (1.0f - 0.45f * replTint) + 0.30f * 0.45f * replTint;
            }

            if (showFullProtein) {
                float rot = (float)i * 1.7f;
                // Stride scales with atom count AND idle/active state.
                // Active or focal particles get dense sampling; idle
                // proteins are downsampled to keep the frame budget.
                int totalAtoms = (int)protMol->atoms.size();
                bool idle = (p.job == JOB_IDLE);
                int stride = 1;
                if (totalAtoms > 1500) stride = idle ? 8 : 4;
                else if (totalAtoms > 500) stride = idle ? 6 : 2;
                else if (totalAtoms > 200) stride = idle ? 3 : 1;
                else stride = idle ? 2 : 1;
                addProteinBackbone(allAtoms, allBonds, *protMol, worldPos, cellR,
                                   p.colorR, p.colorG, p.colorB, rot, stride);
            } else if (showChemicalMolecule) {
                float minHalfExtentFrac = fminf(0.018f, (renderRadius * 1.45f) / fmaxf(cellR, 0.001f));
                float molScale = chemistryScaleFromRealSize(*smallMol, cellR, 18.0f, minHalfExtentFrac);
                int stride = (smallMol->atoms.size() > 96) ? 2 : 1;
                addMoleculeGeometry(allAtoms, allBonds, *smallMol, worldPos, molScale,
                                    renderCr, renderCg, renderCb, (float)i * 1.7f, stride,
                                    fmaxf(renderRadius * 0.22f, cellR * 0.0008f),
                                    0.38f, 0.56f, 0.20f, 0.0f, 0.0f,
                                    p.stateIndex + i * 19);
            } else {
                addSphere(allAtoms, worldPos, renderRadius, renderCr, renderCg, renderCb);
            }

            // DNA backbone connections
            if (p.type == PT_DNA_NODE) {
                int slot = mitosisChromosomeSlot(p);
                simd_float3 tint = chromosomeSlotColor(slot);
                float bondRadius = p.radius * 0.4f;
                float rungRadius = p.radius * 0.25f;
                float bondCr = 0.7f, bondCg = 0.5f, bondCb = 0.2f;
                float rungCr = 0.4f, rungCg = 0.4f, rungCb = 0.5f;
                if (mitoticDNA) {
                    float condense = mitosisCondenseBlend();
                    bondRadius *= 1.0f + condense * 1.25f;
                    rungRadius *= 1.0f + condense * 0.95f;
                    bondCr = tint.x; bondCg = tint.y; bondCb = tint.z;
                    rungCr = 0.75f * tint.x + 0.20f;
                    rungCg = 0.75f * tint.y + 0.20f;
                    rungCb = 0.75f * tint.z + 0.20f;
                }
                if (p.tetherId >= 0 && p.tetherId < (int)gIntraPhys.particles.size()) {
                    auto& other = gIntraPhys.particles[p.tetherId];
                    addStick(allBonds, worldPos, particleToWorld(other), bondRadius,
                             bondCr, bondCg, bondCb);
                }
                if (i + 1 < (int)gIntraPhys.particles.size() &&
                    gIntraPhys.particles[i + 1].type == PT_DNA_NODE &&
                    gIntraPhys.particles[i + 1].stateIndex == p.stateIndex) {
                    addStick(allBonds, worldPos, particleToWorld(gIntraPhys.particles[i + 1]),
                             rungRadius, rungCr, rungCg, rungCb);
                }
            }

            // Vesicles: add cargo sphere inside
            if (p.type == PT_VESICLE_COPII || p.type == PT_VESICLE_SECRETORY) {
                addSphere(allAtoms, worldPos, p.radius * 0.4f, 0.4f, 0.3f, 0.6f);
            }

            // ── mRNA chain backbone ────────────────────────────────
            // Pre-mRNA and mature mRNA beads are tethered to form a
            // polymer. Render the bond between consecutive beads so the
            // transcript looks like a real strand trailing the RNA Pol.
            // Bright pink/magenta for pre-mRNA, cyan for mature mRNA.
            if ((p.type == PT_PRE_MRNA || p.type == PT_MRNA_NODE)
                && p.tetherId >= 0 && p.tetherId < (int)gIntraPhys.particles.size()) {
                simd_float3 op = particleToWorld(gIntraPhys.particles[p.tetherId]);
                float bondR = p.radius * 0.50f;
                float br = (p.type == PT_MRNA_NODE) ? 0.2f : 0.9f;
                float bg = (p.type == PT_MRNA_NODE) ? 0.9f : 0.3f;
                float bb = (p.type == PT_MRNA_NODE) ? 0.9f : 0.7f;
                addStick(allBonds, worldPos, op, bondR, br, bg, bb);
                // Bead slightly larger so chain reads clearly
                addSphere(allAtoms, worldPos, p.radius * 1.4f, br, bg, bb);
            }
        }
    }

    // ── MITOSIS RENDERING (full staged sequence through cytokinesis) ──
    // Runs any time the persistent condensation blend is above zero — this
    // keeps chromosome rods, spindle, and centrosomes on screen while they
    // fade into the interphase nucleus, rather than snapping off the frame
    // cytokinesis ends.
    float chromRenderAlpha = gMitosis.active ? mitosisCondenseBlend() : 0.0f;
    if (renderPrimaryInterior && chromRenderAlpha > 0.02f) {
        float cs = cellSize;
        float spindleAlpha = fmaxf(0.08f, 1.0f - gMitosis.spindleDisassembly * 0.92f);

        // ── Centrosomes (poles) — CLAMPED inside cell membrane ────────
        // Poles are in local space. Clamp distance from origin to 85% of cellR
        // so they NEVER escape the membrane.
        auto clampInside = [&](simd_float3 local, float maxR) -> simd_float3 {
            float d = sqrtf(local.x*local.x + local.y*local.y + local.z*local.z);
            if (d > maxR && d > 0.001f) {
                float s = maxR / d;
                return {local.x*s, local.y*s, local.z*s};
            }
            return local;
        };
        simd_float3 clampedPoleA = clampInside(gMitosis.poleA, cellR * 0.85f);
        simd_float3 clampedPoleB = clampInside(gMitosis.poleB, cellR * 0.85f);

        // ── CENTROSOME LOG — every frame during mitosis ─────────────
        {
            static int cFrame = 0;
            cFrame++;
            const char* centrosomePath = gCentrosomeLogPath.empty()
                ? "/tmp/cellsim_centrosome.log"
                : gCentrosomeLogPath.c_str();
            FILE* cLog = fopen(centrosomePath, (cFrame == 1) ? "w" : "a");
            if (cLog) {
                simd_float3 wA = toWorldLocal(clampedPoleA);
                simd_float3 wB = toWorldLocal(clampedPoleB);
                // Distance from cell center (world)
                float dA = sqrtf((wA.x-cellPos.x)*(wA.x-cellPos.x)+(wA.y-cellPos.y)*(wA.y-cellPos.y)+(wA.z-cellPos.z)*(wA.z-cellPos.z));
                float dB = sqrtf((wB.x-cellPos.x)*(wB.x-cellPos.x)+(wB.y-cellPos.y)*(wB.y-cellPos.y)+(wB.z-cellPos.z)*(wB.z-cellPos.z));
                // Local distance (should be < cellR*0.85)
                float dAlocal = sqrtf(clampedPoleA.x*clampedPoleA.x+clampedPoleA.y*clampedPoleA.y+clampedPoleA.z*clampedPoleA.z);
                float dBlocal = sqrtf(clampedPoleB.x*clampedPoleB.x+clampedPoleB.y*clampedPoleB.y+clampedPoleB.z*clampedPoleB.z);
                float rawA = sqrtf(gMitosis.poleA.x*gMitosis.poleA.x+gMitosis.poleA.y*gMitosis.poleA.y+gMitosis.poleA.z*gMitosis.poleA.z);
                float rawB = sqrtf(gMitosis.poleB.x*gMitosis.poleB.x+gMitosis.poleB.y*gMitosis.poleB.y+gMitosis.poleB.z*gMitosis.poleB.z);
                bool escaped = (dA > cellSize || dB > cellSize);
                if (cFrame % 30 == 0 || escaped) { // Every 0.5s or if escaped
                    fprintf(cLog, "f=%d phase=%d cellR=%.2f cellSize=%.2f poleSep=%.2f "
                            "rawA=%.3f rawB=%.3f clampA=%.3f clampB=%.3f limit=%.3f "
                            "worldDistA=%.3f worldDistB=%.3f %s\n",
                            cFrame, gMitosis.phase, cellR, cellSize, gMitosis.poleSeparation,
                            rawA, rawB, dAlocal, dBlocal, cellR*0.85f,
                            dA, dB, escaped ? "*** ESCAPED ***" : "ok");
                }
                fclose(cLog);
            }
        }

        addSphere(allAtoms, toWorldLocal(clampedPoleA), cs * 0.025f, 0.9f * spindleAlpha, 0.9f * spindleAlpha, 0.3f);
        addSphere(allAtoms, toWorldLocal(clampedPoleB), cs * 0.025f, 0.9f * spindleAlpha, 0.9f * spindleAlpha, 0.3f);
        addSphere(allAtoms, toWorldLocal(clampedPoleA), cs * 0.04f, 0.5f * spindleAlpha, 0.5f * spindleAlpha, 0.15f);
        addSphere(allAtoms, toWorldLocal(clampedPoleB), cs * 0.04f, 0.5f * spindleAlpha, 0.5f * spindleAlpha, 0.15f);

        // ── Astral microtubules — pole → cell cortex ─────────────────
        // Each MT is a ray from the pole in a hemispheric direction; we
        // solve ray-sphere intersection against the cortex so the tip is
        // guaranteed to land ON the membrane, never outside. This is the
        // biological reality — astral MT plus-ends dock into cortical
        // dynein at the inner membrane face. Ref: Kiyomitsu 2012 Dev Cell.
        {
            constexpr int ASTRAL_PER_POLE = 24;
            float membraneR = cellR * 0.90f; // render to just inside membrane
            for (int pole = 0; pole < 2; pole++) {
                simd_float3 polePosLocal = (pole == 0) ? clampedPoleA : clampedPoleB;
                simd_float3 polePos = toWorldLocal(polePosLocal);
                float sign = (pole == 0) ? 1.0f : -1.0f;
                for (int f = 0; f < ASTRAL_PER_POLE; f++) {
                    float u = (float)f / (float)(ASTRAL_PER_POLE - 1);
                    float v = ((f * 7) % ASTRAL_PER_POLE) / (float)(ASTRAL_PER_POLE - 1);
                    float theta = u * (float)M_PI;        // polar angle
                    float phi   = v * 2.0f * (float)M_PI; // azimuth
                    float dirX = sign * cosf(theta);
                    float dirY = sinf(theta) * cosf(phi);
                    float dirZ = sinf(theta) * sinf(phi);

                    // Ray-sphere intersection: find t > 0 such that
                    //   |polePosLocal + t·dir|² = membraneR².
                    // Solves t² + 2·(P·D)·t + (|P|² − R²) = 0.
                    float px = polePosLocal.x, py = polePosLocal.y, pz = polePosLocal.z;
                    float pDotD = px*dirX + py*dirY + pz*dirZ;
                    float pSq   = px*px + py*py + pz*pz;
                    float disc  = pDotD * pDotD - (pSq - membraneR * membraneR);
                    float tHit  = (disc > 0.0f) ? (-pDotD + sqrtf(disc)) : cellR * 0.35f;
                    // Random length jitter (some MTs shorter — dynamic
                    // instability). Clamp to positive so they never start
                    // at the pole and extend backward.
                    float len = fmaxf(cellR * 0.25f, tHit);
                    len *= (0.70f + 0.30f * fabsf(sinf((float)f * 3.17f
                                                        + gMitosis.phaseTimer * 2.3f)));
                    // Re-clamp after jitter so we never exceed the cortex.
                    if (len > tHit) len = tHit;

                    simd_float3 tipLocal = {
                        px + dirX * len,
                        py + dirY * len,
                        pz + dirZ * len
                    };
                    // Belt-and-braces clamp to the cell sphere.
                    simd_float3 originLocal = {0, 0, 0};
                    tipLocal = clampPointToSphere(tipLocal, originLocal, membraneR);
                    simd_float3 tip = toWorldLocal(tipLocal);
                    addStick(allBonds, polePos, tip, cs * 0.0018f,
                             0.30f * spindleAlpha, 0.85f * spindleAlpha, 0.40f * spindleAlpha);
                }
            }
        }

        // ── Interpolar microtubules — pole-to-pole bridging bundles ──
        // Real metaphase spindle is a fusiform (lemon / rugby-ball) shape
        // filling most of the cell. Interpolar bundles are MT fibers that
        // arc from one pole through the equator to the opposite pole;
        // their lateral displacement from the spindle axis is maximum at
        // the equator (the spindle's "belly") and minimum at the poles.
        // Each bundle is a smooth curve, not two straight sticks — this
        // gives the characteristic ellipsoidal spindle silhouette.
        //
        // Each bundle = 3-10 individual MTs crosslinked by PRC1 / EG5 /
        // KIF4A. We render N=24 bundles distributed around the spindle
        // axis as ~8-segment curves so the spindle reads as a continuous
        // elastic ellipse.
        // Ref: Mastronarde 1993 JCB; Kajtez 2016 Nat Commun (HeLa spindle
        // barrel shape measured by light-sheet microscopy).
        if (gMitosis.poleSeparation > 0.2f && spindleAlpha > 0.1f) {
            constexpr int INTERPOLAR_BUNDLES = 28;
            constexpr int BUNDLE_SEGMENTS = 12;
            float bundleRMax = cellR * 0.42f;   // equatorial belly
            float bundleRPole = cellR * 0.04f;  // bundles converge at pole
            float tNow = gMitosis.phaseTimer;
            for (int f = 0; f < INTERPOLAR_BUNDLES; f++) {
                float ang = ((float)f / INTERPOLAR_BUNDLES) * 2.0f * (float)M_PI;
                float cosA = cosf(ang), sinA = sinf(ang);
                // Per-bundle length jitter (±15%) and time-varying wobble
                // so the silhouette never looks frozen — real MTs undergo
                // dynamic instability on 10–100s timescales.
                float bundleSeed = (float)f * 3.71f;
                float lenJitter = 0.85f + 0.15f * sinf(bundleSeed);
                float wobbleT = 1.0f + 0.10f * sinf(tNow * 1.7f + bundleSeed);
                // Random bend: shift the bulge off-center per bundle so
                // not every fiber peaks at t=0.5. Also break perfect
                // antipodal symmetry by offsetting the pole-contact angle.
                float bulgeCenter = 0.42f + 0.16f * sinf(bundleSeed * 1.3f);
                float bulgeWidth  = 0.35f + 0.15f * cosf(bundleSeed * 0.9f);
                simd_float3 prev;
                bool hasPrev = false;
                for (int s = 0; s <= BUNDLE_SEGMENTS; s++) {
                    float t = (float)s / (float)BUNDLE_SEGMENTS;
                    // Skewed bell curve for lateral bulge
                    float u = (t - bulgeCenter) / bulgeWidth;
                    float bell = expf(-u * u * 1.6f);
                    float radial = (bundleRPole + (bundleRMax - bundleRPole)
                                    * bell * lenJitter) * wobbleT;
                    float x = clampedPoleA.x + (clampedPoleB.x - clampedPoleA.x) * t;
                    // Per-segment tangential wobble so the fiber doesn't
                    // lie on a perfect plane — adds 3D texture.
                    float tanWobble = 0.015f * cellR *
                                      sinf(tNow * 3.1f + bundleSeed + t * 4.2f);
                    simd_float3 local = {
                        x + tanWobble * 0.3f,
                        cosA * radial + tanWobble * sinA,
                        sinA * radial - tanWobble * cosA
                    };
                    simd_float3 wp = toWorldLocal(local);
                    if (hasPrev) {
                        addStick(allBonds, prev, wp, cs * 0.0022f,
                                 0.25f * spindleAlpha,
                                 0.72f * spindleAlpha,
                                 0.32f * spindleAlpha);
                    }
                    prev = wp;
                    hasPrev = true;
                }
            }
        }

        // ── Chromosomes ─────────────────────────────────────────────
        // Their thickness / arm-length / brightness are all multiplied by
        // the persistent condensation blend. After cytokinesis this blend
        // decays smoothly toward 0, so the rods shrink and dim continuously
        // rather than vanishing in one frame.
        for (int i = 0; i < MitosisState::NUM_CHROMO; i++) {
            auto& chr = gMitosis.chromosomes[i];
            simd_float3 tint = chromosomeSlotColor(i);

            // Attachment-type recoloring — so the viewer can see which
            // chromosomes are still being corrected by Aurora B. Bright
            // green for amphitelic (correct), yellow for monotelic, orange
            // for syntelic, red for merotelic (lagging risk), white for
            // fully unattached.
            float typeR = tint.x, typeG = tint.y, typeB = tint.z;
            if (gMitosis.phase >= MITO_PROMETAPHASE) {
                switch (chr.attachmentType) {
                    case ATTACH_AMPHITELIC: typeR = 0.25f; typeG = 0.95f; typeB = 0.35f; break;
                    case ATTACH_MONOTELIC:  typeR = 0.95f; typeG = 0.85f; typeB = 0.20f; break;
                    case ATTACH_SYNTELIC:   typeR = 0.98f; typeG = 0.55f; typeB = 0.15f; break;
                    case ATTACH_MEROTELIC:  typeR = 0.95f; typeG = 0.20f; typeB = 0.20f; break;
                    case ATTACH_UNATTACHED: typeR = 0.90f; typeG = 0.90f; typeB = 0.90f; break;
                }
                // Blend with per-slot rainbow so neighbors are still
                // distinguishable even when they share a type.
                float blend = 0.65f;
                typeR = typeR * blend + tint.x * (1.0f - blend);
                typeG = typeG * blend + tint.y * (1.0f - blend);
                typeB = typeB * blend + tint.z * (1.0f - blend);
            }

            float effCondense = chr.condensation * chromRenderAlpha;
            float cr = typeR * chromRenderAlpha;
            float cg = typeG * chromRenderAlpha;
            float cb = typeB * chromRenderAlpha;

            if (effCondense > 0.02f) {
                float thick = cs * (0.008f + 0.010f * effCondense);
                float armLen = cs * (0.030f + 0.040f * effCondense);
                simd_float3 wp = toWorldLocal(chr.position);
                simd_float3 sp = toWorldLocal(chr.sisterPosition);
                // Lagging chromatid visualization: in anaphase, merotelic
                // chromatids stay near the midline instead of cleanly
                // reaching their pole (Cimini 2001 JCB). Pull the sister
                // back toward x=0 proportional to separation progress.
                if (gMitosis.phase >= MITO_ANAPHASE && chr.lagging) {
                    float lagPull = gMitosis.chromatidSeparation * 0.85f;
                    sp.x = sp.x * (1.0f - lagPull);
                    wp.x = wp.x * (1.0f - lagPull * 0.60f);
                    // Red flash for lagging chromatids — they're the
                    // literal cause of aneuploidy.
                    cr = fminf(1.0f, cr + 0.35f);
                    cg *= 0.4f;
                    cb *= 0.4f;
                }
                simd_float3 topA = {wp.x, wp.y+armLen, wp.z};
                simd_float3 botA = {wp.x, wp.y-armLen, wp.z};
                simd_float3 topB = {wp.x+armLen*0.3f, wp.y+armLen, wp.z};
                simd_float3 botB = {wp.x-armLen*0.3f, wp.y-armLen, wp.z};
                addSphere(allAtoms, wp, thick, cr, cg, cb);
                addStick(allBonds, topA, botA, thick*0.7f, cr, cg, cb);
                addStick(allBonds, topB, botB, thick*0.7f, cr*0.8f, cg*0.8f, cb*0.8f);
                addSphere(allAtoms, wp, thick*0.5f, 0.9f, 0.9f, 0.1f);
                simd_float3 sTA={sp.x,sp.y+armLen,sp.z}, sBA={sp.x,sp.y-armLen,sp.z};
                simd_float3 sTB={sp.x+armLen*0.3f,sp.y+armLen,sp.z}, sBB={sp.x-armLen*0.3f,sp.y-armLen,sp.z};
                addSphere(allAtoms, sp, thick, cr, cg, cb);
                addStick(allBonds, sTA, sBA, thick*0.7f, cr, cg, cb);
                addStick(allBonds, sTB, sBB, thick*0.7f, cr*0.8f, cg*0.8f, cb*0.8f);
                addSphere(allAtoms, sp, thick*0.5f, 0.9f, 0.9f, 0.1f);
                if (gMitosis.phase >= MITO_PROMETAPHASE && spindleAlpha > 0.08f) {
                    float attachA = chr.attachmentA * spindleAlpha;
                    float attachB = chr.attachmentB * spindleAlpha;
                    // Kinetochore fibers (k-fibers) — the THICK MT bundles
                    // directly connecting each kinetochore to its pole.
                    // Bright yellow-green to stand out from astral (green)
                    // and interpolar (desaturated). These are the "cables"
                    // that actually pull the chromatids apart at anaphase.
                    // Ref: Kiewisz 2022 eLife EM tomography.
                    auto drawKFiber = [&](simd_float3 kineticoreLocal,
                                          simd_float3 poleLocal,
                                          float attach, float sign, int hashSeed) {
                        if (attach <= 0.05f) return;
                        constexpr int SEGS = 6;
                        // K-fibers are the THICKEST bundles in the spindle
                        // (20-30 MTs). Render 2.5× thicker than interpolar.
                        float thick = cs * 0.0055f * (0.80f + attach * 0.50f);
                        // Each k-fiber bows slightly with a unique axis so
                        // the 46 bundles don't all arc the same way.
                        float bowA = (float)hashSeed * 0.39f;
                        float bowAx = cosf(bowA), bowAz = sinf(bowA);
                        float bowMag = cs * 0.05f * (0.8f + attach * 0.4f);
                        simd_float3 prev;
                        bool hasPrev = false;
                        for (int s = 0; s <= SEGS; s++) {
                            float t = (float)s / (float)SEGS;
                            float bow = sinf(t * (float)M_PI) * bowMag * sign;
                            simd_float3 local = {
                                kineticoreLocal.x + (poleLocal.x - kineticoreLocal.x) * t,
                                kineticoreLocal.y + (poleLocal.y - kineticoreLocal.y) * t + bow * 0.3f,
                                kineticoreLocal.z + (poleLocal.z - kineticoreLocal.z) * t
                                    + bow * 0.3f * bowAz
                            };
                            simd_float3 wp2 = toWorldLocal(local);
                            if (hasPrev) {
                                // Bright yellow-green, scales with attachment
                                addStick(allBonds, prev, wp2, thick,
                                         0.70f * attach + 0.20f,
                                         0.95f * attach,
                                         0.25f);
                            }
                            prev = wp2;
                            hasPrev = true;
                        }
                        // Kinetochore "plate" — bright marker at the attachment
                        simd_float3 kineWorld = toWorldLocal(kineticoreLocal);
                        addSphere(allAtoms, kineWorld, thick * 1.6f,
                                  0.95f, 0.95f * attach, 0.20f);
                    };
                    float bowSign = (i % 2 == 0) ? 1.0f : -1.0f;
                    drawKFiber(chr.position,         clampedPoleA, attachA,  bowSign, i);
                    drawKFiber(chr.sisterPosition,   clampedPoleB, attachB, -bowSign, i + 23);
                }
            } else if (chromRenderAlpha > 0.05f) {
                // Only draw the fallback dot while the condensation blend is
                // still visible, so it fades with the rest of the rod.
                addSphere(allAtoms, toWorldLocal(chr.position),
                          cs*0.008f*chromRenderAlpha, cr, cg, cb);
            }
        }

        // ── Contractile ring (visible ring at equator during cytokinesis) ──
        if (gMitosis.furrowDepth > 0.05f) {
            int ringSegs = 24;
            float ringR = cellR * (1.0f - gMitosis.furrowDepth * 0.85f);
            for (int s = 0; s < ringSegs; s++) {
                float a1 = (float)s / ringSegs * 2.0f * M_PI;
                float a2 = (float)(s + 1) / ringSegs * 2.0f * M_PI;
                simd_float3 p1 = {cellPos.x + cosf(a1) * ringR, cellPos.y, cellPos.z + sinf(a1) * ringR};
                simd_float3 p2 = {cellPos.x + cosf(a2) * ringR, cellPos.y, cellPos.z + sinf(a2) * ringR};
                float ringGlow = 0.45f + 0.55f * gMitosis.contractileRingAssembly;
                addStick(allBonds, p1, p2, cs * 0.004f, 0.9f, 0.35f + 0.25f * ringGlow, 0.5f); // actin ring — pink
            }
        }
    }

    // ATP synthase and hemoglobin are rendered through particle system only
    // (no static placement — reduces clutter and prevents escape)

    if (renderPrimaryInterior) {
        auto addRNAChain = [&](const std::string& seq, simd_float3 anchor, simd_float3 axis,
                               float wiggle, float baseRadius,
                               int focusStart, int focusCount, float alpha) {
            if (seq.empty() || alpha <= 0.01f) return;
            int samples = std::min(42, (int)seq.size());
            float denom = fmaxf((float)(samples - 1), 1.0f);
            simd_float3 prev = anchor;
            for (int s = 0; s < samples; s++) {
                float t = (float)s / denom;
                int seqIdx = (int)roundf(t * (float)(seq.size() - 1));
                char base = seq[seqIdx];
                BaseColor bc = getBaseColor(base, true);
                bool focus = focusStart >= 0 && seqIdx >= focusStart && seqIdx < focusStart + focusCount;
                float wave = time * 0.7f + seqIdx * 0.37f;
                simd_float3 tp = {
                    anchor.x + axis.x * t + sinf(wave) * wiggle * (0.35f + 0.65f * t),
                    anchor.y + axis.y * t + cosf(wave * 1.21f) * wiggle * 0.65f,
                    anchor.z + axis.z * t + sinf(wave * 0.83f + 1.2f) * wiggle
                };
                float radius = baseRadius * (focus ? 1.45f : 1.0f);
                float tint = focus ? 1.15f : 0.92f;
                addSphere(allAtoms, tp, radius, bc.r * tint * alpha, bc.g * tint * alpha, bc.b * tint * alpha);
                if (s > 0) {
                    addStick(allBonds, prev, tp, radius * 0.42f,
                             0.72f * alpha, 0.42f * alpha, 0.14f * alpha);
                }
                prev = tp;
            }
        };

        auto findDNAWorldMarker = [&](int bpIndex) -> simd_float3 {
            for (const auto& q : gIntraPhys.particles) {
                if (!q.active || q.type != PT_DNA_NODE || q.stateIndex != bpIndex) continue;
                return particleToWorld(q);
            }
            return cellPos;
        };

        // ── Nascent pre-mRNA trailing active RNA polymerases ──────────
        for (int slot = 0; slot < (int)primaryRNAPolParticles.size(); slot++) {
            const auto& ts = gCDogma.transcription[slot];
            const auto* tx = gCDogma.getTranscript(ts.transcriptIndex);
            if (!ts.active || !tx || tx->preMRNA.empty()) continue;
            simd_float3 wp = particleToWorld(gIntraPhys.particles[primaryRNAPolParticles[slot]]);
            int tailLen = std::min((int)tx->preMRNA.size(), 42);
            int newestLen = std::min(tailLen, 6);
            std::string tail = tx->preMRNA.substr((int)tx->preMRNA.size() - tailLen, tailLen);
            simd_float3 trailAxis = {0.0f, -cellSize * 0.12f, -cellSize * 0.05f};
            addRNAChain(tail, wp, trailAxis,
                        cellSize * 0.015f, cellSize * 0.0048f,
                        tailLen - newestLen, newestLen, 1.0f);
        }

        // ── Replication forks, proofreading pauses, and mismatch-repair foci ──
        for (int slot = 0; slot < (int)primaryDNAPolParticles.size(); slot++) {
            const auto* fork = gCDogma.getReplicationFork(slot);
            if (!fork) continue;
            simd_float3 wp = particleToWorld(gIntraPhys.particles[primaryDNAPolParticles[slot]]);
            float stress = clampf(fork->stress, 0.0f, 1.0f);
            addSphere(allAtoms, wp, cellSize * 0.0062f,
                      0.45f + 0.35f * stress,
                      0.75f + 0.20f * fork->recruitment,
                      1.0f);
            if (fork->proofreading) {
                addSphere(allAtoms, wp, cellSize * 0.0084f,
                          1.0f, 0.75f, 0.28f);
            }
        }
        for (int i = 0; i < CDOGMA_MAX_REPL_ERRORS; i++) {
            const auto& err = gCDogma.replicationErrors[i];
            if (!err.active) continue;
            simd_float3 marker = findDNAWorldMarker(err.bpPosition);
            float pulse = 0.72f + 0.28f * sinf(time * 7.0f + err.bpPosition * 0.17f);
            float mr = cellSize * (err.mmrDetected ? 0.0065f : 0.0052f);
            float cr = err.mmrDetected ? 1.0f : 0.98f;
            float cg = err.mmrDetected ? 0.72f : 0.28f;
            float cb = err.mmrDetected ? 0.22f : 0.22f;
            addSphere(allAtoms, marker, mr * pulse, cr, cg, cb);
        }

        // ── Mature mRNA pool: capped, spliced, polyadenylated transcripts ──
        for (int txIdx = 0; txIdx < CDOGMA_MAX_TRANSCRIPTS; txIdx++) {
            const auto& tx = gCDogma.transcripts[txIdx];
            if (!tx.active || tx.matureMRNA.empty()) continue;

            int focusBase = -1;
            for (const auto& tr : gCDogma.translation) {
                if (tr.active && tr.transcriptIndex == txIdx) {
                    focusBase = gCDogma.startCodonBase + tr.currentCodon * 3;
                    break;
                }
            }

            float orbital = time * 0.35f + tx.uid * 0.91f;
            float exportBlend = tx.exported ? 1.0f : tx.exportProgress;
            simd_float3 anchor = {
                cellPos.x + cosf(orbital) * cellSize * (0.06f + 0.12f * exportBlend),
                cellPos.y + (txIdx - 2.5f) * cellSize * 0.018f,
                cellPos.z + sinf(orbital * 0.9f) * cellSize * (0.04f + 0.10f * exportBlend)
            };
            simd_float3 axis = tx.exported
                ? simd_float3{cellSize * 0.20f, 0.0f, cellSize * 0.05f}
                : simd_float3{cellSize * 0.06f, cellSize * 0.11f, 0.0f};
            std::string seq = tx.exported ? (tx.matureMRNA + tx.polyATail) : tx.preMRNA;
            addRNAChain(seq, anchor, axis,
                        tx.exported ? cellSize * 0.010f : cellSize * 0.014f,
                        tx.exported ? cellSize * 0.0037f : cellSize * 0.0042f,
                        focusBase, 3, tx.exported ? 0.95f : 0.60f);

            if (tx.capped) {
                addSphere(allAtoms, anchor, cellSize * 0.006f, 0.2f, 0.7f, 1.0f);
            }
        }

        // ── Translation: codon/anticodon decoding + nascent peptides ──
        for (int slot = 0; slot < (int)primaryRibosomeParticles.size(); slot++) {
            const auto& tr = gCDogma.translation[slot];
            if (!tr.active) continue;

            simd_float3 ribW = particleToWorld(gIntraPhys.particles[primaryRibosomeParticles[slot]]);
            const auto* tx = gCDogma.getTranscript(tr.transcriptIndex);
            if (tx && !tx->matureMRNA.empty()) {
                int codonBase = gCDogma.startCodonBase + tr.currentCodon * 3;
                int ribbonStart = std::max(0, codonBase - 9);
                int ribbonEnd = std::min((int)tx->matureMRNA.size(), codonBase + 15);
                if (ribbonEnd > ribbonStart) {
                    simd_float3 ribbonAnchor = {ribW.x - cellSize * 0.12f, ribW.y - cellSize * 0.018f, ribW.z};
                    simd_float3 ribbonAxis = {cellSize * 0.24f, 0.0f, 0.0f};
                    addRNAChain(tx->matureMRNA.substr(ribbonStart, ribbonEnd - ribbonStart),
                                ribbonAnchor, ribbonAxis,
                                cellSize * 0.004f, cellSize * 0.0036f,
                                codonBase - ribbonStart, 3, 1.0f);
                }
            }

            if (tr.codon.size() == 3) {
                for (int b = 0; b < 3; b++) {
                    BaseColor codonColor = getBaseColor(tr.codon[b], true);
                    BaseColor antiColor = getBaseColor(tr.anticodon[b], true);
                    simd_float3 codonPos = {
                        ribW.x + (b - 1) * cellSize * 0.012f,
                        ribW.y - cellSize * 0.010f,
                        ribW.z - cellSize * 0.024f
                    };
                    simd_float3 antiPos = {
                        codonPos.x,
                        codonPos.y - cellSize * (0.016f + 0.014f * (1.0f - tr.codonProgress)),
                        codonPos.z + cellSize * 0.014f
                    };
                    addSphere(allAtoms, codonPos, cellSize * 0.0048f,
                              codonColor.r, codonColor.g, codonColor.b);
                    addSphere(allAtoms, antiPos, cellSize * 0.0044f,
                              antiColor.r, antiColor.g, antiColor.b);
                    addStick(allBonds, codonPos, antiPos, cellSize * 0.0018f,
                             0.85f, 0.8f, 0.4f);
                }
            }

            if (tr.incomingAA != '?' && tr.codonProgress < 0.999f) {
                AminoAcid aa = aminoAcidByLetter(tr.incomingAA);
                simd_float3 cargo = {
                    ribW.x - cellSize * 0.03f,
                    ribW.y - cellSize * (0.07f - 0.05f * tr.codonProgress),
                    ribW.z + cellSize * 0.01f
                };
                addSphere(allAtoms, cargo, cellSize * 0.0065f, aa.r, aa.g, aa.b);
            }

            int chainSamples = std::min(42, (int)tr.nascentPeptide.size());
            if (chainSamples > 0) {
                float denom = fmaxf((float)(chainSamples - 1), 1.0f);
                simd_float3 prev = ribW;
                for (int s = 0; s < chainSamples; s++) {
                    float t = (float)s / denom;
                    int aaIdx = (int)roundf(t * (float)(tr.nascentPeptide.size() - 1));
                    AminoAcid aa = aminoAcidByLetter(tr.nascentPeptide[aaIdx]);
                    simd_float3 cp = {
                        ribW.x + t * cellSize * 0.17f +
                            sinf(t * 9.0f + time * 0.5f + slot) * cellSize * 0.007f,
                        ribW.y - cellSize * 0.018f - t * cellSize * 0.11f,
                        ribW.z + cosf(t * 6.0f + time * 0.4f) * cellSize * 0.010f
                    };
                    addSphere(allAtoms, cp, cellSize * 0.0052f, aa.r, aa.g, aa.b);
                    addStick(allBonds, prev, cp, cellSize * 0.002f, 0.55f, 0.55f, 0.55f);
                    prev = cp;
                }
            }
        }

        // ── Released proteins: folding and maturation after termination ──
        int proteinOrdinal = 0;
        for (const auto& protein : gCDogma.proteins) {
            if (!protein.active || protein.aminoAcids.empty()) continue;
            const MoleculeData* productMol = nullptr;
            if (!protein.structureAsset.empty()) {
                productMol = gProtCache.get(protein.structureAsset);
            } else if ((int)protein.aminoAcids.size() >= 120 && (int)protein.aminoAcids.size() <= 170) {
                productMol = gProtCache.get("hbb_folded.pdb");
            }
            int samples = std::min(36, (int)protein.aminoAcids.size());
            float phase = time * 0.24f + proteinOrdinal * 0.93f;
            float orbit = cellSize * (0.16f + 0.03f * (proteinOrdinal % 4));
            simd_float3 center = {
                cellPos.x + cosf(phase) * orbit,
                cellPos.y + ((proteinOrdinal % 3) - 1) * cellSize * 0.06f,
                cellPos.z + sinf(phase * 0.85f) * orbit
            };

            if (productMol) {
                float matureTint = protein.mature ? 1.12f : 0.96f;
                float foldBlend = fmaxf(0.0f, fminf(protein.foldProgress, 1.0f));
                int stride = (productMol->atoms.size() > 1400) ? 3 : ((productMol->atoms.size() > 700) ? 2 : 1);
                addProteinBackbone(allAtoms, allBonds, *productMol, center, cellSize,
                                   0.88f * matureTint, 0.34f * matureTint, 0.28f * matureTint,
                                   phase * 0.65f, stride, foldBlend);
                if (!protein.folded) {
                    float cloudAlpha = 1.0f - foldBlend;
                    addSphere(allAtoms, center, cellSize * (0.010f + 0.008f * cloudAlpha),
                              0.95f * cloudAlpha, 0.48f * cloudAlpha, 0.24f * cloudAlpha);
                }
                proteinOrdinal++;
                continue;
            }

            float compact = protein.foldProgress;
            float denom = fmaxf((float)(samples - 1), 1.0f);
            simd_float3 prev = center;
            for (int s = 0; s < samples; s++) {
                float t = (float)s / denom;
                int aaIdx = (int)roundf(t * (float)(protein.aminoAcids.size() - 1));
                AminoAcid aa = aminoAcidByLetter(protein.aminoAcids[aaIdx]);
                float helixTheta = t * 6.0f * M_PI + phase;
                float radial = cellSize * (0.05f * (1.0f - compact) + 0.018f * compact);
                float axial = cellSize * (0.15f * (1.0f - compact) + 0.025f * compact);
                simd_float3 pp = {
                    center.x + cosf(helixTheta) * radial,
                    center.y + (t - 0.5f) * axial,
                    center.z + sinf(helixTheta) * radial
                };
                float matureTint = protein.mature ? 1.15f : 0.95f;
                addSphere(allAtoms, pp, cellSize * 0.0048f * (protein.folded ? 1.1f : 1.0f),
                          aa.r * matureTint, aa.g * matureTint, aa.b * matureTint);
                if (s > 0) {
                    addStick(allBonds, prev, pp, cellSize * 0.0018f,
                             0.60f * matureTint, 0.58f * matureTint, 0.50f * matureTint);
                }
                prev = pp;
            }
            proteinOrdinal++;
        }
    }

    // No render-time escape clamp — the physics membrane handles containment.

    // ── Render interiors for non-primary cells ────────────────────────
    // Only simulate/render the nearest MAX_ACTIVE_INTERIORS cells for performance.
    // All cells get DNA/particles but distant cells are paused.
    {
        int nCells = (int)gSim.cells.size();
        constexpr int MAX_ACTIVE_INTERIORS = 3; // Max non-primary cells to simulate

        // Find nearest cells to camera
        simd_float3 camPos = gCamera.getPosition();
        std::vector<std::pair<float, int>> cellDists;
        for (int ci = 1; ci < nCells; ci++) {
            auto& c = gSim.cells[ci];
            if (!c.alive) continue;
            float dx = c.position.x - camPos.x;
            float dy = c.position.y - camPos.y;
            float dz = c.position.z - camPos.z;
            cellDists.push_back({dx*dx + dy*dy + dz*dz, ci});
        }
        std::sort(cellDists.begin(), cellDists.end());

        int activeCount = 0;
        for (auto& [dist, ci] : cellDists) {
            if (!shouldRenderSourceCellIndex(ci)) continue;
            auto& c = gSim.cells[ci];
            auto& interior = gCellInteriors[c.cellUid];
            float cSize = c.radius * c.size;

            // Initialize if needed
            if (!interior.initialized) {
                initCellInterior(interior, cSize);
            }

            if (c.program.cdogmaInitialized) {
                syncGenomeReplicationVisuals(interior.phys, c.program.cdogma, cSize, false);
            }

            interior.postMitoticCondense = fmaxf(0.0f, fminf(1.0f, c.postDivisionRecovery / 6.0f));
            if (interior.postMitoticCondense > 0.02f) {
                restoreInterphaseNuclearLayout(interior.phys, cSize, interior.postMitoticCondense);
            }

            // Only step physics for nearest cells
            bool active = (activeCount < MAX_ACTIVE_INTERIORS);
            if (active) {
                activeCount++;
                float cR = cSize;
                float nucR_ci = cR * 0.45f;
                interior.phys.params.cellRadius = cR;
                interior.phys.params.cellVelocity = {0, 0, 0};
                for (auto& p : interior.phys.particles) {
                    if (!p.active) continue;
                    p.confineCenter = {0, 0, 0};
                    if (isNuclearParticleType(p.type)) {
                        p.confineRadius = nucR_ci;
                    } else {
                        p.confineRadius = cR;
                    }
                }
                driveInteriorActivity(interior.phys, cR, nucR_ci, c.cellUid, time);
                interior.phys.skipCollisions = true; // Always skip for secondary cells
                interior.phys.step(simDt);
            }

            // Render particles → allAtoms/allBonds (local → world)
            simd_float3 cp = c.position;
            for (int pi = 0; pi < (int)interior.phys.particles.size(); pi++) {
                auto& p = interior.phys.particles[pi];
                if (!p.active) continue;
                simd_float3 wp = {p.position.x + cp.x, p.position.y + cp.y, p.position.z + cp.z};
                float renderRadius = p.radius;
                float renderCr = p.colorR, renderCg = p.colorG, renderCb = p.colorB;
                if (p.type == PT_DNA_NODE && interior.postMitoticCondense > 0.02f) {
                    int slot = mitosisChromosomeSlot(p);
                    simd_float3 tint = chromosomeSlotColor(slot);
                    float condense = interior.postMitoticCondense;
                    float slotSpan = fmaxf((float)HBB_LENGTH / (float)MitosisState::NUM_CHROMO, 1.0f);
                    float slotOffset = fmodf((float)fmax(p.stateIndex, 0), slotSpan);
                    float slotLocal = (slotOffset + 0.5f) / slotSpan - 0.5f;
                    float centromere = 1.0f - fminf(fabsf(slotLocal) * 2.2f, 1.0f);
                    renderRadius *= 1.0f + condense * 1.20f + centromere * 0.18f;
                    renderCr = p.colorR * (1.0f - condense * 0.35f) + tint.x * (0.45f + condense * 0.55f);
                    renderCg = p.colorG * (1.0f - condense * 0.35f) + tint.y * (0.45f + condense * 0.55f);
                    renderCb = p.colorB * (1.0f - condense * 0.35f) + tint.z * (0.45f + condense * 0.55f);
                    renderCr = fminf(renderCr + centromere * 0.18f, 1.0f);
                    renderCg = fminf(renderCg + centromere * 0.18f, 1.0f);
                }
                const char* molId = (p.type == PT_MOLECULE) ? particleTagToMolecule(p.tag) : nullptr;
                const MoleculeData* smallMol = molId ? gMolCache.get(molId) : nullptr;
                const char* bgProtFile = particleToProtein(p.type);
                const MoleculeData* bgProtMol = (bgProtFile && p.type != PT_DNA_NODE)
                                                ? gProtCache.get(bgProtFile) : nullptr;
                if (smallMol) {
                    float minHalfExtentFrac = fminf(0.016f, (renderRadius * 1.35f) / fmaxf(c.radius, 0.001f));
                    float molScale = chemistryScaleFromRealSize(*smallMol, c.radius, 18.0f, minHalfExtentFrac);
                    int stride = (smallMol->atoms.size() > 96) ? 2 : 1;
                    addMoleculeGeometry(allAtoms, allBonds, *smallMol, wp, molScale,
                                        renderCr, renderCg, renderCb, (float)pi * 1.37f, stride,
                                        fmaxf(renderRadius * 0.20f, c.radius * 0.0008f),
                                        0.38f, 0.56f, 0.20f, 0.0f, 0.0f,
                                        c.cellUid * 97 + pi * 17);
                } else if (bgProtMol) {
                    // Background proteins render atoms at a heavy stride
                    // so 600+ cells × dozens of proteins stays tractable.
                    int totalAtoms = (int)bgProtMol->atoms.size();
                    int stride = (totalAtoms > 1500) ? 16 : (totalAtoms > 500) ? 8 : 4;
                    addProteinBackbone(allAtoms, allBonds, *bgProtMol, wp, c.radius,
                                       p.colorR, p.colorG, p.colorB,
                                       (float)pi * 1.37f, stride);
                } else {
                    addSphere(allAtoms, wp, renderRadius, renderCr, renderCg, renderCb);
                }

                if (p.type == PT_VESICLE_COPII || p.type == PT_VESICLE_SECRETORY) {
                    addSphere(allAtoms, wp, p.radius * 0.4f, 0.4f, 0.3f, 0.6f);
                }

                // DNA backbone bonds
                if (p.type == PT_DNA_NODE && p.tetherId >= 0 &&
                    p.tetherId < (int)interior.phys.particles.size()) {
                    auto& other = interior.phys.particles[p.tetherId];
                    if (other.active) {
                        simd_float3 owp = {other.position.x + cp.x, other.position.y + cp.y, other.position.z + cp.z};
                        float bondRadius = p.radius * 0.4f;
                        float bondCr = 0.7f, bondCg = 0.5f, bondCb = 0.2f;
                        if (interior.postMitoticCondense > 0.02f) {
                            simd_float3 tint = chromosomeSlotColor(mitosisChromosomeSlot(p));
                            float condense = interior.postMitoticCondense;
                            bondRadius *= 1.0f + condense * 1.10f;
                            bondCr = tint.x; bondCg = tint.y; bondCb = tint.z;
                        }
                        addStick(allBonds, wp, owp, bondRadius, bondCr, bondCg, bondCb);
                    }
                }
                // Base pair rungs
                if (p.type == PT_DNA_NODE && pi + 1 < (int)interior.phys.particles.size()) {
                    auto& next = interior.phys.particles[pi + 1];
                    if (next.active && next.type == PT_DNA_NODE && next.stateIndex == p.stateIndex) {
                        simd_float3 nwp = {next.position.x + cp.x, next.position.y + cp.y, next.position.z + cp.z};
                        float rungRadius = p.radius * 0.25f;
                        float rungCr = 0.4f, rungCg = 0.4f, rungCb = 0.5f;
                        if (interior.postMitoticCondense > 0.02f) {
                            simd_float3 tint = chromosomeSlotColor(mitosisChromosomeSlot(p));
                            float condense = interior.postMitoticCondense;
                            rungRadius *= 1.0f + condense * 0.85f;
                            rungCr = 0.75f * tint.x + 0.20f;
                            rungCg = 0.75f * tint.y + 0.20f;
                            rungCb = 0.75f * tint.z + 0.20f;
                        }
                        addStick(allBonds, wp, nwp, rungRadius, rungCr, rungCg, rungCb);
                    }
                }
            }
        }
    }

    // Print atom/bond count every 60 frames to catch growth
    {
        static int gpuFrame = 0;
        if (++gpuFrame % 60 == 0 || allAtoms.size() > 100000) {
            printf("[GPU] frame=%d atoms=%zu bonds=%zu particles=%d cells=%d\n",
                   gpuFrame, allAtoms.size(), allBonds.size(),
                   (int)gIntraPhys.particles.size(), (int)gSim.cells.size());
        }
    }
    // Safety cap: 800k atoms / 300k bonds fits the ribosome atomic
    // detail (~600k) + membrane mosaic (~12k) with headroom. 1.5M hit
    // a GPU buffer edge case that crashed around frame 115.
    constexpr size_t MAX_GPU_ATOMS = 800000;
    constexpr size_t MAX_GPU_BONDS = 300000;
    if (allAtoms.size() > MAX_GPU_ATOMS) {
        static bool warned = false;
        if (!warned) {
            printf("[PERF WARNING] allAtoms=%zu > cap %zu — truncating. "
                   "Particles=%d cells=%d\n",
                   allAtoms.size(), MAX_GPU_ATOMS,
                   (int)gIntraPhys.particles.size(), (int)gSim.cells.size());
            warned = true;
        }
        allAtoms.resize(MAX_GPU_ATOMS);
    }
    if (allBonds.size() > MAX_GPU_BONDS) {
        allBonds.resize(MAX_GPU_BONDS);
    }

    // Upload to GPU
    if (!allAtoms.empty()) {
        size_t sz = allAtoms.size() * sizeof(GPUAtomInstance);
        if (!gMolRender.atomBuffer || gMolRender.atomBuffer.length < sz)
            gMolRender.atomBuffer = [gCtx.device() newBufferWithLength:std::max(sz, (size_t)256)
                                                              options:MTLResourceStorageModeShared];
        memcpy(gMolRender.atomBuffer.contents, allAtoms.data(), sz);
    }
    if (!allBonds.empty()) {
        size_t sz = allBonds.size() * sizeof(GPUBondInstance);
        if (!gMolRender.bondBuffer || gMolRender.bondBuffer.length < sz)
            gMolRender.bondBuffer = [gCtx.device() newBufferWithLength:std::max(sz, (size_t)256)
                                                              options:MTLResourceStorageModeShared];
        memcpy(gMolRender.bondBuffer.contents, allBonds.data(), sz);
    }

    gMolRender.atomCount = (int)allAtoms.size();
    gMolRender.bondCount = (int)allBonds.size();
}

// ── Switch simulation mode ──────────────────────────────────────────────
static void switchMode(SimMode newMode) {
    gSimMode = newMode;
    gSim.init(newMode);
    initFluidHaze();                // dark-blue translucent fluid backdrop
    initMediumChemicals();          // re-seed the volumetric chemical swarm
    gSelectedCell = 0;
    gSelectedCellUid = gSim.cells.empty() ? -1 : gSim.cells[0].cellUid;
    gPostMitosisPairCameraTimer = 0.0f;
    gPostMitosisPairUidA = -1;
    gPostMitosisPairUidB = -1;
    gIntraInitialized = false;
    gCellInteriors.clear(); // Clear all per-cell particle systems
    gOrganelleMotionStates.clear();
    initializeDogmaStatesFromSimulation();
    gMolRender.clear();
    if (newMode == MODE_SINGLE_CELL)
        gCamera.setSingleCellView();
    else
        gCamera.setColonyView();
}

// Apply gSetup to the running simulation. Replaces the existing cells
// with a freshly-seeded population and bakes the medium with the user's
// chosen initial concentrations. Called when the user clicks "Start" on
// the setup overlay, or whenever they re-open it and re-start.
static void applySetupAndStartSimulation() {
    gSim.timeScale = gSetup.initialTimeScale;
    gSim.init(MODE_SINGLE_CELL);
    // Seed additional cells to match gSetup.initCells. Place each cell
    // at its proper rest-on-floor Y (FLOOR_Y + radius × size × 0.85) so
    // the lower hemisphere doesn't clip through the substrate on the
    // first frame before physics has a chance to settle it.
    while ((int)gSim.cells.size() < gSetup.initCells
           && (int)gSim.cells.size() < MAX_CELLS) {
        SimCell c;
        float r = sqrtf((float)rand()/RAND_MAX) * SCENE_BOUND * 0.45f;
        float a = (float)rand()/RAND_MAX * 2.0f * (float)M_PI;
        // Temporary y; overwritten after init() reads radius/size.
        simd_float3 p = {r*cosf(a), FLOOR_Y, r*sinf(a)};
        c.init(p, (int)gSim.cells.size());
        c.cellUid = gSim.allocateCellUid();
        // Rest cell on the floor — matches updatePhysics anchor.
        c.position.y = FLOOR_Y + c.radius * c.size * 0.85f;
        gSim.cells.push_back(c);
    }
    // Override the medium-field initial concentrations with the user's
    // template values. Existing values are already set from Constants.h
    // defaults; we overwrite only the species the UI exposes.
    for (int i = 0; i < MediumField::CELLS; i++) {
        gSim.nutrients.c[MS_GLUCOSE][i]   = gSetup.glucoseMM;
        gSim.nutrients.c[MS_GLUTAMINE][i] = gSetup.glutamineMM;
        gSim.nutrients.c[MS_AA_POOL][i]   = gSetup.aaPoolMM;
        gSim.nutrients.c[MS_O2][i]        = gSetup.o2MM;
        gSim.nutrients.c[MS_CO2][i]       = gSetup.co2MM;
        gSim.nutrients.c[MS_HPLUS][i]     = gSetup.pH;
        gSim.nutrients.c[MS_GROWTH_F][i]  = gSetup.growthFactorNgML;
    }
    // Reset mass-balance baseline so the checker doesn't flag the new
    // totals as drift.
    for (int s = 0; s < MS_COUNT; s++) {
        double sum = 0.0;
        for (int i = 0; i < MediumField::CELLS; i++) sum += gSim.nutrients.c[s][i];
        gSim.nutrients.totalAtStart[s] = sum;
    }
    initFluidHaze();
    initMediumChemicals();
    gSelectedCell = 0;
    gSelectedCellUid = gSim.cells.empty() ? -1 : gSim.cells[0].cellUid;
    gCamera.setSingleCellView();
}

// Apply one of the named templates to gSetup.
static void applySetupTemplate(int idx) {
    gSetup.templateIdx = idx;
    switch (idx) {
        case 0:   // HeLa Standard (DMEM HG + 10 % FBS, 5 % CO2 atmosphere)
            gSetup.initCells = 5;
            gSetup.glucoseMM = 25.0f;  gSetup.glutamineMM = 4.0f;
            gSetup.aaPoolMM = 5.0f;    gSetup.o2MM = 0.20f;
            gSetup.co2MM = 1.20f;      gSetup.pH = 7.40f;
            gSetup.growthFactorNgML = 50.0f;
            gSetup.initialTimeScale = 20.0f;
            break;
        case 1:   // CTC Fluo-N2DL-HeLa (Mitocheck, 10 % CO2, GlutaMAX + NEAA)
            gSetup.initCells = 125;
            gSetup.glucoseMM = 25.0f;  gSetup.glutamineMM = 6.0f;
            gSetup.aaPoolMM = 7.0f;    gSetup.o2MM = 0.20f;
            gSetup.co2MM = 2.40f;      gSetup.pH = 7.40f;
            gSetup.growthFactorNgML = 50.0f;
            gSetup.initialTimeScale = 60.0f;
            break;
        case 2:   // Hypoxia (3 % O2, low glucose, stressful)
            gSetup.initCells = 10;
            gSetup.glucoseMM = 5.5f;   gSetup.glutamineMM = 2.0f;
            gSetup.aaPoolMM = 3.0f;    gSetup.o2MM = 0.04f;    // 3 % atm
            gSetup.co2MM = 1.20f;      gSetup.pH = 7.20f;
            gSetup.growthFactorNgML = 10.0f;
            gSetup.initialTimeScale = 20.0f;
            break;
        case 3:   // Starvation (low glucose + low AAs — drives apoptosis)
            gSetup.initCells = 5;
            gSetup.glucoseMM = 1.0f;   gSetup.glutamineMM = 0.5f;
            gSetup.aaPoolMM = 1.0f;    gSetup.o2MM = 0.20f;
            gSetup.co2MM = 1.20f;      gSetup.pH = 7.40f;
            gSetup.growthFactorNgML = 5.0f;
            gSetup.initialTimeScale = 20.0f;
            break;
        case 4:   // Rich Medium (overfed — serum-rich IGF/EGF boost)
            gSetup.initCells = 5;
            gSetup.glucoseMM = 35.0f;  gSetup.glutamineMM = 8.0f;
            gSetup.aaPoolMM = 10.0f;   gSetup.o2MM = 0.25f;
            gSetup.co2MM = 1.20f;      gSetup.pH = 7.40f;
            gSetup.growthFactorNgML = 150.0f;
            gSetup.initialTimeScale = 20.0f;
            break;
        case 5:   // Single Cell Detail (classic 1-cell mode)
            gSetup.initCells = 1;
            gSetup.glucoseMM = 25.0f;  gSetup.glutamineMM = 4.0f;
            gSetup.aaPoolMM = 5.0f;    gSetup.o2MM = 0.20f;
            gSetup.co2MM = 1.20f;      gSetup.pH = 7.40f;
            gSetup.growthFactorNgML = 50.0f;
            gSetup.initialTimeScale = 1.0f;
            break;
        case 6:   // Dense Colony (near-confluent 400 cells)
            gSetup.initCells = 400;
            gSetup.glucoseMM = 25.0f;  gSetup.glutamineMM = 4.0f;
            gSetup.aaPoolMM = 5.0f;    gSetup.o2MM = 0.20f;
            gSetup.co2MM = 1.20f;      gSetup.pH = 7.40f;
            gSetup.growthFactorNgML = 50.0f;
            gSetup.initialTimeScale = 60.0f;
            break;
    }
}

// ══════════════════════════════════════════════════════════════════════════
//  Complete-state export — one-click dump of every simulation artifact
//  to a timestamped folder on disk.
// ══════════════════════════════════════════════════════════════════════════

// Last-export status, surfaced in the Export panel.
static std::string gLastExportPath;
static double      gLastExportWallTime = 0.0;

// Per-grid-cell concentration dump. One row per grid cell, one column per species.
static void exportMediumFieldCSV(const std::string& path, const Simulation& sim) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return;
    fprintf(f, "ix,iz,world_x_wu,world_z_wu");
    for (int s = 0; s < MS_COUNT; s++)
        fprintf(f, ",%s_%s", mediumSpeciesName(s), mediumSpeciesUnit(s));
    fprintf(f, ",extracellular_protease\n");
    const int N = MediumField::N;
    for (int iz = 0; iz < N; iz++) {
        for (int ix = 0; ix < N; ix++) {
            int idx = iz * N + ix;
            float gx = (float)ix / (N - 1) * 2 * SCENE_BOUND - SCENE_BOUND;
            float gz = (float)iz / (N - 1) * 2 * SCENE_BOUND - SCENE_BOUND;
            fprintf(f, "%d,%d,%.3f,%.3f", ix, iz, gx, gz);
            for (int s = 0; s < MS_COUNT; s++)
                fprintf(f, ",%.4f", sim.nutrients.c[s][idx]);
            fprintf(f, ",%.4f\n", sim.nutrients.extracellularProteaseLevel[idx]);
        }
    }
    fclose(f);
}

// Dump every apoptotic body's pose, mass, and phase state.
static void exportApoptoticBodiesCSV(const std::string& path) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return;
    fprintf(f, "idx,origin_cell_uid,kind,phase,pos_x,pos_y,pos_z,"
               "radius_wu,radius0_wu,remaining_biomass_bm,remaining_membrane_bm,"
               "remaining_receptor_bm,initial_biomass_bm,age_biosec,"
               "osmotic_swell,membrane_integrity\n");
    for (size_t i = 0; i < gApoBodies.size(); i++) {
        const ApoptoticBody& b = gApoBodies[i];
        const char* phaseN[] = {"DRIFT","SWELL","BURST"};
        fprintf(f, "%zu,%d,%u,%s,%.3f,%.3f,%.3f,"
                   "%.4f,%.4f,%.4f,%.4f,"
                   "%.4f,%.4f,%.2f,"
                   "%.4f,%.4f\n",
                i, b.originCellUid, (unsigned)b.kind,
                (b.phase < 3 ? phaseN[b.phase] : "?"),
                b.position.x, b.position.y, b.position.z,
                b.radius, b.radius0,
                b.remainingBiomass, b.remainingMembrane,
                b.remainingReceptor, b.initialBiomass,
                b.ageBioSec,
                b.osmoticSwell, b.membraneIntegrity);
    }
    fclose(f);
}

// Dump every medium-chemical particle (visual swarm tracked per-state).
static void exportMediumChemicalsCSV(const std::string& path) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return;
    fprintf(f, "idx,species,species_name,state,"
               "pos_x,pos_y,pos_z,home_x,home_y,home_z,"
               "target_cell_idx,state_timer_biosec,phase\n");
    const char* stateN[] = {"FREE","ATTRACTED","BINDING","TRANSPORT","DESPAWN"};
    for (size_t i = 0; i < gMediumChemicals.size(); i++) {
        const MediumChemical& mc = gMediumChemicals[i];
        fprintf(f, "%zu,%d,%s,%s,"
                   "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"
                   "%d,%.2f,%.4f\n",
                i, mc.species, mediumSpeciesName(mc.species),
                (mc.state < 5 ? stateN[mc.state] : "?"),
                mc.position.x, mc.position.y, mc.position.z,
                mc.home.x, mc.home.y, mc.home.z,
                mc.targetCellIdx, mc.stateTimer, mc.phase);
    }
    fclose(f);
}

// Run summary + closed-system mass-balance drift for each species.
static void exportManifest(const std::string& path, const Simulation& sim) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) return;
    time_t now = ::time(nullptr);
    char tbuf[64];
    std::strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S",
                  std::localtime(&now));
    fprintf(f, "CellSim Export Manifest\n");
    fprintf(f, "=======================\n");
    fprintf(f, "wall_time          : %s\n", tbuf);
    fprintf(f, "session            : %s\n", gRunSessionTag.c_str());
    fprintf(f, "mode               : %s\n",
            sim.mode == MODE_SINGLE_CELL ? "SINGLE_CELL" : "COLONY");
    fprintf(f, "bio_seconds        : %.0f\n", sim.bioTime);
    fprintf(f, "bio_hours          : %.3f\n", sim.bioTime / 3600.0f);
    fprintf(f, "time_scale         : %.2f\n", sim.timeScale);
    fprintf(f, "\n-- Population --\n");
    fprintf(f, "cells_alive        : %d\n", sim.statAlive);
    fprintf(f, "cells_prolif       : %d\n", sim.statProlif);
    fprintf(f, "cells_quiescent    : %d\n", sim.statQuiescent);
    fprintf(f, "cells_apoptotic    : %d\n", sim.statApoptotic);
    fprintf(f, "cells_necrotic     : %d\n", sim.statNecrotic);
    fprintf(f, "total_divisions    : %d\n", sim.statDivisions);
    fprintf(f, "total_deaths       : %d\n", sim.statDeaths);
    fprintf(f, "apoptotic_bodies   : %zu\n", gApoBodies.size());
    fprintf(f, "medium_chemicals   : %zu\n", gMediumChemicals.size());
    fprintf(f, "\n-- Phase distribution --\n");
    float total = fmaxf(1.0f, (float)sim.statAlive);
    const char* phaseN[4] = {"G1","S","G2","M"};
    for (int p = 0; p < 4; p++)
        fprintf(f, "phase_%s           : %d (%.1f%%)\n",
                phaseN[p], sim.statPhases[p],
                sim.statPhases[p] / total * 100.0f);
    fprintf(f, "\n-- Medium (dish-mean, closed system) --\n");
    for (int s = 0; s < MS_COUNT; s++)
        fprintf(f, "%-10s         : %.3f %s\n",
                mediumSpeciesName(s), sim.nutrients.mean(s),
                mediumSpeciesUnit(s));
    fprintf(f, "\n-- Mass-balance drift (should be < 1e-3) --\n");
    double drift[MS_COUNT];
    double maxDrift = sim.nutrients.checkBalance(drift);
    for (int s = 0; s < MS_COUNT; s++)
        fprintf(f, "drift_%-10s    : %+.6e\n",
                mediumSpeciesName(s), drift[s]);
    fprintf(f, "MAX_DRIFT          : %.6e %s\n", maxDrift,
            maxDrift < 1e-3 ? "(OK)" : "(BAD — mass conservation violated)");
    fprintf(f, "\n-- Files in this bundle --\n");
    fprintf(f, "population.csv     : time-series (one row per sample)\n");
    fprintf(f, "cells.csv          : per-cell snapshot at export time\n");
    fprintf(f, "medium_field.csv   : 64×64 grid × 12 species concentrations\n");
    fprintf(f, "apo_bodies.csv     : apoptotic-body ledgers\n");
    fprintf(f, "medium_chems.csv   : floating chemical-particle swarm\n");
    fclose(f);
}

// One-click "export everything". Creates a timestamped subfolder inside
// `destDir` and writes every artifact into it. Updates gLastExportPath so
// the UI can surface the destination and a "Reveal in Finder" button.
static void exportAllData(const std::string& destDir) {
    time_t now = ::time(nullptr);
    struct tm* t = localtime(&now);
    char folder[512];
    snprintf(folder, sizeof(folder),
             "%s/cellsim_export_%04d%02d%02d_%02d%02d%02d",
             destDir.c_str(),
             t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
             t->tm_hour, t->tm_min, t->tm_sec);
    if (mkdir(folder, 0755) != 0 && errno != EEXIST) {
        printf("[Export] Failed to create %s (errno=%d)\n", folder, errno);
        gLastExportPath = std::string("ERROR creating ") + folder;
        return;
    }
    std::string base(folder);
    gTS.exportCSV(base + "/population.csv");
    gTS.exportCellSnapshot(base + "/cells.csv", gSim);
    exportMediumFieldCSV(base + "/medium_field.csv", gSim);
    exportApoptoticBodiesCSV(base + "/apo_bodies.csv");
    exportMediumChemicalsCSV(base + "/medium_chems.csv");
    exportManifest(base + "/manifest.txt", gSim);
    gLastExportPath = base;
    gLastExportWallTime = (double)::time(nullptr);
    printf("[Export] Complete bundle → %s\n", base.c_str());
}

// Capture an in-RAM snapshot of everything needed for export, so the
// actual file-write can run on the Cocoa main thread after the native
// folder-picker closes without racing the sim loop. For now the export
// is fast enough (< 1 s on a 1000-cell run) that we just dispatch the
// whole call; `gSim` / `gTS` / bodies vectors are stable during a paused
// UI tick, matching how the existing CSV/Snapshot export works.
static void exportAllDataWithPicker() {
    // Open an NSOpenPanel that lets the user pick a DIRECTORY. Start the
    // picker at the project's exports/ folder (gExportDir); persist the
    // chosen dir back into gExportDir so the next export defaults to
    // the same place.
    dispatch_async(dispatch_get_main_queue(), ^{
        NSOpenPanel* panel = [NSOpenPanel openPanel];
        [panel setTitle:@"Choose export destination folder"];
        [panel setCanChooseDirectories:YES];
        [panel setCanChooseFiles:NO];
        [panel setAllowsMultipleSelection:NO];
        [panel setCanCreateDirectories:YES];
        [panel setPrompt:@"Export Here"];
        [panel setDirectoryURL:[NSURL fileURLWithPath:
            [NSString stringWithUTF8String:gExportDir.c_str()]]];
        if ([panel runModal] == NSModalResponseOK) {
            NSURL* url = [panel URL];
            std::string chosen = [[url path] UTF8String];
            gExportDir = chosen;                   // remember for next time
            exportAllData(chosen);
        }
    });
}

// ── ImGui UI ────────────────────────────────────────────────────────────
static void drawUI() {
    // ── Setup Overlay (on startup / when user presses R) ─────────────
    // Modal-ish panel at the top-center. User picks a preset template
    // (which fills every field below with the matching values) OR
    // manually tunes each field, then clicks Start.
    if (gSetup.showOverlay) {
        ImVec2 ds = ImGui::GetIO().DisplaySize;
        float winW = 560, winH = 520;
        ImGui::SetNextWindowPos(ImVec2((ds.x - winW) * 0.5f, (ds.y - winH) * 0.5f),
                                ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(winW, winH), ImGuiCond_Always);
        ImGui::Begin("Cell Mode — Setup", nullptr,
                     ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize);

        ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f),
                           "Pick a template or manually configure, then Start");
        ImGui::Separator();

        // ── Preset templates (top row buttons) ──
        ImGui::Text("Templates:");
        const char* TPL_NAMES[] = {
            "HeLa Standard",
            "CTC HeLa (Mitocheck)",
            "Hypoxia 3% O2",
            "Starvation",
            "Rich Medium",
            "Single Cell Detail",
            "Dense Colony (400)"
        };
        const int TPL_COUNT = sizeof(TPL_NAMES) / sizeof(TPL_NAMES[0]);
        for (int i = 0; i < TPL_COUNT; i++) {
            bool active = (i == gSetup.templateIdx);
            if (active) ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.25f, 0.55f, 0.85f, 1.0f));
            if (ImGui::Button(TPL_NAMES[i], ImVec2(175, 0))) {
                applySetupTemplate(i);
            }
            if (active) ImGui::PopStyleColor();
            // 3 buttons per row
            if ((i % 3) != 2 && i < TPL_COUNT - 1) ImGui::SameLine();
        }

        ImGui::Separator();
        ImGui::Text("Configuration:");
        ImGui::SliderInt("Initial cells",      &gSetup.initCells,       1, 500);
        ImGui::SliderFloat("Glucose (mM)",     &gSetup.glucoseMM,       0.0f, 50.0f, "%.2f");
        ImGui::SliderFloat("Glutamine (mM)",   &gSetup.glutamineMM,     0.0f, 10.0f, "%.2f");
        ImGui::SliderFloat("AA pool (mM)",     &gSetup.aaPoolMM,        0.0f, 15.0f, "%.2f");
        ImGui::SliderFloat("O2 (mM)",          &gSetup.o2MM,            0.0f,  0.30f, "%.3f");
        ImGui::SliderFloat("CO2 (mM)",         &gSetup.co2MM,           0.0f,  5.0f,  "%.2f");
        ImGui::SliderFloat("pH",               &gSetup.pH,              6.5f,  8.0f,  "%.2f");
        ImGui::SliderFloat("Growth factors",   &gSetup.growthFactorNgML, 0.0f, 200.0f, "%.1f ng/mL");
        ImGui::SliderFloat("Initial timeScale",&gSetup.initialTimeScale, 0.1f, 100.0f, "%.1f×");

        ImGui::Separator();
        if (ImGui::Button("▶  Start simulation", ImVec2(-1, 36))) {
            applySetupAndStartSimulation();
            gSetup.showOverlay = false;
        }
        ImGui::TextDisabled("(press R anytime to re-open this panel)");
        ImGui::End();
        // Don't draw the rest of the UI while setup is open; the user
        // hasn't picked their conditions yet.
        return;
    }

    // ── Cell Mode runtime controls (top-center) ──────────────────────
    // Replaces the old "Mode [Single Cell | Colony]" switcher. There's
    // now ONE cell mode, configured at startup via the Setup overlay.
    {
        float winW = 240;
        ImGui::SetNextWindowPos(ImVec2((ImGui::GetIO().DisplaySize.x - winW) * 0.5f, 10),
                                ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(winW, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Cell Mode", nullptr, ImGuiWindowFlags_NoCollapse);
        if (ImGui::Button("Re-setup (R)", ImVec2(-1, 0))) {
            gSetup.showOverlay = true;
        }
        activeFocusCellIndex();
        ImGui::Text("Cells: %d", (int)gSim.cells.size());
        ImGui::Checkbox("Follow cell (F)", &gFollowCell);
        ImGui::SameLine();
        ImGui::Checkbox("Solo focus", &gSoloFocusCell);
        if (gFollowCell) {
            ImGui::Text("Zoom: %.2fx  (scroll)", gFollowZoom);
            ImGui::SameLine();
            if (ImGui::SmallButton("Reset zoom")) gFollowZoom = 1.0f;
            ImGui::TextColored(ImVec4(0.5f,0.7f,1.0f,1.0f),
                "WASD = jump to nearest cell in that direction");
        }
        if (gSim.cells.size() > 1) {
            ImGui::SetNextItemWidth(120);
            int prevSelected = gSelectedCell;
            ImGui::SliderInt("Cell index", &gSelectedCell, 0, (int)gSim.cells.size()-1);
            if (gSelectedCell != prevSelected) {
                selectCellIndex(gSelectedCell);
            }
        }
        ImGui::End();
    }

    // ── Drug Lab panel (right side, always visible) ──────────────────────
    // Emergent drug system: pick a compound from the library (loaded from
    // data/bioagents/drugs.csv), set a dish-level concentration, apply it
    // uniformly. The backend runs binding-matcher scoring, pharmacokinetics,
    // and per-target modulators — phenotype is discovered, not prescribed.
    {
        float winW = 310;
        float px = ImGui::GetIO().DisplaySize.x - winW - 10;
        float py = 560.0f;   // below Export Data panel
        ImGui::SetNextWindowPos(ImVec2(px, py), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(winW, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Drug Lab (chemistry-driven)");
        ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f),
                           "MOA-free: structure → emergent effect");
        ImGui::Separator();

        int nDrugs = (int)gBioagents.all().size();
        if (nDrugs == 0) {
            ImGui::TextColored(ImVec4(0.9f, 0.5f, 0.3f, 1.0f),
                               "No drugs loaded. Check data/bioagents/drugs.csv");
        } else {
            static int selected = 0;
            if (selected >= nDrugs) selected = 0;
            std::vector<const char*> names;
            names.reserve(nDrugs);
            for (int i = 0; i < nDrugs; i++) {
                names.push_back(gBioagents.all()[i].id.c_str());
            }
            ImGui::Combo("Compound", &selected, names.data(), nDrugs);

            const ChemicalEntity& d = gBioagents.all()[selected];
            ImGui::TextColored(ImVec4(0.7f,0.75f,0.8f,1.0f),
                               "%s  MW=%.0f  logP=%.2f",
                               d.name.c_str(), d.mw, d.logP);
            ImGui::TextColored(ImVec4(0.7f,0.75f,0.8f,1.0f),
                               "TPSA=%.0f  HBD=%d  HBA=%d  rings=%d",
                               d.tpsa, d.hbd, d.hba, d.aromatic_rings);

            static float logConc_uM = 1.0f;  // log10
            ImGui::SliderFloat("log[C] µM", &logConc_uM, -2.0f, 3.0f, "%.1f");
            float conc = powf(10.0f, logConc_uM);
            ImGui::Text("Concentration: %.3f µM", conc);

            ImGui::Separator();
            if (ImGui::Button("APPLY UNIFORMLY", ImVec2(-1, 30))) {
                applyDrugVisuals(d.id, conc);
            }

            // Show top affinities of the MOST-RECENTLY-APPLIED drug.
            if (!gSim.appliedDrugs.empty()) {
                ImGui::Separator();
                const auto& last = gSim.appliedDrugs.back();
                const ChemicalEntity& lastDrug =
                    gBioagents.all()[last.entityIdx];
                ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.2f, 1.0f),
                                   "Active: %s @ %.2f µM",
                                   lastDrug.name.c_str(), last.dishConc_uM);
                ImGui::TextColored(ImVec4(0.6f, 0.75f, 0.9f, 1.0f),
                                   "Top affinities (score > 0.50):");
                for (int t = 0; t < (int)last.affinityPerTarget.size(); t++) {
                    const auto& aff = last.affinityPerTarget[t];
                    if (aff.score <= 0.50f) continue;
                    ImGui::Text("  %-22s  Kd=%.2f mM  s=%.2f",
                                gTargets.idAt(t).c_str(),
                                aff.Kd_mM, aff.score);
                }
            }
        }
        ImGui::End();
    }

    // ── Export Data panel (top-right, always visible) ────────────────────
    // One-click complete dump of every simulation artifact (time-series,
    // per-cell state, medium grid, apoptotic bodies, medium particles, and
    // a manifest with the mass-balance check) to a timestamped folder
    // inside the project's exports/ directory.
    {
        float winW = 280;
        float px = gSimMode == MODE_COLONY
            ? (ImGui::GetIO().DisplaySize.x - winW - 10)      // left of Population Stats? no, above
            : (ImGui::GetIO().DisplaySize.x - winW - 10);
        // Stack below Population Stats in colony mode, at top-right otherwise.
        float py = (gSimMode == MODE_COLONY) ? 300.0f : 10.0f;
        ImGui::SetNextWindowPos(ImVec2(px, py), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(winW, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Export Data");
        ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f),
                           "Dump every artifact to disk");
        ImGui::Separator();
        if (ImGui::Button("EXPORT EVERYTHING...", ImVec2(-1, 36))) {
            exportAllDataWithPicker();
        }
        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.6f, 1.0f),
                           "Bundle: population + cells + medium");
        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.6f, 1.0f),
                           "+ apoptotic bodies + chemicals + manifest");
        ImGui::Separator();
        ImGui::Text("Default folder (last used):");
        ImGui::TextWrapped("%s", gExportDir.c_str());
        if (!gLastExportPath.empty()) {
            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.4f, 0.85f, 0.45f, 1.0f),
                               "Last export:");
            // Show just the last folder name for readability.
            size_t slash = gLastExportPath.rfind('/');
            const char* leaf = slash == std::string::npos
                ? gLastExportPath.c_str()
                : gLastExportPath.c_str() + slash + 1;
            ImGui::TextWrapped("%s", leaf);
            if (ImGui::Button("Reveal in Finder", ImVec2(-1, 0))) {
                std::string cmd = "open \"" + gLastExportPath + "\"";
                (void)system(cmd.c_str());
            }
        }
        ImGui::End();
    }

    // Population panel (top-right) — colony mode only
    if (gSimMode == MODE_COLONY) {
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(270, 0), ImGuiCond_FirstUseEver);
    ImGui::Begin("Population Stats");
    {
        float bioH = gSim.bioTime / 3600;
        float bioM = fmodf(gSim.bioTime / 60, 60);
        ImGui::Text("BIO-TIME: %.0fh %02.0fm", floorf(bioH), bioM);
        ImGui::Separator();
        ImGui::Text("CELLS: %d", gSim.statAlive);
        ImGui::TextColored(ImVec4(0,1,0.6f,1), "  Proliferating: %d", gSim.statProlif);
        ImGui::TextColored(ImVec4(0.3f,0.5f,1,1), "  Quiescent: %d", gSim.statQuiescent);
        ImGui::TextColored(ImVec4(1,0.3f,0.1f,1), "  Apoptotic: %d", gSim.statApoptotic);
        ImGui::Separator();
        ImGui::Text("AVG ATP: %.1f", gSim.statAvgATP);
        ImGui::Text("Divisions: %d  Deaths: %d", gSim.statDivisions, gSim.statDeaths);
        ImGui::Separator();
        ImGui::TextColored(ImVec4(1,0.5f,0,1), "  Necrotic: %d", gSim.statNecrotic);
        ImGui::TextColored(ImVec4(1,0.7f,0.3f,1), "  Glycolytic: %d", gSim.statGlycolytic);
        ImGui::Separator();
        ImGui::Text("PHASE DISTRIBUTION:");
        float total = fmaxf(1, (float)gSim.statAlive);
        ImGui::TextColored(ImVec4(0.3f,0.5f,1,1),  "  G1: %d (%.0f%%)", gSim.statPhases[0], gSim.statPhases[0]/total*100);
        ImGui::TextColored(ImVec4(0.2f,1,0.6f,1),   "  S:  %d (%.0f%%)", gSim.statPhases[1], gSim.statPhases[1]/total*100);
        ImGui::TextColored(ImVec4(1,0.8f,0.2f,1),   "  G2: %d (%.0f%%)", gSim.statPhases[2], gSim.statPhases[2]/total*100);
        ImGui::TextColored(ImVec4(1,0.3f,0.6f,1),   "  M:  %d (%.0f%%)", gSim.statPhases[3], gSim.statPhases[3]/total*100);
    }
    ImGui::End();
    } // end colony-only population panel

    // Culture-medium readout (bottom-left). Closed-system dish: NO source
    // / sink controls — the user reads the running concentrations as the
    // cells consume / produce species. Initial recipe is DMEM HG + 10 % FBS.
    ImGui::SetNextWindowPos(ImVec2(10, ImGui::GetIO().DisplaySize.y - 360), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(310, 0), ImGuiCond_FirstUseEver);
    ImGui::Begin("Culture Medium (closed dish)");
    {
        ImGui::TextColored(ImVec4(0.7f,0.85f,1.0f,1.0f), "DMEM HG + 10%% FBS  (sealed)");
        ImGui::Separator();
        auto row = [](const char* label, float val, float init, const char* unit,
                      bool risingIsBad){
            float ratio = (init > 1e-6f) ? (val / init) : 1.0f;
            // Color: green near 1, amber when depleted/elevated, red when extreme
            float drift = fabsf(ratio - 1.0f);
            ImVec4 col(0.45f, 0.95f, 0.55f, 1.0f);
            if (drift > 0.20f) col = ImVec4(0.95f, 0.85f, 0.30f, 1.0f);
            if (drift > 0.50f) col = ImVec4(0.95f, 0.45f, 0.30f, 1.0f);
            (void)risingIsBad;
            ImGui::Text("%-12s", label);
            ImGui::SameLine(110);
            ImGui::TextColored(col, "%6.2f → %6.2f %s", init, val, unit);
        };
        row("Glucose",   (float)gSim.nutrients.mean(MS_GLUCOSE),   MediumComposition::DMEM_GLUCOSE_MM,   "mM", false);
        row("Glutamine", (float)gSim.nutrients.mean(MS_GLUTAMINE), MediumComposition::DMEM_GLUTAMINE_MM, "mM", false);
        row("Pyruvate",  (float)gSim.nutrients.mean(MS_PYRUVATE),  MediumComposition::DMEM_PYRUVATE_MM,  "mM", false);
        row("AA pool",   (float)gSim.nutrients.mean(MS_AA_POOL),   MediumComposition::DMEM_AA_POOL_MM,   "mM", false);
        row("O2",        (float)gSim.nutrients.mean(MS_O2),        MediumComposition::DMEM_O2_MM,        "mM", false);
        row("CO2",       (float)gSim.nutrients.mean(MS_CO2),       MediumComposition::DMEM_CO2_MM,       "mM", true);
        row("Lactate",   (float)gSim.nutrients.mean(MS_LACTATE),   1.0f,                                  "mM", true);
        row("pH",        (float)gSim.nutrients.mean(MS_HPLUS),     MediumComposition::DMEM_PH,           "",   false);
        row("Growth f.", (float)gSim.nutrients.mean(MS_GROWTH_F),  MediumComposition::DMEM_GROWTH_FACTOR_NG_PER_ML, "ng/mL", false);
        // Mass-balance bar — green if every conserved species drift < 1 %.
        double drifts[MS_COUNT];
        double maxDrift = gSim.nutrients.checkBalance(drifts);
        ImVec4 mbCol = (maxDrift < 0.001) ? ImVec4(0.45f,0.95f,0.55f,1)
                     : (maxDrift < 0.01)  ? ImVec4(0.95f,0.85f,0.30f,1)
                                          : ImVec4(0.95f,0.30f,0.30f,1);
        ImGui::Separator();
        ImGui::TextColored(mbCol, "Mass balance drift: %.4f%%", maxDrift * 100.0);
        ImGui::Separator();
        // Logarithmic speed slider (0.1× to 100×)
        float logSpeed = log10f(gSim.timeScale);
        if (ImGui::SliderFloat("##speed", &logSpeed, -1.0f, 2.0f, "")) {
            gSim.timeScale = powf(10.0f, logSpeed);
        }
        float bioMinPerSec = gSim.timeScale * BIO_MIN_PER_SEC;
        if (bioMinPerSec < 1.0f)
            ImGui::Text("%.1fx | %.0fs bio/s", gSim.timeScale, bioMinPerSec * 60);
        else if (bioMinPerSec >= 60)
            ImGui::Text("%.1fx | %.1fh bio/s", gSim.timeScale, bioMinPerSec / 60);
        else
            ImGui::Text("%.1fx | %.1fm bio/s", gSim.timeScale, bioMinPerSec);

        if (gSim.pendingScaledDt > 0.01f) {
            ImGui::TextColored(ImVec4(0.95f, 0.72f, 0.30f, 1.0f),
                               "CPU backlog: %.2fs queued (no sim dropped)", gSim.pendingScaledDt);
        }

        // Speed presets
        if (ImGui::SmallButton("0.5x")) gSim.timeScale = 0.5f; ImGui::SameLine();
        if (ImGui::SmallButton("1x")) gSim.timeScale = 1.0f; ImGui::SameLine();
        if (ImGui::SmallButton("5x")) gSim.timeScale = 5.0f; ImGui::SameLine();
        if (ImGui::SmallButton("20x")) gSim.timeScale = 20.0f; ImGui::SameLine();
        if (ImGui::SmallButton("50x")) gSim.timeScale = 50.0f;

        if (ImGui::Button(gSim.paused ? ">> RESUME" : "|| PAUSE", ImVec2(-1, 0)))
            gSim.paused = !gSim.paused;

        ImGui::Separator();
        ImGui::TextColored(ImVec4(0.4f,0.6f,0.8f,0.6f), "FPS: %.0f | Cells: %d",
                           ImGui::GetIO().Framerate, gSim.statAlive);

        // Camera sensitivity (collapsible, inside same window)
        if (ImGui::CollapsingHeader("Camera Sensitivity")) {
            ImGui::SliderFloat("Rotate", &gCamera.rotateSensitivity, 0.1f, 3.0f, "%.1f");
            ImGui::SliderFloat("Pan",    &gCamera.panSensitivity,    0.1f, 3.0f, "%.1f");
            ImGui::SliderFloat("Zoom",   &gCamera.zoomSensitivity,   0.1f, 3.0f, "%.1f");
            ImGui::SliderFloat("Move",   &gCamera.moveSensitivity,   0.1f, 3.0f, "%.1f");
            if (ImGui::SmallButton("Reset")) {
                gCamera.rotateSensitivity = gCamera.panSensitivity =
                gCamera.zoomSensitivity = gCamera.moveSensitivity = 1.0f;
            }
        }
        ImGui::TextColored(ImVec4(0.4f,0.6f,0.8f,0.4f), "WASD:move Space:up Shift:down P:pause");
        if (ImGui::Button("Reset View", ImVec2(-1, 24))) {
            if (gSimMode == MODE_SINGLE_CELL)
                gCamera.setSingleCellView();
            else
                gCamera.setColonyView();
        }
    }
    ImGui::End();

    // Time-series plots — colony mode only
    if (gSimMode == MODE_COLONY) {
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(420, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Population Dynamics");
    {
        int n = gTS.count();
        ImGui::Text("Data points: %d", n);

        if (n > 1) {
            const float* t = gTS.time.data();

            if (ImPlot::BeginPlot("Cell Count", ImVec2(-1, 130))) {
                ImPlot::SetupAxes("Bio-time (h)", "Cells");
                ImPlot::PlotLine("Total", t, gTS.population.data(), n);
                ImPlot::PlotLine("Prolif", t, gTS.proliferating.data(), n);
                ImPlot::PlotLine("Quiesc", t, gTS.quiescent.data(), n);
                ImPlot::EndPlot();
            }

            if (ImPlot::BeginPlot("Metabolism", ImVec2(-1, 130))) {
                ImPlot::SetupAxes("Bio-time (h)", "");
                ImPlot::PlotLine("Avg ATP", t, gTS.avgATP.data(), n);
                ImPlot::PlotLine("Avg Stress", t, gTS.avgStress.data(), n);
                ImPlot::PlotLine("Glycolytic %", t, gTS.glycolyticPct.data(), n);
                ImPlot::EndPlot();
            }

            if (ImPlot::BeginPlot("Phase Distribution %", ImVec2(-1, 130))) {
                ImPlot::SetupAxes("Bio-time (h)", "%");
                ImPlot::PlotLine("G1", t, gTS.phaseG1.data(), n);
                ImPlot::PlotLine("S", t, gTS.phaseS.data(), n);
                ImPlot::PlotLine("G2", t, gTS.phaseG2.data(), n);
                ImPlot::PlotLine("M", t, gTS.phaseM.data(), n);
                ImPlot::EndPlot();
            }
        } else {
            ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1), "Waiting for data...");
        }

        ImGui::Separator();
        if (ImGui::Button("EXPORT CSV (Save As...)", ImVec2(-1, 28))) {
            // Default timestamped filename
            time_t now = ::time(nullptr);
            struct tm* t = localtime(&now);
            char fname_buf[128];
            snprintf(fname_buf, sizeof(fname_buf), "cellsim_%04d%02d%02d_%02d%02d%02d.csv",
                     t->tm_year+1900, t->tm_mon+1, t->tm_mday,
                     t->tm_hour, t->tm_min, t->tm_sec);
            NSString* defaultName = [NSString stringWithUTF8String:fname_buf];

            // Pop up native macOS Save dialog
            dispatch_async(dispatch_get_main_queue(), ^{
                NSSavePanel* panel = [NSSavePanel savePanel];
                [panel setTitle:@"Export CellSim Data"];
                [panel setNameFieldStringValue:defaultName];
                [panel setAllowedContentTypes:@[[UTType typeWithFilenameExtension:@"csv"]]];
                // Default to project exports folder
                [panel setDirectoryURL:[NSURL fileURLWithPath:
                    [NSString stringWithUTF8String:gExportDir.c_str()]]];
                [panel setCanCreateDirectories:YES];

                if ([panel runModal] == NSModalResponseOK) {
                    NSURL* url = [panel URL];
                    std::string savePath = [[url path] UTF8String];
                    gTS.exportCSV(savePath);
                    // Remember the chosen directory for next time
                    gExportDir = [[[url path] stringByDeletingLastPathComponent] UTF8String];
                }
            });
        }
        // Per-cell snapshot export
        if (ImGui::Button("EXPORT CELL SNAPSHOT", ImVec2(-1, 24))) {
            time_t now2 = ::time(nullptr);
            struct tm* t2 = localtime(&now2);
            char fname2[128];
            snprintf(fname2, sizeof(fname2), "cellsim_cells_%04d%02d%02d_%02d%02d%02d.csv",
                     t2->tm_year+1900, t2->tm_mon+1, t2->tm_mday,
                     t2->tm_hour, t2->tm_min, t2->tm_sec);
            NSString* defaultName2 = [NSString stringWithUTF8String:fname2];
            dispatch_async(dispatch_get_main_queue(), ^{
                NSSavePanel* panel = [NSSavePanel savePanel];
                [panel setTitle:@"Export Cell Snapshot"];
                [panel setNameFieldStringValue:defaultName2];
                [panel setAllowedContentTypes:@[[UTType typeWithFilenameExtension:@"csv"]]];
                [panel setDirectoryURL:[NSURL fileURLWithPath:
                    [NSString stringWithUTF8String:gExportDir.c_str()]]];
                [panel setCanCreateDirectories:YES];
                if ([panel runModal] == NSModalResponseOK) {
                    std::string savePath = [[[panel URL] path] UTF8String];
                    gTS.exportCellSnapshot(savePath, gSim);
                    gExportDir = [[[[panel URL] path] stringByDeletingLastPathComponent] UTF8String];
                }
            });
        }
        ImGui::TextColored(ImVec4(0.4f,0.6f,0.8f,0.5f), "Dir: %s",
                           gExportDir.substr(gExportDir.rfind('/')+1).c_str());
    }
    ImGui::End();
    } // end colony-only time-series plots

    // Selected cell detail — always show in single cell mode, conditional in colony
    activeFocusCellIndex();
    if (gSelectedCell >= 0 && gSelectedCell < (int)gSim.cells.size()) {
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 350), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(270, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Cell Research Panel");
        {
            auto& c = gSim.cells[gSelectedCell];
            int prevSelected = gSelectedCell;
            ImGui::SliderInt("Cell #", &gSelectedCell, 0, (int)gSim.cells.size()-1);
            if (gSelectedCell != prevSelected) {
                selectCellIndex(gSelectedCell);
            }
            ImGui::Separator();

            const char* phaseNames[] = {"G1","S","G2","M"};
            const char* fateNames[] = {"Undetermined","Proliferating","Quiescent","Apoptotic"};
            ImVec4 phaseColors[] = {{0.3f,0.5f,1,1},{0.2f,1,0.6f,1},{1,0.8f,0.2f,1},{1,0.3f,0.6f,1}};

            ImGui::TextColored(phaseColors[c.phase], "Phase: %s", phaseNames[c.phase]);
            ImGui::Text("Fate: %s", fateNames[c.fate]);
            ImGui::ProgressBar(c.cycleProgress, ImVec2(-1,12), "Cycle");

            ImGui::Separator();
            ImGui::Text("CDK/CYCLIN (Novak-Tyson):");
            ImGui::ProgressBar(c.cdk.CycD/1.5f, ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycD %.3f", c.cdk.CycD);
            ImGui::ProgressBar(c.cdk.Rb,         ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Rb   %.3f", c.cdk.Rb);
            ImGui::ProgressBar(c.cdk.E2F,        ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("E2F  %.3f", c.cdk.E2F);
            ImGui::ProgressBar(c.cdk.CycE/1.5f,  ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycE %.3f", c.cdk.CycE);
            ImGui::ProgressBar(c.cdk.CycA/1.5f,  ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycA %.3f", c.cdk.CycA);
            ImGui::ProgressBar(c.cdk.CycB/1.5f,  ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycB %.3f", c.cdk.CycB);
            ImGui::ProgressBar(c.cdk.p21,        ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("p21  %.3f", c.cdk.p21);

            ImGui::Separator();
            ImGui::Text("METABOLISM:");
            ImGui::ProgressBar(c.ATP/100,      ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("ATP    %.1f", c.ATP);
            ImGui::ProgressBar(c.stress/100,   ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Stress %.1f", c.stress);
            ImGui::ProgressBar(c.ROS/100,      ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("ROS    %.1f", c.ROS);
            ImGui::ProgressBar(c.biomass/2.3f, ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Biomass %.2f", c.biomass);
            ImGui::ProgressBar(c.damageLevel/2,ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Damage %.3f", c.damageLevel);

            ImGui::Separator();
            ImGui::Text("MITOCHONDRIA:");
            ImGui::ProgressBar((c.mitoPotential-40)/180, ImVec2(-1,8), ""); ImGui::SameLine();
            ImGui::Text("Psi %.0f mV", c.mitoPotential);
            ImGui::ProgressBar(c.mitoHealth, ImVec2(-1,8), ""); ImGui::SameLine();
            ImGui::Text("Health %.2f", c.mitoHealth);
            ImGui::TextColored(c.glycolytic ? ImVec4(1,0.5f,0.2f,1) : ImVec4(0.5f,0.9f,0.6f,1),
                "Mode: %s", c.glycolytic ? "GLYCOLYTIC (Warburg)" : "OXIDATIVE");
            if (c.necrotic) ImGui::TextColored(ImVec4(1,0.3f,0,1), "!! NECROTIC (O2 critical)");

            ImGui::Separator();
            ImGui::Text("GENOME:");
            ImGui::Text("  Gen:%d Clone:%d", c.generation, c.cloneId);
            ImGui::Text("  Telomere: %.0f bp %s", c.telomere, c.senescent ? "[SENESCENT]" : "");
            ImGui::Text("  gly:%.2f pro:%.2f ros:%.2f rep:%.2f",
                         c.glycolysisBias, c.prolifBias, c.rosTolerance, c.repairRate);
            ImGui::Text("  Pressure: %.2f  Hypoxia: %.0fs", c.localPressure, c.hypoxiaTimer);
            ImGui::Text("  ATP danger: %.0f/%.0fs", c.atpDangerTimer, ATP_DANGER_DURATION);

            // Drug response readout removed 2026-04-19 — pending rewrite.
        }
        ImGui::End();
    }

    // ══════════════════════════════════════════════════════════════════
    //  MOLECULE VIEWER (single cell mode)
    // ══════════════════════════════════════════════════════════════════
    if (gSimMode == MODE_SINGLE_CELL) {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(250, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Molecules Inside Cell");
        {
            int molCount = gMolCache.loadedCount();

            if (molCount > 0) {
                ImGui::Text("Particles: %d", (int)gIntraPhys.particles.size());
                ImGui::Text("Rendered: %d atoms, %d bonds", gMolRender.atomCount, gMolRender.bondCount);

                ImGui::Separator();
                ImGui::TextColored(ImVec4(0.3f,0.7f,1.0f,1), "DNA: HBB gene (%d bp, real sequence)", HBB_LENGTH);
                int visibleDNA = 0, visibleSisterDNA = 0;
                countActiveDNAVisuals(gIntraPhys, visibleDNA, visibleSisterDNA);
                ImGui::TextColored(ImVec4(0.95f,0.82f,0.35f,1),
                                   "Visible DNA nodes: %d  sister copies: %d",
                                   visibleDNA, visibleSisterDNA);

                // Central dogma status — show REAL per-cell totals from
                // SimCell biology (literature counts, Milo 2013 Bioessays)
                // alongside the currently-engaged teaching-animation
                // subset. "Active" = right now elongating a transcript /
                // splicing / translating. Total = literature count per
                // cell (most are idle at any moment).
                int activeTranscr = gCDogma.countActiveTranscription();
                int activeSplice = gCDogma.countActiveSplicing();
                int activeTransl = gCDogma.countActiveTranslation();
                int activeProteins = gCDogma.countActiveProteins();
                int activeTRNAs = gCDogma.countActiveTRNAs();

                auto fmtSciStr = [](double v) -> std::string {
                    char buf[32];
                    if (v >= 1e12)      snprintf(buf, sizeof(buf), "%.2fT", v / 1e12);
                    else if (v >= 1e9)  snprintf(buf, sizeof(buf), "%.2fG", v / 1e9);
                    else if (v >= 1e6)  snprintf(buf, sizeof(buf), "%.2fM", v / 1e6);
                    else if (v >= 1e3)  snprintf(buf, sizeof(buf), "%.2fk", v / 1e3);
                    else                snprintf(buf, sizeof(buf), "%.0f",  v);
                    return std::string(buf);
                };

                const SimCell* pc = gSim.cells.empty() ? nullptr : &gSim.cells[0];
                double tRibo = pc ? pc->ribosomeCount  : 0.0;
                double tPolII= pc ? pc->rnaPolII       : 0.0;
                double tSpl  = pc ? pc->spliceosomes   : 0.0;
                double tmRNA = pc ? pc->mRNACount      : 0.0;
                double ttRNA = pc ? pc->tRNACount      : 0.0;
                double tProt = pc ? pc->totalProteins  : 0.0;
                double tPore = pc ? pc->nuclearPores   : 0.0;
                double tChap = pc ? pc->chaperones     : 0.0;
                double tProte= pc ? pc->proteasomes    : 0.0;

                ImGui::TextColored(ImVec4(0.6f,0.2f,0.8f,1), "RNA Pol II: %s total (%d engaged)", fmtSciStr(tPolII).c_str(), activeTranscr);
                ImGui::TextColored(ImVec4(0.8f,0.6f,0.2f,1), "Spliceosomes: %s total (%d active)", fmtSciStr(tSpl).c_str(), activeSplice);
                ImGui::TextColored(ImVec4(0.5f,0.4f,0.8f,1), "Ribosomes: %s total (%d translating)", fmtSciStr(tRibo).c_str(), activeTransl);
                ImGui::TextColored(ImVec4(0.9f,0.4f,0.1f,1), "Mature mRNAs: %s (teaching subset: %d)", fmtSciStr(tmRNA).c_str(), gCDogma.activeMRNAs);
                ImGui::TextColored(ImVec4(0.3f,0.85f,0.6f,1), "tRNAs: %s total (%d charged now)", fmtSciStr(ttRNA).c_str(), activeTRNAs);
                ImGui::TextColored(ImVec4(0.9f,0.7f,0.35f,1), "Total protein: %s (teaching: %d)", fmtSciStr(tProt).c_str(), activeProteins);
                ImGui::TextColored(ImVec4(0.4f,0.7f,0.4f,1), "Nuclear pores: %s", fmtSciStr(tPore).c_str());
                ImGui::TextColored(ImVec4(0.7f,0.5f,0.3f,1), "Chaperones: %s", fmtSciStr(tChap).c_str());
                ImGui::TextColored(ImVec4(0.7f,0.7f,0.4f,1), "Proteasomes: %s", fmtSciStr(tProte).c_str());

                bool replicationTelemetryVisible =
                    gCDogma.sPhaseProgramStarted ||
                    gCDogma.replicationProgress > 0.0f ||
                    gCDogma.unresolvedReplicationErrors > 0 ||
                    gCDogma.countActiveReplicationForks() > 0;
                if (replicationTelemetryVisible) {
                    ImGui::TextColored(ImVec4(1.0f,0.3f,0.3f,1),
                                       "DNA REPLICATION: %.0f%%  forks:%d  origins:%d/%d",
                                       gCDogma.replicationProgress * 100,
                                       gCDogma.countActiveReplicationForks(),
                                       gCDogma.countFiredOrigins(),
                                       CDOGMA_MAX_REPL_ORIGINS);
                    ImGui::TextColored(ImVec4(0.9f,0.65f,0.25f,1),
                                       "dNTP supply: %.0f%%  CHK1: %.0f%%  quality: %.0f%%",
                                       gCDogma.dntpAvailability * 100.0f,
                                       gCDogma.chk1Signal * 100.0f,
                                       gCDogma.replicationQuality * 100.0f);
                    ImGui::TextColored(ImVec4(0.65f,0.85f,0.95f,1),
                                       "Fork calibration: %.0f bp/s (~%.2f kb/min literature target)",
                                       CDOGMA_LIT_FORK_SPEED_BP_PER_SEC,
                                       CDOGMA_LIT_FORK_SPEED_BP_PER_SEC * 60.0f / 1000.0f);
                    ImGui::TextColored(ImVec4(0.55f,0.8f,1.0f,1),
                                       "Proofread: %d  MMR fixed: %d  unresolved: %d  escaped: %d",
                                       gCDogma.proofreadCorrections,
                                       gCDogma.mmrCorrections,
                                       gCDogma.unresolvedReplicationErrors,
                                       gCDogma.escapedErrors);
                    ImGui::TextColored(ImVec4(0.8f,0.7f,1.0f,1),
                                       "Dormant origins fired: %d  predicted fork risk: %.2f%%",
                                       gCDogma.dormantOriginsFired,
                                       gCDogma.predictedErrorRisk * 100.0f);
                }

                ImGui::Separator();
                ImGui::TextColored(ImVec4(0.4f,0.6f,0.4f,1), "Physics: Langevin dynamics");

                // List loaded molecules
                auto loadedIds = gMolCache.getLoadedIds();
                if (ImGui::TreeNode("Molecule List")) {
                    for (auto& id : loadedIds) {
                        const MoleculeData* mol = gMolCache.get(id);
                        std::string displayName = id;
                        for (int i = 0; i < MOLECULE_ENTRY_COUNT; i++) {
                            if (id == MOLECULE_ENTRIES[i].id) {
                                displayName = MOLECULE_ENTRIES[i].name; break;
                            }
                        }
                        if (mol) {
                            ImGui::BulletText("%s (%d atoms, %.1f A dia)",
                                              displayName.c_str(),
                                              (int)mol->atoms.size(),
                                              mol->maxExtent * 2.0f);
                        }
                    }
                    ImGui::TreePop();
                }
            } else {
                ImGui::TextColored(ImVec4(1,0.5f,0.2f,1), "No molecules loaded.");
                ImGui::TextWrapped("Run: python3 scripts/generate_library.py");
                if (gMolGen.isReady()) {
                    if (ImGui::Button("Generate Library", ImVec2(-1, 28))) {
                        gMolGen.generateLibrary(gMoleculeDir);
                        gMolCache.init(gMoleculeDir);
                    }
                }
            }
        }
        ImGui::End();
    }

    if (gSimMode == MODE_SINGLE_CELL) {
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, ImGui::GetIO().DisplaySize.y - 210),
                                ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(270, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Mitosis Debug");
        {
            ImGui::Checkbox("Show Overlay", &gShowMitosisDebugOverlay);
            ImGui::Checkbox("Cell Labels", &gDebugOverlayShowLabels);
            ImGui::Checkbox("Nucleus Markers", &gDebugOverlayShowNucleusMarkers);
            ImGui::Checkbox("Collision Rings", &gDebugOverlayShowCollisionRings);
            ImGui::Checkbox("Warnings", &gDebugOverlayShowWarnings);
            ImGui::Separator();
            ImGui::Text("Cells: %d", (int)gSim.cells.size());
            ImGui::Text("Mitosis: %s", gMitosis.active ? "ACTIVE" : "idle");
            if (gMitosis.active) {
                ImGui::Text("Phase: %d", (int)gMitosis.phase);
            }
            ImGui::TextWrapped("Overlay shows cell ID, clone ID, DNA count, nucleus marker, and interior-missing warnings directly over each cell.");
        }
        ImGui::End();
    }

    // Drug Treatment panel removed 2026-04-19 — pending rewrite.
}

// ── Render frame ────────────────────────────────────────────────────────
static void renderFrame(float time, float dt) {
    @autoreleasepool {
        if (gSimMode == MODE_SINGLE_CELL) {
            ensurePrimaryCellBinding();
        }
        persistBoundDogmaState();
        persistBoundMitosisState();

        // Sync visual mitosis flag before simulation step
        gSim.visualMitosisActive = gMitosis.active;
        gSim.primaryReplicationProgress = gCDogma.replicationProgress;
        gSim.primaryReplicationCheckpoint = gCDogma.chk1Signal;
        gSim.primaryReplicationQuality = gCDogma.replicationQuality;
        gSim.primaryReplicationReady = gCDogma.replicationReadyForM();
        stepBackgroundDogmaStates(dt);
        // Skip all simulation + chemical-animation updates while the
        // Setup overlay is open — the user hasn't picked conditions yet,
        // so nothing should be moving / dividing in the background.
        if (!gSetup.showOverlay) {
            gSim.update(dt);
            updateMediumChemicals(dt);
            // Apoptotic-body system: spawn bodies from cells that
            // crossed into FRAGMENTATION/BODIES this tick, then update
            // the drift/decomposition/lysis state of existing bodies,
            // then let nearby live cells engulf them (efferocytosis)
            // and digest lysosomal contents.
            for (SimCell& c : gSim.cells) {
                if (c.apoPhase == Apoptosis::FRAGMENTATION
                    && !c.bodiesSpawned
                    && c.alive) {
                    spawnApoptoticBodies(c);
                    c.bodiesSpawned = true;
                }
            }
            updateApoptoticBodies(dt);
            updateEfferocytosis(dt);
            updateLysosomalDigestion(dt);
            applyDampEventsToNeighbors(dt);
        }
        // Comprehensive science-grade telemetry — opens on first frame,
        // samples once per bio-minute, captures the full population +
        // per-element environment for cross-validation against real
        // HeLa experiments.
        if (!gTelemetryOpened) {
            char sid[64];
            time_t nowT = ::time(nullptr);
            std::strftime(sid, sizeof(sid),
                          "session_%Y%m%d_%H%M%S", std::localtime(&nowT));
            gTelemetryOpened = gTelemetryLog.open(sid);
            if (gTelemetryOpened) {
                printf("[Telemetry] logging to logs/telemetry_population_%s.csv\n", sid);
            }
        }
        if (gTelemetryOpened) {
            gTelemetryLog.tick(gSim, gSim.bioTime, (float)glfwGetTime());
        }

        bool sawBackgroundDivision = false;
        Simulation::DivisionEvent latestBackgroundDivision{};

        // Process division events from background cells (interior splitting)
        for (auto& evt : gSim.pendingDivisions) {
            auto parentIt = gCellInteriors.find(evt.parentCellUid);
            if (parentIt != gCellInteriors.end() && parentIt->second.initialized) {
                // Split the parent's interior between the two daughters
                auto& parentPhys = parentIt->second.phys;
                float dR = parentPhys.params.cellRadius * 0.794f;
                prepareInteriorForPhysicalDivision(parentPhys, parentPhys.params.cellRadius);
                MitosisHalfCenters halfCenters =
                    computeHalfCentersForInterior(parentPhys, parentPhys.params.cellRadius * 0.20f, false);
                CellInterior intA = splitMitosisDaughterInterior(parentPhys, 0, halfCenters.local[0], dR);
                CellInterior intB = splitMitosisDaughterInterior(parentPhys, 1, halfCenters.local[1], dR);
                gCellInteriors[evt.daughterACellUid] = std::move(intA);
                gCellInteriors[evt.daughterBCellUid] = std::move(intB);
                splitOrganelleMotionState(evt.parentCellUid,
                                          evt.daughterACellUid,
                                          evt.daughterBCellUid);
                if (evt.parentCellUid != evt.daughterACellUid &&
                    evt.parentCellUid != evt.daughterBCellUid) {
                    gCellInteriors.erase(evt.parentCellUid);
                }
            }
            // Freshly initialize each daughter's own replication program.
            // The daughters' CellCycleProgram lives on the SimCell itself,
            // so we locate them by UID rather than a side map.
            if (SimCell* dA = findCellByUid(evt.daughterACellUid)) {
                dA->program.cdogmaInitialized = true;
            }
            if (SimCell* dB = findCellByUid(evt.daughterBCellUid)) {
                dB->program.cdogmaInitialized = true;
            }
            if (gSimMode == MODE_SINGLE_CELL) {
                latestBackgroundDivision = evt;
                sawBackgroundDivision = true;
            }
            // If no interior existed, initCellInterior will create one on demand
        }
        gSim.pendingDivisions.clear();

        // Cleanup stale cellUid entries from dead cells
        {
            std::set<int> liveUids;
            for (auto& c : gSim.cells) if (c.alive) liveUids.insert(c.cellUid);
            for (auto it = gCellInteriors.begin(); it != gCellInteriors.end(); ) {
                if (liveUids.find(it->first) == liveUids.end()) it = gCellInteriors.erase(it);
                else ++it;
            }
            for (auto it = gOrganelleMotionStates.begin(); it != gOrganelleMotionStates.end(); ) {
                if (liveUids.find(it->first) == liveUids.end()) it = gOrganelleMotionStates.erase(it);
                else ++it;
            }
            // cell.program lifetimes now tied to SimCell lifetime — no
            // separate map to clean up.
        }

        // Keep the detailed single-cell pipeline bound to one deliberate
        // lineage cell instead of whatever happened to slide into cells[0]
        // after dead-cell erasure or background divisions.
        if (gSimMode == MODE_SINGLE_CELL) {
            ensurePrimaryCellBinding();
        }

        if (gSimMode == MODE_SINGLE_CELL && sawBackgroundDivision) {
            int idxA = findCellIndexByUid(latestBackgroundDivision.daughterACellUid);
            int idxB = findCellIndexByUid(latestBackgroundDivision.daughterBCellUid);
            if (idxA >= 0 && idxB >= 0) {
                const auto& a = gSim.cells[idxA];
                const auto& b = gSim.cells[idxB];
                gPostMitosisPairUidA = a.cellUid;
                gPostMitosisPairUidB = b.cellUid;
                gPostMitosisPairA = a.position;
                gPostMitosisPairB = b.position;
                gPostMitosisPairRadius = fmaxf(a.radius * a.size, b.radius * b.size);
                gPostMitosisPairCameraTimer = fmaxf(gPostMitosisPairCameraTimer, 2.75f);
            }
        }

        // Step the full cell biology engine (all 16 modules)
        if (gSimMode == MODE_SINGLE_CELL && !gSim.cells.empty()) {
            if (!gCellBio.initialized) {
                gCellBio.init();
                syncPrimarySimATPIntoBioEngine();
            }
            // Couple external O2 (from nutrient field) into the metabolism pool.
            // NutrientField O2 is 0–1 (fraction of atmospheric). 100% → ~0.03 mM
            // dissolved O2 at tissue level (Alberts Ch 14).
            float envO2_frac = gSim.nutrients.getO2(gSim.cells[0].position.x,
                                                    gSim.cells[0].position.z);
            gCellBio.metabolism.setExternalO2(envO2_frac * 0.03f);
            gCellBio.step(gSim.lastExecutedScaledDt * BIO_MIN_PER_SEC);
            // Feed ATP from new engine back to old simulation
            gSim.cells[0].ATP = bioATPToSimPercent(gCellBio.metabolism.pool.ATP);
        }

        // Update the primary interior before rebuilding cell instances.
        // Mitosis finalization happens inside uploadCellInterior(), so syncing
        // shells first can leave the render one frame behind and cause a
        // visible disappear/reappear pop at cytokinesis handoff.
        if (gSimMode == MODE_SINGLE_CELL && !gSim.cells.empty()) {
            auto& cell = gSim.cells[0];
            uploadCellInterior(cell.position, cell.radius * cell.size, time, dt);
            persistBoundDogmaState();
            persistBoundMitosisState();
        }

        syncCellInstances();

        // Time-series sampling (every 0.5 real seconds)
        gTS.sampleTimer += dt;
        if (gTS.sampleTimer > 0.5f) {
            gTS.sampleTimer = 0;
            gTS.sample(gSim);
        }

        id<CAMetalDrawable> drawable = [gCtx.metalLayer() nextDrawable];
        if (!drawable) return;

        // DO NOT early-return when cellCount == 0. When the Setup overlay
        // is showing (no sim inited yet), cellCount is 0 but we still
        // need to clear the screen to black and draw ImGui on top.
        // Returning here leaves the drawable unpresented → white screen.
        int cellCount = (int)gCellInstances.size();

        Uniforms uni;
        uni.viewProjection = gCamera.getViewProjection();
        uni.model = matrix_identity_float4x4;
        uni.cameraPos = gCamera.getPosition();
        uni.time = time;
        uni.lightDir = simd_normalize(simd_make_float3(0.5f, 0.9f, 0.6f));
        uni.pad0 = 0;

        id<MTLCommandBuffer> cmd = [gCtx.commandQueue() commandBuffer];
        MTLRenderPassDescriptor* pass = [MTLRenderPassDescriptor renderPassDescriptor];
        pass.colorAttachments[0].texture = drawable.texture;
        pass.colorAttachments[0].loadAction = MTLLoadActionClear;
        pass.colorAttachments[0].storeAction = MTLStoreActionStore;
        pass.colorAttachments[0].clearColor = MTLClearColorMake(0.004, 0.008, 0.031, 1.0);
        pass.depthAttachment.texture = gCtx.depthTexture();
        pass.depthAttachment.loadAction = MTLLoadActionClear;
        pass.depthAttachment.storeAction = MTLStoreActionDontCare;
        pass.depthAttachment.clearDepth = 1.0;

        id<MTLRenderCommandEncoder> enc = [cmd renderCommandEncoderWithDescriptor:pass];

        // While the Setup overlay is showing, don't draw ANY scene
        // geometry — substrate, fluid, organelles, cells, molecules all
        // suppressed. The user sees only the clear-colour background and
        // the ImGui setup panel on top. Simulation isn't running yet.
        if (!gSetup.showOverlay) {

        // 1. Substrate
        {
            [enc setDepthStencilState:gCtx.depthState()];
            const MeshData& sub = gMeshes.substrate();
            Uniforms su = uni;
            simd_float3 fp = {0, FLOOR_Y, 0};
            su.model = mat4_translation(fp);
            [enc setRenderPipelineState:gCtx.substratePipeline()];
            [enc setVertexBuffer:sub.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBytes:&su length:sizeof(Uniforms) atIndex:2];
            [enc setFragmentBytes:&su length:sizeof(Uniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle indexCount:sub.indexCount
                             indexType:MTLIndexTypeUInt32 indexBuffer:sub.indexBuffer indexBufferOffset:0];
        }

        // 1b. Fluid backdrop — drawn FIRST (before organelles/cells/molecules)
        // so cells render fully on top of it instead of being washed out
        // by it. The dish is full of dark medium; everything inside the
        // medium glows ON TOP at its full unfaded color.
        if (gFluidVertexBuffer && gFluidIndexBuffer && gCtx.fluidPipeline()) {
            FluidUniforms fu;
            fu.viewProjection = uni.viewProjection;
            fu.cameraPos      = uni.cameraPos;
            fu.time           = time;
            fu.lightDir       = uni.lightDir;
            fu.floorY         = FLOOR_Y;
            // Very dark navy — deep medium so glowing organelles + molecules pop.
            fu.waterColor     = {0.02f, 0.05f, 0.18f};
            fu.waterAlpha     = 0.65f;     // dense visible body
            fu.radius         = SCENE_BOUND * 0.99f;
            fu.fluidHeight    = MEDIUM_FLUID_HEIGHT;
            fu.waveAmp        = 0.55f;
            fu.waveSpeed      = 0.65f;
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            [enc setRenderPipelineState:gCtx.fluidPipeline()];
            [enc setVertexBuffer:gFluidVertexBuffer offset:0 atIndex:0];
            [enc setVertexBytes:&fu length:sizeof(FluidUniforms) atIndex:1];
            [enc setFragmentBytes:&fu length:sizeof(FluidUniforms) atIndex:1];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                            indexCount:gFluidIndexCount
                             indexType:MTLIndexTypeUInt32
                           indexBuffer:gFluidIndexBuffer
                     indexBufferOffset:0];
            // Restore depth state for subsequent passes
            [enc setDepthStencilState:gCtx.depthState()];
        }

        // 2. GLB Organelles (only for nearby cells — LOD optimization)
        {
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            [enc setRenderPipelineState:gCtx.glbOrganellePipeline()];

            simd_float3 camPos = gCamera.getPosition();
            for (int i = 0; i < cellCount; i++) {
                simd_float3 cp = gCellInstances[i].position;
                float cs = gCellInstances[i].radius;
                float toff = time + (float)i * 7.3f;
                int sourceIdx = renderSourceIndex(i);
                int orgKey = organelleMotionKeyForRenderIndex(i);
                OrganelleConfigs orgCfg = animatedOrganelleConfigsForCell(orgKey, time);

                // LOD: only render organelles for cells within 60 units of camera
                float dx = cp.x-camPos.x, dy = cp.y-camPos.y, dz = cp.z-camPos.z;
                if (dx*dx+dy*dy+dz*dz > 3600) continue;

                // NOTE: the previous build scripted an "organelle duplication
                // and migration" animation during mitosis that reset the
                // render to a hardcoded 3-mito, single-organelle config and
                // pushed two copies toward the poles. That looked fake: the
                // cell never actually re-spawns organelles as a choreographed
                // event. Real cells keep their existing organelles through
                // M-phase; physical partitioning at cytokinesis is what
                // determines which daughter gets which. We simply continue
                // the normal full-count render below (with all 500 mitos,
                // 80 ribosomes, etc.) during every phase of mitosis — no
                // transform, no re-spawn, just simulation continuing. The
                // spindle / chromosome / furrow rendering overlays the
                // mitotic machinery on top.

                // Always render organelles — they scale with cellSize automatically.
                // During MITO_COMPLETE, both daughter cells render organelles (i=0 and i=1).
                drawGLBOrganelle(enc, orgCfg.nucleus, cp, cs, uni, toff);
                drawGLBOrganelle(enc, orgCfg.smoothER, cp, cs, uni, toff+1);
                drawGLBOrganelle(enc, orgCfg.roughER, cp, cs, uni, toff+2);
                drawGLBOrganelle(enc, orgCfg.golgi, cp, cs, uni, toff+3);

                // Render every mitochondrion in the cell at its biological
                // size. HeLa has 383-882 per cell (Posakony 1977 JCB) —
                // they fit because each is only ~2 µm × 0.5 µm ≈ 0.4 µm³,
                // and 500 of them = 200 µm³ ≈ 5% of a 4000 µm³ cell
                // (Alberts Ch 14 — mitochondria occupy ~5-10% of cytoplasm
                // volume). So crowded, but not overlapping.
                //
                // Render only a visible SAMPLE of mitochondria. The real
                // count (hundreds) is tracked as a scalar on SimCell and
                // drives the biology (ATP output, etc.). Rendering all
                // 500 physically-sized mitos blocks the view of DNA,
                // chromosomes, and cytoplasm. 60 is enough to read as a
                // "mito network" without crowding the interior.
                int mitoCount = SimCell::MITO_TARGET_DEFAULT;
                if (sourceIdx >= 0 && sourceIdx < (int)gSim.cells.size()) {
                    mitoCount = gSim.cells[sourceIdx].mitoCount;
                }
                float camDx = cp.x - camPos.x, camDy = cp.y - camPos.y, camDz = cp.z - camPos.z;
                float camDist = sqrtf(camDx*camDx + camDy*camDy + camDz*camDz);
                int maxForDist = 60;
                if (camDist > 25.0f) maxForDist = 25;
                if (camDist > 45.0f) maxForDist = 10;
                int visibleMitos = (int)sqrtf((float)mitoCount) * 2;
                if (visibleMitos > maxForDist) visibleMitos = maxForDist;
                if (visibleMitos < 3) visibleMitos = 3;
                for (int m = 0; m < visibleMitos; m++) {
                    OrganelleConfig cfg = orgCfg.mito[m % 3];
                    // 3D golden-angle sphere distribution. Radius range
                    // 0.55-0.88 × cellR so mitos occupy the CYTOPLASM
                    // shell (outside nucleus at 0.50) right up to the
                    // cortex — real mitochondrial network hugs the
                    // membrane in HeLa. Previous range (0.22-0.72) left
                    // an empty rim between organelles and membrane.
                    const float PHI = 2.39996323f;
                    float phi = (float)m * PHI + (float)(i % 17) * 0.41f;
                    float ySlice = ((float)m + 0.5f) / (float)visibleMitos * 2.0f - 1.0f;
                    float yScale = sqrtf(1.0f - ySlice * ySlice);
                    float rAlt = 0.55f + 0.33f * fmodf((float)m * 0.271828f + (float)(i * 13 % 97) * 0.0103f, 1.0f);
                    float yJitter = 0.85f + 0.15f * sinf((float)m * 1.61f + (float)(i % 23));
                    cfg.position.x = cosf(phi) * yScale * rAlt;
                    cfg.position.y = ySlice * rAlt * 0.70f * yJitter;
                    cfg.position.z = sinf(phi) * yScale * rAlt;
                    cfg.rotation.y += (float)m * 0.27f;
                    drawGLBOrganelle(enc, cfg, cp, cs, uni, toff+4+m);
                }
            }
        }

        // 3. Cell membranes
        {
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            const MeshData& sphere = gMeshes.sphereLOD(2);
            [enc setRenderPipelineState:gCtx.cellPipeline()];
            [enc setVertexBuffer:sphere.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBuffer:gCellBuffer offset:0 atIndex:1];
            [enc setVertexBytes:&uni length:sizeof(Uniforms) atIndex:2];
            [enc setFragmentBytes:&uni length:sizeof(Uniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle indexCount:sphere.indexCount
                             indexType:MTLIndexTypeUInt32 indexBuffer:sphere.indexBuffer
                     indexBufferOffset:0 instanceCount:cellCount];
        }

        // 4. Wireframe
        {
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            const MeshData& sphere = gMeshes.sphereLOD(1);
            [enc setRenderPipelineState:gCtx.wirePipeline()];
            [enc setVertexBuffer:sphere.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBuffer:gCellBuffer offset:0 atIndex:1];
            [enc setVertexBytes:&uni length:sizeof(Uniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeLine indexCount:sphere.indexCount
                             indexType:MTLIndexTypeUInt32 indexBuffer:sphere.indexBuffer
                     indexBufferOffset:0 instanceCount:cellCount];
        }

        // 4a. Pierce-through glow pass for the SELECTED cell only.
        // Depth test always-pass so the organelle outline shows through
        // the cell membrane and any obstructing geometry. Mitochondria
        // intentionally skipped (per user — keeps the existing mito look).
        {
            int focusIdx = activeFocusCellIndex();
            if (focusIdx >= 0 && focusIdx < cellCount) {
                [enc setDepthStencilState:gCtx.depthStateAlways()];
                [enc setRenderPipelineState:gCtx.glbOrganellePipeline()];
                simd_float3 cp = gCellInstances[focusIdx].position;
                float cs = gCellInstances[focusIdx].radius;
                float toff = time + (float)focusIdx * 7.3f;
                int orgKey = organelleMotionKeyForRenderIndex(focusIdx);
                OrganelleConfigs orgCfg = animatedOrganelleConfigsForCell(orgKey, time);
                // Glow boost makes the sinusoidal pulse strong enough to
                // dominate the colour. Nucleus extra-bright per user ask.
                drawOrganelleGlowPass(enc, uni, cp, cs, orgCfg.nucleus,  toff,     4.5f);
                drawOrganelleGlowPass(enc, uni, cp, cs, orgCfg.smoothER, toff + 1, 3.0f);
                drawOrganelleGlowPass(enc, uni, cp, cs, orgCfg.roughER,  toff + 2, 3.0f);
                drawOrganelleGlowPass(enc, uni, cp, cs, orgCfg.golgi,    toff + 3, 3.0f);
            }
        }

        // 5. Molecules (single cell mode — render inside the cell, before membrane)
        if (gSimMode == MODE_SINGLE_CELL && gMolRender.atomCount > 0
            && gCtx.moleculeAtomPipeline()) {
            struct MolUniforms {
                simd_float4x4 viewProjection;
                simd_float3   cameraPos;
                float         time;
                simd_float3   lightDir;
                float         pad0;
                simd_float3   cellCenter;
                float         cellRadius;
            };
            MolUniforms mu;
            mu.viewProjection = uni.viewProjection;
            mu.cameraPos = uni.cameraPos;
            mu.time = time;
            mu.lightDir = uni.lightDir;
            mu.pad0 = 0;
            // The CPU-side intracellular physics already keeps particles inside
            // their own cell. A single global GPU membrane clamp breaks as soon
            // as we render multiple daughter interiors in one atom buffer.
            mu.cellCenter = {0, 0, 0};
            mu.cellRadius = 0.0f;

            // Atoms (instanced spheres) — depth write ON so molecules occlude each other
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            [enc setRenderPipelineState:gCtx.moleculeAtomPipeline()];
            const MeshData& atomSphere = gMeshes.sphereLOD(1);
            [enc setVertexBuffer:atomSphere.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBuffer:gMolRender.atomBuffer offset:0 atIndex:1];
            [enc setVertexBytes:&mu length:sizeof(MolUniforms) atIndex:2];
            [enc setFragmentBytes:&mu length:sizeof(MolUniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                            indexCount:atomSphere.indexCount
                             indexType:MTLIndexTypeUInt32
                           indexBuffer:atomSphere.indexBuffer
                     indexBufferOffset:0
                         instanceCount:gMolRender.atomCount];

            // Bonds (instanced cylinders)
            if (gMolRender.bondCount > 0 && gCtx.moleculeBondPipeline()) {
                [enc setRenderPipelineState:gCtx.moleculeBondPipeline()];
                const MeshData& cyl = gMeshes.cylinder();
                [enc setVertexBuffer:cyl.vertexBuffer offset:0 atIndex:0];
                [enc setVertexBuffer:gMolRender.bondBuffer offset:0 atIndex:1];
                [enc setVertexBytes:&mu length:sizeof(MolUniforms) atIndex:2];
                [enc setFragmentBytes:&mu length:sizeof(MolUniforms) atIndex:2];
                [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                                indexCount:cyl.indexCount
                                 indexType:MTLIndexTypeUInt32
                               indexBuffer:cyl.indexBuffer
                         indexBufferOffset:0
                             instanceCount:gMolRender.bondCount];
            }
        }

        } // end if(!gSetup.showOverlay)

        // 6. ImGui render (always draws — this is the Setup overlay too)
        {
            ImGui_ImplMetal_NewFrame(pass);
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            drawUI();
            drawMitosisDebugOverlay(uni.viewProjection);
            ImGui::Render();
            ImGui_ImplMetal_RenderDrawData(ImGui::GetDrawData(), cmd, enc);
        }

        [enc endEncoding];
        [cmd presentDrawable:drawable];
        [cmd commit];
    }
}

// ── Main ────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    @autoreleasepool {
        printf("╔══════════════════════════════════════════════════════════╗\n");
        printf("║  CellSim — Computational Cell Biology Simulator    ║\n");
        printf("║  Full Simulation · CDK/Cyclin ODE · Fick Diffusion     ║\n");
        printf("║  Fate Decision · Cell Division · ImGui Research UI     ║\n");
        printf("╚══════════════════════════════════════════════════════════╝\n\n");

        NSString* execPath = [[NSBundle mainBundle] executablePath];
        NSString* execDir  = [execPath stringByDeletingLastPathComponent];
        NSFileManager* fm = [NSFileManager defaultManager];
        NSString* candidateA = [execDir stringByDeletingLastPathComponent];
        NSString* candidateB = [candidateA stringByDeletingLastPathComponent];
        NSString* cwdPath = [[fm currentDirectoryPath] stringByStandardizingPath];
        auto looksLikeProjectRoot = ^bool(NSString* path) {
            return [fm fileExistsAtPath:[path stringByAppendingPathComponent:@"assets"]] &&
                   [fm fileExistsAtPath:[path stringByAppendingPathComponent:@"src"]];
        };
        NSString* projectDir = nil;
        if (looksLikeProjectRoot(candidateA)) projectDir = candidateA;
        else if (looksLikeProjectRoot(candidateB)) projectDir = candidateB;
        else if (looksLikeProjectRoot(cwdPath)) projectDir = cwdPath;
        else projectDir = candidateA;
        initializeSessionLogging([projectDir UTF8String]);
        NSString* modelDir = [[execDir stringByDeletingLastPathComponent]
                              stringByAppendingPathComponent:@"assets/models"];
        gModelDir = [modelDir UTF8String];

        printf("[CellSim] Session: %s\n", gRunSessionTag.c_str());
        printf("[CellSim] Launch time: %s\n", gRunLaunchTimeText.c_str());
        printf("[CellSim] Mitosis log: %s\n", gMitosisLogPath.c_str());
        printf("[CellSim] Diag log: %s\n", gDiagLogPath.c_str());
        printf("[CellSim] Centrosome log: %s\n\n", gCentrosomeLogPath.c_str());

        // Export directory = project exports folder
        NSString* exportDir = [[execDir stringByDeletingLastPathComponent]
                               stringByAppendingPathComponent:@"exports"];
        [[NSFileManager defaultManager] createDirectoryAtPath:exportDir
            withIntermediateDirectories:YES attributes:nil error:nil];
        gExportDir = [exportDir UTF8String];

        if (!glfwInit()) return 1;
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        std::string windowTitle = "CellSim — Computational Cell Biology Simulator [" + gRunSessionTag + "]";
        GLFWwindow* window = glfwCreateWindow(INIT_WIDTH, INIT_HEIGHT,
            windowTitle.c_str(), nullptr, nullptr);
        if (!window) { glfwTerminate(); return 1; }

        NSWindow* nsWin = glfwGetCocoaWindow(window);
        CAMetalLayer* metalLayer = [CAMetalLayer layer];
        metalLayer.contentsScale = nsWin.backingScaleFactor;
        nsWin.contentView.layer = metalLayer;
        nsWin.contentView.wantsLayer = YES;

        int fbW, fbH;
        glfwGetFramebufferSize(window, &fbW, &fbH);
        metalLayer.drawableSize = CGSizeMake(fbW, fbH);

        if (!gCtx.init(metalLayer)) { glfwTerminate(); return 1; }
        gCtx.recreateDepthTexture(fbW, fbH);
        if (!gMeshes.init(gCtx.device())) { glfwTerminate(); return 1; }
        gCamera.onResize(fbW, fbH);

        // ImGui setup
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGui::StyleColorsDark();
        ImGuiStyle& style = ImGui::GetStyle();
        style.WindowRounding = 4; style.FrameRounding = 2;
        style.Colors[ImGuiCol_WindowBg] = ImVec4(0.01f, 0.03f, 0.06f, 0.85f);
        style.Colors[ImGuiCol_TitleBg] = ImVec4(0.0f, 0.04f, 0.08f, 0.9f);
        style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.0f, 0.06f, 0.12f, 0.95f);
        style.Colors[ImGuiCol_FrameBg] = ImVec4(0.0f, 0.04f, 0.08f, 0.7f);
        style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.0f, 0.5f, 0.8f, 1.0f);

        // Set our callbacks FIRST, then ImGui chains on top with install_callbacks=true
        glfwSetMouseButtonCallback(window, mouseButtonCB);
        glfwSetCursorPosCallback(window, cursorPosCB);
        glfwSetScrollCallback(window, scrollCB);
        glfwSetFramebufferSizeCallback(window, framebufferCB);

        ImGui_ImplGlfw_InitForOther(window, true);
        ImGui_ImplMetal_Init(gCtx.device());

        // Load models
        loadModels();

        // Init molecule cache
        gProjectRoot = [projectDir UTF8String];
        NSString* molDir = [[execDir stringByDeletingLastPathComponent]
                            stringByAppendingPathComponent:@"assets/molecules"];
        gMoleculeDir = [molDir UTF8String];
        gMolCache.init(gMoleculeDir);
        gMolGen.init(gProjectRoot);

        // Init protein structure cache
        NSString* protDir = [[execDir stringByDeletingLastPathComponent]
                             stringByAppendingPathComponent:@"assets/proteins"];
        gProteinDir = [protDir UTF8String];
        gProtCache.init(gProteinDir);

        // NO simulation init yet — we show the Setup overlay on the very
        // first frame and only build the sim when the user clicks Start.
        // Until then the dish is empty: no cells, no fluid, no molecules.
        applySetupTemplate(0);            // seed default field values only
        gSim.paused = true;
        gSelectedCell = 0;
        gSelectedCellUid = -1;
        gPostMitosisPairCameraTimer = 0.0f;
        gPostMitosisPairUidA = -1;
        gPostMitosisPairUidB = -1;
        gCamera.setSingleCellView();

        int initCount = (gSimMode == MODE_SINGLE_CELL) ? SINGLE_CELL_COUNT : INIT_CELLS;
        printf("[CellSim] Simulation started: %d cell(s) · %d models loaded · Mode: %s\n",
               initCount, (int)gOrganelleMeshes.size(),
               gSimMode == MODE_SINGLE_CELL ? "Single Cell" : "Colony");
        printf("[CellSim] Drag:rotate  Right-drag:pan  Scroll:zoom  ESC:quit\n");

        float lastTime = (float)glfwGetTime();
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
                glfwSetWindowShouldClose(window, true);
            // P to pause (Space is now fly-up)
            static bool pWasDown = false;
            bool pDown = glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS;
            if (pDown && !pWasDown) gSim.paused = !gSim.paused;
            pWasDown = pDown;
            // R to re-open the startup Setup overlay
            static bool rWasDown = false;
            bool rDown = glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS;
            if (rDown && !rWasDown && !ImGui::GetIO().WantCaptureKeyboard) {
                gSetup.showOverlay = true;
            }
            rWasDown = rDown;

            float now = (float)glfwGetTime();
            float dt = fminf(now - lastTime, 0.05f);
            lastTime = now;

            // WASD has two modes:
            //  - Follow mode ON: WASD switches to the nearest cell in that
            //    screen-space direction (edge-triggered, one jump per press).
            //  - Follow mode OFF: WASD flies the free camera.
            if (!ImGui::GetIO().WantCaptureKeyboard) {
                static bool wPrev=false, aPrev=false, sPrev=false, dPrev=false, fPrev=false;
                bool wDown = glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS;
                bool aDown = glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS;
                bool sDown = glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS;
                bool dDown = glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS;
                bool fDown = glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS;

                // F toggles follow mode (edge-triggered)
                if (fDown && !fPrev) gFollowCell = !gFollowCell;
                fPrev = fDown;

                if (gFollowCell && gSimMode == MODE_SINGLE_CELL && gSim.cells.size() > 1) {
                    int curIdx = activeFocusCellIndex();
                    simd_float3 fwd = gCamera.getForward();
                    // Build camera-aligned directions on the ground plane (y constant)
                    simd_float3 fwdPlanar = {fwd.x, 0, fwd.z};
                    float fl = sqrtf(fwdPlanar.x*fwdPlanar.x + fwdPlanar.z*fwdPlanar.z);
                    if (fl > 1e-4f) { fwdPlanar.x/=fl; fwdPlanar.z/=fl; }
                    // Right = forward × world-up (right-handed: +X when facing -Z)
                    simd_float3 rightPlanar = {-fwdPlanar.z, 0, fwdPlanar.x};

                    if (wDown && !wPrev) {
                        int next = nearestCellInDirection(curIdx, fwdPlanar);
                        if (next != curIdx) selectCellIndex(next);
                    }
                    if (sDown && !sPrev) {
                        simd_float3 back = {-fwdPlanar.x, 0, -fwdPlanar.z};
                        int next = nearestCellInDirection(curIdx, back);
                        if (next != curIdx) selectCellIndex(next);
                    }
                    if (dDown && !dPrev) {
                        int next = nearestCellInDirection(curIdx, rightPlanar);
                        if (next != curIdx) selectCellIndex(next);
                    }
                    if (aDown && !aPrev) {
                        simd_float3 left = {-rightPlanar.x, 0, -rightPlanar.z};
                        int next = nearestCellInDirection(curIdx, left);
                        if (next != curIdx) selectCellIndex(next);
                    }
                } else {
                    // Follow mode off: WASD flies the camera normally.
                    gCamera.updateFPS(dt,
                        wDown, aDown, sDown, dDown,
                        glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS,
                        glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS ||
                        glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS);
                }
                wPrev=wDown; aPrev=aDown; sPrev=sDown; dPrev=dDown;
            }

            // In single-cell mode, stay locked to one selected cell so mitosis
            // details remain large enough to inspect frame-by-frame.
            if (gFollowCell && !gSim.cells.empty()) {
                refreshTrackedPostMitosisPairFromSimulation();
                if (trackedPostMitosisPairActive()) {
                    if (gPostMitosisPairCameraTimer > 0.0f) {
                        gPostMitosisPairCameraTimer = fmaxf(0.0f, gPostMitosisPairCameraTimer - dt);
                    }
                    simd_float3 center = {
                        (gPostMitosisPairA.x + gPostMitosisPairB.x) * 0.5f,
                        (gPostMitosisPairA.y + gPostMitosisPairB.y) * 0.5f,
                        (gPostMitosisPairA.z + gPostMitosisPairB.z) * 0.5f
                    };
                    simd_float3 d = {
                        gPostMitosisPairA.x - gPostMitosisPairB.x,
                        gPostMitosisPairA.y - gPostMitosisPairB.y,
                        gPostMitosisPairA.z - gPostMitosisPairB.z
                    };
                    float halfSep = 0.5f * sqrtf(d.x*d.x + d.y*d.y + d.z*d.z);
                    gCamera.followCell(center, gPostMitosisPairRadius + halfSep * 1.25f, gFollowZoom);
                } else if (gSimMode == MODE_SINGLE_CELL && gSim.cells.size() == 2) {
                    const auto& a = gSim.cells[0];
                    const auto& b = gSim.cells[1];
                    simd_float3 center = {
                        (a.position.x + b.position.x) * 0.5f,
                        (a.position.y + b.position.y) * 0.5f,
                        (a.position.z + b.position.z) * 0.5f
                    };
                    simd_float3 d = {
                        a.position.x - b.position.x,
                        a.position.y - b.position.y,
                        a.position.z - b.position.z
                    };
                    float halfSep = 0.5f * sqrtf(d.x*d.x + d.y*d.y + d.z*d.z);
                    float pairRadius = fmaxf(a.radius * a.size, b.radius * b.size);
                    gCamera.followCell(center, pairRadius + halfSep * 1.25f, gFollowZoom);
                } else {
                    int focusIdx = activeFocusCellIndex();
                    auto& c = gSim.cells[focusIdx];
                    gCamera.followCell(c.position, c.radius * c.size, gFollowZoom);
                }
            }

            renderFrame(now, dt);
        }

        ImGui_ImplMetal_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        gCtx.shutdown();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 0;
    }
}
