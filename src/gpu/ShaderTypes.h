#ifndef SHADER_TYPES_H
#define SHADER_TYPES_H

#include <simd/simd.h>

// ══════════════════════════════════════════════════════════════════════════
//  ShaderTypes.h — Shared between C++ and Metal shaders
//  All structs must be plain C types with simd types for GPU compatibility
// ══════════════════════════════════════════════════════════════════════════

struct Vertex {
    simd_float3 position;
    simd_float3 normal;
    simd_float2 uv;
};

struct Uniforms {
    simd_float4x4 viewProjection;
    simd_float4x4 model;
    simd_float3   cameraPos;
    float          time;
    simd_float3   lightDir;
    float          pad0;
};

// Per-instance data for instanced cell rendering
struct CellInstance {
    simd_float3 position;
    float       radius;
    simd_float4 color;       // membrane tint color
    float       glowIntensity;
    uint32_t    lodLevel;
    float       phase;       // 0-3: G1/S/G2/M
    float       pad;
};

// Per-instance data for procedural organelles
struct OrganelleInstance {
    simd_float3 position;
    float       radius;
    simd_float4 color;       // emissive RGBA
    float       glowStrength;
    float       pad[3];
};

// Per-cell transform (GPU physics buffer)
struct CellTransform {
    simd_float3 position;
    float       radius;
    simd_float3 velocity;
    float       size;        // growth scaling factor [0.1 - 2.0]
};

// Per-cell biology state (GPU compute buffer)
struct CellBiology {
    // Ca²⁺ / IP₃ ODE state (Li-Rinzel model)
    float ip3;
    float ca_cyt;
    float ca_er;
    float h_gate;
    float pkc;
    float mapk;

    // CDK/Cyclin ODE state (Novak-Tyson model)
    float CycD;    // Cyclin D / CDK4/6
    float Rb;      // Retinoblastoma protein (active)
    float E2F;     // E2F transcription factor
    float CycE;    // Cyclin E / CDK2
    float CycA;    // Cyclin A / CDK2
    float CycB;    // Cyclin B / CDK1 (MPF)
    float p21;     // CDK inhibitor p21CIP1

    // Metabolism
    float ATP;
    float stress;
    float ROS;
    float mitoPotential;
    float UPR;
    float biomass;
    float misfoldedProtein;

    // Phase / Fate
    uint32_t phase;          // 0=G1, 1=S, 2=G2, 3=M
    uint32_t fate;           // 0=undetermined, 1=prolif, 2=quiescent, 3=apoptotic
    float    cycleTimer;
    float    cycleProgress;
    float    damageLevel;

    // Flags
    uint32_t metabolismMode; // 0=OXIDATIVE, 1=GLYCOLYTIC
    uint32_t flags;          // bitfield: necrotic, senescent, markedForRemoval, etc.
    float    age;
    float    pad;
};

// Per-cell genome (less frequently accessed)
struct CellGenome {
    float    telomereA[23];
    float    telomereB[23];
    uint32_t mutations[23];

    // Heritable traits
    float glycolysisBias;
    float rosTolerance;
    float prolifBias;
    float couplingStrength;
    float o2Efficiency;
    float repairRate;
    float elasticModulus;
    float invasionBias;

    uint32_t generation;
    uint32_t cloneId;
    float    fitnessScore;
    float    pad;
};

// Cell state flags (bitfield for CellBiology.flags)
enum CellFlags : uint32_t {
    CELL_FLAG_NECROTIC          = 1 << 0,
    CELL_FLAG_SENESCENT         = 1 << 1,
    CELL_FLAG_MARKED_FOR_REMOVAL = 1 << 2,
    CELL_FLAG_DIVIDING          = 1 << 3,
    CELL_FLAG_ADHERED           = 1 << 4,
};

// Cell cycle phases
enum CellPhase : uint32_t {
    PHASE_G1 = 0,
    PHASE_S  = 1,
    PHASE_G2 = 2,
    PHASE_M  = 3,
};

// Cell fate states
enum CellFate : uint32_t {
    FATE_UNDETERMINED  = 0,
    FATE_PROLIFERATING = 1,
    FATE_QUIESCENT     = 2,
    FATE_APOPTOTIC     = 3,
};

#endif // SHADER_TYPES_H
