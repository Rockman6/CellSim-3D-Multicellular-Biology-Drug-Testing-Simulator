#include <metal_stdlib>
using namespace metal;

// ══════════════════════════════════════════════════════════════════════════
//  MoleculeRender.metal — Ball-and-stick rendering for molecular structures
//  Bio-luminescent style matching cell aesthetic
//  Atoms: instanced spheres with CPK coloring + glow + subsurface
//  Bonds: instanced cylinders with emissive edges
// ══════════════════════════════════════════════════════════════════════════

struct Vertex {
    float3 position;
    float3 normal;
    float2 uv;
};

struct MolUniforms {
    float4x4 viewProjection;
    float3   cameraPos;
    float    time;
    float3   lightDir;
    float    pad0;
    float3   cellCenter;     // World position of cell center
    float    cellRadius;     // Visual membrane radius
};

// Per-atom instance data
struct AtomInstance {
    float3 position;    // World position
    float  radius;      // Display radius
    float3 color;       // CPK color
    float  pad;
};

// Per-bond instance data
struct BondInstance {
    float3 posA;        // Start atom position
    float  radius;      // Bond cylinder radius
    float3 posB;        // End atom position
    float  pad;
    float3 color;       // Bond color
    float  pad2;
};

struct MolVertexOut {
    float4 clipPos   [[position]];
    float3 worldPos;
    float3 worldNormal;
    float3 color;
    float  time;
};

// ── Atom (sphere) vertex shader ─────────────────────────────────────────
vertex MolVertexOut moleculeAtomVertex(
    uint                    vid    [[vertex_id]],
    uint                    iid    [[instance_id]],
    const device Vertex*    verts  [[buffer(0)]],
    const device AtomInstance* atoms [[buffer(1)]],
    constant MolUniforms&   u      [[buffer(2)]]
) {
    Vertex v = verts[vid];
    AtomInstance atom = atoms[iid];

    float3 worldPos = v.position * atom.radius + atom.position;
    float3 worldNormal = normalize(v.normal);

    // ── MEMBRANE WALL ───────────────────────────────────────────
    if (u.cellRadius > 0.01) {
        float3 toCell = worldPos - u.cellCenter;
        float d = length(toCell);
        float maxDist = u.cellRadius * 0.93;
        if (d > maxDist && d > 0.001) {
            worldPos = u.cellCenter + normalize(toCell) * maxDist;
        }
    }

    MolVertexOut out;
    out.clipPos = u.viewProjection * float4(worldPos, 1.0);
    out.worldPos = worldPos;
    out.worldNormal = worldNormal;
    out.color = atom.color;
    out.time = u.time;
    return out;
}

// ── Bond (cylinder) vertex shader ───────────────────────────────────────
vertex MolVertexOut moleculeBondVertex(
    uint                    vid    [[vertex_id]],
    uint                    iid    [[instance_id]],
    const device Vertex*    verts  [[buffer(0)]],
    const device BondInstance* bonds [[buffer(1)]],
    constant MolUniforms&   u      [[buffer(2)]]
) {
    Vertex v = verts[vid];
    BondInstance bond = bonds[iid];

    float3 dir = bond.posB - bond.posA;
    float len = length(dir);
    if (len < 0.001) len = 0.001;
    float3 up = dir / len;

    float3 arbitrary = abs(up.y) < 0.99 ? float3(0,1,0) : float3(1,0,0);
    float3 right = normalize(cross(up, arbitrary));
    float3 forward = cross(right, up);

    float3 localPos = v.position;
    float3 worldPos = bond.posA + dir * 0.5
                    + right * localPos.x * bond.radius
                    + up * localPos.y * len * 0.5
                    + forward * localPos.z * bond.radius;

    float3 worldNormal = normalize(
        right * v.normal.x + up * v.normal.y + forward * v.normal.z
    );

    // ── MEMBRANE WALL ───────────────────────────────────────────
    if (u.cellRadius > 0.01) {
        float3 toCell = worldPos - u.cellCenter;
        float d = length(toCell);
        float maxDist = u.cellRadius * 0.93;
        if (d > maxDist && d > 0.001) {
            worldPos = u.cellCenter + normalize(toCell) * maxDist;
        }
    }

    MolVertexOut out;
    out.clipPos = u.viewProjection * float4(worldPos, 1.0);
    out.worldPos = worldPos;
    out.worldNormal = worldNormal;
    out.color = bond.color;
    out.time = u.time;
    return out;
}

// ── Fragment shader — bio-luminescent molecular style ───────────────────
fragment float4 moleculeFragment(
    MolVertexOut in [[stage_in]],
    constant MolUniforms& u [[buffer(2)]]
) {
    float3 N = normalize(in.worldNormal);
    float3 V = normalize(u.cameraPos - in.worldPos);
    float3 L = normalize(u.lightDir);
    float3 H = normalize(L + V);

    float NdotL = max(dot(N, L), 0.0);
    float NdotH = max(dot(N, H), 0.0);
    float NdotV = max(dot(N, V), 0.0);

    // Fresnel rim — glowing edges like organelles
    float fresnel = pow(1.0 - NdotV, 2.8);

    // Soft diffuse
    float diffuse = NdotL * 0.35 + 0.15;

    // Glossy specular (wet molecular surface)
    float spec = pow(NdotH, 80.0) * 0.5;

    // Emissive glow — molecules glow softly from within
    float pulse = 0.88 + 0.12 * sin(in.time * 1.8 + in.worldPos.x * 4.0 + in.worldPos.z * 3.0);
    float3 emissive = in.color * 0.4 * pulse;

    // Subsurface scattering approximation — light through translucent molecule
    float sss = pow(max(dot(-L, V), 0.0), 3.0) * 0.25;
    float3 subsurface = in.color * sss;

    // Combine
    float3 finalColor = in.color * diffuse
                      + emissive
                      + float3(0.4, 0.6, 0.9) * spec
                      + in.color * fresnel * 0.4
                      + subsurface;

    // Fog (match cell fog)
    float dist = length(u.cameraPos - in.worldPos);
    float fog = 1.0 - exp(-dist * 0.012);
    finalColor = mix(finalColor, float3(0.004, 0.008, 0.031), fog);

    // Ghostly — base alpha < 10% so membrane mosaic and inner molecules
    // read as a translucent cloud the viewer sees through clearly.
    // Fresnel rim retains a faint sharp edge so each dot stays visible.
    float alpha = 0.06 + fresnel * 0.10 + emissive.r * 0.05;
    alpha = clamp(alpha, 0.0, 0.22);
    alpha *= (1.0 - fog * 0.3);

    return float4(finalColor, alpha);
}
