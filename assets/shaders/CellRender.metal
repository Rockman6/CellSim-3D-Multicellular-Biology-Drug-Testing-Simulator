#include <metal_stdlib>
using namespace metal;

// ══════════════════════════════════════════════════════════════════════════
//  CellRender.metal — Translucent bio-luminescent cell rendering
//  Bio-luminescent visual style:
//    - Transparent membrane (opacity 0.12-0.28)
//    - Interior point light glow (subsurface scattering approximation)
//    - Fresnel rim highlighting
//    - Wireframe overlay
//    - Procedural organelle spheres (nucleus, mito, ER)
// ══════════════════════════════════════════════════════════════════════════

struct Vertex {
    float3 position;
    float3 normal;
    float2 uv;
};

struct Uniforms {
    float4x4 viewProjection;
    float4x4 model;
    float3   cameraPos;
    float    time;
    float3   lightDir;
    float    pad0;
};

struct CellInstance {
    float3   position;
    float    radius;
    float4   color;          // membrane tint color (low alpha)
    float    glowIntensity;  // interior light strength
    uint     lodLevel;       // 0-3
    float    phase;          // 0-3 (G1/S/G2/M) for interior color
    float    pad;
};

// ── Vertex outputs ──────────────────────────────────────────────────────

struct CellVertexOut {
    float4 clipPos   [[position]];
    float3 worldPos;
    float3 localPos;         // position on unit sphere
    float3 normal;
    float4 memColor;         // membrane tint
    float3 interiorColor;    // glow color
    float  glowIntensity;
    float  memRadius;
    float  time;
};

// ── FATE COLOR TABLE ──────────────────────────

// Interior glow colors by phase
constant float3 PHASE_GLOW_COLORS[] = {
    float3(0.35, 0.55, 1.0),   // G1 — cool blue
    float3(0.2,  1.0,  0.6),   // S  — green (DNA synthesis)
    float3(0.9,  0.75, 0.2),   // G2 — amber
    float3(1.0,  0.27, 0.67),  // M  — magenta (mitosis)
};

// Membrane opacity by phase
constant float PHASE_MEM_OPACITY[] = {
    0.17,  // G1 (UNDETERMINED-like)
    0.22,  // S  (PROLIFERATING-like)
    0.17,  // G2
    0.28,  // M  (dividing, more opaque)
};

// ══════════════════════════════════════════════════════════════════════════
//  CELL MEMBRANE — translucent outer shell with bio-luminescent interior
// ══════════════════════════════════════════════════════════════════════════

vertex CellVertexOut cellVertex(
    uint        vid   [[vertex_id]],
    uint        iid   [[instance_id]],
    const device Vertex*       vertices  [[buffer(0)]],
    const device CellInstance* instances [[buffer(1)]],
    constant     Uniforms&     u         [[buffer(2)]]
) {
    Vertex v = vertices[vid];
    CellInstance inst = instances[iid];

    // Organic membrane deformation (spherical harmonic noise)
    float3 p = v.position; // unit sphere
    float seed = float(iid) * 7.3;
    float theta = atan2(p.z, p.x);
    float phi   = acos(clamp(p.y, -1.0f, 1.0f));
    float noise = 0.055 * (
        sin(2.0*theta + seed) * cos(phi + seed*1.3) +
        0.65 * sin(3.0*theta + seed*2.1) * cos(2.0*phi + seed*0.7) +
        0.4  * sin(theta + seed*0.9) * cos(3.0*phi + seed*1.8) +
        0.25 * cos(4.0*theta + seed*3.3) * cos(phi*2.0 + seed*2.5)
    );
    // Gentle breathing animation
    float breathe = 1.0 + 0.015 * sin(u.time * 1.2 + seed * 0.3);
    float r = (1.0 + noise) * breathe;

    float3 deformedPos = normalize(p) * r;
    float3 worldPos = deformedPos * inst.radius + inst.position;

    // Phase-dependent colors
    uint phaseIdx = clamp(uint(inst.phase), 0u, 3u);
    float baseOpacity = PHASE_MEM_OPACITY[phaseIdx];

    CellVertexOut out;
    out.clipPos        = u.viewProjection * float4(worldPos, 1.0);
    out.worldPos       = worldPos;
    out.localPos       = deformedPos;
    out.normal         = normalize(deformedPos); // organic normal
    out.memColor       = float4(inst.color.rgb, baseOpacity);
    out.interiorColor  = PHASE_GLOW_COLORS[phaseIdx];
    out.glowIntensity  = inst.glowIntensity;
    out.memRadius      = inst.radius;
    out.time           = u.time;
    return out;
}

fragment float4 cellFragment(
    CellVertexOut    in [[stage_in]],
    constant Uniforms& u  [[buffer(2)]]
) {
    float3 N = normalize(in.normal);
    float3 V = normalize(u.cameraPos - in.worldPos);
    float3 L = normalize(u.lightDir);
    float3 H = normalize(L + V);

    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float NdotH = max(dot(N, H), 0.0);

    // ── Fresnel (Schlick approximation) ─────────────────────────────
    // Strong Fresnel = cells glow at edges, transparent in center
    float fresnel = pow(1.0 - NdotV, 3.5);

    // ── Membrane surface ────────────────────────────────────────────
    // Very subtle diffuse (membrane is mostly transparent)
    float diffuse = NdotL * 0.15 + 0.05;

    // Specular highlights (wet membrane surface)
    float spec = pow(NdotH, 80.0) * 0.5;

    // ── Interior glow (subsurface scattering approximation) ─────────
    // Light passes through the translucent membrane, revealing the interior
    // Brighter at center (where organelles are), fades at edges
    float centerWeight = pow(NdotV, 0.6); // brighter when looking straight through
    float interiorGlow = in.glowIntensity * centerWeight * 0.7;

    // Pulsing interior (organelle activity simulation)
    float pulse = 0.85 + 0.15 * sin(in.time * 2.0 + in.localPos.x * 3.0 + in.localPos.y * 2.0);
    interiorGlow *= pulse;

    // ── Combine ─────────────────────────────────────────────────────
    float3 memSurface = in.memColor.rgb * diffuse;
    float3 rimGlow    = in.memColor.rgb * fresnel * 0.6;
    float3 interior   = in.interiorColor * interiorGlow;
    float3 specColor  = float3(0.3, 0.6, 0.9) * spec;

    float3 finalColor = memSurface + rimGlow + interior + specColor;

    // ── Alpha: transparent membrane with visible rim ────────────────
    float alpha = in.memColor.a + fresnel * 0.35 + interiorGlow * 0.15;
    alpha = clamp(alpha, 0.0, 0.85);

    // ── Fog ───────────────────────────────────
    float dist = length(u.cameraPos - in.worldPos);
    float fog  = 1.0 - exp(-dist * 0.012);
    finalColor = mix(finalColor, float3(0.004, 0.008, 0.031), fog);
    alpha *= (1.0 - fog * 0.5);

    return float4(finalColor, alpha);
}

// ══════════════════════════════════════════════════════════════════════════
//  WIREFRAME OVERLAY — thin blue grid over membrane
// ══════════════════════════════════════════════════════════════════════════

struct WireVertexOut {
    float4 clipPos [[position]];
    float4 color;
};

vertex WireVertexOut wireVertex(
    uint        vid   [[vertex_id]],
    uint        iid   [[instance_id]],
    const device Vertex*       vertices  [[buffer(0)]],
    const device CellInstance* instances [[buffer(1)]],
    constant     Uniforms&     u         [[buffer(2)]]
) {
    Vertex v = vertices[vid];
    CellInstance inst = instances[iid];

    float3 worldPos = v.position * inst.radius * 1.02 + inst.position;

    WireVertexOut out;
    out.clipPos = u.viewProjection * float4(worldPos, 1.0);

    // Wire color: faint cyan
    float dist = length(u.cameraPos - worldPos);
    float fogFade = exp(-dist * 0.015);
    float wireAlpha = 0.06 * fogFade;
    out.color = float4(0.0, 0.87, 1.0, wireAlpha);
    return out;
}

fragment float4 wireFragment(WireVertexOut in [[stage_in]]) {
    return in.color;
}

// ══════════════════════════════════════════════════════════════════════════
//  PROCEDURAL ORGANELLES — nucleus, mitochondria, ER rendered as spheres
//  Until GLB models are available, we render glowing spheres inside cells
// ══════════════════════════════════════════════════════════════════════════

struct OrganelleInstance {
    float3   position;       // world position
    float    radius;
    float4   color;          // emissive color
    float    glowStrength;
    float    pad[3];
};

struct OrgVertexOut {
    float4 clipPos [[position]];
    float3 worldPos;
    float3 normal;
    float4 color;
    float  glow;
};

vertex OrgVertexOut organelleVertex(
    uint        vid   [[vertex_id]],
    uint        iid   [[instance_id]],
    const device Vertex*            vertices  [[buffer(0)]],
    const device OrganelleInstance*  instances [[buffer(1)]],
    constant     Uniforms&          u         [[buffer(2)]]
) {
    Vertex v = vertices[vid];
    OrganelleInstance inst = instances[iid];

    float3 worldPos = v.position * inst.radius + inst.position;

    OrgVertexOut out;
    out.clipPos  = u.viewProjection * float4(worldPos, 1.0);
    out.worldPos = worldPos;
    out.normal   = v.normal;
    out.color    = inst.color;
    out.glow     = inst.glowStrength;
    return out;
}

fragment float4 organelleFragment(
    OrgVertexOut     in [[stage_in]],
    constant Uniforms& u  [[buffer(2)]]
) {
    float3 N = normalize(in.normal);
    float3 V = normalize(u.cameraPos - in.worldPos);
    float3 L = normalize(u.lightDir);

    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);

    float fresnel = pow(1.0 - NdotV, 2.5);

    float3 emissive = in.color.rgb * in.glow;
    float3 diffuse  = in.color.rgb * NdotL * 0.3;
    float3 rim      = in.color.rgb * fresnel * 0.5;

    float3 finalColor = emissive + diffuse + rim;

    float dist = length(u.cameraPos - in.worldPos);
    float fog  = 1.0 - exp(-dist * 0.012);
    finalColor = mix(finalColor, float3(0.004, 0.008, 0.031), fog);

    float alpha = in.color.a * (0.6 + 0.4 * in.glow) + fresnel * 0.2;
    alpha = clamp(alpha, 0.0, 0.9);

    return float4(finalColor, alpha);
}

// ══════════════════════════════════════════════════════════════════════════
//  SUBSTRATE — floor + grid (same as before but brighter grid)
// ══════════════════════════════════════════════════════════════════════════

struct SubstrateVertexOut {
    float4 position [[position]];
    float3 worldPos;
    float2 uv;
};

vertex SubstrateVertexOut substrateVertex(
    uint         vid  [[vertex_id]],
    const device Vertex*   vertices [[buffer(0)]],
    constant     Uniforms& u        [[buffer(2)]]
) {
    Vertex v = vertices[vid];
    float4 wp = u.model * float4(v.position, 1.0);

    SubstrateVertexOut out;
    out.position = u.viewProjection * wp;
    out.worldPos = wp.xyz;
    out.uv       = v.uv;
    return out;
}

fragment float4 substrateFragment(
    SubstrateVertexOut in [[stage_in]],
    constant Uniforms& u  [[buffer(2)]]
) {
    // Grid pattern
    float2 grid = abs(fract(in.uv * 12.0) - 0.5);
    float  line = min(grid.x, grid.y);
    float  gridAlpha = 1.0 - smoothstep(0.0, 0.03, line);

    // Radial fade (petri dish edge)
    float dist = length(in.worldPos.xz);
    float radialFade = 1.0 - smoothstep(45.0, 55.0, dist);

    // Bright rim at petri dish edge
    float rimGlow = smoothstep(50.0, 54.0, dist) * (1.0 - smoothstep(54.0, 55.0, dist));

    // Base color (dark blue-black)
    float3 baseColor = float3(0.004, 0.012, 0.022);
    float3 gridColor = float3(0.0, 0.08, 0.12);
    float3 rimColor  = float3(0.0, 0.15, 0.25);

    float3 finalColor = mix(baseColor, gridColor, gridAlpha * 0.5 * radialFade)
                       + rimColor * rimGlow * 0.5;

    // Fog
    float camDist = length(u.cameraPos - in.worldPos);
    float fog = 1.0 - exp(-camDist * 0.012);
    finalColor = mix(finalColor, float3(0.004, 0.008, 0.031), fog);

    return float4(finalColor, 1.0);
}
