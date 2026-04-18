#include <metal_stdlib>
using namespace metal;

// ══════════════════════════════════════════════════════════════════════════
//  FluidRender.metal — Continuous translucent culture-medium volume
//  -----------------------------------------------------------------------
//  Renders one cylinder mesh that fills the dish above the floor. The TOP
//  RING of vertices is displaced vertically by a sum-of-sines wave so the
//  surface looks like ocean. Inside the volume, the fragment shader draws
//  a soft water tint with fresnel-driven transparency (more visible at
//  grazing angles, near 100 % transparent looking straight down).
//
//  Vertex layout (fluidVertex):
//    position.xz = unit-disk position (x, z) in [-1, 1]
//    position.y  = vertex tag: 0.0 = bottom ring, 1.0 = top ring
//    normal      = surface normal in object space (computed CPU-side)
//    uv.x        = radial coord in [0, 1] (0 = on axis, 1 = at outer rim)
//    uv.y        = vertex role: 0=side, 1=top-cap, 2=bottom-cap
// ══════════════════════════════════════════════════════════════════════════

struct FluidVertex {
    float3 position;
    float3 normal;
    float2 uv;
};

struct FluidUniforms {
    float4x4 viewProjection;
    float3   cameraPos;
    float    time;
    float3   lightDir;
    float    floorY;
    float3   waterColor;     // base water tint
    float    waterAlpha;     // overall alpha multiplier
    float    radius;         // dish radius (wu)
    float    fluidHeight;    // height above floor (wu)
    float    waveAmp;        // surface displacement amplitude (wu)
    float    waveSpeed;      // angular frequency
};

struct FluidVertexOut {
    float4 clipPos [[position]];
    float3 worldPos;
    float3 worldNormal;
    float  surfaceFrac;      // 0 = bottom, 1 = top
    float  radialFrac;       // 0 = axis, 1 = rim
    float  time;
};

// ── Gerstner trochoidal ocean wave bank ─────────────────────────────
// Reference: Tessendorf 2001 "Simulating Ocean Water"; Mastin/Watterberg
// 1987 Gerstner waves. Each wave shifts surface vertices BOTH vertically
// (height) AND horizontally (toward crest), producing sharp crests and
// rounded troughs — the characteristic "choppy" look of real ocean and
// stirred lab medium. Six wave components at varying directions /
// wavelengths / phases give natural-looking interference.
//
// Steepness Q controls crest sharpness: 0 = pure sine (smooth);
// 1/(ω·A·N) = sharpest before self-intersection. We use moderate Q for
// a lively but stable surface.
struct GerstnerOut {
    float dx;        // horizontal displacement X
    float dz;        // horizontal displacement Z
    float dy;        // vertical displacement Y (above mean surface)
    float3 normal;   // analytic surface normal at displaced point
};

inline GerstnerOut gerstnerWaves(float x, float z, float t,
                                 float amp, float speed) {
    // Wave bank: (dirX, dirZ, wavelength, amplitudeMul, phaseSpeedMul, Q)
    // No address-space qualifier — Metal doesn't allow `constant` on
    // function-local arrays.
    float W[6][6] = {
        { 0.95,  0.31, 14.0, 1.00, 1.00, 0.55 },
        { 0.71,  0.71,  8.5, 0.65, 1.20, 0.45 },
        {-0.45,  0.89,  5.0, 0.40, 1.55, 0.35 },
        { 0.31, -0.95,  3.2, 0.30, 1.85, 0.28 },
        {-0.83, -0.55,  2.0, 0.22, 2.30, 0.22 },
        { 0.55, -0.83,  1.3, 0.16, 2.80, 0.18 }
    };
    GerstnerOut g;
    g.dx = 0.0; g.dz = 0.0; g.dy = 0.0;
    float3 nAcc = float3(0, 1, 0);

    for (int k = 0; k < 6; k++) {
        float dx_k = W[k][0];
        float dz_k = W[k][1];
        float wavelen = W[k][2];
        float A = amp * W[k][3];
        float phaseSpeed = speed * W[k][4];
        float Q = W[k][5];

        float omega = 2.0 * 3.14159265 / wavelen;        // angular wavenumber
        float c = sqrt(9.81 / omega);                    // deep-water phase velocity
        float theta = omega * (dx_k * x + dz_k * z) + c * phaseSpeed * t;
        float cosT = cos(theta);
        float sinT = sin(theta);

        g.dx += Q * A * dx_k * cosT;
        g.dz += Q * A * dz_k * cosT;
        g.dy += A * sinT;

        // Analytic partial derivatives for the normal (Tessendorf eq. 12)
        float WA  = omega * A;
        float qWA = Q * WA;
        nAcc.x -= dx_k * WA * cosT;
        nAcc.z -= dz_k * WA * cosT;
        nAcc.y -= qWA * sinT;
    }

    // Optional micro-chop FBM for fine detail (small amplitude, high freq).
    float chop = sin(0.55 * x + 0.71 * z + 1.7 * t) * 0.10
               + cos(0.83 * x - 0.61 * z + 1.3 * t) * 0.08
               + sin(0.42 * (x + z) + 2.4 * t) * 0.06;
    g.dy += chop * amp * 0.35;

    g.normal = normalize(nAcc);
    return g;
}

vertex FluidVertexOut fluidVertex(
    uint                       vid  [[vertex_id]],
    const device FluidVertex*  verts [[buffer(0)]],
    constant FluidUniforms&    u    [[buffer(1)]]
) {
    FluidVertex v = verts[vid];

    float surfaceTag = v.position.y;     // 0 = bottom, 1 = top
    float xUnit = v.position.x;
    float zUnit = v.position.z;

    // Base world position: scale unit-disk XZ by dish radius, lift Y.
    float3 worldPos;
    worldPos.x = xUnit * u.radius;
    worldPos.z = zUnit * u.radius;

    // Bottom sits 0.20 wu above the floor (so it doesn't z-fight the
    // substrate quad); top sits at floor + fluidHeight.
    float baseY = u.floorY + 0.20;
    float topY  = u.floorY + u.fluidHeight;
    worldPos.y = mix(baseY, topY, surfaceTag);

    // Apply Gerstner wave displacement (XZ + Y) to TOP cap vertices.
    // Rim damping keeps the surface attached to the cylinder wall and
    // prevents wave energy from "escaping" past the dish edge.
    float3 worldNormal = v.normal;
    if (surfaceTag > 0.5) {
        float radial = sqrt(xUnit * xUnit + zUnit * zUnit);
        float rimDamp = 1.0 - smoothstep(0.80, 1.00, radial);
        GerstnerOut g = gerstnerWaves(worldPos.x, worldPos.z, u.time,
                                      u.waveAmp, u.waveSpeed);
        worldPos.x += g.dx * rimDamp;
        worldPos.z += g.dz * rimDamp;
        worldPos.y += g.dy * rimDamp;
        // Use the analytic Gerstner normal on the top cap; on the sides
        // keep the CPU-supplied radial normal.
        if (v.uv.y > 0.5 && v.uv.y < 1.5) {
            worldNormal = mix(float3(0, 1, 0), g.normal, rimDamp);
        }
    }

    FluidVertexOut out;
    out.clipPos = u.viewProjection * float4(worldPos, 1.0);
    out.worldPos = worldPos;
    out.worldNormal = worldNormal;
    out.surfaceFrac = surfaceTag;
    out.radialFrac  = v.uv.x;
    out.time = u.time;
    return out;
}

fragment float4 fluidFragment(
    FluidVertexOut in [[stage_in]],
    constant FluidUniforms& u [[buffer(1)]]
) {
    float3 N = normalize(in.worldNormal);
    float3 V = normalize(u.cameraPos - in.worldPos);
    float3 L = normalize(u.lightDir);

    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);

    // Fresnel: stronger at grazing angles. We deliberately keep the
    // outward normals on the side-walls so when viewed from inside the
    // dish, dot(N,V) is small and the fresnel rim glows on the cylinder
    // walls. Top-cap viewed from above gets a bright specular sheen on
    // the wave crests.
    float fresnel = pow(1.0 - NdotV, 2.5);

    // Specular highlight on wave crests (top cap only, where N has a
    // significant horizontal component due to wave tilt).
    float waveTilt = clamp(1.0 - N.y, 0.0, 1.0);
    float spec = waveTilt * pow(NdotL, 8.0) * 0.6;

    // Soft sub-surface scattering — lighter where light penetrates.
    float sss = clamp(NdotL * 0.4 + 0.3, 0.0, 1.0);

    // Base water color modulated by depth (slightly darker at bottom,
    // lighter near the surface). depth = 1 - surfaceFrac.
    float depth = 1.0 - in.surfaceFrac;
    float3 baseCol = u.waterColor * (1.0 - 0.30 * depth);

    float3 finalColor = baseCol * sss
                      + baseCol * fresnel * 0.8
                      + float3(0.6, 0.85, 1.0) * spec;

    // Alpha: base waterAlpha gives the body its visible tint, fresnel
    // rim and wave-crest specular lift edges and crests further. Top
    // surface gets a small boost so the meniscus reads clearly.
    float topBoost = (in.surfaceFrac > 0.95) ? 0.10 : 0.0;
    float alpha = u.waterAlpha
                + fresnel * 0.18
                + spec * 0.30
                + topBoost;
    alpha = clamp(alpha, 0.0, 0.90);

    // Distance fog so far cylinder walls fade into background.
    float dist = length(u.cameraPos - in.worldPos);
    float fog = 1.0 - exp(-dist * 0.010);
    finalColor = mix(finalColor, float3(0.01, 0.02, 0.05), fog * 0.45);
    alpha *= (1.0 - fog * 0.20);

    return float4(finalColor, alpha);
}
