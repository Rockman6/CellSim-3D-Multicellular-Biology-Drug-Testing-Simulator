#include <metal_stdlib>
using namespace metal;

// ══════════════════════════════════════════════════════════════════════════
//  OrganelleRender.metal — Render actual GLB organelle meshes
//  Each organelle is drawn with a per-draw-call model matrix
//  (not instanced, since each organelle has a unique mesh)
//  Uses emissive bio-luminescent shading
// ══════════════════════════════════════════════════════════════════════════

struct Vertex {
    float3 position;
    float3 normal;
    float2 uv;
};

struct OrgUniforms {
    float4x4 viewProjection;
    float4x4 model;
    float3   cameraPos;
    float    time;
    float3   lightDir;
    float    pad0;
    float3   baseColor;
    float    emissiveIntensity;
    float3   emissiveColor;
    float    pad1;
    float3   cellCenter;     // World position of cell center
    float    cellRadius;     // Visual membrane radius
    float    furrowDepth;    // Cytokinesis pinch (must match CellRender)
    // Pierce-through glow pass (set non-zero only on selected-cell extra
    // pass; vanilla draw uses 0). At runtime: emissive *= (1 + glowBoost
    // * sinPulse), alpha lifted to ~0.85 so the organelle outline reads
    // through any other geometry.
    float    glowBoost;
    float    pad3;
    float    pad4;
};

struct OrgVertexOut {
    float4 clipPos   [[position]];
    float3 worldPos;
    float3 worldNormal;
    float3 baseColor;
    float3 emissiveColor;
    float  emissiveIntensity;
    float  time;
};

vertex OrgVertexOut glbOrganelleVertex(
    uint                vid  [[vertex_id]],
    const device Vertex* verts [[buffer(0)]],
    constant OrgUniforms& u    [[buffer(1)]]
) {
    Vertex v = verts[vid];

    float4 worldPos = u.model * float4(v.position, 1.0);
    float3 worldNormal = normalize((u.model * float4(v.normal, 0.0)).xyz);

    // ── MEMBRANE WALL: clamp vertex inside cell envelope ─────────
    // Follows the same furrow/elongation deformation as the membrane
    // shader so organelles never escape the pinched waist during
    // cytokinesis. Without this an organelle at (0, 0.85R, 0) would
    // pass the sphere-clamp but stick out through the pinched membrane
    // at x=0 where the real envelope contracts to (1-fw)×R.
    if (u.cellRadius > 0.01) {
        float3 toCell = worldPos.xyz - u.cellCenter;
        float fw = clamp(u.furrowDepth, 0.0, 1.0);
        // Compute local direction inside unit sphere.
        float R = u.cellRadius;
        float3 localP = toCell / R;
        // Match CellRender.metal: pinch along y/z at x≈0, elongate x.
        float pinch = 1.0 - fw * exp(-localP.x * localP.x * 6.0);
        pinch = max(pinch, 0.02);
        // If fw > 0, compute the max allowed radial extent in y/z
        // at this x, and clamp.
        float maxLateral = pinch * 0.90;     // 10% safety inside
        float maxAxial   = (1.0 + fw * 0.45) * 0.90;
        // Separate axial (x) and lateral (y,z) components.
        float latLen = length(float2(localP.y, localP.z));
        if (latLen > maxLateral && latLen > 0.001) {
            float k = maxLateral / latLen;
            localP.y *= k;
            localP.z *= k;
        }
        if (fabs(localP.x) > maxAxial) {
            localP.x = sign(localP.x) * maxAxial;
        }
        // Final spherical fallback so we never exceed 90% of R overall.
        float d = length(localP);
        if (d > 0.90 && d > 0.001) {
            localP *= 0.90 / d;
        }
        worldPos.xyz = u.cellCenter + localP * R;
    }

    OrgVertexOut out;
    out.clipPos          = u.viewProjection * worldPos;
    out.worldPos         = worldPos.xyz;
    out.worldNormal      = worldNormal;
    out.baseColor        = u.baseColor;
    out.emissiveColor    = u.emissiveColor;
    out.emissiveIntensity = u.emissiveIntensity;
    out.time             = u.time;
    return out;
}

fragment float4 glbOrganelleFragment(
    OrgVertexOut       in [[stage_in]],
    constant OrgUniforms& u  [[buffer(1)]]
) {
    float3 N = normalize(in.worldNormal);
    float3 V = normalize(u.cameraPos - in.worldPos);
    float3 L = normalize(u.lightDir);
    float3 H = normalize(L + V);

    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float NdotH = max(dot(N, H), 0.0);

    // Fresnel rim
    float fresnel = pow(1.0 - NdotV, 2.5);

    // Diffuse (soft Lambertian)
    float diffuse = NdotL * 0.4 + 0.15;

    // Specular (wet surface)
    float spec = pow(NdotH, 55.0) * 0.35;

    // Emissive glow (subsurface / bioluminescence) — gentle baseline
    // wash matched to the organelle's signature hue. The pierce-through
    // pulse (below) is what makes the selected cell's organelles light up.
    float3 emissive = in.emissiveColor * in.emissiveIntensity;

    // Sinusoidal pulse — slow rhythmic "breathing" of every organelle.
    float pulse = 0.9 + 0.1 * sin(in.time * 1.5 + in.worldPos.x * 2.0 + in.worldPos.z * 1.5);
    emissive *= pulse;

    // ── PIERCE-THROUGH GLOW PASS (only when glowBoost > 0) ────────────
    // The selected-cell pass sets glowBoost to ~3-6 so the organelle
    // outline becomes a strong pulsing emissive that punches through
    // anything in front (depth test disabled by the C++ side). Sinusoidal
    // pulse here is fast and bright — like a highlight halo.
    float glowPulse = 0.55 + 0.45 * sin(in.time * 3.5
                       + in.worldPos.x * 1.5 + in.worldPos.z * 1.2);
    if (u.glowBoost > 0.001) {
        // Use baseColor to drive the glow so each organelle keeps its
        // signature hue (nucleus blue, ER teal, Golgi yellow). Strong
        // edge contribution so the SHAPE (especially nucleus outline)
        // reads clearly.
        float3 glow = in.baseColor * (u.glowBoost * glowPulse)
                    + in.baseColor * fresnel * (u.glowBoost * 1.5);
        emissive += glow;
    }

    // Combine — original look (when glowBoost=0): base + diffuse +
    // emissive + spec + fresnel. Selected-cell pass adds the glow into
    // emissive above, dominating the colour.
    float3 finalColor = in.baseColor * (0.55 + diffuse * 0.40)
                      + emissive
                      + float3(0.3, 0.5, 0.8) * spec
                      + in.baseColor * fresnel * 0.55;

    // Very light fog — don't wash the color out to black.
    float dist = length(u.cameraPos - in.worldPos);
    float fog = 1.0 - exp(-dist * 0.008);
    finalColor = mix(finalColor, float3(0.02, 0.03, 0.08), fog * 0.4);

    // Vanilla pass: original transparent organelle. Selected-cell pass:
    // alpha lifted by glow so the outline pierces other geometry.
    float baseAlpha = 0.02 + fresnel * 0.18 + in.emissiveIntensity * 0.08;
    float glowAlpha = (u.glowBoost > 0.001)
        ? (0.45 + fresnel * 0.50 + glowPulse * 0.20)
        : 0.0;
    float alpha = clamp(max(baseAlpha, glowAlpha), 0.0, 0.95);
    alpha *= (1.0 - fog * 0.2);

    return float4(finalColor, alpha);
}
