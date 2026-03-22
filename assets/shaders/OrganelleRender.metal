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
    float4x4 model;          // per-organelle transform
    float3   cameraPos;
    float    time;
    float3   lightDir;
    float    pad0;
    float3   baseColor;      // organelle base color
    float    emissiveIntensity;
    float3   emissiveColor;
    float    pad1;
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

    // Emissive glow (subsurface / bioluminescence)
    float3 emissive = in.emissiveColor * in.emissiveIntensity;

    // Pulse (organelle activity)
    float pulse = 0.9 + 0.1 * sin(in.time * 1.5 + in.worldPos.x * 2.0 + in.worldPos.z * 1.5);
    emissive *= pulse;

    // Combine
    float3 finalColor = in.baseColor * diffuse
                      + emissive
                      + float3(0.3, 0.5, 0.8) * spec
                      + in.baseColor * fresnel * 0.3;

    // Fog
    float dist = length(u.cameraPos - in.worldPos);
    float fog = 1.0 - exp(-dist * 0.012);
    finalColor = mix(finalColor, float3(0.004, 0.008, 0.031), fog);

    // Translucent
    float alpha = 0.55 + fresnel * 0.25 + in.emissiveIntensity * 0.15;
    alpha = clamp(alpha, 0.0, 0.9);
    alpha *= (1.0 - fog * 0.3);

    return float4(finalColor, alpha);
}
