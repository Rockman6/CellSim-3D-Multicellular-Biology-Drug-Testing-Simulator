#import "MeshLibrary.h"
#include <cmath>
#include <map>

// ══════════════════════════════════════════════════════════════════════════
//  MeshLibrary — Icosphere generation + GPU buffer upload
// ══════════════════════════════════════════════════════════════════════════

bool MeshLibrary::init(id<MTLDevice> device) {
    // Generate spheres at 3 LOD levels
    // LOD 0: 1 subdivision (42 vertices)  — far
    // LOD 1: 2 subdivisions (162 vertices) — medium
    // LOD 2: 3 subdivisions (642 vertices) — near
    for (int i = 0; i < 3; i++) {
        auto mesh = generateIcosphere(i + 1);
        spheres_[i] = upload(device, mesh);
    }

    // Substrate disc (r=55, matching SCENE_BOUND)
    auto disc = generateDisc(55.0f, 64);
    substrate_ = upload(device, disc);

    return true;
}

const MeshData& MeshLibrary::sphereLOD(int level) const {
    return spheres_[std::max(0, std::min(2, level))];
}

MeshData MeshLibrary::upload(id<MTLDevice> device, const RawMesh& mesh) {
    MeshData data;
    data.vertexBuffer = [device newBufferWithBytes:mesh.vertices.data()
                                           length:mesh.vertices.size() * sizeof(Vertex)
                                          options:MTLResourceStorageModeShared];
    data.indexBuffer = [device newBufferWithBytes:mesh.indices.data()
                                          length:mesh.indices.size() * sizeof(uint32_t)
                                         options:MTLResourceStorageModeShared];
    data.indexCount = (uint32_t)mesh.indices.size();
    return data;
}

// ── Icosphere generation ────────────────────────────────────────────────
// Start with icosahedron (12 vertices, 20 faces), subdivide edges

MeshLibrary::RawMesh MeshLibrary::generateIcosphere(int subdivisions) {
    RawMesh mesh;

    // Golden ratio
    const float t = (1.0f + sqrtf(5.0f)) / 2.0f;

    // Icosahedron base vertices (normalized to unit sphere)
    std::vector<simd_float3> positions = {
        simd_normalize(simd_make_float3(-1,  t,  0)),
        simd_normalize(simd_make_float3( 1,  t,  0)),
        simd_normalize(simd_make_float3(-1, -t,  0)),
        simd_normalize(simd_make_float3( 1, -t,  0)),
        simd_normalize(simd_make_float3( 0, -1,  t)),
        simd_normalize(simd_make_float3( 0,  1,  t)),
        simd_normalize(simd_make_float3( 0, -1, -t)),
        simd_normalize(simd_make_float3( 0,  1, -t)),
        simd_normalize(simd_make_float3( t,  0, -1)),
        simd_normalize(simd_make_float3( t,  0,  1)),
        simd_normalize(simd_make_float3(-t,  0, -1)),
        simd_normalize(simd_make_float3(-t,  0,  1)),
    };

    // Icosahedron faces (20 triangles)
    std::vector<uint32_t> indices = {
        0,11,5,  0,5,1,   0,1,7,  0,7,10, 0,10,11,
        1,5,9,   5,11,4,  11,10,2, 10,7,6, 7,1,8,
        3,9,4,   3,4,2,   3,2,6,  3,6,8,  3,8,9,
        4,9,5,   2,4,11,  6,2,10, 8,6,7,  9,8,1,
    };

    // Subdivide
    using EdgeKey = uint64_t;
    auto makeEdgeKey = [](uint32_t a, uint32_t b) -> EdgeKey {
        if (a > b) std::swap(a, b);
        return ((uint64_t)a << 32) | (uint64_t)b;
    };

    for (int s = 0; s < subdivisions; s++) {
        std::map<EdgeKey, uint32_t> midpointCache;
        std::vector<uint32_t> newIndices;

        auto getMidpoint = [&](uint32_t i0, uint32_t i1) -> uint32_t {
            EdgeKey key = makeEdgeKey(i0, i1);
            auto it = midpointCache.find(key);
            if (it != midpointCache.end()) return it->second;

            simd_float3 mid = simd_normalize((positions[i0] + positions[i1]) * 0.5f);
            uint32_t idx = (uint32_t)positions.size();
            positions.push_back(mid);
            midpointCache[key] = idx;
            return idx;
        };

        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t v0 = indices[i], v1 = indices[i+1], v2 = indices[i+2];
            uint32_t a = getMidpoint(v0, v1);
            uint32_t b = getMidpoint(v1, v2);
            uint32_t c = getMidpoint(v2, v0);

            newIndices.insert(newIndices.end(), {v0, a, c});
            newIndices.insert(newIndices.end(), {v1, b, a});
            newIndices.insert(newIndices.end(), {v2, c, b});
            newIndices.insert(newIndices.end(), {a, b, c});
        }
        indices = newIndices;
    }

    // Build Vertex array (position = normal for unit sphere)
    mesh.vertices.reserve(positions.size());
    for (const auto& p : positions) {
        Vertex v;
        v.position = p;
        v.normal = p; // Unit sphere: normal = position
        // Simple spherical UV
        v.uv.x = 0.5f + atan2f(p.z, p.x) / (2.0f * (float)M_PI);
        v.uv.y = 0.5f - asinf(p.y) / (float)M_PI;
        mesh.vertices.push_back(v);
    }
    mesh.indices = indices;

    return mesh;
}

// ── Disc mesh (substrate floor) ─────────────────────────────────────────

MeshLibrary::RawMesh MeshLibrary::generateDisc(float radius, int segments) {
    RawMesh mesh;

    // Center vertex
    Vertex center;
    center.position = {0, 0, 0};
    center.normal = {0, 1, 0};
    center.uv = {0.5f, 0.5f};
    mesh.vertices.push_back(center);

    // Rim vertices
    for (int i = 0; i <= segments; i++) {
        float angle = (float)i / (float)segments * 2.0f * M_PI;
        Vertex v;
        v.position = {radius * cosf(angle), 0, radius * sinf(angle)};
        v.normal = {0, 1, 0};
        v.uv = {0.5f + 0.5f * cosf(angle), 0.5f + 0.5f * sinf(angle)};
        mesh.vertices.push_back(v);
    }

    // Triangle fan
    for (int i = 1; i <= segments; i++) {
        mesh.indices.push_back(0);
        mesh.indices.push_back(i);
        mesh.indices.push_back(i + 1);
    }

    return mesh;
}
