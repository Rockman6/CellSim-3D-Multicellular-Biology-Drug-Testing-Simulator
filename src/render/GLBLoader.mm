#define CGLTF_IMPLEMENTATION
#include <cgltf.h>

#import "GLBLoader.h"
#include <cstdio>
#include <cmath>
#include <cfloat>

// ══════════════════════════════════════════════════════════════════════════
//  GLBLoader — Parse GLB, apply node transforms (scale+rotation),
//  upload pre-transformed vertices to Metal
// ══════════════════════════════════════════════════════════════════════════

// Transform a point by a 4x4 column-major matrix
static simd_float3 transformPoint(const float m[16], simd_float3 p) {
    return simd_make_float3(
        m[0]*p.x + m[4]*p.y + m[8]*p.z  + m[12],
        m[1]*p.x + m[5]*p.y + m[9]*p.z  + m[13],
        m[2]*p.x + m[6]*p.y + m[10]*p.z + m[14]
    );
}

// Transform a normal (no translation, re-normalize)
static simd_float3 transformNormal(const float m[16], simd_float3 n) {
    simd_float3 r = simd_make_float3(
        m[0]*n.x + m[4]*n.y + m[8]*n.z,
        m[1]*n.x + m[5]*n.y + m[9]*n.z,
        m[2]*n.x + m[6]*n.y + m[10]*n.z
    );
    float len = simd_length(r);
    return len > 0.0001f ? r / len : simd_make_float3(0, 1, 0);
}

LoadedMesh GLBLoader::load(id<MTLDevice> device, const std::string& path, uint32_t maxTriangles) {
    LoadedMesh result;
    result.valid = false;

    cgltf_options options = {};
    cgltf_data* data = nullptr;

    if (cgltf_parse_file(&options, path.c_str(), &data) != cgltf_result_success) {
        fprintf(stderr, "[GLBLoader] Failed to parse: %s\n", path.c_str());
        return result;
    }
    if (cgltf_load_buffers(&options, data, path.c_str()) != cgltf_result_success) {
        fprintf(stderr, "[GLBLoader] Failed to load buffers: %s\n", path.c_str());
        cgltf_free(data);
        return result;
    }

    std::vector<Vertex>   allVertices;
    std::vector<uint32_t> allIndices;

    // Process each node that has a mesh — apply the node's world transform
    // This is what Three.js GLTFLoader does automatically
    for (cgltf_size ni = 0; ni < data->nodes_count; ni++) {
        cgltf_node* node = &data->nodes[ni];
        if (!node->mesh) continue;

        // Get node world transform (includes parent hierarchy + node scale/rotation/translation)
        float worldMat[16];
        cgltf_node_transform_world(node, worldMat);

        cgltf_mesh* mesh = node->mesh;
        for (cgltf_size pi = 0; pi < mesh->primitives_count; pi++) {
            cgltf_primitive* prim = &mesh->primitives[pi];
            if (prim->type != cgltf_primitive_type_triangles) continue;

            cgltf_accessor* posAcc = nullptr;
            cgltf_accessor* normAcc = nullptr;
            cgltf_accessor* uvAcc = nullptr;

            for (cgltf_size ai = 0; ai < prim->attributes_count; ai++) {
                if (prim->attributes[ai].type == cgltf_attribute_type_position)
                    posAcc = prim->attributes[ai].data;
                else if (prim->attributes[ai].type == cgltf_attribute_type_normal)
                    normAcc = prim->attributes[ai].data;
                else if (prim->attributes[ai].type == cgltf_attribute_type_texcoord)
                    uvAcc = prim->attributes[ai].data;
            }
            if (!posAcc) continue;

            uint32_t baseVertex = (uint32_t)allVertices.size();

            for (cgltf_size vi = 0; vi < posAcc->count; vi++) {
                Vertex v = {};
                float tmp3[3] = {0, 0, 0};
                float tmp2[2] = {0, 0};

                cgltf_accessor_read_float(posAcc, vi, tmp3, 3);
                simd_float3 rawPos = simd_make_float3(tmp3[0], tmp3[1], tmp3[2]);
                // Apply node world transform to position
                v.position = transformPoint(worldMat, rawPos);

                if (normAcc) {
                    cgltf_accessor_read_float(normAcc, vi, tmp3, 3);
                    simd_float3 rawNorm = simd_make_float3(tmp3[0], tmp3[1], tmp3[2]);
                    v.normal = transformNormal(worldMat, rawNorm);
                } else {
                    float len = simd_length(v.position);
                    v.normal = len > 0.0001f ? v.position / len : simd_make_float3(0, 1, 0);
                }

                if (uvAcc) {
                    cgltf_accessor_read_float(uvAcc, vi, tmp2, 2);
                    v.uv = simd_make_float2(tmp2[0], tmp2[1]);
                }

                allVertices.push_back(v);
            }

            if (prim->indices) {
                for (cgltf_size ii = 0; ii < prim->indices->count; ii++) {
                    uint32_t idx = (uint32_t)cgltf_accessor_read_index(prim->indices, ii);
                    allIndices.push_back(baseVertex + idx);
                }
            } else {
                for (cgltf_size vi = 0; vi < posAcc->count; vi++) {
                    allIndices.push_back(baseVertex + (uint32_t)vi);
                }
            }
        }
    }

    cgltf_free(data);

    if (allVertices.empty() || allIndices.empty()) {
        fprintf(stderr, "[GLBLoader] No geometry found in: %s\n", path.c_str());
        return result;
    }

    // Compute bounding box of the TRANSFORMED vertices
    float minX=FLT_MAX, minY=FLT_MAX, minZ=FLT_MAX;
    float maxX=-FLT_MAX, maxY=-FLT_MAX, maxZ=-FLT_MAX;
    for (const auto& v : allVertices) {
        if (v.position.x < minX) minX = v.position.x;
        if (v.position.x > maxX) maxX = v.position.x;
        if (v.position.y < minY) minY = v.position.y;
        if (v.position.y > maxY) maxY = v.position.y;
        if (v.position.z < minZ) minZ = v.position.z;
        if (v.position.z > maxZ) maxZ = v.position.z;
    }
    float sx = maxX-minX, sy = maxY-minY, sz = maxZ-minZ;
    result.maxDimension = fmaxf(sx, fmaxf(sy, sz));
    if (result.maxDimension < 0.0001f) result.maxDimension = 1.0f;
    result.center = simd_make_float3((minX+maxX)*0.5f, (minY+maxY)*0.5f, (minZ+maxZ)*0.5f);

    // Decimate if too many triangles
    uint32_t totalTris = (uint32_t)(allIndices.size() / 3);
    if (maxTriangles > 0 && totalTris > maxTriangles) {
        uint32_t step = totalTris / maxTriangles;
        if (step < 2) step = 2;
        std::vector<uint32_t> decimated;
        decimated.reserve(maxTriangles * 3);
        for (uint32_t t = 0; t < totalTris && decimated.size() < (size_t)maxTriangles * 3; t += step) {
            decimated.push_back(allIndices[t * 3]);
            decimated.push_back(allIndices[t * 3 + 1]);
            decimated.push_back(allIndices[t * 3 + 2]);
        }
        printf("[GLBLoader] Decimated %s: %u → %zu triangles (step=%u)\n",
               path.c_str(), totalTris, decimated.size() / 3, step);
        allIndices = std::move(decimated);
    }

    printf("[GLBLoader] %s: %zu verts, %zu idx, bbox=%.2f x %.2f x %.2f\n",
           path.c_str(), allVertices.size(), allIndices.size(), sx, sy, sz);

    // Upload to Metal
    result.vertexBuffer = [device newBufferWithBytes:allVertices.data()
                                             length:allVertices.size() * sizeof(Vertex)
                                            options:MTLResourceStorageModeShared];
    result.indexBuffer = [device newBufferWithBytes:allIndices.data()
                                            length:allIndices.size() * sizeof(uint32_t)
                                           options:MTLResourceStorageModeShared];
    result.indexCount = (uint32_t)allIndices.size();
    result.valid = true;

    return result;
}
