#pragma once

#import <Metal/Metal.h>
#include "../gpu/ShaderTypes.h"
#include <vector>
#include <string>
#include <simd/simd.h>

// ══════════════════════════════════════════════════════════════════════════
//  GLBLoader — Parse GLB/glTF binary files into Metal vertex/index buffers
//  Uses cgltf (single-header C library)
// ══════════════════════════════════════════════════════════════════════════

struct LoadedMesh {
    id<MTLBuffer> vertexBuffer;
    id<MTLBuffer> indexBuffer;
    uint32_t      indexCount;
    bool          valid = false;
    // Bounding box info for normalization
    float         maxDimension = 1.0f;   // max of width/height/depth
    simd_float3   center = {0, 0, 0};    // bounding box center
};

class GLBLoader {
public:
    // Load a GLB file, extract all meshes, apply node transforms, upload to GPU
    // maxTriangles: if > 0, decimate meshes with more triangles than this
    static LoadedMesh load(id<MTLDevice> device, const std::string& path, uint32_t maxTriangles = 0);
};

// ══════════════════════════════════════════════════════════════════════════
//  OrganelleConfig — organelle placement, rotation, scale, color
// ══════════════════════════════════════════════════════════════════════════

struct OrganelleConfig {
    std::string filename;
    simd_float3 position;
    simd_float3 rotation; // Euler XYZ radians
    float       scale;
    // Stylize colors (from GLB model styling)
    simd_float3 color;    // base color
    simd_float3 emissive; // emissive color
    float       emissiveIntensity;
};

// All organelle configs with correct transforms
struct OrganelleConfigs {
    OrganelleConfig nucleus;
    OrganelleConfig smoothER;
    OrganelleConfig roughER;
    OrganelleConfig golgi;
    // Mitochondria: 3 instances with different transforms
    OrganelleConfig mito[3];

    static OrganelleConfigs defaults() {
        OrganelleConfigs c;

        // nucleus.glb — center of cell
        c.nucleus = {
            "nucleus.glb",
            {0, 0, 0}, {0, 0, 0}, 0.15f,
            {0.267f, 0.4f, 0.8f},    // 0x4466cc
            {0.0f, 0.067f, 0.2f},    // 0x001133
            0.35f
        };

        // smooth ER.glb
        c.smoothER = {
            "smooth ER.glb",
            {1.0f, 0, 0.3f}, {0, 0, 0}, 0.8f,
            {0.133f, 0.867f, 0.667f}, // 0x22ddaa
            {0.0f, 0.2f, 0.133f},    // 0x003322
            0.12f
        };

        // rough ER.glb
        c.roughER = {
            "rough ER.glb",
            {0, 0, 1.0f}, {0, 9.3f, 0}, 0.5f,
            {0.0f, 0.867f, 0.733f},   // 0x00ddbb
            {0.0f, 0.267f, 0.133f},   // 0x004422
            0.10f
        };

        // golgi apparatus.glb
        c.golgi = {
            "golgi apparatus.glb",
            {0, 0, -1.5f}, {0, 0.8f, 0}, 0.65f,
            {1.0f, 0.8f, 0.267f},     // 0xffcc44
            {0.267f, 0.2f, 0.0f},     // 0x443300
            0.10f
        };

        // mitochondria.glb — 3 instances
        c.mito[0] = {
            "mitochondria.glb",
            {0, 1.0f, 1.5f}, {0, 0, 0}, 0.6f,
            {1.0f, 0.4f, 0.533f},     // 0xff6688
            {0.2f, 0.0f, 0.067f},     // 0x330011
            0.08f
        };
        c.mito[1] = {
            "mitochondria.glb",
            {0, -1.0f, 0}, {0, 1.2f, 0}, 0.56f,
            {1.0f, 0.4f, 0.533f},
            {0.2f, 0.0f, 0.067f},
            0.08f
        };
        c.mito[2] = {
            "mitochondria.glb",
            {-1.0f, 0, -1.0f}, {0, -2.0f, 0.3f}, 0.72f,
            {1.0f, 0.4f, 0.533f},
            {0.2f, 0.0f, 0.067f},
            0.08f
        };

        return c;
    }
};
