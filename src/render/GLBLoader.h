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

        // All organelle positions in unit sphere coordinates (sphere of
        // radius 1 = cell membrane; then scaled by cellR world units).
        // World unit anchor: 1 wu = 5 µm (MICROMETERS_PER_WORLD_UNIT).
        // Target radii in µm / 5 µm/wu = wu; GLB scale = target_wu / (bbox_max/2).
        //
        // Real-world organelle sizes (Alberts MBoC Ch 12-13; BioNumbers):
        //   Nucleus radius    ≈ 5 µm   → 0.50 × cellR
        //   Mitochondrion len ≈ 2 µm   → 0.20 × cellR
        //   Rough ER sheet    ≈ 4–5 µm → 0.45 × cellR
        //   Smooth ER patch   ≈ 2 µm   → 0.20 × cellR
        //   Golgi stack       ≈ 5 µm   → 0.50 × cellR
        //
        // Scale formula: scale = (target_wu × 2) / bbox_max (reported in the
        // GLBLoader decimation log at startup).

        // Original muted base + emissive palette restored. The ADDITIVE
        // pierce-through glow pass for the selected cell is rendered
        // separately on top, so the static look stays calm and the glow
        // is a deliberate "selected" treatment.
        // nucleus.glb (bbox_max=9.40): 0.50 × 2 / 9.40 = 0.106
        c.nucleus = {
            "nucleus.glb",
            {0, 0, 0}, {0, 0, 0}, 0.106f,
            {0.267f, 0.4f, 0.8f},
            {0.0f, 0.067f, 0.2f},
            0.35f
        };

        // smooth ER.glb (bbox_max=1.20)
        c.smoothER = {
            "smooth ER.glb",
            {0.3f, 0, 0.1f}, {0, 0, 0}, 0.35f,
            {0.133f, 0.867f, 0.667f},
            {0.0f, 0.2f, 0.133f},
            0.12f
        };

        // rough ER.glb (bbox_max=3.91)
        c.roughER = {
            "rough ER.glb",
            {0, 0, 0.25f}, {0, 9.3f, 0}, 0.18f,
            {0.0f, 0.867f, 0.733f},
            {0.0f, 0.267f, 0.133f},
            0.10f
        };

        // golgi apparatus.glb (bbox_max=2.41)
        c.golgi = {
            "golgi apparatus.glb",
            {0, 0, -0.3f}, {0, 0.8f, 0}, 0.21f,
            {1.0f, 0.8f, 0.267f},
            {0.267f, 0.2f, 0.0f},
            0.10f
        };

        // mitochondria.glb (bbox_max=1.87): target 1.2 µm length = 0.24 wu.
        // Coral-pink emissive boosted so each mito glows like a power plant.
        c.mito[0] = {
            "mitochondria.glb",
            {0, 0.2f, 0.3f}, {0, 0, 0}, 0.13f,
            {1.0f, 0.45f, 0.55f},
            {1.00f, 0.40f, 0.30f},   // strong coral glow (was 0.60,0.20,0.10)
            0.70f                     // (was 0.15)
        };
        c.mito[1] = {
            "mitochondria.glb",
            {0, -0.2f, 0}, {0, 1.2f, 0}, 0.12f,
            {1.0f, 0.45f, 0.55f},
            {1.00f, 0.40f, 0.30f},
            0.70f
        };
        c.mito[2] = {
            "mitochondria.glb",
            {-0.25f, 0, -0.25f}, {0, -2.0f, 0.3f}, 0.14f,
            {1.0f, 0.45f, 0.55f},
            {1.00f, 0.40f, 0.30f},
            0.70f
        };

        return c;
    }
};
