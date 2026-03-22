#pragma once

#import <Metal/Metal.h>
#include <vector>
#include "../gpu/ShaderTypes.h"

// ══════════════════════════════════════════════════════════════════════════
//  MeshLibrary — Generates shared geometry (icospheres at multiple LODs)
// ══════════════════════════════════════════════════════════════════════════

struct MeshData {
    id<MTLBuffer> vertexBuffer;
    id<MTLBuffer> indexBuffer;
    uint32_t      indexCount;
};

class MeshLibrary {
public:
    bool init(id<MTLDevice> device);

    // LOD sphere meshes
    const MeshData& sphereLOD(int level) const; // 0=low, 1=med, 2=high

    // Substrate disc mesh
    const MeshData& substrate() const { return substrate_; }

private:
    MeshData spheres_[3]; // LOD 0, 1, 2
    MeshData substrate_;

    struct RawMesh {
        std::vector<Vertex>   vertices;
        std::vector<uint32_t> indices;
    };

    RawMesh generateIcosphere(int subdivisions);
    RawMesh generateDisc(float radius, int segments);
    MeshData upload(id<MTLDevice> device, const RawMesh& mesh);
};
