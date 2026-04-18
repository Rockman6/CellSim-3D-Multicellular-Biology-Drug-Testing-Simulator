#pragma once

#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#include <vector>

// ══════════════════════════════════════════════════════════════════════════
//  MetalContext — Manages MTLDevice, command queue, render pipeline states
// ══════════════════════════════════════════════════════════════════════════

class MetalContext {
public:
    bool init(CAMetalLayer* layer);
    void shutdown();

    // Getters
    id<MTLDevice>              device()       const { return device_; }
    id<MTLCommandQueue>        commandQueue() const { return commandQueue_; }
    id<MTLRenderPipelineState> cellPipeline() const { return cellPipeline_; }
    id<MTLRenderPipelineState> wirePipeline() const { return wirePipeline_; }
    id<MTLRenderPipelineState> organellePipeline() const { return organellePipeline_; }
    id<MTLRenderPipelineState> glbOrganellePipeline() const { return glbOrganellePipeline_; }
    id<MTLRenderPipelineState> substratePipeline() const { return substratePipeline_; }
    id<MTLRenderPipelineState> moleculeAtomPipeline() const { return moleculeAtomPipeline_; }
    id<MTLRenderPipelineState> moleculeBondPipeline() const { return moleculeBondPipeline_; }
    id<MTLRenderPipelineState> fluidPipeline()        const { return fluidPipeline_; }
    id<MTLDepthStencilState>   depthState()   const { return depthState_; }
    id<MTLDepthStencilState>   depthStateNoWrite() const { return depthStateNoWrite_; }
    id<MTLDepthStencilState>   depthStateAlways()  const { return depthStateAlways_; }
    CAMetalLayer*              metalLayer()   const { return metalLayer_; }
    id<MTLTexture>             depthTexture() const { return depthTexture_; }

    void recreateDepthTexture(uint32_t width, uint32_t height);

private:
    bool createCellPipeline();
    bool createWirePipeline();
    bool createOrganellePipeline();
    bool createGLBOrganellePipeline();
    bool createSubstratePipeline();
    bool createMoleculeAtomPipeline();
    bool createMoleculeBondPipeline();
    bool createFluidPipeline();

    id<MTLDevice>              device_       = nil;
    id<MTLCommandQueue>        commandQueue_ = nil;
    id<MTLLibrary>             shaderLibrary_ = nil;
    id<MTLRenderPipelineState> cellPipeline_ = nil;
    id<MTLRenderPipelineState> wirePipeline_ = nil;
    id<MTLRenderPipelineState> organellePipeline_ = nil;
    id<MTLRenderPipelineState> glbOrganellePipeline_ = nil;
    id<MTLRenderPipelineState> substratePipeline_ = nil;
    id<MTLRenderPipelineState> moleculeAtomPipeline_ = nil;
    id<MTLRenderPipelineState> moleculeBondPipeline_ = nil;
    id<MTLRenderPipelineState> fluidPipeline_        = nil;
    id<MTLDepthStencilState>   depthState_   = nil;
    id<MTLDepthStencilState>   depthStateNoWrite_ = nil;
    id<MTLDepthStencilState>   depthStateAlways_  = nil;
    CAMetalLayer*              metalLayer_   = nil;
    id<MTLTexture>             depthTexture_ = nil;

    // Runtime-compiled shader libraries (one per .metal file)
    std::vector<id<MTLLibrary>> runtimeLibraries_;

    // Find a function across all libraries
    id<MTLFunction> findFunction(NSString* name);
};
