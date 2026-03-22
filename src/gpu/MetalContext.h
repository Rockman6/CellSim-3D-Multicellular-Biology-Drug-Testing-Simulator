#pragma once

#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>

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
    id<MTLDepthStencilState>   depthState()   const { return depthState_; }
    id<MTLDepthStencilState>   depthStateNoWrite() const { return depthStateNoWrite_; }
    CAMetalLayer*              metalLayer()   const { return metalLayer_; }
    id<MTLTexture>             depthTexture() const { return depthTexture_; }

    void recreateDepthTexture(uint32_t width, uint32_t height);

private:
    bool createCellPipeline();
    bool createWirePipeline();
    bool createOrganellePipeline();
    bool createGLBOrganellePipeline();
    bool createSubstratePipeline();

    id<MTLDevice>              device_       = nil;
    id<MTLCommandQueue>        commandQueue_ = nil;
    id<MTLLibrary>             shaderLibrary_ = nil;
    id<MTLRenderPipelineState> cellPipeline_ = nil;
    id<MTLRenderPipelineState> wirePipeline_ = nil;
    id<MTLRenderPipelineState> organellePipeline_ = nil;
    id<MTLRenderPipelineState> glbOrganellePipeline_ = nil;
    id<MTLRenderPipelineState> substratePipeline_ = nil;
    id<MTLDepthStencilState>   depthState_   = nil;
    id<MTLDepthStencilState>   depthStateNoWrite_ = nil;
    CAMetalLayer*              metalLayer_   = nil;
    id<MTLTexture>             depthTexture_ = nil;
};
