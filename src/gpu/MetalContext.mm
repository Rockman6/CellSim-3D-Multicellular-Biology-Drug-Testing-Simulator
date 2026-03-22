#import "MetalContext.h"
#import <Foundation/Foundation.h>

// ══════════════════════════════════════════════════════════════════════════
//  MetalContext — Device, queue, pipeline states
// ══════════════════════════════════════════════════════════════════════════

static MTLRenderPipelineDescriptor* makeBlendedPipelineDesc(
    id<MTLFunction> vert, id<MTLFunction> frag, bool depthWrite
) {
    MTLRenderPipelineDescriptor* desc = [[MTLRenderPipelineDescriptor alloc] init];
    desc.vertexFunction   = vert;
    desc.fragmentFunction = frag;
    desc.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
    desc.colorAttachments[0].blendingEnabled = YES;
    desc.colorAttachments[0].sourceRGBBlendFactor      = MTLBlendFactorSourceAlpha;
    desc.colorAttachments[0].destinationRGBBlendFactor  = MTLBlendFactorOneMinusSourceAlpha;
    desc.colorAttachments[0].sourceAlphaBlendFactor     = MTLBlendFactorOne;
    desc.colorAttachments[0].destinationAlphaBlendFactor= MTLBlendFactorOneMinusSourceAlpha;
    desc.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float;
    return desc;
}

static MTLRenderPipelineDescriptor* makeAdditivePipelineDesc(
    id<MTLFunction> vert, id<MTLFunction> frag
) {
    MTLRenderPipelineDescriptor* desc = [[MTLRenderPipelineDescriptor alloc] init];
    desc.vertexFunction   = vert;
    desc.fragmentFunction = frag;
    desc.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
    desc.colorAttachments[0].blendingEnabled = YES;
    desc.colorAttachments[0].sourceRGBBlendFactor      = MTLBlendFactorSourceAlpha;
    desc.colorAttachments[0].destinationRGBBlendFactor  = MTLBlendFactorOne; // additive
    desc.colorAttachments[0].sourceAlphaBlendFactor     = MTLBlendFactorOne;
    desc.colorAttachments[0].destinationAlphaBlendFactor= MTLBlendFactorOne;
    desc.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float;
    return desc;
}

bool MetalContext::init(CAMetalLayer* layer) {
    metalLayer_ = layer;
    device_ = MTLCreateSystemDefaultDevice();
    if (!device_) { NSLog(@"[Metal] No device"); return false; }

    commandQueue_ = [device_ newCommandQueue];
    if (!commandQueue_) { NSLog(@"[Metal] No queue"); return false; }

    metalLayer_.device = device_;
    metalLayer_.pixelFormat = MTLPixelFormatBGRA8Unorm;
    metalLayer_.framebufferOnly = YES;

    // Load shader library
    NSError* error = nil;
    NSString* execPath = [[NSBundle mainBundle] executablePath];
    NSString* execDir  = [execPath stringByDeletingLastPathComponent];
    NSString* libPath  = [execDir stringByAppendingPathComponent:@"default.metallib"];

    if ([[NSFileManager defaultManager] fileExistsAtPath:libPath]) {
        shaderLibrary_ = [device_ newLibraryWithURL:[NSURL fileURLWithPath:libPath] error:&error];
    }
    if (!shaderLibrary_) shaderLibrary_ = [device_ newDefaultLibrary];
    if (!shaderLibrary_) { NSLog(@"[Metal] No shaders: %@", error); return false; }

    // Depth states
    {
        MTLDepthStencilDescriptor* d = [[MTLDepthStencilDescriptor alloc] init];
        d.depthCompareFunction = MTLCompareFunctionLess;
        d.depthWriteEnabled = YES;
        depthState_ = [device_ newDepthStencilStateWithDescriptor:d];
    }
    {
        MTLDepthStencilDescriptor* d = [[MTLDepthStencilDescriptor alloc] init];
        d.depthCompareFunction = MTLCompareFunctionLess;
        d.depthWriteEnabled = NO; // for transparent objects
        depthStateNoWrite_ = [device_ newDepthStencilStateWithDescriptor:d];
    }

    if (!createSubstratePipeline()) return false;
    if (!createOrganellePipeline()) return false;
    if (!createGLBOrganellePipeline()) return false;
    if (!createCellPipeline()) return false;
    if (!createWirePipeline()) return false;

    NSLog(@"[Metal] Ready — %@", device_.name);
    return true;
}

bool MetalContext::createCellPipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        [shaderLibrary_ newFunctionWithName:@"cellVertex"],
        [shaderLibrary_ newFunctionWithName:@"cellFragment"], false
    );
    cellPipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!cellPipeline_) { NSLog(@"[Metal] Cell pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createWirePipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        [shaderLibrary_ newFunctionWithName:@"wireVertex"],
        [shaderLibrary_ newFunctionWithName:@"wireFragment"], false
    );
    wirePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!wirePipeline_) { NSLog(@"[Metal] Wire pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createOrganellePipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        [shaderLibrary_ newFunctionWithName:@"organelleVertex"],
        [shaderLibrary_ newFunctionWithName:@"organelleFragment"], false
    );
    organellePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!organellePipeline_) { NSLog(@"[Metal] Organelle pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createGLBOrganellePipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        [shaderLibrary_ newFunctionWithName:@"glbOrganelleVertex"],
        [shaderLibrary_ newFunctionWithName:@"glbOrganelleFragment"], false
    );
    glbOrganellePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!glbOrganellePipeline_) { NSLog(@"[Metal] GLB Organelle pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createSubstratePipeline() {
    NSError* error = nil;
    MTLRenderPipelineDescriptor* desc = [[MTLRenderPipelineDescriptor alloc] init];
    desc.vertexFunction   = [shaderLibrary_ newFunctionWithName:@"substrateVertex"];
    desc.fragmentFunction = [shaderLibrary_ newFunctionWithName:@"substrateFragment"];
    desc.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
    desc.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float;
    substratePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!substratePipeline_) { NSLog(@"[Metal] Substrate pipe: %@", error); return false; }
    return true;
}

void MetalContext::recreateDepthTexture(uint32_t width, uint32_t height) {
    MTLTextureDescriptor* desc = [MTLTextureDescriptor
        texture2DDescriptorWithPixelFormat:MTLPixelFormatDepth32Float
        width:width height:height mipmapped:NO];
    desc.storageMode = MTLStorageModePrivate;
    desc.usage = MTLTextureUsageRenderTarget;
    depthTexture_ = [device_ newTextureWithDescriptor:desc];
}

void MetalContext::shutdown() {
    cellPipeline_ = nil; wirePipeline_ = nil;
    organellePipeline_ = nil; glbOrganellePipeline_ = nil; substratePipeline_ = nil;
    depthState_ = nil; depthStateNoWrite_ = nil;
    depthTexture_ = nil; shaderLibrary_ = nil;
    commandQueue_ = nil; device_ = nil;
}
