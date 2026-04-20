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

    // Try precompiled metallib first
    if ([[NSFileManager defaultManager] fileExistsAtPath:libPath]) {
        NSDictionary* attrs = [[NSFileManager defaultManager] attributesOfItemAtPath:libPath error:nil];
        if ([attrs fileSize] > 0) {
            shaderLibrary_ = [device_ newLibraryWithURL:[NSURL fileURLWithPath:libPath] error:&error];
        }
    }
    if (!shaderLibrary_) shaderLibrary_ = [device_ newDefaultLibrary];

    // Fallback: compile shaders from source at runtime (each file separately)
    if (!shaderLibrary_) {
        NSLog(@"[Metal] No precompiled metallib — compiling shaders from source...");
        NSString* shaderDir = [[execDir stringByDeletingLastPathComponent]
                               stringByAppendingPathComponent:@"assets/shaders"];
        NSArray* shaderFiles = [[NSFileManager defaultManager]
                                contentsOfDirectoryAtPath:shaderDir error:nil];
        MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
        for (NSString* file in shaderFiles) {
            if (![file hasSuffix:@".metal"]) continue;
            NSString* path = [shaderDir stringByAppendingPathComponent:file];
            NSString* src = [NSString stringWithContentsOfFile:path
                                                     encoding:NSUTF8StringEncoding error:nil];
            if (!src) continue;
            NSError* compErr = nil;
            id<MTLLibrary> lib = [device_ newLibraryWithSource:src options:opts error:&compErr];
            if (lib) {
                runtimeLibraries_.push_back(lib);
                NSLog(@"[Metal] Compiled shader: %@", file);
            } else {
                NSLog(@"[Metal] Failed to compile %@: %@", file, compErr);
            }
        }
        if (!runtimeLibraries_.empty()) {
            shaderLibrary_ = runtimeLibraries_[0];
        }
    }
    if (!shaderLibrary_ && runtimeLibraries_.empty()) {
        NSLog(@"[Metal] No shaders: %@", error); return false;
    }

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
        d.depthWriteEnabled = NO;
        depthStateNoWrite_ = [device_ newDepthStencilStateWithDescriptor:d];
    }
    {
        MTLDepthStencilDescriptor* d = [[MTLDepthStencilDescriptor alloc] init];
        d.depthCompareFunction = MTLCompareFunctionAlways;  // pierce-through
        d.depthWriteEnabled = NO;
        depthStateAlways_ = [device_ newDepthStencilStateWithDescriptor:d];
    }

    if (!createSubstratePipeline()) return false;
    if (!createOrganellePipeline()) return false;
    if (!createGLBOrganellePipeline()) return false;
    if (!createCellPipeline()) return false;
    if (!createWirePipeline()) return false;
    createMoleculeAtomPipeline();
    createMoleculeBondPipeline();
    createFluidPipeline();

    NSLog(@"[Metal] Ready — %@", device_.name);
    return true;
}

id<MTLFunction> MetalContext::findFunction(NSString* name) {
    // Try primary library first
    id<MTLFunction> fn = [shaderLibrary_ newFunctionWithName:name];
    if (fn) return fn;
    // Search runtime libraries
    for (auto& lib : runtimeLibraries_) {
        fn = [lib newFunctionWithName:name];
        if (fn) return fn;
    }
    return nil;
}

bool MetalContext::createCellPipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        findFunction(@"cellVertex"),
        findFunction(@"cellFragment"), false
    );
    cellPipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!cellPipeline_) { NSLog(@"[Metal] Cell pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createWirePipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        findFunction(@"wireVertex"),
        findFunction(@"wireFragment"), false
    );
    wirePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!wirePipeline_) { NSLog(@"[Metal] Wire pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createOrganellePipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        findFunction(@"organelleVertex"),
        findFunction(@"organelleFragment"), false
    );
    organellePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!organellePipeline_) { NSLog(@"[Metal] Organelle pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createGLBOrganellePipeline() {
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(
        findFunction(@"glbOrganelleVertex"),
        findFunction(@"glbOrganelleFragment"), false
    );
    glbOrganellePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!glbOrganellePipeline_) { NSLog(@"[Metal] GLB Organelle pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createSubstratePipeline() {
    NSError* error = nil;
    MTLRenderPipelineDescriptor* desc = [[MTLRenderPipelineDescriptor alloc] init];
    desc.vertexFunction   = findFunction(@"substrateVertex");
    desc.fragmentFunction = findFunction(@"substrateFragment");
    desc.colorAttachments[0].pixelFormat = MTLPixelFormatBGRA8Unorm;
    desc.depthAttachmentPixelFormat = MTLPixelFormatDepth32Float;
    substratePipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!substratePipeline_) { NSLog(@"[Metal] Substrate pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createMoleculeAtomPipeline() {
    id<MTLFunction> vert = findFunction(@"moleculeAtomVertex");
    id<MTLFunction> frag = findFunction(@"moleculeFragment");
    if (!vert || !frag) {
        NSLog(@"[Metal] Molecule atom shader functions not found (optional)");
        return false;
    }
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(vert, frag, false);
    moleculeAtomPipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!moleculeAtomPipeline_) { NSLog(@"[Metal] Molecule atom pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createMoleculeBondPipeline() {
    id<MTLFunction> vert = findFunction(@"moleculeBondVertex");
    id<MTLFunction> frag = findFunction(@"moleculeFragment");
    if (!vert || !frag) {
        NSLog(@"[Metal] Molecule bond shader functions not found (optional)");
        return false;
    }
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(vert, frag, false);
    moleculeBondPipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!moleculeBondPipeline_) { NSLog(@"[Metal] Molecule bond pipe: %@", error); return false; }
    return true;
}

bool MetalContext::createFluidPipeline() {
    id<MTLFunction> vert = findFunction(@"fluidVertex");
    id<MTLFunction> frag = findFunction(@"fluidFragment");
    if (!vert || !frag) {
        NSLog(@"[Metal] Fluid shader functions not found (optional)");
        return false;
    }
    NSError* error = nil;
    auto desc = makeBlendedPipelineDesc(vert, frag, false);
    fluidPipeline_ = [device_ newRenderPipelineStateWithDescriptor:desc error:&error];
    if (!fluidPipeline_) { NSLog(@"[Metal] Fluid pipe: %@", error); return false; }
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
    moleculeAtomPipeline_ = nil; moleculeBondPipeline_ = nil;
    depthState_ = nil; depthStateNoWrite_ = nil;
    depthTexture_ = nil; shaderLibrary_ = nil;
    runtimeLibraries_.clear();
    commandQueue_ = nil; device_ = nil;
}
