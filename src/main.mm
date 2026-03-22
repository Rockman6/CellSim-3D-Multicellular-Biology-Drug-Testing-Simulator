#import <Cocoa/Cocoa.h>
#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#import <UniformTypeIdentifiers/UniformTypeIdentifiers.h>

#define GLFW_INCLUDE_NONE
#define GLFW_EXPOSE_NATIVE_COCOA
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_metal.h"
#include "implot.h"

#include "gpu/MetalContext.h"
#include "gpu/ShaderTypes.h"
#include "render/Camera.h"
#include "render/MeshLibrary.h"
#include "render/GLBLoader.h"
#include "simulation/Simulation.h"

#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <map>

// ══════════════════════════════════════════════════════════════════════════
//  CellSim — Full Simulation + ImGui Research UI
// ══════════════════════════════════════════════════════════════════════════

static const int INIT_WIDTH  = 1600;
static const int INIT_HEIGHT = 1000;

static Camera       gCamera;
static MetalContext  gCtx;
static MeshLibrary  gMeshes;
static Simulation   gSim;
static std::string  gModelDir;

// Rendering
static std::vector<CellInstance> gCellInstances;
static id<MTLBuffer> gCellBuffer = nil;
static std::map<std::string, LoadedMesh> gOrganelleMeshes;
static int gSelectedCell = 0;

// Time-series recording (growing vectors, no circular buffer issues)
struct TimeSeriesData {
    std::vector<float> time;        // bio-hours
    std::vector<float> population;
    std::vector<float> proliferating;
    std::vector<float> quiescent;
    std::vector<float> apoptotic;
    std::vector<float> necrotic;
    std::vector<float> avgATP;
    std::vector<float> avgStress;
    std::vector<float> glycolyticPct;
    std::vector<float> phaseG1, phaseS, phaseG2, phaseM;
    std::vector<float> divisions;   // cumulative
    std::vector<float> deaths;      // cumulative
    std::vector<float> viability;   // % viable vs pre-drug count
    std::vector<float> avgDrugDamage;
    float sampleTimer = 0;

    void sample(const Simulation& sim) {
        float bioH = sim.bioTime / 3600.0f;
        time.push_back(bioH);
        population.push_back((float)sim.statAlive);
        proliferating.push_back((float)sim.statProlif);
        quiescent.push_back((float)sim.statQuiescent);
        apoptotic.push_back((float)sim.statApoptotic);
        necrotic.push_back((float)sim.statNecrotic);
        avgATP.push_back(sim.statAvgATP);
        // Compute avg stress
        float sumStress = 0;
        for (auto& c : sim.cells) if (c.alive) sumStress += c.stress;
        avgStress.push_back(sim.statAlive > 0 ? sumStress / sim.statAlive : 0);
        glycolyticPct.push_back(sim.statAlive > 0 ? (float)sim.statGlycolytic / sim.statAlive * 100 : 0);
        float total = fmaxf(1, (float)sim.statAlive);
        phaseG1.push_back(sim.statPhases[0] / total * 100);
        phaseS.push_back(sim.statPhases[1] / total * 100);
        phaseG2.push_back(sim.statPhases[2] / total * 100);
        phaseM.push_back(sim.statPhases[3] / total * 100);
        divisions.push_back((float)sim.statDivisions);
        deaths.push_back((float)sim.statDeaths);
        viability.push_back(sim.statViability);
        avgDrugDamage.push_back(sim.statAvgDrugDamage);
    }

    int count() const { return (int)time.size(); }

    void exportCSV(const std::string& path) const {
        FILE* f = fopen(path.c_str(), "w");
        if (!f) { printf("[Export] Failed to write %s\n", path.c_str()); return; }
        fprintf(f, "bio_time_h,population,proliferating,quiescent,apoptotic,necrotic,"
                   "avg_ATP,avg_stress,glycolytic_pct,"
                   "phase_G1_pct,phase_S_pct,phase_G2_pct,phase_M_pct,"
                   "cumulative_divisions,cumulative_deaths,"
                   "viability_pct,avg_drug_damage\n");
        for (int i = 0; i < count(); i++) {
            fprintf(f, "%.4f,%g,%g,%g,%g,%g,%.2f,%.2f,%.1f,%.1f,%.1f,%.1f,%.1f,%g,%g,%.1f,%.4f\n",
                    time[i], population[i], proliferating[i], quiescent[i],
                    apoptotic[i], necrotic[i], avgATP[i], avgStress[i], glycolyticPct[i],
                    phaseG1[i], phaseS[i], phaseG2[i], phaseM[i],
                    divisions[i], deaths[i], viability[i], avgDrugDamage[i]);
        }
        fclose(f);
        printf("[Export] Saved %d rows to %s\n", count(), path.c_str());
    }


    // Per-cell snapshot export
    void exportCellSnapshot(const std::string& path, const Simulation& sim) const {
        FILE* f = fopen(path.c_str(), "w");
        if (!f) { printf("[Export] Failed to write %s\n", path.c_str()); return; }
        fprintf(f, "cell_id,pos_x,pos_z,phase,fate,ATP,stress,ROS,biomass,damage,"
                   "telomere,generation,clone_id,glycolytic,mito_health,mito_potential,"
                   "CycD,Rb,E2F,CycE,CycA,CycB,p21,pressure,necrotic,senescent\n");
        for (int i = 0; i < (int)sim.cells.size(); i++) {
            const auto& c = sim.cells[i];
            if (!c.alive) continue;
            const char* phaseN[] = {"G1","S","G2","M"};
            const char* fateN[] = {"undetermined","proliferating","quiescent","apoptotic"};
            fprintf(f, "%d,%.2f,%.2f,%s,%s,%.1f,%.1f,%.1f,%.3f,%.4f,"
                       "%.0f,%d,%d,%d,%.3f,%.1f,"
                       "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%d,%d\n",
                    i, c.position.x, c.position.z,
                    phaseN[c.phase], fateN[c.fate],
                    c.ATP, c.stress, c.ROS, c.biomass, c.damageLevel,
                    c.telomere, c.generation, c.cloneId,
                    c.glycolytic ? 1 : 0, c.mitoHealth, c.mitoPotential,
                    c.cdk.CycD, c.cdk.Rb, c.cdk.E2F, c.cdk.CycE, c.cdk.CycA, c.cdk.CycB, c.cdk.p21,
                    c.localPressure, c.necrotic ? 1 : 0, c.senescent ? 1 : 0);
        }
        fclose(f);
        printf("[Export] Cell snapshot: %d cells to %s\n", sim.statAlive, path.c_str());
    }
};
static TimeSeriesData gTS;
static std::string gExportDir;

// Per-organelle uniforms
struct OrgUniforms {
    simd_float4x4 viewProjection;
    simd_float4x4 model;
    simd_float3   cameraPos;
    float          time;
    simd_float3   lightDir;
    float          pad0;
    simd_float3   baseColor;
    float          emissiveIntensity;
    simd_float3   emissiveColor;
    float          pad1;
};

// ── Matrix helpers ──────────────────────────────────────────────────────
static simd_float4x4 mat4_translation(simd_float3 t) {
    return (simd_float4x4){{
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {t.x,t.y,t.z,1}
    }};
}
static simd_float4x4 mat4_scale(float s) {
    return (simd_float4x4){{
        {s,0,0,0}, {0,s,0,0}, {0,0,s,0}, {0,0,0,1}
    }};
}
static simd_float4x4 mat4_rotateY(float a) {
    float c=cosf(a),s=sinf(a);
    return (simd_float4x4){{{c,0,s,0},{0,1,0,0},{-s,0,c,0},{0,0,0,1}}};
}
static simd_float4x4 mat4_rotateX(float a) {
    float c=cosf(a),s=sinf(a);
    return (simd_float4x4){{{1,0,0,0},{0,c,-s,0},{0,s,c,0},{0,0,0,1}}};
}
static simd_float4x4 mat4_rotateZ(float a) {
    float c=cosf(a),s=sinf(a);
    return (simd_float4x4){{{c,-s,0,0},{s,c,0,0},{0,0,1,0},{0,0,0,1}}};
}

// ── Organelle model matrix ──────────────────────────────────────────────
static simd_float4x4 organelleModelMatrix(
    simd_float3 cellPos, float cellSize,
    const OrganelleConfig& cfg, float timeOffset
) {
    float t = timeOffset;
    simd_float3 jitter = {sinf(t*0.3f)*0.06f, cosf(t*0.5f)*0.04f, sinf(t*0.4f)*0.06f};
    simd_float3 localPos = cfg.position + jitter;

    simd_float4x4 m = mat4_translation(cellPos);
    m = simd_mul(m, mat4_scale(cellSize));
    m = simd_mul(m, mat4_translation(localPos));
    m = simd_mul(m, mat4_rotateX(cfg.rotation.x));
    m = simd_mul(m, mat4_rotateY(cfg.rotation.y));
    m = simd_mul(m, mat4_rotateZ(cfg.rotation.z));
    m = simd_mul(m, mat4_scale(cfg.scale));
    return m;
}

// ── Draw one GLB organelle ──────────────────────────────────────────────
static void drawGLBOrganelle(
    id<MTLRenderCommandEncoder> enc, const OrganelleConfig& cfg,
    simd_float3 cellPos, float cellSize, const Uniforms& baseUni, float toff
) {
    auto it = gOrganelleMeshes.find(cfg.filename);
    if (it == gOrganelleMeshes.end() || !it->second.valid) return;
    const LoadedMesh& mesh = it->second;

    OrgUniforms ou;
    ou.viewProjection = baseUni.viewProjection;
    ou.model = organelleModelMatrix(cellPos, cellSize, cfg, toff);
    ou.cameraPos = baseUni.cameraPos;
    ou.time = baseUni.time;
    ou.lightDir = baseUni.lightDir;
    ou.pad0 = 0;
    ou.baseColor = cfg.color;
    ou.emissiveIntensity = cfg.emissiveIntensity;
    ou.emissiveColor = cfg.emissive;
    ou.pad1 = 0;

    [enc setVertexBuffer:mesh.vertexBuffer offset:0 atIndex:0];
    [enc setVertexBytes:&ou length:sizeof(OrgUniforms) atIndex:1];
    [enc setFragmentBytes:&ou length:sizeof(OrgUniforms) atIndex:1];
    [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle
                    indexCount:mesh.indexCount indexType:MTLIndexTypeUInt32
                   indexBuffer:mesh.indexBuffer indexBufferOffset:0];
}

// ── Sync simulation → rendering instances ───────────────────────────────
static void syncCellInstances() {
    int n = (int)gSim.cells.size();
    gCellInstances.resize(n);
    for (int i = 0; i < n; i++) {
        auto& c = gSim.cells[i];
        auto& inst = gCellInstances[i];
        inst.position = c.position;
        inst.radius = c.radius * c.size;
        inst.phase = (float)c.phase;
        inst.lodLevel = 2;
        inst.pad = c.size; // store for organelle transform

        // Phase + fate dependent colors
        float glow = 1.0f;
        switch (c.fate) {
            case SIM_FATE_PROLIF:
                inst.color = {0.0f, 0.13f, 0.20f, 0.22f}; glow = 2.5f; break;
            case SIM_FATE_QUIESCENT:
                inst.color = {0.0f, 0.06f, 0.13f, 0.12f}; glow = 0.7f; break;
            case SIM_FATE_APOPTOTIC:
                inst.color = {0.15f, 0.02f, 0.0f, 0.28f}; glow = 1.8f; break;
            default: // UNDETERMINED — color by phase
                switch (c.phase) {
                    case 0: inst.color = {0.0f,0.08f,0.15f,0.17f}; glow=1.0f; break;
                    case 1: inst.color = {0.0f,0.10f,0.18f,0.20f}; glow=1.8f; break;
                    case 2: inst.color = {0.06f,0.04f,0.0f,0.17f}; glow=1.5f; break;
                    case 3: inst.color = {0.15f,0.02f,0.0f,0.25f}; glow=2.0f; break;
                }
        }
        // Override for necrotic cells (orange-red glow, swollen)
        if (c.necrotic) {
            inst.color = {0.25f, 0.04f, 0.0f, 0.35f};
            glow = 1.2f;
        }
        inst.glowIntensity = glow;
    }

    // Resize buffer if needed
    size_t bufSize = n * sizeof(CellInstance);
    if (!gCellBuffer || gCellBuffer.length < bufSize) {
        gCellBuffer = [gCtx.device() newBufferWithLength:std::max(bufSize, (size_t)256)
                                                options:MTLResourceStorageModeShared];
    }
    if (n > 0) memcpy(gCellBuffer.contents, gCellInstances.data(), bufSize);
}

// ── GLFW Callbacks ──────────────────────────────────────────────────────
static void mouseButtonCB(GLFWwindow* w, int b, int a, int m) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    double x, y; glfwGetCursorPos(w, &x, &y);
    gCamera.onMouseButton(b, a, x, y);
}
static void cursorPosCB(GLFWwindow* w, double x, double y) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    gCamera.onMouseMove(x, y);
}
static void scrollCB(GLFWwindow* w, double xo, double yo) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    gCamera.onScroll(-yo * 3.0);
}
static void framebufferCB(GLFWwindow* w, int width, int height) {
    gCamera.onResize(width, height);
    gCtx.metalLayer().drawableSize = CGSizeMake(width, height);
    gCtx.recreateDepthTexture(width, height);
}

// ── Load models ─────────────────────────────────────────────────────────
static void loadModels() {
    struct Spec { const char* f; uint32_t maxTris; };
    Spec specs[] = {
        {"nucleus.glb",0}, {"smooth ER.glb",0}, {"rough ER.glb",30000},
        {"golgi apparatus.glb",0}, {"mitochondria.glb",0}
    };
    for (auto& s : specs) {
        std::string path = gModelDir + "/" + s.f;
        auto mesh = GLBLoader::load(gCtx.device(), path, s.maxTris);
        if (mesh.valid) gOrganelleMeshes[s.f] = mesh;
    }
}

// ── ImGui UI ────────────────────────────────────────────────────────────
static void drawUI() {
    // Population panel (top-right)
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(270, 0), ImGuiCond_FirstUseEver);
    ImGui::Begin("Population Stats");
    {
        float bioH = gSim.bioTime / 3600;
        float bioM = fmodf(gSim.bioTime / 60, 60);
        ImGui::Text("BIO-TIME: %.0fh %02.0fm", floorf(bioH), bioM);
        ImGui::Separator();
        ImGui::Text("CELLS: %d", gSim.statAlive);
        ImGui::TextColored(ImVec4(0,1,0.6f,1), "  Proliferating: %d", gSim.statProlif);
        ImGui::TextColored(ImVec4(0.3f,0.5f,1,1), "  Quiescent: %d", gSim.statQuiescent);
        ImGui::TextColored(ImVec4(1,0.3f,0.1f,1), "  Apoptotic: %d", gSim.statApoptotic);
        ImGui::Separator();
        ImGui::Text("AVG ATP: %.1f", gSim.statAvgATP);
        ImGui::Text("Divisions: %d  Deaths: %d", gSim.statDivisions, gSim.statDeaths);
        ImGui::Separator();
        ImGui::TextColored(ImVec4(1,0.5f,0,1), "  Necrotic: %d", gSim.statNecrotic);
        ImGui::TextColored(ImVec4(1,0.7f,0.3f,1), "  Glycolytic: %d", gSim.statGlycolytic);
        ImGui::Separator();
        ImGui::Text("PHASE DISTRIBUTION:");
        float total = fmaxf(1, (float)gSim.statAlive);
        ImGui::TextColored(ImVec4(0.3f,0.5f,1,1),  "  G1: %d (%.0f%%)", gSim.statPhases[0], gSim.statPhases[0]/total*100);
        ImGui::TextColored(ImVec4(0.2f,1,0.6f,1),   "  S:  %d (%.0f%%)", gSim.statPhases[1], gSim.statPhases[1]/total*100);
        ImGui::TextColored(ImVec4(1,0.8f,0.2f,1),   "  G2: %d (%.0f%%)", gSim.statPhases[2], gSim.statPhases[2]/total*100);
        ImGui::TextColored(ImVec4(1,0.3f,0.6f,1),   "  M:  %d (%.0f%%)", gSim.statPhases[3], gSim.statPhases[3]/total*100);
    }
    ImGui::End();

    // Environment controls (bottom-left)
    ImGui::SetNextWindowPos(ImVec2(10, ImGui::GetIO().DisplaySize.y - 160), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(250, 0), ImGuiCond_FirstUseEver);
    ImGui::Begin("Environment Controls");
    {
        float o2pct = gSim.envO2 * 100;
        float glupct = gSim.envGlucose * 100;
        if (ImGui::SliderFloat("O2 %", &o2pct, 0, 100)) gSim.envO2 = o2pct / 100;
        if (ImGui::SliderFloat("Glucose %", &glupct, 0, 100)) gSim.envGlucose = glupct / 100;
        ImGui::Separator();
        // Logarithmic speed slider (0.1× to 100×)
        float logSpeed = log10f(gSim.timeScale);
        if (ImGui::SliderFloat("##speed", &logSpeed, -1.0f, 2.0f, "")) {
            gSim.timeScale = powf(10.0f, logSpeed);
        }
        float bioMinPerSec = gSim.timeScale * BIO_MIN_PER_SEC;
        if (bioMinPerSec < 1.0f)
            ImGui::Text("%.1fx | %.0fs bio/s", gSim.timeScale, bioMinPerSec * 60);
        else if (bioMinPerSec >= 60)
            ImGui::Text("%.1fx | %.1fh bio/s", gSim.timeScale, bioMinPerSec / 60);
        else
            ImGui::Text("%.1fx | %.1fm bio/s", gSim.timeScale, bioMinPerSec);

        // Speed presets
        if (ImGui::SmallButton("0.5x")) gSim.timeScale = 0.5f; ImGui::SameLine();
        if (ImGui::SmallButton("1x")) gSim.timeScale = 1.0f; ImGui::SameLine();
        if (ImGui::SmallButton("5x")) gSim.timeScale = 5.0f; ImGui::SameLine();
        if (ImGui::SmallButton("20x")) gSim.timeScale = 20.0f; ImGui::SameLine();
        if (ImGui::SmallButton("50x")) gSim.timeScale = 50.0f;

        if (ImGui::Button(gSim.paused ? ">> RESUME" : "|| PAUSE", ImVec2(-1, 0)))
            gSim.paused = !gSim.paused;

        ImGui::Separator();
        ImGui::TextColored(ImVec4(0.4f,0.6f,0.8f,0.6f), "FPS: %.0f | Cells: %d",
                           ImGui::GetIO().Framerate, gSim.statAlive);
    }
    ImGui::End();

    // Time-series plots
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(420, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Population Dynamics");
    {
        int n = gTS.count();
        ImGui::Text("Data points: %d", n);

        if (n > 1) {
            const float* t = gTS.time.data();

            if (ImPlot::BeginPlot("Cell Count", ImVec2(-1, 130))) {
                ImPlot::SetupAxes("Bio-time (h)", "Cells");
                ImPlot::PlotLine("Total", t, gTS.population.data(), n);
                ImPlot::PlotLine("Prolif", t, gTS.proliferating.data(), n);
                ImPlot::PlotLine("Quiesc", t, gTS.quiescent.data(), n);
                ImPlot::EndPlot();
            }

            if (ImPlot::BeginPlot("Metabolism", ImVec2(-1, 130))) {
                ImPlot::SetupAxes("Bio-time (h)", "");
                ImPlot::PlotLine("Avg ATP", t, gTS.avgATP.data(), n);
                ImPlot::PlotLine("Avg Stress", t, gTS.avgStress.data(), n);
                ImPlot::PlotLine("Glycolytic %", t, gTS.glycolyticPct.data(), n);
                ImPlot::EndPlot();
            }

            if (ImPlot::BeginPlot("Phase Distribution %", ImVec2(-1, 130))) {
                ImPlot::SetupAxes("Bio-time (h)", "%");
                ImPlot::PlotLine("G1", t, gTS.phaseG1.data(), n);
                ImPlot::PlotLine("S", t, gTS.phaseS.data(), n);
                ImPlot::PlotLine("G2", t, gTS.phaseG2.data(), n);
                ImPlot::PlotLine("M", t, gTS.phaseM.data(), n);
                ImPlot::EndPlot();
            }
        } else {
            ImGui::TextColored(ImVec4(0.5f,0.5f,0.5f,1), "Waiting for data...");
        }

        ImGui::Separator();
        if (ImGui::Button("EXPORT CSV (Save As...)", ImVec2(-1, 28))) {
            // Default timestamped filename
            time_t now = ::time(nullptr);
            struct tm* t = localtime(&now);
            char fname_buf[128];
            snprintf(fname_buf, sizeof(fname_buf), "cellsim_%04d%02d%02d_%02d%02d%02d.csv",
                     t->tm_year+1900, t->tm_mon+1, t->tm_mday,
                     t->tm_hour, t->tm_min, t->tm_sec);
            NSString* defaultName = [NSString stringWithUTF8String:fname_buf];

            // Pop up native macOS Save dialog
            dispatch_async(dispatch_get_main_queue(), ^{
                NSSavePanel* panel = [NSSavePanel savePanel];
                [panel setTitle:@"Export CellSim Data"];
                [panel setNameFieldStringValue:defaultName];
                [panel setAllowedContentTypes:@[[UTType typeWithFilenameExtension:@"csv"]]];
                // Default to project exports folder
                [panel setDirectoryURL:[NSURL fileURLWithPath:
                    [NSString stringWithUTF8String:gExportDir.c_str()]]];
                [panel setCanCreateDirectories:YES];

                if ([panel runModal] == NSModalResponseOK) {
                    NSURL* url = [panel URL];
                    std::string savePath = [[url path] UTF8String];
                    gTS.exportCSV(savePath);
                    // Remember the chosen directory for next time
                    gExportDir = [[[url path] stringByDeletingLastPathComponent] UTF8String];
                }
            });
        }
        // Per-cell snapshot export
        if (ImGui::Button("EXPORT CELL SNAPSHOT", ImVec2(-1, 24))) {
            time_t now2 = ::time(nullptr);
            struct tm* t2 = localtime(&now2);
            char fname2[128];
            snprintf(fname2, sizeof(fname2), "cellsim_cells_%04d%02d%02d_%02d%02d%02d.csv",
                     t2->tm_year+1900, t2->tm_mon+1, t2->tm_mday,
                     t2->tm_hour, t2->tm_min, t2->tm_sec);
            NSString* defaultName2 = [NSString stringWithUTF8String:fname2];
            dispatch_async(dispatch_get_main_queue(), ^{
                NSSavePanel* panel = [NSSavePanel savePanel];
                [panel setTitle:@"Export Cell Snapshot"];
                [panel setNameFieldStringValue:defaultName2];
                [panel setAllowedContentTypes:@[[UTType typeWithFilenameExtension:@"csv"]]];
                [panel setDirectoryURL:[NSURL fileURLWithPath:
                    [NSString stringWithUTF8String:gExportDir.c_str()]]];
                [panel setCanCreateDirectories:YES];
                if ([panel runModal] == NSModalResponseOK) {
                    std::string savePath = [[[panel URL] path] UTF8String];
                    gTS.exportCellSnapshot(savePath, gSim);
                    gExportDir = [[[[panel URL] path] stringByDeletingLastPathComponent] UTF8String];
                }
            });
        }
        ImGui::TextColored(ImVec4(0.4f,0.6f,0.8f,0.5f), "Dir: %s",
                           gExportDir.substr(gExportDir.rfind('/')+1).c_str());
    }
    ImGui::End();

    // Selected cell detail
    if (gSelectedCell >= 0 && gSelectedCell < (int)gSim.cells.size()) {
        ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 280, 350), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(270, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Cell Research Panel");
        {
            auto& c = gSim.cells[gSelectedCell];
            ImGui::SliderInt("Cell #", &gSelectedCell, 0, (int)gSim.cells.size()-1);
            ImGui::Separator();

            const char* phaseNames[] = {"G1","S","G2","M"};
            const char* fateNames[] = {"Undetermined","Proliferating","Quiescent","Apoptotic"};
            ImVec4 phaseColors[] = {{0.3f,0.5f,1,1},{0.2f,1,0.6f,1},{1,0.8f,0.2f,1},{1,0.3f,0.6f,1}};

            ImGui::TextColored(phaseColors[c.phase], "Phase: %s", phaseNames[c.phase]);
            ImGui::Text("Fate: %s", fateNames[c.fate]);
            ImGui::ProgressBar(c.cycleProgress, ImVec2(-1,12), "Cycle");

            ImGui::Separator();
            ImGui::Text("CDK/CYCLIN (Novak-Tyson):");
            ImGui::ProgressBar(c.cdk.CycD/1.5f, ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycD %.3f", c.cdk.CycD);
            ImGui::ProgressBar(c.cdk.Rb,         ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Rb   %.3f", c.cdk.Rb);
            ImGui::ProgressBar(c.cdk.E2F,        ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("E2F  %.3f", c.cdk.E2F);
            ImGui::ProgressBar(c.cdk.CycE/1.5f,  ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycE %.3f", c.cdk.CycE);
            ImGui::ProgressBar(c.cdk.CycA/1.5f,  ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycA %.3f", c.cdk.CycA);
            ImGui::ProgressBar(c.cdk.CycB/1.5f,  ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("CycB %.3f", c.cdk.CycB);
            ImGui::ProgressBar(c.cdk.p21,        ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("p21  %.3f", c.cdk.p21);

            ImGui::Separator();
            ImGui::Text("METABOLISM:");
            ImGui::ProgressBar(c.ATP/100,      ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("ATP    %.1f", c.ATP);
            ImGui::ProgressBar(c.stress/100,   ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Stress %.1f", c.stress);
            ImGui::ProgressBar(c.ROS/100,      ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("ROS    %.1f", c.ROS);
            ImGui::ProgressBar(c.biomass/2.3f, ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Biomass %.2f", c.biomass);
            ImGui::ProgressBar(c.damageLevel/2,ImVec2(-1,8), ""); ImGui::SameLine(); ImGui::Text("Damage %.3f", c.damageLevel);

            ImGui::Separator();
            ImGui::Text("MITOCHONDRIA:");
            ImGui::ProgressBar((c.mitoPotential-40)/180, ImVec2(-1,8), ""); ImGui::SameLine();
            ImGui::Text("Psi %.0f mV", c.mitoPotential);
            ImGui::ProgressBar(c.mitoHealth, ImVec2(-1,8), ""); ImGui::SameLine();
            ImGui::Text("Health %.2f", c.mitoHealth);
            ImGui::TextColored(c.glycolytic ? ImVec4(1,0.5f,0.2f,1) : ImVec4(0.5f,0.9f,0.6f,1),
                "Mode: %s", c.glycolytic ? "GLYCOLYTIC (Warburg)" : "OXIDATIVE");
            if (c.necrotic) ImGui::TextColored(ImVec4(1,0.3f,0,1), "!! NECROTIC (O2 critical)");

            ImGui::Separator();
            ImGui::Text("GENOME:");
            ImGui::Text("  Gen:%d Clone:%d", c.generation, c.cloneId);
            ImGui::Text("  Telomere: %.0f bp %s", c.telomere, c.senescent ? "[SENESCENT]" : "");
            ImGui::Text("  gly:%.2f pro:%.2f ros:%.2f rep:%.2f",
                         c.glycolysisBias, c.prolifBias, c.rosTolerance, c.repairRate);
            ImGui::Text("  Pressure: %.2f  Hypoxia: %.0fs", c.localPressure, c.hypoxiaTimer);
            ImGui::Text("  ATP danger: %.0f/%.0fs", c.atpDangerTimer, ATP_DANGER_DURATION);

            // Drug response for selected cell
            if (gSim.activeDrugIdx > 0) {
                ImGui::Separator();
                ImGui::TextColored(ImVec4(1,0.6f,0,1), "DRUG RESPONSE:");
                ImGui::Text("  Internal: %.3f µM", c.drugInternal);
                ImGui::ProgressBar(c.drugDamage, ImVec2(-1,8), ""); ImGui::SameLine();
                ImGui::Text("Damage %.3f", c.drugDamage);
                ImGui::Text("  Resistant: %s", c.drugResistant ? "YES (MDR)" : "no");
            }
        }
        ImGui::End();
    }

    // ══════════════════════════════════════════════════════════════════
    //  DRUG TREATMENT PANEL
    // ══════════════════════════════════════════════════════════════════
    ImGui::SetNextWindowPos(ImVec2(10, ImGui::GetIO().DisplaySize.y - 350), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(280, 0), ImGuiCond_FirstUseEver);
    ImGui::Begin("Drug Treatment");
    {
        // Drug selector
        static int selectedDrug = 0;
        const char* drugNames[DRUG_COUNT];
        for (int i = 0; i < DRUG_COUNT; i++) drugNames[i] = DRUG_LIBRARY[i].name;
        ImGui::Combo("Drug", &selectedDrug, drugNames, DRUG_COUNT);

        if (selectedDrug > 0) {
            const Drug& d = DRUG_LIBRARY[selectedDrug];
            ImGui::TextColored(ImVec4(0.6f,0.8f,1,0.7f), "EC50: %.3f µM  Hill: %.1f", d.EC50, d.hillCoeff);
            const char* moaNames[] = {"Anti-Prolif","Pro-Apoptosis","DNA Damage","Mito Toxin"};
            if (d.mechanism >= 0) ImGui::Text("MOA: %s", moaNames[d.mechanism]);
            if (d.mechanism2 >= 0) ImGui::Text("MOA2: %s", moaNames[d.mechanism2]);
        }

        // Concentration slider (log scale)
        static float logConc = 0; // log10(µM)
        ImGui::SliderFloat("log[C] µM", &logConc, -3.0f, 2.0f, "%.1f");
        float conc = powf(10.0f, logConc);
        ImGui::Text("Concentration: %.4f µM", conc);

        ImGui::Separator();

        // Application buttons
        if (ImGui::Button("Apply Uniform", ImVec2(-1, 28))) {
            gSim.applyDrugUniform(selectedDrug, conc);
        }
        ImGui::TextColored(ImVec4(0.5f,0.7f,0.9f,0.5f), "Click plate to inject (right-click)");

        if (ImGui::Button("WASH OUT", ImVec2(-1, 24))) {
            gSim.washOutDrug();
        }

        ImGui::Separator();

        // Drug status
        if (gSim.activeDrugIdx > 0) {
            ImGui::TextColored(ImVec4(1,0.8f,0,1), "Active: %s", DRUG_LIBRARY[gSim.activeDrugIdx].name);
            ImGui::Text("Viability: %.1f%%", gSim.statViability);
            ImGui::ProgressBar(gSim.statViability / 100.0f, ImVec2(-1, 12), "");
            ImGui::Text("Avg Drug Damage: %.3f", gSim.statAvgDrugDamage);
        } else {
            ImGui::TextColored(ImVec4(0.4f,0.4f,0.4f,1), "No drug active");
        }
    }
    ImGui::End();
}

// ── Render frame ────────────────────────────────────────────────────────
static void renderFrame(float time, float dt) {
    @autoreleasepool {
        // Update simulation
        gSim.update(dt);
        syncCellInstances();

        // Time-series sampling (every 0.5 real seconds)
        gTS.sampleTimer += dt;
        if (gTS.sampleTimer > 0.5f) {
            gTS.sampleTimer = 0;
            gTS.sample(gSim);
        }

        id<CAMetalDrawable> drawable = [gCtx.metalLayer() nextDrawable];
        if (!drawable) return;

        int cellCount = (int)gCellInstances.size();
        if (cellCount == 0) return;

        Uniforms uni;
        uni.viewProjection = gCamera.getViewProjection();
        uni.model = matrix_identity_float4x4;
        uni.cameraPos = gCamera.getPosition();
        uni.time = time;
        uni.lightDir = simd_normalize(simd_make_float3(0.5f, 0.9f, 0.6f));
        uni.pad0 = 0;

        id<MTLCommandBuffer> cmd = [gCtx.commandQueue() commandBuffer];
        MTLRenderPassDescriptor* pass = [MTLRenderPassDescriptor renderPassDescriptor];
        pass.colorAttachments[0].texture = drawable.texture;
        pass.colorAttachments[0].loadAction = MTLLoadActionClear;
        pass.colorAttachments[0].storeAction = MTLStoreActionStore;
        pass.colorAttachments[0].clearColor = MTLClearColorMake(0.004, 0.008, 0.031, 1.0);
        pass.depthAttachment.texture = gCtx.depthTexture();
        pass.depthAttachment.loadAction = MTLLoadActionClear;
        pass.depthAttachment.storeAction = MTLStoreActionDontCare;
        pass.depthAttachment.clearDepth = 1.0;

        id<MTLRenderCommandEncoder> enc = [cmd renderCommandEncoderWithDescriptor:pass];

        // 1. Substrate
        {
            [enc setDepthStencilState:gCtx.depthState()];
            const MeshData& sub = gMeshes.substrate();
            Uniforms su = uni;
            simd_float3 fp = {0, FLOOR_Y, 0};
            su.model = mat4_translation(fp);
            [enc setRenderPipelineState:gCtx.substratePipeline()];
            [enc setVertexBuffer:sub.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBytes:&su length:sizeof(Uniforms) atIndex:2];
            [enc setFragmentBytes:&su length:sizeof(Uniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle indexCount:sub.indexCount
                             indexType:MTLIndexTypeUInt32 indexBuffer:sub.indexBuffer indexBufferOffset:0];
        }

        // 2. GLB Organelles (only for nearby cells — LOD optimization)
        {
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            [enc setRenderPipelineState:gCtx.glbOrganellePipeline()];
            OrganelleConfigs orgCfg = OrganelleConfigs::defaults();

            simd_float3 camPos = gCamera.getPosition();
            for (int i = 0; i < cellCount; i++) {
                simd_float3 cp = gCellInstances[i].position;
                float cs = gCellInstances[i].pad;
                float toff = time + (float)i * 7.3f;

                // LOD: only render organelles for cells within 60 units of camera
                float dx = cp.x-camPos.x, dy = cp.y-camPos.y, dz = cp.z-camPos.z;
                if (dx*dx+dy*dy+dz*dz > 3600) continue;

                drawGLBOrganelle(enc, orgCfg.nucleus, cp, cs, uni, toff);
                drawGLBOrganelle(enc, orgCfg.smoothER, cp, cs, uni, toff+1);
                drawGLBOrganelle(enc, orgCfg.roughER, cp, cs, uni, toff+2);
                drawGLBOrganelle(enc, orgCfg.golgi, cp, cs, uni, toff+3);
                for (int m = 0; m < 3; m++)
                    drawGLBOrganelle(enc, orgCfg.mito[m], cp, cs, uni, toff+4+m);
            }
        }

        // 3. Cell membranes
        {
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            const MeshData& sphere = gMeshes.sphereLOD(2);
            [enc setRenderPipelineState:gCtx.cellPipeline()];
            [enc setVertexBuffer:sphere.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBuffer:gCellBuffer offset:0 atIndex:1];
            [enc setVertexBytes:&uni length:sizeof(Uniforms) atIndex:2];
            [enc setFragmentBytes:&uni length:sizeof(Uniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeTriangle indexCount:sphere.indexCount
                             indexType:MTLIndexTypeUInt32 indexBuffer:sphere.indexBuffer
                     indexBufferOffset:0 instanceCount:cellCount];
        }

        // 4. Wireframe
        {
            [enc setDepthStencilState:gCtx.depthStateNoWrite()];
            const MeshData& sphere = gMeshes.sphereLOD(1);
            [enc setRenderPipelineState:gCtx.wirePipeline()];
            [enc setVertexBuffer:sphere.vertexBuffer offset:0 atIndex:0];
            [enc setVertexBuffer:gCellBuffer offset:0 atIndex:1];
            [enc setVertexBytes:&uni length:sizeof(Uniforms) atIndex:2];
            [enc drawIndexedPrimitives:MTLPrimitiveTypeLine indexCount:sphere.indexCount
                             indexType:MTLIndexTypeUInt32 indexBuffer:sphere.indexBuffer
                     indexBufferOffset:0 instanceCount:cellCount];
        }

        // 5. ImGui render
        {
            ImGui_ImplMetal_NewFrame(pass);
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            drawUI();
            ImGui::Render();
            ImGui_ImplMetal_RenderDrawData(ImGui::GetDrawData(), cmd, enc);
        }

        [enc endEncoding];
        [cmd presentDrawable:drawable];
        [cmd commit];
    }
}

// ── Main ────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    @autoreleasepool {
        printf("╔══════════════════════════════════════════════════════════╗\n");
        printf("║  CellSim — Computational Cell Biology Simulator    ║\n");
        printf("║  Full Simulation · CDK/Cyclin ODE · Fick Diffusion     ║\n");
        printf("║  Fate Decision · Cell Division · ImGui Research UI     ║\n");
        printf("╚══════════════════════════════════════════════════════════╝\n\n");

        NSString* execPath = [[NSBundle mainBundle] executablePath];
        NSString* execDir  = [execPath stringByDeletingLastPathComponent];
        NSString* modelDir = [[execDir stringByDeletingLastPathComponent]
                              stringByAppendingPathComponent:@"assets/models"];
        gModelDir = [modelDir UTF8String];

        // Export directory = project exports folder
        NSString* exportDir = [[execDir stringByDeletingLastPathComponent]
                               stringByAppendingPathComponent:@"exports"];
        [[NSFileManager defaultManager] createDirectoryAtPath:exportDir
            withIntermediateDirectories:YES attributes:nil error:nil];
        gExportDir = [exportDir UTF8String];

        if (!glfwInit()) return 1;
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        GLFWwindow* window = glfwCreateWindow(INIT_WIDTH, INIT_HEIGHT,
            "CellSim — Computational Cell Biology Simulator", nullptr, nullptr);
        if (!window) { glfwTerminate(); return 1; }

        NSWindow* nsWin = glfwGetCocoaWindow(window);
        CAMetalLayer* metalLayer = [CAMetalLayer layer];
        metalLayer.contentsScale = nsWin.backingScaleFactor;
        nsWin.contentView.layer = metalLayer;
        nsWin.contentView.wantsLayer = YES;

        int fbW, fbH;
        glfwGetFramebufferSize(window, &fbW, &fbH);
        metalLayer.drawableSize = CGSizeMake(fbW, fbH);

        if (!gCtx.init(metalLayer)) { glfwTerminate(); return 1; }
        gCtx.recreateDepthTexture(fbW, fbH);
        if (!gMeshes.init(gCtx.device())) { glfwTerminate(); return 1; }
        gCamera.onResize(fbW, fbH);

        // ImGui setup
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGui::StyleColorsDark();
        ImGuiStyle& style = ImGui::GetStyle();
        style.WindowRounding = 4; style.FrameRounding = 2;
        style.Colors[ImGuiCol_WindowBg] = ImVec4(0.01f, 0.03f, 0.06f, 0.85f);
        style.Colors[ImGuiCol_TitleBg] = ImVec4(0.0f, 0.04f, 0.08f, 0.9f);
        style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.0f, 0.06f, 0.12f, 0.95f);
        style.Colors[ImGuiCol_FrameBg] = ImVec4(0.0f, 0.04f, 0.08f, 0.7f);
        style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.0f, 0.5f, 0.8f, 1.0f);

        // Set our callbacks FIRST, then ImGui chains on top with install_callbacks=true
        glfwSetMouseButtonCallback(window, mouseButtonCB);
        glfwSetCursorPosCallback(window, cursorPosCB);
        glfwSetScrollCallback(window, scrollCB);
        glfwSetFramebufferSizeCallback(window, framebufferCB);

        ImGui_ImplGlfw_InitForOther(window, true);
        ImGui_ImplMetal_Init(gCtx.device());

        // Load models
        loadModels();

        // Init simulation
        gSim.init();

        printf("[CellSim] Simulation started: %d cells · %d models loaded\n",
               INIT_CELLS, (int)gOrganelleMeshes.size());
        printf("[CellSim] Drag:rotate  Right-drag:pan  Scroll:zoom  ESC:quit\n");

        float lastTime = (float)glfwGetTime();
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
                glfwSetWindowShouldClose(window, true);
            if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
                gSim.paused = !gSim.paused;

            float now = (float)glfwGetTime();
            float dt = fminf(now - lastTime, 0.05f);
            lastTime = now;
            renderFrame(now, dt);
        }

        ImGui_ImplMetal_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        gCtx.shutdown();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 0;
    }
}
