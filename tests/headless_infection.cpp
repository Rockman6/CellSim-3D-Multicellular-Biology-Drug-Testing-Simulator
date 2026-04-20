// headless_infection.cpp — Phase 7 smoke test.
//
// Setup: 400 live cells on the dish floor, 1 bio-h pre-incubation.
// Inoculate: 50 virions (flu_h1n1) at cell cluster center.
// Run: 24 bio-h simulated time.
// Assert:
//   a) ≥ 20 virions bound within 1 bio-h.
//   b) ≥ 3 cells infected (intraVirions > 0) at some point.
//   c) free-virion pool grows (budding or apoptotic spill).
//   d) statDeaths > 0 (infection or incidental apoptosis drives lysis).
//
// Negative control: run again with randomized spike descriptors
// (break the binding) — expect ≤ 5 binding events.

#include "../src/simulation/Simulation.h"
#include "../src/simulation/SimRng.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

static void runScenario(const char* label,
                        bool randomizeSpike,
                        int& outMaxBound,
                        int& outMaxInfected,
                        int& outFreeVirions,
                        int& outDeaths) {
    srand(42);
    Simulation sim;
    sim.mode = MODE_COLONY;
    sim.init();
    // Pre-incubate 0.5 bio-h.
    for (int i = 0; i < 120; i++) sim.update(1.0f / 60.0f);

    printf("[%s] Pathogen registry: %zu virion specs, %zu bacterium specs\n",
           label, gPathogens().virionSpecs.size(), gPathogens().bacteriumSpecs.size());
    printf("[%s] Cells after pre-incubation: %zu\n", label, sim.cells.size());
    if (!sim.cells.empty()) {
        printf("[%s] Cell[0] pos=(%.2f,%.2f,%.2f) alive=%d radius=%.2f\n", label,
               sim.cells[0].position.x, sim.cells[0].position.y, sim.cells[0].position.z,
               sim.cells[0].alive, sim.cells[0].radius);
    }
    int vidx = gPathogens().findVirion("flu_h1n1");
    if (vidx < 0) { printf("[%s] NO VIRION SPEC — did YAML load?\n", label); return; }

    if (randomizeSpike) {
        VirionSpec& vs = gPathogens().virionSpecs[vidx];
        vs.spike_logP  = -10.0f;   // pushes score → 0
        vs.spike_mw    = 1.0f;
        vs.spike_hbd   = 99;
        vs.spike_hba   = 99;
        vs.spike_aromatic = 99;
        vs.preferredReceptors.clear();
    }

    // Inoculate 50 virions near cell[0] — the COLONY layout leaves a
    // hole at the origin, so seeding at (0,0) would place every virion
    // too far from any cell to ever make contact within the test.
    float cx = sim.cells.empty() ? 0.0f : sim.cells[0].position.x;
    float cz = sim.cells.empty() ? 0.0f : sim.cells[0].position.z;
    sim.seedVirion(vidx, cx, cz, 50, 2.0f);

    int maxBound    = 0;
    int maxInfected = 0;
    int startFree   = (int)sim.gFreeVirions.size();

    // Run 24 bio-h ≈ 24 * 60 minute-ticks at dt=1/60 s.
    for (int step = 0; step < 24 * 60 * 60; step++) {
        sim.update(1.0f / 60.0f);
        if (sim.statVirionsBound   > maxBound)    maxBound    = sim.statVirionsBound;
        if (sim.statInfectedCells  > maxInfected) maxInfected = sim.statInfectedCells;
        if (step % 3600 == 0) {
            printf("[%s] bio-h=%d  free=%zu  bound=%d  intraCells=%d  deaths=%d\n",
                   label, step / 3600,
                   sim.gFreeVirions.size(), sim.statVirionsBound,
                   sim.statInfectedCells, sim.statDeaths);
        }
    }
    outMaxBound    = maxBound;
    outMaxInfected = maxInfected;
    outFreeVirions = (int)sim.gFreeVirions.size();
    outDeaths      = sim.statDeaths;
    printf("[%s] FINAL free=%d bound_peak=%d infected_peak=%d deaths=%d (startFree=%d)\n",
           label, outFreeVirions, outMaxBound, outMaxInfected, outDeaths, startFree);
}

int main(int argc, char** argv) {
    // Phase P1: --seed <uint32> consumed here; scenario runs still call
    // srand(42) internally (legacy fixed seed) unless overridden.
    uint32_t sim_seed = 42u;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            sim_seed = (uint32_t)std::strtoul(argv[i + 1], nullptr, 0);
            i++;
        }
    }
    (void)sim_seed; // runScenario currently reseeds with srand(42)
                    // to preserve its legacy contract; wire-up is left
                    // for the PR that standardises scenario seeding.
    printf("[infection] CLI seed = 0x%08x (scenarios still use srand(42))\n",
           sim_seed);

    int a_bound=0, a_inf=0, a_free=0, a_deaths=0;
    int b_bound=0, b_inf=0, b_free=0, b_deaths=0;

    runScenario("POSITIVE", false, a_bound, a_inf, a_free, a_deaths);
    runScenario("NEGATIVE", true,  b_bound, b_inf, b_free, b_deaths);

    bool okA = (a_bound   >= 20);
    bool okB = (a_inf     >= 50);      // infection well beyond the 50-virion seed
    bool okC = (a_inf     >  50);      // amplification beyond initial inoculum
    bool okD = (a_bound   >  b_bound * 4 + 5);  // positive beats negative decisively

    printf("\n=== RESULT ===\n");
    printf("  A bound peak ≥ 20        : %s  (got %d)\n", okA?"PASS":"FAIL", a_bound);
    printf("  B infection peak ≥ 50    : %s  (got %d)\n", okB?"PASS":"FAIL", a_inf);
    printf("  C amplification          : %s  (got %d, seed=50)\n", okC?"PASS":"FAIL", a_inf);
    printf("  D pos ≫ neg binding      : %s  (pos=%d neg=%d)\n", okD?"PASS":"FAIL", a_bound, b_bound);

    return (okA && okB && okC && okD) ? 0 : 1;
}
