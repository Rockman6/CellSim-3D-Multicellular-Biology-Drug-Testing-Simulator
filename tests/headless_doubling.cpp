// ══════════════════════════════════════════════════════════════════════════
//  headless_doubling.cpp — biology validator (no GPU, no rendering)
//  -----------------------------------------------------------------------
//  Runs the Simulation core for a configurable bio-time span and writes a
//  CSV of cell-count and phase distribution vs bio-time. This lets us
//  verify cells DO divide and population doubles even when the dev-box
//  GPU edge case prevents the interactive build from running long enough.
//
//  Usage:
//      ./build/cellsim_headless [bio_hours] [bio_seconds_per_step]
//      defaults: 72 bio-hours, dt = 60 bio-seconds per step
//
//  Output:
//      logs/headless_doubling.csv
//      Columns: wall_sec, bio_sec, bio_hours, cell_count,
//               g1, s, g2, m, divisions, deaths
// ══════════════════════════════════════════════════════════════════════════

#include "../src/simulation/Simulation.h"
#include "../src/simulation/TelemetryLog.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static double secondsSince(std::chrono::high_resolution_clock::time_point t0) {
    using namespace std::chrono;
    return duration<double>(high_resolution_clock::now() - t0).count();
}

int main(int argc, char** argv) {
    // ── CLI args ─────────────────────────────────────────────────────────
    float bio_hours_target = (argc >= 2) ? (float)atof(argv[1]) : 72.0f;
    float bio_dt_step      = (argc >= 3) ? (float)atof(argv[2]) :  60.0f;

    printf("[headless] bio-hours target = %.1f, bio_dt step = %.1f s\n",
           bio_hours_target, bio_dt_step);

    // ── Init Simulation ──────────────────────────────────────────────────
    Simulation sim;
    sim.mode = MODE_COLONY;          // need population, not single-cell
    sim.timeScale = 1.0f;            // wall-clock irrelevant; we tick directly
    sim.init(MODE_COLONY);
    printf("[headless] init: %d cells\n", (int)sim.cells.size());

    // Comprehensive telemetry log: writes logs/telemetry_population_*.csv
    // and logs/telemetry_division_*.csv with the full state per sample.
    TelemetryLog tlog;
    tlog.populationSampleBiosec = 60.0f;     // sample once per bio-minute
    char sessionId[64];
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::strftime(sessionId, sizeof(sessionId),
                  "headless_%Y%m%d_%H%M%S", std::localtime(&now));
    if (!tlog.open(sessionId)) {
        fprintf(stderr, "[headless] failed to open telemetry log\n");
    } else {
        printf("[headless] telemetry session id = %s\n", sessionId);
    }

    // To advance bio-time directly we call sim.update with a real_dt that
    // produces the desired bio-step under bioDt(). bioDt = real_dt × ts ×
    // BIO_TIME_SCALE → real_dt = bio_dt / (ts × BIO_TIME_SCALE).
    float real_dt_per_step = bio_dt_step / (sim.timeScale * BIO_TIME_SCALE);

    // ── CSV setup ────────────────────────────────────────────────────────
    FILE* csv = fopen("logs/headless_doubling.csv", "w");
    if (!csv) { fprintf(stderr, "[headless] cannot open csv\n"); return 1; }
    fprintf(csv, "wall_sec,bio_sec,bio_hours,cell_count,"
                 "g1,s,g2,m,divisions,deaths,"
                 "cdk_CycD,cdk_CycE,cdk_CycA,cdk_CycB,cdk_Rb,cdk_p21,"
                 "biomass,ATP,damage,replProgress\n");

    auto wall0 = std::chrono::high_resolution_clock::now();
    int prevDivisions = sim.statDivisions;
    int prevDeaths    = sim.statDeaths;
    int prevAlive     = (int)sim.cells.size();

    // ── Tick loop ────────────────────────────────────────────────────────
    int    sample_every = 60;                     // sample once per N steps
    int    step = 0;
    double bio_sec = 0.0;
    double bio_hours = 0.0;

    printf("[headless] running %.0f bio-h at %.0f s/step → ~%d steps total\n",
           bio_hours_target, bio_dt_step,
           (int)(bio_hours_target * 3600.0 / bio_dt_step));

    while (bio_hours < bio_hours_target) {
        sim.update(real_dt_per_step);
        bio_sec   += bio_dt_step;
        bio_hours  = bio_sec / 3600.0;
        step++;
        // Telemetry tick (self-throttled)
        tlog.tick(sim, (float)bio_sec, (float)secondsSince(wall0));

        if (step % sample_every == 0 || (int)sim.cells.size() != prevAlive) {
            int g1=0, s=0, g2=0, m=0, alive=0;
            for (auto& c : sim.cells) {
                if (!c.alive) continue;
                alive++;
                int p = c.phase;
                if      (p == 0) g1++;
                else if (p == 1) s++;
                else if (p == 2) g2++;
                else if (p == 3) m++;
            }
            double wall = secondsSince(wall0);
            // Average CDK + biomass + ATP across alive cells
            float avgCycD=0, avgCycE=0, avgCycA=0, avgCycB=0;
            float avgRb=0, avgp21=0, avgBio=0, avgATP=0, avgDmg=0, avgRepl=0;
            for (auto& c : sim.cells) {
                if (!c.alive) continue;
                avgCycD += c.cdk.CycD; avgCycE += c.cdk.CycE;
                avgCycA += c.cdk.CycA; avgCycB += c.cdk.CycB;
                avgRb   += c.cdk.Rb;   avgp21  += c.cdk.p21;
                avgBio  += c.biomass;  avgATP  += c.ATP; avgDmg += c.damageLevel;
                avgRepl += (float)c.program.cdogma.replicationProgress;
            }
            float n = alive > 0 ? (float)alive : 1.0f;
            fprintf(csv,
                "%.2f,%.0f,%.3f,%d,%d,%d,%d,%d,%d,%d,"
                "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%.3f,%.3f\n",
                wall, bio_sec, bio_hours, alive,
                g1, s, g2, m,
                sim.statDivisions, sim.statDeaths,
                avgCycD/n, avgCycE/n, avgCycA/n, avgCycB/n,
                avgRb/n, avgp21/n,
                avgBio/n, avgATP/n, avgDmg/n, avgRepl/n);
            fflush(csv);

            if (step % (sample_every * 10) == 0 || alive != prevAlive) {
                printf("[headless] bioh=%5.2f cells=%4d  G1=%d S=%d G2=%d M=%d  "
                       "div+%d  death+%d  (wall %.1fs)\n",
                       bio_hours, alive, g1, s, g2, m,
                       sim.statDivisions - prevDivisions,
                       sim.statDeaths    - prevDeaths,
                       wall);
                prevDivisions = sim.statDivisions;
                prevDeaths    = sim.statDeaths;
                prevAlive     = alive;
            }
        }

        // Progress dot every ~600 steps for liveness while running.
        if (step % 600 == 0) { printf("."); fflush(stdout); }
    }

    fclose(csv);
    tlog.close();
    double wall_total = secondsSince(wall0);
    printf("\n[headless] DONE. %.1f bio-h in %.1f wall-s → CSV at "
           "logs/headless_doubling.csv\n",
           bio_hours_target, wall_total);
    printf("[headless] final cells=%d  total divisions=%d  total deaths=%d\n",
           (int)sim.cells.size(), sim.statDivisions, sim.statDeaths);
    return 0;
}
