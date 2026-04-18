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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

static double secondsSince(std::chrono::high_resolution_clock::time_point t0) {
    using namespace std::chrono;
    return duration<double>(high_resolution_clock::now() - t0).count();
}

int main(int argc, char** argv) {
    // ── CLI args ─────────────────────────────────────────────────────────
    float bio_hours_target = (argc >= 2) ? (float)atof(argv[1]) : 72.0f;
    float bio_dt_step      = (argc >= 3) ? (float)atof(argv[2]) :  60.0f;
    int   init_cells_req   = (argc >= 4) ? atoi(argv[3])        :  -1;
    const char* ref_csv    = (argc >= 5) ? argv[4]
        : "data/reference/growth_curves/ctc_hela_cellcount_seq02.csv";

    printf("[headless] bio-hours target = %.1f, bio_dt step = %.1f s\n",
           bio_hours_target, bio_dt_step);
    printf("[headless] reference CSV = %s\n", ref_csv);

    // ── Init Simulation ──────────────────────────────────────────────────
    Simulation sim;
    sim.mode = MODE_COLONY;          // need population, not single-cell
    sim.timeScale = 1.0f;            // wall-clock irrelevant; we tick directly
    sim.init(MODE_COLONY);
    // Optionally seed extra cells to match reference experiment's initial
    // density (CTC seq02 starts at 125 cells; seq01 at 43).
    if (init_cells_req > 0) {
        while ((int)sim.cells.size() < init_cells_req) {
            SimCell c;
            float r = sqrtf((float)rand()/RAND_MAX) * SCENE_BOUND * 0.50f;
            float a = (float)rand()/RAND_MAX * 2.0f * (float)M_PI;
            simd_float3 p = {r*cosf(a), FLOOR_Y, r*sinf(a)};
            c.init(p, (int)sim.cells.size());
            c.cellUid = sim.allocateCellUid();
            sim.cells.push_back(c);
        }
    }
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

    // ── Compare our count-vs-time curve against a reference experiment ──
    // Reads a CSV with (frame,hours,cell_count,...) format (CTC HeLa
    // tracking datasets use this schema) and writes a combined CSV with
    // both sim and reference columns aligned on the reference's time
    // grid. This lets the user plot both curves on one graph to visually
    // check the growth kinetics match HeLa wet-lab data.
    {
        // Re-read our own CSV into memory to interpolate.
        FILE* simF = fopen("logs/headless_doubling.csv", "r");
        FILE* refF = fopen(ref_csv, "r");
        FILE* outF = fopen("logs/compare_vs_reference.csv", "w");
        if (!refF) {
            printf("[headless] reference CSV not found, skipping comparison\n");
        } else if (!simF || !outF) {
            printf("[headless] could not open compare output\n");
            if (simF) fclose(simF); if (refF) fclose(refF); if (outF) fclose(outF);
        } else {
            // Load sim (bio_hours, cell_count) pairs
            std::vector<std::pair<double,int>> sim_pts;
            char buf[2048]; (void)fgets(buf, sizeof(buf), simF);  // header
            while (fgets(buf, sizeof(buf), simF)) {
                double wall_s, bio_s, bio_h; int alive;
                if (sscanf(buf, "%lf,%lf,%lf,%d", &wall_s, &bio_s, &bio_h, &alive) == 4) {
                    sim_pts.push_back({bio_h, alive});
                }
            }
            fclose(simF);

            fprintf(outF, "hours,sim_cells,ref_cells,ratio,abs_err,rel_err_pct\n");
            // Read reference, interpolate sim at each ref time
            int matched = 0;
            double sum_rel_err = 0.0;
            double sum_abs_err = 0.0;
            while (fgets(buf, sizeof(buf), refF)) {
                if (buf[0] == '#' || buf[0] == 'f') continue;  // skip comments/header
                int frame; double h; int ref_cells; int divs;
                if (sscanf(buf, "%d,%lf,%d,%d", &frame, &h, &ref_cells, &divs) < 3) continue;
                // Interpolate sim at time h
                int sim_cells = 0;
                if (!sim_pts.empty()) {
                    if (h <= sim_pts.front().first) sim_cells = sim_pts.front().second;
                    else if (h >= sim_pts.back().first) sim_cells = sim_pts.back().second;
                    else {
                        for (size_t i = 1; i < sim_pts.size(); i++) {
                            if (sim_pts[i].first >= h) {
                                double t0 = sim_pts[i-1].first, t1 = sim_pts[i].first;
                                int    c0 = sim_pts[i-1].second, c1 = sim_pts[i].second;
                                double u = (h - t0) / (t1 - t0 + 1e-9);
                                sim_cells = (int)(c0 + (c1 - c0) * u + 0.5);
                                break;
                            }
                        }
                    }
                }
                double ratio   = (ref_cells > 0) ? (double)sim_cells / ref_cells : 0.0;
                double abs_err = (double)(sim_cells - ref_cells);
                double rel_err = (ref_cells > 0) ? 100.0 * abs_err / ref_cells : 0.0;
                fprintf(outF, "%.2f,%d,%d,%.4f,%.1f,%.2f\n",
                        h, sim_cells, ref_cells, ratio, abs_err, rel_err);
                matched++;
                sum_rel_err += std::fabs(rel_err);
                sum_abs_err += std::fabs(abs_err);
            }
            fclose(refF); fclose(outF);
            if (matched > 0) {
                printf("[headless] compare_vs_reference.csv: %d rows matched\n", matched);
                printf("[headless] mean |relative error| = %.1f%% vs reference\n",
                       sum_rel_err / matched);
                printf("[headless] mean |absolute error| = %.1f cells\n",
                       sum_abs_err / matched);
            }
        }
    }
    return 0;
}
