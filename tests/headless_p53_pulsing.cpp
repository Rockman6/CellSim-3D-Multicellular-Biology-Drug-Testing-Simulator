// ══════════════════════════════════════════════════════════════════════════
//  headless_p53_pulsing.cpp — Phase G2 biology validator.
//  ────────────────────────────────────────────────────────────────────────
//  Forces sustained DNA damage in a small dish and records per-tick
//  p53, MDM2, ATM time series. Writes logs/p53_pulsing.csv. Exits
//  non-zero if:
//    (a) no-damage control: p53 stays within 25 % of baseline
//        (i.e. the mechanism is specific, not constitutively on).
//    (b) with damageLevel = 0.30, p53 rises to ≥ 1.5× baseline
//        (ATM-gated MDM2-escape is working).
//    (c) with damageLevel = 0.30, MDM2 rises to ≥ 1.5× baseline
//        (p53-driven MDM2 transcription closes the loop).
//    (d) p53(damage=0.30) ≫ p53(damage=0) — specificity to stress.
//
//  This is the first PR in the CellSim history that tests a specific
//  published biological *mechanism* (ATM → p53 stabilisation → MDM2
//  transactivation) rather than population counts. The initial 2-var
//  (p53, MDM2) ODE settles to a stable elevated p53 under sustained
//  damage, not the oscillatory regime of Purvis 2012 — reproducing
//  such pulses requires an MDM2-mRNA intermediate (Geva-Zatorsky 2006
//  3-variable model), which is the target of the next PR (G3).
//
//  References:
//    Purvis TJ et al. 2012 Science 336:1440 — p53 pulsing.
//    Lahav G et al. 2004 Nat Genet 36:147  — single-cell p53 dynamics.
//    Geva-Zatorsky N et al. 2006 Mol Syst Biol 2:0033 — oscillation.
//    Lev Bar-Or R et al. 2000 PNAS 97:11250 — delay-based oscillator model.
//
//  Usage:
//      ./build/cellsim_p53_pulsing [--seed 0xNN]
// ══════════════════════════════════════════════════════════════════════════

#include "../src/simulation/Simulation.h"
#include "../src/simulation/SimRng.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

int main(int argc, char** argv) {
    uint32_t sim_seed = simrng::DEFAULT_SEED;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            sim_seed = (uint32_t)std::strtoul(argv[i + 1], nullptr, 0);
            i++;
        }
    }
    simrng::seed(sim_seed);
    printf("[p53] seed = 0x%08x\n", sim_seed);

    // ── Scenario runner ─────────────────────────────────────────────────
    struct ScenarioResult {
        float peak_p53   = 0.0f;
        float peak_mdm2  = 0.0f;
        float last_p53   = 0.0f;
        int   peak_count = 0;      // number of local p53 maxima detected
        float mean_period_h = 0.0f; // mean interval between peaks (0 if <2 peaks)
    };
    auto runScenario = [&](const char* label, float forced_damage,
                           ScenarioResult& out) {
        simrng::seed(sim_seed);
        Simulation sim;
        sim.mode = MODE_COLONY;
        sim.timeScale = 1.0f;
        sim.init(MODE_COLONY);
        while ((int)sim.cells.size() > 10) sim.cells.pop_back();

        float real_dt_per_step = 60.0f / (sim.timeScale * BIO_TIME_SCALE);
        const float bio_hours_target = 24.0f;
        int   sample_stride = 6; // sample every 6 bio-min

        FILE* csv = nullptr;
        {
            char path[128];
            std::snprintf(path, sizeof(path),
                          "logs/p53_pulsing_%s.csv", label);
            csv = fopen(path, "w");
        }
        if (csv) {
            fprintf(csv, "bio_hours,cell_id,damageLevel,"
                         "ATM_active,p53_protein,MDM2_mRNA,MDM2_protein,cdk_p21\n");
        }

        std::vector<float> p53_trace;
        std::vector<float> t_trace;

        double bio_sec = 0.0;
        int step = 0;
        while (bio_sec / 3600.0 < bio_hours_target) {
            sim.update(real_dt_per_step);
            bio_sec += 60.0;
            step++;

            // Pin damageLevel each tick so the assay probes steady-state
            // stress response, not damage-repair dynamics.
            for (auto& c : sim.cells) {
                if (c.alive) c.damageLevel = forced_damage;
            }

            if (step % sample_stride == 0 && !sim.cells.empty()) {
                for (size_t i = 0; i < sim.cells.size(); i++) {
                    const auto& c = sim.cells[i];
                    if (!c.alive) continue;
                    if (csv) {
                        fprintf(csv,
                                "%.3f,%zu,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f\n",
                                bio_sec / 3600.0, i,
                                c.damageLevel, c.ATM_active,
                                c.p53_protein, c.MDM2_mRNA,
                                c.MDM2_protein, c.cdk.p21);
                    }
                }
                const auto& c0 = sim.cells[0];
                if (c0.alive) {
                    out.peak_p53  = std::max(out.peak_p53,  c0.p53_protein);
                    out.peak_mdm2 = std::max(out.peak_mdm2, c0.MDM2_protein);
                    out.last_p53  = c0.p53_protein;
                    p53_trace.push_back(c0.p53_protein);
                    t_trace.push_back((float)(bio_sec / 3600.0));
                }
            }
        }
        if (csv) fclose(csv);

        // ── Peak detection on cell-0 p53 trace ────────────────────────
        // Skip first 2 h (initial ATM-activation transient). Use a
        // 3-point local-max rule with peak-height gate = 1.3× baseline.
        // Enforce a minimum inter-peak separation of 2.5 h so noise
        // isn't counted as twin peaks.
        const float baseline_p53  = 0.089f;
        const float peak_min      = baseline_p53 * 1.3f;
        const float min_gap_h     = 2.5f;
        std::vector<float> peak_times;
        size_t i0 = 0;
        while (i0 < t_trace.size() && t_trace[i0] < 2.0f) i0++;
        for (size_t i = i0 + 1; i + 1 < p53_trace.size(); i++) {
            float v = p53_trace[i];
            if (v <= p53_trace[i - 1] || v <= p53_trace[i + 1]) continue;
            if (v < peak_min) continue;
            if (!peak_times.empty() &&
                t_trace[i] - peak_times.back() < min_gap_h) continue;
            peak_times.push_back(t_trace[i]);
        }
        out.peak_count = (int)peak_times.size();
        if (peak_times.size() >= 2) {
            out.mean_period_h =
                (peak_times.back() - peak_times.front()) /
                (float)(peak_times.size() - 1);
        }
    };

    // ── Three-scenario assay: control, mid damage, high damage ─────────
    ScenarioResult ctl, dmg_mid, dmg_high;
    runScenario("control",  0.00f, ctl);
    runScenario("damage",   0.30f, dmg_mid);
    runScenario("damage35", 0.35f, dmg_high);

    float ctl_peak_p53 = ctl.peak_p53,    ctl_peak_mdm2 = ctl.peak_mdm2,
          ctl_last_p53 = ctl.last_p53;
    float dmg_peak_p53 = dmg_mid.peak_p53, dmg_peak_mdm2 = dmg_mid.peak_mdm2,
          dmg_last_p53 = dmg_mid.last_p53;

    // Baselines observed at init (analytical steady state of the ODE).
    const float baseline_p53  = 0.089f;
    const float baseline_mdm2 = 0.21f;

    // ── Assertions ─────────────────────────────────────────────────────
    // G2 gates: mechanism is present and specific.
    bool okA = (ctl_peak_p53 > baseline_p53 * 0.75f) &&
               (ctl_peak_p53 < baseline_p53 * 1.25f);
    bool okB = (dmg_peak_p53  >= baseline_p53 * 1.5f);
    bool okC = (dmg_peak_mdm2 >= baseline_mdm2 * 1.5f);
    bool okD = (dmg_peak_p53  >  ctl_peak_p53 * 1.4f);
    // G3 gates: Purvis-2012-style oscillation under sustained damage.
    bool okE = (dmg_mid.peak_count >= 2);
    bool okF = (dmg_mid.mean_period_h >= 3.0f &&
                dmg_mid.mean_period_h <= 9.0f);
    bool okG = (dmg_high.mean_period_h > 0.0f &&
                dmg_mid.mean_period_h > 0.0f &&
                dmg_high.mean_period_h <= dmg_mid.mean_period_h * 1.05f);
    bool okH = (ctl.peak_count == 0);

    printf("\n=== p53 STRESS-RESPONSE + PULSING RESULT (seed=0x%08x) ===\n",
           sim_seed);
    printf("  baseline p53      : %.4f\n", baseline_p53);
    printf("  baseline MDM2     : %.4f\n", baseline_mdm2);
    printf("  CONTROL (dmg=0):   p53 peak %.4f  MDM2 peak %.4f  last p53 %.4f\n",
           ctl_peak_p53, ctl_peak_mdm2, ctl_last_p53);
    printf("  DAMAGE (dmg=0.30): p53 peak %.4f  MDM2 peak %.4f  last p53 %.4f\n",
           dmg_peak_p53, dmg_peak_mdm2, dmg_last_p53);
    printf("  (a) control p53 within ±25%% of basal  : %s\n", okA ? "PASS" : "FAIL");
    printf("  (b) damage  p53 ≥ 1.5× basal          : %s\n", okB ? "PASS" : "FAIL");
    printf("  (c) damage MDM2 ≥ 1.5× basal          : %s\n", okC ? "PASS" : "FAIL");
    printf("  (d) p53(damage) > 1.4× p53(control)   : %s\n", okD ? "PASS" : "FAIL");

    printf("\n--- G3 Purvis-2012-style oscillation ---\n");
    printf("  control peaks=%d  damage(0.30) peaks=%d period=%.2fh  "
           "damage(0.35) peaks=%d period=%.2fh\n",
           ctl.peak_count, dmg_mid.peak_count, dmg_mid.mean_period_h,
           dmg_high.peak_count, dmg_high.mean_period_h);
    printf("  (e) ≥2 p53 peaks under damage 0.30    : %s\n", okE ? "PASS" : "FAIL");
    printf("  (f) mean period ∈ [3, 9] h            : %s\n", okF ? "PASS" : "FAIL");
    printf("  (g) period(0.35) ≤ period(0.30) (+5%%): %s\n", okG ? "PASS" : "FAIL");
    printf("  (h) control has 0 peaks               : %s\n", okH ? "PASS" : "FAIL");

    bool all = okA && okB && okC && okD && okE && okF && okG && okH;
    printf("=> %s\n", all ? "PASS" : "FAIL");
    printf("CSVs: logs/p53_pulsing_{control,damage,damage35}.csv\n");
    return all ? 0 : 1;
}
