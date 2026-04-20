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
    auto runScenario = [&](const char* label, float forced_damage,
                           float& out_peak_p53, float& out_peak_mdm2,
                           float& out_last_p53) {
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
                         "ATM_active,p53_protein,MDM2_protein,cdk_p21\n");
        }

        float peak_p53 = 0.0f, peak_mdm2 = 0.0f, last_p53 = 0.0f;
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
                        fprintf(csv, "%.3f,%zu,%.4f,%.4f,%.5f,%.5f,%.4f\n",
                                bio_sec / 3600.0, i,
                                c.damageLevel, c.ATM_active,
                                c.p53_protein, c.MDM2_protein, c.cdk.p21);
                    }
                }
                const auto& c0 = sim.cells[0];
                if (c0.alive) {
                    peak_p53  = std::max(peak_p53,  c0.p53_protein);
                    peak_mdm2 = std::max(peak_mdm2, c0.MDM2_protein);
                    last_p53  = c0.p53_protein;
                }
            }
        }
        if (csv) fclose(csv);
        out_peak_p53  = peak_p53;
        out_peak_mdm2 = peak_mdm2;
        out_last_p53  = last_p53;
    };

    // ── Two-scenario assay: no damage (control) vs sustained damage ────
    float ctl_peak_p53 = 0, ctl_peak_mdm2 = 0, ctl_last_p53 = 0;
    float dmg_peak_p53 = 0, dmg_peak_mdm2 = 0, dmg_last_p53 = 0;
    runScenario("control", 0.00f, ctl_peak_p53, ctl_peak_mdm2, ctl_last_p53);
    runScenario("damage",  0.30f, dmg_peak_p53, dmg_peak_mdm2, dmg_last_p53);

    // Baselines observed at init (analytical steady state of the ODE).
    const float baseline_p53  = 0.089f;
    const float baseline_mdm2 = 0.21f;

    // ── Assertions ─────────────────────────────────────────────────────
    bool okA = (ctl_peak_p53 > baseline_p53 * 0.75f) &&
               (ctl_peak_p53 < baseline_p53 * 1.25f);  // ±25 % window
    bool okB = (dmg_peak_p53  >= baseline_p53 * 1.5f);
    bool okC = (dmg_peak_mdm2 >= baseline_mdm2 * 1.5f);
    bool okD = (dmg_peak_p53  >  ctl_peak_p53 * 1.4f); // specificity

    printf("\n=== p53 STRESS-RESPONSE RESULT (seed=0x%08x) ===\n", sim_seed);
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

    bool all = okA && okB && okC && okD;
    printf("=> %s\n", all ? "PASS" : "FAIL");
    printf("CSVs: logs/p53_pulsing_control.csv, logs/p53_pulsing_damage.csv\n");
    return all ? 0 : 1;
}
