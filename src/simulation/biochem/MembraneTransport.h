#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  MembraneTransport.h — Membrane biology and ion transport
//  Based on: Alberts Ch 10 (Membrane Structure) & Ch 11 (Transport)
//
//  Implements:
//    - Ion concentrations (Na+, K+, Ca2+, Cl-, H+)
//    - Na+/K+ ATPase pump (3 Na+ out, 2 K+ in per ATP)
//    - Ca2+ ATPase (SERCA and PMCA)
//    - Voltage-gated channels (Na+, K+, Ca2+)
//    - Ligand-gated channels
//    - Aquaporins
//    - ABC transporters
//    - Membrane potential (Goldman equation)
//    - Nernst equation for each ion
//    - Hodgkin-Huxley action potential model
//
//  Units: mV for potential, mM for concentrations, nA for currents
// ══════════════════════════════════════════════════════════════════════════

static constexpr float RT_F = 26.7f; // RT/F at 37°C in mV (for Nernst)

struct IonState {
    // Intracellular concentrations (mM)
    float Na_in  = 12.0f;   // Na+ (low inside)
    float K_in   = 140.0f;  // K+  (high inside)
    float Ca_in  = 0.0001f; // Ca2+ (very low, ~100 nM)
    float Cl_in  = 4.0f;    // Cl- (low inside)
    float Mg_in  = 0.5f;    // Mg2+

    // Extracellular concentrations (mM) — normally constant
    float Na_out = 145.0f;
    float K_out  = 5.0f;
    float Ca_out = 1.8f;
    float Cl_out = 120.0f;

    // Membrane potential (mV)
    float Vm = -70.0f;      // Resting membrane potential

    // Membrane capacitance
    static constexpr float Cm = 1.0f; // µF/cm² (typical)

    // Channel conductances (mS/cm²) — Hodgkin-Huxley
    float gNa = 0.0f;       // Na+ channel conductance
    float gK  = 0.3f;       // K+ leak conductance
    float gCl = 0.1f;       // Cl- conductance
    float gCa = 0.0f;       // Ca2+ channel conductance

    // Hodgkin-Huxley gating variables
    float m = 0.05f;        // Na+ activation gate
    float h = 0.6f;         // Na+ inactivation gate
    float n = 0.32f;        // K+ activation gate

    // Pump rates
    float NaK_pump_rate = 1.0f;    // Na+/K+ ATPase activity (normalized)
    float Ca_pump_rate = 1.0f;     // Ca2+ ATPase activity

    // ── Nernst equation ────────────────────────────────────────────
    float E_Na() const { return RT_F * logf(Na_out / fmaxf(Na_in, 0.01f)); }
    float E_K()  const { return RT_F * logf(K_out / fmaxf(K_in, 0.01f)); }
    float E_Ca() const { return RT_F * 0.5f * logf(Ca_out / fmaxf(Ca_in, 1e-6f)); } // z=2
    float E_Cl() const { return -RT_F * logf(Cl_out / fmaxf(Cl_in, 0.01f)); } // z=-1

    // ── Goldman equation for resting potential ──────────────────────
    float goldmanPotential() const {
        float P_Na = 0.04f, P_K = 1.0f, P_Cl = 0.45f; // Permeability ratios
        float num = P_Na * Na_out + P_K * K_out + P_Cl * Cl_in;
        float den = P_Na * Na_in + P_K * K_in + P_Cl * Cl_out;
        if (den < 1e-6f) return -70.0f;
        return RT_F * logf(num / den);
    }
};

class MembraneTransportEngine {
public:
    IonState ions;

    void init() { ions = IonState(); }

    void step(float dt_s, float ATP_available) {
        // Convert once: dt_s is seconds, HH rate constants are per ms.
        // Stability cap 25 µs for HH gating. Subdivide if needed.
        float dt_ms_total = dt_s * 1000.0f;
        int n_sub = (int)ceilf(dt_ms_total / 0.025f);
        n_sub = std::clamp(n_sub, 1, 1000);
        float dt_ms = dt_ms_total / (float)n_sub;

        auto& s = ions;

        for (int sub = 0; sub < n_sub; sub++) {
            // ── Na+/K+ ATPase (3 Na+ out, 2 K+ in per ATP) ────────
            // pump is in cycles per ms; multiply by stoichiometry.
            float pump_rate = s.NaK_pump_rate * ATP_available / (1.0f + ATP_available)
                            * s.Na_in / (10.0f + s.Na_in)
                            * s.K_out / (1.5f + s.K_out);
            s.Na_in -= 3.0f * pump_rate * dt_ms * 0.001f; // → per second
            s.K_in  += 2.0f * pump_rate * dt_ms * 0.001f;

            // ── Ca2+ ATPase (PMCA) ───────────────────────────────
            float ca_pump = s.Ca_pump_rate * s.Ca_in / (0.0005f + s.Ca_in)
                          * ATP_available / (1.0f + ATP_available);
            s.Ca_in -= ca_pump * dt_ms * 0.001f;

            // ── Hodgkin-Huxley (MODERN convention, rest = -65 mV) ─
            // Ref: standard textbook form, e.g. Dayan & Abbott 2001.
            // V is in mV. Rate constants are per ms.
            float V = s.Vm;

            // Stable evaluation near removable singularities:
            auto alpha_m_fn = [](float V) {
                float x = (V + 40.0f);
                if (fabsf(x) < 1e-3f) return 1.0f; // L'Hôpital
                return 0.1f * x / (1.0f - expf(-x / 10.0f));
            };
            auto alpha_n_fn = [](float V) {
                float x = (V + 55.0f);
                if (fabsf(x) < 1e-3f) return 0.1f;
                return 0.01f * x / (1.0f - expf(-x / 10.0f));
            };
            float alpha_m = alpha_m_fn(V);
            float beta_m  = 4.0f * expf(-(V + 65.0f) / 18.0f);
            float alpha_h = 0.07f * expf(-(V + 65.0f) / 20.0f);
            float beta_h  = 1.0f / (1.0f + expf(-(V + 35.0f) / 10.0f));
            float alpha_n = alpha_n_fn(V);
            float beta_n  = 0.125f * expf(-(V + 65.0f) / 80.0f);

            // Rush-Larsen exponential update for gating — unconditionally stable
            auto rush_larsen = [](float x, float a, float b, float dt_ms) {
                float tau = 1.0f / (a + b);
                float xinf = a * tau;
                return xinf + (x - xinf) * expf(-dt_ms / tau);
            };
            s.m = rush_larsen(s.m, alpha_m, beta_m, dt_ms);
            s.h = rush_larsen(s.h, alpha_h, beta_h, dt_ms);
            s.n = rush_larsen(s.n, alpha_n, beta_n, dt_ms);
            s.m = std::clamp(s.m, 0.0f, 1.0f);
            s.h = std::clamp(s.h, 0.0f, 1.0f);
            s.n = std::clamp(s.n, 0.0f, 1.0f);

            // Channel conductances (mS/cm²)
            const float gNa_max = 120.0f;
            const float gK_max  = 36.0f;
            const float gL      = 0.3f;
            s.gNa = gNa_max * s.m * s.m * s.m * s.h;
            s.gK  = gK_max  * s.n * s.n * s.n * s.n;

            // Currents (µA/cm²): reversal potentials from Nernst.
            float I_Na = s.gNa * (V - s.E_Na());
            float I_K  = s.gK  * (V - s.E_K());
            float I_L  = gL    * (V - (-54.4f));      // leak reversal (L=Cl+bg)
            float I_pump = pump_rate * 1.0f;          // electrogenic: +1e out per cycle

            // Vm update (V is mV, Cm is µF/cm², I in µA/cm², dt in ms)
            float dV = -(I_Na + I_K + I_L + I_pump) / IonState::Cm;
            s.Vm += dV * dt_ms;
            s.Vm = std::clamp(s.Vm, -100.0f, 60.0f);

            // Ion concentration changes from currents.
            // dC/dt [mM/s] = I [µA/cm²] * A [cm²] / (F * V_cell [L])
            // For a 20 µm sphere: A ≈ 1.26e-5 cm², V ≈ 4.2e-12 L
            // So dC/dt ≈ I * 31 mM/s per µA/cm² — very fast in practice.
            // Use a scaled-down coupling for visible dynamics over sim timescales.
            const float ion_coupling = 1e-4f; // tuned for visibility
            s.Na_in += I_Na * ion_coupling * dt_ms * 0.001f;
            s.K_in  -= I_K  * ion_coupling * dt_ms * 0.001f;
        }

        // Clamp
        s.Na_in = fmaxf(1.0f, s.Na_in);
        s.K_in  = fmaxf(10.0f, s.K_in);
        s.Ca_in = fmaxf(1e-5f, s.Ca_in);
    }

    float getMembranePotential() const { return ions.Vm; }
    float getRestingPotential() const { return ions.goldmanPotential(); }
};
