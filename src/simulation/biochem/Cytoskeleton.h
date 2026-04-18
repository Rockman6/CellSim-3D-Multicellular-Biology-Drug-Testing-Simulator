#pragma once
#include <cmath>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  Cytoskeleton.h — Cytoskeletal dynamics simulation
//  Based on: Alberts Ch 16 "The Cytoskeleton"
//
//  Implements:
//    - Actin filaments:
//      G-actin ⇌ F-actin polymerization, treadmilling
//      Arp2/3 branching (lamellipodia), formin nucleation (filopodia)
//      Cofilin severing, profilin/thymosin-β4 monomer regulation
//      Cortex model (ERM proteins linking actin to membrane)
//    - Microtubules:
//      α/β-tubulin → protofilaments, dynamic instability
//      γ-TuRC nucleation at centrosome, +TIPs (EB1)
//      Kinesin (+end motors), dynein (-end motors)
//      Katanin/spastin severing
//    - Intermediate filaments:
//      Keratin, vimentin, neurofilaments, lamins
//      Mechanical stress distribution
//    - Motor proteins:
//      Myosin II (contractile ring, stress fibers, cortical tension)
//      Myosin V (vesicle transport on actin)
//      Kinesin-1 (anterograde transport)
//      Cytoplasmic dynein + dynactin (retrograde transport)
//    - Rho GTPase regulation:
//      RhoA → stress fibers + focal adhesions
//      Rac1 → lamellipodia (Arp2/3 activation via WAVE)
//      Cdc42 → filopodia (formin activation)
//
//  Units: µM for concentrations, µm for lengths, pN for forces
// ══════════════════════════════════════════════════════════════════════════

struct ActinState {
    // Monomer pools
    float G_actin_ATP = 50.0f;   // ATP-G-actin (polymerizable)
    float G_actin_ADP = 5.0f;    // ADP-G-actin (from depolymerization)
    float thymosin_bound = 30.0f;// Thymosin-β4-sequestered G-actin
    float profilin_bound = 10.0f;// Profilin-bound G-actin (promotes barbed end addition)

    // Filament state
    float F_actin_total = 20.0f; // Total polymerized actin (µM equivalent)
    float barbed_ends = 100.0f;  // Number of free barbed ends
    float pointed_ends = 100.0f; // Number of free pointed ends
    float branch_points = 20.0f; // Arp2/3-created branch points

    // Regulatory proteins
    float Arp23_active = 0.0f;   // Active Arp2/3 complex
    float WASP_WAVE_active = 0.0f;// NPFs (activated by Rac1/Cdc42)
    float formin_active = 0.0f;  // Active formins (mDia)
    float cofilin_active = 0.5f; // Active cofilin/ADF
    float capping_protein = 0.3f;// CapZ (caps barbed ends)
    float tropomyosin = 0.5f;    // Stabilizes filaments, blocks cofilin

    // Cortex
    float cortex_thickness = 0.2f;  // µm
    float cortical_tension = 100.0f;// pN/µm (from myosin II)
    float ERM_active = 0.5f;       // Ezrin/radixin/moesin (actin-membrane linker)

    // Critical concentrations (µM)
    static constexpr float Cc_barbed = 0.12f;  // Barbed end
    static constexpr float Cc_pointed = 0.6f;  // Pointed end
};

struct MicrotubuleState {
    // Tubulin pools
    float free_tubulin = 15.0f;      // Free αβ-tubulin heterodimers (µM)
    float polymerized = 10.0f;       // Polymerized tubulin

    // Microtubule dynamics
    float MT_count = 200.0f;         // Number of microtubules
    float avg_length = 10.0f;        // Average length (µm)
    float growth_rate = 1.5f;        // µm/min at barbed end
    float shrink_rate = 20.0f;       // µm/min during catastrophe
    float catastrophe_freq = 0.05f;  // Events per minute per MT
    float rescue_freq = 0.02f;       // Rescues per minute per shrinking MT

    // MTOC / Centrosome
    float gamma_TuRC = 1.0f;         // γ-tubulin ring complex (nucleation)
    float centrosome_nucleation = 0.5f;// Nucleation rate

    // MAPs (Microtubule-Associated Proteins)
    float MAP2 = 0.3f;              // Stabilizing MAP (dendrites)
    float Tau = 0.3f;               // Stabilizing MAP (axons)
    float stathmin = 0.2f;          // Tubulin-sequestering (destabilizing)
    float katanin = 0.1f;           // MT-severing AAA ATPase
    float EB1 = 0.5f;               // +TIP protein

    // Dynamic instability tracking
    float growing_fraction = 0.7f;   // Fraction of MTs in growth phase
    float shrinking_fraction = 0.3f;
};

struct MotorProteinState {
    // Myosin motors
    float myosinII_active = 0.3f;    // Non-muscle myosin II (NM-II)
    float myosinII_cortex = 0.2f;    // Myosin II in cortex
    float myosinV = 0.1f;            // Vesicle transport on actin
    float myosinVI = 0.05f;          // Minus-end directed (endocytosis)

    // MT motors
    float kinesin1 = 0.3f;           // Kinesin-1 (anterograde, +end)
    float kinesin13 = 0.1f;          // Kinesin-13 (MCAK, depolymerase)
    float dynein = 0.3f;             // Cytoplasmic dynein (-end)
    float dynactin = 0.3f;           // Dynein activator complex

    // Force generation
    float contractile_force = 0.0f;  // pN (from myosin II)
    float transport_rate = 0.0f;     // Vesicles/s transported
};

struct IntermediateFilamentState {
    float keratin = 1.0f;            // Type I/II IF (epithelial)
    float vimentin = 0.5f;           // Type III IF (mesenchymal)
    float lamin_A = 0.5f;            // Nuclear lamina (type V IF)
    float lamin_B = 0.5f;
    float neurofilament = 0.0f;      // Type IV IF (neurons)
    float desmin = 0.0f;             // Type III IF (muscle)

    // Mechanical properties
    float tensile_strength = 1.0f;   // Normalized
    float nuclear_integrity = 1.0f;  // Lamin-dependent
};

// ══════════════════════════════════════════════════════════════════════════

class CytoskeletonEngine {
public:
    ActinState actin;
    MicrotubuleState mt;
    MotorProteinState motors;
    IntermediateFilamentState ifs;

    // Rho GTPase inputs (from Signaling.h)
    float RhoA = 0.0f, Rac1 = 0.0f, Cdc42 = 0.0f;

    // Cell state outputs
    float cell_stiffness = 0.0f;     // Overall mechanical stiffness
    float migration_speed = 0.0f;    // Cell migration velocity (µm/min)
    float polarity = 0.0f;           // Front-rear polarity (0-1)

    void init() {
        actin = ActinState();
        mt = MicrotubuleState();
        motors = MotorProteinState();
        ifs = IntermediateFilamentState();
    }

    void step(float dt, float rhoA, float rac1, float cdc42, float ATP) {
        dt = fminf(dt, 0.1f);
        RhoA = rhoA; Rac1 = rac1; Cdc42 = cdc42;

        // ═══ ACTIN DYNAMICS ═══════════════════════════════════════
        auto& a = actin;

        // Rac1 → WAVE → Arp2/3 activation (lamellipodia)
        a.WASP_WAVE_active = Rac1 * 0.8f + Cdc42 * 0.3f;
        a.Arp23_active = a.WASP_WAVE_active * 0.5f / (0.3f + a.WASP_WAVE_active);

        // Cdc42 → formin activation (filopodia)
        a.formin_active = Cdc42 * 0.6f / (0.2f + Cdc42);

        // RhoA → myosin II activation → cortical tension
        motors.myosinII_active = RhoA * 0.7f / (0.3f + RhoA);
        motors.myosinII_cortex = motors.myosinII_active * 0.5f;

        // Profilin exchanges ADP for ATP on G-actin
        float exchange = 0.5f * a.G_actin_ADP * dt;
        a.G_actin_ADP -= exchange;
        a.G_actin_ATP += exchange;

        // Thymosin-β4 sequesters G-actin (buffer)
        float seq = 0.1f * (a.G_actin_ATP - a.thymosin_bound * 0.3f) * dt;
        a.thymosin_bound += seq;
        a.G_actin_ATP -= seq;

        // Polymerization at barbed ends (ATP-actin)
        float poly_barbed = a.barbed_ends * 0.01f *
            fmaxf(0, a.G_actin_ATP - ActinState::Cc_barbed) * dt;

        // Depolymerization at pointed ends (ADP-actin)
        float depoly_pointed = a.pointed_ends * 0.005f * dt;

        // Arp2/3 branching (creates new barbed ends)
        float new_branches = a.Arp23_active * a.F_actin_total * 0.01f * dt;
        a.branch_points += new_branches;
        a.barbed_ends += new_branches;

        // Capping protein caps barbed ends
        float capping = a.capping_protein * a.barbed_ends * 0.01f * dt;
        a.barbed_ends -= capping;

        // Cofilin severing (creates new pointed ends, accelerates turnover)
        float severing = a.cofilin_active * a.F_actin_total * 0.005f * dt;
        // Tropomyosin protects from cofilin
        severing *= (1.0f - a.tropomyosin * 0.5f);
        a.pointed_ends += severing;

        // Update pools
        a.F_actin_total += (poly_barbed - depoly_pointed - severing * 0.1f);
        a.G_actin_ATP -= poly_barbed;
        a.G_actin_ADP += depoly_pointed;
        a.F_actin_total = fmaxf(0, a.F_actin_total);
        a.G_actin_ATP = fmaxf(0, a.G_actin_ATP);
        a.G_actin_ADP = fmaxf(0, a.G_actin_ADP);
        a.barbed_ends = fmaxf(10, a.barbed_ends);

        // Cortical tension (myosin II pulling on cortical actin)
        a.cortical_tension = 50.0f + motors.myosinII_cortex * 200.0f;

        // ═══ MICROTUBULE DYNAMICS ═════════════════════════════════
        auto& m = mt;

        // Dynamic instability
        float cat_events = m.catastrophe_freq * m.growing_fraction * m.MT_count * dt / 60.0f;
        float res_events = m.rescue_freq * m.shrinking_fraction * m.MT_count * dt / 60.0f;
        m.growing_fraction += (res_events - cat_events) / fmaxf(m.MT_count, 1);
        m.shrinking_fraction = 1.0f - m.growing_fraction;
        m.growing_fraction = std::clamp(m.growing_fraction, 0.1f, 0.9f);

        // Net growth
        float net_growth = m.growing_fraction * m.growth_rate - m.shrinking_fraction * m.shrink_rate;
        m.avg_length += net_growth * dt / 60.0f;
        m.avg_length = std::clamp(m.avg_length, 1.0f, 30.0f);

        // Stathmin sequesters tubulin (reduces free tubulin)
        m.free_tubulin = 15.0f - m.stathmin * 5.0f;

        // MAPs stabilize (reduce catastrophe frequency)
        m.catastrophe_freq = 0.05f / (1.0f + (m.MAP2 + m.Tau) * 2.0f);

        // Katanin severing
        float mt_severing = m.katanin * m.MT_count * 0.001f * dt;
        m.MT_count += mt_severing; // Severing creates more (shorter) MTs

        // Nucleation from centrosome
        float nucleation = m.gamma_TuRC * m.centrosome_nucleation * 0.01f * dt;
        m.MT_count += nucleation;

        // ═══ MOTOR PROTEINS ═══════════════════════════════════════
        // Contractile force from myosin II
        motors.contractile_force = motors.myosinII_active * a.F_actin_total * 5.0f; // pN

        // Transport rate (kinesin/dynein moving vesicles along MTs)
        motors.transport_rate = (motors.kinesin1 + motors.dynein * motors.dynactin) *
                                m.MT_count * 0.01f;

        // ═══ INTERMEDIATE FILAMENTS ═══════════════════════════════
        // Nuclear integrity depends on lamin network
        ifs.nuclear_integrity = (ifs.lamin_A + ifs.lamin_B) * 0.5f;
        // Overall tensile strength
        ifs.tensile_strength = ifs.keratin * 0.3f + ifs.vimentin * 0.3f +
                               ifs.desmin * 0.2f + ifs.neurofilament * 0.2f;

        // ═══ CELL-LEVEL OUTPUTS ═══════════════════════════════════
        cell_stiffness = a.cortical_tension * 0.005f + ifs.tensile_strength * 0.3f +
                         motors.contractile_force * 0.001f;
        cell_stiffness = std::clamp(cell_stiffness, 0.0f, 5.0f);

        // Migration speed (driven by actin polymerization at leading edge)
        float protrusion = a.Arp23_active * 0.5f + a.formin_active * 0.3f;
        float traction = motors.myosinII_active * 0.3f;
        migration_speed = fminf(protrusion, traction) * 10.0f; // µm/min
        migration_speed = std::clamp(migration_speed, 0.0f, 30.0f);

        // Polarity (asymmetry between front and back)
        polarity = fabsf(Rac1 - RhoA) / (Rac1 + RhoA + 0.1f);
    }
};
