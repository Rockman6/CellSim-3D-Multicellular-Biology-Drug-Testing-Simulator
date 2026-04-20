#include "TargetLibrary.h"
#include "../Simulation.h"
#include <cmath>

TargetLibrary gTargets;

// ══════════════════════════════════════════════════════════════════════════
//  Target modulators — each one reads a drug-bound occupancy 0..1 and
//  perturbs an existing SimCell biology field. No new MOA code paths,
//  only rate / pool adjustments that the sim already consumes.
//
//  All of these are reversible: when occupancy drops to 0 (drug washes
//  out) the affected field returns to normal via the native dynamics.
// ══════════════════════════════════════════════════════════════════════════

static inline float softMul(float a, float factor) {
    if (factor < 0.0f) factor = 0.0f;
    if (factor > 1.0f) factor = 1.0f;
    return a * factor;
}

// TGT_CDK1 — Cyclin B blockade. Inhibits M-entry.
static void modCdkB(SimCell& c, float occ, float /*dt*/) {
    c.cdk.CycB = softMul(c.cdk.CycB, 1.0f - 0.95f * occ);
}

// TGT_CDK2 — CyclinE + CyclinA blockade. Inhibits S-entry.
static void modCdkE(SimCell& c, float occ, float /*dt*/) {
    c.cdk.CycE = softMul(c.cdk.CycE, 1.0f - 0.95f * occ);
    c.cdk.CycA = softMul(c.cdk.CycA, 1.0f - 0.90f * occ);
}

// TGT_RB — Rb lock. Prevents Rb phosphorylation → G1 arrest via p21 rise.
static void modRb(SimCell& c, float occ, float dt_biosec) {
    c.cdk.Rb  = fminf(1.0f, c.cdk.Rb + 0.02f * occ * dt_biosec / 3600.0f);
    c.cdk.p21 = fminf(1.0f, c.cdk.p21 + 0.10f * occ * dt_biosec / 3600.0f);
}

// TGT_TOPO2 — TopII poison. Stalls replication → DNA damage rises.
// Rates 4× stronger than first pass — TopII poisons like doxorubicin
// produce DSBs at a real rate of ~1 lesion/min/cell in culture
// (Tewey 1984; our prior 0.02/bio-h was ~50× too slow).
static void modTopo2(SimCell& c, float occ, float dt_biosec) {
    c.program.cdogma.replicationActive = (occ < 0.3f);
    c.damageLevel = fminf(2.0f, c.damageLevel + 0.10f * occ * dt_biosec / 3600.0f);
    c.ROS         = fminf(100.0f, c.ROS + 2.5f * occ * dt_biosec / 3600.0f);
    c.apo.state.p53_active = fminf(1.0f, c.apo.state.p53_active + 0.05f * occ * dt_biosec / 3600.0f);
}

// TGT_TUBULIN — spindle disruption. Mitosis SAC fails → apoptosis.
static void modTubulin(SimCell& c, float occ, float dt_biosec) {
    // Force mitosis SAC to accumulate errors by injecting a transient
    // damage signal that propagates to the apoptosis engine.
    if (c.phase == 3) {
        c.damageLevel = fminf(2.0f, c.damageLevel + 0.03f * occ * dt_biosec / 3600.0f);
    }
}

// TGT_BCL2 — anti-apoptotic Bcl-2 inhibited. MOMP threshold drops.
static void modBcl2(SimCell& c, float occ, float /*dt*/) {
    c.apo.state.Bcl2 = softMul(c.apo.state.Bcl2, 1.0f - 0.90f * occ);
    c.apo.state.BclXL = softMul(c.apo.state.BclXL, 1.0f - 0.50f * occ);
}

// TGT_DNA_INTERCALATION — DNA intercalator / crosslinker. Lesions
// accumulate; p53 activation follows; cell dies by intrinsic pathway.
// When damage saturates we directly inject Bax activation (p53-mediated
// BAX transcription bypasses Bcl-2 buffering) so caspase-3 commits on
// physiological timescale rather than waiting for slow Bcl-2 neutralisation.
static void modDnaIntercalate(SimCell& c, float occ, float dt_biosec) {
    c.damageLevel = fminf(2.0f, c.damageLevel + 0.20f * occ * dt_biosec / 3600.0f);
    c.apo.state.p53_active = fminf(1.0f, c.apo.state.p53_active + 0.10f * occ * dt_biosec / 3600.0f);
    c.apo.state.Puma = fminf(1.0f, c.apo.state.Puma + 0.03f * occ * dt_biosec / 3600.0f);
    // Saturating DNA damage → direct Bax transcription (p53 → BAX promoter).
    if (c.damageLevel > 1.2f && occ > 0.15f) {
        c.apo.state.Bax_active = fminf(1.0f,
            c.apo.state.Bax_active + 0.04f * occ * dt_biosec / 3600.0f);
        c.apo.state.MOMP_pores = fminf(1.0f,
            c.apo.state.MOMP_pores + 0.02f * occ * dt_biosec / 3600.0f);
    }
}

// TGT_HDAC — chromatin opens. Gene expression reshape (gentle stress).
static void modHdac(SimCell& c, float occ, float /*dt*/) {
    c.stress = fminf(100.0f, c.stress + 2.0f * occ);
}

// TGT_KINASE_EGFR — proliferation signalling blocked.
static void modKinaseEgfr(SimCell& c, float occ, float /*dt*/) {
    c.cdk.CycD = softMul(c.cdk.CycD, 1.0f - 0.80f * occ);
}

// TGT_PROTEASOME — proteasome inhibition. Misfolded proteins accumulate → ER stress.
static void modProteasome(SimCell& c, float occ, float dt_biosec) {
    c.stress = fminf(100.0f, c.stress + 1.0f * occ * dt_biosec / 3600.0f);
    c.apo.state.Puma = fminf(1.0f, c.apo.state.Puma + 0.01f * occ * dt_biosec / 3600.0f);
}

// TGT_RIBOSOME — translation shutdown. Biomass synthesis drops.
static void modRibosome(SimCell& c, float occ, float /*dt*/) {
    c.biomass = softMul(c.biomass, 1.0f - 0.0005f * occ);  // gentle continuous drop
}

// TGT_POL2 — transcription blocked. Cycle protein supply cut.
static void modPol2(SimCell& c, float occ, float /*dt*/) {
    c.cdk.CycD = softMul(c.cdk.CycD, 1.0f - 0.50f * occ);
    c.cdk.CycE = softMul(c.cdk.CycE, 1.0f - 0.30f * occ);
}

// TGT_MITO_COMPLEX_I — mitochondrial ETC inhibition. ΔΨm drops → MOMP bias.
// Rotenone / Complex-I poisons collapse ΔΨm over ~15 min at saturating
// concentrations (Ramsay 1997). Rate boosted 3× so mitoPotential crosses
// MITO_K = 80 threshold within ~1 bio-h at full occupancy.
static void modMitoCmplxI(SimCell& c, float occ, float dt_biosec) {
    c.mitoPotential = fmaxf(40.0f, c.mitoPotential - 60.0f * occ * dt_biosec / 3600.0f);
    c.ROS = fminf(100.0f, c.ROS + 1.5f * occ * dt_biosec / 3600.0f);
}

void registerBuiltinTargets() {
    auto reg = [](const char* id, const char* kind,
                   float logpmin, float logpmax,
                   float mwmin, float mwmax,
                   int hbd, int hba, int arom,
                   void (*mod)(SimCell&, float, float)) {
        TargetProfile p;
        p.id = id;
        p.pref_logP_min = logpmin; p.pref_logP_max = logpmax;
        p.pref_mw_min = mwmin; p.pref_mw_max = mwmax;
        p.pref_hbd = hbd; p.pref_hba = hba; p.pref_aromatic = arom;
        gTargets.registerTarget(id, kind, p, mod);
    };
    reg("TGT_CDK1",               "enzyme",      1.0f, 4.5f, 250, 500, 1, 4, 2, modCdkB);
    reg("TGT_CDK2",               "enzyme",      1.0f, 4.5f, 250, 500, 1, 4, 2, modCdkE);
    reg("TGT_RB",                 "structural",  1.0f, 4.0f, 300, 600, 2, 5, 1, modRb);
    reg("TGT_TOPO2",              "enzyme",     -1.0f, 3.5f, 400, 800, 2, 8, 3, modTopo2);
    reg("TGT_TUBULIN",            "structural",  3.0f, 6.0f, 500, 900, 2, 10, 2, modTubulin);
    reg("TGT_BCL2",               "structural",  2.0f, 5.5f, 400, 800, 3, 6, 2, modBcl2);
    reg("TGT_DNA_INTERCALATION",  "nucleic_acid",-2.0f, 2.0f, 150, 500, 2, 6, 1, modDnaIntercalate);
    reg("TGT_HDAC",               "enzyme",      1.0f, 4.0f, 250, 500, 1, 5, 1, modHdac);
    reg("TGT_KINASE_EGFR",        "enzyme",      2.0f, 5.0f, 350, 600, 1, 5, 2, modKinaseEgfr);
    reg("TGT_PROTEASOME",         "enzyme",      1.0f, 4.0f, 300, 700, 2, 6, 1, modProteasome);
    reg("TGT_RIBOSOME",           "enzyme",      0.0f, 3.5f, 400, 1000,3, 8, 1, modRibosome);
    reg("TGT_POL2",               "enzyme",      1.0f, 4.0f, 300, 700, 2, 6, 1, modPol2);
    reg("TGT_MITO_COMPLEX_I",     "enzyme",      2.0f, 5.0f, 250, 600, 1, 4, 1, modMitoCmplxI);
}
