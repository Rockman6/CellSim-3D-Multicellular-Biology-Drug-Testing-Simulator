#pragma once

// ══════════════════════════════════════════════════════════════════════════
//  CellBiologyEngine.h — Master integration of ALL biochemistry modules
//  Based on: Alberts "Molecular Biology of the Cell" 7th Ed, Ch 1-24
//
//  This is the unified engine that connects all 16 simulation modules:
//    1.  Metabolism (glycolysis, TCA, OxPhos, fermentation)
//    2.  Signaling (GPCR, MAPK, PI3K-Akt, Wnt, Notch, Hedgehog, NF-κB...)
//    3.  DNARepair (MMR, BER, NER, NHEJ, HR, DDR, telomerase)
//    4.  GeneRegulation (chromatin, histone code, methylation, miRNA)
//    5.  Proteome (folding, chaperones, proteasome, PTMs)
//    6.  MembraneTransport (ion channels, pumps, Hodgkin-Huxley)
//    7.  VesicularTraffic (COPII/COPI/clathrin, SNAREs, autophagy, UPR)
//    8.  Cytoskeleton (actin, microtubules, IFs, motors)
//    9.  Apoptosis (Bcl-2, caspases, death receptors)
//   10.  CellJunctions (cadherins, integrins, tight/gap junctions, ECM)
//   11.  Cancer (oncogenes, tumor suppressors, hallmarks)
//   12.  GenomeModel (23 chromosomes, telomeres, TADs)
//   13.  Development (morphogens, Hox, germ layers, PAR)
//   14.  StemCells (niche, self-renewal, iPSC)
//   15.  Pathogens (bacteria, viruses, microbiota)
//   16.  Immunity (innate TLR/complement/NK, adaptive B/T cells)
//
//  All modules communicate through shared state variables.
//  The step() function runs all modules in the correct order
//  with proper cross-talk between pathways.
// ══════════════════════════════════════════════════════════════════════════

#include "Metabolism.h"
#include "Signaling.h"
#include "DNARepair.h"
#include "GeneRegulation.h"
#include "Proteome.h"
#include "MembraneTransport.h"
#include "VesicularTraffic.h"
#include "Cytoskeleton.h"
#include "Apoptosis.h"
#include "CellJunctions.h"
#include "Cancer.h"
#include "GenomeModel.h"
#include "Development.h"
#include "StemCells.h"
#include "Pathogens.h"
#include "Immunity.h"

class CellBiologyEngine {
public:
    // All 16 subsystems
    MetabolismEngine     metabolism;
    SignalingEngine       signaling;
    DNARepairEngine      dnaRepair;
    GeneRegulationEngine geneRegulation;
    ProteomeEngine       proteome;
    MembraneTransportEngine membrane;
    VesicularTrafficEngine  vesicles;
    CytoskeletonEngine   cytoskeleton;
    ApoptosisEngine      apoptosis;
    CellJunctionsEngine  junctions;
    CancerEngine         cancer;
    GenomeEngine         genome;
    DevelopmentEngine    development;
    StemCellEngine       stemCells;
    PathogenEngine       pathogens;
    ImmunityEngine       immunity;

    bool initialized = false;
    float bioTime = 0.0f; // Biological time in seconds

    void init() {
        metabolism.init();
        signaling.init();
        dnaRepair.init();
        geneRegulation.init();
        proteome.init();
        membrane.init();
        vesicles.init();
        cytoskeleton.init();
        apoptosis.init();
        junctions.init();
        cancer.init();
        genome.init(true); // female by default
        development.init();
        stemCells.init(TISSUE_GENERIC);
        pathogens.init();
        immunity.init();
        initialized = true;
    }

    // ── Main simulation step ────────────────────────────────────────
    // Runs all 16 modules with proper inter-module communication
    void step(float dt) {
        if (!initialized) init();
        dt = fminf(dt, 0.1f); // Max 100ms per step
        bioTime += dt;

        // ════════════════════════════════════════════════════════════
        //  ORDER MATTERS: upstream modules feed downstream ones
        // ════════════════════════════════════════════════════════════

        // 1. METABOLISM — produces ATP, NADH, drives everything
        metabolism.step(dt);
        float ATP = metabolism.pool.ATP;
        float energyCharge = metabolism.pool.energyCharge();

        // 2. MEMBRANE TRANSPORT — ion gradients, membrane potential
        membrane.step(dt, ATP);

        // 3. SIGNALING — integrates external signals, drives decisions
        signaling.step(dt);
        float prolifSignal = signaling.state.proliferation_signal;
        float survivalSignal = signaling.state.survival_signal;
        float growthSignal = signaling.state.growth_signal;

        // 4. GENE REGULATION — transcription/translation driven by signals
        geneRegulation.step(dt,
            signaling.state.ERK_active,      // MAPK signal
            signaling.state.Akt_active,      // PI3K-Akt signal
            dnaRepair.ddr.p53_P,             // p53 from DNA damage
            signaling.state.TCF_active,      // Wnt signal
            dnaRepair.getDamageLevel()        // DNA damage level
        );

        // 5. PROTEOME — protein folding, degradation
        proteome.step(dt, 310.0f, ATP,
            metabolism.pool.amino_acid_pool);

        // 6. DNA REPAIR — detect/repair damage, checkpoint control
        // ROS estimated from NADH leak and O2 availability (simplified)
        float ROS = 0.1f * metabolism.pool.NADH * metabolism.pool.O2;
        int cellPhase = 0; // Would come from cell cycle
        dnaRepair.step(dt, ROS, 0.0f, 0.0f, cellPhase);

        // 7. VESICULAR TRAFFIC — protein sorting, secretion, autophagy
        vesicles.step(dt, ATP, energyCharge,
            signaling.state.mTORC1_active);

        // 8. CYTOSKELETON — actin/MT dynamics, driven by Rho GTPases
        cytoskeleton.step(dt,
            signaling.state.RhoA_GTP,
            signaling.state.Rac1_GTP,
            signaling.state.Cdc42_GTP,
            ATP);

        // 9. APOPTOSIS — death decision based on damage + survival signals
        apoptosis.applySurvival(
            signaling.state.Akt_active,
            signaling.state.NFkB_nuclear);
        apoptosis.applyDNAdamage(dnaRepair.ddr.p53_P);
        apoptosis.step(dt);

        // 10. CELL JUNCTIONS — adhesion (uses Ca2+, tension)
        junctions.step(dt,
            1.8f, // extracellular Ca2+ (mM)
            cytoskeleton.motors.contractile_force * 0.01f,
            metabolism.pool.pH_cytosol(),
            signaling.state.Ca_cytosol);

        // 11. CANCER — mutation tracking, hallmark computation
        // (mutations added externally; cancer.onDivision() called at cell division)

        // 12. GENOME — chromosome state, replication tracking
        // (genome.replicateGenome() called during S-phase)

        // 13. DEVELOPMENT — positional identity
        development.step(dt);

        // 14. STEM CELLS — self-renewal vs differentiation
        stemCells.step(dt);

        // 15. PATHOGENS — infection progression
        pathogens.step(dt,
            immunity.innate.phagocytosis_rate + immunity.adaptive.cytotoxic_killing,
            immunity.innate.antiviral_state);

        // 16. IMMUNITY — immune response
        immunity.step(dt);

        // ════════════════════════════════════════════════════════════
        //  CROSS-TALK: feedback between modules
        // ════════════════════════════════════════════════════════════

        // Metabolism ← Signaling: growth signal increases ATP consumption
        metabolism.atpConsumptionRate = 0.3f + growthSignal * 0.5f + prolifSignal * 0.3f;

        // Signaling ← Metabolism: low ATP activates AMPK → inhibits mTOR
        if (energyCharge < 0.8f) {
            signaling.state.mTORC1_active *= energyCharge;
        }

        // Gene Regulation ← Cancer: mutations alter pathway parameters
        auto effects = cancer.getEffects();
        if (effects.ras_constitutive > 0) {
            signaling.state.Ras_GTP = fmaxf(signaling.state.Ras_GTP, effects.ras_constitutive);
        }
        if (effects.p53_functional < 1.0f) {
            dnaRepair.ddr.p53_level *= effects.p53_functional;
        }

        // Apoptosis ← Pathogen: viral infection can trigger/inhibit apoptosis
        if (pathogens.isInfected()) {
            if (pathogens.pathogen.apoptosis_inhibition > 0) {
                apoptosis.state.Bcl2 += pathogens.pathogen.apoptosis_inhibition * 0.1f * dt;
            }
        }

        // Immunity ← Pathogen: pathogen presence activates immune response
        if (pathogens.isInfected()) {
            immunity.innate.stress_ligands = pathogens.getPathogenLoad() * 0.1f;
            immunity.innate.MHC_I_surface = 1.0f - pathogens.pathogen.MHC_I_downregulation;
        }

        // UPR ← Proteome: misfolded proteins cause ER stress
        vesicles.upr.ER_stress = proteome.misfolded_fraction * 2.0f;

        // Autophagy ← Metabolism: starvation induces autophagy
        if (energyCharge < 0.7f) {
            vesicles.autophagy.mTORC1_inhibition = 1.0f - energyCharge;
        }
    }

    // ── Compatibility: get ATP as percentage (for existing Simulation.h) ──
    float getATPpercent() const { return metabolism.getATPpercent(); }
    float getStress() const { return 100.0f * (1.0f - metabolism.pool.energyCharge()); }
    float getROS() const { return metabolism.pool.NADH * 0.1f; } // Simplified ROS proxy
    bool isApoptotic() const { return apoptosis.isApoptotic(); }
    float getDamageLevel() const { return dnaRepair.getDamageLevel(); }
    float getMigrationSpeed() const { return cytoskeleton.migration_speed; }
    float getCellStiffness() const { return cytoskeleton.cell_stiffness; }

    void printStatus() const {
        printf("╔══ CellBiology Engine Status ══╗\n");
        metabolism.printStatus();
        signaling.printStatus();
        printf("  Genome: %d chromosomes, telomere avg=%.0f bp\n",
               (int)genome.chromosomes.size(), genome.averageTelomereLength());
        printf("  Damage: %.3f, p53=%.2f, p21=%.2f\n",
               dnaRepair.getDamageLevel(), dnaRepair.ddr.p53_P, dnaRepair.ddr.p21_level);
        printf("  Apoptosis: %s (progress=%.1f%%)\n",
               apoptosis.isApoptotic() ? "COMMITTED" : "alive",
               apoptosis.progress() * 100);
        printf("  Vesicles: %d active, ER stress=%.2f, autophagy=%.2f\n",
               vesicles.getVesicleCount(), vesicles.getERstress(), vesicles.getAutophagyLevel());
        printf("  Cytoskeleton: stiffness=%.2f, migration=%.1f µm/min\n",
               cytoskeleton.cell_stiffness, cytoskeleton.migration_speed);
        printf("  Cancer: %d mutations, malignancy=%.2f, transformed=%s\n",
               (int)cancer.mutations.size(), cancer.hallmarks.overallMalignancy(),
               cancer.is_transformed ? "YES" : "no");
        printf("  Infection: %s (load=%.2f)\n",
               pathogens.isInfected() ? "INFECTED" : "clean", pathogens.getPathogenLoad());
        printf("╚══════════════════════════════╝\n");
    }
};
