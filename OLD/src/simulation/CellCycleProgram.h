#pragma once

// ══════════════════════════════════════════════════════════════════════════
//  CellCycleProgram.h — per-cell biological machinery container
//
//  Every SimCell owns one of these. Before this struct existed the simulator
//  had a single global CentralDogmaState and a single global MitosisState
//  tied to cells[0], so only the focused cell had real replication and a
//  real mitosis state machine. Background cells divided on a wall-clock
//  timer with no notion of "my genome is replicated" or "my SAC is clear."
//
//  CellCycleProgram bundles the per-cell replication + mitosis state so
//  every lineage can run the same machinery. The render pipeline still
//  focuses on one chosen cell — it just dereferences the focused SimCell's
//  program instead of a privileged global.
//
//  Memory: each instance is a few KB (HBB_LENGTH=763 bool array dominates).
//  At MAX_CELLS=2000 this is ~few MB total, well within budget.
//
//  References cited in the division plan at
//  ~/plans/make-this-plan-into-virtual-key.md:
//    Alberts MBoC 7e Ch 4-5, 17-18
//    Diffley 2011 Cold Spring Harb (origin licensing)
//    Musacchio 2015 Curr Biol (SAC)
//    Bartek 2004 Nat Rev MCB (replication stress)
//    Kunkel 2004 JBC (post-MMR error rate)
//    Albeck 2008 PLoS Biol (apoptosis kinetics)
// ══════════════════════════════════════════════════════════════════════════

#include "CentralDogma.h"
#include <array>
#include <cstdio>
#include <cstring>

struct CellCycleAuditEntry {
    float time = 0.0f;
    int cellPhase = 0;
    MitosisPhase mitosisPhase = MITO_NONE;
    float replicationProgress = 0.0f;
    float replicationQuality = 1.0f;
    char event[24] = {0};
    char reason[64] = {0};
};

struct CellCycleProgram {
    // Replication / transcription / translation — existing struct reused.
    // For background cells only the replication half needs to run; the
    // transcription/translation portion idles until the cell is focused.
    CentralDogmaState cdogma;

    // Per-cell mitosis state machine (chromatin, spindle, SAC, cytokinesis).
    MitosisState mitosis;

    // Tetraploid flag — set if cytokinesis failed and the cell re-entered
    // G1 with 4N DNA. Next cycle should trigger p53-dependent arrest.
    bool tetraploid = false;

    // Whether cdogma has been initialized (codingMRNA built, origins licensed).
    // Starts false for background cells so we can skip their init until
    // they need to run replication or become focused.
    bool cdogmaInitialized = false;

    // Lightweight in-memory audit ring. This gives each cell its own short
    // checkpoint / mitosis history even before we promote the data to a
    // dedicated on-disk per-cell timeline.
    std::array<CellCycleAuditEntry, 64> audit{};
    int auditHead = 0;

    void ensureCDogmaInitialized() {
        if (cdogmaInitialized) return;
        cdogma.init();
        cdogmaInitialized = true;
    }

    void logEvent(float t, int cellPhase, const char* event, const char* reason = "") {
        CellCycleAuditEntry& entry = audit[(size_t)auditHead % audit.size()];
        entry.time = t;
        entry.cellPhase = cellPhase;
        entry.mitosisPhase = mitosis.phase;
        entry.replicationProgress = cdogma.replicationProgress;
        entry.replicationQuality = cdogma.replicationQuality;
        std::snprintf(entry.event, sizeof(entry.event), "%s", event ? event : "");
        std::snprintf(entry.reason, sizeof(entry.reason), "%s", reason ? reason : "");
        auditHead = (auditHead + 1) % (int)audit.size();
    }

    // Reset everything to a freshly post-mitotic state (for use in
    // deriveDaughter after a successful division).
    void resetForNewCycle() {
        // Only reset replication program — keep transcription state intact
        // so the focused-cell teaching animations stay continuous.
        cdogma.resetReplicationProgram();
        // Clear mitosis state
        mitosis.active = false;
        mitosis.phase = MITO_NONE;
        mitosis.phaseTimer = 0;
        mitosis.totalProgress = 0;
        mitosis.chromatinCondensation = 0;
        mitosis.nucleusAlpha = 1.0f;
        mitosis.nuclearEnvelopeBreakdown = 0.0f;
        mitosis.poleSeparation = 0;
        mitosis.spindleAssembly = 0.0f;
        mitosis.spindleDisassembly = 0.0f;
        mitosis.kinetochoreAttachment = 0.0f;
        mitosis.metaphaseAlignment = 0.0f;
        mitosis.spindleCheckpoint = 0.0f;
        mitosis.chromatidSeparation = 0.0f;
        mitosis.checkpointHold = 0.0f;
        mitosis.furrowDepth = 0;
        mitosis.cellElongation = 0;
        mitosis.contractileRingAssembly = 0.0f;
        mitosis.organelleDuplication = 0;
        mitosis.organelleMigration = 0;
        mitosis.particlesDuplicated = false;
        mitosis.dnaCheckpointPassed = false;
        mitosis.dnaDuplicationProgress = 0.0f;
        mitosis.postDivisionTimer = 0;
        mitosis.fadeAlpha = 1.0f;
        tetraploid = false;
    }
};
