#pragma once
#include <simd/simd.h>
#include <string>
#include <vector>
#include <cstdint>

// ══════════════════════════════════════════════════════════════════════════
//  Virion.h — Individual virus particle in the medium (or inside a host).
//
//  Physics-first discipline (Plan Part 7 §44):
//    * No `tropism = "CD4"` field. A virion binds a host because its spike
//      descriptor vector is compatible with one of the host's surface
//      receptor descriptor vectors (BindingMatcher-style score).
//    * The virion is a physical entity with position, velocity, stage
//      (FREE → BOUND → ENTERING → UNCOATING → REPLICATING → ASSEMBLING →
//      BUDDING → SPILLED). Transitions are triggered by physical
//      conditions (bound-spike count, cytoplasmic pH, elapsed dwell).
//    * Geometry is generated procedurally from Caspar-Klug triangulation
//      (CapsidBuilder consumes T number, protomer id, spike count).
//
//  A virion carries no MOA flag. When it enters a cell, it competes for
//  ribosomes via existing translationRateMul; when enough copies assemble,
//  egress pressure rises; on host apoptosis the intracellular pool spills
//  into MediumField and re-infects neighbours.
// ══════════════════════════════════════════════════════════════════════════

enum class VirionShape : uint8_t {
    ICOSAHEDRAL = 0,   // most small-genome viruses
    HELICAL     = 1,   // TMV, rabies
    COMPLEX     = 2,   // pox, T4 bacteriophage
    PLEOMORPHIC = 3,   // enveloped, variable (flu, coronavirus)
};

enum class GenomeKind : uint8_t {
    SSRNA_POS = 0,   // positive-sense single-strand RNA
    SSRNA_NEG = 1,   // negative-sense
    DSRNA     = 2,
    SSDNA     = 3,
    DSDNA     = 4,
    RETRO     = 5,   // ssRNA + reverse-transcriptase
};

enum class VirionStage : uint8_t {
    FREE        = 0,   // drifting in medium
    BOUND       = 1,   // spike-receptor engaged on host surface
    ENTERING    = 2,   // fusion or endocytic uptake underway
    UNCOATING   = 3,   // capsid disassembly inside cytosol
    REPLICATING = 4,   // hijacking host machinery
    ASSEMBLING  = 5,   // new virions forming
    BUDDING     = 6,   // continuous export through host membrane
    SPILLED     = 7,   // released by lysis / apoptosis (visual only)
    CLEARED     = 8,   // engulfed by phagocyte or decayed
};

// Virion "species" — the catalogue entry loaded from data/pathogens/virions.yaml.
// One VirionSpec, many Virion instances.
struct VirionSpec {
    std::string  id;                  // "flu_h1n1", "hbv_01"
    std::string  displayName;
    VirionShape  shape = VirionShape::ICOSAHEDRAL;
    int          T_number = 1;        // Caspar-Klug triangulation
    bool         enveloped = false;
    float        radius_nm = 50.0f;   // outer radius for rendering + diffusion
    float        hydrodynamicRadius_nm = 55.0f; // Stokes-Einstein
    GenomeKind   genomeKind = GenomeKind::SSRNA_POS;
    int          genomeLength_nt = 10000;
    int          spikesPerVirion = 80;

    // Spike descriptor vector — used by BindingMatcher against host-cell
    // receptor descriptors. 5-dim Lipinski-like: logP, mw, hbd, hba, aromatic.
    // Spike protein = averaged pharmacophore over its binding face.
    float spike_logP = 2.0f;
    float spike_mw   = 60000.0f;      // Da for full spike protein
    int   spike_hbd  = 6;
    int   spike_hba  = 14;
    int   spike_aromatic = 3;

    // Intracellular kinetics
    float uncoatDwell_bioSec    = 1800.0f;    // 30 bio-min
    float replicationRate_per_s = 1.2e-2f;    // genome copies / bio-s
    float assemblyThreshold     = 200.0f;     // genomes before assembly starts
    float budRate_per_s         = 2.0e-3f;    // virions exported / bio-s
    float lyticYield            = 1000.0f;    // burst size at host death

    // Receptor affinity preferences — which host receptor types are
    // structurally compatible. Populated from YAML. Host-cell receptor
    // mosaic is scored per spike against this.
    std::vector<std::string> preferredReceptors;  // e.g. "PT_RECEPTOR_RTK"

    // Cytotoxicity — how much stress per replicating virion, accumulated
    // into the host's existing drug_pro_apop trigger slot.
    float cytotoxicity_per_copy = 5e-5f;
};

// Live virion instance.
struct Virion {
    int specIdx = -1;                 // index into VirionRegistry::specs
    simd_float3 position = {0,0,0};
    simd_float3 velocity = {0,0,0};
    VirionStage stage    = VirionStage::FREE;
    int hostCellIdx      = -1;        // -1 if free in medium
    float stageTimer_s   = 0.0f;

    // Intracellular state (valid when stage ∈ {REPLICATING, ASSEMBLING, BUDDING})
    float genomeCopies   = 0.0f;
    float assembled      = 0.0f;       // complete new virions ready to bud/burst
    float boundSpikes    = 0.0f;       // spike-receptor bindings on host surface

    // Lightweight ID so we can track specific virions across frames.
    uint32_t uid = 0;
};
