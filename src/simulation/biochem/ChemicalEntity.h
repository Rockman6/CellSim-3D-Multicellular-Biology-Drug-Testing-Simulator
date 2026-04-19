#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Forward declarations to avoid heavy includes.
struct MoleculeData;
struct ProteinData;

// ══════════════════════════════════════════════════════════════════════════
//  ChemicalEntity.h — Universal "anything with a structure" record.
//
//  Every object in the simulator that carries a chemical structure —
//  metabolites (ATP, glucose), organelles / enzymes (ribosome, hexokinase),
//  receptors, drugs, viruses, bacteria, membrane lipids — subclasses
//  this single type. The physics / reaction / binding engines dispatch
//  on `ChemicalEntity` without branching on kind; kind only selects
//  which fields are populated.
//
//  See docs/CHEMICAL_ENTITY.md (plan Part Five §24) for the design
//  rationale. This header is a declaration-only skeleton — fields are
//  filled in by Tier 0-2 Python workers writing JSON to
//  data/bioagents/chem/<id>/ and loaded at simulator init.
// ══════════════════════════════════════════════════════════════════════════

enum class ChemicalEntityKind : uint8_t {
    UNKNOWN         = 0,
    METABOLITE      = 1,   // ATP, glucose, water, ions, AA pool, …
    DRUG            = 2,   // small-molecule therapeutic
    VIRUS           = 3,   // nucleocapsid + optional envelope + genome
    BACTERIUM       = 4,   // whole-organism agent
    ORGANELLE       = 5,   // membrane-bounded cellular compartment
    ENZYME          = 6,   // catalytic protein (may be an organelle component)
    RECEPTOR        = 7,   // membrane-embedded signalling protein
    MEMBRANE_LIPID  = 8,   // PC/PE/PS/cholesterol/sphingolipid
    NUCLEIC_ACID    = 9,   // DNA / RNA (chromosomes, mRNA, viral genomes)
    COFACTOR        = 10,  // NAD+/NADH, FAD/FADH2, CoA, heme, ATP-as-cofactor
    PEPTIDE         = 11,  // short peptides / hormones
    ANTIBODY        = 12,  // IgG / IgA / IgM class immunoglobulin
    ION             = 13,  // Na+, K+, Cl−, Ca²⁺, Mg²⁺, H+
};

// Geometric + chemical pharmacophore point used by BindingMatcher.
// Populated from either (a) a PDB binding-pocket analysis for targets,
// or (b) RDKit feature-map detection for drugs.
struct PharmacophoreFeature {
    enum Kind : uint8_t {
        HBOND_DONOR     = 0,
        HBOND_ACCEPTOR  = 1,
        HYDROPHOBIC     = 2,
        AROMATIC        = 3,
        POS_CHARGE      = 4,
        NEG_CHARGE      = 5,
        METAL_CHELATOR  = 6,
    };
    Kind  kind;
    float x, y, z;     // Å relative to molecule centroid
    float radius;      // tolerance Å
};

// Binding kinetics between this entity and another (both identified by
// ChemicalEntity.id). Populated by BindingMatcher (Tier 0 descriptor +
// fingerprint) and progressively refined by Tier 3 docking / Tier 4 MD.
struct BindingAffinity {
    float Kd_mM        = 1e6f;   // dissociation constant
    float k_on_per_uM  = 0.0f;   // on-rate (1/(µM·s))
    float k_off_per_s  = 0.0f;   // off-rate (1/s)
    float score        = 0.0f;   // 0..1 summary shown in UI
    uint8_t tier       = 0;      // highest tier that produced this value
};

// How a bioagent enters a host cell. Irrelevant for intracellular
// ChemicalEntities (they're already inside) — default is NONE.
enum class EntryPathway : uint8_t {
    NONE                 = 0,
    PASSIVE_DIFFUSION    = 1,   // small lipophilic drugs
    RECEPTOR_ENDOCYTOSIS = 2,   // virus spike → host receptor
    MEMBRANE_FUSION      = 3,   // enveloped virus → host membrane
    PHAGOCYTOSIS         = 4,   // bacterium → phagocyte
    INVASION             = 5,   // pathogen actively invading
};

// The universal chemistry record. Tier 0 populates descriptor + FP;
// Tier 1 adds partial charges & LJ; Tier 1b adds polarizability;
// Tier 2 adds HOMO/LUMO; Tier 3 fills affinities; Tier 5 fills
// reactsIn / catalyses. Empty fields are free — default-constructed.
struct ChemicalEntity {
    // ── Identity ──
    std::string        id;       // filesystem-safe key
    std::string        name;     // display name
    ChemicalEntityKind kind = ChemicalEntityKind::UNKNOWN;

    // ── Structure ──
    const MoleculeData* mol     = nullptr;   // ball-and-stick geometry
    const ProteinData*  protein = nullptr;   // macromolecular geometry
    std::vector<PharmacophoreFeature> pharmacophores;

    // ── Chemistry (Tier 0-2 descriptors) ──
    float mw                = 0.0f;   // Da
    float logP              = 0.0f;
    float logS              = 0.0f;
    float tpsa              = 0.0f;   // Å² polar surface area
    int   hbd               = 0;
    int   hba               = 0;
    int   rotatable_bonds   = 0;
    int   aromatic_rings    = 0;
    int   formal_charge     = 0;
    float polarizability    = 0.0f;   // Å³ (xTB Tier 2)
    float homo_eV           = 0.0f;
    float lumo_eV           = 0.0f;
    float dipole_debye      = 0.0f;

    // ── Tier 1 force-field parameters (per-atom, same length as mol->atoms) ──
    std::vector<float> partialCharges;
    std::vector<float> ljSigma_nm;
    std::vector<float> ljEps_kJ_per_mol;

    // ── Tier 1b polarizable multipole (AMOEBA, optional) ──
    std::vector<float> inducedDipoleMag;

    // ── Morgan ECFP4 fingerprint, 2048 bits packed into 256 bytes ──
    std::array<uint8_t, 256> morganFp{};

    // ── Binding affinities to other entities, keyed by other.id ──
    std::vector<std::pair<std::string, BindingAffinity>> affinities;

    // ── Reactions in which this entity participates ──
    std::vector<std::string> reactsIn;    // rule ids where this is substrate
    std::vector<std::string> catalyses;   // rule ids where this is enzyme

    // ── Bioagent-specific (null / empty for intracellular entities) ──
    EntryPathway entry = EntryPathway::NONE;
    std::vector<std::string> cargo;             // viral genome ids, toxin ids
    std::vector<std::string> modulatesTargets;  // shortcut list of TGT_* ids

    bool valid() const { return !id.empty(); }
};
