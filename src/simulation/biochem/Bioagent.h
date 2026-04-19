#pragma once

#include "ChemicalEntity.h"
#include <string>
#include <unordered_map>
#include <vector>

// ══════════════════════════════════════════════════════════════════════════
//  Bioagent.h — Registry + lookup for every ChemicalEntity in the sim.
//
//  A single `BioagentRegistry` holds the entire universe of
//  ChemicalEntities the simulator knows about: metabolites (ATP, glucose,
//  …), organelles (ribosome, hexokinase, …), drugs, viruses, bacteria,
//  receptors, lipids. Registration is a one-shot at simulator init,
//  populated from data/bioagents/*.csv and data/bioagents/chem/*.json.
//
//  Access is by string id (e.g., "atp", "cisplatin", "TGT_CDK1"). Sim
//  hot-loop code reads ChemicalEntity pointers out of the registry and
//  dispatches physics on them without branching on kind.
//
//  This is intentionally a declaration-only skeleton (Phase M1-M3).
//  Loading logic arrives in Phase W; consumers arrive in Phase X.
// ══════════════════════════════════════════════════════════════════════════

class BioagentRegistry {
public:
    // Load a CSV row into a new entity. No-op if id already exists.
    ChemicalEntity* create(const std::string& id, ChemicalEntityKind kind);

    // Return nullptr if not found.
    const ChemicalEntity* get(const std::string& id) const;
    ChemicalEntity*       getMutable(const std::string& id);

    // Iterate all entities, optionally filtered by kind.
    const std::vector<ChemicalEntity>& all() const { return entities_; }

    // Resolve a binding affinity between two entity ids — returns the
    // highest-tier cached value, or a default (Kd = 1 M, tier 0) if none.
    BindingAffinity affinity(const std::string& aId,
                             const std::string& bId) const;

    int count() const { return (int)entities_.size(); }

    // Phase-W population: read data/bioagents/*.csv + per-entity JSON
    // caches and build the registry. Currently a stub — returns 0.
    int loadFromDisk(const std::string& dataDir);

private:
    std::vector<ChemicalEntity>             entities_;
    std::unordered_map<std::string, int>    indexById_;
};

// The canonical global registry. Owned by main.mm / Simulation init.
extern BioagentRegistry gBioagents;
