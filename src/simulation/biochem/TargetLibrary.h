#pragma once
#include "BindingMatcher.h"
#include <string>
#include <vector>
#include <unordered_map>

// Forward declaration to avoid heavy include.
struct SimCell;

// ══════════════════════════════════════════════════════════════════════════
//  TargetLibrary.h — Registry of druggable cellular components + their
//  function modulators.
//
//  A "target" is any molecular component in the cell that a drug can
//  bind. When occupancy > 0, the modulator callback adjusts the cell's
//  existing biology fields. The adjustments only ever change RATES or
//  POOLS the existing code already reads — never add new outcome logic.
//  This is what keeps the MOA-free discipline: we never write
//  "drug kills cell"; we write "drug bound to CDK1 → CDK1 inactive →
//  existing cycle gate fails → existing apoptosis engine triggers."
// ══════════════════════════════════════════════════════════════════════════

struct TargetEntry {
    TargetProfile profile;             // preferred drug descriptor ranges
    std::string   kind;                // enzyme / structural / receptor / ...
    int           baseCopyCount = 100000;
    // Modulator: (cell, occupancy ∈ [0..1], dt_biosec) → void
    void (*modulator)(SimCell& c, float occ, float dt_biosec) = nullptr;
};

class TargetLibrary {
public:
    void registerTarget(const std::string& id,
                        const std::string& kind,
                        const TargetProfile& profile,
                        void (*modulator)(SimCell&, float, float)) {
        int idx = (int)entries_.size();
        TargetEntry e; e.profile = profile; e.kind = kind;
        e.modulator = modulator;
        entries_.push_back(e);
        index_[id] = idx;
        ids_.push_back(id);
    }

    int count() const { return (int)entries_.size(); }
    const std::string& idAt(int i) const { return ids_[i]; }
    const TargetEntry& at(int i) const { return entries_[i]; }
    int indexOf(const std::string& id) const {
        auto it = index_.find(id);
        return (it == index_.end()) ? -1 : it->second;
    }

private:
    std::vector<TargetEntry>              entries_;
    std::vector<std::string>              ids_;
    std::unordered_map<std::string, int>  index_;
};

extern TargetLibrary gTargets;
