#pragma once

#include "PDBParser.h"
#include "MoleculeLibrary.h"
#include <dirent.h>
#include <map>
#include <string>

// ══════════════════════════════════════════════════════════════════════════
//  MoleculeCache — Load and cache molecular structures from disk
//  Checks assets/molecules/ for pre-generated .sdf/.pdb files
// ══════════════════════════════════════════════════════════════════════════

class MoleculeCache {
public:
    void init(const std::string& moleculeDir) {
        moleculeDir_ = moleculeDir;
        printf("[MoleculeCache] Directory: %s\n", moleculeDir.c_str());
        loadAll();
        // Also scan any drug SDFs produced by the chemistry worker
        // (data/bioagents/chem/<id>/ → assets/drugs/<id>.sdf). These
        // entries aren't in MOLECULE_ENTRIES but still need to be
        // renderable once applyDrug() spawns particles for them.
        std::string drugDir = moleculeDir + "/../drugs";
        loadDirectoryGlob(drugDir);
    }

    // Get a cached molecule by ID (e.g., "atp", "glucose")
    const MoleculeData* get(const std::string& id) const {
        auto it = cache_.find(id);
        if (it != cache_.end() && it->second.valid()) return &it->second;
        return nullptr;
    }

    // Get all loaded molecule IDs
    std::vector<std::string> getLoadedIds() const {
        std::vector<std::string> ids;
        for (auto& [id, mol] : cache_) {
            if (mol.valid()) ids.push_back(id);
        }
        return ids;
    }

    int loadedCount() const {
        int n = 0;
        for (auto& [id, mol] : cache_) if (mol.valid()) n++;
        return n;
    }

private:
    // Scan a directory for *.sdf files and index each by its stem.
    // Skips entries already in the cache (so the MOLECULE_ENTRIES
    // registry remains the primary source of truth).
    void loadDirectoryGlob(const std::string& dir) {
        DIR* d = opendir(dir.c_str());
        if (!d) return;
        int added = 0;
        struct dirent* entry;
        while ((entry = readdir(d)) != nullptr) {
            std::string name = entry->d_name;
            if (name.size() < 5) continue;
            if (name.substr(name.size() - 4) != ".sdf") continue;
            std::string id = name.substr(0, name.size() - 4);
            if (cache_.count(id)) continue;  // already loaded
            std::string path = dir + "/" + name;
            MoleculeData mol = PDBParser::parseSDF(path);
            if (mol.valid()) {
                cache_[id] = std::move(mol);
                added++;
            }
        }
        closedir(d);
        if (added > 0) {
            printf("[MoleculeCache] +%d drug SDFs from %s\n", added, dir.c_str());
        }
    }

    void loadAll() {
        for (int i = 0; i < MOLECULE_ENTRY_COUNT; i++) {
            const auto& entry = MOLECULE_ENTRIES[i];
            std::string id = entry.id;

            // Try .sdf first, then .pdb
            std::string sdfPath = moleculeDir_ + "/" + id + ".sdf";
            std::string pdbPath = moleculeDir_ + "/" + id + ".pdb";

            MoleculeData mol;
            FILE* f = fopen(sdfPath.c_str(), "r");
            if (f) {
                fclose(f);
                mol = PDBParser::parseSDF(sdfPath);
            } else {
                f = fopen(pdbPath.c_str(), "r");
                if (f) {
                    fclose(f);
                    mol = PDBParser::parsePDB(pdbPath);
                }
            }

            if (mol.valid()) {
                cache_[id] = std::move(mol);
            }
        }
        printf("[MoleculeCache] Loaded %d molecules\n", loadedCount());
    }

    std::string moleculeDir_;
    std::map<std::string, MoleculeData> cache_;
};
