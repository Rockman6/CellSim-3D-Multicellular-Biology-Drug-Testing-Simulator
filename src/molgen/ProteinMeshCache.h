#pragma once

#include "PDBParser.h"
#include <map>
#include <string>

// ══════════════════════════════════════════════════════════════════════════
//  ProteinMeshCache — Load real protein backbone structures from PDB files
//  Each protein is stored as a MoleculeData (atoms + bonds) that can be
//  rendered as backbone trace spheres/sticks at any position and scale.
// ══════════════════════════════════════════════════════════════════════════

struct ProteinEntry {
    const char* filename;   // e.g., "tRNA-Phe.pdb"
    const char* name;       // Display name
    float defaultScale;     // Scale factor when rendering inside cell
    float colorR, colorG, colorB;  // Base tint color
};

static const ProteinEntry PROTEIN_ENTRIES[] = {
    {"tRNA-Phe.pdb",        "tRNA",           0.003f,  0.3f, 0.8f, 0.4f},
    {"ribosome_80S.pdb",    "Ribosome 80S",   0.0008f, 0.55f, 0.45f, 0.8f},
    {"rna_pol_II.pdb",      "RNA Pol II",     0.0012f, 0.6f, 0.2f, 0.8f},
    {"spliceosome.pdb",     "Spliceosome",    0.0006f, 0.8f, 0.6f, 0.2f},
    {"hemoglobin.pdb",      "Hemoglobin",     0.004f,  0.9f, 0.3f, 0.3f},
    {"hbb_folded.pdb",      "Hemoglobin beta",0.0035f, 0.88f, 0.35f, 0.35f},
    {"hsp70_chaperone.pdb", "Hsp70 Chaperone",0.001f,  0.7f, 0.5f, 0.3f},
    {"dna_polymerase.pdb",  "DNA Polymerase", 0.003f,  0.3f, 0.7f, 0.9f},
    {"atp_synthase.pdb",    "ATP Synthase",   0.0008f, 0.9f, 0.8f, 0.2f},
};
static const int PROTEIN_ENTRY_COUNT = sizeof(PROTEIN_ENTRIES) / sizeof(ProteinEntry);

class ProteinMeshCache {
public:
    void init(const std::string& proteinDir) {
        proteinDir_ = proteinDir;
        loadAll();
    }

    const MoleculeData* get(const std::string& filename) const {
        auto it = cache_.find(filename);
        if (it != cache_.end() && it->second.valid()) return &it->second;
        return nullptr;
    }

    const ProteinEntry* getEntry(const std::string& filename) const {
        for (int i = 0; i < PROTEIN_ENTRY_COUNT; i++) {
            if (filename == PROTEIN_ENTRIES[i].filename) return &PROTEIN_ENTRIES[i];
        }
        return nullptr;
    }

    int loadedCount() const {
        int n = 0;
        for (auto& [k, v] : cache_) if (v.valid()) n++;
        return n;
    }

    std::vector<std::string> getLoadedFiles() const {
        std::vector<std::string> files;
        for (auto& [k, v] : cache_) if (v.valid()) files.push_back(k);
        return files;
    }

private:
    void loadAll() {
        for (int i = 0; i < PROTEIN_ENTRY_COUNT; i++) {
            std::string path = proteinDir_ + "/" + PROTEIN_ENTRIES[i].filename;
            MoleculeData mol = PDBParser::parsePDB(path);
            if (mol.valid()) {
                cache_[PROTEIN_ENTRIES[i].filename] = std::move(mol);
            }
        }
        printf("[ProteinMeshCache] Loaded %d protein structures\n", loadedCount());
    }

    std::string proteinDir_;
    std::map<std::string, MoleculeData> cache_;
};
