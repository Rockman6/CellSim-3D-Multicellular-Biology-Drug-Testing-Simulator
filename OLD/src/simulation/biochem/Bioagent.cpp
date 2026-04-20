#include "Bioagent.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cstring>

// Canonical global registry. Populated at sim init (Phase W).
BioagentRegistry gBioagents;

// Forward declare the target-registration helper.
void registerBuiltinTargets();

ChemicalEntity* BioagentRegistry::create(const std::string& id,
                                          ChemicalEntityKind kind) {
    auto it = indexById_.find(id);
    if (it != indexById_.end()) return &entities_[it->second];
    ChemicalEntity e;
    e.id = id;
    e.kind = kind;
    int idx = (int)entities_.size();
    entities_.push_back(std::move(e));
    indexById_[id] = idx;
    return &entities_[idx];
}

const ChemicalEntity* BioagentRegistry::get(const std::string& id) const {
    auto it = indexById_.find(id);
    return (it == indexById_.end()) ? nullptr : &entities_[it->second];
}

ChemicalEntity* BioagentRegistry::getMutable(const std::string& id) {
    auto it = indexById_.find(id);
    return (it == indexById_.end()) ? nullptr : &entities_[it->second];
}

BindingAffinity BioagentRegistry::affinity(const std::string& aId,
                                            const std::string& bId) const {
    const ChemicalEntity* a = get(aId);
    if (!a) return {};
    for (const auto& [other, aff] : a->affinities) {
        if (other == bId) return aff;
    }
    // Try reverse direction — affinities are symmetric.
    const ChemicalEntity* b = get(bId);
    if (!b) return {};
    for (const auto& [other, aff] : b->affinities) {
        if (other == aId) return aff;
    }
    return {};  // default — Kd = 1 M, tier 0
}

// Minimal single-value extractor for a JSON scalar field.
static float jsonFloat(const std::string& body, const std::string& key,
                        float fallback = 0.0f) {
    size_t p = body.find("\"" + key + "\"");
    if (p == std::string::npos) return fallback;
    p = body.find(':', p);
    if (p == std::string::npos) return fallback;
    return (float)atof(body.c_str() + p + 1);
}
static int jsonInt(const std::string& body, const std::string& key,
                    int fallback = 0) {
    return (int)jsonFloat(body, key, (float)fallback);
}
static std::string jsonString(const std::string& body, const std::string& key) {
    size_t p = body.find("\"" + key + "\"");
    if (p == std::string::npos) return "";
    p = body.find(':', p);
    if (p == std::string::npos) return "";
    p = body.find('"', p);
    if (p == std::string::npos) return "";
    size_t q = body.find('"', p + 1);
    if (q == std::string::npos) return "";
    return body.substr(p + 1, q - p - 1);
}

static std::string readFile(const std::string& path) {
    std::ifstream f(path);
    if (!f) return {};
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}

static std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == ',') { out.push_back(cur); cur.clear(); }
        else         { cur += c; }
    }
    out.push_back(cur);
    return out;
}

int BioagentRegistry::loadFromDisk(const std::string& dataDir) {
    // Register the 13 built-in targets before processing drug CSV so
    // binding affinities can be computed immediately.
    registerBuiltinTargets();

    int loaded = 0;
    // ── drugs.csv ─────────────────────────────────────────────────
    std::ifstream df(dataDir + "/bioagents/drugs.csv");
    std::string line;
    while (df && std::getline(df, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("id,", 0) == 0) continue;
        auto cols = splitCSV(line);
        if (cols.size() < 4) continue;
        const std::string id     = cols[0];
        const std::string name   = cols[1];
        const std::string kind   = cols[2];
        const std::string smiles = cols[3];

        ChemicalEntity* e = create(id, ChemicalEntityKind::DRUG);
        e->name   = name;
        e->entry  = EntryPathway::PASSIVE_DIFFUSION;

        // Pull Tier-0 descriptor cache if present.
        std::string tier0 = readFile(dataDir + "/bioagents/chem/"
                                     + id + "/tier0_descriptors.json");
        if (!tier0.empty()) {
            e->mw              = jsonFloat(tier0, "mw");
            e->logP            = jsonFloat(tier0, "logP");
            e->tpsa            = jsonFloat(tier0, "tpsa");
            e->hbd             = jsonInt  (tier0, "hbd");
            e->hba             = jsonInt  (tier0, "hba");
            e->rotatable_bonds = jsonInt  (tier0, "rotatable_bonds");
            e->aromatic_rings  = jsonInt  (tier0, "aromatic_rings");
            e->formal_charge   = jsonInt  (tier0, "formal_charge");
        }
        loaded++;
    }
    printf("[BioagentRegistry] Loaded %d drugs from disk.\n", loaded);
    return loaded;
}
