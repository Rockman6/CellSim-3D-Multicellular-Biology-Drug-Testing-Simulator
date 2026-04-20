#pragma once
#include "Virion.h"
#include "Bacterium.h"
#include <cstdio>
#include <cstring>

// ══════════════════════════════════════════════════════════════════════════
//  PathogenRegistry.h — catalogue of virion & bacterium species loaded
//  from data/pathogens/*.yaml. Indexed by id; Virion.specIdx /
//  Bacterium.specIdx point into the corresponding vector.
//
//  YAML is a flat CSV-like dialect (one block per species, key: value).
//  We keep the parser trivial — no libyaml dependency.
// ══════════════════════════════════════════════════════════════════════════

class PathogenRegistry {
public:
    std::vector<VirionSpec>     virionSpecs;
    std::vector<BacteriumSpec>  bacteriumSpecs;

    int findVirion(const std::string& id) const {
        for (int i = 0; i < (int)virionSpecs.size(); i++)
            if (virionSpecs[i].id == id) return i;
        return -1;
    }
    int findBacterium(const std::string& id) const {
        for (int i = 0; i < (int)bacteriumSpecs.size(); i++)
            if (bacteriumSpecs[i].id == id) return i;
        return -1;
    }

    // Minimal key:value block parser. A new entry starts when the line
    // begins with "- id:" (YAML list dash). Everything until the next
    // dash-id line is one entry. Comments (#) are ignored. The value is
    // trimmed; trailing inline comments are stripped.
    bool loadVirionsYaml(const char* path) {
        FILE* f = fopen(path, "r");
        if (!f) return false;
        char line[1024];
        VirionSpec cur;
        bool inEntry = false;
        auto commit = [&]() {
            if (inEntry && !cur.id.empty()) virionSpecs.push_back(cur);
            cur = VirionSpec();
            inEntry = false;
        };
        while (fgets(line, sizeof(line), f)) {
            std::string L = trimComment(line);
            if (L.empty()) continue;
            std::string key, val;
            bool dashStart = false;
            parseLine(L, key, val, dashStart);
            if (dashStart || key == "id") {
                if (key == "id" && dashStart) { commit(); cur.id = val; inEntry = true; continue; }
                if (key == "id") { commit(); cur.id = val; inEntry = true; continue; }
            }
            if (!inEntry) continue;
            if      (key == "displayName")           cur.displayName = val;
            else if (key == "shape")                  cur.shape = parseShape(val);
            else if (key == "T_number")               cur.T_number = atoi(val.c_str());
            else if (key == "enveloped")              cur.enveloped = (val == "true" || val == "1");
            else if (key == "radius_nm")              cur.radius_nm = atof(val.c_str());
            else if (key == "hydrodynamicRadius_nm")  cur.hydrodynamicRadius_nm = atof(val.c_str());
            else if (key == "genomeKind")             cur.genomeKind = parseGenomeKind(val);
            else if (key == "genomeLength_nt")        cur.genomeLength_nt = atoi(val.c_str());
            else if (key == "spikesPerVirion")        cur.spikesPerVirion = atoi(val.c_str());
            else if (key == "spike_logP")             cur.spike_logP = atof(val.c_str());
            else if (key == "spike_mw")               cur.spike_mw = atof(val.c_str());
            else if (key == "spike_hbd")              cur.spike_hbd = atoi(val.c_str());
            else if (key == "spike_hba")              cur.spike_hba = atoi(val.c_str());
            else if (key == "spike_aromatic")         cur.spike_aromatic = atoi(val.c_str());
            else if (key == "uncoatDwell_bioSec")     cur.uncoatDwell_bioSec = atof(val.c_str());
            else if (key == "replicationRate_per_s")  cur.replicationRate_per_s = atof(val.c_str());
            else if (key == "assemblyThreshold")      cur.assemblyThreshold = atof(val.c_str());
            else if (key == "budRate_per_s")          cur.budRate_per_s = atof(val.c_str());
            else if (key == "lyticYield")             cur.lyticYield = atof(val.c_str());
            else if (key == "cytotoxicity_per_copy")  cur.cytotoxicity_per_copy = atof(val.c_str());
            else if (key == "preferredReceptors")     splitList(val, cur.preferredReceptors);
        }
        commit();
        fclose(f);
        return true;
    }

    bool loadBacteriaYaml(const char* path) {
        FILE* f = fopen(path, "r");
        if (!f) return false;
        char line[1024];
        BacteriumSpec cur;
        bool inEntry = false;
        auto commit = [&]() {
            if (inEntry && !cur.id.empty()) bacteriumSpecs.push_back(cur);
            cur = BacteriumSpec();
            inEntry = false;
        };
        while (fgets(line, sizeof(line), f)) {
            std::string L = trimComment(line);
            if (L.empty()) continue;
            std::string key, val;
            bool dashStart = false;
            parseLine(L, key, val, dashStart);
            if (key == "id") { commit(); cur.id = val; inEntry = true; continue; }
            if (!inEntry) continue;
            if      (key == "displayName")         cur.displayName = val;
            else if (key == "shape")                cur.shape = parseBacShape(val);
            else if (key == "gram")                 cur.gram = parseGram(val);
            else if (key == "length_um")            cur.length_um = atof(val.c_str());
            else if (key == "width_um")             cur.width_um = atof(val.c_str());
            else if (key == "peptidoglycan_nm")     cur.peptidoglycan_nm = atof(val.c_str());
            else if (key == "has_outer_membrane")   cur.has_outer_membrane = (val == "true" || val == "1");
            else if (key == "has_lps_fringe")       cur.has_lps_fringe = (val == "true" || val == "1");
            else if (key == "has_capsule")          cur.has_capsule = (val == "true" || val == "1");
            else if (key == "flagella_count")       cur.flagella_count = atoi(val.c_str());
            else if (key == "pili_count")           cur.pili_count = atoi(val.c_str());
            else if (key == "adhesin_logP")         cur.adhesin_logP = atof(val.c_str());
            else if (key == "adhesin_mw")           cur.adhesin_mw = atof(val.c_str());
            else if (key == "adhesin_hbd")          cur.adhesin_hbd = atoi(val.c_str());
            else if (key == "adhesin_hba")          cur.adhesin_hba = atoi(val.c_str());
            else if (key == "adhesin_aromatic")     cur.adhesin_aromatic = atoi(val.c_str());
            else if (key == "doublingTime_bioSec")  cur.doublingTime_bioSec = atof(val.c_str());
            else if (key == "toxinRate_per_s")      cur.toxinRate_per_s = atof(val.c_str());
            else if (key == "toxinCytotoxicity")    cur.toxinCytotoxicity = atof(val.c_str());
            else if (key == "biomass_bm")           cur.biomass_bm = atof(val.c_str());
            else if (key == "preferredReceptors")   splitList(val, cur.preferredReceptors);
        }
        commit();
        fclose(f);
        return true;
    }

private:
    static std::string trimComment(const char* raw) {
        std::string s(raw);
        auto h = s.find('#');
        if (h != std::string::npos) s = s.substr(0, h);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r' || s.back() == ' ' || s.back() == '\t'))
            s.pop_back();
        size_t i = 0;
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) i++;
        return s.substr(i);
    }
    static void parseLine(const std::string& L, std::string& key, std::string& val, bool& dashStart) {
        std::string s = L;
        dashStart = false;
        if (!s.empty() && s[0] == '-') {
            dashStart = true;
            size_t i = 1;
            while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) i++;
            s = s.substr(i);
        }
        auto c = s.find(':');
        if (c == std::string::npos) { key.clear(); val.clear(); return; }
        key = s.substr(0, c);
        std::string v = s.substr(c + 1);
        size_t i = 0;
        while (i < v.size() && (v[i] == ' ' || v[i] == '\t' || v[i] == '"')) i++;
        v = v.substr(i);
        while (!v.empty() && (v.back() == ' ' || v.back() == '\t' || v.back() == '"')) v.pop_back();
        val = v;
        while (!key.empty() && (key.back() == ' ' || key.back() == '\t')) key.pop_back();
    }
    static void splitList(const std::string& s, std::vector<std::string>& out) {
        out.clear();
        std::string t = s;
        if (!t.empty() && t.front() == '[') t = t.substr(1);
        if (!t.empty() && t.back() == ']')  t.pop_back();
        size_t i = 0;
        while (i <= t.size()) {
            size_t j = t.find(',', i);
            if (j == std::string::npos) j = t.size();
            std::string tok = t.substr(i, j - i);
            size_t a = 0, b = tok.size();
            while (a < b && (tok[a] == ' ' || tok[a] == '\'' || tok[a] == '"')) a++;
            while (b > a && (tok[b-1] == ' ' || tok[b-1] == '\'' || tok[b-1] == '"')) b--;
            if (b > a) out.push_back(tok.substr(a, b - a));
            i = j + 1;
        }
    }
    static VirionShape parseShape(const std::string& v) {
        if (v == "HELICAL")     return VirionShape::HELICAL;
        if (v == "COMPLEX")     return VirionShape::COMPLEX;
        if (v == "PLEOMORPHIC") return VirionShape::PLEOMORPHIC;
        return VirionShape::ICOSAHEDRAL;
    }
    static GenomeKind parseGenomeKind(const std::string& v) {
        if (v == "SSRNA_NEG") return GenomeKind::SSRNA_NEG;
        if (v == "DSRNA")     return GenomeKind::DSRNA;
        if (v == "SSDNA")     return GenomeKind::SSDNA;
        if (v == "DSDNA")     return GenomeKind::DSDNA;
        if (v == "RETRO")     return GenomeKind::RETRO;
        return GenomeKind::SSRNA_POS;
    }
    static BacteriumShape parseBacShape(const std::string& v) {
        if (v == "COCCUS")     return BacteriumShape::COCCUS;
        if (v == "VIBRIO")     return BacteriumShape::VIBRIO;
        if (v == "SPIRILLUM")  return BacteriumShape::SPIRILLUM;
        if (v == "SPIROCHETE") return BacteriumShape::SPIROCHETE;
        return BacteriumShape::ROD;
    }
    static GramType parseGram(const std::string& v) {
        if (v == "GRAM_POSITIVE") return GramType::GRAM_POSITIVE;
        if (v == "ACID_FAST")     return GramType::ACID_FAST;
        return GramType::GRAM_NEGATIVE;
    }
};

inline PathogenRegistry& gPathogens() {
    static PathogenRegistry reg;
    return reg;
}
