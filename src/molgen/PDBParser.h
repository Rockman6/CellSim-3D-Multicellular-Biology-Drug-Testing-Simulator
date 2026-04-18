#pragma once

#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <simd/simd.h>

// ══════════════════════════════════════════════════════════════════════════
//  PDBParser — Parse PDB and SDF files into atom positions for rendering
//  Produces atom positions + element types for ball-and-stick visualization
// ══════════════════════════════════════════════════════════════════════════

struct Atom {
    simd_float3 position;
    int         element;    // Atomic number
    float       radius;     // Van der Waals radius (Angstroms)
    float       color[3];   // CPK coloring
};

struct Bond {
    int atomA, atomB;
    int order;  // 1=single, 2=double, 3=triple
};

struct MoleculeData {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    simd_float3 center;     // Centroid
    float       maxExtent;  // Max distance from center

    bool valid() const { return !atoms.empty(); }

    void computeBounds() {
        if (atoms.empty()) { center = {0,0,0}; maxExtent = 1; return; }
        simd_float3 sum = {0,0,0};
        for (auto& a : atoms) {
            sum.x += a.position.x;
            sum.y += a.position.y;
            sum.z += a.position.z;
        }
        float n = (float)atoms.size();
        center = {sum.x/n, sum.y/n, sum.z/n};
        maxExtent = 0;
        for (auto& a : atoms) {
            float dx = a.position.x - center.x;
            float dy = a.position.y - center.y;
            float dz = a.position.z - center.z;
            float d = sqrtf(dx*dx + dy*dy + dz*dz);
            if (d > maxExtent) maxExtent = d;
        }
        if (maxExtent < 0.1f) maxExtent = 1.0f;
    }
};

// CPK coloring and Van der Waals radii by element
static void getElementProperties(const char* symbol, int& atomicNum, float& vdwRadius, float color[3]) {
    // Default
    atomicNum = 6; vdwRadius = 1.7f;
    color[0] = 0.5f; color[1] = 0.5f; color[2] = 0.5f;

    if (!symbol) return;
    // Skip leading whitespace
    while (*symbol == ' ') symbol++;

    if (strncmp(symbol, "C", 1) == 0 && (symbol[1] == 0 || symbol[1] == ' ')) {
        atomicNum = 6; vdwRadius = 1.70f;
        color[0] = 0.3f; color[1] = 0.3f; color[2] = 0.3f;
    } else if (strncmp(symbol, "N", 1) == 0 && (symbol[1] == 0 || symbol[1] == ' ')) {
        atomicNum = 7; vdwRadius = 1.55f;
        color[0] = 0.2f; color[1] = 0.3f; color[2] = 0.9f;
    } else if (strncmp(symbol, "O", 1) == 0 && (symbol[1] == 0 || symbol[1] == ' ')) {
        atomicNum = 8; vdwRadius = 1.52f;
        color[0] = 0.9f; color[1] = 0.1f; color[2] = 0.1f;
    } else if (strncmp(symbol, "H", 1) == 0 && (symbol[1] == 0 || symbol[1] == ' ')) {
        atomicNum = 1; vdwRadius = 1.20f;
        color[0] = 0.9f; color[1] = 0.9f; color[2] = 0.9f;
    } else if (strncmp(symbol, "S", 1) == 0 && (symbol[1] == 0 || symbol[1] == ' ')) {
        atomicNum = 16; vdwRadius = 1.80f;
        color[0] = 0.9f; color[1] = 0.8f; color[2] = 0.2f;
    } else if (strncmp(symbol, "P", 1) == 0 && (symbol[1] == 0 || symbol[1] == ' ')) {
        atomicNum = 15; vdwRadius = 1.80f;
        color[0] = 0.9f; color[1] = 0.5f; color[2] = 0.0f;
    } else if (strncmp(symbol, "Fe", 2) == 0) {
        atomicNum = 26; vdwRadius = 2.00f;
        color[0] = 0.7f; color[1] = 0.4f; color[2] = 0.1f;
    } else if (strncmp(symbol, "Ca", 2) == 0) {
        atomicNum = 20; vdwRadius = 2.31f;
        color[0] = 0.2f; color[1] = 0.8f; color[2] = 0.2f;
    } else if (strncmp(symbol, "Mg", 2) == 0) {
        atomicNum = 12; vdwRadius = 1.73f;
        color[0] = 0.1f; color[1] = 0.7f; color[2] = 0.1f;
    } else if (strncmp(symbol, "Cl", 2) == 0) {
        atomicNum = 17; vdwRadius = 1.75f;
        color[0] = 0.1f; color[1] = 0.9f; color[2] = 0.1f;
    } else if (strncmp(symbol, "Na", 2) == 0) {
        atomicNum = 11; vdwRadius = 2.27f;
        color[0] = 0.6f; color[1] = 0.3f; color[2] = 0.9f;
    } else if (strncmp(symbol, "Zn", 2) == 0) {
        atomicNum = 30; vdwRadius = 1.39f;
        color[0] = 0.5f; color[1] = 0.5f; color[2] = 0.7f;
    }
}

// Parse element symbol from PDB atom line (columns 77-78)
static void parseElementFromPDB(const char* line, char elem[3]) {
    elem[0] = elem[1] = elem[2] = 0;
    int len = (int)strlen(line);
    if (len >= 78) {
        elem[0] = line[76];
        elem[1] = line[77];
        elem[2] = 0;
        // Trim
        if (elem[0] == ' ') { elem[0] = elem[1]; elem[1] = 0; }
        if (elem[1] == ' ') elem[1] = 0;
    } else {
        // Fallback: extract from atom name (columns 13-16)
        if (len >= 16) {
            char name[5] = {line[12], line[13], line[14], line[15], 0};
            // First non-space character is usually the element
            for (int i = 0; i < 4; i++) {
                if (name[i] != ' ') {
                    elem[0] = name[i];
                    if (i+1 < 4 && name[i+1] >= 'a' && name[i+1] <= 'z')
                        elem[1] = name[i+1];
                    break;
                }
            }
        }
    }
}

class PDBParser {
public:
    // Parse a PDB file
    static MoleculeData parsePDB(const std::string& path) {
        MoleculeData mol;
        FILE* f = fopen(path.c_str(), "r");
        if (!f) {
            printf("[PDBParser] Failed to open: %s\n", path.c_str());
            return mol;
        }

        char line[256];
        while (fgets(line, sizeof(line), f)) {
            if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0) {
                Atom atom;
                // PDB format: x at 30-38, y at 38-46, z at 46-54
                char xstr[9], ystr[9], zstr[9];
                strncpy(xstr, line+30, 8); xstr[8] = 0;
                strncpy(ystr, line+38, 8); ystr[8] = 0;
                strncpy(zstr, line+46, 8); zstr[8] = 0;
                atom.position = {(float)atof(xstr), (float)atof(ystr), (float)atof(zstr)};

                char elem[3];
                parseElementFromPDB(line, elem);
                getElementProperties(elem, atom.element, atom.radius, atom.color);

                mol.atoms.push_back(atom);
            }
            else if (strncmp(line, "CONECT", 6) == 0) {
                // Parse connectivity
                int serial = atoi(line + 6);
                for (int col = 11; col < 31 && col + 5 <= (int)strlen(line); col += 5) {
                    char buf[6];
                    strncpy(buf, line + col, 5); buf[5] = 0;
                    int bonded = atoi(buf);
                    if (bonded > 0 && serial > 0 && serial <= (int)mol.atoms.size()
                        && bonded <= (int)mol.atoms.size() && serial < bonded) {
                        Bond b;
                        b.atomA = serial - 1;
                        b.atomB = bonded - 1;
                        b.order = 1;
                        mol.bonds.push_back(b);
                    }
                }
            }
        }
        fclose(f);

        // If no CONECT records, infer bonds from distance
        if (mol.bonds.empty()) inferBonds(mol);

        mol.computeBounds();
        printf("[PDBParser] Loaded %s: %d atoms, %d bonds\n",
               path.c_str(), (int)mol.atoms.size(), (int)mol.bonds.size());
        return mol;
    }

    // Parse an SDF/MOL file
    static MoleculeData parseSDF(const std::string& path) {
        MoleculeData mol;
        FILE* f = fopen(path.c_str(), "r");
        if (!f) {
            printf("[PDBParser] Failed to open: %s\n", path.c_str());
            return mol;
        }

        char line[256];
        // Skip first 3 header lines
        for (int i = 0; i < 3; i++) {
            if (!fgets(line, sizeof(line), f)) { fclose(f); return mol; }
        }

        // Counts line: aaabbblll... (atoms, bonds, ...)
        if (!fgets(line, sizeof(line), f)) { fclose(f); return mol; }
        int nAtoms = 0, nBonds = 0;
        sscanf(line, "%d %d", &nAtoms, &nBonds);

        // Atom block
        for (int i = 0; i < nAtoms; i++) {
            if (!fgets(line, sizeof(line), f)) break;
            Atom atom;
            char elemStr[4] = {};
            float px, py, pz;
            sscanf(line, "%f %f %f %3s", &px, &py, &pz, elemStr);
            atom.position = {px, py, pz};
            // Trim element string
            for (int j = 0; j < 3; j++) if (elemStr[j] == ' ') elemStr[j] = 0;
            getElementProperties(elemStr, atom.element, atom.radius, atom.color);
            mol.atoms.push_back(atom);
        }

        // Bond block
        for (int i = 0; i < nBonds; i++) {
            if (!fgets(line, sizeof(line), f)) break;
            Bond bond;
            int a, b, order;
            sscanf(line, "%d %d %d", &a, &b, &order);
            bond.atomA = a - 1;  // 1-indexed → 0-indexed
            bond.atomB = b - 1;
            bond.order = order;
            if (bond.atomA >= 0 && bond.atomB >= 0
                && bond.atomA < (int)mol.atoms.size()
                && bond.atomB < (int)mol.atoms.size()) {
                mol.bonds.push_back(bond);
            }
        }

        fclose(f);
        mol.computeBounds();
        printf("[PDBParser] Loaded %s: %d atoms, %d bonds\n",
               path.c_str(), (int)mol.atoms.size(), (int)mol.bonds.size());
        return mol;
    }

    // Auto-detect format by extension
    static MoleculeData parse(const std::string& path) {
        if (path.size() > 4) {
            std::string ext = path.substr(path.size() - 4);
            if (ext == ".pdb" || ext == ".PDB") return parsePDB(path);
            if (ext == ".sdf" || ext == ".SDF" || ext == ".mol" || ext == ".MOL")
                return parseSDF(path);
        }
        printf("[PDBParser] Unknown format: %s\n", path.c_str());
        return MoleculeData{};
    }

private:
    // Infer bonds from interatomic distances (for PDB files without CONECT)
    static void inferBonds(MoleculeData& mol) {
        int n = (int)mol.atoms.size();
        // Limit bond inference to reasonable-sized molecules
        int limit = (n > 5000) ? 5000 : n;
        for (int i = 0; i < limit; i++) {
            for (int j = i + 1; j < limit; j++) {
                float dx = mol.atoms[i].position.x - mol.atoms[j].position.x;
                float dy = mol.atoms[i].position.y - mol.atoms[j].position.y;
                float dz = mol.atoms[i].position.z - mol.atoms[j].position.z;
                float dist = sqrtf(dx*dx + dy*dy + dz*dz);
                // Typical covalent bond: < sum of covalent radii + tolerance
                float maxBond = (mol.atoms[i].radius + mol.atoms[j].radius) * 0.55f + 0.3f;
                if (dist < maxBond && dist > 0.4f) {
                    mol.bonds.push_back({i, j, 1});
                }
            }
        }
    }
};
