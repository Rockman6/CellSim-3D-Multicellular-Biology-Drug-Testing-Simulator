#pragma once

#include <string>
#include <cstdio>
#include <cstdlib>

// ══════════════════════════════════════════════════════════════════════════
//  MolGenRunner — Run Python scripts for molecular structure generation
//  Wraps RDKit (small molecules) and ESMFold (proteins)
// ══════════════════════════════════════════════════════════════════════════

class MolGenRunner {
public:
    void init(const std::string& projectRoot) {
        projectRoot_ = projectRoot;
        venvPython_ = projectRoot + "/.molgen_venv/bin/python";
        scriptsDir_ = projectRoot + "/scripts";

        // Check if venv exists
        FILE* f = fopen(venvPython_.c_str(), "r");
        if (f) {
            fclose(f);
            venvReady_ = true;
            printf("[MolGenRunner] Venv ready at %s\n", venvPython_.c_str());
        } else {
            printf("[MolGenRunner] Venv not found. Run: python3 scripts/setup_molgen.py\n");
            venvReady_ = false;
        }
    }

    bool isReady() const { return venvReady_; }

    // Generate a small molecule from SMILES → .sdf file
    // Returns true on success
    bool generateMolecule(const std::string& smiles, const std::string& outputPath) {
        if (!venvReady_) return false;
        std::string cmd = venvPython_ + " " + scriptsDir_ + "/generate_molecule.py"
                          + " \"" + smiles + "\" \"" + outputPath + "\"";
        printf("[MolGenRunner] Running: %s\n", cmd.c_str());
        int ret = system(cmd.c_str());
        return ret == 0;
    }

    // Generate a protein from amino acid sequence → .pdb file
    // Note: ESMFold can be slow on CPU (~minutes for long sequences)
    bool generateProtein(const std::string& sequence, const std::string& outputPath) {
        if (!venvReady_) return false;
        std::string cmd = venvPython_ + " " + scriptsDir_ + "/generate_protein.py"
                          + " \"" + sequence + "\" \"" + outputPath + "\"";
        printf("[MolGenRunner] Running: %s\n", cmd.c_str());
        int ret = system(cmd.c_str());
        return ret == 0;
    }

    // Generate the entire predefined library
    bool generateLibrary(const std::string& outputDir) {
        if (!venvReady_) return false;
        std::string cmd = venvPython_ + " " + scriptsDir_ + "/generate_library.py"
                          + " \"" + outputDir + "\"";
        printf("[MolGenRunner] Running: %s\n", cmd.c_str());
        int ret = system(cmd.c_str());
        return ret == 0;
    }

private:
    std::string projectRoot_;
    std::string venvPython_;
    std::string scriptsDir_;
    bool venvReady_ = false;
};
