#pragma once

#include <string>
#include <vector>

// ══════════════════════════════════════════════════════════════════════════
//  MoleculeLibrary — Predefined cellular molecules with SMILES/sequences
//  Used for on-demand generation via RDKit/ESMFold
// ══════════════════════════════════════════════════════════════════════════

enum MoleculeType {
    MOL_SMALL_MOLECULE = 0,  // Generate via RDKit
    MOL_PROTEIN        = 1,  // Generate via ESMFold
};

struct MoleculeEntry {
    const char* id;           // Filename stem (e.g., "atp")
    const char* name;         // Display name
    MoleculeType type;
    const char* data;         // SMILES for small molecules, AA sequence for proteins
    float baseColor[3];       // Default rendering color (RGB 0-1)
    float scale;              // Rendering scale multiplier
};

// Common cellular molecules — matches scripts/generate_library.py
static const MoleculeEntry MOLECULE_ENTRIES[] = {
    // Energy molecules
    {"atp",        "ATP",               MOL_SMALL_MOLECULE,
     "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
     {0.2f, 0.8f, 0.3f}, 0.3f},
    {"adp",        "ADP",               MOL_SMALL_MOLECULE,
     "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
     {0.3f, 0.7f, 0.2f}, 0.3f},
    {"glucose",    "Glucose",           MOL_SMALL_MOLECULE,
     "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
     {0.9f, 0.7f, 0.1f}, 0.4f},
    {"nadh",       "NADH",              MOL_SMALL_MOLECULE,
     "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]2OC(n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",
     {0.1f, 0.3f, 0.9f}, 0.3f},
    {"water",      "Water",             MOL_SMALL_MOLECULE,
     "O",
     {0.8f, 0.85f, 0.9f}, 0.25f},
    {"calcium_ion","Calcium ion",       MOL_SMALL_MOLECULE,
     "[Ca+2]",
     {0.3f, 0.7f, 0.7f}, 0.35f},

    // Membrane lipids
    {"dppc",       "Phospholipid",      MOL_SMALL_MOLECULE,
     "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC",
     {0.8f, 0.6f, 0.2f}, 0.2f},
    {"cholesterol","Cholesterol",       MOL_SMALL_MOLECULE,
     "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
     {0.9f, 0.85f, 0.3f}, 0.25f},

    // Nucleotides
    {"deoxyadenosine", "Deoxyadenosine", MOL_SMALL_MOLECULE,
     "c1nc(c2c(n1)n(cn2)[C@@H]3C[C@@H]([C@H](O3)CO)O)N",
     {0.3f, 0.5f, 0.9f}, 0.35f},

    // Amino acids
    {"glycine",    "Glycine",           MOL_SMALL_MOLECULE,
     "NCC(=O)O",
     {0.85f, 0.35f, 0.35f}, 0.45f},
    {"tryptophan", "Tryptophan",        MOL_SMALL_MOLECULE,
     "c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N",
     {0.7f, 0.3f, 0.7f}, 0.5f},

    // Signaling
    {"camp",       "Cyclic AMP",        MOL_SMALL_MOLECULE,
     "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N",
     {0.5f, 0.9f, 0.5f}, 0.35f},

    // Key cellular proteins (short representative sequences for demo)
    {"actin",      "Actin",             MOL_PROTEIN,
     "MDSEVAALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEH",
     {0.9f, 0.3f, 0.3f}, 0.05f},
    {"tubulin",    "Tubulin",           MOL_PROTEIN,
     "MRECISIHVGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGKHVPRAVFVDLEPTV",
     {0.3f, 0.9f, 0.5f}, 0.05f},
};

static const int MOLECULE_ENTRY_COUNT = sizeof(MOLECULE_ENTRIES) / sizeof(MoleculeEntry);
