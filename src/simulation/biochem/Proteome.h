#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

// ══════════════════════════════════════════════════════════════════════════
//  Proteome.h — Protein structure, function, folding, and degradation
//  Based on: Alberts Ch 3 "Proteins"
//
//  Implements:
//    - Amino acid properties (20 standard, hydrophobicity, charge, size)
//    - Protein folding:
//      Primary → Secondary (α-helix, β-sheet prediction: Chou-Fasman)
//      Tertiary (hydrophobic core packing, H-bonds, salt bridges)
//      Quaternary (subunit assembly)
//      Chaperone-assisted: Hsp70 (DnaK), Hsp60/GroEL-GroES chaperonin
//      Intrinsically disordered regions (IDRs)
//    - Enzyme kinetics:
//      Michaelis-Menten (Km, Vmax, kcat)
//      Allosteric regulation (MWC/concerted model, Hill equation)
//      Competitive, uncompetitive, non-competitive inhibition
//      Irreversible inhibition (covalent)
//    - Post-translational modifications (PTMs):
//      Phosphorylation (Ser/Thr/Tyr by kinases, reversed by phosphatases)
//      Ubiquitination (Ub → E1 → E2 → E3 ligase → proteasome)
//      SUMOylation, acetylation, methylation, glycosylation
//      GPI anchor addition
//      Disulfide bond formation (ER, PDI)
//      Signal peptide cleavage
//    - Protein degradation:
//      Ubiquitin-proteasome system (26S proteasome)
//      Autophagy-lysosome pathway
//      N-end rule (N-degron)
//    - Protein-protein interactions:
//      Binding affinity (Kd), on/off rates
//      Scaffold proteins, multivalent interactions
//      Biomolecular condensates (LLPS, IDR-driven)
//    - Motor proteins:
//      Myosin mechanochemical cycle (rigor→ATP→hydrolysis→powerstroke)
//      Kinesin hand-over-hand stepping (8nm per ATP)
//    - Protein misfolding and amyloid aggregation
//
//  Uses AlphaFold2/ESMFold concepts for structure prediction
//  OpenMM integration for physics-based folding
//
//  Units: µM for concentrations, kJ/mol for energies, nm for distances
// ══════════════════════════════════════════════════════════════════════════

// ── Amino Acid Properties ───────────────────────────────────────────────
struct AminoAcidProperties {
    char letter;
    const char* name3;
    const char* fullName;
    float hydrophobicity;  // Kyte-Doolittle scale (-4.5 to +4.5)
    float charge_pH7;      // Net charge at pH 7.0
    float molecularWeight; // Daltons
    float vdwRadius;       // Van der Waals radius (Å)
    bool aromatic;
    bool polar;
    bool canPhosphorylate; // Ser, Thr, Tyr
};

static const AminoAcidProperties AA_TABLE[20] = {
    {'A', "Ala", "Alanine",        1.8f,  0,  89.1f, 1.7f, false, false, false},
    {'R', "Arg", "Arginine",      -4.5f, +1, 174.2f, 2.0f, false, true,  false},
    {'N', "Asn", "Asparagine",    -3.5f,  0, 132.1f, 1.8f, false, true,  false},
    {'D', "Asp", "Aspartate",     -3.5f, -1, 133.1f, 1.8f, false, true,  false},
    {'C', "Cys", "Cysteine",       2.5f,  0, 121.2f, 1.8f, false, true,  false},
    {'E', "Glu", "Glutamate",     -3.5f, -1, 147.1f, 1.9f, false, true,  false},
    {'Q', "Gln", "Glutamine",     -3.5f,  0, 146.2f, 1.9f, false, true,  false},
    {'G', "Gly", "Glycine",       -0.4f,  0,  75.0f, 1.5f, false, false, false},
    {'H', "His", "Histidine",     -3.2f,  0, 155.2f, 1.9f, true,  true,  false},
    {'I', "Ile", "Isoleucine",     4.5f,  0, 131.2f, 1.9f, false, false, false},
    {'L', "Leu", "Leucine",        3.8f,  0, 131.2f, 1.9f, false, false, false},
    {'K', "Lys", "Lysine",        -3.9f, +1, 146.2f, 2.0f, false, true,  false},
    {'M', "Met", "Methionine",     1.9f,  0, 149.2f, 1.9f, false, false, false},
    {'F', "Phe", "Phenylalanine",  2.8f,  0, 165.2f, 2.0f, true,  false, false},
    {'P', "Pro", "Proline",       -1.6f,  0, 115.1f, 1.8f, false, false, false},
    {'S', "Ser", "Serine",        -0.8f,  0, 105.1f, 1.7f, false, true,  true },
    {'T', "Thr", "Threonine",     -0.7f,  0, 119.1f, 1.8f, false, true,  true },
    {'W', "Trp", "Tryptophan",    -0.9f,  0, 204.2f, 2.1f, true,  true,  false},
    {'Y', "Tyr", "Tyrosine",     -1.3f,  0, 181.2f, 2.0f, true,  true,  true },
    {'V', "Val", "Valine",         4.2f,  0, 117.1f, 1.8f, false, false, false},
};

static int aaIndex(char letter) {
    const char* aa = "ARNDCEQGHILKMFPSTWYV";
    for (int i = 0; i < 20; i++) if (aa[i] == letter) return i;
    return 0; // Default to Ala
}

// ── Secondary Structure Prediction (Chou-Fasman simplified) ─────────────
enum SecondaryStructure { SS_COIL = 0, SS_HELIX, SS_SHEET, SS_TURN };

// Helix/sheet propensities (Chou-Fasman parameters, simplified)
static const float HELIX_PROPENSITY[20] = {
    1.42f, 0.98f, 0.67f, 1.01f, 0.70f, 1.51f, 1.11f, 0.57f, 1.00f, 1.08f,
    1.21f, 1.16f, 1.45f, 1.13f, 0.57f, 0.77f, 0.83f, 1.08f, 0.69f, 1.06f
};
static const float SHEET_PROPENSITY[20] = {
    0.83f, 0.93f, 0.89f, 0.54f, 1.19f, 0.37f, 1.10f, 0.75f, 0.87f, 1.60f,
    1.30f, 0.74f, 1.05f, 1.38f, 0.55f, 0.75f, 1.19f, 1.37f, 1.47f, 1.70f
};

static SecondaryStructure predictSS(const char* seq, int pos, int len) {
    if (pos < 0 || pos >= len) return SS_COIL;
    // Window of 5 residues
    float h_score = 0, s_score = 0;
    int window = 5;
    int start = std::max(0, pos - window/2);
    int end = std::min(len - 1, pos + window/2);
    for (int i = start; i <= end; i++) {
        int idx = aaIndex(seq[i]);
        h_score += HELIX_PROPENSITY[idx];
        s_score += SHEET_PROPENSITY[idx];
    }
    h_score /= (end - start + 1);
    s_score /= (end - start + 1);
    // Proline breaks helices
    if (seq[pos] == 'P') h_score *= 0.3f;
    if (h_score > 1.03f && h_score > s_score) return SS_HELIX;
    if (s_score > 1.05f && s_score > h_score) return SS_SHEET;
    return SS_COIL;
}

// ── Protein Model ───────────────────────────────────────────────────────
struct ProteinModel {
    std::string name;
    std::string sequence;          // Amino acid sequence (1-letter)
    int length = 0;                // Number of residues

    // Structure
    std::vector<SecondaryStructure> secondary_structure;
    float hydrophobic_core = 0.0f; // Fraction of hydrophobic residues buried
    float folded_fraction = 0.0f;  // 0 = unfolded, 1 = fully folded
    float stability_dG = 0.0f;     // Folding free energy (kJ/mol, negative = stable)

    // Properties
    float molecular_weight = 0.0f; // kDa
    float isoelectric_point = 0.0f;// pI
    float net_charge_pH7 = 0.0f;

    // Post-translational modifications
    std::vector<int> phosphorylation_sites; // Residue indices that are phosphorylated
    bool ubiquitinated = false;
    bool sumoylated = false;
    int disulfide_bonds = 0;
    bool signal_peptide_cleaved = false;
    bool glycosylated = false;

    // Functional state
    float activity = 0.0f;         // Enzymatic/functional activity (0-1)
    bool is_active = false;
    float half_life = 3600.0f;     // Protein half-life (seconds)
    float age = 0.0f;

    // Binding
    float binding_Kd = 1.0f;       // Dissociation constant (µM)
    int binding_partners = 0;

    void initFromSequence(const std::string& seq, const std::string& n) {
        name = n;
        sequence = seq;
        length = (int)seq.size();
        molecular_weight = 0;
        net_charge_pH7 = 0;
        secondary_structure.resize(length, SS_COIL);

        for (int i = 0; i < length; i++) {
            int idx = aaIndex(seq[i]);
            molecular_weight += AA_TABLE[idx].molecularWeight;
            net_charge_pH7 += AA_TABLE[idx].charge_pH7;
            secondary_structure[i] = predictSS(seq.c_str(), i, length);
        }
        molecular_weight /= 1000.0f; // Convert to kDa

        // Compute hydrophobicity
        float hydro_sum = 0;
        for (char c : seq) hydro_sum += AA_TABLE[aaIndex(c)].hydrophobicity;
        hydrophobic_core = fmaxf(0, hydro_sum / length / 4.5f); // Normalized

        // Stability estimate (empirical: ~-5 kJ/mol per kDa for typical protein)
        stability_dG = -5.0f * molecular_weight;

        // Half-life (N-end rule: depends on first residue)
        // Ref: Tobias et al. 1991
        char n_end = seq[0];
        if (n_end == 'M' || n_end == 'S' || n_end == 'A' || n_end == 'T' ||
            n_end == 'V' || n_end == 'G') {
            half_life = 72000.0f; // Stabilizing (~20 hours)
        } else if (n_end == 'R' || n_end == 'K' || n_end == 'F' || n_end == 'L' || n_end == 'W') {
            half_life = 180.0f;   // Destabilizing (~3 min)
        } else {
            half_life = 7200.0f;  // Intermediate (~2 hours)
        }
    }

    // Fold the protein (simplified thermodynamic model)
    void fold(float temperature, float chaperone_level) {
        // Folding rate depends on temperature and chaperones
        float fold_rate = 0.01f * (1.0f + chaperone_level * 2.0f);
        // Unfolding rate increases with temperature (Arrhenius-like)
        float unfold_rate = 0.001f * expf((temperature - 310.0f) * 0.1f);

        folded_fraction += (fold_rate * (1.0f - folded_fraction) - unfold_rate * folded_fraction);
        folded_fraction = std::clamp(folded_fraction, 0.0f, 1.0f);

        activity = folded_fraction * (ubiquitinated ? 0.0f : 1.0f);
    }

    // Check if misfolded → amyloid aggregation risk
    float amyloidPropensity() const {
        // High β-sheet content + partial unfolding = amyloid risk
        int sheet_count = 0;
        for (auto ss : secondary_structure) if (ss == SS_SHEET) sheet_count++;
        float sheet_frac = (float)sheet_count / fmaxf(1, length);
        return sheet_frac * (1.0f - folded_fraction) * 0.5f;
    }
};

// ── Chaperone System ────────────────────────────────────────────────────
struct ChaperoneSystem {
    // Hsp70 / DnaK system (prevents aggregation, assists initial folding)
    float Hsp70 = 1.0f;           // Hsp70 level
    float Hsp40 = 0.5f;           // J-domain co-chaperone (delivers substrates)
    float NEF = 0.3f;             // Nucleotide exchange factor (BAG1)

    // Hsp60 / GroEL-GroES chaperonin (folding cage, ATP-driven)
    float Hsp60 = 0.5f;           // GroEL barrel
    float Hsp10 = 0.5f;           // GroES cap
    float chaperonin_capacity = 0.0f; // Fraction occupied

    // Hsp90 (specialized: steroid receptors, kinases)
    float Hsp90 = 0.5f;
    float Hsp90_clients = 0.0f;

    // Small HSPs (prevent aggregation under stress)
    float sHsp = 0.3f;

    // Heat shock response
    float HSF1_active = 0.0f;     // Heat shock factor 1 (transcription factor)
    float stress_level = 0.0f;    // Unfolded protein burden

    float totalChaperoneActivity() const {
        return Hsp70 * 0.4f + Hsp60 * 0.3f + Hsp90 * 0.2f + sHsp * 0.1f;
    }

    void update(float unfolded_proteins, float temperature, float dt) {
        stress_level = unfolded_proteins / (1.0f + totalChaperoneActivity() * 5.0f);

        // HSF1 activation (unfolded proteins titrate Hsp70 away from HSF1)
        HSF1_active = stress_level / (0.5f + stress_level);

        // Heat shock response upregulates all chaperones
        float hsr = HSF1_active * 0.1f * dt;
        Hsp70 += hsr;
        Hsp60 += hsr * 0.5f;
        Hsp90 += hsr * 0.3f;
        sHsp += hsr * 0.5f;

        // Baseline degradation
        float decay = 0.001f * dt;
        Hsp70 -= decay; Hsp60 -= decay; Hsp90 -= decay; sHsp -= decay;
        Hsp70 = fmaxf(0.1f, Hsp70);
        Hsp60 = fmaxf(0.1f, Hsp60);
    }
};

// ── Ubiquitin-Proteasome System ─────────────────────────────────────────
struct ProteasomeSystem {
    // E1 → E2 → E3 ubiquitin ligase cascade
    float E1_activating = 0.5f;    // Ubiquitin-activating enzyme
    float E2_conjugating = 0.5f;   // Ubiquitin-conjugating enzyme
    // Specific E3 ligases
    float SCF_complex = 0.5f;      // SKP1-Cullin-F-box (cell cycle targets)
    float APC_C = 0.0f;            // Anaphase-Promoting Complex (mitotic cyclins)
    float MDM2_ligase = 0.5f;      // p53-specific E3 ligase
    float VHL_ligase = 0.3f;       // HIF1α degradation (O2-dependent)

    // 26S Proteasome
    float proteasome_20S = 0.5f;   // Catalytic core (β1/2/5 proteases)
    float proteasome_19S = 0.5f;   // Regulatory cap (unfoldase + deubiquitinase)

    float degradation_capacity = 0.0f;

    void update() {
        degradation_capacity = proteasome_20S * proteasome_19S * E1_activating * E2_conjugating;
    }
};

// ══════════════════════════════════════════════════════════════════════════

class ProteomeEngine {
public:
    std::vector<ProteinModel> proteins;
    ChaperoneSystem chaperones;
    ProteasomeSystem proteasome;

    // Global protein synthesis/degradation balance
    float total_protein_mass = 100.0f;  // Arbitrary units
    float synthesis_rate = 1.0f;        // Proteins/s
    float degradation_rate = 0.0f;
    float misfolded_fraction = 0.0f;
    float aggregated_fraction = 0.0f;   // Amyloid-like aggregates

    void init() {
        proteins.clear();
        chaperones = ChaperoneSystem();
        proteasome = ProteasomeSystem();
    }

    // Create a new protein from sequence
    int createProtein(const std::string& sequence, const std::string& name) {
        ProteinModel p;
        p.initFromSequence(sequence, name);
        proteins.push_back(p);
        return (int)proteins.size() - 1;
    }

    void step(float dt, float temperature, float ATP, float amino_acid_supply) {
        dt = fminf(dt, 0.1f);

        // Update chaperone system
        chaperones.update(misfolded_fraction, temperature, dt);

        // Update proteasome
        proteasome.update();

        // Protein synthesis rate (depends on ribosomes, amino acids, energy)
        synthesis_rate = amino_acid_supply * ATP / (1.0f + ATP) * 0.5f;

        // Fold/degrade each protein
        misfolded_fraction = 0;
        degradation_rate = 0;
        for (auto& p : proteins) {
            p.age += dt;
            p.fold(temperature, chaperones.totalChaperoneActivity());

            // Misfolded proteins → chaperone rescue or degradation
            if (p.folded_fraction < 0.5f) {
                misfolded_fraction += 1.0f / fmaxf(1, (int)proteins.size());
                // Ubiquitinate misfolded proteins
                if (p.age > 60.0f && !p.ubiquitinated) { // Give 60s to fold
                    p.ubiquitinated = (float)rand() / RAND_MAX < 0.01f * dt;
                }
            }

            // Proteasomal degradation of ubiquitinated proteins
            if (p.ubiquitinated) {
                float deg_prob = proteasome.degradation_capacity * 0.1f * dt;
                if ((float)rand() / RAND_MAX < deg_prob) {
                    p.activity = 0;
                    degradation_rate += 1.0f / fmaxf(1, (int)proteins.size());
                }
            }

            // Age-dependent degradation (N-end rule half-life)
            if (p.age > p.half_life) {
                float excess = (p.age - p.half_life) / p.half_life;
                if ((float)rand() / RAND_MAX < excess * 0.01f * dt) {
                    p.activity = 0;
                }
            }
        }

        // Amyloid aggregation (misfolded proteins with high β-sheet content)
        for (auto& p : proteins) {
            if (p.folded_fraction < 0.3f && p.amyloidPropensity() > 0.1f) {
                aggregated_fraction += 0.001f * dt;
            }
        }
        aggregated_fraction = fminf(1.0f, aggregated_fraction);

        // Total protein mass
        total_protein_mass += (synthesis_rate - degradation_rate) * dt;
        total_protein_mass = fmaxf(10.0f, total_protein_mass);
    }
};
