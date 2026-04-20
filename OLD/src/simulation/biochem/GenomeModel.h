#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  GenomeModel.h — Full genome representation
//  Based on: Alberts Ch 4 "DNA, Chromosomes, and Genomes"
//
//  Implements:
//    - 23 chromosome pairs (diploid human genome)
//    - Chromatin architecture:
//      Nucleosome (147bp + histone octamer H2A/H2B/H3/H4)
//      30nm fiber (H1 linker histone)
//      Chromosome loops (50-200kb, cohesin/condensin)
//      Chromosome territories
//      TADs (Topologically Associating Domains, CTCF/cohesin boundaries)
//    - Euchromatin vs heterochromatin (constitutive vs facultative)
//    - Histone variants (H3.3, CENP-A, H2A.X, macroH2A)
//    - Telomere structure (TTAGGG repeats + shelterin complex)
//    - Centromere (α-satellite DNA + CENP-A)
//    - Gene density per chromosome
//    - Ploidy tracking (diploid → tetraploid during S-phase)
//    - Transposable elements (LINE, SINE/Alu, ~45% of genome)
//    - SNPs and CNVs (genetic variation)
//
//  Simplified representation: each chromosome stored as metadata
//  (not full 3.2 billion bp sequence — that would be ~800MB RAM per cell)
//
//  Units: bp for lengths, Mb for large regions
// ══════════════════════════════════════════════════════════════════════════

// ── Chromosome data (human genome, hg38) ────────────────────────────────
struct ChromosomeInfo {
    int number;              // 1-22 + 23(X) or 24(Y)
    const char* name;
    long long size_bp;       // Size in base pairs
    int gene_count;          // Estimated protein-coding genes
    float gc_content;        // GC content fraction
    float gene_density;      // Genes per Mb
    float repeat_fraction;   // Fraction that is repetitive DNA
};

static const ChromosomeInfo HUMAN_CHROMOSOMES[24] = {
    { 1, "chr1",  248956422, 2058, 0.417f, 8.27f, 0.447f},
    { 2, "chr2",  242193529, 1309, 0.402f, 5.40f, 0.468f},
    { 3, "chr3",  198295559, 1078, 0.397f, 5.44f, 0.459f},
    { 4, "chr4",  190214555,  752, 0.384f, 3.95f, 0.480f},
    { 5, "chr5",  181538259,  876, 0.395f, 4.83f, 0.466f},
    { 6, "chr6",  170805979, 1048, 0.397f, 6.14f, 0.453f},
    { 7, "chr7",  159345973,  989, 0.407f, 6.21f, 0.472f},
    { 8, "chr8",  145138636,  677, 0.401f, 4.66f, 0.468f},
    { 9, "chr9",  138394717,  786, 0.414f, 5.68f, 0.436f},
    {10, "chr10", 133797422,  733, 0.415f, 5.48f, 0.427f},
    {11, "chr11", 135086622, 1298, 0.415f, 9.61f, 0.441f},
    {12, "chr12", 133275309, 1034, 0.408f, 7.76f, 0.441f},
    {13, "chr13", 114364328,  327, 0.385f, 2.86f, 0.505f},
    {14, "chr14", 107043718,  830, 0.409f, 7.76f, 0.453f},
    {15, "chr15", 101991189,  613, 0.420f, 6.01f, 0.437f},
    {16, "chr16",  90338345,  873, 0.448f, 9.66f, 0.410f},
    {17, "chr17",  83257441, 1197, 0.455f, 14.37f, 0.375f},
    {18, "chr18",  80373285,  270, 0.397f, 3.36f, 0.473f},
    {19, "chr19",  58617616, 1472, 0.481f, 25.11f, 0.377f},
    {20, "chr20",  64444167,  544, 0.441f, 8.44f, 0.430f},
    {21, "chr21",  46709983,  234, 0.408f, 5.01f, 0.452f},
    {22, "chr22",  50818468,  488, 0.480f, 9.60f, 0.421f},
    {23, "chrX",  156040895,  842, 0.394f, 5.40f, 0.563f},
    {24, "chrY",   57227415,   63, 0.395f, 1.10f, 0.638f},
};

// ── Chromosome State ────────────────────────────────────────────────────
struct GenomeChromosomeState {
    int id;                          // Chromosome number
    bool maternal;                   // Maternal or paternal copy

    // Structural state
    float condensation = 0.0f;       // 0 = interphase, 1 = metaphase
    float euchromatin_fraction = 0.7f;// Fraction in open chromatin
    float heterochromatin_fraction = 0.3f;

    // Telomere
    float telomere_length = 10000.0f;// bp (5000-15000 normal range)
    float shelterin_bound = 1.0f;    // Shelterin complex occupancy

    // Centromere
    float CENP_A_loading = 1.0f;     // CENP-A (centromeric H3 variant)
    bool kinetochore_assembled = false;
    bool spindle_attached = false;

    // Replication
    float replication_progress = 0.0f;// 0 = unreplicated, 1 = fully replicated
    int active_origins = 0;
    bool fully_replicated = false;

    // Cohesin (holds sister chromatids)
    float cohesin_arm = 1.0f;        // Arm cohesin (removed in anaphase I of meiosis)
    float cohesin_centromere = 1.0f;  // Centromeric cohesin (protected by Shugoshin)

    // Damage
    int unrepaired_DSBs = 0;
    float gammaH2AX_foci = 0.0f;

    // TAD structure (simplified — number of TADs)
    int num_TADs = 0;
    float CTCF_binding = 1.0f;       // CTCF boundary protein
    float cohesin_loading = 1.0f;    // Loop extrusion activity

    // Histone modifications (chromosome-wide average)
    float H3K4me3_avg = 0.3f;        // Active promoters
    float H3K27me3_avg = 0.2f;       // Polycomb repression
    float H3K9me3_avg = 0.1f;        // Constitutive heterochromatin
    float H3K36me3_avg = 0.2f;       // Actively transcribed genes
    float acetylation_avg = 0.4f;    // General acetylation level
};

// ── Nuclear Organization ────────────────────────────────────────────────
struct NuclearOrganization {
    // Chromosome territories (each chr occupies a distinct nuclear volume)
    // Gene-rich chromosomes tend to be interior, gene-poor at periphery
    float territory_positions[24] = {}; // 0 = center, 1 = periphery

    // Nuclear lamina association
    float LAD_fraction = 0.3f;          // Lamina-Associated Domains (repressed)

    // Nucleolus (rRNA gene transcription — chr 13, 14, 15, 21, 22)
    float nucleolus_activity = 1.0f;

    // Nuclear pore complexes
    int NPC_count = 3000;               // ~3000-5000 per nucleus

    // Nuclear envelope
    float lamin_A = 1.0f;
    float lamin_B = 1.0f;
    float envelope_integrity = 1.0f;

    void initTerritories() {
        for (int i = 0; i < 24; i++) {
            // Gene-rich chromosomes (high gene density) tend toward interior
            float density = HUMAN_CHROMOSOMES[i].gene_density;
            territory_positions[i] = 1.0f - fminf(1.0f, density / 25.0f);
        }
    }
};

// ══════════════════════════════════════════════════════════════════════════

class GenomeEngine {
public:
    // 46 chromosomes (23 pairs: maternal + paternal)
    std::vector<GenomeChromosomeState> chromosomes;
    NuclearOrganization nucleus;

    // Genome-wide stats
    long long total_bp = 0;
    int total_genes = 0;
    float ploidy = 2.0f;                // 2n = diploid, 4n = after S-phase
    bool s_phase_complete = false;

    // Transposable elements (genome-wide)
    float LINE_fraction = 0.20f;        // Long Interspersed Nuclear Elements
    float SINE_fraction = 0.13f;        // Short (Alu elements)
    float LTR_fraction = 0.08f;         // LTR retrotransposons
    float DNA_transposon_fraction = 0.03f;

    // Epigenome summary
    float global_methylation = 0.7f;    // CpG methylation (genome-wide)
    float X_inactivation = 0.0f;        // 0 = both active, 1 = one inactivated (females)
    float imprinted_genes_active = 0.5f;// Imprinted loci properly regulated

    void init(bool is_female = true) {
        chromosomes.clear();
        total_bp = 0;
        total_genes = 0;

        // Create 46 chromosomes (23 maternal + 23 paternal)
        for (int copy = 0; copy < 2; copy++) {
            for (int i = 0; i < 23; i++) {
                GenomeChromosomeState c;
                c.id = i + 1;
                c.maternal = (copy == 0);

                // Use actual chromosome data
                auto& info = HUMAN_CHROMOSOMES[i];
                c.telomere_length = 8000.0f + (float)(rand() % 6000); // 8000-14000 bp
                c.euchromatin_fraction = 1.0f - info.repeat_fraction;
                c.num_TADs = (int)(info.size_bp / 500000); // ~1 TAD per 500kb

                total_bp += info.size_bp;
                total_genes += info.gene_count;

                // Second X in females
                if (i == 22 && copy == 1 && is_female) {
                    // X inactivation (Xist-mediated)
                    X_inactivation = 1.0f;
                    c.H3K27me3_avg = 0.8f; // Heavily repressed
                    c.euchromatin_fraction = 0.1f;
                }

                // Y chromosome (only in males)
                if (i == 23 && is_female) continue; // Skip Y for females

                chromosomes.push_back(c);
            }
        }

        nucleus.initTerritories();
        s_phase_complete = false;
        ploidy = 2.0f;
    }

    // S-phase replication
    void replicateGenome(float progress) {
        for (auto& c : chromosomes) {
            // Euchromatin replicates early, heterochromatin late
            float timing = c.euchromatin_fraction; // Early-replicating fraction
            if (progress > (1.0f - timing)) {
                c.replication_progress = fminf(1.0f,
                    (progress - (1.0f - timing)) / timing);
            }
            c.fully_replicated = c.replication_progress >= 0.99f;
        }
        s_phase_complete = progress >= 0.99f;
        ploidy = 2.0f + 2.0f * progress; // Increases from 2n to 4n
    }

    // Telomere shortening (called after each complete cell division)
    void shortenTelomeres() {
        for (auto& c : chromosomes) {
            float loss = 50.0f + (float)(rand() % 50); // 50-100 bp per division
            c.telomere_length -= loss;
        }
    }

    // Check if any telomere is critically short
    bool hasCriticallyShortTelomere() const {
        for (auto& c : chromosomes) {
            if (c.telomere_length < 4000.0f) return true;
        }
        return false;
    }

    float averageTelomereLength() const {
        float sum = 0;
        for (auto& c : chromosomes) sum += c.telomere_length;
        return sum / fmaxf(1, (int)chromosomes.size());
    }

    // Mitotic condensation
    void condenseChromosomes(float level) {
        for (auto& c : chromosomes) {
            c.condensation = level;
            if (level > 0.5f) {
                c.kinetochore_assembled = true;
                c.CENP_A_loading = 1.0f;
            }
        }
    }
};
