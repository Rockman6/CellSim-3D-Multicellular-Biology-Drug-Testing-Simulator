#pragma once
#include <cstring>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <cstdint>
#include <simd/simd.h>

// Local clampf so this header is independent of include order. Simulation.h
// also defines one (same body); both are `static`, so they don't collide.
#ifndef CENTRAL_DOGMA_CLAMPF_DEFINED
#define CENTRAL_DOGMA_CLAMPF_DEFINED
static float clampf(float v, float lo, float hi) { return fmaxf(lo, fminf(hi, v)); }
#endif

// ══════════════════════════════════════════════════════════════════════════
//  CentralDogma.h — Real gene sequence, codon table, intron/exon map
//  Gene: Human beta-globin (HBB), GenBank NM_000518 / genomic with introns
//  Implements: DNA replication, transcription, RNA splicing, translation
// ══════════════════════════════════════════════════════════════════════════

// ── HBB genomic sequence (1558 bp full length) ─────────────────────
// Source: Ensembl REST API, GRCh38 chr11:5225464-5227021 (-strand).
// Gene: HBB (beta-globin), NCBI Gene ID 3043. Includes 3 exons + 2 introns.
// Approximate structure (0-indexed inclusive):
//   Exon 1:   0-89    (90 bp)   — 5'UTR + first half of coding sequence
//   Intron 1: 90-218  (129 bp)  — retained for splicing visualization
//   Exon 2:  219-441  (223 bp)
//   Intron 2: 442-1293 (852 bp) — full length, was truncated to 120 bp before
//   Exon 3: 1294-1557 (264 bp)  — includes 3'UTR
// NOTE: The existing transcription/replication machinery operates on this
// entire sequence. Error-rate hotspots (CpG, trinucleotide repeats) in
// hotspotMultiplier() work at base-level so no changes needed there.
static const char HBB_SEQUENCE[] =
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG"
    "AACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG"
    "TTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAG"
    "ACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATT"
    "TTCCCACCCTTAGG"
    "CTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCA"
    "CTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGG"
    "TGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTG"
    "AGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGG"
    "GTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCAT"
    "GTCATAGGAAGGGGATAAGTAACAGGGTACAGTTTAGAATGGGAAACAGACGAATGAT"
    "TGCATCAGTGTGGAAGTCTCAGGATCGTTTTAGTTTCTTTTATTTGCTGTTCATAACA"
    "ATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTCTTCTCCGCAATTTTTACT"
    "ATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATT"
    "AAGTAACTTAAAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATA"
    "TGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTG"
    "ATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATGTGTACAC"
    "ATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTT"
    "TAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGG"
    "GCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATA"
    "ATTTCTGGGTTAAGGCAATAGCAATATCTCTGCATATAAATATTTCTGCATATAAATT"
    "GTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCT"
    "GCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTT"
    "TGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAG"
    "CTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCC"
    "CACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCA"
    "CAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTC"
    "CCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCT"
    "GCCTAATAAAAAACATTTATTTTCATTGCAA";

static const int HBB_LENGTH = 1558;
static float hotspotMultiplier(int bpIdx);

static int homopolymerRunLengthAt(const char* seq, int len, int idx) {
    if (idx < 0 || idx >= len) return 1;
    char base = seq[idx];
    int run = 1;
    for (int i = idx - 1; i >= 0 && seq[i] == base; --i) run++;
    for (int i = idx + 1; i < len && seq[i] == base; ++i) run++;
    return run;
}

static float localGCFraction(const char* seq, int len, int center, int radius) {
    int lo = std::max(0, center - radius);
    int hi = std::min(len - 1, center + radius);
    int total = 0;
    int gc = 0;
    for (int i = lo; i <= hi; i++) {
        total++;
        char b = seq[i];
        if (b == 'G' || b == 'g' || b == 'C' || b == 'c') gc++;
    }
    return total > 0 ? (float)gc / (float)total : 0.5f;
}

// ── Intron/Exon boundaries ──────────────────────────────────────────────
struct GeneRegion {
    int start, end;  // 0-indexed, inclusive
    bool isExon;
};

// Approximate exon/intron boundaries for the full 1558 bp HBB sequence.
// Used by isExon() for coloring (introns render dim) and by splicing state.
static const GeneRegion HBB_REGIONS[] = {
    {0,     89,  true},   // Exon 1 (90 bp, includes ATG start)
    {90,    218, false},  // Intron 1 (129 bp)
    {219,   441, true},   // Exon 2 (223 bp)
    {442,   1293, false}, // Intron 2 (852 bp, full length)
    {1294,  1557, true},  // Exon 3 (264 bp, includes 3'UTR + polyA)
};
static const int HBB_REGION_COUNT = 5;

static bool isExon(int bpIndex) {
    for (int i = 0; i < HBB_REGION_COUNT; i++) {
        if (bpIndex >= HBB_REGIONS[i].start && bpIndex <= HBB_REGIONS[i].end)
            return HBB_REGIONS[i].isExon;
    }
    return false;
}

// ── Base colors (educational standard) ──────────────────────────────────
struct BaseColor { float r, g, b; };

static BaseColor getBaseColor(char base, bool exon) {
    float dim = exon ? 1.0f : 0.3f;  // Introns rendered dim
    switch (base) {
        case 'A': case 'a': return {0.2f*dim, 0.4f*dim, 0.95f*dim}; // Blue — Adenine
        case 'T': case 't': return {0.2f*dim, 0.85f*dim, 0.3f*dim}; // Green — Thymine
        case 'G': case 'g': return {0.95f*dim, 0.2f*dim, 0.2f*dim}; // Red — Guanine
        case 'C': case 'c': return {0.95f*dim, 0.85f*dim, 0.15f*dim}; // Yellow — Cytosine
        case 'U': case 'u': return {0.9f*dim, 0.5f*dim, 0.1f*dim};  // Orange — Uracil (RNA)
        default: return {0.5f, 0.5f, 0.5f};
    }
}

static char complement(char base) {
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'G': case 'g': return 'C';
        case 'C': case 'c': return 'G';
        default: return 'N';
    }
}

static char transcribe(char dnaBase) {
    // Stored HBB sequence is the coding strand (starts with ATG).
    // Transcription therefore mirrors the coding strand with T -> U.
    switch (dnaBase) {
        case 'A': case 'a': return 'A';
        case 'T': case 't': return 'U';
        case 'G': case 'g': return 'G';
        case 'C': case 'c': return 'C';
        default: return 'N';
    }
}

// ── Standard Genetic Code ───────────────────────────────────────────────
// 64 codons → amino acid (single letter) + 3-letter name + color
struct AminoAcid {
    char letter;
    const char* name3;
    float r, g, b;  // Rendering color
};

static int codonIndex(char b1, char b2, char b3) {
    auto idx = [](char c) -> int {
        switch (c) {
            case 'U': case 'u': return 0;
            case 'C': case 'c': return 1;
            case 'A': case 'a': return 2;
            case 'G': case 'g': return 3;
            default: return 0;
        }
    };
    return idx(b1) * 16 + idx(b2) * 4 + idx(b3);
}

// Standard genetic code table (UUU=0, UUC=1, ..., GGG=63)
static const AminoAcid CODON_TABLE[64] = {
    // UUx
    {'F', "Phe", 0.8f, 0.6f, 0.2f}, {'F', "Phe", 0.8f, 0.6f, 0.2f},
    {'L', "Leu", 0.2f, 0.7f, 0.4f}, {'L', "Leu", 0.2f, 0.7f, 0.4f},
    // UCx
    {'S', "Ser", 0.9f, 0.8f, 0.3f}, {'S', "Ser", 0.9f, 0.8f, 0.3f},
    {'S', "Ser", 0.9f, 0.8f, 0.3f}, {'S', "Ser", 0.9f, 0.8f, 0.3f},
    // UAx
    {'Y', "Tyr", 0.6f, 0.3f, 0.8f}, {'Y', "Tyr", 0.6f, 0.3f, 0.8f},
    {'*', "Stop",0.9f, 0.1f, 0.1f}, {'*', "Stop",0.9f, 0.1f, 0.1f},
    // UGx
    {'C', "Cys", 0.9f, 0.9f, 0.2f}, {'C', "Cys", 0.9f, 0.9f, 0.2f},
    {'*', "Stop",0.9f, 0.1f, 0.1f}, {'W', "Trp", 0.5f, 0.2f, 0.7f},
    // CUx
    {'L', "Leu", 0.2f, 0.7f, 0.4f}, {'L', "Leu", 0.2f, 0.7f, 0.4f},
    {'L', "Leu", 0.2f, 0.7f, 0.4f}, {'L', "Leu", 0.2f, 0.7f, 0.4f},
    // CCx
    {'P', "Pro", 0.3f, 0.8f, 0.8f}, {'P', "Pro", 0.3f, 0.8f, 0.8f},
    {'P', "Pro", 0.3f, 0.8f, 0.8f}, {'P', "Pro", 0.3f, 0.8f, 0.8f},
    // CAx
    {'H', "His", 0.4f, 0.5f, 0.9f}, {'H', "His", 0.4f, 0.5f, 0.9f},
    {'Q', "Gln", 0.7f, 0.4f, 0.6f}, {'Q', "Gln", 0.7f, 0.4f, 0.6f},
    // CGx
    {'R', "Arg", 0.2f, 0.4f, 0.9f}, {'R', "Arg", 0.2f, 0.4f, 0.9f},
    {'R', "Arg", 0.2f, 0.4f, 0.9f}, {'R', "Arg", 0.2f, 0.4f, 0.9f},
    // AUx
    {'I', "Ile", 0.3f, 0.6f, 0.3f}, {'I', "Ile", 0.3f, 0.6f, 0.3f},
    {'I', "Ile", 0.3f, 0.6f, 0.3f}, {'M', "Met", 0.2f, 0.9f, 0.5f}, // AUG = Start
    // ACx
    {'T', "Thr", 0.7f, 0.5f, 0.3f}, {'T', "Thr", 0.7f, 0.5f, 0.3f},
    {'T', "Thr", 0.7f, 0.5f, 0.3f}, {'T', "Thr", 0.7f, 0.5f, 0.3f},
    // AAx
    {'N', "Asn", 0.8f, 0.4f, 0.5f}, {'N', "Asn", 0.8f, 0.4f, 0.5f},
    {'K', "Lys", 0.3f, 0.3f, 0.9f}, {'K', "Lys", 0.3f, 0.3f, 0.9f},
    // AGx
    {'S', "Ser", 0.9f, 0.8f, 0.3f}, {'S', "Ser", 0.9f, 0.8f, 0.3f},
    {'R', "Arg", 0.2f, 0.4f, 0.9f}, {'R', "Arg", 0.2f, 0.4f, 0.9f},
    // GUx
    {'V', "Val", 0.5f, 0.8f, 0.2f}, {'V', "Val", 0.5f, 0.8f, 0.2f},
    {'V', "Val", 0.5f, 0.8f, 0.2f}, {'V', "Val", 0.5f, 0.8f, 0.2f},
    // GCx
    {'A', "Ala", 0.6f, 0.6f, 0.6f}, {'A', "Ala", 0.6f, 0.6f, 0.6f},
    {'A', "Ala", 0.6f, 0.6f, 0.6f}, {'A', "Ala", 0.6f, 0.6f, 0.6f},
    // GAx
    {'D', "Asp", 0.9f, 0.3f, 0.4f}, {'D', "Asp", 0.9f, 0.3f, 0.4f},
    {'E', "Glu", 0.8f, 0.2f, 0.5f}, {'E', "Glu", 0.8f, 0.2f, 0.5f},
    // GGx
    {'G', "Gly", 0.7f, 0.7f, 0.7f}, {'G', "Gly", 0.7f, 0.7f, 0.7f},
    {'G', "Gly", 0.7f, 0.7f, 0.7f}, {'G', "Gly", 0.7f, 0.7f, 0.7f},
};

static AminoAcid translateCodon(char b1, char b2, char b3) {
    return CODON_TABLE[codonIndex(b1, b2, b3)];
}

static char rnaComplement(char base) {
    switch (base) {
        case 'A': case 'a': return 'U';
        case 'U': case 'u': return 'A';
        case 'G': case 'g': return 'C';
        case 'C': case 'c': return 'G';
        default: return 'N';
    }
}

static AminoAcid aminoAcidByLetter(char letter) {
    for (const auto& aa : CODON_TABLE) {
        if (aa.letter == letter && aa.letter != '*') return aa;
    }
    return {'?', "Unk", 0.7f, 0.7f, 0.7f};
}

static std::string anticodonForCodon(char b1, char b2, char b3) {
    std::string anti = "NNN";
    anti[0] = rnaComplement(b3);
    anti[1] = rnaComplement(b2);
    anti[2] = rnaComplement(b1);
    return anti;
}

static std::string buildHBBCodingMRNA() {
    std::string mrna;
    mrna.reserve(HBB_LENGTH);
    for (int bp = 0; bp < HBB_LENGTH; bp++) {
        if (!isExon(bp)) continue;
        mrna.push_back(transcribe(HBB_SEQUENCE[bp]));
    }
    return mrna;
}

static int findStartCodon(const std::string& mrna) {
    for (int i = 0; i + 2 < (int)mrna.size(); i++) {
        if (mrna[i] == 'A' && mrna[i + 1] == 'U' && mrna[i + 2] == 'G') return i;
    }
    return -1;
}

static int findStopCodon(const std::string& mrna, int startBase) {
    if (startBase < 0) return -1;
    for (int i = startBase; i + 2 < (int)mrna.size(); i += 3) {
        if (translateCodon(mrna[i], mrna[i + 1], mrna[i + 2]).letter == '*') return i;
    }
    return -1;
}

// ── Viral / extra gene template (Phase 8A) ──────────────────────────
//
// Part 8 §55 Stage 7: when a virion uncoats inside a cell, its genome
// becomes an additional template the host's RNA Pol II (or viral
// RdRp) transcribes alongside HBB. Each GeneTemplate is one gene
// added to the cell's gene registry. The existing HBB pipeline is
// unaffected — HBB is still referenced directly via HBB_SEQUENCE
// when transcript.geneIndex == -1 (the default).
//
// This is the load-bearing change that lets viral replication go
// through real transcription/translation rather than the toy
// `dG = rate × intra × dt` autocatalytic ODE.
enum class GenePolymerase : uint8_t {
    HOST_POLII   = 0,   // default host nuclear Pol II
    VIRAL_RdRp   = 1,   // virion-carried RNA-dependent RNA polymerase
    VIRAL_RT     = 2,   // retrovirus reverse transcriptase (cDNA path)
    VIRAL_POL    = 3,   // DNA-virus own polymerase (HBV, adeno)
};

struct GeneTemplate {
    std::string     id;                // e.g. "flu_HA_seg4"
    std::string     dnaSeq;            // DNA as sense strand
    int             startCodon    = -1; // position within mRNA (set after build)
    int             stopCodon     = -1;
    std::string     codingMRNA;        // cached mRNA built once
    float           promoterStrength = 1.0f; // 0..N; weighted random start
    GenePolymerase  polymerase    = GenePolymerase::HOST_POLII;
    std::string     proteinName;       // tag attached to released peptides
    int             copiesMade    = 0; // monitoring / assembly counter
};

// Forward declarations for mitosis (full definition below)
enum MitosisPhase {
    MITO_NONE = 0, MITO_PROPHASE, MITO_PROMETAPHASE, MITO_METAPHASE,
    MITO_ANAPHASE, MITO_TELOPHASE, MITO_CYTOKINESIS, MITO_COMPLETE,
};

// ── Central Dogma State Machine ─────────────────────────────────────────
// Tracks the progression of molecular biology processes in real-time

enum CDogmaPhase {
    CD_IDLE = 0,
    CD_REPLICATION,     // S-phase: DNA being copied
    CD_TRANSCRIPTION,   // RNA Pol II reading DNA → pre-mRNA
    CD_SPLICING,        // Spliceosome removing introns
    CD_EXPORT,          // mRNA exiting nucleus through pore
    CD_TRANSLATION,     // Ribosome reading mRNA → polypeptide
    CD_FOLDING,         // ER chaperone folding
    CD_TRANSPORT,       // Vesicle ER→Golgi→membrane
};

struct TranscriptionState {
    float rnaPolPosition;  // 0-1 along gene
    int currentBP;         // base pair index being read
    bool active;
    float timer;
    int transcriptIndex = -1;
    // -1 = transcribing HBB (host default).
    // >= 0 = index into CentralDogmaState::viralGenes (Phase 8A).
    int geneIndex = -1;
};

struct SplicingState {
    float progress;        // 0-1
    bool intron1Removed;
    bool intron2Removed;
    bool complete;
    bool active = false;
    int transcriptIndex = -1;
};

struct TranslationState {
    float ribosomePosition; // 0-1 along mRNA
    int currentCodon;       // codon index being read
    int polypeptideLength;  // amino acids added so far
    bool active;
    float timer;
    int transcriptIndex = -1;
    int trnaIndex = -1;
    float codonProgress = 0.0f;
    bool waitingForTRNA = true;
    bool stopReached = false;
    char incomingAA = '?';
    std::string codon;
    std::string anticodon;
    std::string nascentPeptide;
};

static constexpr int CDOGMA_MAX_TRANSCRIPTS = 6;
static constexpr int CDOGMA_MAX_TRNA_POOL = 12;
static constexpr int CDOGMA_MAX_PROTEINS = 12;

struct TranscriptState {
    bool active = false;
    int uid = 0;
    int sourcePolymerase = -1;
    int genomicBases = 0;
    int ribosomeLoad = 0;
    int spliceosomeIndex = -1;
    float transcriptionProgress = 0.0f;
    float processingProgress = 0.0f;
    float exportProgress = 0.0f;
    float age = 0.0f;
    bool transcriptionComplete = false;
    bool capped = false;
    bool intron1Removed = false;
    bool intron2Removed = false;
    bool spliced = false;
    bool polyadenylated = false;
    bool exported = false;
    std::string preMRNA;
    std::string matureMRNA;
    std::string polyATail;
    // Phase 8A — which gene template produced this transcript.
    // -1 = HBB host default; >= 0 = index into viralGenes.
    int  geneIndex        = -1;
    // Per-transcript start codon, used by translation instead of the
    // cell-wide `startCodonBase` when geneIndex >= 0.
    int  transcriptStartCodon = -1;
    int  transcriptStopCodon  = -1;
    // Tag propagated to ProteinProductState.proteinName on release.
    std::string geneId;
    std::string proteinTag;
};

struct ChargedTRNAState {
    bool active = false;
    bool charged = false;
    int ribosomeIndex = -1;
    float shuttle = 0.0f;
    std::string anticodon;
    char aminoAcid = '?';
};

struct ProteinProductState {
    bool active = false;
    int sourceTranscript = -1;
    float releaseAge = 0.0f;
    float foldProgress = 0.0f;
    float modificationProgress = 0.0f;
    bool folded = false;
    bool mature = false;
    std::string aminoAcids;
    std::string structureAsset;
    // Phase 8A — gene + protein tag so downstream code can route
    // viral proteins into the per-cell `viralProteinCount` ledger
    // (Part 8 §55 Stage 8).
    std::string geneId;
    std::string proteinName;  // e.g. "flu_HA"
    int geneIndex = -1;       // -1 = HBB host protein
};

static constexpr int CDOGMA_MAX_REPL_ORIGINS = 3;
static constexpr int CDOGMA_MAX_REPL_FORKS = 4;
static constexpr int CDOGMA_MAX_REPL_ERRORS = 24;
static constexpr float CDOGMA_LIT_FORK_SPEED_BP_PER_SEC = 24.0f; // ~1.44 kb/min, aligned with mammalian DNA-combing/single-molecule fork-rate measurements.

struct ReplicationOriginState {
    int baseIndex = 0;
    bool licensed = false;
    bool fired = false;
    bool dormant = false;
    bool rescued = false;
    bool passivelyReplicated = false;
    float activation = 0.0f;
    float fireTimer = 0.0f;
};

struct ReplicationForkState {
    bool active = false;
    int originIndex = -1;
    int direction = 1; // -1 left, +1 right
    float basePosition = 0.0f;
    float previousBasePosition = 0.0f;
    float velocityBpPerSec = 0.0f;
    float recruitment = 0.0f;
    float stress = 0.0f;
    float pauseTimer = 0.0f;
    bool proofreading = false;
    bool stalled = false;
};

struct ReplicationMismatchState {
    bool active = false;
    int bpPosition = 0;
    int forkIndex = -1;
    char correctBase = 'N';
    char wrongBase = 'N';
    float sequenceRisk = 0.0f;
    float severity = 0.0f;
    float age = 0.0f;
    bool proofreadAttempted = false;
    bool mmrDetected = false;
    bool repaired = false;
    bool repairedByProofreading = false;
    bool repairedByMMR = false;
    bool escaped = false;
};

struct CentralDogmaState {
    // Multiple concurrent processes
    TranscriptionState transcription[3];  // Up to 3 RNA Pol II
    SplicingState splicing[2];            // Up to 2 spliceosomes
    TranslationState translation[6];      // Up to 6 ribosomes translating
    TranscriptState transcripts[CDOGMA_MAX_TRANSCRIPTS];
    ChargedTRNAState trnaPool[CDOGMA_MAX_TRNA_POOL];
    ProteinProductState proteins[CDOGMA_MAX_PROTEINS];
    int activeMRNAs;                      // mature mRNAs in cytoplasm

    // Phase 8A — extra gene templates available for transcription.
    // Appended at runtime when a virion uncoats (Part 8 §55 Stage 7)
    // or a drug induces heterologous expression. HBB remains the
    // implicit "default" gene consulted when no viral gene is loaded
    // or when the random-weighted pick selects it (promoterStrength
    // of HBB is 1.0 by convention).
    std::vector<GeneTemplate> viralGenes;
    // Aggregate viral-protein counts per proteinName, harvested by
    // Simulation::updatePathogens to drive stoichiometric assembly.
    // This replaces Part 7's float-valued `assembledVirions[]`.
    // Phase P1: swapped from std::unordered_map to std::map so future
    // telemetry / diagnostics that iterate this ledger produce deterministic
    // sorted-key order. N is tiny (<100 viral proteins at high MOI); the
    // log-N insert cost is unmeasurable vs amortised O(1).
    std::map<std::string, int> viralProteinCount;

    float replicationProgress;   // 0-1 during S-phase
    bool replicationActive;
    int nextTranscriptUid;
    std::string codingMRNA;
    int startCodonBase;
    int stopCodonBase;
    ReplicationOriginState replicationOrigins[CDOGMA_MAX_REPL_ORIGINS];
    ReplicationForkState replicationForks[CDOGMA_MAX_REPL_FORKS];
    ReplicationMismatchState replicationErrors[CDOGMA_MAX_REPL_ERRORS];
    bool replicatedBase[HBB_LENGTH];
    bool sPhaseProgramStarted = false;
    float polymeraseRecruitment = 0.0f;
    float dntpAvailability = 1.0f;
    float replicationStress = 0.0f;
    float chk1Signal = 0.0f;
    float replicationQuality = 1.0f;
    float predictedErrorRisk = 0.0f;
    int proofreadCorrections = 0;
    int mmrCorrections = 0;
    int escapedErrors = 0;
    int dormantOriginsFired = 0;
    int unresolvedReplicationErrors = 0;
    int totalReplicationErrors = 0;

    static float logistic(float x) {
        return 1.0f / (1.0f + expf(-x));
    }

    void clearMismatchState(ReplicationMismatchState& err) {
        err = {};
        err.correctBase = 'N';
        err.wrongBase = 'N';
        err.forkIndex = -1;
    }

    void resetReplicationProgram() {
        replicationProgress = 0.0f;
        replicationActive = false;
        sPhaseProgramStarted = false;
        polymeraseRecruitment = 0.0f;
        dntpAvailability = 0.90f;
        replicationStress = 0.0f;
        chk1Signal = 0.0f;
        replicationQuality = 1.0f;
        predictedErrorRisk = 0.0f;
        proofreadCorrections = 0;
        mmrCorrections = 0;
        escapedErrors = 0;
        dormantOriginsFired = 0;
        unresolvedReplicationErrors = 0;
        totalReplicationErrors = 0;
        std::fill(replicatedBase, replicatedBase + HBB_LENGTH, false);
        for (auto& origin : replicationOrigins) origin = {};
        for (auto& fork : replicationForks) fork = {};
        for (auto& err : replicationErrors) clearMismatchState(err);

        replicationOrigins[0].baseIndex = HBB_LENGTH / 2;
        replicationOrigins[0].licensed = true;
        replicationOrigins[0].dormant = false;

        replicationOrigins[1].baseIndex = HBB_LENGTH / 4;
        replicationOrigins[1].licensed = true;
        replicationOrigins[1].dormant = true;

        replicationOrigins[2].baseIndex = (HBB_LENGTH * 3) / 4;
        replicationOrigins[2].licensed = true;
        replicationOrigins[2].dormant = true;
    }

    int countActiveReplicationForks() const {
        int count = 0;
        for (const auto& fork : replicationForks) if (fork.active) count++;
        return count;
    }

    int countFiredOrigins() const {
        int count = 0;
        for (const auto& origin : replicationOrigins) if (origin.fired) count++;
        return count;
    }

    int findFreeForkSlot() const {
        for (int i = 0; i < CDOGMA_MAX_REPL_FORKS; i++) {
            if (!replicationForks[i].active) return i;
        }
        return -1;
    }

    int findFreeReplicationErrorSlot() const {
        for (int i = 0; i < CDOGMA_MAX_REPL_ERRORS; i++) {
            if (!replicationErrors[i].active) return i;
        }
        return -1;
    }

    void markBaseReplicated(int bpIndex) {
        if (bpIndex < 0 || bpIndex >= HBB_LENGTH) return;
        replicatedBase[bpIndex] = true;
    }

    void updateReplicationCoverage() {
        int covered = 0;
        for (int i = 0; i < HBB_LENGTH; i++) {
            if (replicatedBase[i]) covered++;
        }
        replicationProgress = (float)covered / (float)HBB_LENGTH;
    }

    float predictReplicationErrorRisk(int bpIndex, float forkStress) const {
        int run = homopolymerRunLengthAt(HBB_SEQUENCE, HBB_LENGTH, bpIndex);
        float gcFrac = localGCFraction(HBB_SEQUENCE, HBB_LENGTH, bpIndex, 8);
        bool cpgHotspot =
            (bpIndex + 1 < HBB_LENGTH &&
             (HBB_SEQUENCE[bpIndex] == 'C' || HBB_SEQUENCE[bpIndex] == 'c') &&
             (HBB_SEQUENCE[bpIndex + 1] == 'G' || HBB_SEQUENCE[bpIndex + 1] == 'g')) ||
            (bpIndex > 0 &&
             (HBB_SEQUENCE[bpIndex - 1] == 'C' || HBB_SEQUENCE[bpIndex - 1] == 'c') &&
             (HBB_SEQUENCE[bpIndex] == 'G' || HBB_SEQUENCE[bpIndex] == 'g'));
        bool nearBoundary = false;
        for (int i = 0; i < HBB_REGION_COUNT; i++) {
            if (abs(bpIndex - HBB_REGIONS[i].start) <= 2 || abs(bpIndex - HBB_REGIONS[i].end) <= 2) {
                nearBoundary = true;
                break;
            }
        }

        float logit =
            -4.6f +
            0.48f * (float)std::max(run - 1, 0) +
            0.32f * hotspotMultiplier(bpIndex) +
            0.55f * forkStress +
            0.38f * (1.0f - dntpAvailability) +
            0.22f * fabsf(gcFrac - 0.5f) * 2.0f +
            (cpgHotspot ? 0.55f : 0.0f) +
            (nearBoundary ? 0.18f : 0.0f);
        return logistic(logit);
    }

    char transitionBase(char base) const {
        switch (base) {
            case 'A': return 'G';
            case 'G': return 'A';
            case 'C': return 'T';
            case 'T': return 'C';
            default: return 'N';
        }
    }

    char chooseWrongBase(char correctBase, int bpIndex) const {
        bool cpgHotspot =
            (bpIndex + 1 < HBB_LENGTH &&
             correctBase == 'C' &&
             (HBB_SEQUENCE[bpIndex + 1] == 'G' || HBB_SEQUENCE[bpIndex + 1] == 'g'));
        if (cpgHotspot) return 'T';
        if (((float)rand() / RAND_MAX) < 0.72f) return transitionBase(correctBase);
        const char bases[] = {'A', 'T', 'G', 'C'};
        char wrong = correctBase;
        while (wrong == correctBase) wrong = bases[rand() % 4];
        return wrong;
    }

    void enqueueReplicationError(int bpIndex, int forkIndex, float risk, float forkStress) {
        int slot = findFreeReplicationErrorSlot();
        if (slot < 0 || bpIndex < 0 || bpIndex >= HBB_LENGTH) return;
        auto& err = replicationErrors[slot];
        clearMismatchState(err);
        err.active = true;
        err.bpPosition = bpIndex;
        err.forkIndex = forkIndex;
        err.correctBase = HBB_SEQUENCE[bpIndex];
        err.wrongBase = chooseWrongBase(err.correctBase, bpIndex);
        err.sequenceRisk = risk;
        err.severity = clampf(0.25f + risk * 0.9f + forkStress * 0.45f, 0.1f, 1.0f);
        totalReplicationErrors++;
    }

    void fireReplicationOrigin(int originIndex, bool rescue) {
        if (originIndex < 0 || originIndex >= CDOGMA_MAX_REPL_ORIGINS) return;
        auto& origin = replicationOrigins[originIndex];
        if (!origin.licensed || origin.fired || origin.passivelyReplicated) return;

        origin.fired = true;
        origin.rescued = rescue;
        origin.activation = 0.15f;
        origin.fireTimer = 0.0f;
        if (rescue) dormantOriginsFired++;
        markBaseReplicated(origin.baseIndex);

        for (int dir = 0; dir < 2; dir++) {
            int slot = findFreeForkSlot();
            if (slot < 0) break;
            auto& fork = replicationForks[slot];
            fork = {};
            fork.active = true;
            fork.originIndex = originIndex;
            fork.direction = (dir == 0) ? -1 : 1;
            fork.basePosition = (float)origin.baseIndex;
            fork.previousBasePosition = fork.basePosition;
            fork.velocityBpPerSec = 0.0f;
            fork.recruitment = rescue ? 0.45f : 0.22f;
            fork.stress = rescue ? 0.10f : 0.04f;
            fork.pauseTimer = rescue ? 0.10f : 0.0f;
            fork.proofreading = false;
            fork.stalled = false;
        }
    }

    const ReplicationForkState* getReplicationFork(int slot) const {
        if (slot < 0 || slot >= CDOGMA_MAX_REPL_FORKS) return nullptr;
        if (!replicationForks[slot].active) return nullptr;
        return &replicationForks[slot];
    }

    int replicationForkBase(int slot) const {
        const auto* fork = getReplicationFork(slot);
        if (!fork) return -1;
        return std::max(0, std::min(HBB_LENGTH - 1, (int)roundf(fork->basePosition)));
    }

    bool isBaseReplicated(int bpIndex) const {
        if (bpIndex < 0 || bpIndex >= HBB_LENGTH) return false;
        return replicatedBase[bpIndex];
    }

    int replicatedBaseCount() const {
        int count = 0;
        for (int i = 0; i < HBB_LENGTH; i++) {
            if (replicatedBase[i]) count++;
        }
        return count;
    }

    bool replicationReadyForM() const {
        return replicationProgress >= 0.995f &&
               countActiveReplicationForks() == 0 &&
               unresolvedReplicationErrors == 0 &&
               chk1Signal < 0.35f &&
               replicationQuality > 0.82f;
    }

    // ── Phase 8A helpers ─────────────────────────────────────────────
    // Route every sequence read through here so HBB is still a direct
    // pointer-to-literal access (same perf) while any registered viral
    // gene goes through its std::string storage.
    inline const char* geneSequence(int geneIndex, int& outLen) const {
        if (geneIndex < 0 || geneIndex >= (int)viralGenes.size()) {
            outLen = HBB_LENGTH;
            return HBB_SEQUENCE;
        }
        outLen = (int)viralGenes[geneIndex].dnaSeq.size();
        return viralGenes[geneIndex].dnaSeq.c_str();
    }
    inline int geneLength(int geneIndex) const {
        if (geneIndex < 0 || geneIndex >= (int)viralGenes.size())
            return HBB_LENGTH;
        return (int)viralGenes[geneIndex].dnaSeq.size();
    }
    // Pick a gene to start transcribing using promoter-weighted random:
    // HBB_weight = 1.0 by convention. Returns geneIndex (-1 for HBB).
    int pickGeneForTranscription() {
        if (viralGenes.empty()) return -1;
        float total = 1.0f; // HBB weight
        for (const auto& g : viralGenes) total += g.promoterStrength;
        float r = ((float)rand() / (float)RAND_MAX) * total;
        if (r < 1.0f) return -1;
        r -= 1.0f;
        for (int i = 0; i < (int)viralGenes.size(); i++) {
            r -= viralGenes[i].promoterStrength;
            if (r <= 0.0f) return i;
        }
        return (int)viralGenes.size() - 1;
    }
    // Register a new gene template; returns its index.
    int registerViralGene(const GeneTemplate& gene) {
        viralGenes.push_back(gene);
        // Pre-build the mRNA if start/stop codon unknown.
        GeneTemplate& g = viralGenes.back();
        if (g.codingMRNA.empty()) {
            g.codingMRNA.reserve(g.dnaSeq.size());
            for (char c : g.dnaSeq) {
                // Keep viral mRNA simple: no splicing — prokaryote-like
                // one-exon; works for flu / covid / any segmented genome.
                g.codingMRNA.push_back(transcribe(c));
            }
        }
        if (g.startCodon < 0) g.startCodon = findStartCodon(g.codingMRNA);
        if (g.stopCodon < 0 && g.startCodon >= 0)
            g.stopCodon = findStopCodon(g.codingMRNA, g.startCodon);
        return (int)viralGenes.size() - 1;
    }

    void clearTranscriptionState(TranscriptionState& ts) {
        ts = {};
        ts.transcriptIndex = -1;
    }

    void clearSplicingState(SplicingState& sp) {
        sp = {};
        sp.active = false;
        sp.complete = false;
        sp.transcriptIndex = -1;
    }

    void clearTranslationState(TranslationState& tr) {
        tr = {};
        tr.transcriptIndex = -1;
        tr.trnaIndex = -1;
        tr.waitingForTRNA = true;
        tr.incomingAA = '?';
    }

    void clearTranscriptState(TranscriptState& tx) {
        tx = {};
    }

    void clearTRNAState(ChargedTRNAState& trna) {
        trna = {};
        trna.ribosomeIndex = -1;
        trna.aminoAcid = '?';
    }

    void clearProteinState(ProteinProductState& protein) {
        protein = {};
        protein.sourceTranscript = -1;
    }

    TranscriptState* getTranscript(int idx) {
        if (idx < 0 || idx >= CDOGMA_MAX_TRANSCRIPTS) return nullptr;
        if (!transcripts[idx].active) return nullptr;
        return &transcripts[idx];
    }

    const TranscriptState* getTranscript(int idx) const {
        if (idx < 0 || idx >= CDOGMA_MAX_TRANSCRIPTS) return nullptr;
        if (!transcripts[idx].active) return nullptr;
        return &transcripts[idx];
    }

    int totalCodingCodons() const {
        if (startCodonBase < 0 || stopCodonBase <= startCodonBase) return 0;
        return (stopCodonBase - startCodonBase) / 3;
    }

    int countActiveTranscription() const {
        int count = 0;
        for (const auto& ts : transcription) if (ts.active) count++;
        return count;
    }

    int countActiveSplicing() const {
        int count = 0;
        for (const auto& sp : splicing) if (sp.active && !sp.complete) count++;
        return count;
    }

    int countActiveTranslation() const {
        int count = 0;
        for (const auto& tr : translation) if (tr.active) count++;
        return count;
    }

    int countActiveProteins() const {
        int count = 0;
        for (const auto& p : proteins) if (p.active) count++;
        return count;
    }

    int countActiveTRNAs() const {
        int count = 0;
        for (const auto& t : trnaPool) if (t.active) count++;
        return count;
    }

    int findFreeTranscriptSlot() const {
        for (int i = 0; i < CDOGMA_MAX_TRANSCRIPTS; i++) {
            if (!transcripts[i].active) return i;
        }
        return -1;
    }

    int findFreeSpliceosome() const {
        for (int i = 0; i < 2; i++) {
            if (!splicing[i].active) return i;
        }
        return -1;
    }

    int findFreeTRNA() const {
        for (int i = 0; i < CDOGMA_MAX_TRNA_POOL; i++) {
            if (!trnaPool[i].active) return i;
        }
        return -1;
    }

    int findFreeProteinSlot() const {
        for (int i = 0; i < CDOGMA_MAX_PROTEINS; i++) {
            if (!proteins[i].active) return i;
        }
        return -1;
    }

    int findTranscriptForTranslation() const {
        int bestIdx = -1;
        int bestLoad = 9999;
        for (int i = 0; i < CDOGMA_MAX_TRANSCRIPTS; i++) {
            const auto& tx = transcripts[i];
            if (!tx.active || !tx.exported || tx.matureMRNA.empty()) continue;
            if (tx.ribosomeLoad >= 3) continue;
            if (tx.ribosomeLoad < bestLoad) {
                bestLoad = tx.ribosomeLoad;
                bestIdx = i;
            }
        }
        return bestIdx;
    }

    void bindSpliceosomeIfNeeded(int transcriptIndex) {
        auto* tx = getTranscript(transcriptIndex);
        if (!tx || tx->spliced || !tx->transcriptionComplete || tx->spliceosomeIndex >= 0) return;
        int slot = findFreeSpliceosome();
        if (slot < 0) return;
        auto& sp = splicing[slot];
        clearSplicingState(sp);
        sp.active = true;
        sp.transcriptIndex = transcriptIndex;
        tx->spliceosomeIndex = slot;
    }

    bool startTranscription(int polymeraseIndex) {
        int txSlot = findFreeTranscriptSlot();
        if (txSlot < 0) return false;

        clearTranscriptState(transcripts[txSlot]);
        auto& tx = transcripts[txSlot];
        tx.active = true;
        tx.uid = nextTranscriptUid++;
        tx.sourcePolymerase = polymeraseIndex;
        // Phase 8A: pick which gene this Pol II transcribes. Promoter-
        // weighted random across host HBB + registered viral genes.
        tx.geneIndex = pickGeneForTranscription();
        if (tx.geneIndex >= 0 && tx.geneIndex < (int)viralGenes.size()) {
            const GeneTemplate& g = viralGenes[tx.geneIndex];
            tx.geneId = g.id;
            tx.proteinTag = g.proteinName.empty() ? g.id : g.proteinName;
            tx.transcriptStartCodon = g.startCodon;
            tx.transcriptStopCodon  = g.stopCodon;
        } else {
            tx.geneIndex = -1;
            tx.geneId = "HBB";
            tx.proteinTag = "HBB";
            tx.transcriptStartCodon = startCodonBase;
            tx.transcriptStopCodon  = stopCodonBase;
        }

        auto& pol = transcription[polymeraseIndex];
        clearTranscriptionState(pol);
        pol.active = true;
        pol.transcriptIndex = txSlot;
        pol.geneIndex = tx.geneIndex;
        pol.timer = 0.0f;
        pol.rnaPolPosition = 0.0f;
        pol.currentBP = 0;
        return true;
    }

    void primeCurrentCodon(TranslationState& tr) {
        tr.codon.clear();
        tr.anticodon.clear();
        tr.incomingAA = '?';
        tr.stopReached = false;
        tr.codonProgress = 0.0f;
        tr.waitingForTRNA = true;

        const auto* tx = getTranscript(tr.transcriptIndex);
        if (!tx || tx->matureMRNA.empty()) {
            tr.stopReached = true;
            return;
        }
        // Phase 8A: pick per-transcript start codon (viral mRNAs have
        // their own ORF); fall back to cell-wide HBB start for host.
        int txStart = (tx->transcriptStartCodon >= 0)
                      ? tx->transcriptStartCodon : startCodonBase;
        int txStop  = (tx->transcriptStopCodon  >= 0)
                      ? tx->transcriptStopCodon  : stopCodonBase;
        if (txStart < 0) {
            tr.stopReached = true;
            return;
        }

        int base = txStart + tr.currentCodon * 3;
        if (base + 2 >= (int)tx->matureMRNA.size()) {
            tr.stopReached = true;
            return;
        }

        tr.codon = tx->matureMRNA.substr(base, 3);
        AminoAcid aa = translateCodon(tr.codon[0], tr.codon[1], tr.codon[2]);
        if (aa.letter == '*') {
            tr.stopReached = true;
            return;
        }
        tr.anticodon = anticodonForCodon(tr.codon[0], tr.codon[1], tr.codon[2]);
        tr.incomingAA = aa.letter;
        int totalCodons = (txStop > txStart) ? ((txStop - txStart) / 3)
                                              : totalCodingCodons();
        tr.ribosomePosition = (totalCodons > 0) ? ((float)(tr.currentCodon + 1) / (float)totalCodons) : 0.0f;
    }

    bool startTranslation(int ribosomeIndex, int transcriptIndex) {
        auto* tx = getTranscript(transcriptIndex);
        if (!tx || !tx->exported || tx->matureMRNA.empty()) return false;

        auto& tr = translation[ribosomeIndex];
        clearTranslationState(tr);
        tr.active = true;
        tr.transcriptIndex = transcriptIndex;
        tr.timer = 0.0f;
        primeCurrentCodon(tr);
        if (tr.stopReached) {
            clearTranslationState(tr);
            return false;
        }
        tx->ribosomeLoad++;
        return true;
    }

    void releaseTRNA(int trnaIndex) {
        if (trnaIndex < 0 || trnaIndex >= CDOGMA_MAX_TRNA_POOL) return;
        clearTRNAState(trnaPool[trnaIndex]);
    }

    void finishTranslation(int ribosomeIndex) {
        if (ribosomeIndex < 0 || ribosomeIndex >= 6) return;
        auto& tr = translation[ribosomeIndex];
        if (tr.trnaIndex >= 0) releaseTRNA(tr.trnaIndex);
        auto* tx = getTranscript(tr.transcriptIndex);
        if (tx && tx->ribosomeLoad > 0) tx->ribosomeLoad--;
        clearTranslationState(tr);
    }

    void releaseProtein(const TranslationState& tr) {
        if (tr.nascentPeptide.empty()) return;
        int proteinSlot = findFreeProteinSlot();
        if (proteinSlot < 0) return;
        auto& protein = proteins[proteinSlot];
        clearProteinState(protein);
        protein.active = true;
        protein.sourceTranscript = tr.transcriptIndex;
        protein.aminoAcids = tr.nascentPeptide;
        // Phase 8A: propagate gene origin to the protein so callers
        // can route viral peptides into viralProteinCount without
        // sequence-matching the polypeptide.
        const auto* tx = getTranscript(tr.transcriptIndex);
        if (tx) {
            protein.geneIndex = tx->geneIndex;
            protein.geneId = tx->geneId;
            protein.proteinName = tx->proteinTag;
            if (protein.geneIndex >= 0) {
                viralProteinCount[tx->proteinTag]++;
                if (protein.geneIndex < (int)viralGenes.size()) {
                    viralGenes[protein.geneIndex].copiesMade++;
                }
            }
        }
        if ((int)tr.nascentPeptide.size() >= 120 && (int)tr.nascentPeptide.size() <= 170) {
            protein.structureAsset = "hbb_folded.pdb";
        }
    }

    void stepReplicationMismatchRepair(float dt) {
        unresolvedReplicationErrors = 0;
        for (auto& err : replicationErrors) {
            if (!err.active) continue;
            err.age += dt;

            if (!err.proofreadAttempted && err.age >= 0.12f + err.severity * 0.22f) {
                err.proofreadAttempted = true;
                float proofreadChance =
                    clampf(0.985f - err.severity * 0.18f - replicationStress * 0.10f, 0.72f, 0.995f);
                if (((float)rand() / RAND_MAX) < proofreadChance) {
                    err.repaired = true;
                    err.repairedByProofreading = true;
                    proofreadCorrections++;
                    err.active = false;
                    continue;
                }
            }

            if (!err.mmrDetected && err.age >= 0.45f) {
                float detectRate = clampf(0.55f + err.sequenceRisk * 1.8f + replicationStress * 0.4f, 0.15f, 0.98f);
                if (((float)rand() / RAND_MAX) < detectRate * dt * 1.8f) {
                    err.mmrDetected = true;
                }
            }

            if (err.mmrDetected && !err.repaired && err.age >= 1.10f) {
                float repairChance =
                    clampf(0.90f - err.severity * 0.15f - replicationStress * 0.08f, 0.45f, 0.985f);
                if (((float)rand() / RAND_MAX) < repairChance * dt * 1.5f) {
                    err.repaired = true;
                    err.repairedByMMR = true;
                    mmrCorrections++;
                    err.active = false;
                    continue;
                }
            }

            if (err.age >= 4.0f) {
                err.escaped = true;
                escapedErrors++;
                err.active = false;
                continue;
            }

            unresolvedReplicationErrors++;
        }
    }

    void stepReplicationProgram(float dt, int cellPhase) {
        bool inSPhase = (cellPhase == 1);
        bool allowLateCompletion = sPhaseProgramStarted &&
            (replicationProgress < 0.999f ||
             countActiveReplicationForks() > 0 ||
             unresolvedReplicationErrors > 0 ||
             chk1Signal > 0.35f);
        bool activeReplicationWindow = inSPhase || allowLateCompletion;

        if (inSPhase && !sPhaseProgramStarted) {
            resetReplicationProgram();
            sPhaseProgramStarted = true;
            polymeraseRecruitment = 0.18f;
            fireReplicationOrigin(0, false);
        }

        if (!sPhaseProgramStarted) {
            replicationActive = false;
            replicationQuality = 1.0f;
            predictedErrorRisk = 0.0f;
            unresolvedReplicationErrors = 0;
            return;
        }

        for (auto& origin : replicationOrigins) {
            if (replicatedBase[origin.baseIndex]) origin.passivelyReplicated = true;
            if (origin.fired) {
                origin.fireTimer += dt;
                origin.activation = fminf(1.0f, origin.activation + dt * 1.8f);
            }
        }

        int activeForksBefore = countActiveReplicationForks();
        dntpAvailability = clampf(
            dntpAvailability + dt * (0.10f - activeForksBefore * 0.035f - chk1Signal * 0.04f),
            0.62f, 1.05f);
        polymeraseRecruitment = clampf(
            polymeraseRecruitment + dt * (activeReplicationWindow ? 0.42f : -0.18f),
            0.0f, 1.0f);

        if (activeReplicationWindow && activeForksBefore == 0 && replicationProgress < 0.995f) {
            int uncoveredBase = -1;
            for (int bp = 0; bp < HBB_LENGTH; bp++) {
                if (!replicatedBase[bp]) {
                    uncoveredBase = bp;
                    break;
                }
            }
            if (uncoveredBase >= 0) {
                int rescueOrigin = -1;
                for (int i = 1; i < CDOGMA_MAX_REPL_ORIGINS; i++) {
                    if (replicationOrigins[i].licensed && !replicationOrigins[i].fired) {
                        rescueOrigin = i;
                        break;
                    }
                }
                if (rescueOrigin < 0) rescueOrigin = CDOGMA_MAX_REPL_ORIGINS - 1;

                auto& origin = replicationOrigins[rescueOrigin];
                origin.baseIndex = uncoveredBase;
                origin.passivelyReplicated = false;
                origin.fired = false;
                origin.rescued = false;
                origin.activation = 0.0f;
                origin.fireTimer = 0.0f;
                fireReplicationOrigin(rescueOrigin, true);
            }
        }

        if (activeReplicationWindow) {
            for (int i = 1; i < CDOGMA_MAX_REPL_ORIGINS; i++) {
                auto& origin = replicationOrigins[i];
                if (!origin.licensed || origin.fired || origin.passivelyReplicated) continue;

                float uncoveredBias = replicatedBase[origin.baseIndex] ? 0.0f : 0.42f;
                float rescueBias = replicationStress * 0.95f + (1.0f - replicationProgress) * 0.20f + uncoveredBias;
                float chk1Brake = chk1Signal * 0.28f;
                float fireChance = clampf((rescueBias - chk1Brake) * dt * 0.55f, 0.0f, 0.18f);
                if (((float)rand() / RAND_MAX) < fireChance) {
                    fireReplicationOrigin(i, true);
                }
            }
        }

        float latestRisk = 0.0f;
        float stressAccum = 0.0f;
        int stressCount = 0;

        for (int forkIdx = 0; forkIdx < CDOGMA_MAX_REPL_FORKS; forkIdx++) {
            auto& fork = replicationForks[forkIdx];
            if (!fork.active) continue;

            fork.previousBasePosition = fork.basePosition;
            if (fork.pauseTimer > 0.0f) {
                fork.pauseTimer = fmaxf(0.0f, fork.pauseTimer - dt);
                fork.proofreading = true;
            } else {
                fork.proofreading = false;
            }

            float targetRecruitment = polymeraseRecruitment;
            if (replicationOrigins[fork.originIndex].rescued) targetRecruitment = fmaxf(targetRecruitment, 0.70f);
            fork.recruitment = fminf(targetRecruitment, fork.recruitment + dt * 0.8f);

            float baseForkSpeed = CDOGMA_LIT_FORK_SPEED_BP_PER_SEC;
            float throttle = (0.55f + 0.45f * fork.recruitment) *
                             (0.72f + 0.28f * dntpAvailability) *
                             (1.0f - chk1Signal * 0.30f);
            if (fork.pauseTimer > 0.0f) throttle *= 0.10f;
            fork.velocityBpPerSec = baseForkSpeed * throttle;

            float stochasticStress = (((float)rand() / RAND_MAX) - 0.5f) * 0.06f;
            fork.stress = clampf(
                fork.stress + dt * ((1.0f - dntpAvailability) * 0.55f + chk1Signal * 0.20f - 0.18f) + stochasticStress,
                0.0f, 1.0f);
            fork.stalled = fork.stress > 0.72f;
            if (fork.stalled && fork.pauseTimer <= 0.0f) {
                fork.pauseTimer = 0.20f + fork.stress * 0.22f;
                fork.proofreading = true;
            }

            if (activeReplicationWindow) {
                fork.basePosition += (float)fork.direction * fork.velocityBpPerSec * dt;
            }
            fork.basePosition = clampf(fork.basePosition, 0.0f, (float)(HBB_LENGTH - 1));

            int start = (int)floorf(std::min(fork.previousBasePosition, fork.basePosition));
            int end = (int)floorf(std::max(fork.previousBasePosition, fork.basePosition));
            for (int bp = start; bp <= end; bp++) {
                if (bp < 0 || bp >= HBB_LENGTH) continue;
                bool newlyReplicated = !replicatedBase[bp];
                markBaseReplicated(bp);
                if (!newlyReplicated) continue;

                float risk = predictReplicationErrorRisk(bp, fork.stress);
                latestRisk = fmaxf(latestRisk, risk);
                float baseErrorRate = 5.0e-4f * (1.0f + risk * 4.0f);
                if (((float)rand() / RAND_MAX) < baseErrorRate) {
                    enqueueReplicationError(bp, forkIdx, risk, fork.stress);
                    fork.pauseTimer = fmaxf(fork.pauseTimer, 0.12f + risk * 0.35f);
                    fork.proofreading = true;
                }
            }

            int nextBase = (int)roundf(fork.basePosition) + fork.direction;
            if (nextBase < 0 || nextBase >= HBB_LENGTH ||
                (nextBase >= 0 && nextBase < HBB_LENGTH && replicatedBase[nextBase])) {
                fork.active = false;
                fork.velocityBpPerSec = 0.0f;
                fork.proofreading = false;
            } else {
                stressAccum += fork.stress;
                stressCount++;
            }
        }

        stepReplicationMismatchRepair(dt);
        updateReplicationCoverage();
        predictedErrorRisk = latestRisk;

        float avgForkStress = stressCount > 0 ? (stressAccum / (float)stressCount) : 0.0f;
        float unresolvedPenalty = unresolvedReplicationErrors * 0.08f;
        float escapedPenalty = escapedErrors * 0.10f;
        replicationStress = clampf(avgForkStress * 0.7f + unresolvedPenalty + escapedPenalty, 0.0f, 1.0f);
        chk1Signal = clampf(replicationStress * 0.9f + (replicationProgress < 0.98f ? 0.10f : 0.0f), 0.0f, 1.0f);
        replicationQuality = clampf(1.0f - unresolvedPenalty * 0.8f - escapedPenalty - avgForkStress * 0.18f, 0.0f, 1.0f);
        replicationActive = countActiveReplicationForks() > 0;

        // If all forks have terminated and >99% of bases are replicated,
        // force-complete the remaining few bases. In real biology, residual
        // unreplicated gaps are resolved by post-replicative repair before
        // G2 entry. Without this, 3 stuck bases keep the cell in S-phase
        // indefinitely.
        if (countActiveReplicationForks() == 0 && replicationProgress >= 0.99f && replicationProgress < 1.0f) {
            for (int i = 0; i < HBB_LENGTH; i++) {
                replicatedBase[i] = true;
            }
            int covered = 0;
            for (int i = 0; i < HBB_LENGTH; i++) if (replicatedBase[i]) covered++;
            replicationProgress = (float)covered / (float)HBB_LENGTH;
        }

        if (replicationProgress >= 0.999f && countActiveReplicationForks() == 0 && unresolvedReplicationErrors == 0) {
            chk1Signal = fmaxf(0.0f, chk1Signal - dt * 0.8f);
            replicationQuality = fmaxf(replicationQuality, 0.90f);
        }
    }

    void init() {
        activeMRNAs = 0;
        nextTranscriptUid = 1;
        codingMRNA = buildHBBCodingMRNA();
        startCodonBase = findStartCodon(codingMRNA);
        stopCodonBase = findStopCodon(codingMRNA, startCodonBase);
        // Phase 8A: viralGenes + protein counter reset on every init
        // so a fresh cell has no residual infection.
        viralGenes.clear();
        viralProteinCount.clear();

        for (auto& ts : transcription) clearTranscriptionState(ts);
        for (auto& sp : splicing) clearSplicingState(sp);
        for (auto& tr : translation) clearTranslationState(tr);
        for (auto& tx : transcripts) clearTranscriptState(tx);
        for (auto& trna : trnaPool) clearTRNAState(trna);
        for (auto& protein : proteins) clearProteinState(protein);
        resetReplicationProgram();

        // Start one real transcription event immediately so the nucleus is active.
        startTranscription(0);
        transcription[1].timer = 2.0f;
        transcription[2].timer = 4.0f;
    }

    void update(float dt, int cellPhase) {
        bool synthesisAllowed = (cellPhase == 0 || cellPhase == 2);
        bool mitoticPause = (cellPhase == 3);

        // ── DNA Replication — origin firing, fork progression, proofreading,
        //    mismatch repair, and CHK1-mediated rescue origin control ───────
        stepReplicationProgram(dt, cellPhase);

        // ── Transcription (real pre-mRNA growth on the HBB genomic strand) ───
        for (int i = 0; i < 3; i++) {
            auto& ts = transcription[i];
            if (!ts.active) {
                ts.timer += dt;
                if (!mitoticPause && synthesisAllowed && ts.timer >= (3.5f + i * 2.0f)) {
                    if (startTranscription(i)) ts.timer = 0.0f;
                }
                continue;
            }
            if (mitoticPause) continue;

            auto* tx = getTranscript(ts.transcriptIndex);
            if (!tx) {
                clearTranscriptionState(ts);
                continue;
            }

            // Phase 8A: route the sequence read via the transcript's
            // gene template. HBB (geneIndex == -1) still uses the
            // embedded HBB_SEQUENCE literal; viral genes use their
            // registered string. No-op change for HBB.
            int seqLen = 0;
            const char* seqData = geneSequence(tx->geneIndex, seqLen);
            ts.rnaPolPosition = fminf(1.0f, ts.rnaPolPosition + dt * 0.22f);
            ts.currentBP = std::min(seqLen, (int)floorf(ts.rnaPolPosition * (float)seqLen));
            while (tx->genomicBases < ts.currentBP && tx->genomicBases < seqLen) {
                tx->preMRNA.push_back(transcribe(seqData[tx->genomicBases]));
                tx->genomicBases++;
            }
            tx->transcriptionProgress = ts.rnaPolPosition;
            if (!tx->capped && tx->genomicBases >= 18) tx->capped = true;

            if (ts.currentBP >= seqLen || ts.rnaPolPosition >= 1.0f) {
                tx->transcriptionComplete = true;
                // Viral genes skip splicing (prokaryote-like /
                // monocistronic RNA virus): jump straight to mature.
                if (tx->geneIndex >= 0 && tx->geneIndex < (int)viralGenes.size()) {
                    tx->spliced = true;
                    tx->polyadenylated = true;
                    tx->matureMRNA = viralGenes[tx->geneIndex].codingMRNA;
                    tx->polyATail.assign(48, 'A');
                    tx->intron1Removed = true;
                    tx->intron2Removed = true;
                    tx->processingProgress = 1.0f;
                } else {
                    bindSpliceosomeIfNeeded(ts.transcriptIndex);
                }
                clearTranscriptionState(ts);
                ts.timer = 0.0f;
            }
        }

        // ── Splicing / capping / poly(A) maturation ───────────────────
        for (int i = 0; i < CDOGMA_MAX_TRANSCRIPTS; i++) {
            auto& tx = transcripts[i];
            if (!tx.active || !tx.transcriptionComplete || tx.spliced) continue;
            bindSpliceosomeIfNeeded(i);
        }
        for (auto& sp : splicing) {
            if (!sp.active || sp.complete || mitoticPause) continue;
            auto* tx = getTranscript(sp.transcriptIndex);
            if (!tx) {
                clearSplicingState(sp);
                continue;
            }

            sp.progress = fminf(1.0f, sp.progress + dt * 0.45f);
            tx->processingProgress = sp.progress;
            if (sp.progress >= 0.25f) {
                sp.intron1Removed = true;
                tx->intron1Removed = true;
            }
            if (sp.progress >= 0.70f) {
                sp.intron2Removed = true;
                tx->intron2Removed = true;
            }
            if (sp.progress >= 1.0f) {
                sp.complete = true;
                tx->spliced = true;
                tx->polyadenylated = true;
                tx->matureMRNA = codingMRNA;
                tx->polyATail.assign(48, 'A');
                tx->spliceosomeIndex = -1;
                clearSplicingState(sp);
            }
        }

        // ── mRNA export through pores ─────────────────────────────────
        for (auto& tx : transcripts) {
            if (!tx.active) continue;
            tx.age += dt;
            if (tx.spliced && !tx.exported && !mitoticPause) {
                tx.exportProgress = fminf(1.0f, tx.exportProgress + dt * 0.60f);
                if (tx.exportProgress >= 1.0f) tx.exported = true;
            }
        }

        // ── Translation initiation / elongation / termination ─────────
        for (int i = 0; i < 6; i++) {
            auto& tr = translation[i];
            if (!tr.active) {
                tr.timer += dt;
                if (!mitoticPause) {
                    int txIdx = findTranscriptForTranslation();
                    if (txIdx >= 0 && tr.timer >= (0.8f + i * 0.25f)) {
                        startTranslation(i, txIdx);
                    }
                }
                continue;
            }
            if (mitoticPause) continue;

            tr.timer += dt;
            auto* tx = getTranscript(tr.transcriptIndex);
            if (!tx || !tx->exported) {
                finishTranslation(i);
                continue;
            }

            if (tr.stopReached) {
                releaseProtein(tr);
                finishTranslation(i);
                continue;
            }

            if (tr.trnaIndex < 0 && !tr.codon.empty() && tr.incomingAA != '?') {
                int cargoSlot = findFreeTRNA();
                if (cargoSlot >= 0) {
                    auto& cargo = trnaPool[cargoSlot];
                    clearTRNAState(cargo);
                    cargo.active = true;
                    cargo.charged = true;
                    cargo.ribosomeIndex = i;
                    cargo.anticodon = tr.anticodon;
                    cargo.aminoAcid = tr.incomingAA;
                    tr.trnaIndex = cargoSlot;
                }
            }

            if (tr.trnaIndex >= 0) {
                auto& cargo = trnaPool[tr.trnaIndex];
                cargo.shuttle = fminf(1.0f, cargo.shuttle + dt * 8.0f);
                tr.codonProgress = cargo.shuttle;
                if (cargo.shuttle >= 1.0f) {
                    tr.nascentPeptide.push_back(cargo.aminoAcid);
                    tr.polypeptideLength = (int)tr.nascentPeptide.size();
                    releaseTRNA(tr.trnaIndex);
                    tr.trnaIndex = -1;
                    tr.currentCodon++;
                    primeCurrentCodon(tr);
                    if (tr.stopReached) {
                        releaseProtein(tr);
                        finishTranslation(i);
                    }
                }
            }
        }

        // ── Protein folding / simple post-translational maturation ────
        for (auto& protein : proteins) {
            if (!protein.active) continue;
            protein.releaseAge += dt;
            if (!protein.folded) {
                protein.foldProgress = fminf(1.0f, protein.foldProgress + dt * 0.30f);
                protein.folded = protein.foldProgress >= 0.999f;
            } else if (!protein.mature) {
                protein.modificationProgress = fminf(1.0f, protein.modificationProgress + dt * 0.18f);
                protein.mature = protein.modificationProgress >= 0.999f;
            }
            if (protein.releaseAge > 45.0f) {
                clearProteinState(protein);
            }
        }

        // ── Transcript turnover keeps the cytoplasm from filling forever ──
        for (auto& tx : transcripts) {
            if (!tx.active) continue;
            if (tx.exported && tx.ribosomeLoad == 0 && tx.age > 30.0f) {
                clearTranscriptState(tx);
            }
        }

        activeMRNAs = 0;
        for (const auto& tx : transcripts) {
            if (tx.active && tx.exported) activeMRNAs++;
        }
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  Mitosis State Machine — detailed sub-phases
//  Division axis: X (left/right), matching shader cleavage furrow at x=0
// ══════════════════════════════════════════════════════════════════════════

// Kinetochore–microtubule attachment types, Alberts MBoC 7e Ch 17;
// Musacchio 2015 Curr Biol; Cimini 2008 JCS.
enum KinetochoreAttachmentType {
    ATTACH_UNATTACHED = 0,  // neither sister bound — max MCC signal
    ATTACH_MONOTELIC  = 1,  // one sister bound to one pole, other free
    ATTACH_SYNTELIC   = 2,  // both sisters bound to the SAME pole → same daughter
    ATTACH_MEROTELIC  = 3,  // one sister bound to BOTH poles → lagging chromatid
    ATTACH_AMPHITELIC = 4,  // correct: sisters to opposite poles, high tension
};

struct ChromosomeState {
    simd_float3 position;
    simd_float3 sisterPosition;
    simd_float3 axis;
    float condensation;
    float congress;
    float attachmentA;
    float attachmentB;
    float tension;
    float centromereStretch;
    bool biOriented;
    bool separated;
    float hue;

    // ── SAC / error-correction state ────────────────────────────────
    int   attachmentType      = ATTACH_UNATTACHED;
    int   correctionAttempts  = 0;     // times Aurora B has destabilized
    bool  misSegregated       = false; // ended up in wrong daughter
    bool  lagging             = false; // merotelic → lagging in anaphase
};

// ── DNA Replication Errors ──────────────────────────────────────────────
struct ReplicationError {
    int bpPosition;
    char wrongBase;
    char correctBase;
    bool detected;
    bool repaired;
    float timer;
};

static float hotspotMultiplier(int bpIdx) {
    if (bpIdx + 2 >= HBB_LENGTH) return 1.0f;
    char a = HBB_SEQUENCE[bpIdx], b = HBB_SEQUENCE[bpIdx+1];
    // CpG site
    if ((a == 'C' || a == 'c') && (b == 'G' || b == 'g')) return 2.5f;
    // Dinucleotide repeat
    if (a == b) return 3.0f;
    // Trinucleotide repeat
    if (bpIdx + 2 < HBB_LENGTH && a == HBB_SEQUENCE[bpIdx+2]) return 4.0f;
    return 1.0f;
}

struct MitosisState {
    bool active = false;
    MitosisPhase phase = MITO_NONE;
    float phaseTimer = 0;
    float totalProgress = 0;

    // Phase durations in "speed units" (dt passed to update). Tuned so
    // the sum (~5 speed units) matches the real HeLa ~60 bio-min M phase
    // when speed = dt from headless/simulation loop (dt ≈ 0.333 per
    // bio-min → 5 units ≈ 15 bio-min at current tick rate, leaving room
    // for SAC correction overhead to bring real M toward 60 min).
    // Prior sum was 11.5 → cells piled up in M at ~50% of the population.
    static constexpr float PROPHASE_DUR     = 0.80f;
    static constexpr float PROMETAPHASE_DUR = 0.60f;
    static constexpr float METAPHASE_DUR    = 0.80f;
    static constexpr float ANAPHASE_DUR     = 0.40f;
    static constexpr float TELOPHASE_DUR    = 0.60f;
    static constexpr float CYTOKINESIS_DUR  = 0.80f;

    // Chromatin / nucleus
    float chromatinCondensation = 0;
    float nucleusAlpha = 1.0f;
    float nuclearEnvelopeBreakdown = 0.0f;

    // Spindle — poles separate along X
    simd_float3 poleA, poleB;
    float poleSeparation = 0;
    float spindleAssembly = 0.0f;
    float spindleDisassembly = 0.0f;

    // Chromosomes — human diploid karyotype: 22 autosome pairs + XX/XY = 46.
    // After S-phase each chromosome has 2 sister chromatids joined at the
    // centromere, so metaphase plate shows 46 X-shaped chromosomes =
    // 92 chromatids total. Division requires all 46 bi-oriented and
    // aligned (SAC checkpoint).
    static constexpr int NUM_CHROMO = 46;
    ChromosomeState chromosomes[NUM_CHROMO];
    float kinetochoreAttachment = 0.0f;
    float metaphaseAlignment = 0.0f;
    float spindleCheckpoint = 0.0f;
    float chromatidSeparation = 0.0f;
    float checkpointHold = 0.0f;

    // ── SAC biochemistry (Musacchio 2015 Curr Biol; Lara-Gonzalez 2021) ──
    // MCC = Mad2 + BubR1 + Bub3 + Cdc20. Produced by every unattached
    // kinetochore. Inhibits APC/C-Cdc20 → stabilizes Securin + CycB
    // → anaphase cannot begin while ANY kinetochore is still signaling.
    //
    // Rule: mccLevel = fraction of chromosomes producing wait-anaphase
    // signal (not amphitelic and not silenced). Anaphase entry requires
    // mccLevel < 0.02 (effectively all 46 kinetochores correctly attached).
    float mccLevel = 1.0f;            // 1 = full inhibition, 0 = silenced

    // Aurora B (CPC: Aurora B + INCENP + survivin + borealin). Senses low
    // inter-sister tension and phosphorylates outer kinetochore proteins
    // (Ndc80, KNL1) to release incorrect (syntelic/merotelic) attachments.
    // Active at ~0.85 in HeLa; mutations reduce it → chromosome instability.
    float auroraBActivity = 0.85f;

    // SAC persistence tracking. Typical HeLa metaphase lasts ~30 min; if
    // the SAC fails to silence within ~2 h the cell undergoes "mitotic
    // slippage" — CycB gradually leaks, cell exits M without cytokinesis,
    // resulting in a tetraploid G1 cell (often apoptosis-prone via p53).
    // Ref: Brito & Rieder 2006 Curr Biol; Gascoigne & Taylor 2008 Cancer Cell.
    float slippageTimer = 0.0f;
    bool  mitoticSlippage = false;

    // Error accounting (per-division audit fields).
    int unattachedKinetochores = NUM_CHROMO;  // current count, drives mccLevel
    int misSegregationCount    = 0;           // how many chromatids ended up wrong
    int laggingChromatidCount  = 0;           // merotelic → lagging at anaphase

    // Cytokinesis
    float furrowDepth = 0;
    float cellElongation = 0;
    float contractileRingAssembly = 0.0f;

    // ── Organelle duplication ───────────────────────────────────────
    float organelleDuplication = 0; // 0=single, 1=fully duplicated
    float organelleMigration = 0;   // 0=center, 1=at poles
    // X offset for left/right organelle copies — needs to be large so they reach poles
    float organelleOffsetY() const { return organelleMigration * poleSeparation * 0.8f; }

    // ── Particle duplication ───────────────────────────────────────
    bool particlesDuplicated = false; // True after S-phase duplication
    bool dnaCheckpointPassed = false;
    float dnaDuplicationProgress = 0.0f;

    // ── Post-division: show 2 cells for 5 real seconds ─────────────
    float postDivisionTimer = 0;     // Real-time seconds
    int survivingDaughter = 0;       // 0 or 1 (random)
    float fadeAlpha = 1.0f;          // Fading daughter alpha
    simd_float3 daughterA, daughterB; // Positions of 2 daughter cells

    // ── Replication errors (generated during S-phase) ──────────────
    static constexpr int MAX_ERRORS = 16;
    ReplicationError errors[MAX_ERRORS];
    int errorCount = 0;
    float replicationQuality = 1.0f; // 0-1, used by checkpoint engine

    void generateReplicationErrors() {
        errorCount = 0;
        float baseRate = 1e-3f; // Exaggerated for visibility (real: 1e-7)
        for (int bp = 0; bp < HBB_LENGTH && errorCount < MAX_ERRORS; bp++) {
            float rate = baseRate * hotspotMultiplier(bp);
            if (((float)rand() / RAND_MAX) < rate) {
                auto& e = errors[errorCount++];
                e.bpPosition = bp;
                e.correctBase = HBB_SEQUENCE[bp];
                // Pick a wrong base
                const char bases[] = "ATGC";
                do { e.wrongBase = bases[rand() % 4]; } while (e.wrongBase == e.correctBase);
                e.detected = false;
                e.repaired = false;
                e.timer = 0;
            }
        }
        // Calculate quality score
        replicationQuality = 1.0f - (float)errorCount / 10.0f;
        if (replicationQuality < 0) replicationQuality = 0;
    }

    void updateMMR(float dt) {
        // Mismatch repair: each error gets detected over time, then repaired with 99% probability
        for (int i = 0; i < errorCount; i++) {
            auto& e = errors[i];
            if (e.repaired) continue;
            e.timer += dt;
            if (!e.detected && e.timer > 0.5f + ((float)rand()/RAND_MAX) * 2.0f) {
                e.detected = true;
            }
            if (e.detected && !e.repaired && e.timer > 1.5f) {
                e.repaired = ((float)rand() / RAND_MAX) < 0.99f;
                if (!e.repaired) e.timer = 0; // Try again later (but likely won't fix)
            }
        }
        // Update quality
        int unrepaired = 0;
        for (int i = 0; i < errorCount; i++)
            if (!errors[i].repaired) unrepaired++;
        replicationQuality = 1.0f - (float)unrepaired / 5.0f;
        if (replicationQuality < 0) replicationQuality = 0;
    }

    void start(simd_float3 cellCenter, float cellR) {
        storedCellR = cellR;
        active = true;
        phase = MITO_PROPHASE;
        phaseTimer = 0;
        totalProgress = 0;
        chromatinCondensation = 0;
        nucleusAlpha = 1.0f;
        nuclearEnvelopeBreakdown = 0.0f;
        poleSeparation = 0;
        spindleAssembly = 0.0f;
        spindleDisassembly = 0.0f;
        kinetochoreAttachment = 0.0f;
        metaphaseAlignment = 0.0f;
        spindleCheckpoint = 0.0f;
        chromatidSeparation = 0.0f;
        checkpointHold = 0.0f;
        furrowDepth = 0;
        cellElongation = 0;
        contractileRingAssembly = 0.0f;
        organelleDuplication = 0;
        organelleMigration = 0;
        particlesDuplicated = false;
        dnaCheckpointPassed = false;
        dnaDuplicationProgress = 0.0f;
        postDivisionTimer = 0;
        fadeAlpha = 1.0f;
        survivingDaughter = rand() % 2;

        poleA = cellCenter;
        poleB = cellCenter;
        daughterA = cellCenter;
        daughterB = cellCenter;

        for (int i = 0; i < NUM_CHROMO; i++) {
            auto& c = chromosomes[i];
            float a = (float)i / NUM_CHROMO * 2.0f * M_PI;
            float r = cellR * 0.08f;
            // LOCAL space (origin = cell center)
            c.position = {cosf(a) * r, sinf(a * 0.7f) * r * 0.5f, sinf(a) * r};
            c.sisterPosition = c.position;
            c.axis = normalizeOr({0.0f, 1.0f, 0.35f * sinf(a)}, {0, 1, 0});
            c.condensation = 0;
            c.congress = 0.0f;
            c.attachmentA = 0.0f;
            c.attachmentB = 0.0f;
            c.tension = 0.0f;
            c.centromereStretch = 0.0f;
            c.biOriented = false;
            c.separated = false;
            c.hue = (float)i / NUM_CHROMO;
            // SAC state — reset per division.
            c.attachmentType = ATTACH_UNATTACHED;
            c.correctionAttempts = 0;
            c.misSegregated = false;
            c.lagging = false;
        }

        // Reset SAC biochemistry.
        mccLevel = 1.0f;
        auroraBActivity = 0.85f;
        slippageTimer = 0.0f;
        mitoticSlippage = false;
        unattachedKinetochores = NUM_CHROMO;
        misSegregationCount = 0;
        laggingChromatidCount = 0;

        // Generate replication errors from S-phase
        generateReplicationErrors();
    }

    // ── Initial kinetochore attachment distribution ─────────────────
    // Called on entry to prometaphase when the nuclear envelope breaks
    // down and microtubules first contact kinetochores.
    //
    // Literature baseline (normal cells): ~1% mis-segregation/division
    // → ~1 in 100 chromosomes ends up wrong at ANAPHASE after correction.
    // Pre-correction attachment distribution in HeLa (Cimini 2001 JCB):
    //   ~25% already amphitelic, ~25% monotelic, ~15% syntelic,
    //   ~20% merotelic, ~15% still unattached.
    // Aurora B corrects the non-amphi types during prometaphase/metaphase
    // (see auroraBCorrectionTick below).
    void initializeKinetochoreAttachments() {
        auto roll = []() { return (float)rand() / (float)RAND_MAX; };
        unattachedKinetochores = 0;
        for (int i = 0; i < NUM_CHROMO; i++) {
            float r = roll();
            int t;
            if      (r < 0.25f) t = ATTACH_AMPHITELIC;
            else if (r < 0.50f) t = ATTACH_MONOTELIC;
            else if (r < 0.65f) t = ATTACH_SYNTELIC;
            else if (r < 0.85f) t = ATTACH_MEROTELIC;
            else                t = ATTACH_UNATTACHED;
            chromosomes[i].attachmentType = t;
            chromosomes[i].correctionAttempts = 0;
            chromosomes[i].misSegregated = false;
            chromosomes[i].lagging = false;
            if (t != ATTACH_AMPHITELIC) unattachedKinetochores++;
        }
    }

    // ── Aurora B / CPC error correction ─────────────────────────────
    // Aurora B phosphorylates outer kinetochore proteins (Ndc80, KNL1)
    // when tension is low. This destabilizes incorrect attachments
    // (syntelic, merotelic, monotelic) so they can be re-made. Each
    // destabilization has a chance of resolving into the correct
    // amphitelic state. The correction rate scales with auroraBActivity.
    // Source: Lampson & Cheeseman 2011 Trends Cell Biol; Krenn 2015.
    void auroraBCorrectionTick(float dt) {
        auto roll = []() { return (float)rand() / (float)RAND_MAX; };
        int remaining = 0;
        for (int i = 0; i < NUM_CHROMO; i++) {
            auto& c = chromosomes[i];
            if (c.attachmentType == ATTACH_AMPHITELIC) continue;
            // Destabilize rate depends on attachment type (least tension
            // → fastest correction). Values tuned so an HeLa cell with
            // auroraBActivity=0.85 converges to < 1 mis-segregation per
            // division within the compressed sim-time metaphase window.
            // These are ~3× the biological per-second rates because the
            // simulator compresses ~30 min of metaphase into a few sim
            // seconds; slower values would leave cells stuck in SAC.
            float k = 0.0f;
            switch (c.attachmentType) {
                case ATTACH_UNATTACHED: k = 8.00f; break;  // was 3.50
                case ATTACH_MONOTELIC:  k = 6.50f; break;  // was 2.80
                case ATTACH_SYNTELIC:   k = 5.00f; break;  // was 2.20
                case ATTACH_MEROTELIC:  k = 3.00f; break;  // was 1.30
                default: break;
            }
            float p_correct = k * auroraBActivity * dt;
            if (roll() < p_correct) {
                c.correctionAttempts++;
                // Re-roll the attachment after destabilization. With 46
                // chromosomes each needing to simultaneously converge to
                // AMPHITELIC, the success probability per retry must be
                // high (0.85) so cells don't spend metaphase waiting for
                // the last stuck chromosome. Real Aurora B is extremely
                // effective; HeLa mis-segregation rate is ~1 % per
                // division (Thompson & Compton 2008 JCB).
                float r = roll();
                if      (r < 0.85f) c.attachmentType = ATTACH_AMPHITELIC;
                else if (r < 0.92f) c.attachmentType = ATTACH_MONOTELIC;
                else if (r < 0.96f) c.attachmentType = ATTACH_SYNTELIC;
                else                c.attachmentType = ATTACH_MEROTELIC;
            }
            if (chromosomes[i].attachmentType != ATTACH_AMPHITELIC) remaining++;
        }
        unattachedKinetochores = remaining;
        // MCC level proportional to unattached fraction. Floor lowered
        // from 0.15 → 0.03 so SAC silences when only 1–2 chromosomes are
        // still mis-attached (real HeLa SAC has non-zero basal signaling
        // but isn't strict enough to block indefinitely). This prevents
        // the "pile up in M" that was dominating population dynamics.
        float frac = (float)remaining / (float)NUM_CHROMO;
        mccLevel = fmaxf(remaining > 0 ? 0.03f : 0.0f, frac);
    }

    float storedCellR = 2.0f; // Set at start(), tracks actual cell radius

    static float saturate(float v) {
        return fmaxf(0.0f, fminf(1.0f, v));
    }

    static float smooth(float v) {
        v = saturate(v);
        return v * v * (3.0f - 2.0f * v);
    }

    static simd_float3 normalizeOr(simd_float3 v, simd_float3 fallback) {
        float d = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
        if (d < 1.0e-5f) return fallback;
        float invD = 1.0f / d;
        return {v.x * invD, v.y * invD, v.z * invD};
    }

    static simd_float3 lerp3(simd_float3 a, simd_float3 b, float t) {
        return {
            a.x + (b.x - a.x) * t,
            a.y + (b.y - a.y) * t,
            a.z + (b.z - a.z) * t
        };
    }

    // Poles along X, LOCAL space. In real HeLa metaphase the centrosomes
    // sit ~80% of the way from cell center to membrane, giving a spindle
    // that nearly fills the cell. At anaphase the poles continue to push
    // out until they're essentially at the cortex. Coefficient 1.10 ×
    // max poleSeparation 1.2 = 1.32 × cellR; clampInside() at the render
    // site caps this at 0.85 × cellR so poles never pierce the membrane.
    // Ref: Dumont & Mitchison 2009; Kajtez 2016 Nat Commun (HeLa spindle
    // geometry by light-sheet microscopy).
    void updatePoles() {
        float sep = poleSeparation * storedCellR * 1.10f;
        poleA = {sep, 0, 0};
        poleB = {-sep, 0, 0};
    }

    void update(float speed, simd_float3 cellCenter, float actualCellR) {
        if (!active) return;
        phaseTimer += speed;
        totalProgress += speed;
        storedCellR = actualCellR; // Track growing cell
        float cellR = actualCellR; // Use REAL radius, not hardcoded
        (void)cellCenter;

        auto chromosomeAxisFor = [&](int idx, float twist) -> simd_float3 {
            float a = ((float)idx / (float)NUM_CHROMO) * 2.0f * M_PI + twist;
            simd_float3 raw = {
                0.10f * sinf(a * 1.3f),
                1.0f,
                0.45f * cosf(a * 0.9f)
            };
            return normalizeOr(raw, {0, 1, 0});
        };
        // Distribute 46 chromosomes across the metaphase plate as a
        // sunflower / golden-ratio disk tiling (Vogel 1979) so they fill
        // the equatorial cross-section evenly rather than collapsing into
        // a single ring. Real metaphase plate in HeLa has chromosomes
        // scattered across a disk ~8 µm wide (Jaqaman 2010 JCB).
        auto platePosFor = [&](int idx, float radius) -> simd_float3 {
            const float PHI = 2.39996323f; // golden angle (radians)
            float a = (float)idx * PHI;
            // sqrt(v) gives uniform disk density
            float v = ((float)idx + 0.5f) / (float)NUM_CHROMO;
            float r = sqrtf(v) * radius;
            // Tiny time-varying wobble so chromosomes jitter with
            // cytoplasmic flow instead of sitting perfectly still.
            float wobble = sinf(totalProgress * 1.3f + (float)idx * 0.73f) * radius * 0.04f;
            return {
                wobble,
                sinf(a) * r,
                cosf(a) * r
            };
        };
        auto chromosomeCenter = [&](const ChromosomeState& c) -> simd_float3 {
            return {
                (c.position.x + c.sisterPosition.x) * 0.5f,
                (c.position.y + c.sisterPosition.y) * 0.5f,
                (c.position.z + c.sisterPosition.z) * 0.5f
            };
        };
        auto setSisterPair = [&](ChromosomeState& c, simd_float3 center, simd_float3 axis, float gap) {
            simd_float3 dir = normalizeOr(axis, {0, 1, 0});
            simd_float3 halfGap = {dir.x * gap * 0.5f, dir.y * gap * 0.5f, dir.z * gap * 0.5f};
            c.axis = dir;
            c.position = {center.x + halfGap.x, center.y + halfGap.y, center.z + halfGap.z};
            c.sisterPosition = {center.x - halfGap.x, center.y - halfGap.y, center.z - halfGap.z};
        };

        switch (phase) {
        case MITO_PROPHASE: {
            float p = smooth(phaseTimer / PROPHASE_DUR);
            chromatinCondensation = 0.72f * p;
            nucleusAlpha = 1.0f - 0.18f * p;
            nuclearEnvelopeBreakdown = 0.0f;
            poleSeparation = p * 0.42f;
            spindleAssembly = p * 0.30f;
            spindleDisassembly = 0.0f;
            organelleDuplication = p;
            organelleMigration = 0.0f;
            kinetochoreAttachment = 0.0f;
            metaphaseAlignment = 0.0f;
            spindleCheckpoint = 0.0f;
            chromatidSeparation = 0.0f;
            checkpointHold = 0.0f;
            contractileRingAssembly = 0.0f;
            furrowDepth = 0.0f;
            cellElongation = 0.02f * p;
            updatePoles();
            for (int i = 0; i < NUM_CHROMO; i++) {
                auto& c = chromosomes[i];
                float a = ((float)i / (float)NUM_CHROMO) * 2.0f * M_PI + totalProgress * 0.18f;
                float orbitR = cellR * (0.08f + 0.02f * sinf(totalProgress * 0.9f + i));
                simd_float3 center = {
                    cosf(a) * orbitR,
                    sinf(a * 1.2f) * orbitR * 0.45f,
                    sinf(a) * orbitR
                };
                c.condensation = chromatinCondensation;
                c.congress = 0.0f;
                c.attachmentA = 0.0f;
                c.attachmentB = 0.0f;
                c.tension = 0.0f;
                c.centromereStretch = 0.0f;
                c.biOriented = false;
                c.separated = false;
                setSisterPair(c, center, chromosomeAxisFor(i, totalProgress * 0.25f),
                              cellR * (0.015f + 0.020f * p));
            }
            updateMMR(speed);
            if (p >= 1.0f && dnaCheckpointPassed) {
                phase = MITO_PROMETAPHASE;
                phaseTimer = 0;
                // Nuclear envelope breakdown exposes kinetochores to
                // microtubules — initial attachment types are rolled now.
                initializeKinetochoreAttachments();
            } else if (p >= 1.0f && !dnaCheckpointPassed) {
                // Hold at the prophase checkpoint until both daughter genomes exist.
                phaseTimer = PROPHASE_DUR;
                chromatinCondensation = 1.0f;
                nucleusAlpha = fmaxf(nucleusAlpha, 0.05f);
            }
            break;
        }
        case MITO_PROMETAPHASE: {
            float pRaw = saturate(phaseTimer / PROMETAPHASE_DUR);
            float p = smooth(pRaw);
            chromatinCondensation = 0.72f + 0.28f * p;
            nucleusAlpha = 1.0f - p;
            nuclearEnvelopeBreakdown = p;
            poleSeparation = 0.42f + p * 0.34f;
            spindleAssembly = 0.30f + 0.70f * p;
            spindleDisassembly = 0.0f;
            organelleMigration = p * 0.35f;
            chromatidSeparation = 0.0f;
            contractileRingAssembly = 0.0f;
            furrowDepth = 0.0f;
            cellElongation = 0.03f + 0.03f * p;
            updatePoles();
            float attachSum = 0.0f;
            float alignSum = 0.0f;
            for (int i = 0; i < NUM_CHROMO; i++) {
                auto& c = chromosomes[i];
                // Per-chromosome stochastic stagger — each chromosome hits
                // its MT search-and-capture at a slightly different time.
                // Scale the index offset by 1/NUM_CHROMO so the stagger
                // spans [0, ~0.4] regardless of chromosome count (the old
                // formula used `i * 0.05` which starved most chromosomes
                // of attachment once NUM_CHROMO grew past ~10).
                float frac = (float)i / (float)NUM_CHROMO;
                float biasA = saturate(pRaw * 1.25f - frac * 0.35f);
                float biasB = saturate(pRaw * 1.20f - (1.0f - frac) * 0.35f);
                c.attachmentA = fmaxf(c.attachmentA, 0.18f + 0.82f * biasA);
                c.attachmentB = fmaxf(c.attachmentB, 0.15f + 0.85f * biasB);
                c.biOriented = (c.attachmentA > 0.40f && c.attachmentB > 0.40f);
                c.congress = fminf(1.0f, c.congress + speed * (0.24f + 0.60f * (c.attachmentA + c.attachmentB) * 0.5f));
                // Plate radius widens as alignment progresses — chromosomes
                // spread across the full equatorial disk (real HeLa plate
                // diameter ≈ 8 µm = 0.40 × cellR at cellR = 2 wu).
                simd_float3 targetCenter = platePosFor(i, cellR * 0.38f);
                simd_float3 center = lerp3(chromosomeCenter(c), targetCenter, 0.08f + 0.08f * c.congress);
                c.tension = saturate((fminf(c.attachmentA, c.attachmentB) - 0.25f) / 0.65f);
                c.centromereStretch = c.tension;
                setSisterPair(c, center, chromosomeAxisFor(i, 0.35f + totalProgress * 0.12f),
                              cellR * (0.022f + 0.028f * c.tension));
                c.condensation = chromatinCondensation;
                attachSum += fminf(c.attachmentA, c.attachmentB);
                alignSum += 1.0f - saturate(fabsf(center.x) / fmaxf(cellR * 0.10f, 0.001f));
            }
            kinetochoreAttachment = attachSum / (float)NUM_CHROMO;
            metaphaseAlignment = alignSum / (float)NUM_CHROMO;
            spindleCheckpoint = fminf(kinetochoreAttachment, metaphaseAlignment);
            // Run Aurora B error correction each tick during prometaphase.
            // This progressively converts syntelic/merotelic/monotelic
            // attachments into amphitelic ones, reducing mccLevel.
            auroraBCorrectionTick(speed);
            if (p >= 1.0f && spindleCheckpoint > 0.55f) {
                phase = MITO_METAPHASE;
                phaseTimer = 0;
            }
            break;
        }
        case MITO_METAPHASE: {
            float pRaw = saturate(phaseTimer / METAPHASE_DUR);
            float p = smooth(pRaw);
            chromatinCondensation = 1.0f;
            nucleusAlpha = 0.0f;
            nuclearEnvelopeBreakdown = 1.0f;
            poleSeparation = 0.76f;
            spindleAssembly = 1.0f;
            spindleDisassembly = 0.0f;
            organelleMigration = 0.35f + p * 0.25f;
            chromatidSeparation = 0.08f + 0.10f * p;
            contractileRingAssembly = 0.0f;
            furrowDepth = 0.0f;
            cellElongation = 0.06f;
            updatePoles();
            float attachSum = 0.0f;
            float alignSum = 0.0f;
            float tensionSum = 0.0f;
            for (int i = 0; i < NUM_CHROMO; i++) {
                auto& c = chromosomes[i];
                simd_float3 platePos = platePosFor(i, cellR * 0.40f);
                simd_float3 center = lerp3(chromosomeCenter(c), platePos, 0.12f + 0.10f * p);
                c.attachmentA = fminf(1.0f, c.attachmentA + speed * 0.90f);
                c.attachmentB = fminf(1.0f, c.attachmentB + speed * 0.90f);
                c.tension = fminf(1.0f, c.tension + speed * 0.85f);
                c.centromereStretch = 0.25f + 0.75f * c.tension;
                c.biOriented = true;
                c.congress = 1.0f;
                c.separated = false;
                setSisterPair(c, center, chromosomeAxisFor(i, 0.70f),
                              cellR * (0.030f + 0.035f * c.tension));
                attachSum += fminf(c.attachmentA, c.attachmentB);
                alignSum += 1.0f - saturate(fabsf(center.x) / fmaxf(cellR * 0.06f, 0.001f));
                tensionSum += c.tension;
            }
            kinetochoreAttachment = attachSum / (float)NUM_CHROMO;
            metaphaseAlignment = alignSum / (float)NUM_CHROMO;
            float avgTension = tensionSum / (float)NUM_CHROMO;
            spindleCheckpoint = fminf(kinetochoreAttachment, fminf(metaphaseAlignment, avgTension));

            // Aurora B keeps correcting during metaphase too — any
            // residual syntelic/merotelic attachments are reworked until
            // the MCC signal falls below threshold.
            auroraBCorrectionTick(speed);

            // Mitotic slippage: if the SAC stays active for very long
            // (CycB slowly leaks despite APC/C inhibition), the cell
            // exits M-phase without proper division → tetraploid G1.
            // Real HeLa: metaphase ~30 min, slippage after ~2-6 h.
            // We compress that ratio into sim time: slippage after the
            // cell has been stuck at mccLevel > 0.05 for ~25× the normal
            // metaphase duration. Most normal cycles should SAC-silence
            // first via Aurora B correction.
            if (mccLevel > 0.05f) slippageTimer += speed;
            else slippageTimer = 0.0f;
            // Slippage: tightened from 25× → 4× metaphase duration so
            // cells that can't silence SAC don't pile up for hours. HeLa
            // biologically slips in 2–6 h of stalled metaphase; at our
            // compressed sim timescale 4× METAPHASE_DUR (~8 bio-min) is
            // the practical window before we let the cell slip.
            if (slippageTimer > METAPHASE_DUR * 4.0f) {
                mitoticSlippage = true;
            }

            if (spindleCheckpoint > 0.96f) checkpointHold += speed;
            else checkpointHold = fmaxf(0.0f, checkpointHold - speed * 0.5f);

            // Anaphase gate — TWO paths out of metaphase:
            //  1. SAC silenced properly: mccLevel low enough AND alignment
            //     good AND minimum metaphase duration met → clean division.
            //  2. Mitotic slippage: cell gives up after extended wait.
            //     Division still proceeds but daughters inherit errors.
            // mccLevel < 0.10 is the practical threshold — real cells have
            // residual weak MCC signal even at anaphase onset (Dick &
            // Gerlich 2013 Nat Cell Biol).
            // SAC silencing thresholds loosened to match the lower MCC
            // floor (0.03). Real metaphase-to-anaphase transition triggers
            // when MCC sig drops below the APC activation threshold —
            // does NOT require perfect attachment of every chromosome.
            bool sacSilenced = (mccLevel < 0.05f &&
                                spindleCheckpoint > 0.70f &&
                                phaseTimer > METAPHASE_DUR * 0.25f);
            if (sacSilenced || mitoticSlippage) {
                phase = MITO_ANAPHASE;
                phaseTimer = 0;
                // Tally segregation errors: any chromosome still not
                // amphitelic at anaphase entry segregates incorrectly.
                //  - syntelic → both sisters go to same daughter
                //  - merotelic → one sister is a "lagging chromatid"
                //    (Cimini 2001 JCB), may end up in either daughter
                //    or trapped in a micronucleus
                //  - monotelic/unattached → random partition (50/50)
                misSegregationCount = 0;
                laggingChromatidCount = 0;
                for (int i = 0; i < NUM_CHROMO; i++) {
                    auto& c = chromosomes[i];
                    if (c.attachmentType == ATTACH_AMPHITELIC) continue;
                    c.misSegregated = true;
                    misSegregationCount++;
                    if (c.attachmentType == ATTACH_MEROTELIC) {
                        c.lagging = true;
                        laggingChromatidCount++;
                    }
                }
            }
            break;
        }
        case MITO_ANAPHASE: {
            float p = smooth(phaseTimer / ANAPHASE_DUR);
            chromatinCondensation = 1.0f;
            nucleusAlpha = 0.0f;
            nuclearEnvelopeBreakdown = 1.0f;
            poleSeparation = 0.76f + p * 0.44f;
            spindleAssembly = 1.0f;
            spindleDisassembly = 0.0f;
            chromatidSeparation = p;
            cellElongation = 0.08f + p * 0.28f;
            furrowDepth = 0.05f + p * 0.18f;
            contractileRingAssembly = 0.25f + p * 0.30f;
            organelleMigration = 0.60f + p * 0.30f;
            updatePoles();
            for (int i = 0; i < NUM_CHROMO; i++) {
                auto& c = chromosomes[i];
                c.separated = true;
                c.biOriented = true;
                c.tension = fmaxf(0.25f, 1.0f - p * 0.70f);
                c.centromereStretch = 1.0f - p * 0.65f;
                float angle = ((float)i / (float)NUM_CHROMO) * 2.0f * M_PI;
                float poleR = cellR * 0.10f;
                simd_float3 tA = {poleA.x + cosf(angle) * poleR * 0.35f,
                                  sinf(angle) * poleR * 0.45f,
                                  poleA.z + cosf(angle * 0.7f) * poleR};
                simd_float3 tB = {poleB.x + cosf(angle) * poleR * 0.35f,
                                  sinf(angle) * poleR * 0.45f,
                                  poleB.z + cosf(angle * 0.7f) * poleR};
                float spd = 0.10f + p * 0.12f;
                c.position = lerp3(c.position, tA, spd);
                c.sisterPosition = lerp3(c.sisterPosition, tB, spd);
            }
            if (p >= 1.0f) { phase = MITO_TELOPHASE; phaseTimer = 0; }
            break;
        }
        case MITO_TELOPHASE: {
            float p = smooth(phaseTimer / TELOPHASE_DUR);
            chromatinCondensation = 1.0f - p * 0.65f;
            nucleusAlpha = 0.15f + 0.85f * p;
            nuclearEnvelopeBreakdown = 1.0f - p;
            poleSeparation = 1.20f - p * 0.45f;
            spindleAssembly = 1.0f - p * 0.55f;
            spindleDisassembly = p;
            organelleMigration = 0.90f + p * 0.10f;
            furrowDepth = 0.23f + p * 0.37f;
            contractileRingAssembly = 0.55f + p * 0.25f;
            cellElongation = 0.36f - p * 0.10f;
            updatePoles();
            for (int i = 0; i < NUM_CHROMO; i++) {
                auto& c = chromosomes[i];
                float angle = ((float)i / (float)NUM_CHROMO) * 2.0f * M_PI;
                float nucR = cellR * 0.09f;
                simd_float3 tA = {poleA.x * 0.88f,
                                  sinf(angle) * nucR * 0.45f,
                                  cosf(angle) * nucR};
                simd_float3 tB = {poleB.x * 0.88f,
                                  sinf(angle) * nucR * 0.45f,
                                  cosf(angle) * nucR};
                c.position = lerp3(c.position, tA, 0.09f);
                c.sisterPosition = lerp3(c.sisterPosition, tB, 0.09f);
                c.condensation = chromatinCondensation;
                c.tension = fmaxf(0.0f, c.tension - speed * 0.5f);
                c.centromereStretch = fmaxf(0.0f, c.centromereStretch - speed * 0.6f);
            }
            if (p >= 1.0f) { phase = MITO_CYTOKINESIS; phaseTimer = 0; }
            break;
        }
        case MITO_CYTOKINESIS: {
            float p = smooth(phaseTimer / CYTOKINESIS_DUR);
            chromatinCondensation = 0.35f * (1.0f - p);
            nucleusAlpha = 1.0f;
            nuclearEnvelopeBreakdown = 0.0f;
            furrowDepth = 0.60f + p * 0.40f;
            cellElongation = 0.26f - p * 0.08f;
            organelleMigration = 1.0f;
            spindleAssembly = 0.45f * (1.0f - p);
            spindleDisassembly = 1.0f;
            contractileRingAssembly = 1.0f;
            chromatidSeparation = 1.0f;
            poleSeparation = 0.75f * (1.0f - p);
            updatePoles();
            for (int i = 0; i < NUM_CHROMO; i++) {
                auto& c = chromosomes[i];
                c.condensation = chromatinCondensation;
                c.tension = 0.0f;
                c.centromereStretch = 0.0f;
            }
            if (p >= 1.0f) {
                float daughterOffsetX = cellR * 0.35f;
                daughterA = { daughterOffsetX, 0, 0 };
                daughterB = { -daughterOffsetX, 0, 0 };
                phase = MITO_COMPLETE;
                phaseTimer = 0;
                postDivisionTimer = 0;
            }
            break;
        }
        case MITO_COMPLETE: {
            // Post-division: phaseTimer accumulates at mitosis speed, but we
            // want 5 REAL seconds. Since speed ≈ 0.0033/frame at 60fps,
            // and we call update() each frame, use a separate real-time counter.
            // postDivisionTimer is incremented externally with real dt in main.mm
            spindleAssembly = 0.0f;
            spindleDisassembly = 1.0f;
            contractileRingAssembly = 1.0f;
            if (postDivisionTimer > 3.5f) {
                fadeAlpha = fmaxf(0, 1.0f - (postDivisionTimer - 3.5f) * 1.5f);
            }
            break;
        }
        default: break;
        }
    }

    bool postDivisionComplete() const {
        // Finalize atomically the frame cytokinesis ends. No artificial delay:
        // in real biology the two cells exist the instant the furrow closes.
        return phase == MITO_COMPLETE;
    }
};
