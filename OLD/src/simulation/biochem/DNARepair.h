#pragma once
#include <cmath>
#include <vector>
#include <algorithm>

// ══════════════════════════════════════════════════════════════════════════
//  DNARepair.h — DNA replication, damage detection, and repair pathways
//  Based on: Alberts Ch 5 "DNA Replication, Repair, and Recombination"
//
//  Implements:
//    - Replication fork machinery:
//      Helicase (MCM2-7), SSB (RPA), Primase, DNA Pol δ/ε,
//      Sliding clamp (PCNA), Clamp loader (RFC),
//      Leading strand (continuous), Lagging strand (Okazaki fragments),
//      Ligase (nick sealing), RNase H1/FEN1 (primer removal)
//    - Proofreading (3'→5' exonuclease, error rate 10^-7)
//    - Mismatch Repair (MMR): MutSα(MSH2-MSH6), MutLα(MLH1-PMS2)
//    - Base Excision Repair (BER): glycosylase → AP endo → Pol β → ligase
//    - Nucleotide Excision Repair (NER): global genome + transcription-coupled
//    - Double-Strand Break Repair:
//      NHEJ: Ku70/80 → DNA-PKcs → Artemis → XRCC4-Ligase IV
//      HR: MRN → CtIP → RPA → Rad51 → strand invasion → resolution
//    - DNA Damage Response (DDR):
//      ATM (DSBs) / ATR (stalled forks) → Chk1/Chk2 → p53/p21 → arrest
//    - Telomerase (TERT + TERC, telomere maintenance)
//    - Origin licensing: ORC → Cdc6 → Cdt1 → MCM loading (pre-RC)
//
//  Damage types modeled:
//    - Spontaneous depurination/deamination (~10,000/cell/day)
//    - Oxidative damage (8-oxoguanine, from ROS)
//    - UV damage (cyclobutane pyrimidine dimers, 6-4 photoproducts)
//    - Alkylation (O6-methylguanine)
//    - Single-strand breaks
//    - Double-strand breaks (most dangerous)
//    - Replication errors (mismatches, insertions, deletions)
//
//  Units: counts (number of lesions), rates per second
// ══════════════════════════════════════════════════════════════════════════

// ── Damage Types ────────────────────────────────────────────────────────
enum DamageType {
    DMG_NONE = 0,
    DMG_DEPURINATION,      // AP site (abasic)
    DMG_DEAMINATION,       // C→U or 5mC→T
    DMG_OXIDATION,         // 8-oxoG, thymine glycol
    DMG_ALKYLATION,        // O6-methylguanine
    DMG_UV_CPD,            // Cyclobutane pyrimidine dimer
    DMG_UV_64PP,           // 6-4 photoproduct
    DMG_MISMATCH,          // Replication error (wrong base)
    DMG_INSERTION,         // Small insertion
    DMG_DELETION,          // Small deletion
    DMG_SSB,               // Single-strand break
    DMG_DSB,               // Double-strand break
    DMG_CROSSLINK,         // Interstrand crosslink
    DMG_COUNT
};

struct DNALesion {
    DamageType type;
    int position;          // Base pair position in genome
    int chromosome;        // Which chromosome (0-22 + X/Y)
    float severity;        // 0-1 (affects repair priority)
    bool detected;         // Has repair machinery found it?
    bool repaired;         // Successfully repaired?
    float repair_timer;    // Time since detection
};

// ── Replication Fork State ──────────────────────────────────────────────
struct ReplicationFork {
    bool active = false;
    int origin_position = 0;           // Origin of replication
    int fork_position_leading = 0;     // Leading strand position
    int fork_position_lagging = 0;     // Lagging strand position
    float speed = 1000.0f;             // bp/s (human: ~1000 bp/s per fork)
    int okazaki_count = 0;             // Number of Okazaki fragments made
    int errors_made = 0;               // Replication errors
    bool stalled = false;              // Fork stalled (damage ahead)

    // Proteins at the fork
    float helicase_MCM = 1.0f;         // MCM2-7 helicase activity
    float primase = 1.0f;              // Primase activity
    float polDelta = 1.0f;             // DNA Pol δ (lagging)
    float polEpsilon = 1.0f;           // DNA Pol ε (leading)
    float PCNA = 1.0f;                 // Sliding clamp loaded
    float RFC = 1.0f;                  // Clamp loader
    float RPA = 1.0f;                  // Single-strand binding protein
    float ligase = 1.0f;               // DNA Ligase I
    float FEN1 = 1.0f;                 // Flap endonuclease 1
    float RNaseH = 1.0f;               // RNase H1 (primer removal)
};

// ── DNA Damage Response Proteins ────────────────────────────────────────
struct DDRState {
    // Sensor kinases
    float ATM_active = 0.0f;           // Activated by DSBs
    float ATR_active = 0.0f;           // Activated by stalled forks/ssDNA

    // Effector kinases
    float Chk1_active = 0.0f;          // Downstream of ATR
    float Chk2_active = 0.0f;          // Downstream of ATM

    // CDC25 phosphatases (cell cycle drivers)
    float CDC25A = 1.0f;               // Activates CDK2 (S-phase entry)
    float CDC25C = 1.0f;               // Activates CDK1 (M-phase entry)

    // p53 pathway
    float p53_level = 0.1f;            // p53 protein level
    float p53_P = 0.0f;                // Phosphorylated (stabilized) p53
    float MDM2 = 0.5f;                 // MDM2 (ubiquitin ligase, degrades p53)
    float p21_level = 0.0f;            // p21/CDKN1A (CDK inhibitor)

    // γ-H2AX foci (DSB markers)
    float gammaH2AX = 0.0f;            // Phosphorylated H2AX (DSB marker)

    // Checkpoint status
    bool G1S_arrested = false;
    bool intraS_arrested = false;
    bool G2M_arrested = false;
};

// ── Repair Pathway Activities ───────────────────────────────────────────
struct RepairCapacity {
    // Mismatch Repair (MMR)
    float MSH2 = 1.0f;                // MutSα component (mismatch recognition)
    float MSH6 = 1.0f;                // MutSα component
    float MLH1 = 1.0f;                // MutLα component (strand discrimination)
    float PMS2 = 1.0f;                // MutLα component
    float MMR_rate = 0.99f;           // Fraction of mismatches repaired

    // Base Excision Repair (BER)
    float OGG1 = 1.0f;                // 8-oxoG glycosylase
    float UNG = 1.0f;                 // Uracil-DNA glycosylase
    float APE1 = 1.0f;                // AP endonuclease 1
    float PolBeta = 1.0f;             // DNA Pol β (gap filling)
    float XRCC1 = 1.0f;              // Scaffold protein
    float BER_rate = 0.95f;

    // Nucleotide Excision Repair (NER)
    float XPC = 1.0f;                 // Global genome NER sensor
    float CSB = 1.0f;                 // Transcription-coupled NER
    float XPA = 1.0f;                 // Damage verification
    float XPF_ERCC1 = 1.0f;          // 5' incision
    float XPG = 1.0f;                 // 3' incision
    float NER_rate = 0.90f;

    // Non-Homologous End Joining (NHEJ)
    float Ku70 = 1.0f;               // DSB end binding
    float Ku80 = 1.0f;
    float DNA_PKcs = 1.0f;           // Kinase
    float Artemis = 1.0f;            // End processing nuclease
    float XRCC4 = 1.0f;              // Ligase scaffold
    float LigaseIV = 1.0f;           // Ligase
    float NHEJ_rate = 0.85f;         // Error-prone!

    // Homologous Recombination (HR)
    float MRE11 = 1.0f;              // MRN complex (end resection)
    float RAD50 = 1.0f;
    float NBS1 = 1.0f;
    float CtIP = 1.0f;               // End resection
    float EXO1 = 1.0f;               // Extended resection
    float RPA_repair = 1.0f;         // ssDNA protection
    float RAD51 = 1.0f;              // Strand invasion (RecA homolog)
    float BRCA1 = 1.0f;              // HR pathway choice
    float BRCA2 = 1.0f;              // RAD51 loader
    float HR_rate = 0.99f;           // High fidelity (needs sister chromatid)

    // Direct reversal
    float MGMT = 0.5f;               // O6-methylguanine-DNA methyltransferase
};

// ── Telomere State ──────────────────────────────────────────────────────
struct TelomereState {
    float length = 10000.0f;          // Base pairs (human: 5000-15000)
    float loss_per_division = 50.0f;  // bp lost per replication (end problem)
    float telomerase_TERT = 0.0f;     // Telomerase reverse transcriptase
    float telomerase_TERC = 0.0f;     // Telomerase RNA component
    float shelterin = 1.0f;           // Shelterin complex (protects ends)
    bool critically_short = false;    // Triggers senescence/crisis

    void replicateOnce() {
        float added = telomerase_TERT * telomerase_TERC * 100.0f; // Telomerase adds ~100bp
        length = length - loss_per_division + added;
        if (length < 4000.0f) critically_short = true;
    }
};

// ── Origin Licensing ────────────────────────────────────────────────────
struct OriginLicensing {
    float ORC = 1.0f;                 // Origin Recognition Complex (bound to origins)
    float Cdc6 = 0.0f;               // Loaded in late M/G1
    float Cdt1 = 0.0f;               // MCM loader
    float MCM_loaded = 0.0f;         // MCM2-7 helicase loaded (pre-RC formed)
    float geminin = 0.0f;            // CDK-dependent inhibitor of Cdt1
    bool licensed = false;            // Origins licensed for replication?

    void licenseOrigins(float cdk_activity) {
        // Licensing only occurs in G1 (low CDK)
        if (cdk_activity < 0.2f) {
            Cdc6 += 0.1f * (1.0f - Cdc6);
            Cdt1 += 0.1f * (1.0f - Cdt1) * (1.0f - geminin);
            MCM_loaded += 0.05f * Cdc6 * Cdt1 * ORC * (1.0f - MCM_loaded);
        } else {
            // High CDK prevents re-licensing (prevents re-replication)
            geminin = cdk_activity;
            Cdt1 *= 0.95f;
        }
        licensed = MCM_loaded > 0.8f;
    }
};

// ══════════════════════════════════════════════════════════════════════════
//  DNA REPAIR ENGINE
// ══════════════════════════════════════════════════════════════════════════

class DNARepairEngine {
public:
    // Current damage state
    std::vector<DNALesion> lesions;
    int total_damage_count[DMG_COUNT] = {};
    int repaired_count[DMG_COUNT] = {};

    // Sub-systems
    DDRState ddr;
    RepairCapacity repair;
    TelomereState telomeres;
    OriginLicensing origins;
    std::vector<ReplicationFork> forks;

    // Summary stats
    float total_damage = 0;
    float total_unrepaired = 0;
    float mutation_rate = 0;     // Permanent mutations per division
    float genome_integrity = 1.0f; // 0-1

    void init() {
        lesions.clear();
        ddr = DDRState();
        repair = RepairCapacity();
        telomeres = TelomereState();
        origins = OriginLicensing();
        forks.clear();
        memset(total_damage_count, 0, sizeof(total_damage_count));
        memset(repaired_count, 0, sizeof(repaired_count));
    }

    // ── Generate spontaneous damage ─────────────────────────────────
    // Ref: Lindahl 1993 Nature 362:709. Rates calibrated to lesions/cell/day.
    // Converted to lesions/second = X / 86400.
    void generateSpontaneousDamage(float dt, float ROS_level, float UV_exposure) {
        // Depurination: ~5000/day = 0.058/s
        if (randFloat() < 0.058f * dt) addLesion(DMG_DEPURINATION);
        // Deamination C→U: ~100-500/day ≈ 0.0035/s
        if (randFloat() < 0.0035f * dt) addLesion(DMG_DEAMINATION);
        // Oxidative (8-oxoG): ~1000/day under normal ROS, scales with ROS_level
        if (randFloat() < 0.0116f * ROS_level * dt) addLesion(DMG_OXIDATION);
        // UV (CPD): dose-dependent, input is photons/cell/s
        if (UV_exposure > 0 && randFloat() < UV_exposure * dt) addLesion(DMG_UV_CPD);
        // SSBs: ~10,000/day ≈ 0.116/s
        if (randFloat() < 0.116f * dt) addLesion(DMG_SSB);
        // DSBs: ~10-50/day ≈ 0.0003/s (rare but catastrophic)
        if (randFloat() < 0.0003f * dt) addLesion(DMG_DSB);
    }

    // ── Run repair pathways ─────────────────────────────────────────
    void step(float dt, float ROS_level, float UV_exposure, float cdk_activity,
              int cell_phase) {
        generateSpontaneousDamage(dt, ROS_level, UV_exposure);

        // ── Detect and repair lesions ───────────────────────────────
        int dsb_count = 0, ssb_count = 0, base_damage = 0;
        for (auto& lesion : lesions) {
            if (lesion.repaired) continue;
            lesion.repair_timer += dt;

            switch (lesion.type) {
            case DMG_DEPURINATION:
            case DMG_DEAMINATION:
            case DMG_OXIDATION:
            case DMG_ALKYLATION:
                // BER pathway
                if (!lesion.detected && randFloat() < repair.BER_rate * dt * 0.5f) {
                    lesion.detected = true;
                }
                if (lesion.detected && randFloat() < repair.PolBeta * 0.3f * dt) {
                    lesion.repaired = true;
                    repaired_count[lesion.type]++;
                }
                base_damage++;
                break;

            case DMG_UV_CPD:
            case DMG_UV_64PP:
                // NER pathway
                if (!lesion.detected && randFloat() < repair.XPC * 0.2f * dt) {
                    lesion.detected = true;
                }
                if (lesion.detected && randFloat() < repair.NER_rate * 0.1f * dt) {
                    lesion.repaired = true;
                    repaired_count[lesion.type]++;
                }
                base_damage++;
                break;

            case DMG_MISMATCH:
            case DMG_INSERTION:
            case DMG_DELETION:
                // MMR pathway
                if (!lesion.detected && randFloat() < repair.MSH2 * repair.MSH6 * 0.5f * dt) {
                    lesion.detected = true;
                }
                if (lesion.detected && randFloat() < repair.MMR_rate * 0.3f * dt) {
                    lesion.repaired = true;
                    repaired_count[lesion.type]++;
                }
                break;

            case DMG_SSB:
                ssb_count++;
                // Quick repair via BER/SSBR
                if (randFloat() < 0.5f * dt) {
                    lesion.repaired = true;
                    repaired_count[DMG_SSB]++;
                }
                break;

            case DMG_DSB:
                dsb_count++;
                // Choose pathway based on cell cycle phase
                if (cell_phase == 1 || cell_phase == 2) {
                    // S/G2: prefer HR (accurate, uses sister chromatid)
                    if (!lesion.detected) {
                        lesion.detected = true;
                    }
                    if (randFloat() < repair.RAD51 * repair.BRCA2 * repair.HR_rate * 0.05f * dt) {
                        lesion.repaired = true;
                        repaired_count[DMG_DSB]++;
                    }
                } else {
                    // G1: use NHEJ (fast but error-prone)
                    if (randFloat() < repair.Ku70 * repair.NHEJ_rate * 0.1f * dt) {
                        lesion.repaired = true;
                        repaired_count[DMG_DSB]++;
                        // NHEJ may cause small deletion (mutagenic)
                        if (randFloat() < 0.1f) mutation_rate += 0.001f;
                    }
                }
                break;

            default:
                break;
            }
        }

        // Remove repaired lesions periodically
        if (lesions.size() > 100) {
            lesions.erase(std::remove_if(lesions.begin(), lesions.end(),
                [](const DNALesion& l) { return l.repaired; }), lesions.end());
        }

        // ── DNA Damage Response (DDR) ──────────────────────────────
        // ATM activated by DSBs
        ddr.ATM_active = fminf(1.0f, dsb_count * 0.2f);
        // ATR activated by stalled forks or extensive ssDNA
        ddr.ATR_active = fminf(1.0f, ssb_count * 0.01f + base_damage * 0.001f);

        // Chk1/Chk2 activation
        ddr.Chk1_active += (0.5f * ddr.ATR_active - 0.1f * ddr.Chk1_active) * dt;
        ddr.Chk2_active += (0.5f * ddr.ATM_active - 0.1f * ddr.Chk2_active) * dt;

        // Chk1/2 inhibit CDC25 → cell cycle arrest
        ddr.CDC25A = 1.0f / (1.0f + ddr.Chk1_active * 5.0f + ddr.Chk2_active * 3.0f);
        ddr.CDC25C = 1.0f / (1.0f + ddr.Chk1_active * 3.0f + ddr.Chk2_active * 5.0f);

        // p53 stabilization
        // ATM/Chk2 phosphorylate p53 (prevents MDM2 binding → stabilizes p53)
        float p53_stabilization = ddr.ATM_active * 0.5f + ddr.Chk2_active * 0.3f;
        ddr.MDM2 = 0.5f / (1.0f + p53_stabilization * 5.0f); // MDM2 activity reduced
        float p53_deg = ddr.MDM2 * ddr.p53_level; // MDM2 ubiquitinates p53
        float p53_syn = 0.1f + p53_stabilization * 0.5f;
        ddr.p53_level += (p53_syn - p53_deg) * dt;
        ddr.p53_P = ddr.p53_level * p53_stabilization / (0.5f + p53_stabilization);

        // p53 induces p21 (CDK inhibitor → G1/S arrest)
        ddr.p21_level += (0.5f * ddr.p53_P - 0.05f * ddr.p21_level) * dt;

        // γ-H2AX (DSB marker)
        ddr.gammaH2AX = fminf(1.0f, dsb_count * 0.1f);

        // Checkpoint status
        ddr.G1S_arrested = ddr.p21_level > 0.5f || ddr.CDC25A < 0.3f;
        ddr.intraS_arrested = ddr.Chk1_active > 0.5f;
        ddr.G2M_arrested = ddr.CDC25C < 0.3f;

        // ── Origin licensing ───────────────────────────────────────
        origins.licenseOrigins(cdk_activity);

        // ── Genome integrity score ─────────────────────────────────
        total_damage = 0;
        total_unrepaired = 0;
        for (auto& l : lesions) {
            if (!l.repaired) {
                total_unrepaired++;
                total_damage += l.severity;
            }
        }
        genome_integrity = 1.0f / (1.0f + total_unrepaired * 0.01f);
    }

    // ── Accessors ───────────────────────────────────────────────────
    float getDamageLevel() const { return total_unrepaired * 0.01f; }
    float getP53level() const { return ddr.p53_P; }
    float getP21level() const { return ddr.p21_level; }
    bool isCheckpointArrested(int phase) const {
        if (phase == 0) return ddr.G1S_arrested;
        if (phase == 1) return ddr.intraS_arrested;
        if (phase == 2) return ddr.G2M_arrested;
        return false;
    }

private:
    void addLesion(DamageType type) {
        DNALesion l;
        l.type = type;
        l.position = rand() % 3200000000; // Random genome position
        l.chromosome = rand() % 23;
        l.severity = (type == DMG_DSB) ? 1.0f : (type == DMG_SSB ? 0.3f : 0.1f);
        l.detected = false;
        l.repaired = false;
        l.repair_timer = 0;
        lesions.push_back(l);
        total_damage_count[type]++;
    }

    float randFloat() { return (float)rand() / (float)RAND_MAX; }
};
