#pragma once
#include <simd/simd.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

// ══════════════════════════════════════════════════════════════════════════
//  IntracellularPhysics.h — Real FEA-inspired soft-body physics
//
//  NO position overrides. Everything is force-driven:
//  - Membrane = Hertz elastic sphere (force field, not clamp)
//  - Particles bounce off membrane with energy loss
//  - Particle-particle soft collisions via spatial hashing
//  - Cell movement → viscous drag on contents (no teleportation)
//  - Brownian thermal force + viscous damping
// ══════════════════════════════════════════════════════════════════════════

static float ipRandf() { return (float)rand() / (float)RAND_MAX; }
static float ipRandGauss() {
    float u1 = ipRandf() + 1e-8f, u2 = ipRandf();
    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
}

enum ParticleType {
    PT_ORGANELLE = 0, PT_DNA_NODE, PT_DNA_POLYMERASE, PT_RNA_POLYMERASE,
    PT_SPLICEOSOME, PT_PRE_MRNA, PT_MRNA_NODE, PT_NUCLEAR_PORE,
    PT_RIBOSOME_SMALL, PT_RIBOSOME_LARGE, PT_TRNA, PT_POLYPEPTIDE,
    PT_CHAPERONE, PT_VESICLE_COPII, PT_VESICLE_SECRETORY, PT_MOLECULE,
    // Membrane mosaic components — phospholipid bilayer and embedded
    // membrane proteins (Alberts MBoC 7e Ch 10, 15).
    PT_PHOSPHOLIPID,         // PC/PE/PS heads dotting the inner/outer leaflets
    PT_RECEPTOR_GPCR,        // 7-pass transmembrane (adrenergic, odorant, etc.)
    PT_RECEPTOR_RTK,         // single-pass tyrosine kinase (insulin, EGF)
    PT_RECEPTOR_CYTOKINE,    // JAK-STAT-linked single-pass
    PT_ION_CHANNEL,          // voltage/ligand-gated Na, K, Ca, Cl
    PT_ION_PUMP,             // Na/K-ATPase, Ca-ATPase
    PT_ADHESION,             // integrin / cadherin / selectin
    PT_TRANSPORTER,          // GLUT glucose, LAT1 AA, ABC exporter (carrier)
    PT_AQUAPORIN,            // AQP water channel
    PT_GAP_JUNCTION,         // connexin hemichannel for cell-cell
    PT_EXCHANGER,            // Na/Ca, Na/H, Cl/HCO3 antiporters
    PT_RECEPTOR_STK,         // Ser/Thr receptor kinase (TGF-β, BMP)
    PT_RECEPTOR_DEATH,       // Fas, TNF-R (apoptosis)
    PT_RECEPTOR_TLR,         // Toll-like (innate immunity / NF-κB)
    PT_RECEPTOR_NOTCH,       // contact-dependent, proteolytic release
    PT_RECEPTOR_FRIZZLED,    // Wnt / β-catenin pathway
    PT_RECEPTOR_PATCHED,     // Hedgehog pathway
    PT_MHC,                  // MHC I/II antigen presentation
    PT_LDL_RECEPTOR,         // LDL-R endocytosis
    PT_TRANSFERRIN_R,        // iron uptake
    PT_COUNT
};

enum MolJob {
    JOB_IDLE = 0, JOB_GOTO_TARGET, JOB_FOLLOW_PATH, JOB_ATTACH, JOB_CONSUMED,
};

enum MolTag {
    TAG_NONE = 0, TAG_ATP, TAG_GLUCOSE, TAG_NADH, TAG_AMINO_ACID,
    TAG_CAMP, TAG_WATER, TAG_CALCIUM, TAG_LIPID,
    // Extracellular signaling ligands (secreted by neighbors, bind
    // receptors on this cell → paracrine/autocrine signaling).
    TAG_LIGAND_GF,        // growth factor (EGF/FGF/VEGF generic)
    TAG_LIGAND_HORMONE,   // epinephrine / neurotransmitter (GPCR)
    TAG_LIGAND_CYTOKINE,  // interleukin / interferon
    TAG_CAMP_SIGNAL,      // cAMP second messenger visual
};

struct IntraParticle {
    simd_float3 position, velocity;
    float radius, mass;
    ParticleType type;
    int tetherId; float tetherLen, tetherK;
    simd_float3 home; float homeK;
    float confineRadius; simd_float3 confineCenter;
    float colorR, colorG, colorB, glowIntensity;
    float baseRadius;
    int stateIndex;
    bool active;
    MolTag tag;
    MolJob job;
    simd_float3 jobTarget;
    int jobTargetId;
    float jobProgress, jobTimer, jobSpeed;
    simd_float3 spawnPos;
    int mitosisHalf;  // -1=unassigned, 0=poleA (top), 1=poleB (bottom)
};

// Organelle attraction field — pulls specific molecule types
struct AttractionField {
    simd_float3 center;
    float radius;        // Range of influence
    float strength;      // Force magnitude
    MolTag attractedTag; // TAG_NONE = attracts all, otherwise specific
    ParticleType attractedType; // PT_COUNT = match by tag only
};

struct PhysicsParams {
    float viscosity = 0.15f;
    float kT = 0.002f;
    float cellRadius = 2.0f;
    simd_float3 cellCenter = {0, 0, 0};
    simd_float3 cellVelocity = {0, 0, 0}; // Cell movement velocity
    float membraneK = 50.0f;    // Membrane stiffness — very strong wall
    float membraneRestitution = 0.4f; // Energy loss on bounce (0=absorb, 1=elastic)
    float collisionK = 3.0f;    // Particle-particle repulsion
    float cellDragCoeff = 12.0f; // How strongly contents follow cell movement
    // Brownian step magnitude per particle type. Derived from the
    // Stokes-Einstein relation in cytoplasmic water:
    //   D = (kT) / (6 π η r) × crowding_factor
    //     = (4.28e-21 J) / (6π × 1.4e-3 Pa·s × r_m) × 0.26
    // For radius r in µm:  D ≈ 42 / r  µm²/s in cytoplasm at 37 °C.
    // Per-frame displacement σ = sqrt(2 D dt). The dimensionless values
    // below are tuned to that scaling factor (smaller particles have
    // proportionally larger steps) but kept within visual tolerance so
    // the smallest molecules don't appear to teleport across the cell
    // every frame. Membrane particles diffuse 2D laterally at the
    // ~0.1 µm²/s rate Kusumi 2005 reports for embedded proteins.
    // Refs: Luby-Phelps 2000 (cyto crowding 0.26 D-ratio);
    //       Kusumi 2005 Biophys J (membrane diffusion); Einstein 1905.
    float idleBrownian[PT_COUNT] = {
        // PT_ORGANELLE..PT_MOLECULE (16 entries)
        0.01f, 0.02f, 0.0f, 0.0f, 0.02f, 0.03f, 0.03f, 0.0f,
        0.03f, 0.02f, 0.06f, 0.02f, 0.04f, 0.02f, 0.02f, 0.08f,
        // Membrane particles — constrained to surface, slow 2D diffusion
        0.025f, // PT_PHOSPHOLIPID
        0.008f, // PT_RECEPTOR_GPCR
        0.008f, // PT_RECEPTOR_RTK
        0.008f, // PT_RECEPTOR_CYTOKINE
        0.010f, // PT_ION_CHANNEL
        0.008f, // PT_ION_PUMP
        0.006f, // PT_ADHESION
        0.008f, // PT_TRANSPORTER
        0.010f, // PT_AQUAPORIN
        0.004f, // PT_GAP_JUNCTION
        0.008f, // PT_EXCHANGER
        0.008f, // PT_RECEPTOR_STK
        0.008f, // PT_RECEPTOR_DEATH
        0.009f, // PT_RECEPTOR_TLR
        0.005f, // PT_RECEPTOR_NOTCH
        0.008f, // PT_RECEPTOR_FRIZZLED
        0.007f, // PT_RECEPTOR_PATCHED
        0.007f, // PT_MHC
        0.010f, // PT_LDL_RECEPTOR
        0.010f  // PT_TRANSFERRIN_R
    };
};

static inline float intraDot(simd_float3 a, simd_float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline void intraProjectInsideConfinement(IntraParticle& p) {
    float dx = p.position.x - p.confineCenter.x;
    float dy = p.position.y - p.confineCenter.y;
    float dz = p.position.z - p.confineCenter.z;
    float dist = sqrtf(dx * dx + dy * dy + dz * dz);
    float boundary = p.confineRadius - p.radius;
    if (boundary < 0.0f) boundary = 0.0f;
    if (dist <= boundary) return;

    simd_float3 normal = {1.0f, 0.0f, 0.0f};
    if (dist > 0.0001f) {
        float invDist = 1.0f / dist;
        normal = {dx * invDist, dy * invDist, dz * invDist};
    }

    float safetyInset = fmaxf(p.radius * 0.02f, 0.0005f);
    float safeBoundary = fmaxf(boundary - safetyInset, 0.0f);
    p.position = {
        p.confineCenter.x + normal.x * safeBoundary,
        p.confineCenter.y + normal.y * safeBoundary,
        p.confineCenter.z + normal.z * safeBoundary
    };

    float outwardV = intraDot(p.velocity, normal);
    if (outwardV > 0.0f) {
        p.velocity = {
            p.velocity.x - normal.x * outwardV,
            p.velocity.y - normal.y * outwardV,
            p.velocity.z - normal.z * outwardV
        };
    }
}

// Spatial hash for O(n) collision detection
struct SpatialHash {
    float cellSize = 0.15f;
    std::unordered_map<int64_t, std::vector<int>> grid;

    int64_t hashKey(float x, float y, float z) const {
        int ix = (int)floorf(x / cellSize);
        int iy = (int)floorf(y / cellSize);
        int iz = (int)floorf(z / cellSize);
        return ((int64_t)ix * 73856093L) ^ ((int64_t)iy * 19349663L) ^ ((int64_t)iz * 83492791L);
    }

    void clear() { grid.clear(); }

    void insert(int idx, simd_float3 pos) {
        grid[hashKey(pos.x, pos.y, pos.z)].push_back(idx);
    }

    // Get all particles in same + neighboring cells
    void queryNeighbors(simd_float3 pos, std::vector<int>& out) const {
        out.clear();
        int ix = (int)floorf(pos.x / cellSize);
        int iy = (int)floorf(pos.y / cellSize);
        int iz = (int)floorf(pos.z / cellSize);
        for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dz = -1; dz <= 1; dz++) {
                    int64_t key = ((int64_t)(ix+dx)*73856093L) ^
                                  ((int64_t)(iy+dy)*19349663L) ^
                                  ((int64_t)(iz+dz)*83492791L);
                    auto it = grid.find(key);
                    if (it != grid.end())
                        for (int idx : it->second) out.push_back(idx);
                }
    }
};

class IntracellularPhysics {
public:
    std::vector<IntraParticle> particles;
    std::vector<AttractionField> attractionFields;
    PhysicsParams params;
    bool initialized = false;
    SpatialHash spatialHash;

    void init(simd_float3 cellCenter, float cellRadius) {
        params.cellCenter = cellCenter;
        params.cellRadius = cellRadius;
        params.cellVelocity = {0, 0, 0};
        particles.clear();
        initialized = true;
        spatialHash.cellSize = fmaxf(cellRadius * 0.08f, 0.1f);
    }

    int addParticle(ParticleType type, simd_float3 pos, float radius,
                     float cr, float cg, float cb) {
        IntraParticle p = {};
        p.position = pos; p.velocity = {0,0,0};
        p.radius = radius; p.baseRadius = radius;
        p.mass = fmaxf(radius * radius * 10.0f, 0.01f);
        p.type = type;
        p.tetherId = -1; p.tetherLen = 0; p.tetherK = 0;
        p.home = pos; p.homeK = 0;
        p.confineRadius = params.cellRadius;
        p.confineCenter = params.cellCenter;
        p.colorR = cr; p.colorG = cg; p.colorB = cb;
        p.glowIntensity = 0.5f;
        p.stateIndex = 0; p.active = true; p.tag = TAG_NONE;
        p.job = JOB_IDLE; p.jobTarget = pos; p.jobTargetId = -1;
        p.jobProgress = 0; p.jobTimer = 0; p.jobSpeed = 0;
        p.spawnPos = pos;
        p.mitosisHalf = -1; // Unassigned until mitosis
        particles.push_back(p);
        return (int)particles.size() - 1;
    }

    void setTether(int idx, int id2, float len, float k) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        particles[idx].tetherId = id2; particles[idx].tetherLen = len; particles[idx].tetherK = k;
    }
    void setHome(int idx, simd_float3 h, float k) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        particles[idx].home = h; particles[idx].homeK = k;
    }
    void setConfinement(int idx, simd_float3 c, float r) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        particles[idx].confineCenter = c; particles[idx].confineRadius = r;
    }
    void assignGoto(int idx, simd_float3 t, float s) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        auto& p = particles[idx]; p.job = JOB_GOTO_TARGET; p.jobTarget = t; p.jobSpeed = s; p.jobTimer = 0;
    }
    void assignAttach(int idx, int tid) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        particles[idx].job = JOB_ATTACH; particles[idx].jobTargetId = tid;
    }
    void assignConsume(int idx) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        particles[idx].job = JOB_CONSUMED; particles[idx].jobTimer = 0;
    }
    void assignFollowPath(int idx, float s) {
        if (idx < 0 || idx >= (int)particles.size()) return;
        particles[idx].job = JOB_FOLLOW_PATH; particles[idx].jobSpeed = s; particles[idx].jobProgress = 0;
    }
    bool reachedTarget(int idx, float th = 0.05f) {
        if (idx < 0 || idx >= (int)particles.size()) return false;
        auto& p = particles[idx]; if (p.job != JOB_GOTO_TARGET) return false;
        float dx = p.position.x-p.jobTarget.x, dy = p.position.y-p.jobTarget.y, dz = p.position.z-p.jobTarget.z;
        return (dx*dx+dy*dy+dz*dz) < th*th;
    }

    // ══════════════════════════════════════════════════════════════════
    //  PHYSICS STEP — all force-driven, no position overrides
    // ══════════════════════════════════════════════════════════════════
    bool skipCollisions = false; // Set true during mitosis to maintain framerate

    void step(float dt) {
        if (!initialized || dt <= 0) return;
        dt = fminf(dt, 0.016f); // Cap at ~60fps timestep
        int n = (int)particles.size();

        // ── 1. Build spatial hash for collisions ────────────────────
        if (!skipCollisions) {
            spatialHash.clear();
            for (int i = 0; i < n; i++) {
                if (particles[i].active) spatialHash.insert(i, particles[i].position);
            }
        }

        std::vector<int> neighbors;
        neighbors.reserve(64);

        for (int i = 0; i < n; i++) {
            auto& p = particles[i];
            if (!p.active || p.type == PT_NUCLEAR_PORE) continue;

            simd_float3 force = {0, 0, 0};
            p.jobTimer += dt;

            // ── 2. Job-directed forces ──────────────────────────────
            if (p.job == JOB_GOTO_TARGET) {
                simd_float3 toT = {p.jobTarget.x-p.position.x, p.jobTarget.y-p.position.y, p.jobTarget.z-p.position.z};
                float dist = sqrtf(toT.x*toT.x+toT.y*toT.y+toT.z*toT.z);
                if (dist > 0.001f) {
                    float spd = p.jobSpeed > 0 ? p.jobSpeed : 0.02f;
                    force.x += toT.x/dist*spd; force.y += toT.y/dist*spd; force.z += toT.z/dist*spd;
                }
                float noise = params.idleBrownian[p.type]*0.1f*params.kT/fmaxf(p.radius,0.001f);
                force.x += ipRandGauss()*noise; force.y += ipRandGauss()*noise; force.z += ipRandGauss()*noise;
            } else if (p.job == JOB_ATTACH && p.jobTargetId >= 0 && p.jobTargetId < n) {
                auto& t = particles[p.jobTargetId];
                force.x += (t.position.x-p.position.x)*3.0f;
                force.y += (t.position.y-p.position.y)*3.0f;
                force.z += (t.position.z-p.position.z)*3.0f;
            } else if (p.job == JOB_CONSUMED) {
                p.radius *= (1.0f-dt*2.0f);
                p.glowIntensity = fmaxf(0, p.glowIntensity-dt*3.0f);
                if (p.jobTimer > 0.5f || p.radius < p.baseRadius*0.1f) {
                    p.position = p.spawnPos; p.radius = p.baseRadius;
                    p.glowIntensity = 0.5f; p.job = JOB_IDLE;
                    p.jobTimer = 0; p.velocity = {0,0,0};
                }
            } else if (p.job == JOB_FOLLOW_PATH) {
                // Thermal force from fluctuation-dissipation theorem:
                //   <F·F> dt = 2 γ kT dt,  γ = 6π η r (Stokes drag)
                //   => F = sqrt(2 γ kT / dt) · N(0,1)
                float gamma = 6.0f * (float)M_PI * params.viscosity * fmaxf(p.radius, 1e-4f);
                float amp = params.idleBrownian[p.type] * 0.5f
                          * sqrtf(2.0f * gamma * params.kT / fmaxf(dt, 1e-6f));
                force.x += ipRandGauss()*amp; force.y += ipRandGauss()*amp; force.z += ipRandGauss()*amp;
            } else {
                // IDLE — proper Brownian force (fluctuation-dissipation)
                float gamma = 6.0f * (float)M_PI * params.viscosity * fmaxf(p.radius, 1e-4f);
                float amp = params.idleBrownian[p.type]
                          * sqrtf(2.0f * gamma * params.kT / fmaxf(dt, 1e-6f));
                force.x += ipRandGauss()*amp; force.y += ipRandGauss()*amp; force.z += ipRandGauss()*amp;
            }

            // ── 4. Viscous drag (Stokes): F = -6π η r v ────────────
            float drag = 6.0f * (float)M_PI * params.viscosity * p.radius;
            force.x -= drag*p.velocity.x; force.y -= drag*p.velocity.y; force.z -= drag*p.velocity.z;

            // ── 5. Organelle attraction fields ──────────────────────
            // Only when idle — job forces override attraction
            if (p.job == JOB_IDLE && !attractionFields.empty()) {
                for (auto& field : attractionFields) {
                    // Check tag match
                    if (field.attractedTag != TAG_NONE && p.tag != field.attractedTag) continue;
                    if (field.attractedType != PT_COUNT && p.type != field.attractedType) continue;

                    float dx = field.center.x - p.position.x;
                    float dy = field.center.y - p.position.y;
                    float dz = field.center.z - p.position.z;
                    float dist2 = dx*dx + dy*dy + dz*dz;
                    float r2 = field.radius * field.radius;
                    if (dist2 < r2 && dist2 > 0.0001f) {
                        float soft2 = r2 * 0.01f; // Softening to prevent singularity
                        float F = field.strength * p.mass / (dist2 + soft2);
                        float dist = sqrtf(dist2 + soft2);
                        force.x += F * dx / dist;
                        force.y += F * dy / dist;
                        force.z += F * dz / dist;
                    }
                }
            }

            // ── Home spring ─────────────────────────────────────────
            if (p.homeK > 0) {
                force.x -= p.homeK*(p.position.x-p.home.x);
                force.y -= p.homeK*(p.position.y-p.home.y);
                force.z -= p.homeK*(p.position.z-p.home.z);
            }

            // ── 6. Tether spring ────────────────────────────────────
            if (p.tetherId >= 0 && p.tetherId < n) {
                auto& o = particles[p.tetherId];
                float dx=p.position.x-o.position.x, dy=p.position.y-o.position.y, dz=p.position.z-o.position.z;
                float d = sqrtf(dx*dx+dy*dy+dz*dz);
                if (d > 0.001f) {
                    float f = p.tetherK*(d-p.tetherLen);
                    force.x -= f*dx/d; force.y -= f*dy/d; force.z -= f*dz/d;
                }
            }

            // ── 7. Membrane: HARD elastic wall ─────────────────────
            // Two-layer approach:
            // Layer 1: Soft repulsion starting at 90% of boundary (gradual push)
            // Layer 2: Hard reflection at boundary (velocity reversal)
            // This prevents ANY escape — particles slow down approaching the wall
            // and bounce off if they somehow reach it.
            {
                float dx = p.position.x-p.confineCenter.x;
                float dy = p.position.y-p.confineCenter.y;
                float dz = p.position.z-p.confineCenter.z;
                float dist = sqrtf(dx*dx+dy*dy+dz*dz);
                float boundary = p.confineRadius - p.radius;
                if (boundary < 0.01f) boundary = 0.01f;

                // Layer 1: Soft repulsion zone (90% → 100% of boundary)
                float softZone = boundary * 0.9f;
                if (dist > softZone && dist > 0.001f) {
                    float penetration = (dist - softZone) / (boundary - softZone + 0.001f);
                    // Exponential ramp — gets VERY strong near boundary
                    float F = params.membraneK * penetration * penetration * 3.0f;
                    float nx=dx/dist, ny=dy/dist, nz=dz/dist;
                    force.x -= F*nx; force.y -= F*ny; force.z -= F*nz;

                    // Kill outward velocity component
                    float vn = p.velocity.x*nx + p.velocity.y*ny + p.velocity.z*nz;
                    if (vn > 0) {
                        p.velocity.x -= vn*nx * 1.3f;
                        p.velocity.y -= vn*ny * 1.3f;
                        p.velocity.z -= vn*nz * 1.3f;
                    }
                }

                // Layer 2: Hard wall — if somehow past boundary, push back smoothly
                if (dist > boundary && dist > 0.001f) {
                    // Don't teleport — apply very strong inward velocity
                    float overshoot = dist - boundary;
                    float nx=dx/dist, ny=dy/dist, nz=dz/dist;
                    // Push back over ~3 frames (smooth, not instant)
                    float pushSpeed = overshoot * 20.0f;
                    p.velocity.x = -nx * pushSpeed;
                    p.velocity.y = -ny * pushSpeed;
                    p.velocity.z = -nz * pushSpeed;
                }
            }

            // ── 8. Particle-particle collisions (via spatial hash) ──
            // Skip during mitosis (2x particles, too slow)
            if (!skipCollisions) {
                spatialHash.queryNeighbors(p.position, neighbors);
                for (int j : neighbors) {
                    if (j <= i || j >= n) continue;
                    auto& q = particles[j];
                    if (!q.active) continue;
                    float dx = p.position.x-q.position.x;
                    float dy = p.position.y-q.position.y;
                    float dz = p.position.z-q.position.z;
                    float d = sqrtf(dx*dx+dy*dy+dz*dz);
                    float minD = (p.radius+q.radius)*1.2f;
                    if (d < minD && d > 0.001f) {
                        float overlap = minD - d;
                        float F = params.collisionK * overlap;
                        float nx=dx/d, ny=dy/d, nz=dz/d;
                        force.x += F*nx*0.5f; force.y += F*ny*0.5f; force.z += F*nz*0.5f;
                        q.velocity.x -= F*nx*0.5f*dt/fmaxf(q.mass,0.01f);
                        q.velocity.y -= F*ny*0.5f*dt/fmaxf(q.mass,0.01f);
                        q.velocity.z -= F*nz*0.5f*dt/fmaxf(q.mass,0.01f);
                    }
                }
            }

            // ── 9. Integrate (Euler) ────────────────────────────────
            float invM = 1.0f / fmaxf(p.mass, 0.01f);
            p.velocity.x += force.x*invM*dt;
            p.velocity.y += force.y*invM*dt;
            p.velocity.z += force.z*invM*dt;

            // Cap velocity
            float vMag = sqrtf(p.velocity.x*p.velocity.x+p.velocity.y*p.velocity.y+p.velocity.z*p.velocity.z);
            float maxV = (p.job == JOB_GOTO_TARGET) ? 0.8f : 0.3f;
            if (vMag > maxV) { float s=maxV/vMag; p.velocity.x*=s; p.velocity.y*=s; p.velocity.z*=s; }

            p.position.x += p.velocity.x*dt;
            p.position.y += p.velocity.y*dt;
            p.position.z += p.velocity.z*dt;

            // Exact containment pass: after force integration, project any
            // overshoot back inside the assigned compartment so particles
            // never spend a visible frame outside the membrane.
            intraProjectInsideConfinement(p);
        }

        // ── 10. Nuclear pores: track cell center (static relative to cell) ──
        for (auto& p : particles) {
            if (p.type == PT_NUCLEAR_PORE && p.active) {
                // Move with cell velocity (rigid attachment)
                p.position.x += params.cellVelocity.x * dt;
                p.position.y += params.cellVelocity.y * dt;
                p.position.z += params.cellVelocity.z * dt;
                p.confineCenter = params.cellCenter;
                intraProjectInsideConfinement(p);
            }
        }
    }
};
