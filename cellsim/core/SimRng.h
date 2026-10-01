#pragma once
// ══════════════════════════════════════════════════════════════════════════
//  SimRng.h — single, deterministic, seedable RNG entry point.
// ──────────────────────────────────────────────────────────────────────────
//  Phase P1 (2026-04-20): empirical audit showed the current simulator is
//  bit-deterministic in CSV content modulo the wall_sec column. That
//  determinism is accidental — it relies on libc's C99 §7.20.2.2 rule that
//  `rand()` starts from seed 1 when `srand()` is never called. Any future
//  code that calls `std::random_device`, `chrono`-seeded RNGs, or threads
//  would shatter it silently.
//
//  This header is the single, named entry point for every future RNG draw.
//  Style rule for new code (enforced by review, not yet by CI):
//      Do NOT call rand() or RAND_MAX directly in simulation code.
//      Call simrng::uniform() / simrng::gauss() / simrng::uniformInt().
//
//  Implementation today: wraps libc rand()/srand() so the 72 pre-existing
//  rand() call sites keep their exact behaviour — zero biology drift, zero
//  calibration shift, 5.1 % CTC-HeLa metric preserved. Those 72 sites will
//  migrate to simrng::* in a later PR; once migrated, the implementation
//  can swap to std::mt19937 in a single line without touching call sites.
//
//  Seeding: callers (main.mm, tests/headless_*.cpp) should invoke
//      simrng::seed(user_seed_or_default);
//  once at startup. Default seed `DEFAULT_SEED` = 1 matches libc's implicit
//  initial state so unseeded runs reproduce today's bit streams exactly.
//
//  Thread safety: libc rand() is NOT thread-safe, and simulation hot paths
//  are single-threaded (headless_doubling, headless_infection). The
//  interactive CellSim binary may call into simrng from a render thread for
//  cosmetic particle placement; that is intentionally non-claimed
//  behaviour and does not participate in any bit-identical gate.
// ══════════════════════════════════════════════════════════════════════════

#include <cstdlib>
#include <cmath>
#include <cstdint>

namespace simrng {

// Matches libc's implicit initial `srand(1)` per C99 §7.20.2.2.
// Headless binaries that want reproducibility without CLI args inherit this.
inline constexpr uint32_t DEFAULT_SEED = 1u;

// Set the global seed. Call once at process startup; subsequent calls are
// allowed (tests may reseed between scenarios). Accepts a uint32 so callers
// can pass hex constants like 0xCE11517E cleanly.
inline void seed(uint32_t s) { std::srand(static_cast<unsigned int>(s)); }

// Uniform float in [0, 1]. Intentionally inclusive at both ends, matching
// the historical `(float)rand()/(float)RAND_MAX` behaviour.
inline float uniform() {
    return static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
}

// Uniform float in [lo, hi].
inline float uniform(float lo, float hi) {
    return lo + (hi - lo) * uniform();
}

// Uniform int in [0, n) using C's rand() modulo. N must be small relative
// to RAND_MAX (on macOS RAND_MAX = 0x7FFFFFFF so this is fine for any
// realistic simulation index).
inline int uniformInt(int n) {
    return (n > 0) ? (std::rand() % n) : 0;
}

// Box-Muller standard-normal draw. Matches the ipRandGauss() pattern that
// lives in IntracellularPhysics.h today.
inline float gauss() {
    float u1 = uniform() + 1e-8f;
    float u2 = uniform();
    return std::sqrt(-2.0f * std::log(u1)) *
           std::cos(2.0f * static_cast<float>(M_PI) * u2);
}

} // namespace simrng
