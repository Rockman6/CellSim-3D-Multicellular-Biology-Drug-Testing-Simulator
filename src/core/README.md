# src/core — Physics-neutral utilities

Reused from the OLD/ tree because they are not biology-coupled:

- `SimRng.h` — seedable deterministic RNG wrapper. Default seed 1
  matches C99 implicit `srand(1)`. Basis for reproducible ensemble
  MD, conformal splits, and blind-bench seed pinning.
- `TelemetryLog.h` — deterministic CSV logger; `wall_sec` emitted
  as `0.0` so cross-machine output is bit-identical.
- `Constants.h` — physical constants (R, T, k_B, N_A, Å, ps, …).

No changes vs OLD/; identical headers are kept here so the new tree
builds without any dependency on `OLD/src/`.
