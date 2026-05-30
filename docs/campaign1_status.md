# Campaign 1 — Status (2026-05-30)

One-page summary of where each Campaign-1 layer stands and what
the realistic timeline is to Stage 1 completion on CPU-only.

## Per-layer status

| Layer | Status | Evidence | Open work |
|---|---|---|---|
| **1.1 Chem foundation** (RDKit + OpenFF + AM1-BCC) | ✅ DONE | 10k ChEMBL round-trip, ADMET + AM1-BCC cache, profile dashboard | — |
| **1.2 Classical MD** (OpenMM + ff14SB + TIP3P) | ✅ DONE (MVP) | 1 ps ubiquitin Cα RMSD 0.74 Å; 100 ns ubiquitin gate awaits GPU | full 100 ns gate (GPU) |
| **1.3a Docking** (Vina + Meeko + PoseBusters + fpocket) | ✅ DONE | mini-bench 2/3 canonical re-dock; off-target panel; CYP3A4 SoM; strain-gate; triage rules + viewer | — |
| **1.3b FEP Milestone A** (hydration ΔG) | ✅ **PASS** (today) | FreeSolv-12 MAE **1.42 kcal/mol** (gate ≤ 1.5), Pearson r +0.913, GHMC 99%, tag `milestone-a-pilot-3` | 12/12 closure on 2 polar fails (structural, doc'd; not blocking) |
| **1.3c FEP Milestone B** (binding ΔG) | ⏳ in flight | streptavidin bench RUNNING on CPU (PID 75132, started 2026-05-30 22:15, ETA ~20 h); EGFR queued | both benches → fep-report → tarball |
| **1.4 Quantum** (xTB GFN2 + PySCF DFT) | ✅ DONE (MVP) | 10/10 sanity + 3/3 CYP3A4 SoM; 2/3 on literature validation (aspirin ✓, midazolam ✓, diazepam doc'd-fail) | optional: extend literature set to 20 |
| **1.5 Coarse-grained Martini 3** | ⏳ scaffold only | `src/cg/{bilayer,protein_cg}.py` are NotImplementedError stubs (~90 LOC total) | full bilayer builder + martinize2 wrapper + OpenMM-Martini MD driver + area-per-lipid validation |
| **1.6 UQ** (Sobol + MC + MAPIE) | ✅ DONE | MC-dock seeds, Sobol sensitivity, conformal quantiles, streptavidin/trypsin/EGFR calibration bundles | — |
| **1.7 Blind harness** (PDBBind / CASF / PoseBusters / ChEMBL) | ⏳ partial | PoseBusters integrated, fpocket integrated, 3-cocrystal mini-bench shipped; `benchmarks/{pdbbind,casf}/` are empty dirs | PDBBind 500-compound CPU subset (full 5k needs GPU); CASF-2016 ranking gate |

## CPU-only timeline to Stage 1 complete

**This week (5 days CPU dedicated):**
- Streptavidin binding bench: ~20 h
- EGFR binding bench: ~80 h
- Milestone B verdict + tarball + prof email: ~Friday 2026-06-05

**Layer 1.5 Martini 3 (4-6 weeks focused dev, no GPU needed):**
- bilayer builder (martinize2 wrapper + cube/tetrahedral packing)
- protein → CG converter (martinize2 + elastic-network)
- OpenMM-Martini MD driver (1-10 µs runs on CPU are feasible at CG resolution)
- area-per-lipid validation against POPC literature (40-65 Å²/lipid)
- viewer: lipid-head colour-coded bilayer with area gauge
- Target: **mid-July 2026**

**Layer 1.7 Blind harness (4-6 weeks, CPU-bound subset version):**
- PDBBind refined-set fetcher + 500-compound subset processor
- CASF-2016 scoring + ranking + screening + docking-power gates
- ChEMBL kinase panels (already partial — has the EGFR rank-order)
- automated red-team slot scaffold (quarterly external adversarial)
- Target: **end of August 2026**

**Stage 1 (Campaign 1) complete: end August / early September 2026.**

If GPU access returns mid-stream, Layer 1.7 PDBBind full 5k subset
becomes feasible and pulls the date forward by 2-4 weeks.

## What's currently running

- PID 75132: `bash scripts/run_binding_streptavidin_gpu.sh` (despite the GPU in the name, runs on CPU because that's what's available; auto-selects fastest available backend)
- Output: `run/fep/streptavidin_20260530_221529/`
- Expected completion: 2026-05-31 ~18:00 (Sun afternoon)

## Tarballs to expect

1. `run/fep/streptavidin_20260530_221529.tar.gz` — Sun afternoon
2. `run/fep/egfr_XXXX.tar.gz` — Fri 2026-06-05

Both get processed via `bash scripts/finalize_run.sh <run-dir>`
(generates report.md + table.csv + parity.png + csv_tldr + prof
email draft + tarball, one command).

## Open questions

- **Are we shipping Milestone B even if EGFR Kendall τ < 0.6?**
  The interesting result either way. PASS = headline "FEP rescues
  Vina's anti-correlated docking". FAIL with diagnosis = honest
  intermediate result.
- **Layer 1.5 priority vs Layer 1.7:** can be done in parallel
  (no compute conflict) but bandwidth-limited. My vote: start 1.7
  first because the data plumbing is mostly orthogonal to Milestone
  B running, then 1.5 once binding completes and frees the mental
  bandwidth for the heavier CG implementation.
