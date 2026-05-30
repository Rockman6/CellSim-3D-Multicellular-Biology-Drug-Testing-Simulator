# Milestone B run — binding ΔG (streptavidin + EGFR)

For the friend running on the M5 Max. Same workflow shape as
Milestone A (`run_freesolv_m5max.sh`), now for binding instead
of hydration. Code is on tag `milestone-a-pilot-3` and forward
(any commit on `main` since `8ab37a1` is binding-correct).

## What changed since Milestone A

- Sign convention on binding was always correct (verified locally
  by a methane-in-weak-pocket cycle check). Hydration's c461053
  sign bug was specific to the hydration cycle and never touched
  `compute_absolute_binding_dg`.
- `--resume` + atomic CSV writes are in both binding bench paths
  (streptavidin + EGFR). Same crash recovery as freesolv.
- `OUT_DIR` env override now consistent across all three run
  scripts. Resume a crashed run with:

    ```
    OUT_DIR=run/fep/streptavidin_20260523_103045 \
        bash scripts/run_binding_streptavidin_gpu.sh
    ```

## What to run

1. **Streptavidin (smaller, faster — start here):**

   ```
   git checkout milestone-a-pilot-3  # or `main` for latest
   conda activate cellsim
   bash scripts/run_binding_streptavidin_gpu.sh
   ```

   4 compounds (biotin, desthiobiotin, 2-iminobiotin,
   biotin_methyl_ester). Validate already gives an OpenCL ETA of
   ~15 min on M-series. On the friend's pilot-1 M5 Max where
   OpenCL was 4-6× wall-time of the estimate, real wall is
   probably 1-2h. If the pipeline is correct, biotin should land
   close to its published −18.30 kcal/mol.

2. **EGFR (the headline test):**

   ```
   bash scripts/run_binding_egfr_gpu.sh
   ```

   6 compounds. Validate flags lapatinib (11 rotatable bonds) and
   suggests `--production-steps 100000` if you want it to converge
   cleanly — the script defaults are conservative; bump there if
   the report flags low GHMC on lapatinib.

   The HEADLINE for this set is Kendall τ ≥ 0.6 — Vina gave
   τ ≈ −0.49 on the same series (anti-correlated! actively
   worse than guessing), so a positive τ is the rescue claim. If
   the FEP run gets τ ≥ 0.6, we have a directly publishable
   "FEP rescues docking" result.

## When the run finishes

Same flow as Milestone A:

```
# Local verdict (auto-runs inside the script; check the markdown)
cat run/fep/streptavidin_*/report.md

# One-line Slack/standup status
python scripts/csv_tldr.py run/fep/streptavidin_*/streptavidin_results.csv

# Prof email draft
python scripts/fill_prof_email.py run/fep/streptavidin_* \
    --hardware "Apple M5 Max (40-core GPU)" \
    --platform "Metal / OpenCL" \
    --out prof_streptavidin.txt

# Tarball for Henry
tar -czf streptavidin_$(date +%Y%m%d).tar.gz run/fep/streptavidin_*
```

## Sanity expectations (so you know what "right" looks like)

| compound | expt ΔG_bind | rough physical sign |
|---|---:|---|
| biotin / streptavidin | −18.30 | very strong binder (negative) |
| desthiobiotin / strep | −13.20 | strong binder (negative) |
| 2-iminobiotin / strep | −10.80 | strong binder (negative) |
| biotin_methyl_ester / strep | −14.10 | strong binder (negative) |
| erlotinib / EGFR | −11.86 | strong binder (negative) |
| gefitinib / EGFR | −10.20 | strong binder (negative) |
| AG-1478 / EGFR | −11.62 | strong binder (negative) |
| lapatinib / EGFR | −10.88 | strong binder (negative) |
| 4-anilinoquinazoline / EGFR | −8.17 | medium binder (negative) |
| tyrphostin_AG-494 / EGFR | −8.11 | medium binder (negative) |

If ANY compound predicts ≥ 0 (positive ΔG_bind), the pipeline says
"non-binder" for a known binder — that's the same gate as the
hydration sign check, and it'd be a red flag worth reporting back
before sending the prof email.

## If the friend hits a problem

- `cellsim doctor` failed → check `OUT_DIR/doctor.log`, the
  message names the specific install issue
- MBAR overlap failure → the biologist-actionable translator
  fires with concrete `--n-windows N --n-production-steps M`
  bumps. Often means tight binder needs deeper sampling.
- Lid-close / power-loss → just re-run with `OUT_DIR=<that-dir>
  bash <same-script>` and `--resume` picks up where it died.

## Why Milestone B matters

Milestone A (hydration FEP) showed the alchemical pipeline is
end-to-end correct on small-molecule solvation. Milestone B
extends to **protein-ligand binding**, which is the actual drug-
discovery quantity. The DDM (double-decoupling method) used here
is the standard reference protocol (Hamelberg-Gilson 2004), and
the openmmtools alchemical primitives that make it work are
exactly the same ones the validated Milestone A pipeline uses.
If Milestone A passes (it does, MAE 1.42), Milestone B is the
high-value proof point that the chemistry layer can be trusted
for binding screens.
