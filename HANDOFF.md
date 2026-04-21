# CellSim FEP Handoff — Apple-silicon Mac runner

You are running CellSim's Milestone A accuracy gate on behalf of
a friend. This is one command, ~4–6 hours wall time on an M5 Max,
and produces one tarball to send back.

This README is written for you **and** your Claude Code. Either
can follow it; Claude will probably run it end-to-end with minimal
human intervention.

## What you're running

Absolute hydration free energy (ΔG_hyd) of 12 small molecules
using alchemical free-energy perturbation (FEP). The 12 molecules
span +2.0 kcal/mol (methane, hydrophobic) to −9.7 kcal/mol
(acetamide, hydrophilic).

**Target**: predicted ΔG_hyd within 1.5 kcal/mol MAE of the
published FreeSolv experimental values.

The pipeline is:
  `SMILES → OpenFF Sage 2.1.0 + AM1-BCC charges → OpenMM TIP3P
   solvated System → openmmtools alchemical decoupling →
   MCMC sampling with GHMC integrator + MBAR free-energy estimate`

## Hardware requirements

- **Apple silicon Mac** (M1/M2/M3/M4/M5, any Pro/Max/Ultra
  variant). Intel Mac also works but slower.
- **32 GB RAM recommended**, 16 GB will work but may swap.
- **~5 GB free disk** in the CellSim repo (data + output).
- **Plugged in**, display-sleep OK but full-sleep NOT OK.
  Set before running:
  ```bash
  sudo pmset -c disablesleep 1
  ```
  Re-enable after:
  ```bash
  sudo pmset -c disablesleep 0
  ```

## One-time setup (first run only)

1. Install mamba/miniforge if you don't have it:
   ```bash
   brew install miniforge
   ```

2. Clone the release:
   ```bash
   git clone https://github.com/Rockman6/CellSim-3D-Multicellular-Biology-Drug-Testing-Simulator.git CellSim
   cd CellSim
   git checkout <release-tag>   # Henry will tell you which tag
   ```

3. Create the conda environment:
   ```bash
   mamba env create -f environment.yml
   conda activate cellsim
   ```
   This installs ~60 packages and takes 10–20 minutes.

4. Install `terminal-notifier` so the script can ping you when
   it finishes:
   ```bash
   brew install terminal-notifier
   ```

5. Verify the install:
   ```bash
   ./scripts/cellsim doctor
   ```
   You should see "**39/39 checks passed**". If fewer pass, send
   Henry the output — do **not** proceed.

## The run

One command:

```bash
./scripts/run_freesolv_m5max.sh
```

The script will:
1. Activate the `cellsim` conda env.
2. Print the OpenMM platform list (Metal / OpenCL / CPU) so you
   can see which backend is used.
3. Run `cellsim doctor` and abort if it fails.
4. Execute the 12-compound FreeSolv gate at production parameters:
   - 11 alchemical windows per leg
   - 50 ps equilibration + 50 ps production per window
   - 100 samples per window (2 ps stride)
   - 2 legs per compound (vacuum + TIP3P-solvated)
   - = ~2.2 ns simulated MD per compound × 12 = 26 ns total
5. Write everything to `run/fep/<timestamp>/` and bundle into
   `freesolv_m5max_<timestamp>.tar.gz`.
6. **Auto-run `cellsim fep-report`** on the result — writes
   `report.md` + `parity.png` into the output dir and prints the
   PASS/FAIL verdict on your terminal *before* you send. So you
   see immediately whether the gate cleared; Henry gets the
   verdict pre-computed inside the tarball.
7. Fire a macOS notification ("CellSim FreeSolv FEP: PASS/FAIL")
   when complete.

Expected wall time:
- **M5 Max (40-core GPU)**: 4–6 hours
- **M2/M3 Max**: 6–9 hours
- **M1 Pro/Max**: 8–12 hours
- **Intel Mac (CPU only)**: 24–48 hours — not recommended

You can close the Terminal window during the run — the process
continues. Use `tail -f run/fep/*/run.log` in another Terminal to
watch progress.

## What to send back

When the script finishes, send **one file**:

```
freesolv_m5max_<timestamp>.tar.gz
```

(in the CellSim repo root). Send via:
- **WeChat file transfer** (works for files up to ~100 MB; tarball
  will be ~1–5 MB)
- **AirDrop** (if you're close geographically)
- **Email / Google Drive / Dropbox**
- **GitHub Gist** (if public data is OK — it is; FreeSolv
  numbers are public benchmark data)

The tarball contains:

| File | What it is | Why Henry wants it |
|---|---|---|
| `env.log` | OpenMM/Python/platform versions | Reproducibility — must match exactly |
| `doctor.log` | cellsim doctor output | Catches env drift between machines |
| `run.log` | Full pipeline stdout | Contains per-compound ΔG, per-window GHMC acceptance, wall times, any NaN warnings |
| `freesolv_results.csv` | 12-row CSV of predictions | The numbers Henry analyses |
| `report.md` | Auto-generated PASS/FAIL verdict + per-compound table | Henry can read it in a text editor; doesn't need to re-run the analyser |
| `parity.png` | Predicted vs experimental scatter w/ ±1.5 kcal/mol gate band | Visual sanity for the professor debrief |
| `table.csv` | Normalised CSV (adds abs_residual, sign_correct, flags) | Consumable by the Campaign-2 prior-emitter |

## What Henry specifically wants in the CSV

Each row of `freesolv_results.csv` has:

```
name, smiles, dG_expt_kcalmol, dG_pred_kcalmol, uncertainty_kcalmol,
residual_kcalmol, wall_seconds, ok, reason
```

Specifically he'll compute:
- **MAE** = mean |residual| across all 12 → gate is ≤ 1.5 kcal/mol.
- **Per-compound sign correctness** — methane should give ≈ +2 kcal/mol
  (positive), not −2 (negative); the CPU pilot gave −2 due to
  under-sampling, so seeing it flip to positive on M5 Max is the
  load-bearing finding.
- **GHMC acceptance rate** — must be ≥ 70% per window (≥ 75%
  preferred). Low acceptance means the ΔG is unreliable.
- **Wall time per compound** — tells Henry the M5 Max throughput
  so he can plan the next round (EGFR kinase relative FEP).

## Troubleshooting

### `cellsim doctor` fails

Copy the terminal output, send it to Henry. Do NOT try to fix it.

### Script errors immediately with "python: command not found"

You didn't activate the env. Run:
```bash
conda activate cellsim
```

### Script errors with "Particle coordinate is NaN"

**Send Henry the log anyway.** That failure mode is itself
scientifically valuable. The CPU version sometimes hits this at
shorter MD lengths; Henry needs to know if it shows up on
M5 Max too or if longer sampling dissolves it.

### Run times out after 6 hours and hasn't finished

That's fine. The tarball won't exist (script never reached the
`tar czf` step), so you'll need to build one yourself from
whatever's in `run/fep/<stamp>/`:

```bash
# Auto-compute the PASS/FAIL verdict on the partial data + bundle
STAMP=$(ls -1 run/fep/ | tail -1)
./scripts/cellsim fep-report "run/fep/${STAMP}" \
    --yaml benchmarks/fep/freesolv_12.yaml \
    --out-dir "run/fep/${STAMP}"
tar czf "freesolv_m5max_${STAMP}_partial.tar.gz" \
    -C run/fep "${STAMP}/"
```

The analyser will automatically flag it as "FAIL (partial run)"
in the `report.md` header — no manual annotation needed. Send
`freesolv_m5max_<stamp>_partial.tar.gz`.

### Machine overheats / fan roars / battery drains

Normal for sustained GPU work. Machine is plugged in and
physically safe. Continue.

### You want to reduce the run to ~2 hours

Edit `scripts/run_freesolv_m5max.sh`, change:
```
--production-steps 25000
```
to:
```
--production-steps 10000
```

This gives a less-converged result (still useful as a first pass).
Tell Henry you used shorter production when you send the tarball.

## What this actually tests (Henry's framing)

The CellSim project got a long critique from Henry's professor
about whether the simulator can actually replace wet-lab
experiments — specifically whether "physics-derived priors" from
FEP are good enough to build cell biology on top of. Methane
hydration is the simplest FEP test: if we can't get methane right,
nothing else matters.

Your run is **the M5 Max-accessible test** of whether longer MD
sampling dissolves the entropy-capture issue Henry saw on CPU.
Professor threshold:
- Acceptance ≥ 75%
- ΔG_hyd(methane) = +2.0 ± 0.5 kcal/mol
- Uncertainty < 0.6 kcal/mol

Hit those and Milestone A clears, Campaign 2 (real cell biology)
opens.

## Questions

Ask Henry directly. He will give Claude explicit instructions
too — if your Claude is unsure about something, let it paste
the question into WeChat to Henry and wait.

Thanks for lending the machine. Good luck.
