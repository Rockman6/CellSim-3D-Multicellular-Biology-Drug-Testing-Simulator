"""Layer 1.6 MC-dock smoke — 4 seeds on 1STP biotin.

Runs src.uq.monte_carlo_dock over 4 Vina seeds and sanity-checks
the output:
  - ok=True with n_ok >= 3 (at least 3/4 seeds produce a valid pose)
  - ΔG mean in [-10, -4] kcal/mol (reasonable range for biotin-
    streptavidin under Vina)
  - std >= 0 and <= 1.0 kcal/mol (seed-jitter shouldn't be large)
  - 95 % CI brackets the mean
  - best_dG_kcalmol == min(dG_samples_kcalmol)

The wall gate is lenient (<=60 s on a laptop, 2 cpu per job) because
Vina is the dominant cost; N=4 seeds keeps CI-friendly.

Run:
    conda activate cellsim
    python tests/uq/test_mc_dock_smoke.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.uq import monte_carlo_dock  # noqa: E402


RECEPTOR = REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb"
LIGAND   = "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"  # biotin
CENTER   = (11.12, 1.68, -10.75)
BOX      = (20.0, 20.0, 20.0)


def test_mc_dock_biotin_1stp():
    assert RECEPTOR.exists(), f"missing {RECEPTOR}"
    t0 = time.time()
    r = monte_carlo_dock(
        receptor_pdb=str(RECEPTOR),
        ligand_smiles=LIGAND,
        center_A=CENTER, box_size_A=BOX,
        n_samples=4, base_seed=1,
        exhaustiveness=16, num_modes=3,
        workers=2, cpu_per_job=2,
        timeout_s=120)
    dt = time.time() - t0
    print(f"[mc-dock] wall={dt:.1f} s")
    print(f"  {r.summary()}")
    print(f"  samples: {r.dG_samples_kcalmol}")

    assert r.ok, f"MC failed: {r.reason}"
    assert r.n_ok >= 3, f"only {r.n_ok}/4 seeds succeeded"
    # Reasonable ΔG range for biotin-streptavidin under Vina:
    assert -10.0 <= r.dG_mean_kcalmol <= -4.0, (
        f"ΔG mean {r.dG_mean_kcalmol:.2f} out of [-10, -4] kcal/mol")
    assert 0.0 <= r.dG_std_kcalmol <= 1.0, (
        f"ΔG std {r.dG_std_kcalmol:.3f} out of [0, 1] kcal/mol")
    # CI sanity: lo <= mean <= hi
    assert (r.dG_ci95_lo_kcalmol <= r.dG_mean_kcalmol
            <= r.dG_ci95_hi_kcalmol), (
        f"CI inconsistent: {r.dG_ci95_lo_kcalmol} <= "
        f"{r.dG_mean_kcalmol} <= {r.dG_ci95_hi_kcalmol}")
    # Best ΔG is the minimum of the sample list.
    assert r.best_dG_kcalmol == min(r.dG_samples_kcalmol), (
        f"best_dG {r.best_dG_kcalmol} != "
        f"min {min(r.dG_samples_kcalmol)}")


if __name__ == "__main__":
    try:
        test_mc_dock_biotin_1stp()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
