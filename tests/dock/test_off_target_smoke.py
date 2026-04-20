"""Off-target smoke — biotin vs 3 receptors.

Gates (small cache-aware run):
  - n_receptors = 3, all entries present in result
  - at least 2 of 3 produce a valid ΔG (the third may fail if
    fpocket can't find a druggable pocket)
  - streptavidin (the biological target for biotin) lands at rank 1
    with ΔG < -6.5 kcal/mol
  - selectivity ΔΔG (best vs 2nd-best) > 0 (strongest target beats
    off-targets)
  - CSV writer produces a parseable file with the expected columns

Run:
    conda activate cellsim
    python tests/dock/test_off_target_smoke.py
"""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cache import Cache  # noqa: E402
from src.dock import off_target_screen  # noqa: E402
from src.dock.off_target import write_csv  # noqa: E402


LIGAND = "OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12"   # biotin

RECEPTORS = [
    ("streptavidin", str(REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb")),
    ("trypsin",      str(REPO_ROOT / "benchmarks" / "dock" / "3ptb.pdb")),
    ("ubiquitin",    str(REPO_ROOT / "benchmarks" / "md" / "1ubq.pdb")),
]


def test_off_target_biotin():
    for _, p in RECEPTORS:
        assert Path(p).exists(), f"missing {p}"

    with tempfile.TemporaryDirectory(prefix="cellsim-off-") as tmp:
        cache = Cache(Path(tmp) / "c.sqlite")
        r = off_target_screen(
            LIGAND, RECEPTORS,
            exhaustiveness=16, num_modes=3, seed=1, cpu_per_job=2,
            cache=cache)
        print(r.summary())

        assert len(r.entries) == 3
        ok = [e for e in r.entries if e.ok]
        assert len(ok) >= 2, f"only {len(ok)}/3 docked"

        ranked = r.sorted_by_affinity()
        top = ranked[0]
        assert top.ok, "top entry not ok"
        assert top.name == "streptavidin", (
            f"expected streptavidin at rank 1, got {top.name}")
        assert top.dG_kcalmol < -6.5, (
            f"streptavidin ΔG {top.dG_kcalmol} weaker than -6.5")

        sel = r.selectivity_kcalmol()
        assert sel is not None and sel > 0, (
            f"selectivity ΔΔG = {sel}; expected > 0")
        print(f"[off] selectivity ΔΔG = {sel:+.2f} kcal/mol")

        # CSV writer.
        csv_path = Path(tmp) / "out.csv"
        write_csv(r, csv_path)
        assert csv_path.exists()
        # utf-8-sig so the BOM cellsim writes is stripped before
        # DictReader sees the header.
        with csv_path.open(encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        # 3 rows (one per receptor), with the expected header keys.
        assert len(rows) == 3
        for k in ["rank", "name", "receptor_pdb", "ok",
                  "dG_kcalmol", "Kd_nM"]:
            assert k in rows[0], f"CSV missing column '{k}'"
        cache.close()


if __name__ == "__main__":
    try:
        test_off_target_biotin()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
