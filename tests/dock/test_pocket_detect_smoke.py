"""Layer 1.3 pocket-detection smoke — fpocket on 1STP.

Verifies that fpocket correctly identifies streptavidin's biotin
binding pocket as the top druggable pocket:
  - Top pocket drug score >= 0.5 (well-druggable cutoff)
  - Top pocket center within 3 Å of the known biotin centroid
    (11.12, 1.68, -10.75) computed from the crystal BTN HETATM
    block and documented in benchmarks/dock/README.md
  - Suggested search box is at least 18 Å per axis and at most 30 Å

Run:
    conda activate cellsim
    python tests/dock/test_pocket_detect_smoke.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.dock import detect_pockets  # noqa: E402

RECEPTOR = REPO_ROOT / "benchmarks" / "dock" / "1stp.pdb"

# Known biotin centroid (from extract_hetatm_ligand + mean over xyz).
KNOWN_BIOTIN_CENTER = (11.12, 1.68, -10.75)
CENTER_TOLERANCE_A = 3.0


def test_fpocket_1stp():
    assert RECEPTOR.exists(), f"missing {RECEPTOR}"
    pockets = detect_pockets(RECEPTOR, top_k=5)
    assert pockets, "fpocket produced no pockets (is fpocket in PATH?)"
    top = pockets[0]
    print(f"[pocket] detected {len(pockets)} pockets; top:")
    for p in pockets[:3]:
        print(f"  {p.summary()}")

    assert top.drug_score is not None and top.drug_score >= 0.5, (
        f"top pocket drug score {top.drug_score} < 0.5")

    assert top.center_A is not None, "no center_A on top pocket"
    dx = top.center_A[0] - KNOWN_BIOTIN_CENTER[0]
    dy = top.center_A[1] - KNOWN_BIOTIN_CENTER[1]
    dz = top.center_A[2] - KNOWN_BIOTIN_CENTER[2]
    d = (dx * dx + dy * dy + dz * dz) ** 0.5
    print(f"[pocket] top center distance to known biotin site: "
          f"{d:.2f} Å")
    assert d < CENTER_TOLERANCE_A, (
        f"top pocket center off by {d:.2f} Å > "
        f"{CENTER_TOLERANCE_A} Å tolerance")

    assert top.suggested_box_A is not None, "no suggested_box_A"
    for dim in top.suggested_box_A:
        assert 18.0 <= dim <= 30.0, (
            f"suggested box dim {dim} not in [18, 30]")


if __name__ == "__main__":
    try:
        test_fpocket_1stp()
        print("PASS")
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"ERROR: {e}")
        sys.exit(2)
