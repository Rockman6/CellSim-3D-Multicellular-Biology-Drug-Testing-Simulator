"""Layer 1.1 smoke test: parametrise 10 canonical drugs.

Exit-criterion surrogate for Layer 1.1 (the full criterion is 10 k
ChEMBL compounds with ≥ 99 % success). This 10-compound file keeps
the gate runnable in seconds on any dev box.

Run:
    pytest tests/chem/test_parametrize_smoke.py -v
    # or without pytest:
    python tests/chem/test_parametrize_smoke.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Let pytest and `python …` both find the repo root.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.chem import parametrize_smiles  # noqa: E402

SMI_FILE = REPO_ROOT / "benchmarks" / "chembl" / "smoke_10.smi"


def load_smiles_file(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(None, 1)
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else smi
        entries.append((smi, name))
    return entries


def _rdkit_ok(r) -> bool:
    """RDKit-only success: embed produced 3D coords + elements + bonds."""
    return (r.positions_nm is not None
            and r.elements is not None
            and r.bonds is not None
            and r.n_atoms is not None and r.n_atoms > 0)


def test_smoke_pass_rate():
    """Two-tier gate:

    - Full (OpenFF present): ≥ 8/10 parameterise to an OpenMM system.
    - Degraded (OpenFF absent): ≥ 9/10 RDKit-embed + 3D coords. This
      keeps CI meaningful on laptops without conda while the full
      gate is reserved for the conda-configured runner.
    """
    entries = load_smiles_file(SMI_FILE)
    results = [parametrize_smiles(smi) for smi, _ in entries]
    full_ok = sum(1 for r in results if r.ok)
    rdkit_ok = sum(1 for r in results if _rdkit_ok(r))

    report_lines = []
    for (smi, name), r in zip(entries, results):
        if r.ok:
            tag = "FULL"
            detail = f"formula={r.formula} atoms={r.n_atoms}"
        elif _rdkit_ok(r):
            tag = "RDK "
            detail = (f"formula={r.formula} atoms={r.n_atoms}  "
                      f"(openff unavailable)")
        else:
            tag = "FAIL"
            detail = f"reason: {r.reason}"
        report_lines.append(f"  {tag}  {name:25s}  {detail}")
    report = "\n".join(report_lines)
    print("\nLayer 1.1 smoke:\n" + report
          + f"\n  full  pass-rate: {full_ok}/{len(entries)}"
          + f"\n  rdkit pass-rate: {rdkit_ok}/{len(entries)}")

    # Pick the tier by whether OpenFF is actually installed.
    try:
        import openff.toolkit  # noqa: F401
        have_openff = True
    except ImportError:
        have_openff = False

    if have_openff:
        assert full_ok >= 8, (
            f"full pass-rate {full_ok}/{len(entries)} < 8; details:\n"
            + report)
    else:
        assert rdkit_ok >= 9, (
            f"rdkit-only pass-rate {rdkit_ok}/{len(entries)} < 9; "
            f"details:\n" + report)


if __name__ == "__main__":
    try:
        test_smoke_pass_rate()
        print("PASS")
        sys.exit(0)
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
