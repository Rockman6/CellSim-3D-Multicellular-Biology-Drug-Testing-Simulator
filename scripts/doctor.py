#!/usr/bin/env python3
"""`cellsim doctor` — end-to-end install health check.

Walks the Campaign-1 pipeline in the smallest possible way and
tells the biologist exactly what's broken, so a half-set-up env
doesn't turn into a day of mysterious errors.

Steps (each one < 1 s):
    1. Conda env activated? AMBERHOME set?
    2. rdkit / openff.toolkit / openmm / pdbfixer / vina / meeko /
       posebusters / xtb / salib / fpocket / pyscf importable?
    3. Vina binary on PATH?
    4. xTB binary on PATH?
    5. fpocket binary on PATH?
    6. bundled cocrystal 1STP loadable?
    7. RDKit can parse and embed aspirin?
    8. compound_hash deterministic?
    9. src.cache round-trip works?
"""

from __future__ import annotations

import importlib
import os
import shutil
import sys
from pathlib import Path


GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
DIM = "\033[90m"
RESET = "\033[0m"
CHECK = GREEN + "✓" + RESET
CROSS = RED + "✗" + RESET
WARN = YELLOW + "!" + RESET


def tag_pass(msg: str) -> None:
    print(f"  {CHECK}  {msg}")


def tag_fail(msg: str) -> None:
    print(f"  {CROSS}  {msg}")


def tag_warn(msg: str) -> None:
    print(f"  {WARN}  {msg}")


def tag_info(msg: str) -> None:
    print(f"  {DIM}·  {msg}{RESET}")


def main() -> int:
    ok_count = 0
    total = 0

    REPO = Path(__file__).resolve().parent.parent

    def check(label: str, cond: bool, hint: str = "") -> None:
        nonlocal ok_count, total
        total += 1
        if cond:
            tag_pass(label)
            ok_count += 1
        else:
            tag_fail(label)
            if hint:
                tag_info(hint)

    print()
    print("  cellsim doctor")
    print("  " + "─" * 60)

    # --- 1. env activation ----------------------------------------
    print(f"\n  1. conda env")
    prefix = os.environ.get("CONDA_PREFIX", "")
    check("CONDA_PREFIX set",
          bool(prefix),
          "run 'conda activate cellsim' first")
    check("AMBERHOME set (needed for AM1-BCC charges)",
          bool(os.environ.get("AMBERHOME")),
          "if CONDA_PREFIX is set but AMBERHOME is not, your env "
          "may be under-provisioned — re-run 'conda env create -f "
          "environment.yml'")
    if prefix and "cellsim" in prefix:
        tag_info(f"prefix: {prefix}")

    # --- 2. Python imports ----------------------------------------
    print(f"\n  2. Python packages")
    for name, module in [
        ("rdkit", "rdkit"),
        ("openmm", "openmm"),
        ("openff.toolkit", "openff.toolkit"),
        ("pdbfixer", "pdbfixer"),
        ("vina", "vina"),
        ("meeko", "meeko"),
        ("posebusters", "posebusters"),
        ("xtb-python", "xtb"),
        ("SALib", "SALib"),
        ("pyscf", "pyscf"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
    ]:
        try:
            m = importlib.import_module(module)
            ver = getattr(m, "__version__", "?")
            check(f"{name:<18s} {ver}", True)
        except ImportError as e:
            check(f"{name:<18s} import failed",
                  False, f"pip/conda install {name}")

    # --- 3. CLI binaries ------------------------------------------
    print(f"\n  3. CLI binaries")
    for exe in ["vina", "xtb", "fpocket", "antechamber",
                "mk_prepare_receptor.py"]:
        path = shutil.which(exe)
        if path:
            tag_pass(f"{exe:<30s} {path}")
            ok_count += 1
        else:
            tag_fail(f"{exe:<30s} not on PATH")
        total += 1

    # --- 4. Bundled benchmark files -------------------------------
    print(f"\n  4. Bundled benchmarks")
    for rel in [
        "benchmarks/chembl/smoke_10.smi",
        "benchmarks/md/1ubq.pdb",
        "benchmarks/dock/1stp.pdb",
        "benchmarks/dock/1m17.pdb",
        "benchmarks/dock/3ptb.pdb",
        "benchmarks/dock/mini_bench.yaml",
    ]:
        p = REPO / rel
        check(f"{rel:<40s} ({p.stat().st_size if p.exists() else 0} B)",
              p.exists(),
              f"file missing — check repo layout")

    # --- 5. Core smoke calls --------------------------------------
    print(f"\n  5. Core smoke calls")
    sys.path.insert(0, str(REPO))

    # 5a RDKit parses aspirin
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O"))
        ok = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == 0
        check("RDKit parse + embed aspirin", ok,
              "RDKit installed but couldn't embed — try rebuilding env")
    except Exception as e:
        check(f"RDKit parse + embed aspirin  [{e}]", False)

    # 5b compound_hash deterministic
    try:
        from src.cache import compound_hash
        h1 = compound_hash("CC(=O)O")
        h2 = compound_hash("OC(C)=O")
        check(f"compound_hash stable under reorder (h={h1[:8] if h1 else 'none'})",
              h1 is not None and h1 == h2,
              "hashing produced None or inconsistent values")
    except Exception as e:
        check(f"compound_hash  [{e}]", False)

    # 5c Cache round-trip
    try:
        import tempfile
        from src.cache import Cache
        with tempfile.TemporaryDirectory() as tmp:
            c = Cache(Path(tmp) / "c.sqlite")
            c.put("k", "m", {"v": 42})
            hit = c.get("k", "m")
            check("Cache put/get round-trip",
                  hit is not None and hit.value == {"v": 42})
            c.close()
    except Exception as e:
        check(f"Cache round-trip  [{e}]", False)

    # --- report --------------------------------------------------
    print()
    color = GREEN if ok_count == total else (
        YELLOW if ok_count > total * 0.7 else RED)
    print(f"  {color}{ok_count}/{total} checks passed{RESET}")
    if ok_count < total:
        print(f"  {YELLOW}→ run 'conda activate cellsim' then "
              f"'cellsim doctor' again{RESET}")
        print(f"  {YELLOW}→ or re-create env: "
              f"'mamba env create -f environment.yml'{RESET}")
    else:
        print(f"  → ready to run: "
              f"'{GREEN}cellsim dock --smi benchmarks/dock/1stp_batch_5.smi "
              f"--receptor benchmarks/dock/1stp.pdb{RESET}'")
    print()
    return 0 if ok_count == total else 1


if __name__ == "__main__":
    sys.exit(main())
