#!/usr/bin/env python3
"""
setup_molgen.py — Create a Python virtualenv with ESMFold and RDKit
for molecular structure generation in CellSim.

Usage:
    python3 scripts/setup_molgen.py
"""

import subprocess
import sys
import os
import shutil

VENV_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".molgen_venv")
VENV_DIR = os.path.normpath(VENV_DIR)

def find_python():
    """Find a compatible Python (3.9-3.12) for rdkit-pypi."""
    for ver in ["3.12", "3.11", "3.10", "3.9"]:
        path = shutil.which(f"python{ver}")
        if path:
            return path
    # Fallback to current
    return sys.executable

def main():
    python = find_python()
    print(f"[setup_molgen] Using Python: {python}")
    print(f"[setup_molgen] Creating virtualenv at {VENV_DIR}")
    subprocess.check_call([python, "-m", "venv", VENV_DIR])

    pip = os.path.join(VENV_DIR, "bin", "pip")

    print("[setup_molgen] Installing RDKit...")
    subprocess.check_call([pip, "install", "rdkit-pypi", "numpy<2"])

    print("[setup_molgen] Installing ESMFold dependencies...")
    subprocess.check_call([pip, "install", "torch", "fair-esm"])

    print("[setup_molgen] Installing additional utilities...")
    subprocess.check_call([pip, "install", "numpy", "biopython"])

    print(f"[setup_molgen] Done. Virtualenv ready at {VENV_DIR}")
    print(f"[setup_molgen] Generate molecules with:")
    print(f"  {VENV_DIR}/bin/python scripts/generate_molecule.py <SMILES>")
    print(f"  {VENV_DIR}/bin/python scripts/generate_protein.py <sequence>")

if __name__ == "__main__":
    main()
