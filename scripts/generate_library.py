#!/usr/bin/env python3
"""
generate_library.py — Pre-generate the standard molecule library for CellSim.
Generates 3D structures for common cellular molecules using RDKit.

Usage:
    python generate_library.py [output_dir]
    Default output: assets/molecules/
"""

import sys
import os
import importlib

# Common cellular molecules with SMILES
MOLECULE_LIBRARY = {
    # Energy molecules
    "atp": {
        "smiles": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
        "name": "Adenosine Triphosphate (ATP)"
    },
    "adp": {
        "smiles": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
        "name": "Adenosine Diphosphate (ADP)"
    },
    "nadh": {
        "smiles": "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@@H](COP(=O)([O-])OP(=O)([O-])OC[C@H]2OC(n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",
        "name": "NADH"
    },
    "glucose": {
        "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "name": "D-Glucose"
    },

    # Nucleotides
    "deoxyadenosine": {
        "smiles": "c1nc(c2c(n1)n(cn2)[C@@H]3C[C@@H]([C@H](O3)CO)O)N",
        "name": "Deoxyadenosine"
    },
    "thymidine": {
        "smiles": "Cc1cn([C@@H]2C[C@@H]([C@H](O2)CO)O)c(=O)[nH]c1=O",
        "name": "Thymidine"
    },

    # Lipids (membrane components)
    "dppc": {
        "smiles": "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC",
        "name": "DPPC (Phosphatidylcholine)"
    },
    "cholesterol": {
        "smiles": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
        "name": "Cholesterol"
    },

    # Amino acids (representative set)
    "alanine": {
        "smiles": "C[C@@H](C(=O)O)N",
        "name": "L-Alanine"
    },
    "glycine": {
        "smiles": "NCC(=O)O",
        "name": "Glycine"
    },
    "tryptophan": {
        "smiles": "c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N",
        "name": "L-Tryptophan"
    },

    # Signaling molecules
    "calcium_ion": {
        "smiles": "[Ca+2]",
        "name": "Calcium Ion (Ca2+)"
    },
    "camp": {
        "smiles": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N",
        "name": "Cyclic AMP"
    },

    # ROS
    "water": {
        "smiles": "O",
        "name": "Water"
    },
}


def generate_all(output_dir: str):
    # Import generate_molecule from sibling script
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from generate_molecule import generate_3d

    os.makedirs(output_dir, exist_ok=True)
    success = 0
    failed = 0

    for mol_id, info in MOLECULE_LIBRARY.items():
        out_path = os.path.join(output_dir, f"{mol_id}.sdf")
        try:
            generate_3d(info["smiles"], out_path)
            success += 1
        except Exception as e:
            print(f"[generate_library] FAILED: {mol_id} ({info['name']}): {e}")
            failed += 1

    print(f"\n[generate_library] Done: {success} generated, {failed} failed")
    print(f"[generate_library] Output: {output_dir}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_dir = os.path.join(script_dir, "..", "assets", "molecules")
    output_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.normpath(default_dir)
    generate_all(output_dir)


if __name__ == "__main__":
    main()
