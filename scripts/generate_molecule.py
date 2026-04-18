#!/usr/bin/env python3
"""
generate_molecule.py — Generate 3D structure from SMILES using RDKit.
Outputs a .sdf file with 3D coordinates.

Usage:
    python generate_molecule.py <SMILES> <output_path>
    python generate_molecule.py "C10H12N5O13P3" assets/molecules/atp.sdf
"""

import sys
import os

def generate_3d(smiles: str, output_path: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[generate_molecule] ERROR: Invalid SMILES: {smiles}", file=sys.stderr)
        sys.exit(1)

    mol = Chem.AddHs(mol)
    try:
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    except TypeError:
        # Fallback for older rdkit versions
        result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    if result == -1:
        # Last resort: random coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)

    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    writer = Chem.SDWriter(output_path)
    writer.write(mol)
    writer.close()

    n_atoms = mol.GetNumAtoms()
    mw = Descriptors.MolWt(mol)
    print(f"[generate_molecule] {smiles} -> {output_path} ({n_atoms} atoms, MW={mw:.1f})")

def main():
    if len(sys.argv) < 3:
        print("Usage: generate_molecule.py <SMILES> <output.sdf>")
        sys.exit(1)
    generate_3d(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
