#!/usr/bin/env python3
"""
generate_protein.py — Predict protein 3D structure from amino acid sequence using ESMFold.
Outputs a .pdb file.

Usage:
    python generate_protein.py <sequence> <output_path>
    python generate_protein.py "MKFLILLFNILCLFPVLAADNHGVS" assets/molecules/protein.pdb
"""

import sys
import os
import torch

def predict_structure(sequence: str, output_path: str):
    import esm

    print(f"[generate_protein] Loading ESMFold model...")
    model = esm.pretrained.esmfold_v1()
    model = model.eval()

    # Use CPU on Mac (MPS support is limited for ESMFold)
    device = torch.device("cpu")
    if torch.cuda.is_available():
        device = torch.device("cuda")
    model = model.to(device)

    print(f"[generate_protein] Predicting structure for {len(sequence)} residues...")
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write(output)

    print(f"[generate_protein] {len(sequence)} residues -> {output_path}")

def main():
    if len(sys.argv) < 3:
        print("Usage: generate_protein.py <amino_acid_sequence> <output.pdb>")
        sys.exit(1)
    predict_structure(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
