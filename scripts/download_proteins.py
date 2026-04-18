#!/usr/bin/env python3
"""
download_proteins.py — Download real protein structures from RCSB PDB
and extract Cα backbone traces for efficient rendering.

Outputs simplified PDB files with only Cα atoms (1 per residue).
"""

import os
import sys
import urllib.request

# PDB IDs and descriptions for the central dogma simulation
PROTEINS = {
    "1EHZ": {"name": "tRNA-Phe", "desc": "Yeast tRNA-Phenylalanine", "atom_filter": "P"},  # nucleic acid → use P atoms
    "4UG0": {"name": "ribosome_80S", "desc": "Human 80S ribosome", "atom_filter": "CA"},
    "1I6H": {"name": "rna_pol_II", "desc": "RNA Polymerase II", "atom_filter": "CA"},
    "5LJ3": {"name": "spliceosome", "desc": "Spliceosome U4/U6.U5 tri-snRNP", "atom_filter": "CA"},
    "1A3N": {"name": "hemoglobin", "desc": "Human hemoglobin (HBB product)", "atom_filter": "CA"},
    "5E84": {"name": "hsp70_chaperone", "desc": "Hsp70 chaperone (BiP)", "atom_filter": "CA"},
    "3KTQ": {"name": "dna_polymerase", "desc": "DNA Polymerase delta", "atom_filter": "CA"},
    "5ARA": {"name": "atp_synthase", "desc": "Mitochondrial ATP synthase", "atom_filter": "CA"},
}

def download_pdb(pdb_id: str, output_dir: str) -> str:
    """Download a PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    raw_path = os.path.join(output_dir, f"{pdb_id}_raw.pdb")

    if os.path.exists(raw_path):
        print(f"  [cached] {raw_path}")
        return raw_path

    print(f"  Downloading {url} ...")
    try:
        urllib.request.urlretrieve(url, raw_path)
        print(f"  Saved to {raw_path}")
    except Exception as e:
        # Try mmCIF format as fallback for large structures
        cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        print(f"  PDB failed, trying CIF: {cif_url}")
        try:
            cif_path = os.path.join(output_dir, f"{pdb_id}_raw.cif")
            urllib.request.urlretrieve(cif_url, cif_path)
            # Convert CIF to simplified PDB (just extract coords)
            return extract_from_cif(cif_path, raw_path)
        except Exception as e2:
            print(f"  ERROR: Could not download {pdb_id}: {e2}")
            return None

    return raw_path


def extract_from_cif(cif_path: str, pdb_path: str) -> str:
    """Extract ATOM records from mmCIF into PDB format."""
    with open(cif_path, 'r') as f:
        lines = f.readlines()

    atoms = []
    in_atom = False
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # mmCIF ATOM lines are space-separated
            parts = line.split()
            if len(parts) >= 15:
                try:
                    atom_name = parts[3]
                    res_name = parts[5]
                    chain = parts[6] if len(parts[6]) == 1 else parts[6][0]
                    x = float(parts[10])
                    y = float(parts[11])
                    z = float(parts[12])
                    element = parts[2] if len(parts) > 2 else atom_name[0]
                    atoms.append((atom_name, res_name, chain, x, y, z, element))
                except (ValueError, IndexError):
                    continue

    # Write as PDB
    with open(pdb_path, 'w') as f:
        for i, (aname, rname, chain, x, y, z, elem) in enumerate(atoms):
            f.write(f"ATOM  {i+1:5d}  {aname:<3s} {rname:>3s} {chain}{i//100:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n")
    print(f"  Converted CIF → {pdb_path} ({len(atoms)} atoms)")
    return pdb_path


def extract_backbone(raw_pdb_path: str, output_path: str, atom_filter: str = "CA") -> int:
    """Extract only backbone atoms (Cα for proteins, P for nucleic acids)."""
    if not raw_pdb_path or not os.path.exists(raw_pdb_path):
        return 0

    atoms = []
    with open(raw_pdb_path, 'r') as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != atom_filter:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = line[22:26].strip()
                element = line[76:78].strip() if len(line) > 77 else atom_filter[0]
                atoms.append((atom_name, res_name, chain, res_num, x, y, z, element))
            except (ValueError, IndexError):
                continue

    # Write simplified PDB
    with open(output_path, 'w') as f:
        f.write(f"REMARK  Backbone trace ({atom_filter} atoms only)\n")
        for i, (aname, rname, chain, rnum, x, y, z, elem) in enumerate(atoms):
            f.write(f"ATOM  {i+1:5d}  {aname:<3s} {rname:>3s} {chain}{rnum:>4s}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n")
        # Add CONECT records for backbone chain
        for i in range(1, len(atoms)):
            # Connect consecutive residues in same chain
            if atoms[i][2] == atoms[i-1][2]:  # same chain
                f.write(f"CONECT{i:5d}{i+1:5d}\n")
    return len(atoms)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.normpath(os.path.join(script_dir, "..", "assets", "proteins"))
    raw_dir = os.path.join(output_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    print("═══ Downloading protein structures from RCSB PDB ═══\n")

    for pdb_id, info in PROTEINS.items():
        print(f"[{pdb_id}] {info['name']} — {info['desc']}")

        # Download raw PDB
        raw_path = download_pdb(pdb_id, raw_dir)

        # Extract backbone
        out_path = os.path.join(output_dir, f"{info['name']}.pdb")
        n_atoms = extract_backbone(raw_path, out_path, info['atom_filter'])

        if n_atoms > 0:
            print(f"  ✓ Backbone: {n_atoms} atoms → {out_path}")
        else:
            print(f"  ✗ No atoms extracted")
        print()

    print("═══ Done ═══")
    print(f"Output: {output_dir}")


if __name__ == "__main__":
    main()
