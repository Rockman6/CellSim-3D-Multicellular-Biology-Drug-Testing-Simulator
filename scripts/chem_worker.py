#!/usr/bin/env python3
"""chem_worker.py — Background compute queue for the bioagent chemistry stack.

Watches `data/bioagents/chem/_pending/` for tickets produced by the sim UI.
Each ticket is a JSON file naming a drug id + requested tier (0..4).
The worker consumes tiers in order, writing each tier's output into
`data/bioagents/chem/<drug_id>/tier<N>_*.json`. The sim polls the output
directory and picks up new results without blocking.

Usage:
    python scripts/chem_worker.py [--data-dir data/bioagents]
                                  [--gpu mps|cuda|cpu]
                                  [--interval 2]

Tiers and the Python packages each requires:
    T0: rdkit               (descriptors, Morgan FP, 3D pharmacophore)
    T1: openff-toolkit      (partial charges + LJ)
    T2: xtb-python          (HOMO/LUMO, dipole, reaction energies)
    T3: gnina binary        (CNN-rescored docking)
    T4: openmm + mace-torch (binding MD)

Any missing package makes that tier fall through silently — the sim keeps
running on whatever tiers are available. Tier 0 always works because RDKit
is already a project dependency.

This is a minimal stub — real tier implementations land in later phases.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time
from pathlib import Path

# ── Tier availability probes (lazy-imported) ──────────────────────────
def _probe(name: str) -> bool:
    try:
        __import__(name)
        return True
    except ImportError:
        return False

TIERS = {
    0: {"name": "RDKit descriptors",        "probe": "rdkit"},
    1: {"name": "OpenFF parametrisation",   "probe": "openff.toolkit"},
    2: {"name": "xTB semi-empirical QM",    "probe": "xtb"},
    3: {"name": "GNINA docking",            "probe": "gnina"},
    4: {"name": "OpenMM + MACE MD",         "probe": "openmm"},
}


def ensure_dirs(data_dir: Path):
    (data_dir / "chem" / "_pending").mkdir(parents=True, exist_ok=True)
    (data_dir / "chem" / "_done").mkdir(parents=True, exist_ok=True)


def consume_ticket(ticket_path: Path, data_dir: Path) -> None:
    """Process one JSON ticket. Ticket schema:
        {"drug_id": "cisplatin", "smiles": "N.N.Cl[Pt]Cl",
         "tiers": [0, 1, 2, 3], "targets": ["TGT_CDK1", ...]}
    """
    try:
        ticket = json.loads(ticket_path.read_text())
    except Exception as e:
        print(f"[chem_worker] Ticket {ticket_path.name} unreadable: {e}",
              file=sys.stderr)
        return

    drug_id = ticket["drug_id"]
    smiles  = ticket.get("smiles", "")
    tiers   = sorted(ticket.get("tiers", [0]))
    out_dir = data_dir / "chem" / drug_id
    out_dir.mkdir(parents=True, exist_ok=True)

    for tier in tiers:
        handler = HANDLERS.get(tier)
        if handler is None:
            continue
        try:
            handler(drug_id, smiles, ticket, out_dir)
        except Exception as e:
            print(f"[chem_worker] Tier {tier} failed for {drug_id}: {e}",
                  file=sys.stderr)

    # Move the ticket to _done so we don't reprocess.
    done = data_dir / "chem" / "_done" / ticket_path.name
    ticket_path.rename(done)


# ── Tier 0 — RDKit descriptors ────────────────────────────────────────
def tier0_rdkit(drug_id: str, smiles: str, ticket: dict, out: Path) -> None:
    """Compute descriptor bundle + 2D/3D fingerprint + pharmacophore."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
        from rdkit.Chem import rdPartialCharges
    except ImportError:
        print("[chem_worker] rdkit unavailable; skipping tier 0",
              file=sys.stderr)
        return
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    molH = Chem.AddHs(mol)
    AllChem.EmbedMolecule(molH, randomSeed=42, useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(molH, maxIters=500)
    except Exception:
        pass
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    fp_hex = fp.ToBitString()  # long string; caller can pack

    data = {
        "drug_id": drug_id,
        "smiles":  smiles,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw":      float(Descriptors.MolWt(mol)),
        "logP":    float(Descriptors.MolLogP(mol)),
        "tpsa":    float(Descriptors.TPSA(mol)),
        "hbd":     int(Descriptors.NumHDonors(mol)),
        "hba":     int(Descriptors.NumHAcceptors(mol)),
        "rotatable_bonds": int(Descriptors.NumRotatableBonds(mol)),
        "aromatic_rings":  int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "formal_charge":   int(Chem.GetFormalCharge(mol)),
        "num_atoms":       int(mol.GetNumAtoms()),
        "morgan_fp_2048_bits": fp_hex,
    }
    (out / "tier0_descriptors.json").write_text(json.dumps(data, indent=2))

    # Also write the 3D SDF for the rest of the pipeline to consume.
    assets = Path(__file__).parent.parent / "assets" / "drugs"
    assets.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(assets / f"{drug_id}.sdf"))
    writer.write(molH)
    writer.close()
    print(f"[chem_worker] tier0 done: {drug_id}")


# Later tiers are stubs for now — implemented in later phases.
def tier1_openff(drug_id, smiles, ticket, out):
    pass

def tier2_xtb(drug_id, smiles, ticket, out):
    pass

def tier3_gnina(drug_id, smiles, ticket, out):
    pass

def tier4_openmm(drug_id, smiles, ticket, out):
    pass


HANDLERS = {
    0: tier0_rdkit,
    1: tier1_openff,
    2: tier2_xtb,
    3: tier3_gnina,
    4: tier4_openmm,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/bioagents")
    ap.add_argument("--interval", type=float, default=2.0,
                    help="polling interval in seconds")
    ap.add_argument("--once", action="store_true",
                    help="process current pending tickets once and exit")
    args = ap.parse_args()

    data_dir = Path(args.data_dir).resolve()
    ensure_dirs(data_dir)

    # Report tier availability.
    print("[chem_worker] Bioagent chemistry worker starting.")
    print(f"[chem_worker] Data dir: {data_dir}")
    print("[chem_worker] Tier availability:")
    for t, info in TIERS.items():
        ok = _probe(info["probe"])
        print(f"  T{t} {info['name']:30s}  "
              f"{'OK' if ok else 'MISSING (fall-through)'}")

    pending = data_dir / "chem" / "_pending"
    while True:
        for ticket in sorted(pending.glob("*.json")):
            consume_ticket(ticket, data_dir)
        if args.once:
            break
        time.sleep(args.interval)
    return 0


if __name__ == "__main__":
    sys.exit(main())
