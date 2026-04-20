"""Non-AI ADMET descriptors for drug-likeness triage.

The 10 properties every med-chem grad student is expected to know
about a candidate compound, all computed from RDKit-native
descriptors — no learned model anywhere.

Reported properties:

  Physicochemical:
    - MW         molecular weight (Da)
    - logP       Crippen octanol-water partition coefficient
    - tpsa       topological polar surface area (Å²; Ertl 2000)
    - hba        H-bond acceptor count (Lipinski definition)
    - hbd        H-bond donor count (Lipinski definition)
    - rotb       rotatable bond count
    - aromatic_rings
    - heavy_atoms
    - formal_charge
  Drug-likeness summary:
    - ro5_violations  Lipinski rule-of-five violations (0-4)
    - ro5_pass        True iff MW<=500 & logP<=5 & HBA<=10 & HBD<=5
    - qed             Bickerton 2012 QED score ∈ [0, 1]
  Solubility (empirical Delaney 2004 ESOL equation):
    - logS_ESOL   log of molar solubility
    - solubility_classification {highly soluble, soluble,
                                 moderately soluble, slightly
                                 soluble, insoluble}

All formulas are published literature expressions evaluated over
RDKit features. No neural network, no training set. Suitable for
the Campaign-1 triage stage and for sanity-checking compounds that
come out of Vina / xTB pipelines.

References:
  - Wildman & Crippen 1999 (J Chem Inf Comput Sci 39:868) — logP.
  - Ertl, Rohde, Selzer 2000 (J Med Chem 43:3714) — TPSA.
  - Lipinski, Lombardo, Dominy, Feeney 1997 (Adv Drug Deliv
    Rev 23:3) — rule of five.
  - Bickerton et al 2012 (Nat Chem 4:90) — QED.
  - Delaney 2004 (J Chem Inf Comput Sci 44:1000) — ESOL.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional


@dataclass
class AdmetResult:
    """Envelope for a non-AI ADMET descriptor bundle."""

    smiles: str
    ok: bool
    reason: str = ""
    formula: Optional[str] = None
    inchi_key: Optional[str] = None

    # Physicochemistry
    MW: Optional[float] = None
    logP: Optional[float] = None
    tpsa: Optional[float] = None
    hba: Optional[int] = None
    hbd: Optional[int] = None
    rotb: Optional[int] = None
    aromatic_rings: Optional[int] = None
    heavy_atoms: Optional[int] = None
    formal_charge: Optional[int] = None

    # Drug-likeness
    ro5_violations: Optional[int] = None
    ro5_pass: Optional[bool] = None
    qed: Optional[float] = None

    # Solubility
    logS_ESOL: Optional[float] = None
    solubility_class: Optional[str] = None

    tool_versions: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] ADMET {self.smiles}  {self.reason}"
        ro5 = ("✓" if self.ro5_pass else f"✗ (×{self.ro5_violations})")
        return (f"[OK]   ADMET {self.formula or self.smiles}  "
                f"MW={self.MW:.1f}  logP={self.logP:+.2f}  "
                f"TPSA={self.tpsa:.1f}  HBA={self.hba}  HBD={self.hbd}  "
                f"rotb={self.rotb}  Ro5 {ro5}  QED={self.qed:.2f}  "
                f"logS={self.logS_ESOL:+.2f} ({self.solubility_class})")


def _classify_solubility(log_s: float) -> str:
    """Delaney-style qualitative classification.

    logS ranges (mol/L):  highly > -1 > soluble > -2 > moderately >
    -4 > slightly > -6 > insoluble.
    """
    if log_s > -1.0:
        return "highly soluble"
    if log_s > -2.0:
        return "soluble"
    if log_s > -4.0:
        return "moderately soluble"
    if log_s > -6.0:
        return "slightly soluble"
    return "insoluble"


def _tool_versions() -> dict:
    try:
        import rdkit
        return {"rdkit": rdkit.__version__}
    except Exception:
        return {}


def compute_admet(smiles: str) -> AdmetResult:
    """Compute all ADMET descriptors for one SMILES."""
    result = AdmetResult(smiles=smiles, ok=False,
                         tool_versions=_tool_versions())

    try:
        from rdkit import Chem
        from rdkit.Chem import (
            AllChem, Crippen, Descriptors, Lipinski,
            rdMolDescriptors, QED,
        )
    except ImportError as e:
        result.reason = f"rdkit import failed: {e}"
        return result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result.reason = "RDKit could not parse SMILES"
        return result

    mol_h = Chem.AddHs(mol)

    try:
        result.formula = rdMolDescriptors.CalcMolFormula(mol)
        result.inchi_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))

        result.MW = float(Descriptors.MolWt(mol))
        result.logP = float(Crippen.MolLogP(mol))
        result.tpsa = float(rdMolDescriptors.CalcTPSA(mol))
        result.hba = int(Lipinski.NumHAcceptors(mol))
        result.hbd = int(Lipinski.NumHDonors(mol))
        result.rotb = int(Descriptors.NumRotatableBonds(mol))
        result.aromatic_rings = int(
            rdMolDescriptors.CalcNumAromaticRings(mol))
        result.heavy_atoms = int(mol.GetNumHeavyAtoms())
        result.formal_charge = int(Chem.GetFormalCharge(mol))

        # Lipinski Rule of Five: MW<=500, logP<=5, HBA<=10, HBD<=5
        violations = 0
        if result.MW > 500.0:     violations += 1
        if result.logP > 5.0:     violations += 1
        if result.hba > 10:       violations += 1
        if result.hbd > 5:        violations += 1
        result.ro5_violations = violations
        result.ro5_pass = (violations == 0)

        # QED (Bickerton 2012)
        try:
            result.qed = float(QED.qed(mol))
        except Exception:
            result.qed = None

        # ESOL (Delaney 2004):
        #   logS = 0.16 - 0.63*logP - 0.0062*MW + 0.066*rotb
        #          - 0.74*aromatic_atom_proportion
        n_atoms = max(1, mol.GetNumHeavyAtoms())
        n_aro = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        aromatic_proportion = n_aro / n_atoms
        log_s = (0.16
                 - 0.63 * result.logP
                 - 0.0062 * result.MW
                 + 0.066 * result.rotb
                 - 0.74 * aromatic_proportion)
        result.logS_ESOL = float(log_s)
        result.solubility_class = _classify_solubility(log_s)

    except Exception as e:
        result.reason = f"descriptor compute failed: {e}"
        return result

    result.ok = True
    return result


if __name__ == "__main__":
    import argparse
    import sys

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    r = compute_admet(args.smiles)
    if args.json:
        import json
        print(json.dumps(r.as_dict(), indent=2, default=str))
    else:
        print(r.summary())
    sys.exit(0 if r.ok else 1)
