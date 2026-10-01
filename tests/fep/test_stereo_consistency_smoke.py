"""Layer 1.3 FEP — stereo consistency across legs (BUG_AUDIT.md #6).

The binding complex builders pass allow_undefined_stereo=True, but the
solvent / hydration path (`ligand_hydration_fep`, `_build_alchemical_
legs`) did not. A SMILES with an undefined stereocentre therefore built
the complex leg fine and then raised UndefinedStereochemistryError on
the solvent leg — an asymmetric failure that killed any ΔG_bind /
ΔG_hyd on such an input mid-run.

We use alanine (`CC(N)C(=O)O`, one undefined tetrahedral centre), which
raises without the flag and builds with it. These tests do a real
OpenFF build (~tens of seconds) so they live in the slow smoke tier.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

UNDEF_STEREO_SMILES = "CC(N)C(=O)O"   # alanine — undefined stereocentre


def test_build_alchemical_legs_accepts_undefined_stereo():
    """`_build_alchemical_legs` is the production solvent/vacuum leg
    for both binding (smirnoff path) and hydration. It must not raise
    on an undefined-stereo ligand."""
    from src.fep import _build_alchemical_legs
    out = _build_alchemical_legs(UNDEF_STEREO_SMILES)
    # (vac_alch, solv_alch, vac_top, solv_top, vac_pos, solv_pos, n)
    assert len(out) == 7
    vac_alch, solv_alch = out[0], out[1]
    assert vac_alch.getNumParticles() > 0
    assert solv_alch.getNumParticles() > vac_alch.getNumParticles(), (
        "solvated system should have more atoms than vacuum")


def test_hydration_scaffold_accepts_undefined_stereo():
    """The Phase-1 hydration scaffold must also accept undefined
    stereo (same fix site, __init__.py:219)."""
    from src.fep import ligand_hydration_fep
    r = ligand_hydration_fep(UNDEF_STEREO_SMILES)
    assert r.ok, f"scaffold failed on undefined-stereo ligand: {r.reason}"
    assert "both_legs" in (r.phase or "") or r.phase == "scaffolded", (
        f"unexpected phase {r.phase!r}")


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
