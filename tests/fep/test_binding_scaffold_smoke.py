"""Layer 1.3 FEP Milestone B binding scaffold smoke.

Pins the Hamelberg–Gilson analytical restraint correction against
its closed-form Gaussian integral and verifies the scaffold builder
wires a valid alchemical system end-to-end (PDBFixer + Interchange
+ AbsoluteAlchemicalFactory + CustomCentroidBondForce) on a small
receptor (ubiquitin + methane, ~14k atoms).

Why both checks:
- Formula drift: if someone edits `_harmonic_restraint_free_energy_
  kcalmol` constants (kT, 1 M standard-state volume) the computed
  ΔG_R silently changes — the analytical check traps that.
- API drift: the builder depends on openff-toolkit Topology.from_pdb
  + Interchange.from_smirnoff + openmmtools.alchemy. Any of those
  APIs breaks and this gate catches it before production runs.

Kept out of the headless validators tree on purpose: this test
needs a working conda env (openff, openmmtools, pdbfixer) and is
slow enough (~20 s on CPU) that it doesn't belong in a fast-loop
unit suite. Run it before cutting a release.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.fep.binding import (  # noqa: E402
    _harmonic_restraint_free_energy_kcalmol,
    _build_complex_alchemical_system,
    _build_complex_alchemical_system_amber14,
    _build_solvent_alchemical_system_amber14,
    _prepare_protein_topology,
)


def _nonbonded_cutoff_nm(system):
    """Cutoff (nm) of the system's NonbondedForce, or None."""
    from openmm import NonbondedForce, unit
    for f in system.getForces():
        if isinstance(f, NonbondedForce):
            return f.getCutoffDistance().value_in_unit(unit.nanometer)
    return None


def test_amber14_solvent_leg_matches_complex_leg_nonbonded():
    """DDM cancellation requires the solvent leg and the complex leg
    to decouple the ligand from the SAME water model at the SAME
    cutoff. Before BUG_AUDIT.md #5 the solvent leg ran through the
    pure-SMIRNOFF builder (tip3p.offxml, 0.9 nm switched cutoff) while
    the amber14 complex leg used amber14/tip3pfb at 1.0 nm PME. Pin the
    cutoffs equal so that mismatch can't silently return."""
    # padding_nm=0.8 deliberately: the builder must floor the solvent
    # box so the 1.0 nm cutoff stays ≤ half the box. A small ligand at
    # tight padding otherwise trips OpenMM's "cutoff > half box" at
    # sampler minimisation (caught end-to-end, not by scaffold build).
    solv_alch, top, _pos, n_lig = (
        _build_solvent_alchemical_system_amber14("C", padding_nm=0.8))
    assert n_lig == 5, f"methane should have 5 atoms; got {n_lig}"
    assert solv_alch.getNumParticles() > 100, (
        "solvent leg should be a solvated box, not bare methane")
    solv_cut = _nonbonded_cutoff_nm(solv_alch)
    assert solv_cut is not None and abs(solv_cut - 1.0) < 1e-6, (
        f"solvent-leg cutoff {solv_cut} nm != 1.0 nm")
    # Every box vector must exceed 2×cutoff or PME rejects the system.
    import numpy as np
    from openmm import unit as _u
    box = np.asarray(
        solv_alch.getDefaultPeriodicBoxVectors().value_in_unit(
            _u.nanometer) if hasattr(
            solv_alch.getDefaultPeriodicBoxVectors(), "value_in_unit")
        else [[v.x, v.y, v.z] for v in
              solv_alch.getDefaultPeriodicBoxVectors()])
    min_span = min(box[i][i] for i in range(3))
    assert min_span > 2.0 * solv_cut, (
        f"solvent box min span {min_span:.2f} nm <= 2×cutoff "
        f"{2*solv_cut:.2f} nm — PME will reject this at sampling time")

    pdb = REPO_ROOT / "benchmarks/md/1ubq.pdb"
    cx = _build_complex_alchemical_system_amber14("C", pdb, padding_nm=0.8)
    cx_cut = _nonbonded_cutoff_nm(cx["alch_system"])
    assert cx_cut is not None and abs(solv_cut - cx_cut) < 1e-6, (
        f"leg cutoff mismatch: solvent {solv_cut} vs complex "
        f"{cx_cut} nm — DDM water contribution won't cancel")


def test_harmonic_restraint_correction_matches_closed_form():
    """Cross-check `_harmonic_restraint_free_energy_kcalmol` against
    the Gaussian integral written out by hand. The helper returns the
    confinement / standard-state correction to ADD to ΔG_bind:
    +kT·ln(V_std/V_harm) > 0 (BUG_AUDIT.md #1/#3)."""
    T_K = 298.15
    k_kJ_per_nm2 = 4184.0           # 10 kcal/mol/Å²
    kT_kJ = 0.0083144626 * T_K      # 2.479 kJ/mol
    # Volume of an isotropic 3D Gaussian: (2π·kT/k)^{3/2}
    v_harm_nm3 = (2.0 * math.pi * kT_kJ / k_kJ_per_nm2) ** 1.5
    v_std_nm3 = 1.66054             # 1 M per-molecule volume
    # POSITIVE confinement work — this is the additive correction, not
    # the (negative) release free energy.
    dG_kJ = kT_kJ * math.log(v_std_nm3 / v_harm_nm3)
    expected_kcalmol = dG_kJ / 4.184

    got = _harmonic_restraint_free_energy_kcalmol(k_kJ_per_nm2)
    assert got > 0, (
        f"standard-state correction must be POSITIVE (confinement "
        f"penalty); got {got:.4f}. A negative value is the reverted "
        "BUG_AUDIT #1/#3 sign error.")
    assert abs(got - expected_kcalmol) < 1e-6, (
        f"correction drift: got {got:.6f}, "
        f"expected {expected_kcalmol:.6f}")


def test_correction_scales_monotonically_with_k():
    """Tighter spring → smaller restraint volume → LARGER positive
    standard-state correction. If this breaks, the sign convention in
    the formula has flipped (BUG_AUDIT.md #1/#3)."""
    tight = _harmonic_restraint_free_energy_kcalmol(4184.0)
    loose = _harmonic_restraint_free_energy_kcalmol(418.4)
    # Both positive (confinement cost). Tighter spring confines the
    # ligand to a smaller volume → larger entropy penalty.
    assert tight > 0 and loose > 0
    assert tight > loose, (
        f"expected tight ({tight:.2f}) MORE positive than loose "
        f"({loose:.2f}); the spring-k monotonicity is broken")


def test_restraint_correction_matches_openmmtools_ssc():
    """Ground-truth the standard-state correction against openmmtools'
    HarmonicRestraintForce.compute_standard_state_correction. The DDM
    cycle adds the CONFINEMENT work, which is the negative of
    openmmtools' (release-direction) SSC. This is the load-bearing
    regression for BUG_AUDIT.md #1/#3 — it pins the sign to an
    independent, YANK-validated reference rather than our own algebra."""
    from openmm import System, NonbondedForce, unit
    from openmmtools import forces, states

    T_K = 298.15
    k_kJ_per_nm2 = 4184.0
    kT_kcal = 0.0083144626 * T_K / 4.184

    sys_ = System()
    sys_.addParticle(12.0 * unit.dalton)
    sys_.addParticle(12.0 * unit.dalton)
    nb = NonbondedForce()
    nb.addParticle(0.0, 0.1, 0.0)
    nb.addParticle(0.0, 0.1, 0.0)
    sys_.addForce(nb)
    restraint = forces.HarmonicRestraintForce(
        spring_constant=(k_kJ_per_nm2
                         * unit.kilojoule_per_mole / unit.nanometer**2),
        restrained_atom_indices1=[0],
        restrained_atom_indices2=[1])
    sys_.addForce(restraint)
    ts = states.ThermodynamicState(
        system=sys_, temperature=T_K * unit.kelvin)
    ssc_kcal = restraint.compute_standard_state_correction(ts) * kT_kcal

    ours = _harmonic_restraint_free_energy_kcalmol(k_kJ_per_nm2)
    # openmmtools reports the release direction (negative); our
    # additive correction is its negation (positive).
    assert ours > 0 > ssc_kcal, (
        f"sign convention broken: ours={ours:+.4f} (want > 0), "
        f"openmmtools SSC={ssc_kcal:+.4f} (want < 0)")
    assert abs(ours - (-ssc_kcal)) < 0.05, (
        f"helper {ours:+.4f} != -(openmmtools SSC) {-ssc_kcal:+.4f} "
        "kcal/mol")


def test_restraint_correction_is_r0_aware():
    """BUG_AUDIT.md #9/#10: the standard-state correction must use the
    r0-aware restraint volume. At r0=0 it equals the Gaussian value;
    for r0>0 the accessible shell is larger so the positive correction
    shrinks monotonically."""
    from src.fep.binding import (
        _harmonic_restraint_free_energy_kcalmol as corr,
        _harmonic_restraint_volume_nm3 as vol,
    )
    k = 4184.0
    # r0=0 reduces to the isotropic Gaussian volume (2π kT/k)^1.5.
    kT_kJ = 0.0083144626 * 298.15
    v_gauss = (2.0 * math.pi * kT_kJ / k) ** 1.5
    assert abs(vol(k, r0_nm=0.0) - v_gauss) < 1e-9

    c0 = corr(k, r0_nm=0.0)
    c02 = corr(k, r0_nm=0.2)
    c05 = corr(k, r0_nm=0.5)
    # All positive (confinement), shrinking as the shell grows.
    assert c0 > c02 > c05 > 0, (c0, c02, c05)
    # Ballpark from the verified closed form.
    assert abs(c0 - 5.27) < 0.05
    assert abs(c05 - 1.28) < 0.05


def test_restraint_r0_from_geometry_is_placement_to_anchor_distance():
    from src.fep.binding import _restraint_r0_from_geometry
    import numpy as np
    # Three anchor atoms with centroid at origin; placement 0.3 nm away.
    pos = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    anchor = [0, 1, 2]              # centroid = (0,0,0)
    placement = np.array([0.3, 0.0, 0.0])
    r0 = _restraint_r0_from_geometry(pos, anchor, placement)
    assert abs(r0 - 0.3) < 1e-9


def test_placed_pose_sits_near_restraint_minimum():
    """BUG_AUDIT.md #7: with r0 = placement→anchor distance the docked
    ligand starts at (not ~125 kcal/mol above) the restraint minimum.
    Verify the built pose's ligand-centroid→anchor-centroid distance
    equals the restraint r0, so the restraint energy there is ~0 and
    never exceeds what the old r0=0 choice would have cost."""
    import numpy as np
    from openmm import unit as u
    pdb = REPO_ROOT / "benchmarks/md/1ubq.pdb"
    cx = _build_complex_alchemical_system_amber14("C", pdb, padding_nm=0.8)
    pos = np.asarray(cx["positions"].value_in_unit(u.nanometer))
    lig_c = pos[cx["ligand_indices"]].mean(axis=0)
    anch_c = pos[cx["anchor_indices"]].mean(axis=0)
    r = float(np.linalg.norm(lig_c - anch_c))
    r0 = cx["restraint_r0_nm"]
    k = cx["restraint_k_kJ_per_nm2"]
    e_fixed = 0.5 * k * (r - r0) ** 2 / 4.184        # kcal/mol
    e_r0_zero = 0.5 * k * r ** 2 / 4.184             # the old-bug energy
    assert e_fixed < 5.0, (
        f"placed-pose restraint energy {e_fixed:.2f} kcal/mol too high "
        f"(r={r:.3f}, r0={r0:.3f}); ligand not at restraint minimum")
    assert e_fixed <= e_r0_zero + 1e-6, (
        "geometry-derived r0 should never cost more than r0=0")


def test_complex_alchemical_builder_on_ubiquitin():
    """Scaffold build ubiquitin + methane complex — verifies the
    PDBFixer → Topology.from_pdb → Interchange → alchemical factory
    + restraint pipeline."""
    pdb = REPO_ROOT / "benchmarks/md/1ubq.pdb"
    out = _build_complex_alchemical_system(
        "C", pdb, padding_nm=0.8)
    assert out["n_ligand_atoms"] == 5, (
        f"methane should have 5 atoms; got {out['n_ligand_atoms']}")
    # Ubiquitin has 76 residues; after PDBFixer hydrogenation we
    # expect on the order of ~1200 protein atoms (ubiquitin canonical
    # is 1231 atoms per the PDBFixer output).
    assert 1000 <= out["n_protein_atoms"] <= 1500, (
        f"unexpected protein atom count {out['n_protein_atoms']}")
    assert out["n_total_atoms"] > 10_000, (
        "solvated box should be > 10k atoms at 0.8 nm padding")
    # Ligand indices contiguous at end of protein.
    assert out["ligand_indices"] == list(
        range(out["n_protein_atoms"],
              out["n_protein_atoms"] + out["n_ligand_atoms"]))
    assert len(out["anchor_indices"]) >= 5, (
        "restraint anchor should have at least 5 Cα atoms even on "
        "a small protein — fallthrough to the nearest-10 pick "
        "isn't kicking in")


def test_amber14_builder_on_streptavidin():
    """Phase-2 regression guard. Before the terminal-missing-residue
    filter + amber14 builder landed (commit 7f1b388), streptavidin
    1stp scaffold OOM'd at 2 GB / 30 min. After: builds in < 30 s
    with a solvated system < 100 k atoms. If either number regresses
    past those thresholds, something reverted."""
    import time
    pdb = REPO_ROOT / "benchmarks/dock/1stp.pdb"
    t0 = time.time()
    out = _build_complex_alchemical_system_amber14(
        "C", pdb, padding_nm=0.8)
    elapsed = time.time() - t0
    assert elapsed < 30.0, (
        f"amber14 builder on 1stp took {elapsed:.1f}s — "
        "terminal-missing-residue filter may have regressed")
    assert out["n_total_atoms"] < 100_000, (
        f"1stp solvated system has {out['n_total_atoms']} atoms; "
        "expected < 100k. If > 500k, the terminal-residue bloat "
        "returned.")
    assert out["n_ligand_atoms"] == 5
    assert 1500 < out["n_protein_atoms"] < 2500, (
        f"1stp post-PDBFixer protein should be ~1 744 atoms; "
        f"got {out['n_protein_atoms']}. If > 2 200 the terminal "
        "residues are being re-added.")
    assert "amber14" in out["ff_stack"]


def test_terminal_missing_residue_filter():
    """Unit-level check: _prepare_protein_topology must NOT
    introduce atoms > 5 nm from the main protein body on 1stp.
    That was the exact failure mode before the filter landed."""
    import numpy as np
    from openmm import unit as u
    pdb = REPO_ROOT / "benchmarks/dock/1stp.pdb"
    tmp_path, top_omm, positions = _prepare_protein_topology(pdb)
    import os
    try:
        os.unlink(tmp_path)
    except OSError:
        pass
    pos_nm = np.asarray(
        positions.value_in_unit(u.nanometer), dtype=float)
    extent = pos_nm.max(axis=0) - pos_nm.min(axis=0)
    # Streptavidin monomer is ~3.5 × 3.5 × 5 nm — any dimension
    # > 6 nm means a mis-fitted terminal residue snuck through.
    assert max(extent) < 6.0, (
        f"1stp post-PDBFixer extent {extent} nm — terminal "
        "mis-fit has returned, solvation box will bloat 50×")


def test_compute_absolute_binding_dg_force_field_path_flag():
    """Verify the force_field_path kwarg dispatches to the right
    builder. Both should scaffold_both_legs ok; timing differs."""
    from src.fep.binding import compute_absolute_binding_dg
    pdb = REPO_ROOT / "benchmarks/md/1ubq.pdb"

    r_amber = compute_absolute_binding_dg(
        "C", pdb, padding_nm=0.8,
        force_field_path="amber14",
        sample=False)
    assert r_amber.ok
    assert r_amber.phase == "scaffolded_both_legs"

    r_smirn = compute_absolute_binding_dg(
        "C", pdb, padding_nm=0.8,
        force_field_path="smirnoff",
        sample=False)
    assert r_smirn.ok
    assert r_smirn.phase == "scaffolded_both_legs"

    # Bogus value errors cleanly, not crashes.
    r_bad = compute_absolute_binding_dg(
        "C", pdb, padding_nm=0.8,
        force_field_path="not_a_thing",
        sample=False)
    assert not r_bad.ok
    assert "unknown force_field_path" in r_bad.reason


if __name__ == "__main__":
    funcs = [
        ("restraint correction closed-form",
         test_harmonic_restraint_correction_matches_closed_form),
        ("restraint k monotonicity",
         test_correction_scales_monotonically_with_k),
        ("restraint correction vs openmmtools SSC",
         test_restraint_correction_matches_openmmtools_ssc),
        ("amber14 solvent leg nonbonded matches complex leg",
         test_amber14_solvent_leg_matches_complex_leg_nonbonded),
        ("complex builder on 1ubq + methane",
         test_complex_alchemical_builder_on_ubiquitin),
        ("amber14 builder on 1stp (Phase-2 regression guard)",
         test_amber14_builder_on_streptavidin),
        ("terminal-missing-residue filter on 1stp",
         test_terminal_missing_residue_filter),
        ("compute_absolute_binding_dg force_field_path flag",
         test_compute_absolute_binding_dg_force_field_path_flag),
    ]
    fails = []
    try:
        for label, f in funcs:
            try:
                f()
                print(f"[PASS] {label}")
            except AssertionError as e:
                print(f"[FAIL] {label}: {e}")
                fails.append(label)
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"[ERROR] {label}: {e}")
                fails.append(label)
        if fails:
            sys.exit(1)
        print(f"{len(funcs)}/{len(funcs)} PASS")
    except AssertionError as e:
        print(f"[FAIL] {e}")
        sys.exit(1)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[ERROR] {e}")
        sys.exit(2)
