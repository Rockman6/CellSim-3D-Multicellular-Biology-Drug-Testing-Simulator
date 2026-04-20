"""OpenMM ligand-in-vacuum MD driver.

Consumes Layer 1.1's `ParametrizeResult`, rebuilds the OpenMM System
(parametrize.py intentionally did not preserve System objects across
its Python boundary — they are not picklable in all OpenMM versions
and we want the Layer 1.2 driver to own its integrator / platform /
options anyway), minimises energy, and runs NVT Langevin dynamics.

Telemetry per stored frame:
    - time (ps)
    - potential energy (kJ/mol)
    - kinetic energy (kJ/mol)
    - temperature (K)
    - RMSD vs frame 0 (Å, all-atom)
    - positions (Å)

Gate-quality checks on the output:
    - no NaN / inf anywhere
    - final temperature within ±50 K of setpoint (crude equilibration
      check)
    - no dramatic geometric blow-up (RMSD stays < 10 Å vs frame 0)
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Optional
import logging
import math
import sys
from pathlib import Path

# Layer-1.1 imports
_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[2]))

from src.chem.parametrize import ParametrizeResult, parametrize_smiles  # noqa: E402

logger = logging.getLogger(__name__)


@dataclass
class TrajectoryResult:
    """Envelope around a (possibly failed) MD run."""

    smiles: str
    ok: bool
    reason: str = ""
    formula: Optional[str] = None
    n_atoms: Optional[int] = None
    n_frames: int = 0
    dt_ps: Optional[float] = None
    total_ps: Optional[float] = None
    temperature_setpoint_K: Optional[float] = None
    platform: Optional[str] = None
    ff_version: Optional[str] = None
    charge_method: Optional[str] = None
    tool_versions: dict = field(default_factory=dict)
    # Per-frame scalars
    times_ps: list = field(default_factory=list)
    potential_energies_kJmol: list = field(default_factory=list)
    kinetic_energies_kJmol: list = field(default_factory=list)
    temperatures_K: list = field(default_factory=list)
    rmsd_A: list = field(default_factory=list)
    # Per-frame positions (Å), shape (n_frames, n_atoms, 3)
    positions_A: list = field(default_factory=list)
    # Topology for the viewer
    elements: list = field(default_factory=list)
    bonds: list = field(default_factory=list)

    def as_dict(self) -> dict:
        d = asdict(self)
        # Trim big arrays in summary dumps if caller asks.
        return d

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] {self.smiles}  {self.reason}"
        t_final = self.temperatures_K[-1] if self.temperatures_K else float("nan")
        e_final = (self.potential_energies_kJmol[-1]
                   if self.potential_energies_kJmol else float("nan"))
        rmsd_final = self.rmsd_A[-1] if self.rmsd_A else float("nan")
        return (f"[OK]   {self.smiles}  {self.formula}  "
                f"atoms={self.n_atoms}  frames={self.n_frames}  "
                f"total={self.total_ps:.2f} ps  T_final={t_final:.1f} K  "
                f"E_pot={e_final:.2f} kJ/mol  RMSD={rmsd_final:.3f} Å  "
                f"platform={self.platform}")


def _tool_versions() -> dict:
    versions = {}
    try:
        import openmm

        versions["openmm"] = openmm.__version__
    except Exception:
        pass
    try:
        import openff.toolkit

        versions["openff-toolkit"] = openff.toolkit.__version__
    except Exception:
        pass
    return versions


def _pick_platform(requested: Optional[str]) -> Any:
    """Return the best available OpenMM Platform.

    Preference order when no request: CUDA > OpenCL > CPU. Metal is
    selected only when explicitly requested because OpenMM's Metal
    backend is still experimental as of 8.x.
    """
    import openmm

    available = {openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())}
    if requested:
        if requested in available:
            return openmm.Platform.getPlatformByName(requested)
        raise RuntimeError(
            f"requested platform '{requested}' not available; "
            f"available: {sorted(available)}")
    for name in ("CUDA", "OpenCL", "CPU"):
        if name in available:
            return openmm.Platform.getPlatformByName(name)
    # Pick anything
    return openmm.Platform.getPlatform(0)


def _rmsd_all_atom(p0_A: "np.ndarray", p_A: "np.ndarray") -> float:
    """Simple no-alignment all-atom RMSD (centroid-subtracted).

    For a small ligand in vacuum this is a sufficient "did the
    structure blow up?" signal without pulling in an alignment
    dependency like MDAnalysis.
    """
    import numpy as np

    c0 = p0_A.mean(axis=0)
    c = p_A.mean(axis=0)
    d = (p_A - c) - (p0_A - c0)
    return float(np.sqrt((d * d).sum(axis=1).mean()))


def simulate_ligand(
    source: str | ParametrizeResult,
    *,
    n_steps: int = 5000,
    dt_fs: float = 2.0,
    save_every: int = 100,
    temperature_K: float = 300.0,
    friction_per_ps: float = 1.0,
    minimize: bool = True,
    platform: Optional[str] = None,
    random_seed: int = 1,
    ff_name: str = "openff-2.1.0.offxml",
    charge_method: str = "am1bcc",
) -> TrajectoryResult:
    """Run vacuum NVT Langevin MD on a small molecule.

    `source` is either a SMILES string (we parametrise internally) or
    a pre-computed `ParametrizeResult` (preferred — so Layer 1.1 and
    Layer 1.2 share one parametrisation per compound).

    Defaults: 5 000 steps × 2 fs = 10 ps, 100 frames. Small enough
    to finish in seconds for ligands on CPU, big enough to prove
    temperature equilibration.
    """
    import numpy as np

    versions = _tool_versions()

    if isinstance(source, str):
        param = parametrize_smiles(
            source, charge_method=charge_method, ff_name=ff_name,
            random_seed=random_seed)
    else:
        param = source

    if not param.ok:
        return TrajectoryResult(
            smiles=param.smiles, ok=False,
            reason=f"parametrise failed: {param.reason}",
            formula=param.formula, n_atoms=param.n_atoms,
            tool_versions=versions,
            elements=param.elements or [],
            bonds=param.bonds or [])

    # Rebuild OpenMM System from the OpenFF Molecule. We need an
    # OpenFF Molecule because ParametrizeResult intentionally only
    # stores primitive data; reparametrise here (fast — RDKit cached
    # the conformer already via the result's positions).
    try:
        from openff.toolkit import ForceField, Molecule
        import openmm
        from openmm import unit
        from openmm.app import Simulation
    except ImportError as e:
        return TrajectoryResult(
            smiles=param.smiles, ok=False,
            reason=f"openmm/openff import failed: {e}",
            tool_versions=versions,
            elements=param.elements or [],
            bonds=param.bonds or [])

    try:
        off_mol = Molecule.from_smiles(
            param.smiles, allow_undefined_stereo=True)
        off_mol.generate_conformers(n_conformers=1)
        off_mol.assign_partial_charges(partial_charge_method=charge_method)
        ff = ForceField(ff_name)
        system = ff.create_openmm_system(off_mol.to_topology())
        topology = off_mol.to_topology().to_openmm()
    except Exception as e:
        return TrajectoryResult(
            smiles=param.smiles, ok=False,
            reason=f"system build failed: {str(e).splitlines()[0]}",
            formula=param.formula, n_atoms=param.n_atoms,
            charge_method=charge_method, ff_version=ff_name,
            tool_versions=versions,
            elements=param.elements or [],
            bonds=param.bonds or [])

    # Initial positions (nm).
    pos_nm = np.array([[p[0], p[1], p[2]] for p in param.positions_nm])
    try:
        plat = _pick_platform(platform)
    except Exception as e:
        return TrajectoryResult(
            smiles=param.smiles, ok=False,
            reason=f"platform init failed: {e}",
            formula=param.formula, n_atoms=param.n_atoms,
            charge_method=charge_method, ff_version=ff_name,
            tool_versions=versions,
            elements=param.elements or [],
            bonds=param.bonds or [])

    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction_per_ps / unit.picosecond,
        dt_fs * unit.femtoseconds)
    integrator.setRandomNumberSeed(int(random_seed))

    sim = Simulation(topology, system, integrator, plat)
    sim.context.setPositions(pos_nm * unit.nanometer)

    try:
        if minimize:
            sim.minimizeEnergy(maxIterations=500)
        sim.context.setVelocitiesToTemperature(
            temperature_K * unit.kelvin, int(random_seed))
    except Exception as e:
        return TrajectoryResult(
            smiles=param.smiles, ok=False,
            reason=f"minimize / init failed: {str(e).splitlines()[0]}",
            formula=param.formula, n_atoms=param.n_atoms,
            charge_method=charge_method, ff_version=ff_name,
            platform=plat.getName(),
            tool_versions=versions,
            elements=param.elements or [],
            bonds=param.bonds or [])

    # Collect frame 0 state.
    times_ps: list[float] = []
    pes: list[float] = []
    kes: list[float] = []
    temps: list[float] = []
    rmsds: list[float] = []
    frames_A: list[np.ndarray] = []

    def _record(step: int):
        state = sim.context.getState(getEnergy=True, getPositions=True)
        pe_kJ = state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)
        ke_kJ = state.getKineticEnergy().value_in_unit(
            unit.kilojoule_per_mole)
        n_dof = max(1, 3 * system.getNumParticles() - system.getNumConstraints())
        t_K = (2.0 * ke_kJ * 1000.0
               / (n_dof * 8.314462618))  # ke*1000 J/mol, R = 8.314 J/(mol·K)
        pos_A = np.array(
            state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
        t_ps = step * dt_fs * 1e-3
        times_ps.append(t_ps)
        pes.append(float(pe_kJ))
        kes.append(float(ke_kJ))
        temps.append(float(t_K))
        frames_A.append(pos_A)
        rmsds.append(_rmsd_all_atom(frames_A[0], pos_A))

    _record(0)
    steps_run = 0
    while steps_run < n_steps:
        block = min(save_every, n_steps - steps_run)
        try:
            sim.step(block)
        except Exception as e:
            return TrajectoryResult(
                smiles=param.smiles, ok=False,
                reason=f"integration step failed at {steps_run}: "
                f"{str(e).splitlines()[0]}",
                formula=param.formula, n_atoms=param.n_atoms,
                charge_method=charge_method, ff_version=ff_name,
                platform=plat.getName(),
                dt_ps=dt_fs * 1e-3,
                total_ps=steps_run * dt_fs * 1e-3,
                tool_versions=versions,
                times_ps=times_ps, potential_energies_kJmol=pes,
                kinetic_energies_kJmol=kes, temperatures_K=temps,
                rmsd_A=rmsds, positions_A=[f.tolist() for f in frames_A],
                elements=[a.symbol for a in off_mol.atoms],
                bonds=[(b.atom1_index, b.atom2_index, float(b.bond_order))
                       for b in off_mol.bonds],
                n_frames=len(times_ps))
        steps_run += block
        _record(steps_run)

    # Gate checks
    reason = ""
    ok = True
    if any(math.isnan(v) or math.isinf(v) for v in pes + temps + rmsds):
        ok = False
        reason = "NaN / inf in telemetry"
    elif temps and abs(temps[-1] - temperature_K) > 50.0:
        ok = False
        reason = (f"final T {temps[-1]:.1f} K > 50 K from setpoint "
                  f"{temperature_K}")
    elif rmsds and rmsds[-1] > 10.0:
        ok = False
        reason = f"RMSD blow-up {rmsds[-1]:.2f} Å"

    return TrajectoryResult(
        smiles=param.smiles, ok=ok, reason=reason,
        formula=param.formula,
        n_atoms=off_mol.n_atoms,
        n_frames=len(times_ps),
        dt_ps=dt_fs * 1e-3,
        total_ps=n_steps * dt_fs * 1e-3,
        temperature_setpoint_K=temperature_K,
        platform=plat.getName(),
        ff_version=ff_name,
        charge_method=charge_method,
        tool_versions=versions,
        times_ps=times_ps,
        potential_energies_kJmol=pes,
        kinetic_energies_kJmol=kes,
        temperatures_K=temps,
        rmsd_A=rmsds,
        positions_A=[f.tolist() for f in frames_A],
        elements=[a.symbol for a in off_mol.atoms],
        bonds=[(b.atom1_index, b.atom2_index, float(b.bond_order))
               for b in off_mol.bonds],
    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("smiles")
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--dt", type=float, default=2.0, help="fs")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--platform", default=None)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    r = simulate_ligand(
        args.smiles,
        n_steps=args.steps, dt_fs=args.dt,
        temperature_K=args.temp,
        platform=args.platform, random_seed=args.seed)
    print(r.summary())
    sys.exit(0 if r.ok else 1)
