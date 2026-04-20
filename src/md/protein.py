"""Protein loader: PDB -> minimised, solvated OpenMM system.

Pipeline:
    PDB -> PDBFixer (missing residues / atoms / terminals; protonate
           at pH 7.0) -> OpenMM Modeller (solvate in TIP3P cubic box
           with NaCl 150 mM) -> AMBER14 / TIP3P-FB force fields
           -> openmm.System
    then energy minimisation in the returned `_initial_frame`.

Gate-quality contract:
    - reason is set when anything fails; `ok` is False in that case.
    - `_initial_frame` holds a minimised frame ready for a Layer-1.2
      production simulation.
    - Structure + atom count + box info are stored on the result so
      downstream gates can assert on them without re-loading.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional
import logging
import sys

logger = logging.getLogger(__name__)


@dataclass
class ProteinSystemResult:
    """Envelope for a loaded + minimised protein system."""

    source: str
    ok: bool
    reason: str = ""
    pdb_id: Optional[str] = None
    n_residues_protein: Optional[int] = None
    n_atoms_protein: Optional[int] = None
    n_atoms_total: Optional[int] = None
    n_waters: Optional[int] = None
    n_ions_na: int = 0
    n_ions_cl: int = 0
    box_nm: Optional[tuple] = None
    ff_protein: Optional[str] = None
    ff_water: Optional[str] = None
    ph: Optional[float] = None
    padding_nm: Optional[float] = None
    platform: Optional[str] = None
    e_initial_kJmol: Optional[float] = None
    e_minimised_kJmol: Optional[float] = None
    minimised_positions_nm: list = field(default_factory=list)
    tool_versions: dict = field(default_factory=dict)
    # Private handles (not serialised):
    _system: Any = None
    _topology: Any = None
    _modeller: Any = None

    def as_dict(self) -> dict:
        d = asdict(self)
        for k in ("_system", "_topology", "_modeller"):
            d.pop(k, None)
        return d

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] {self.source}  {self.reason}"
        return (f"[OK]   {self.source}  ff={self.ff_protein}  "
                f"protein_atoms={self.n_atoms_protein}  "
                f"total_atoms={self.n_atoms_total}  "
                f"waters={self.n_waters}  "
                f"Na+={self.n_ions_na}  Cl-={self.n_ions_cl}  "
                f"box={self.box_nm}  "
                f"E_min={self.e_minimised_kJmol:.1f} kJ/mol  "
                f"dE={self.e_initial_kJmol - self.e_minimised_kJmol:.1f} "
                f"kJ/mol")


def _tool_versions() -> dict:
    versions: dict = {}
    for mod_name, attr in (("openmm", "__version__"),
                           ("pdbfixer", "__version__"),
                           ("openff.toolkit", "__version__")):
        try:
            mod = __import__(mod_name)
            if "." in mod_name:
                for part in mod_name.split(".")[1:]:
                    mod = getattr(mod, part)
            versions[mod_name] = getattr(mod, attr, "?")
        except Exception:
            pass
    return versions


def _pick_platform(requested: Optional[str]) -> Any:
    import openmm

    available = {openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())}
    if requested:
        if requested in available:
            return openmm.Platform.getPlatformByName(requested)
        raise RuntimeError(
            f"requested platform '{requested}' not available; "
            f"have {sorted(available)}")
    for name in ("CUDA", "OpenCL", "CPU"):
        if name in available:
            return openmm.Platform.getPlatformByName(name)
    return openmm.Platform.getPlatform(0)


def load_protein_pdb(
    pdb_path: str | Path,
    *,
    ff_protein: str = "amber14-all.xml",
    ff_water: str = "amber14/tip3pfb.xml",
    ph: float = 7.0,
    padding_nm: float = 1.0,
    ionic_strength_M: float = 0.15,
    minimise_max_iter: int = 2000,
    platform: Optional[str] = None,
) -> ProteinSystemResult:
    """Load a PDB, fix it, solvate it, minimise it."""
    versions = _tool_versions()
    pdb_path = Path(pdb_path)
    source = str(pdb_path)

    try:
        from pdbfixer import PDBFixer
        import openmm
        from openmm import unit
        from openmm.app import ForceField, Modeller, PDBFile, Simulation
    except ImportError as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"openmm/pdbfixer import failed: {e}",
            tool_versions=versions)

    if not pdb_path.exists():
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"file not found: {pdb_path}",
            tool_versions=versions)

    try:
        fixer = PDBFixer(filename=str(pdb_path))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=ph)
    except Exception as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"PDBFixer failed: {str(e).splitlines()[0]}",
            tool_versions=versions)

    # Snapshot protein-only atom / residue counts before solvation.
    n_res_protein = len(list(fixer.topology.residues()))
    n_atoms_protein = fixer.topology.getNumAtoms()

    try:
        ff = ForceField(ff_protein, ff_water)
    except Exception as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"ForceField load failed: {str(e).splitlines()[0]}",
            ff_protein=ff_protein, ff_water=ff_water,
            tool_versions=versions)

    modeller = Modeller(fixer.topology, fixer.positions)
    try:
        modeller.addSolvent(
            ff,
            padding=padding_nm * unit.nanometer,
            ionicStrength=ionic_strength_M * unit.molar,
            model="tip3p")
    except Exception as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"addSolvent failed: {str(e).splitlines()[0]}",
            ff_protein=ff_protein, ff_water=ff_water,
            n_residues_protein=n_res_protein,
            n_atoms_protein=n_atoms_protein,
            tool_versions=versions)

    try:
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=openmm.app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=openmm.app.HBonds,
            rigidWater=True)
    except Exception as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"createSystem failed: {str(e).splitlines()[0]}",
            ff_protein=ff_protein, ff_water=ff_water,
            n_residues_protein=n_res_protein,
            n_atoms_protein=n_atoms_protein,
            tool_versions=versions)

    try:
        plat = _pick_platform(platform)
    except Exception as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"platform init failed: {e}",
            ff_protein=ff_protein, ff_water=ff_water,
            tool_versions=versions)

    # Minimise. Use a dummy Langevin integrator because OpenMM
    # requires an integrator to build a Simulation, but we only
    # actually call minimizeEnergy + getState here.
    integrator = openmm.LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtoseconds)
    sim = Simulation(modeller.topology, system, integrator, plat)
    sim.context.setPositions(modeller.positions)

    try:
        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
        e0_kJ = float(e0.value_in_unit(unit.kilojoule_per_mole))
        sim.minimizeEnergy(maxIterations=minimise_max_iter)
        state = sim.context.getState(getEnergy=True, getPositions=True)
        e1_kJ = float(state.getPotentialEnergy()
                      .value_in_unit(unit.kilojoule_per_mole))
    except Exception as e:
        return ProteinSystemResult(
            source=source, ok=False,
            reason=f"minimize failed: {str(e).splitlines()[0]}",
            ff_protein=ff_protein, ff_water=ff_water,
            n_residues_protein=n_res_protein,
            n_atoms_protein=n_atoms_protein,
            platform=plat.getName(),
            tool_versions=versions)

    # Count species after solvation.
    residues = list(modeller.topology.residues())
    n_waters = sum(1 for r in residues if r.name in ("HOH", "WAT", "TIP3"))
    n_na = sum(1 for r in residues if r.name in ("NA", "NA+"))
    n_cl = sum(1 for r in residues if r.name in ("CL", "CL-"))

    box_vecs = modeller.topology.getPeriodicBoxVectors()
    if box_vecs is not None:
        box = tuple(float(box_vecs[i][i].value_in_unit(unit.nanometer))
                    for i in range(3))
    else:
        box = None

    pos_nm = [list(p.value_in_unit(unit.nanometer))
              for p in state.getPositions(asNumpy=False)]

    return ProteinSystemResult(
        source=source, ok=True,
        pdb_id=pdb_path.stem.upper(),
        n_residues_protein=n_res_protein,
        n_atoms_protein=n_atoms_protein,
        n_atoms_total=modeller.topology.getNumAtoms(),
        n_waters=n_waters,
        n_ions_na=n_na, n_ions_cl=n_cl,
        box_nm=box,
        ff_protein=ff_protein, ff_water=ff_water,
        ph=ph, padding_nm=padding_nm,
        platform=plat.getName(),
        e_initial_kJmol=e0_kJ, e_minimised_kJmol=e1_kJ,
        minimised_positions_nm=pos_nm,
        tool_versions=versions,
        _system=system, _topology=modeller.topology,
        _modeller=modeller,
    )


def short_protein_md(
    result: ProteinSystemResult,
    *,
    n_steps: int = 1000,
    dt_fs: float = 2.0,
    save_every: int = 100,
    temperature_K: float = 300.0,
    friction_per_ps: float = 1.0,
    random_seed: int = 1,
) -> "ProteinTrajectoryResult":
    """Short equilibration run on a loaded ProteinSystemResult.

    Reports backbone-Cα RMSD vs the minimised start. Intended for
    spot-checks, not the full 100 ns ubiquitin gate (that lives on a
    GPU pipeline).
    """
    import math

    import numpy as np
    import openmm
    from openmm import unit
    from openmm.app import Simulation

    if not result.ok:
        return ProteinTrajectoryResult(
            source=result.source, ok=False,
            reason=f"cannot run MD: protein load failed ({result.reason})")

    platform = openmm.Platform.getPlatformByName(result.platform)
    integrator = openmm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin,
        friction_per_ps / unit.picosecond,
        dt_fs * unit.femtoseconds)
    integrator.setRandomNumberSeed(int(random_seed))
    sim = Simulation(result._topology, result._system, integrator, platform)
    sim.context.setPositions(
        np.array(result.minimised_positions_nm) * unit.nanometer)
    sim.context.setVelocitiesToTemperature(
        temperature_K * unit.kelvin, int(random_seed))

    # Identify Cα indices once.
    ca_indices = [a.index for a in result._topology.atoms() if a.name == "CA"]
    all_indices = list(range(result.n_atoms_total or
                             sum(1 for _ in result._topology.atoms())))

    def _record(step: int):
        state = sim.context.getState(getEnergy=True, getPositions=True)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        n_dof = max(1, 3 * result._system.getNumParticles()
                    - result._system.getNumConstraints())
        t_K = (2.0 * ke * 1000.0 / (n_dof * 8.314462618))
        pos_A = np.array(state.getPositions(asNumpy=True)
                         .value_in_unit(unit.angstrom))
        return step * dt_fs * 1e-3, pe, ke, t_K, pos_A

    times_ps: list = []
    pes: list = []
    kes: list = []
    temps: list = []
    rmsds_ca: list = []
    rmsds_all: list = []

    t0, pe0, ke0, T0, pos0 = _record(0)
    pos0_ca = pos0[ca_indices]
    times_ps.append(t0)
    pes.append(pe0); kes.append(ke0); temps.append(T0)
    rmsds_ca.append(0.0); rmsds_all.append(0.0)

    steps_run = 0
    while steps_run < n_steps:
        block = min(save_every, n_steps - steps_run)
        sim.step(block)
        steps_run += block
        t, pe, ke, T, pos = _record(steps_run)
        times_ps.append(t); pes.append(pe); kes.append(ke); temps.append(T)
        # Cα RMSD (no alignment — just centroid-subtracted)
        pc = pos[ca_indices]
        c0 = pos0_ca.mean(axis=0); c = pc.mean(axis=0)
        d = (pc - c) - (pos0_ca - c0)
        rmsds_ca.append(float(np.sqrt((d * d).sum(axis=1).mean())))
        # All-atom RMSD (centroid-subtracted; coarse signal)
        c0a = pos0.mean(axis=0); ca = pos.mean(axis=0)
        da = (pos - ca) - (pos0 - c0a)
        rmsds_all.append(float(np.sqrt((da * da).sum(axis=1).mean())))

    ok = True
    reason = ""
    if any(math.isnan(v) or math.isinf(v) for v in pes + temps + rmsds_ca):
        ok = False
        reason = "NaN / inf in telemetry"
    elif temps and abs(temps[-1] - temperature_K) > 50.0:
        ok = False
        reason = (f"final T {temps[-1]:.1f} K > 50 K from setpoint "
                  f"{temperature_K}")
    elif rmsds_ca and rmsds_ca[-1] > 5.0:
        ok = False
        reason = f"Cα RMSD blow-up {rmsds_ca[-1]:.2f} Å"

    return ProteinTrajectoryResult(
        source=result.source, ok=ok, reason=reason,
        pdb_id=result.pdb_id,
        n_atoms=result.n_atoms_total,
        n_residues_protein=result.n_residues_protein,
        n_frames=len(times_ps),
        dt_ps=dt_fs * 1e-3,
        total_ps=n_steps * dt_fs * 1e-3,
        temperature_setpoint_K=temperature_K,
        platform=result.platform,
        times_ps=times_ps,
        potential_energies_kJmol=pes,
        kinetic_energies_kJmol=kes,
        temperatures_K=temps,
        rmsd_ca_A=rmsds_ca,
        rmsd_all_A=rmsds_all,
    )


@dataclass
class ProteinTrajectoryResult:
    source: str
    ok: bool
    reason: str = ""
    pdb_id: Optional[str] = None
    n_atoms: Optional[int] = None
    n_residues_protein: Optional[int] = None
    n_frames: int = 0
    dt_ps: Optional[float] = None
    total_ps: Optional[float] = None
    temperature_setpoint_K: Optional[float] = None
    platform: Optional[str] = None
    times_ps: list = field(default_factory=list)
    potential_energies_kJmol: list = field(default_factory=list)
    kinetic_energies_kJmol: list = field(default_factory=list)
    temperatures_K: list = field(default_factory=list)
    rmsd_ca_A: list = field(default_factory=list)
    rmsd_all_A: list = field(default_factory=list)

    def summary(self) -> str:
        if not self.ok:
            return f"[FAIL] {self.pdb_id or self.source}  {self.reason}"
        tf = self.temperatures_K[-1] if self.temperatures_K else float("nan")
        rf = self.rmsd_ca_A[-1] if self.rmsd_ca_A else float("nan")
        pef = (self.potential_energies_kJmol[-1]
               if self.potential_energies_kJmol else float("nan"))
        return (f"[OK]   {self.pdb_id or self.source}  atoms={self.n_atoms}  "
                f"frames={self.n_frames}  total={self.total_ps:.2f} ps  "
                f"T_final={tf:.1f} K  RMSD_Cα={rf:.3f} Å  "
                f"E_pot={pef:.1f} kJ/mol  platform={self.platform}")


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb")
    ap.add_argument("--md-steps", type=int, default=0,
                    help="also run N MD steps after minimisation")
    ap.add_argument("--padding", type=float, default=1.0)
    ap.add_argument("--platform", default=None)
    args = ap.parse_args()

    r = load_protein_pdb(args.pdb, padding_nm=args.padding,
                         platform=args.platform)
    print(r.summary())
    if r.ok and args.md_steps > 0:
        tr = short_protein_md(r, n_steps=args.md_steps)
        print(tr.summary())
    sys.exit(0 if r.ok else 1)
