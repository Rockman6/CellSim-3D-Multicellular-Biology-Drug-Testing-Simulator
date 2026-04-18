#!/usr/bin/env python3
"""
fold_protein.py — Physics-based protein folding using OpenMM

Implements REAL molecular dynamics:
  - Builds peptide chain from amino acid sequence (via PDBFixer)
  - AMBER ff14SB force field (bonded: bonds, angles, torsions; nonbonded: LJ, Coulomb)
  - Implicit solvent (GBn2 model for hydrophobic effect)
  - Langevin dynamics at 300K
  - Energy minimization + MD folding

This is the physics described in Alberts Ch 3:
  - Peptide bonds (C-N, planar)
  - Hydrogen bonds → α-helix, β-sheet
  - Hydrophobic effect → core packing
  - Van der Waals → close packing
  - Electrostatics → salt bridges

Usage:
    python fold_protein.py <sequence> <output.pdb> [--steps N]
"""

import sys, os, time, tempfile
from io import StringIO

def sequence_to_pdb(sequence):
    """Build an extended peptide chain PDB from 1-letter amino acid sequence."""
    AA3 = {
        'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','E':'GLU','Q':'GLN',
        'G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS','M':'MET','F':'PHE',
        'P':'PRO','S':'SER','T':'THR','W':'TRP','Y':'TYR','V':'VAL'
    }
    import math
    lines = []
    a = 1
    for i, aa in enumerate(sequence):
        r = AA3.get(aa.upper(), 'ALA')
        # Extended chain: phi=-180, psi=180 → 3.8Å between Cα atoms
        # Place backbone N-Cα-C along a line with proper geometry
        x_base = i * 3.8
        # N
        lines.append(f"ATOM  {a:5d}  N   {r} A{i+1:4d}    {x_base:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           N")
        a += 1
        # CA
        lines.append(f"ATOM  {a:5d}  CA  {r} A{i+1:4d}    {x_base+1.458:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
        a += 1
        # C
        lines.append(f"ATOM  {a:5d}  C   {r} A{i+1:4d}    {x_base+2.400:8.3f}{1.200:8.3f}{0.0:8.3f}  1.00  0.00           C")
        a += 1
        # O
        lines.append(f"ATOM  {a:5d}  O   {r} A{i+1:4d}    {x_base+2.400:8.3f}{2.400:8.3f}{0.0:8.3f}  1.00  0.00           O")
        a += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def fold(sequence, output_path, num_steps=10000):
    """Fold a protein using OpenMM molecular dynamics with AMBER force field."""
    import openmm
    from openmm import app, unit, LangevinMiddleIntegrator, Platform
    from pdbfixer import PDBFixer

    n = len(sequence)
    print(f"[fold] Sequence ({n} aa): {sequence[:40]}{'...' if n > 40 else ''}")

    # Step 1: Build extended chain PDB
    print("[fold] Building extended peptide chain...")
    pdb_str = sequence_to_pdb(sequence)
    tmp = tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False)
    tmp.write(pdb_str)
    tmp.close()

    # Step 2: Use PDBFixer to add missing atoms (side chains, hydrogens)
    print("[fold] Adding side chains and hydrogens with PDBFixer...")
    try:
        fixer = PDBFixer(filename=tmp.name)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)  # pH 7.0

        # Save fixed PDB
        fixed_path = tmp.name + ".fixed.pdb"
        with open(fixed_path, 'w') as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
        print(f"[fold] Fixed structure: {fixer.topology.getNumAtoms()} atoms")
    except Exception as e:
        print(f"[fold] PDBFixer failed ({e}), using minimal backbone")
        fixed_path = tmp.name
        fixer = None

    # Step 3: Set up OpenMM simulation with AMBER force field
    print("[fold] Setting up AMBER ff14SB + implicit solvent (GBn2)...")
    try:
        pdb = app.PDBFile(fixed_path)
        forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            hydrogenMass=1.5 * unit.amu  # Hydrogen mass repartitioning for larger timestep
        )
    except Exception as e:
        print(f"[fold] Force field setup failed ({e})")
        print("[fold] Saving unfolded structure...")
        os.rename(tmp.name, output_path)
        return False

    # Step 4: Langevin integrator (thermostat at 300K)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,        # Temperature
        1.0 / unit.picosecond,    # Friction coefficient
        4.0 * unit.femtoseconds   # Timestep (4fs with HMR)
    )

    # Use best available platform
    platform_name = 'CPU'
    for name in ['OpenCL', 'CPU']:
        try:
            Platform.getPlatformByName(name)
            platform_name = name
            break
        except:
            continue
    platform = Platform.getPlatformByName(platform_name)
    print(f"[fold] Platform: {platform_name}")

    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    # Step 5: Energy minimization (remove steric clashes)
    print("[fold] Energy minimization...")
    t0 = time.time()
    simulation.minimizeEnergy(maxIterations=2000)
    state = simulation.context.getState(getEnergy=True)
    e_min = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"[fold] Minimized energy: {e_min:.1f} kJ/mol ({time.time()-t0:.1f}s)")

    # Step 6: MD folding simulation
    print(f"[fold] Running {num_steps} steps of MD at 300K...")
    t0 = time.time()

    # Report every 10% of steps
    report_interval = max(num_steps // 10, 100)
    for step in range(0, num_steps, report_interval):
        simulation.step(min(report_interval, num_steps - step))
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
        elapsed = time.time() - t0
        pct = min(100, (step + report_interval) * 100 // num_steps)
        print(f"  [{pct:3d}%] PE={pe:10.1f} KE={ke:8.1f} kJ/mol  ({elapsed:.1f}s)")

    # Step 7: Final energy minimization to settle
    print("[fold] Final minimization...")
    simulation.minimizeEnergy(maxIterations=1000)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    e_final = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"[fold] Final energy: {e_final:.1f} kJ/mol")

    # Step 8: Save folded structure as PDB
    positions = state.getPositions()
    with open(output_path, 'w') as f:
        app.PDBFile.writeFile(pdb.topology, positions, f)

    # Count atoms
    n_atoms = pdb.topology.getNumAtoms()
    print(f"[fold] Saved {output_path}: {n_atoms} atoms, {n} residues")
    print(f"[fold] Total time: {time.time()-t0:.1f}s")

    # Cleanup
    os.unlink(tmp.name)
    if os.path.exists(fixed_path):
        os.unlink(fixed_path)

    return True


def main():
    if len(sys.argv) < 3:
        print("Usage: fold_protein.py <sequence> <output.pdb> [--steps N]")
        print("Example: fold_protein.py MVHLTPEEKSAVTALWGKVNV folded.pdb --steps 50000")
        sys.exit(1)

    sequence = sys.argv[1]
    output_path = sys.argv[2]
    num_steps = 10000
    for i, arg in enumerate(sys.argv):
        if arg == '--steps' and i + 1 < len(sys.argv):
            num_steps = int(sys.argv[i + 1])

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    fold(sequence, output_path, num_steps)


if __name__ == "__main__":
    main()
