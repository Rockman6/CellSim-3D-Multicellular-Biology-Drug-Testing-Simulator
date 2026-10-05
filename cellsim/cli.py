"""`cellsim` console entry point.

Python twin of `scripts/cellsim` (the bash dispatcher), so the command
exists after `pip install -e .` without a repository checkout on PATH.
Module-backed subcommands run via `runpy`; the handful that live under
`scripts/` are located relative to the repository root and run the same
way, which means they work from a checkout and report clearly when the
package was installed without one.
"""
from __future__ import annotations

import importlib.metadata
import runpy
import sys
import warnings
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# subcommand -> module run as __main__
MODULES: dict[str, str] = {
    "dock": "cellsim.dock.batch",
    "dock-one": "cellsim.dock.viewer",
    "off-target": "cellsim.dock.off_target",
    "triage-png": "cellsim.dock.triage_viewer",
    "shortlist": "cellsim.dock.shortlist",
    "to-md": "cellsim.dock.to_md",
    "redock": "cellsim.dock.vina",
    "pocket": "cellsim.dock.pocket_detect",
    "cyp-inhibit": "cellsim.dock.cyp_inhibition",
    "profile": "cellsim.chem.profile",
    "admet": "cellsim.chem.admet",
    "parametrize": "cellsim.chem.parametrize",
    "xtb": "cellsim.quantum.xtb",
    "dft": "cellsim.quantum.dft",
    "som": "cellsim.quantum.metabolism",
    "som-validate": "cellsim.quantum.som_validation",
    "fep-hyd": "cellsim.fep",
    "fep-freesolv": "cellsim.fep.freesolv_validate",
    "fep-report": "cellsim.fep.report",
    "fep-binding": "cellsim.fep.binding",
    "md-ligand": "cellsim.md.simulate",
    "md-protein": "cellsim.md.protein",
    "md-view": "cellsim.md.protein_viewer",
    "uq-mc": "cellsim.uq.dock_mc",
    "uq-sobol": "cellsim.uq.sobol",
    "cell-sim": "cellsim.cell.engine",
    "cell-stream": "cellsim.cell.stream",
    "dish": "cellsim.cell.dish",
    "plate": "cellsim.plate",
}

# subcommand -> script under <repo>/scripts
SCRIPTS: dict[str, str] = {
    "doctor": "doctor.py",
    "bench-all": "bench_all_yamls.py",
    "calib": "calib_cli.py",
    "bench": "bench_cli.py",
    "fetch-pdb": "fetch_pdb.py",
    "fetch-chembl": "fetch_chembl_sample.py",
    "validate-gdsc": "validate_gdsc.py",
}

USAGE = """usage: cellsim <subcommand> [args...]

  Cell drug-response engine (Phase 1)
    cell-sim      dose-response of one cell line to one drug -> viability + IC50
    cell-stream   one population under a dosing schedule -> per-tick JSON Lines for a UI
    dish          cells in space (monolayer or spheroid) -> per-cell JSON Lines for a UI
    plate         your plate-reader export -> viability, IC50 with a CI, GR metrics
    validate-gdsc fit one potency per drug on p53-wt lines, predict the rest vs GDSC

  Screening and triage
    dock          batch docking screen -> ranked CSV ± MC error bars ± profile PNGs
    dock-one      single-compound dock + viewer
    off-target    one compound vs N receptors (selectivity)
    triage-png    triage-breakdown dashboard from a batch CSV
    shortlist     keep only follow_up + review rows of a batch CSV
    to-md         batch CSV -> markdown table
    redock        single-pair re-dock with crystal RMSD
    pocket        fpocket binding-site detection
    cyp-inhibit   CYP3A4 inhibition / DDI-risk screen

  Compound properties
    profile       one-page drug profile dashboard
    admet         Lipinski / QED / logS / TPSA descriptors
    parametrize   SMILES -> parametrised OpenMM system
    xtb           xTB GFN2 single-point + HOMO/LUMO + charges
    dft           PySCF B3LYP/def2-SVP single-point
    som           CYP3A4 site-of-metabolism (advisory)
    som-validate  SoM predictor vs literature sites

  Experimental free energy (not a validated path; see docs/VALIDATION.md)
    fep-hyd fep-freesolv fep-report fep-binding

  MD and uncertainty
    md-ligand md-protein md-view uq-mc uq-sobol

  Repository tools (need a checkout)
    doctor        install sanity check
    bench-all     validate every benchmarks/fep/*.yaml
    calib         run one calibration-set YAML + scatter PNG
    bench         run all bundled benchmarks
    fetch-pdb     download a PDB by ID
    fetch-chembl  drug-like SMILES subset from ChEMBL

  help / version

Every subcommand accepts --help for its own arguments.
"""

_VERSION_DISTS = [
    ("rdkit", "rdkit"), ("openmm", "openmm"), ("openmmtools", "openmmtools"),
    ("openmmforcefields", "openmmforcefields"), ("openff-toolkit", "openff-toolkit"),
    ("openff-interchange", "openff-interchange"), ("pdbfixer", "pdbfixer"),
    ("pymbar", "pymbar"), ("vina", "vina"), ("meeko", "meeko"),
    ("posebusters", "posebusters"), ("xtb", "xtb-python"), ("pyscf", "pyscf"),
    ("salib", "SALib"),
]


def _version() -> int:
    try:
        own = importlib.metadata.version("cellsim")
    except importlib.metadata.PackageNotFoundError:
        own = "(not installed; running from checkout)"
    print(f"cellsim   {own}")
    print(f"python    {sys.version.split()[0]}")
    for label, dist in _VERSION_DISTS:
        try:
            print(f"{label:<20s}{importlib.metadata.version(dist)}")
        except importlib.metadata.PackageNotFoundError:
            print(f"{label:<20s}MISSING")
    return 0


def _run_module(module: str, argv: list[str]) -> int:
    sys.argv = [f"cellsim {argv[0]}" if argv else "cellsim", *argv[1:]]
    # Several packages import their submodules in __init__ (e.g.
    # cellsim.chem imports admet), so runpy warns that the module is
    # already in sys.modules. Harmless here; the module is re-executed
    # as __main__ exactly as `python -m` would.
    warnings.filterwarnings(
        "ignore", message=".*found in sys.modules after import of package.*",
        category=RuntimeWarning)
    try:
        runpy.run_module(module, run_name="__main__", alter_sys=True)
    except SystemExit as e:  # propagate the module's exit status
        return int(e.code or 0) if not isinstance(e.code, str) else 1
    return 0


def _run_script(script: str, argv: list[str]) -> int:
    path = REPO_ROOT / "scripts" / script
    if not path.exists():
        print(
            f"cellsim: '{argv[0]}' needs a repository checkout "
            f"(expected {path}). Clone the repo and run from inside it.",
            file=sys.stderr,
        )
        return 2
    sys.argv = [str(path), *argv[1:]]
    try:
        runpy.run_path(str(path), run_name="__main__")
    except SystemExit as e:
        return int(e.code or 0) if not isinstance(e.code, str) else 1
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(USAGE)
        return 0
    cmd = argv[0]
    if cmd in ("version", "--version", "-V"):
        return _version()
    if cmd in MODULES:
        return _run_module(MODULES[cmd], argv)
    if cmd in SCRIPTS:
        return _run_script(SCRIPTS[cmd], argv)
    print(f"cellsim: unknown subcommand '{cmd}'. Run 'cellsim help' for a list.",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
