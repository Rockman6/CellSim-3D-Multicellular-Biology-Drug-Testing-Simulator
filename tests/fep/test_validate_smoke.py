"""fep-binding validate dry-run regression tests.

The validator is the biologist's last cheap check before committing
to a multi-hour sampled run. If it silently breaks — lets a typo'd
SMILES through, or misses a missing receptor — the user burns GPU
hours on a run doomed from parse time.

These tests run in < 1 s because the validator doesn't import
OpenFF / OpenMM.
"""

from __future__ import annotations

import io
import sys
import tempfile
import textwrap
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


def _run(yaml_text: str, extra_args=()):
    """Invoke the validate CLI on a temporary YAML and return the
    exit code. Captures stdout so the test harness stays quiet."""
    from src.fep.binding import main

    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False) as tmp:
        tmp.write(textwrap.dedent(yaml_text))
        tmp_path = tmp.name
    try:
        # Silence the validator's stdout by redirecting during the
        # call. Tests still inspect exit code and any trailing
        # state we expose.
        old = sys.stdout
        sys.stdout = io.StringIO()
        try:
            rc = main(["validate", tmp_path, *extra_args])
            captured = sys.stdout.getvalue()
        finally:
            sys.stdout = old
    finally:
        Path(tmp_path).unlink(missing_ok=True)
    return rc, captured


def test_valid_binding_yaml_passes():
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: biotin
        smiles: "O=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -18.3
    """
    rc, out = _run(yaml)
    assert rc == 0, f"valid YAML should pass; got {rc}, out:\n{out}"
    assert "PASS" in out
    assert "biotin" in out


def test_valid_hydration_yaml_passes_without_receptor():
    yaml = """
    entries:
      - name: methane
        smiles: "C"
        dG_hydration_kcalmol: 2.0
      - name: ethanol
        smiles: "CCO"
        dG_hydration_kcalmol: -5.01
    """
    rc, out = _run(yaml)
    assert rc == 0
    assert "hydration" in out
    assert "n/a (hydration YAML)" in out


def test_typod_smiles_fails():
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: bad_smiles
        smiles: "X(not a real smiles"
        dG_bind_kcalmol: -10.0
    """
    rc, out = _run(yaml)
    assert rc != 0, "typo'd SMILES should fail validation"
    assert "RDKit" in out or "parse" in out.lower()


def test_missing_receptor_fails_on_binding_yaml():
    yaml = """
    receptor:
      pdb_path: /path/that/does/not/exist.pdb

    entries:
      - name: biotin
        smiles: "O=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -18.3
    """
    rc, out = _run(yaml)
    assert rc != 0
    assert "MISSING" in out or "not found" in out.lower()


def test_undefined_stereo_bond_flagged():
    """Exocyclic imino (C=N) with unpinned E/Z should be flagged —
    the openff builder would otherwise pick arbitrarily.
    Precision: must NOT also flag amide C=O (biotin) or urea C=N
    ring-constrained geometry."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: iminobiotin_unpinned
        smiles: "N=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -10.8
      - name: biotin_clean
        smiles: "O=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -18.3
    """
    rc, out = _run(yaml)
    assert rc != 0, "unpinned imino should fail validation"
    assert "iminobiotin_unpinned" in out
    assert "undefined E/Z" in out
    # Precision: biotin must NOT be flagged (its C=O is terminal).
    lines = [line for line in out.splitlines()
             if "biotin_clean" in line and "undefined" in line]
    assert not lines, (
        f"biotin_clean was false-positive-flagged; "
        f"amide C=O is terminal, not E/Z-ambiguous: {lines}")


def test_pinned_stereo_passes():
    """Pinning E/Z with / \\ should let the SMILES through."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: iminobiotin_pinned
        smiles: "[H]/N=C1/N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -10.8
    """
    rc, out = _run(yaml)
    assert rc == 0, (
        f"pinned-stereo SMILES should pass; got rc={rc}\n{out}")
    assert "PASS" in out


def test_charged_ligand_flagged():
    """Formal-charge != 0 needs a Rocklin/Warren correction we
    haven't wired. Validator must flag so the biologist knows."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: acetate
        smiles: "CC(=O)[O-]"
        dG_bind_kcalmol: -5.0
    """
    rc, out = _run(yaml)
    assert rc != 0, "charged ligand should be flagged as an issue"
    assert "formal charge" in out.lower() or "rocklin" in out.lower()


def test_duplicate_names_flagged():
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: biotin
        smiles: "C"
        dG_bind_kcalmol: -10.0
      - name: biotin
        smiles: "CC"
        dG_bind_kcalmol: -11.0
    """
    rc, out = _run(yaml)
    assert rc != 0, "duplicate names should be flagged"
    assert "duplicate" in out.lower()


def test_empty_entries_list_fails():
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries: []
    """
    rc, out = _run(yaml)
    assert rc != 0
    assert "no entries" in out.lower()


def test_bundled_streptavidin_yaml_validates():
    """End-to-end on the committed binding benchmark — this is the
    CI gate that stops a SMILES typo from landing on main."""
    from src.fep.binding import main

    old = sys.stdout
    sys.stdout = io.StringIO()
    try:
        rc = main([
            "validate",
            "benchmarks/fep/binding_streptavidin.yaml"])
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old
    assert rc == 0, (
        f"bundled streptavidin YAML regressed: rc={rc}\n{out}")
    assert "PASS" in out
    assert "biotin" in out and "desthiobiotin" in out


def test_bundled_freesolv_yaml_validates():
    from src.fep.binding import main

    old = sys.stdout
    sys.stdout = io.StringIO()
    try:
        rc = main([
            "validate",
            "benchmarks/fep/freesolv_12.yaml"])
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old
    assert rc == 0
    assert "PASS" in out


if __name__ == "__main__":
    funcs = [
        test_valid_binding_yaml_passes,
        test_valid_hydration_yaml_passes_without_receptor,
        test_typod_smiles_fails,
        test_missing_receptor_fails_on_binding_yaml,
        test_undefined_stereo_bond_flagged,
        test_pinned_stereo_passes,
        test_charged_ligand_flagged,
        test_duplicate_names_flagged,
        test_empty_entries_list_fails,
        test_bundled_streptavidin_yaml_validates,
        test_bundled_freesolv_yaml_validates,
    ]
    fails = []
    for f in funcs:
        try:
            f()
            print(f"[PASS] {f.__name__}")
        except AssertionError as e:
            print(f"[FAIL] {f.__name__}: {e}")
            fails.append(f.__name__)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"[ERROR] {f.__name__}: {e}")
            fails.append(f.__name__)
    print(f"{len(funcs) - len(fails)}/{len(funcs)} PASS")
    sys.exit(0 if not fails else 1)
