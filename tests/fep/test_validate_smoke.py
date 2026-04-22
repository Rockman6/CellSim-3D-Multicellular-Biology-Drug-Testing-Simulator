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


def test_warning_does_not_fail_exit():
    """Rotatable-bond > 10 should emit a warning in stdout AND
    leave exit code 0 — the biologist gets informed, not blocked.
    Uses a deliberately flexible aliphatic chain (14 rot bonds);
    lapatinib (11 rot) is the real-world motivator.
    """
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: c20_chain
        smiles: "CCCCCCCCCCCCCCCCCCCC"
        dG_bind_kcalmol: -8.0
    """
    rc, out = _run(yaml)
    assert rc == 0, (
        f"advisory warning must not fail exit; got rc={rc}\n{out}")
    assert "warning" in out.lower()
    assert "rotatable bonds > 10" in out
    assert ("PASS (with warnings)" in out
            or "PASS — ready" in out), (
        f"verdict line should still say PASS; got:\n{out[-300:]}")


def test_error_does_fail_exit():
    """A real error (unpinned stereo) must continue to hard-FAIL
    even when mixed with warnings on other entries."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: flexible
        smiles: "CCCCCCCCCCCCCCCCCCCC"
        dG_bind_kcalmol: -8.0
      - name: unpinned_imino
        smiles: "N=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -10.0
    """
    rc, out = _run(yaml)
    assert rc != 0, "any error should fail, even with warnings"
    # Both printed: warning AND error sections.
    assert "warning" in out.lower()
    assert "issue(s)" in out
    assert "FAIL" in out


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
    """Formal-charge != 0 carries a PBC self-interaction error in
    absolute ΔG (Rocklin 2013). Validator emits a WARNING (not
    a hard error) so biologists studying ionic series aren't
    blocked — the error cancels in same-charge ΔΔG.
    """
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: acetate
        smiles: "CC(=O)[O-]"
        dG_bind_kcalmol: -5.0
    """
    rc, out = _run(yaml)
    # Warning, not hard fail — proceed.
    assert rc == 0, (
        f"charged ligand should warn but pass; got rc={rc}\n{out}")
    assert "formal charge" in out.lower()
    assert "rocklin" in out.lower()
    assert "PBC self-interaction error" in out
    assert "ΔΔG within a same-charge" in out, (
        "biologist actionability: must mention ΔΔG-cancels caveat")


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


def test_binding_site_center_unit_mismatch_warns():
    """Biologist copy-pastes a coordinate from the docking YAML
    (Ångströms) into a binding FEP YAML (nm). Protein bounding
    boxes are 3-10 nm — any axis > 10 nm is almost certainly a
    factor-of-10 unit error. Validator warns (not FAIL) with
    the suggested nm conversion."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1m17.pdb
      binding_site_center_nm: [22.01, 0.25, 52.79]

    entries:
      - name: erlotinib
        smiles: "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"
        dG_bind_kcalmol: -11.86
    """
    rc, out = _run(yaml)
    assert rc == 0, (
        "unit warning must not fail exit (advisory only); "
        f"got rc={rc}\n{out}")
    assert "likely Å, not nm" in out
    assert "[2.201, 0.025, 5.279]" in out, (
        f"suggested nm conversion should be printed; got:\n{out}")


def test_yaml_entry_key_typo_flagged():
    """Typos in per-entry keys (dG_bond_kcalmol instead of
    dG_bind_kcalmol) should warn with a specific 'did you mean X'
    suggestion. The bench driver silently ignores unknown keys so
    undetected typos would produce a run with missing data.

    Two scenarios:
      A. Typo on dG_bind_kcalmol alongside valid smiles + missing
         ΔG → warning + 'missing expt ΔG' error (rc != 0).
      B. Only an extra custom-looking key (no load-bearing typo) →
         warning + clean pass."""
    # Scenario A: typo on the expt ΔG key.
    yaml_a = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: biotin
        smiles: "O=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bond_kcalmol: -18.3
    """
    rc, out = _run(yaml_a)
    # The typo'd key would silently drop expt ΔG — validator
    # detects 'missing expt ΔG' as a hard error AND flags the
    # typo.
    assert rc != 0
    assert "unknown entry key 'dG_bond_kcalmol'" in out
    assert "did you mean 'dG_bind_kcalmol'" in out

    # Scenario B: author extends schema with a custom field — no
    # load-bearing typo, just a warning.
    yaml_b = """
    receptor:
      pdb_path: benchmarks/dock/1stp.pdb

    entries:
      - name: biotin
        smiles: "O=C1N[C@@H]2[C@H](SC[C@@H]2CCCCC(=O)O)N1"
        dG_bind_kcalmol: -18.3
        assay_id: "Green1975"
    """
    rc, out = _run(yaml_b)
    assert rc == 0, (
        f"custom-looking extra key should warn but not fail; "
        f"got rc={rc}\n{out}")
    assert "unknown entry key 'assay_id'" in out


def test_binding_site_center_outside_protein_bbox_warns():
    """Correct nm units but coordinate far outside the protein —
    possible wrong-receptor copy-paste. Validator warns (not FAIL)
    with the bounding-box numbers so biologist can sanity-check."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1m17.pdb
      binding_site_center_nm: [8.0, 8.0, 8.0]

    entries:
      - name: erlotinib
        smiles: "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"
        dG_bind_kcalmol: -11.86
    """
    rc, out = _run(yaml)
    assert rc == 0, "off-box warning must not fail exit"
    assert "outside the protein's atom-coord bounding box" in out
    # Must mention which axis is off.
    assert "axis" in out


def test_binding_site_center_correct_nm_no_warn():
    """Sanity: valid nm coordinates should NOT trigger the unit
    warning. Use the bundled binding_egfr.yaml values."""
    yaml = """
    receptor:
      pdb_path: benchmarks/dock/1m17.pdb
      binding_site_center_nm: [2.201, 0.025, 5.279]

    entries:
      - name: erlotinib
        smiles: "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"
        dG_bind_kcalmol: -11.86
    """
    rc, out = _run(yaml)
    assert rc == 0
    assert "likely Å, not nm" not in out, (
        f"valid nm coordinate should NOT warn; got:\n{out}")


def test_hydration_wall_time_estimate_matches_m5max_observation():
    """Pin the hydration estimator against the M5 Max ground-truth
    (friend's in-flight FreeSolv-12 run at production parameters
    takes 4-6 h GPU). Our estimator must land in [3 h, 10 h] for
    the bundled freesolv_12.yaml.

    If the hydration estimator reverts to the binding per-atom-
    step model (which gave 4.5 min → 60× underestimate against
    real hardware), this test fires. Similarly catches any
    accidental unit drift (s vs min vs h).
    """
    import io as _io
    import re
    from src.fep.binding import main

    old = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        main(["validate", "benchmarks/fep/freesolv_12.yaml"])
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old
    m = re.search(r"sampled \(GPU\) : ([\d.]+)\s*(s|min|h)", out)
    assert m, f"GPU estimate missing from hydration output:\n{out}"
    val, unit = m.groups()
    gpu_s = float(val) * {"s": 1, "min": 60, "h": 3600}[unit]
    assert 3 * 3600 < gpu_s < 10 * 3600, (
        f"hydration GPU estimate {gpu_s/3600:.1f} h outside "
        "[3 h, 10 h] — M5 Max ground truth is 4-6 h")


def test_wall_time_estimate_within_reasonable_range():
    """Pin the wall-time estimator constants against silent drift.

    Calibrated 2026-04-22 from the 1ubq + methane sampled smoke
    (82.5 s / 1200 steps / 15k-atom complex → 5e-7 s per atom-
    step CPU, 20× faster GPU). For EGFR (6 compounds on 1m17,
    ~110k-atom complex × 6 = 55k atom-steps × 550k steps):

      - CPU ≈ 100 × 0.05 = ~12 h (within [5, 50] h)
      - GPU ≈ 12 h / 20 = ~35 min (within [10, 180] min)

    Wide bands tolerate reasonable refinement but catch a 10×
    reversal. If someone drops an exponent on atom-count scaling
    or forgets the legs/windows factor, this fires.
    """
    from src.fep.binding import main

    old = sys.stdout
    sys.stdout = io.StringIO()
    try:
        main(["validate", "benchmarks/fep/binding_egfr.yaml"])
        out = sys.stdout.getvalue()
    finally:
        sys.stdout = old

    assert "estimated wall time" in out
    # Parse the 'sampled (CPU)' line: expect hours.
    import re
    cpu_m = re.search(
        r"sampled \(CPU\) : ([\d.]+)\s*(s|min|h)", out)
    gpu_m = re.search(
        r"sampled \(GPU\) : ([\d.]+)\s*(s|min|h)", out)
    assert cpu_m, f"CPU estimate missing from output:\n{out}"
    assert gpu_m, f"GPU estimate missing from output:\n{out}"

    def _s(val, unit):
        return (float(val) * {"s": 1, "min": 60, "h": 3600}[unit])

    cpu_s = _s(*cpu_m.groups())
    gpu_s = _s(*gpu_m.groups())
    assert 5 * 3600 < cpu_s < 50 * 3600, (
        f"CPU estimate {cpu_s:.0f}s ({cpu_s/3600:.1f}h) "
        "outside [5h, 50h] — constants drifted")
    assert 10 * 60 < gpu_s < 180 * 60, (
        f"GPU estimate {gpu_s:.0f}s ({gpu_s/60:.1f}min) "
        "outside [10min, 180min] — constants drifted")
    # CPU should always be considerably slower than GPU.
    assert cpu_s > 5 * gpu_s, (
        "CPU should be >> GPU; got CPU={cpu_s}, GPU={gpu_s}")


if __name__ == "__main__":
    funcs = [
        test_valid_binding_yaml_passes,
        test_valid_hydration_yaml_passes_without_receptor,
        test_typod_smiles_fails,
        test_missing_receptor_fails_on_binding_yaml,
        test_undefined_stereo_bond_flagged,
        test_warning_does_not_fail_exit,
        test_error_does_fail_exit,
        test_pinned_stereo_passes,
        test_charged_ligand_flagged,
        test_duplicate_names_flagged,
        test_empty_entries_list_fails,
        test_bundled_streptavidin_yaml_validates,
        test_bundled_freesolv_yaml_validates,
        test_wall_time_estimate_within_reasonable_range,
        test_hydration_wall_time_estimate_matches_m5max_observation,
        test_binding_site_center_unit_mismatch_warns,
        test_binding_site_center_outside_protein_bbox_warns,
        test_binding_site_center_correct_nm_no_warn,
        test_yaml_entry_key_typo_flagged,
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
