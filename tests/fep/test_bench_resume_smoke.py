"""--resume regression test.

Reads a prior bench CSV, skips the compounds already present,
only scaffolds the missing ones. Validates the path a biologist
takes after a crashed multi-hour run.

This runs scaffold-only (no --sample) so it's CI-cheap
(~30 s on the 3-compound ubiquitin test YAML).
"""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _write_yaml(path: Path, smiles_list):
    lines = [
        "receptor:",
        "  pdb_path: benchmarks/md/1ubq.pdb",
        "  binding_site_center_nm: null",
        "",
        "entries:",
    ]
    for name, smi in smiles_list:
        lines.append(f"  - name: {name}")
        lines.append(f"    smiles: \"{smi}\"")
        lines.append("    dG_bind_kcalmol: -2.0")
    path.write_text("\n".join(lines))


def _seed_csv(path: Path, completed_names):
    cols = ["name", "smiles", "dG_expt_kcalmol", "dG_pred_kcalmol",
            "uncertainty_kcalmol", "residual_kcalmol",
            "ghmc_accept_mean", "ghmc_accept_min",
            "wall_seconds", "ok", "reason"]
    with path.open("w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=cols)
        w.writeheader()
        for name in completed_names:
            w.writerow({
                "name": name, "smiles": "X",  # any non-empty SMILES
                "dG_expt_kcalmol": -2.0,
                "dG_pred_kcalmol": -1.85,
                "uncertainty_kcalmol": 0.4,
                "residual_kcalmol": 0.15,
                "ghmc_accept_mean": 0.78,
                "ghmc_accept_min": 0.72,
                "wall_seconds": 1200.0,
                "ok": True, "reason": "",
            })


def test_resume_skips_completed_compounds():
    """3-entry YAML, CSV pre-seeded with 1 complete row, expect
    --resume to skip that one and scaffold the other 2."""
    with tempfile.TemporaryDirectory(prefix="cellsim_resume_") as tmp:
        tmp = Path(tmp)
        yaml = tmp / "resume.yaml"
        csv_path = tmp / "resume.csv"
        _write_yaml(yaml, [
            ("methane", "C"),
            ("ethane", "CC"),
            ("propane", "CCC"),
        ])
        _seed_csv(csv_path, completed_names=["methane"])

        r = subprocess.run(
            ["python", "-m", "src.fep.binding", "bench",
             str(yaml),
             "--padding", "0.6",
             "--out-csv", str(csv_path),
             "--resume"],
            cwd=REPO_ROOT,
            capture_output=True, text=True,
        )
        assert r.returncode == 0, (
            f"bench --resume failed rc={r.returncode}\n"
            f"stdout: {r.stdout}\nstderr: {r.stderr}")

        out = r.stdout + r.stderr
        assert "resumed:    1 compounds" in out, (
            f"resume banner missing:\n{out[-600:]}")
        assert "(skipped, resumed from CSV)" in out
        # Final CSV should have 3 rows + 1 header.
        with csv_path.open("r", encoding="utf-8-sig") as fi:
            rows = list(csv.DictReader(fi))
        assert len(rows) == 3, (
            f"expected 3 rows (1 resumed + 2 scaffolded); "
            f"got {len(rows)}: {[r['name'] for r in rows]}")

        # Methane row must preserve the PRIOR prediction (-1.85),
        # not be re-scaffolded (which would wipe it to None).
        methane = next(r for r in rows if r["name"] == "methane")
        assert methane["dG_pred_kcalmol"] == "-1.85", (
            f"resumed methane row was overwritten: {methane}")
        # Ethane/propane should be freshly scaffolded (no pred).
        ethane = next(r for r in rows if r["name"] == "ethane")
        assert not ethane["dG_pred_kcalmol"], (
            f"freshly-scaffolded ethane should have empty pred: "
            f"{ethane}")


def test_ctrl_c_exits_cleanly_with_partial_csv():
    """Simulate Ctrl-C during bench: monkey-patch
    compute_absolute_binding_dg to raise KeyboardInterrupt on the
    2nd call. Expected behaviour:
      - bench returns exit 130 (SIGINT convention)
      - output contains 'INTERRUPTED'
      - CSV has the 1 completed row (crash-proof incremental
        write from commit 9732e0b preserves it)
      - output suggests --resume as the restart path
    """
    import tempfile
    with tempfile.TemporaryDirectory(
            prefix="cellsim_ctrlc_") as tmp:
        tmp = Path(tmp)
        yaml = tmp / "ctrlc.yaml"
        csv_path = tmp / "ctrlc.csv"
        _write_yaml(yaml, [
            ("methane", "C"),
            ("ethane", "CC"),
            ("propane", "CCC"),
        ])
        # Use a subprocess that raises on 2nd compound via a
        # wrapper python -c script. Monkey-patching main()
        # from here is fraught because bench imports eagerly.
        harness = f"""
import sys
sys.path.insert(0, "{REPO_ROOT}")
import src.fep.binding as binding
_real = binding.compute_absolute_binding_dg
_counter = {{"n": 0}}
def mock(*a, **kw):
    _counter["n"] += 1
    if _counter["n"] == 2:
        raise KeyboardInterrupt()
    return _real(*a, **kw)
binding.compute_absolute_binding_dg = mock
sys.exit(binding.main([
    "bench", "{yaml}", "--padding", "0.6",
    "--out-csv", "{csv_path}",
]))
"""
        r = subprocess.run(
            ["python", "-c", harness],
            cwd=REPO_ROOT,
            capture_output=True, text=True,
        )
        assert r.returncode == 130, (
            f"Ctrl-C should exit 130; got {r.returncode}\n"
            f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}")
        out = r.stdout + r.stderr
        assert "INTERRUPTED" in out
        assert "--resume" in out, (
            f"partial-run message should suggest --resume\n{out}")
        # CSV should have one completed row (methane).
        assert csv_path.exists()
        with csv_path.open("r", encoding="utf-8-sig") as fi:
            rows = list(csv.DictReader(fi))
        assert len(rows) == 1
        assert rows[0]["name"] == "methane"


def test_resume_without_flag_overwrites():
    """Without --resume, an existing CSV is overwritten and every
    compound is re-scaffolded. Regression guard on default."""
    with tempfile.TemporaryDirectory(prefix="cellsim_resume_") as tmp:
        tmp = Path(tmp)
        yaml = tmp / "noresume.yaml"
        csv_path = tmp / "noresume.csv"
        _write_yaml(yaml, [("methane", "C")])
        _seed_csv(csv_path, completed_names=["methane"])

        r = subprocess.run(
            ["python", "-m", "src.fep.binding", "bench",
             str(yaml),
             "--padding", "0.6",
             "--out-csv", str(csv_path)],   # no --resume
            cwd=REPO_ROOT,
            capture_output=True, text=True,
        )
        assert r.returncode == 0
        assert "resumed:" not in r.stdout  # no resume banner
        # Methane was re-scaffolded, so pred should now be empty.
        with csv_path.open("r", encoding="utf-8-sig") as fi:
            rows = list(csv.DictReader(fi))
        assert len(rows) == 1
        assert not rows[0]["dG_pred_kcalmol"], (
            "without --resume, the seeded prior pred should be "
            f"overwritten; got {rows[0]}")


if __name__ == "__main__":
    funcs = [
        test_resume_skips_completed_compounds,
        test_resume_without_flag_overwrites,
        test_ctrl_c_exits_cleanly_with_partial_csv,
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
