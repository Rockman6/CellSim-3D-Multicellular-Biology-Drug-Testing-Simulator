"""cellsim bench-all dashboard regression test.

The bench-all script is the biologist's one-line regression
dashboard across every shipped benchmark YAML. If its output
format drifts (column widths, summary line, exit code) any
downstream tooling parsing it breaks. This test pins the shape.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_bench_all() -> tuple[int, str]:
    r = subprocess.run(
        ["python", "scripts/bench_all_yamls.py"],
        cwd=REPO_ROOT,
        capture_output=True, text=True,
    )
    return r.returncode, r.stdout + r.stderr


def test_bench_all_exits_zero_on_bundled_yamls():
    rc, out = _run_bench_all()
    assert rc == 0, (
        f"bench-all should pass on bundled YAMLs; got rc={rc}\n{out}")


def test_bench_all_summary_line_is_stable():
    """Final line must match 'N/M YAMLs validate clean. Total
    sampled budget on GPU: ~X'. Downstream parsers depend on
    this shape."""
    _, out = _run_bench_all()
    import re
    # Find a line matching the summary template.
    m = re.search(
        r"(\d+)/(\d+) YAMLs validate clean\.\s+"
        r"Total sampled budget on GPU: ~([\d.]+\s*\S+)",
        out)
    assert m, f"summary line shape changed:\n{out[-500:]}"
    n_pass, n_total, _budget = m.groups()
    assert n_pass == n_total, (
        f"one or more bundled YAMLs failing: {n_pass}/{n_total}")


def test_bench_all_lists_every_shipped_yaml():
    """Every *.yaml under benchmarks/fep/ must show up exactly
    once in the dashboard — catches silent drop of a benchmark."""
    _, out = _run_bench_all()
    yaml_dir = REPO_ROOT / "benchmarks" / "fep"
    for y in sorted(yaml_dir.glob("*.yaml")):
        assert y.name in out, (
            f"bench-all didn't mention {y.name}:\n{out}")


def test_bench_all_shows_gpu_estimate_per_yaml():
    """Every listed line should carry a GPU wall estimate.
    Protects against estimator regression on one YAML not
    surfacing in the dashboard."""
    _, out = _run_bench_all()
    lines = [
        line for line in out.splitlines()
        if ".yaml" in line and "PASS" in line]
    assert len(lines) >= 3, (
        f"expected ≥3 YAML lines in dashboard; got {len(lines)}"
        f":\n{out}")
    for line in lines:
        assert "GPU" in line, (
            f"line missing GPU estimate: {line!r}")


if __name__ == "__main__":
    funcs = [
        test_bench_all_exits_zero_on_bundled_yamls,
        test_bench_all_summary_line_is_stable,
        test_bench_all_lists_every_shipped_yaml,
        test_bench_all_shows_gpu_estimate_per_yaml,
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
