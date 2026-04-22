#!/usr/bin/env python3
"""Biologist regression dashboard: run `cellsim fep-binding
validate` on every YAML under benchmarks/fep/ and emit a
one-line-per-YAML summary plus an aggregate pass count.

Zero MD. Zero GPU. < 10 s total. Intended to be run:
  - locally when pulling a fresh clone, to confirm every
    shipped benchmark still parses cleanly against the
    current force fields;
  - in CI (from scripts/ for reuse) to block a PR that
    breaks a YAML.

Output example:

    [bench-all-yamls]
      binding_egfr.yaml         PASS (with warnings)  6 entries  → 38 min GPU
      binding_streptavidin.yaml PASS                  4 entries  → 15 min GPU
      freesolv_12.yaml          PASS                 12 entries  → 6 h   GPU

      3/3 YAMLs validate clean.  Total sampled budget on GPU: ~7 h

Exit 0 if all validate clean (errors are hard; warnings are ok).
Exit 1 if any YAML has hard errors.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
YAML_DIR = REPO_ROOT / "benchmarks" / "fep"


def _parse_minutes_hours(text: str) -> float:
    """Grab the 'sampled (GPU) : X h|min|s' line and return seconds."""
    m = re.search(r"sampled \(GPU\) : ([\d.]+)\s*(s|min|h)", text)
    if not m:
        return 0.0
    val, unit = m.groups()
    return float(val) * {"s": 1, "min": 60, "h": 3600}[unit]


def _fmt_wall(s: float) -> str:
    if s < 90:
        return f"{s:>3.0f} s"
    if s < 3600:
        return f"{s/60:>3.0f} min"
    return f"{s/3600:>3.1f} h"


def main() -> int:
    yamls = sorted(YAML_DIR.glob("*.yaml"))
    if not yamls:
        print(f"bench-all-yamls: no *.yaml files under {YAML_DIR}",
              file=sys.stderr)
        return 2

    print("[bench-all-yamls]")
    n_hard_fail = 0
    total_gpu_s = 0.0
    lines: list[str] = []

    for y in yamls:
        r = subprocess.run(
            ["python", "-m", "src.fep.binding", "validate", str(y)],
            cwd=REPO_ROOT,
            capture_output=True, text=True,
        )
        out = r.stdout + r.stderr
        n_entries_m = re.search(r"entries:\s+(\d+)", out)
        n_entries = int(n_entries_m.group(1)) if n_entries_m else 0
        if r.returncode == 0:
            warn = "(with warnings)" in out
            verdict = "PASS (with warnings)" if warn else "PASS"
        else:
            verdict = "FAIL"
            n_hard_fail += 1
        gpu_s = _parse_minutes_hours(out)
        total_gpu_s += gpu_s
        wall_tag = (
            f" → {_fmt_wall(gpu_s)} GPU" if gpu_s > 0 else "")
        lines.append(
            f"  {y.name:<28s} {verdict:<22s} "
            f"{n_entries:>2d} entries{wall_tag}")

    for line in lines:
        print(line)
    print()
    print(f"  {len(yamls) - n_hard_fail}/{len(yamls)} "
          "YAMLs validate clean.  "
          f"Total sampled budget on GPU: ~{_fmt_wall(total_gpu_s)}")
    return 0 if n_hard_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
