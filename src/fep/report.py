"""FreeSolv tarball analyser — Milestone A gate verdict in 30 seconds.

Consumes a `freesolv_m5max_<stamp>.tar.gz` (or an already-extracted
run directory, or just a raw `freesolv_results.csv`) and emits a
biologist-readable pass/fail report:

    - Overall gate verdict: MAE ≤ 1.5 kcal/mol?
    - Per-compound residual table ranked worst → best.
    - Sign-correctness audit (methane MUST be positive; acetamide
      MUST be strongly negative; etc.).
    - GHMC acceptance audit (parsed from run.log when present;
      flags anything below 0.70 per the prof's checklist).
    - Wall-time distribution and throughput (tells us what
      Campaign-2 dosage-response screens will cost).
    - Environment / doctor provenance copied verbatim so anyone
      reading the report can reproduce the run from fresh.

Used two ways:
  1. Locally when the friend's tarball arrives: one command, know
     within 30 s whether Milestone A cleared.
  2. In CI on every FEP change: drops a markdown summary next to
     the CSV and non-zero-exits on regression.

This is reader-facing: biologists who don't read source will only
see the markdown output, so wording of the verdict + per-compound
flags has to be self-explanatory. No jargon-only notes; every
flag has a one-line "what this means for your experiment".
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import re
import shutil
import sys
import tarfile
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# FreeSolv gate thresholds — hard-coded because they're literally
# the Milestone A exit criteria Henry + the prof agreed on; they
# belong in code, not config.
GATE_MAE_KCALMOL = 1.5
GATE_GHMC_ACCEPT_MIN = 0.70
# Per-compound sign-correctness: ΔG_hyd sign must match experiment
# for at least METHANE (hydrophobic → positive) and ACETAMIDE
# (strongly polar → strongly negative). Getting methane wrong on a
# long GPU run would be the professor's "still undersampled" finding
# and is a strict stop-the-campaign signal.
_SIGN_CRITICAL_NAMES = {"methane", "acetamide"}


@dataclass
class CompoundRow:
    name: str
    smiles: str
    dG_expt_kcalmol: float
    dG_pred_kcalmol: Optional[float] = None
    uncertainty_kcalmol: Optional[float] = None
    residual_kcalmol: Optional[float] = None
    wall_seconds: Optional[float] = None
    ok: bool = False
    reason: str = ""
    # Parsed from run.log when available — absent on legacy runs.
    ghmc_accept_mean: Optional[float] = None
    ghmc_accept_min: Optional[float] = None
    sign_correct: Optional[bool] = None
    # |residual| ≤ uncertainty (statistically consistent with
    # experiment). None when residual/σ is missing (scaffold-only).
    within_sigma: Optional[bool] = None
    # Warnings attached to this row: e.g. "abs residual > 3 kcal/mol",
    # "GHMC mean < 0.70", "NaN in uncertainty".
    flags: list = field(default_factory=list)


@dataclass
class ReportResult:
    source: str
    n_total: int = 0
    n_ok: int = 0
    mae_kcalmol: Optional[float] = None
    rmse_kcalmol: Optional[float] = None
    pearson_r: Optional[float] = None
    spearman_rho: Optional[float] = None
    kendall_tau: Optional[float] = None
    # MAE gate this report was scored against; stored so the
    # markdown + parity PNG render the right threshold.
    mae_gate_kcalmol: float = GATE_MAE_KCALMOL
    # Expected row count (from --expected flag). When set and
    # the CSV has fewer rows, the run is flagged as 'partial'
    # and cannot pass even if those rows all look fine.
    expected_rows: Optional[int] = None
    is_partial: bool = False
    # 'hydration' | 'binding' — drives the sign-critical rule and
    # the markdown header wording.
    yaml_kind: str = "hydration"
    # Aggregate GHMC stats (across all windows × compounds that
    # report data); None if no GHMC info in log.
    ghmc_accept_overall_mean: Optional[float] = None
    ghmc_accept_overall_min: Optional[float] = None
    ghmc_accept_compound_count: int = 0
    wall_total_seconds: Optional[float] = None
    wall_mean_per_compound_seconds: Optional[float] = None
    wall_max_per_compound_seconds: Optional[float] = None
    rows: list = field(default_factory=list)
    # Gate verdicts — True/False/None (None = couldn't evaluate).
    pass_mae: Optional[bool] = None
    pass_ghmc: Optional[bool] = None
    pass_sign_critical: Optional[bool] = None
    pass_overall: Optional[bool] = None
    # Provenance: env + doctor output verbatim (truncated) so the
    # markdown report includes "this is exactly what ran".
    env_log: Optional[str] = None
    doctor_tail: Optional[str] = None
    platform_line: Optional[str] = None
    git_commit: Optional[str] = None


def _autodiscover_latest() -> Optional[Path]:
    """Find the most recent run artefact in the repo: either a
    tarball at the repo root or a run/fep/<stamp>/ directory. The
    biologist workflow is "I just ran the gate, now what?" — they
    shouldn't have to remember a timestamp.

    Returns the path (tarball or directory) with the newest mtime,
    or None if nothing's there.
    """
    cwd = Path.cwd()
    candidates: list[tuple[float, Path]] = []

    for tb in cwd.glob("freesolv_m5max_*.tar.gz"):
        candidates.append((tb.stat().st_mtime, tb))
    for tb in cwd.glob("freesolv_*.tar.gz"):
        candidates.append((tb.stat().st_mtime, tb))

    run_root = cwd / "run" / "fep"
    if run_root.is_dir():
        for sub in run_root.iterdir():
            csv = sub / "freesolv_results.csv"
            if csv.is_file():
                candidates.append((sub.stat().st_mtime, sub))

    if not candidates:
        return None
    candidates.sort(reverse=True)
    return candidates[0][1]


def _resolve_inputs(path: Path) -> tuple[Path, Path]:
    """Resolve (run_dir, csv_path). Accepts:

    - a raw CSV (e.g. `bench_mini_results.csv`): csv_path = path;
      run_dir = path.parent (used for optional env.log / run.log
      siblings if present).
    - a directory: find any `*results*.csv` under it (prefers
      `freesolv_results.csv` for back-compat with Milestone A
      tarballs).
    - a tarball: extract + recurse.

    Raises ValueError with actionable messages — biologists running
    this shouldn't have to read stack traces.
    """
    if path.is_file():
        if path.suffix == ".csv":
            return path.parent, path
        if path.suffixes[-2:] == [".tar", ".gz"] or path.suffix == ".tgz":
            tmp = Path(tempfile.mkdtemp(prefix="fep_report_"))
            with tarfile.open(path, "r:gz") as tf:
                tf.extractall(tmp)
            return _resolve_inputs(tmp)
    if path.is_dir():
        cands = list(path.rglob("freesolv_results.csv"))
        if cands:
            return cands[0].parent, cands[0]
        cands = sorted(path.rglob("*results*.csv"))
        if cands:
            return cands[0].parent, cands[0]
        raise ValueError(
            f"directory {path} has no *results*.csv under it; "
            f"top-level: "
            f"{sorted(p.name for p in path.iterdir())[:20]}")
    raise ValueError(f"no such file or directory: {path}")


def _read_csv(csv_path: Path) -> list[CompoundRow]:
    rows: list[CompoundRow] = []
    with csv_path.open("r", encoding="utf-8-sig") as fi:
        reader = csv.DictReader(fi)
        for r in reader:
            def _f(key):
                v = r.get(key, "")
                if v in (None, "", "None"):
                    return None
                try:
                    return float(v)
                except ValueError:
                    return None

            def _b(key):
                v = (r.get(key, "") or "").strip().lower()
                return v in ("true", "1", "yes", "ok")

            row = CompoundRow(
                name=r.get("name", "").strip(),
                smiles=r.get("smiles", "").strip(),
                dG_expt_kcalmol=float(
                    r.get("dG_expt_kcalmol") or "nan"),
                dG_pred_kcalmol=_f("dG_pred_kcalmol"),
                uncertainty_kcalmol=_f("uncertainty_kcalmol"),
                residual_kcalmol=_f("residual_kcalmol"),
                wall_seconds=_f("wall_seconds"),
                ok=_b("ok"),
                reason=r.get("reason", "").strip(),
                ghmc_accept_mean=_f("ghmc_accept_mean"),
                ghmc_accept_min=_f("ghmc_accept_min"),
            )
            rows.append(row)
    return rows


_GHMC_LINE_RE = re.compile(
    r"GHMC accept\s+mean=(?P<mean>[\d.]+)%\s+min=(?P<min>[\d.]+)%")
_COMPOUND_HEADER_RE = re.compile(r"^\[freesolv\]\s+(\S+)\s+(.+)$")
_PLATFORM_RE = re.compile(
    r"FEP sampling platform:\s+(\S+)",
    re.IGNORECASE)
_GIT_COMMIT_RE = re.compile(r"git commit:\s+([0-9a-f]{7,40})", re.I)


def _parse_run_log(log_path: Path, rows: list[CompoundRow]) -> dict:
    """Best-effort enrich the CompoundRow list with per-compound GHMC
    acceptance data from `run.log` (if it was printed there). Returns
    aggregate metadata: platform, git commit, overall GHMC stats.

    The run.log is created by `tee` on top of `run_freesolv_validation`
    output; per-compound boundaries are the `[freesolv] {name} {smiles}`
    header lines. GHMC lines only appear if a caller printed
    AlchemicalSamplingResult.summary() — the current
    `run_freesolv_validation` does NOT, so on the first M5 Max run
    this section will just be empty and we report "no GHMC data in
    log (legacy format)". The parser still works when Phase-2 wires
    per-compound acceptance through the CSV.

    Provenance fallback: if run.log doesn't carry `git commit:` (the
    shell header block tees into env.log on some pre-2026-04-21 runs
    and into run.log after), also scan env.log and doctor.log in the
    same directory. Biologists care that *some* log captured the
    commit — not which one.
    """
    meta: dict = {"platform": None, "git_commit": None,
                  "ghmc_means": [], "ghmc_mins": []}
    if not log_path.exists():
        # Fall through to the sibling-log scan below so we at least
        # pick up git_commit / platform from env.log when present.
        return _scan_sibling_logs_for_provenance(
            log_path.parent, meta)

    name_to_row = {r.name: r for r in rows}
    current: Optional[CompoundRow] = None
    try:
        txt = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return meta
    for line in txt.splitlines():
        m = _PLATFORM_RE.search(line)
        if m and not meta["platform"]:
            meta["platform"] = m.group(1).strip()
            continue
        m = _GIT_COMMIT_RE.search(line)
        if m and not meta["git_commit"]:
            meta["git_commit"] = m.group(1).strip()
            continue
        m = _COMPOUND_HEADER_RE.match(line.strip())
        if m:
            current = name_to_row.get(m.group(1))
            continue
        m = _GHMC_LINE_RE.search(line)
        if m and current is not None:
            mean = float(m.group("mean")) / 100.0
            mn = float(m.group("min")) / 100.0
            # If multiple legs per compound emit two lines, keep the
            # worst (more conservative audit).
            if (current.ghmc_accept_mean is None
                    or mean < current.ghmc_accept_mean):
                current.ghmc_accept_mean = mean
            if (current.ghmc_accept_min is None
                    or mn < current.ghmc_accept_min):
                current.ghmc_accept_min = mn
            meta["ghmc_means"].append(mean)
            meta["ghmc_mins"].append(mn)
    # If provenance still missing, look in sibling logs.
    if not meta["git_commit"] or not meta["platform"]:
        meta = _scan_sibling_logs_for_provenance(log_path.parent, meta)
    return meta


def _scan_sibling_logs_for_provenance(
    dir_path: Path, meta: dict,
) -> dict:
    """Scan env.log and doctor.log next to run.log for git_commit
    and platform. These files are written by the shell wrapper
    scripts; the header-tee pattern can land the provenance in any
    of them depending on the script version."""
    for sibling in ("env.log", "doctor.log", "header.log"):
        p = dir_path / sibling
        if not p.exists():
            continue
        try:
            t = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if not meta["git_commit"]:
            m = _GIT_COMMIT_RE.search(t)
            if m:
                meta["git_commit"] = m.group(1).strip()
        if not meta["platform"]:
            m = _PLATFORM_RE.search(t)
            if m:
                meta["platform"] = m.group(1).strip()
        if meta["git_commit"] and meta["platform"]:
            break
    return meta


def _spearman(x: list[float], y: list[float]) -> Optional[float]:
    """Rank-order Pearson — no scipy dep, keeps the CLI lean."""
    if len(x) < 3:
        return None
    def _ranks(a):
        order = sorted(range(len(a)), key=lambda i: a[i])
        ranks = [0.0] * len(a)
        for rank, idx in enumerate(order):
            ranks[idx] = float(rank)
        return ranks
    rx = _ranks(x)
    ry = _ranks(y)
    return _pearson(rx, ry)


def _pearson(x: list[float], y: list[float]) -> Optional[float]:
    if len(x) < 3:
        return None
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    num = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    dx = math.sqrt(sum((xi - mx) ** 2 for xi in x))
    dy = math.sqrt(sum((yi - my) ** 2 for yi in y))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def _kendall_tau(x: list[float], y: list[float]) -> Optional[float]:
    """Plain Kendall τ-b — O(n²) but n=12 in our case, so fine."""
    n = len(x)
    if n < 3:
        return None
    concordant = discordant = tx = ty = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            if dx == 0 and dy == 0:
                tx += 1; ty += 1
            elif dx == 0:
                tx += 1
            elif dy == 0:
                ty += 1
            elif (dx > 0) == (dy > 0):
                concordant += 1
            else:
                discordant += 1
    n_pairs = n * (n - 1) / 2
    denom_x = math.sqrt((n_pairs - tx) * (n_pairs - ty))
    if denom_x == 0:
        return None
    return (concordant - discordant) / denom_x


def analyse(path: str | Path,
            *,
            mae_gate_kcalmol: float = GATE_MAE_KCALMOL,
            expected_rows: Optional[int] = None,
            yaml_kind: str = "hydration",
            ) -> ReportResult:
    """Top-level: path in, ReportResult out.

    `mae_gate_kcalmol` lets a caller override the default
    hydration-era 1.5 kcal/mol gate — binding runs (streptavidin
    etc.) use 2.0 per BENCHMARKS.md Milestone B scope.

    `expected_rows` is the number of compounds the YAML asked for.
    If the CSV has fewer rows, the run is flagged as 'partial'
    and cannot PASS regardless of MAE — a 6/12 run that happens to
    have good MAE on its 6 rows is NOT a Milestone A pass; the
    other 6 compounds might have been the hard ones. Matches the
    professor's 'do not pass on incomplete data' rule.

    `yaml_kind` selects the sign-critical rule:
      'hydration' → methane must be ≥0 and acetamide must be <0
          (rule applies per the specific FreeSolv benchmark).
      'binding'   → every compound must have ΔG_bind < 0; an
          entry in a binding YAML that the sampler predicts as
          non-binding is a loud 'pipeline broke' signal.
    """
    path = Path(path)
    run_dir, csv_path = _resolve_inputs(path)

    rows = _read_csv(csv_path)
    meta = _parse_run_log(run_dir / "run.log", rows)

    # Per-compound flags + sign check.
    for r in rows:
        if not r.ok:
            continue
        if (r.dG_pred_kcalmol is not None
                and not math.isnan(r.dG_expt_kcalmol)):
            r.sign_correct = (
                (r.dG_pred_kcalmol >= 0) == (r.dG_expt_kcalmol >= 0)
                or abs(r.dG_expt_kcalmol) < 0.3)   # near-zero ok either way
        # within_sigma: |residual| ≤ uncertainty. Single source of
        # truth — markdown + table.csv + JSON all read from the
        # row attribute now instead of recomputing.
        if (r.residual_kcalmol is not None
                and r.uncertainty_kcalmol is not None
                and r.uncertainty_kcalmol > 0):
            r.within_sigma = (
                abs(r.residual_kcalmol) <= r.uncertainty_kcalmol)
        if (r.residual_kcalmol is not None
                and abs(r.residual_kcalmol) > 3.0):
            r.flags.append(
                f"|residual| {abs(r.residual_kcalmol):.2f} > 3 "
                f"kcal/mol (outlier)")
        if r.sign_correct is False:
            r.flags.append(
                "sign WRONG — pred and expt have opposite sign "
                "(loud signal: sampling not converged or FF failure)")
        if (r.ghmc_accept_mean is not None
                and r.ghmc_accept_mean < GATE_GHMC_ACCEPT_MIN):
            r.flags.append(
                f"GHMC mean acceptance {r.ghmc_accept_mean:.0%} < "
                f"{GATE_GHMC_ACCEPT_MIN:.0%} threshold — trust the "
                "ΔG less; timestep likely too large")

    # `n_ok` counts rows that ran without errors (including scaffold-
    # only ones where dG_pred is None); `ok_rows` is the stricter
    # subset used for statistics (MAE/RMSE/Pearson/etc.), which
    # requires an actual predicted ΔG.
    ok_rows = [r for r in rows
               if r.ok and r.dG_pred_kcalmol is not None]
    n_total = len(rows)
    n_ok = sum(1 for r in rows if r.ok)

    expts = [r.dG_expt_kcalmol for r in ok_rows]
    preds = [r.dG_pred_kcalmol for r in ok_rows]
    resids = [r.residual_kcalmol for r in ok_rows
              if r.residual_kcalmol is not None]
    walls = [r.wall_seconds for r in rows
             if r.wall_seconds is not None]

    mae = (sum(abs(v) for v in resids) / len(resids)
           if resids else None)
    rmse = (math.sqrt(sum(v ** 2 for v in resids) / len(resids))
            if resids else None)

    result = ReportResult(
        source=str(path),
        n_total=n_total, n_ok=n_ok,
        mae_kcalmol=mae, rmse_kcalmol=rmse,
        mae_gate_kcalmol=mae_gate_kcalmol,
        expected_rows=expected_rows,
        is_partial=(
            expected_rows is not None and n_total < expected_rows),
        yaml_kind=yaml_kind,
        pearson_r=_pearson(preds, expts),
        spearman_rho=_spearman(preds, expts),
        kendall_tau=_kendall_tau(preds, expts),
        rows=rows,
        wall_total_seconds=sum(walls) if walls else None,
        wall_mean_per_compound_seconds=(
            sum(walls) / len(walls) if walls else None),
        wall_max_per_compound_seconds=max(walls) if walls else None,
        platform_line=meta.get("platform"),
        git_commit=meta.get("git_commit"),
    )

    # Aggregate GHMC (only populated if Phase-2 CSVs carry the cols,
    # OR the legacy log had GHMC lines — future runs will).
    means = [r.ghmc_accept_mean for r in rows
             if r.ghmc_accept_mean is not None]
    mins = [r.ghmc_accept_min for r in rows
            if r.ghmc_accept_min is not None]
    if means:
        result.ghmc_accept_overall_mean = sum(means) / len(means)
        result.ghmc_accept_compound_count = len(means)
    if mins:
        result.ghmc_accept_overall_min = min(mins)

    # Gate verdicts. Tri-state:
    #   True   = gate passed
    #   False  = gate evaluable and failed
    #   None   = gate not evaluable from this CSV (scaffold-only
    #            run, missing column, no critical compounds etc.)
    if mae is not None:
        result.pass_mae = mae <= mae_gate_kcalmol
    else:
        result.pass_mae = None
    if means:
        result.pass_ghmc = all(
            m >= GATE_GHMC_ACCEPT_MIN for m in means)
    else:
        # No acceptance data → can't evaluate; warn but don't fail.
        result.pass_ghmc = None

    # Sign-critical check. Rule depends on YAML kind:
    if yaml_kind == "binding":
        # Every binder must predict ΔG_bind < 0. A positive
        # prediction on a compound labelled as a binder is a loud
        # 'sampler broke / restraint escaped / force field wrong'
        # signal and must fail the gate regardless of MAE.
        bind_rows = [
            r for r in rows
            if r.ok and r.dG_pred_kcalmol is not None]
        if bind_rows:
            result.pass_sign_critical = all(
                r.dG_pred_kcalmol < 0 for r in bind_rows)
            # Tag per-row flag for the markdown table.
            for r in bind_rows:
                if r.dG_pred_kcalmol >= 0:
                    r.sign_correct = False
                    if "predicted non-binder" not in " ".join(r.flags):
                        r.flags.append(
                            "predicted ΔG_bind ≥ 0 (non-binder) — "
                            "unexpected for a curated binding set")
                else:
                    r.sign_correct = True
        else:
            result.pass_sign_critical = None
    else:
        # Hydration: methane + acetamide must have correct sign.
        # Only count rows that have a prediction — scaffold-only
        # rows have sign_correct=None and would incorrectly fail.
        sign_critical = [r for r in rows
                         if r.name in _SIGN_CRITICAL_NAMES
                         and r.ok
                         and r.sign_correct is not None]
        if sign_critical:
            result.pass_sign_critical = all(
                r.sign_correct for r in sign_critical)
        else:
            result.pass_sign_critical = None

    # Overall verdict tri-state:
    #   True  = at least one gate evaluable and all evaluable pass
    #   False = at least one evaluable gate failed, OR partial run
    #           (< expected rows present) regardless of per-row MAE
    #   None  = nothing evaluable (scaffold-only / empty CSV)
    checks = [g for g in (result.pass_mae,
                          result.pass_ghmc,
                          result.pass_sign_critical)
              if g is not None]
    if result.is_partial:
        # A 6/12 run can't pass Milestone A even if the 6 look good —
        # the other 6 might be the hard ones. Force FAIL-or-None.
        if not checks:
            result.pass_overall = None
        else:
            result.pass_overall = False
    elif not checks:
        result.pass_overall = None
    else:
        result.pass_overall = all(checks)

    # Provenance — optional logs.
    env_path = run_dir / "env.log"
    if env_path.exists():
        result.env_log = env_path.read_text(
            encoding="utf-8", errors="replace")[:4000]
    doctor_path = run_dir / "doctor.log"
    if doctor_path.exists():
        txt = doctor_path.read_text(
            encoding="utf-8", errors="replace").splitlines()
        result.doctor_tail = "\n".join(txt[-15:])

    return result


def format_markdown(r: ReportResult) -> str:
    """Biologist-readable summary. Must stand alone as a lab-
    notebook entry — no external references required. The order
    matters: verdict first, then per-compound detail, then
    reproducibility metadata."""
    lines: list[str] = []
    verdict_icon = {
        True: "PASS", False: "FAIL", None: "inconclusive"}
    overall_icon = (
        "PASS" if r.pass_overall
        else "FAIL" if r.pass_overall is False
        else "inconclusive")
    if r.is_partial:
        overall_icon = "FAIL (partial run)" if r.pass_overall is False else "partial"
    # Title reflects the YAML kind the run was scored against so
    # biologists reading the markdown know immediately whether
    # they're looking at a hydration or binding report, not the
    # legacy FreeSolv-era generic title.
    if r.yaml_kind == "binding":
        _title = "Binding FEP report"
    elif r.yaml_kind == "hydration":
        _title = "Hydration FEP report"
    else:
        _title = "FEP report"
    lines.append(f"# {_title} — {overall_icon}")
    lines.append("")
    lines.append(f"- source: `{r.source}`")
    if r.platform_line:
        lines.append(f"- platform: **{r.platform_line}**")
    if r.git_commit:
        lines.append(f"- git commit: `{r.git_commit}`")
    if r.expected_rows is not None:
        completeness = (
            f" ({r.n_total}/{r.expected_rows} expected"
            + (", **PARTIAL — Milestone A cannot pass**"
               if r.is_partial else "")
            + ")")
    else:
        completeness = ""
    lines.append(
        f"- compounds: {r.n_ok} ok / {r.n_total} total"
        f"{completeness}")
    if r.wall_total_seconds is not None:
        lines.append(
            f"- wall time: total {r.wall_total_seconds/60:.1f} min  "
            f"(mean {r.wall_mean_per_compound_seconds/60:.1f} min / "
            f"max {r.wall_max_per_compound_seconds/60:.1f} min per "
            "compound)")
    lines.append("")

    lines.append("## Gate verdict")
    lines.append("")
    mae_str = (f"{r.mae_kcalmol:.2f}" if r.mae_kcalmol is not None
               else "n/a")
    lines.append(
        f"- **MAE gate** (≤ {r.mae_gate_kcalmol} kcal/mol): "
        f"{verdict_icon[r.pass_mae]}  (MAE = {mae_str} kcal/mol)")
    if r.pass_ghmc is None:
        lines.append(
            "- **GHMC acceptance gate** "
            f"(≥ {GATE_GHMC_ACCEPT_MIN:.0%} per compound): "
            "inconclusive (no per-compound acceptance data in "
            "log — this is the legacy CSV schema; Phase-2 runs "
            "will surface it)")
    else:
        gmean = (f"{r.ghmc_accept_overall_mean:.0%}"
                 if r.ghmc_accept_overall_mean is not None
                 else "n/a")
        gmin = (f"{r.ghmc_accept_overall_min:.0%}"
                if r.ghmc_accept_overall_min is not None
                else "n/a")
        lines.append(
            "- **GHMC acceptance gate** "
            f"(≥ {GATE_GHMC_ACCEPT_MIN:.0%} per compound): "
            f"{verdict_icon[r.pass_ghmc]}  "
            f"(overall mean {gmean}, worst {gmin}, "
            f"{r.ghmc_accept_compound_count} compounds report)")
    if r.yaml_kind == "binding":
        sign_label = "every ΔG_bind < 0"
    else:
        sign_label = "/".join(sorted(_SIGN_CRITICAL_NAMES))
    if r.pass_sign_critical is None:
        lines.append(
            "- **Sign-critical check** "
            f"({sign_label}): "
            "inconclusive (no rows for critical compounds)")
    else:
        lines.append(
            "- **Sign-critical check** "
            f"({sign_label}): "
            f"{verdict_icon[r.pass_sign_critical]}")
    lines.append("")

    if r.mae_kcalmol is not None or r.pearson_r is not None:
        lines.append("## Aggregate accuracy")
        lines.append("")
        if r.mae_kcalmol is not None:
            lines.append(f"- MAE  = {r.mae_kcalmol:.3f} kcal/mol")
        if r.rmse_kcalmol is not None:
            lines.append(f"- RMSE = {r.rmse_kcalmol:.3f} kcal/mol")
        if r.pearson_r is not None:
            lines.append(f"- Pearson r    = {r.pearson_r:+.3f}")
        if r.spearman_rho is not None:
            lines.append(f"- Spearman ρ   = {r.spearman_rho:+.3f}")
        if r.kendall_tau is not None:
            lines.append(f"- Kendall τ    = {r.kendall_tau:+.3f}")
        # Within-σ summary: how many predictions agree with
        # experiment inside the MBAR error bar? If it's most of
        # them, the FF is accurate. If it's 0/N, systematic
        # under- or over-binding.
        within_rows = [
            rr for rr in r.rows
            if rr.ok and rr.within_sigma is not None]
        if within_rows:
            within_ok = sum(
                1 for rr in within_rows if rr.within_sigma)
            lines.append(
                f"- within σ     = {within_ok}/{len(within_rows)} "
                f"compounds ({100*within_ok/len(within_rows):.0f}%) "
                "agree with experiment inside MBAR uncertainty")
        lines.append("")

    # Per-compound table — sorted by |residual| descending so the
    # worst offenders are at the top.
    lines.append("## Per-compound results")
    lines.append("")
    # 'within-σ' column: '✓' if |residual| ≤ uncertainty (prediction
    # statistically consistent with expt), '—' if we lack data.
    # Biologists skim this column to spot which residuals are
    # actually bigger than the error bar.
    hdr = ("| name | smiles | expt | pred | ± | residual | "
           "|resid| | within σ | wall (s) | flags |")
    sep = ("|------|--------|------|------|---|----------|"
           "-----:|---------:|---------:|-------|")
    lines.append(hdr)
    lines.append(sep)

    def sort_key(rr: CompoundRow) -> float:
        if rr.residual_kcalmol is None:
            return -float("inf")
        return abs(rr.residual_kcalmol)

    for rr in sorted(r.rows, key=sort_key, reverse=True):
        if not rr.ok:
            lines.append(
                f"| {rr.name} | `{rr.smiles}` | "
                f"{rr.dG_expt_kcalmol:+.2f} | FAIL | — | — | — | "
                f"— | {(rr.wall_seconds or 0):.0f} | "
                f"{rr.reason[:60]} |")
            continue
        # Scaffold-only row: ok=True but pred is None. Render as
        # "scaffolded" so bench pre-flight runs don't 500 the
        # formatter.
        if rr.dG_pred_kcalmol is None:
            lines.append(
                f"| {rr.name} | `{rr.smiles}` | "
                f"{rr.dG_expt_kcalmol:+.2f} | scaffolded | — | "
                f"— | — | — | "
                f"{(rr.wall_seconds or 0):.0f} | "
                "(no sample) |")
            continue
        resid = rr.residual_kcalmol
        resid_abs = abs(resid) if resid is not None else None
        unc = rr.uncertainty_kcalmol
        if rr.within_sigma is None:
            within_sigma = "—"
        else:
            within_sigma = "✓" if rr.within_sigma else ""
        flag_str = "; ".join(rr.flags) if rr.flags else ""
        if rr.ghmc_accept_mean is not None:
            flag_str = (flag_str + "; " if flag_str else "") + (
                f"GHMC {rr.ghmc_accept_mean:.0%}")
        if rr.sign_correct is False:
            flag_str = ("SIGN WRONG"
                        + ("; " + flag_str if flag_str else ""))
        lines.append(
            f"| {rr.name} | `{rr.smiles}` | "
            f"{rr.dG_expt_kcalmol:+.2f} | "
            f"{rr.dG_pred_kcalmol:+.2f} | "
            f"{(unc or 0):.2f} | "
            f"{(resid if resid is not None else 0):+.2f} | "
            f"{(resid_abs if resid_abs is not None else 0):.2f} | "
            f"{within_sigma} | "
            f"{(rr.wall_seconds or 0):.0f} | "
            f"{flag_str} |")
    lines.append("")

    # Notable flags section — pulls sign errors + GHMC failures
    # to the top so the reader doesn't have to scan the table.
    notable = [rr for rr in r.rows
               if rr.ok and (rr.flags or rr.sign_correct is False)]
    if notable:
        lines.append("## Notable flags")
        lines.append("")
        for rr in notable:
            all_flags = list(rr.flags)
            if rr.sign_correct is False:
                all_flags.insert(0, "SIGN WRONG")
            lines.append(f"- **{rr.name}**: " + "; ".join(all_flags))
        lines.append("")

    if r.env_log:
        lines.append("## Environment (env.log)")
        lines.append("")
        lines.append("```")
        lines.append(r.env_log.strip())
        lines.append("```")
        lines.append("")

    if r.doctor_tail:
        lines.append("## cellsim doctor (last 15 lines)")
        lines.append("")
        lines.append("```")
        lines.append(r.doctor_tail.strip())
        lines.append("```")
        lines.append("")

    lines.append(
        "_Report generated by `cellsim fep-report`. "
        "Source: [src/fep/report.py](src/fep/report.py)._")
    return "\n".join(lines)


def _write_parity_png(r: ReportResult, out_path: Path) -> bool:
    """Scatter plot: predicted vs experimental with y=x reference
    and per-point error bars. Returns True on success.

    matplotlib is already in the cellsim env (via seaborn deps of
    pandas), so this should always succeed. If not, we skip and
    note in the report — the markdown stands on its own.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        logger.warning("matplotlib unavailable: %s", e)
        return False
    ok = [rr for rr in r.rows
          if rr.ok and rr.dG_pred_kcalmol is not None]
    if not ok:
        return False
    x = [rr.dG_expt_kcalmol for rr in ok]
    y = [rr.dG_pred_kcalmol for rr in ok]
    err = [rr.uncertainty_kcalmol or 0.0 for rr in ok]

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    # y = x reference band at ±1.5 (the MAE gate)
    lo = min(min(x), min(y)) - 1
    hi = max(max(x), max(y)) + 1
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.5)
    ax.fill_between(
        [lo, hi],
        [lo - r.mae_gate_kcalmol, hi - r.mae_gate_kcalmol],
        [lo + r.mae_gate_kcalmol, hi + r.mae_gate_kcalmol],
        color="gray", alpha=0.15,
        label=f"±{r.mae_gate_kcalmol} kcal/mol gate")
    ax.errorbar(x, y, yerr=err, fmt="o", capsize=3, lw=1)
    for rr, xi, yi in zip(ok, x, y):
        ax.annotate(rr.name, (xi, yi),
                    xytext=(4, 4), textcoords="offset points",
                    fontsize=7, alpha=0.75)
    # Axis + title labels mirror the kind so a biologist reading
    # a binding parity PNG isn't confused by 'ΔG_hyd' axis labels.
    if r.yaml_kind == "binding":
        label = "ΔG_bind"
        title_prefix = "Binding FEP parity"
    else:
        label = "ΔG_hyd"
        title_prefix = "Hydration FEP parity"
    ax.set_xlabel(f"experimental {label} (kcal/mol)")
    ax.set_ylabel(f"predicted {label} (kcal/mol)")
    title_bits = [f"MAE = {r.mae_kcalmol:.2f}" if r.mae_kcalmol else ""]
    if r.pearson_r is not None:
        title_bits.append(f"r = {r.pearson_r:+.2f}")
    ax.set_title(
        f"{title_prefix} — " + ", ".join(t for t in title_bits if t))
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", "box")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


def _write_table_csv(r: ReportResult, out_path: Path) -> None:
    """Normalised CSV (same schema as input plus flag + sign cols).
    Used by anyone consuming the report programmatically — e.g.
    Campaign-2's prior-emitter may want the residual column but
    not care about free text `reason`."""
    cols = ["name", "smiles", "dG_expt_kcalmol", "dG_pred_kcalmol",
            "uncertainty_kcalmol", "residual_kcalmol", "abs_residual",
            "within_sigma", "sign_correct",
            "ghmc_accept_mean", "ghmc_accept_min",
            "wall_seconds", "ok", "reason", "flags"]
    with out_path.open(
            "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=cols)
        w.writeheader()
        for rr in r.rows:
            abs_r = (abs(rr.residual_kcalmol)
                     if rr.residual_kcalmol is not None else None)
            # within_sigma is now populated once in analyse()
            # and read from the row attribute here.
            w.writerow({
                "name": rr.name,
                "smiles": rr.smiles,
                "dG_expt_kcalmol": rr.dG_expt_kcalmol,
                "dG_pred_kcalmol": rr.dG_pred_kcalmol,
                "uncertainty_kcalmol": rr.uncertainty_kcalmol,
                "residual_kcalmol": rr.residual_kcalmol,
                "abs_residual": abs_r,
                "within_sigma": rr.within_sigma,
                "sign_correct": rr.sign_correct,
                "ghmc_accept_mean": rr.ghmc_accept_mean,
                "ghmc_accept_min": rr.ghmc_accept_min,
                "wall_seconds": rr.wall_seconds,
                "ok": rr.ok,
                "reason": rr.reason,
                "flags": "; ".join(rr.flags),
            })


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Analyse a FEP run (hydration or binding) — "
                    "tarball, run directory, or raw CSV — and emit "
                    "a pass/fail markdown report + parity PNG + "
                    "normalised table. --yaml auto-infers kind "
                    "(gate, sign-rule, expected row count). Exit 0 "
                    "on gate pass OR inconclusive, 1 on real FAIL.")
    ap.add_argument(
        "path", nargs="?", default=None,
        help="freesolv_m5max_<stamp>.tar.gz, the extracted "
             "directory, or freesolv_results.csv directly. "
             "Omit to auto-discover the most recent run under "
             "run/fep/ or any freesolv_*.tar.gz at the repo root.")
    ap.add_argument(
        "--yaml", dest="yaml_path", default=None,
        help="source benchmark YAML — auto-sets --expected (entry "
             "count), --mae-gate (1.5 hydration / 2.0 binding), "
             "and the sign-critical rule (methane+acetamide for "
             "hydration; every ΔG < 0 for binding). Explicit "
             "--mae-gate still wins. Cannot pass with --expected.")
    ap.add_argument(
        "--out-dir", type=Path, default=None,
        help="write report.md + table.csv + parity.png here")
    ap.add_argument(
        "--json", action="store_true",
        help="emit the ReportResult as JSON instead of markdown")
    ap.add_argument(
        "--quiet", action="store_true",
        help="suppress markdown stdout (info/chatter); --json "
             "output always prints; exit code is always set")
    ap.add_argument(
        "--mae-gate", type=float, default=GATE_MAE_KCALMOL,
        help=f"MAE pass threshold in kcal/mol "
             f"(default {GATE_MAE_KCALMOL}; use 2.0 for binding-FEP "
             f"benchmarks per BENCHMARKS.md Milestone B)")
    ap.add_argument(
        "--expected", type=int, default=None,
        help="expected row count from the source YAML. When the "
             "CSV has fewer rows (run was killed early) the report "
             "is flagged 'partial' and cannot pass — a 6/12 run "
             "that happens to look good on its 6 is not a pass.")
    args = ap.parse_args(argv)

    # Auto-discover: if no path supplied, find the most recent run.
    if args.path is None:
        discovered = _autodiscover_latest()
        if discovered is None:
            print("fep-report: no tarball or run/fep/<stamp>/ found. "
                  "Run `cellsim fep-freesolv ...` first, or pass an "
                  "explicit path.", file=sys.stderr)
            return 2
        if not args.quiet:
            print(f"auto-discovered: {discovered}", file=sys.stderr)
        args.path = str(discovered)

    # Auto-derive --expected + --mae-gate + yaml_kind from --yaml.
    # YAML with dG_hydration_kcalmol → FreeSolv-style, gate 1.5.
    # YAML with dG_bind_kcalmol      → binding-style, gate 2.0
    # (per BENCHMARKS.md Milestone B scope — protein FEP carries
    # additional error sources so the absolute-MAE gate is looser).
    # Only auto-sets --mae-gate if the user didn't override it on
    # the CLI — their choice wins.
    user_set_mae_gate = any(
        s == "--mae-gate" or s.startswith("--mae-gate=")
        for s in (argv or sys.argv[1:]))
    yaml_kind = "hydration"  # default
    if args.yaml_path is not None:
        if args.expected is not None:
            print("fep-report: pass either --yaml or --expected, "
                  "not both", file=sys.stderr)
            return 2
        try:
            import yaml as _yaml
            data = _yaml.safe_load(Path(args.yaml_path).read_text())
            entries = (data or {}).get("entries") or []
            args.expected = len(entries)
            is_binding = any(
                "dG_bind_kcalmol" in (e or {}) for e in entries)
            is_hydration = any(
                "dG_hydration_kcalmol" in (e or {}) for e in entries)
            if is_binding and not is_hydration:
                yaml_kind = "binding"
                if not user_set_mae_gate:
                    args.mae_gate = 2.0
            elif is_hydration and not is_binding:
                yaml_kind = "hydration"
                if not user_set_mae_gate:
                    args.mae_gate = 1.5
            if not args.quiet:
                print(f"yaml {args.yaml_path}: "
                      f"{args.expected} expected rows, "
                      f"gate = {args.mae_gate} kcal/mol",
                      file=sys.stderr)
        except Exception as e:
            print(f"fep-report: cannot read --yaml {args.yaml_path}: "
                  f"{e}", file=sys.stderr)
            return 2

    try:
        r = analyse(args.path,
                    mae_gate_kcalmol=args.mae_gate,
                    expected_rows=args.expected,
                    yaml_kind=yaml_kind)
    except Exception as e:
        print(f"fep-report: failed to parse {args.path}: {e}",
              file=sys.stderr)
        return 2

    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        md = format_markdown(r)
        (args.out_dir / "report.md").write_text(md, encoding="utf-8")
        _write_table_csv(r, args.out_dir / "table.csv")
        png_ok = _write_parity_png(
            r, args.out_dir / "parity.png")
        if not args.quiet:
            print(f"wrote {args.out_dir / 'report.md'}")
            print(f"wrote {args.out_dir / 'table.csv'}")
            if png_ok:
                print(f"wrote {args.out_dir / 'parity.png'}")

    # --json emits machine-readable output even in --quiet mode;
    # --quiet was meant to silence chatter (info/progress), not
    # suppress the primary output format the caller asked for.
    if args.json:
        d = asdict(r)
        d["rows"] = [asdict(rr) for rr in r.rows]
        print(json.dumps(d, indent=2, default=str))
    elif not args.quiet:
        print(format_markdown(r))

    # Exit: 0 for pass OR inconclusive (scaffold-only / empty CSV
    # shouldn't fail CI); 1 for a real evaluable failure.
    return 0 if r.pass_overall is not False else 1


if __name__ == "__main__":
    sys.exit(main())
