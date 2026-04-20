"""Triage-breakdown PNG for a batch docking run.

Given the records list from `src.dock.batch.run_batch` (or a CSV),
render a single-page dashboard that a wet-lab user can paste into
a meeting:

  1. Verdict breakdown (horizontal stacked bar):
     follow_up | review | deprioritise | drop
     Widths are counts. Coloured by the same palette used in the
     stdout rank table.

  2. ΔG histogram, faceted by verdict: shows where the follow-up
     compounds sit on the affinity axis versus the drops.

  3. Pose-trust column (stacked) — good / acceptable / suspicious /
     reject — so you can see at a glance how much of the hit list
     you should actually trust.

  4. Safety filter pass rates: Ro5, hERG, Ames. Bar chart with
     count of compounds that PASS each. Reveals which filter is
     gating the campaign.

Non-AI; just matplotlib over the batch CSV.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]

TRIAGE_ORDER = ["follow_up", "review", "deprioritise", "drop"]
TRIAGE_COLOURS = {
    "follow_up":    "#2b8a3e",   # green
    "review":       "#f0b429",   # amber
    "deprioritise": "#868e96",   # grey
    "drop":         "#c92a2a",   # red
}
STRAIN_ORDER = ["good", "acceptable", "suspicious", "reject"]
STRAIN_COLOURS = {
    "good":       "#2b8a3e",
    "acceptable": "#74c0fc",
    "suspicious": "#f0b429",
    "reject":     "#c92a2a",
}


def _load_records(csv_path: Path) -> list[dict]:
    records: list[dict] = []
    with csv_path.open() as fp:
        for row in csv.DictReader(fp):
            records.append(row)
    return records


def _parse_float(s) -> Optional[float]:
    if s in (None, "", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def render_triage_dashboard(
    records: list[dict],
    save: Path,
    *,
    title: str = "CellSim batch screen — triage breakdown",
) -> None:
    import matplotlib.pyplot as plt

    # Count verdicts (defensive against missing columns in older CSVs).
    verdict_counts = {v: 0 for v in TRIAGE_ORDER}
    strain_counts = {s: 0 for s in STRAIN_ORDER}
    per_verdict_dGs: dict[str, list[float]] = {
        v: [] for v in TRIAGE_ORDER}

    n_total = 0
    n_ok = 0
    ro5_pass = 0
    herg_pass = 0      # low + medium risk counted as pass
    mut_pass = 0       # same convention
    ro5_known = 0
    herg_known = 0
    mut_known = 0

    for r in records:
        n_total += 1
        if str(r.get("ok", "")).lower() != "true":
            continue
        n_ok += 1
        v = r.get("triage") or "drop"
        if v not in verdict_counts:
            v = "drop"
        verdict_counts[v] += 1
        s = r.get("strain_band") or ""
        if s in strain_counts:
            strain_counts[s] += 1
        dG = _parse_float(r.get("dG_kcalmol"))
        if dG is not None:
            per_verdict_dGs[v].append(dG)
        # Safety filters.
        ro5 = r.get("ro5_pass")
        if ro5 not in (None, "", "None"):
            ro5_known += 1
            if str(ro5).lower() == "true":
                ro5_pass += 1
        herg = r.get("herg_risk")
        if herg:
            herg_known += 1
            if herg != "high":
                herg_pass += 1
        mut = r.get("mutagenic_risk")
        if mut:
            mut_known += 1
            if mut != "high":
                mut_pass += 1

    # Figure layout: 2 x 2 grid.
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5))
    fig.suptitle(
        f"{title}\n"
        f"{n_ok} of {n_total} compounds docked OK",
        fontsize=13)

    # --- (1) Verdict stacked horizontal bar ---------------------
    ax = axes[0, 0]
    left = 0
    for v in TRIAGE_ORDER:
        c = verdict_counts[v]
        if c == 0:
            continue
        ax.barh(0, c, left=left, color=TRIAGE_COLOURS[v],
                edgecolor="white", linewidth=1.0, label=v)
        ax.text(left + c / 2, 0, f"{v}\n{c}",
                ha="center", va="center", fontsize=9,
                color="white", fontweight="bold")
        left += c
    ax.set_yticks([])
    ax.set_xlim(0, max(n_ok, 1))
    ax.set_title("Verdict breakdown (stacked, n-compounds)",
                  fontsize=11)
    ax.set_xlabel("count")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    # --- (2) ΔG facet scatter by verdict ------------------------
    ax = axes[0, 1]
    y_off = 0
    ytick_pos = []
    ytick_lbl = []
    for v in TRIAGE_ORDER:
        dGs = per_verdict_dGs[v]
        if not dGs:
            continue
        ax.scatter(dGs, [y_off] * len(dGs),
                   color=TRIAGE_COLOURS[v], s=40,
                   edgecolor="white", linewidth=0.7, alpha=0.9)
        ytick_pos.append(y_off)
        ytick_lbl.append(f"{v}  (n={len(dGs)})")
        y_off += 1
    ax.set_yticks(ytick_pos)
    ax.set_yticklabels(ytick_lbl)
    ax.invert_yaxis()
    ax.axvline(-6.0, color="#c92a2a", linestyle="--", alpha=0.4,
                linewidth=1)
    ax.axvline(-7.5, color="#f0b429", linestyle="--", alpha=0.4,
                linewidth=1)
    ax.set_xlabel("ΔG (kcal/mol)  — lower is tighter")
    ax.set_title("ΔG by verdict", fontsize=11)
    ax.grid(axis="x", linestyle=":", alpha=0.4)

    # --- (3) Strain-band bar ------------------------------------
    ax = axes[1, 0]
    strain_total = sum(strain_counts.values())
    if strain_total == 0:
        ax.text(0.5, 0.5, "strain disabled / no strain data",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="#868e96")
        ax.axis("off")
    else:
        bars = [strain_counts[s] for s in STRAIN_ORDER]
        colors = [STRAIN_COLOURS[s] for s in STRAIN_ORDER]
        ax.bar(STRAIN_ORDER, bars, color=colors, edgecolor="white")
        for i, b in enumerate(bars):
            if b == 0:
                continue
            ax.text(i, b, str(b), ha="center",
                    va="bottom", fontsize=9, fontweight="bold")
        ax.set_ylabel("count")
        ax.set_title(
            "Pose-trust (UFF-ensemble strain band)", fontsize=11)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    # --- (4) Safety pass-rates ---------------------------------
    ax = axes[1, 1]
    labels = ["Ro5", "hERG", "Ames"]
    passes = [ro5_pass, herg_pass, mut_pass]
    known = [ro5_known, herg_known, mut_known]
    xs = range(len(labels))
    bar_known = ax.bar(xs, known, color="#dee2e6",
                        edgecolor="white", label="compounds evaluated")
    bar_pass = ax.bar(xs, passes, color="#2b8a3e",
                       edgecolor="white", label="pass")
    for i, (p, k) in enumerate(zip(passes, known)):
        if k == 0:
            continue
        ax.text(i, k, f"{p}/{k}\n({100 * p / k:.0f}%)",
                ha="center", va="bottom", fontsize=9,
                fontweight="bold")
    ax.set_xticks(list(xs))
    ax.set_xticklabels(labels)
    ax.set_ylabel("count")
    ax.set_title("Safety-filter pass rates", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path,
                    help="batch output CSV (from `cellsim dock`)")
    ap.add_argument("--out", type=Path, required=True,
                    help="output PNG path")
    ap.add_argument("--title", default="CellSim batch screen — "
                                         "triage breakdown")
    args = ap.parse_args(argv)

    records = _load_records(args.csv)
    if not records:
        print(f"no records in {args.csv}", file=sys.stderr)
        return 2
    render_triage_dashboard(records, args.out, title=args.title)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
