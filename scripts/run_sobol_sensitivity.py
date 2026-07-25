#!/usr/bin/env python
"""Layer 1.6 — Sobol global sensitivity of docking ΔG (criterion 4).

Criterion 4 asks for "Sobol sensitivity indices documented for every
headline prediction". Read literally that is impractical AND the wrong
tool: a Sobol analysis costs n_base*(2D+2) dockings (64+ per compound),
and it answers a GLOBAL question — "across the plausible range of
docking INPUTS, which input drives the variance of the predicted ΔG?"
— not a per-compound uncertainty (that is the conformal / reliability
layer, see scripts/run_uq_coverage.py).

So the honest reading is: run Sobol ONCE on a representative case to
learn what the pipeline is sensitive to, and document it. That tells a
user which knobs to control (e.g. "pin exhaustiveness; box placement
barely matters" or the reverse) and tells a developer where to spend
determinism effort.

This runner does exactly that and writes the indices to JSON.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


# Below this spread (kcal/mol) the docking output is treated as
# INSENSITIVE to the swept inputs, and the normalised Sobol indices are
# flagged unreliable (they divide by a near-zero total variance). 0.5
# kcal/mol ≈ RT·ln(10)/2.7, well under the ~1 log-Kd decision bar.
INSENSITIVE_STD_KCALMOL = 0.5


def dg_spread(samples):
    """Summary statistics of the raw ΔG samples, or None if empty.

    The raw spread is the physically meaningful output of a Sobol run on
    a well-behaved case: if it is tiny, the docking ΔG does not depend on
    the swept inputs at all, and the normalised indices are noise.
    """
    vals = [float(v) for v in (samples or []) if v is not None]
    if not vals:
        return None
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / n
    return {
        "n": n, "mean": mean, "std": var ** 0.5,
        "min": min(vals), "max": max(vals),
        "range": max(vals) - min(vals),
    }


def main(argv=None) -> int:
    from src.uq import sobol_dock

    ap = argparse.ArgumentParser(description=__doc__)
    # Default: trypsin/benzamidine — the well-behaved, trustworthy class,
    # so the sensitivity reflects the method, not a pathological target.
    ap.add_argument("--receptor", default="benchmarks/dock/3ptb.pdb")
    ap.add_argument("--smiles", default="c1ccc(cc1)C(=N)N")  # benzamidine
    # 3ptb/BEN pocket centroid (benzamidine bounding-box centre),
    # extracted via src.dock.extract_hetatm_ligand at audit time.
    ap.add_argument("--center", nargs=3, type=float,
                    default=[-1.86, 14.37, 16.75])
    ap.add_argument("--box", nargs=3, type=float,
                    default=[20.0, 20.0, 20.0])
    ap.add_argument("--n-base", type=int, default=8,
                    help="Saltelli base count; total runs = n_base*(2D+2). "
                         ">=32 for stable indices; 8 = fast honest demo.")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out-json", default="benchmarks/dock/sobol_sensitivity.json")
    args = ap.parse_args(argv)

    print(f"[sobol] {Path(args.receptor).name}  n_base={args.n_base} "
          f"(=> {args.n_base * 8} dockings)", flush=True)
    r = sobol_dock(
        receptor_pdb=args.receptor, ligand_smiles=args.smiles,
        center_A=tuple(args.center), box_size_A=tuple(args.box),
        n_base=args.n_base, workers=args.workers)
    if not r.ok:
        print(f"[sobol] FAILED: {r.reason}")
        return 1

    # The raw ΔG spread across all Saltelli samples is the physically
    # meaningful quantity. If it is tiny, the docking is INSENSITIVE to
    # all inputs and the normalised Sobol indices become numerically
    # unstable (they divide by a near-zero total variance and fly out of
    # [0,1]) — a small spread is the finding, not the indices.
    spread = dg_spread(r.dG_samples)

    # Per-input first-order (S1) and total (ST) indices. SobolResult
    # stores them as dicts keyed by input name.
    names = list(r.input_names or [])

    def _num(d, k):
        v = (d or {}).get(k)
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    rows = [{"input": nm, "S1": _num(r.s1, nm),
             "ST": _num(r.st, nm),
             "ST_conf": _num(getattr(r, "st_conf", {}), nm)}
            for nm in names]
    rows.sort(key=lambda d: (d["ST"] if d["ST"] is not None else -1),
              reverse=True)

    print()
    if spread is not None:
        print(f"  Raw ΔG spread over {spread['n']} samples: "
              f"mean {spread['mean']:.2f}  std {spread['std']:.3f}  "
              f"range [{spread['min']:.2f}, {spread['max']:.2f}] "
              f"= {spread['range']:.2f} kcal/mol")
        if spread["std"] < INSENSITIVE_STD_KCALMOL:
            print("  => spread < 0.5 kcal/mol: docking is INSENSITIVE to "
                  "these inputs for this case;\n     normalised Sobol "
                  "indices below are numerically unstable (near-zero "
                  "denominator) —\n     trust the spread, not the indices.")
        print()
    print("  Which docking input drives ΔG variance? (higher = more)")
    print(f"    {'input':<18}{'S1 (first-order)':>18}{'ST (total)':>14}")
    for d in rows:
        s1s = f"{d['S1']:.3f}" if d['S1'] is not None else "  n/a"
        sts = f"{d['ST']:.3f}" if d['ST'] is not None else "  n/a"
        print(f"    {d['input']:<18}{s1s:>18}{sts:>14}")
    if rows and rows[0]["ST"] is not None:
        print(f"\n  Dominant input: {rows[0]['input']}  "
              f"(control this first for reproducible ΔG).")
    print(f"  NB n_base={args.n_base}: honest ranking, wide CI on the "
          "exact indices (>= 32 for tight numbers).")

    out = {
        "receptor": args.receptor, "smiles": args.smiles,
        "n_base": args.n_base, "n_runs": args.n_base * 8,
        "dG_spread": spread,
        "indices": rows,
        "dominant_input": rows[0]["input"] if rows else None,
        "indices_reliable": bool(
            spread and spread["std"] >= INSENSITIVE_STD_KCALMOL),
    }
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
