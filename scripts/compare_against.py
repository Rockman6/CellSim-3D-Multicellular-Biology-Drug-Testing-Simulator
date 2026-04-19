#!/usr/bin/env python3
"""Compare a cellsim_headless output CSV against one reference dataset.

Usage:
    ./scripts/compare_against.py <sim_csv> <ref_csv> <ref_schema>

ref_schema ∈ {'viability', 'growth_curve', 'single_cell', 'ctc_cellcount'}

Outputs a per-timepoint table (sim vs ref at reference's time grid) and a
summary row: mean |relative error|, mean |absolute error|, doubling-time
comparison.
"""
import csv
import math
import sys


def read_sim(path):
    """Return [(bio_hours, cell_count), ...]"""
    out = []
    with open(path) as f:
        r = csv.reader(f)
        header = next(r)
        # bio_hours is col 2, cell_count is col 3 (headless_doubling.csv schema)
        h_i = header.index("bio_hours") if "bio_hours" in header else 2
        c_i = header.index("cell_count") if "cell_count" in header else 3
        for row in r:
            if not row or row[0].startswith("#"):
                continue
            try:
                out.append((float(row[h_i]), int(row[c_i])))
            except (ValueError, IndexError):
                continue
    return out


def read_ref(path, schema):
    """Return [(hours, ref_value), ...] per schema."""
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(("frame,", "hours,", "study,", "condition,",
                                "metabolite,", "cell_type,")):
                continue
            parts = line.split(",")
            try:
                if schema == "ctc_cellcount":
                    # frame,hours,cell_count,divisions_this_frame
                    out.append((float(parts[1]), int(parts[2])))
                elif schema == "viability":
                    # hours,relative_proliferation,std_dev,seeded_cells,notes
                    out.append((float(parts[0]), float(parts[1])))
                elif schema == "growth_curve":
                    # hours,cells,phase,notes
                    out.append((float(parts[0]), int(parts[1])))
                elif schema == "single_cell":
                    # hours,expected_cells,min_cells,max_cells,divisions,notes
                    out.append((float(parts[0]), int(parts[1])))
                else:
                    raise ValueError(f"unknown schema: {schema}")
            except (ValueError, IndexError):
                continue
    return out


def interpolate(pts, x):
    """Piecewise-linear interpolation on a sorted list of (x,y)."""
    if not pts:
        return 0.0
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]
    for i in range(1, len(pts)):
        if pts[i][0] >= x:
            t0, y0 = pts[i - 1]
            t1, y1 = pts[i]
            u = (x - t0) / (t1 - t0 + 1e-9)
            return y0 + (y1 - y0) * u
    return pts[-1][1]


def doubling_time(pts, h0=6.0, h1=None):
    """Estimate doubling time from log-phase slope between h0 and h1."""
    if len(pts) < 2:
        return None
    if h1 is None:
        # Use the interval before any plateau — first 1/3 of the curve.
        h1 = max(h0 + 1, pts[0][0] + (pts[-1][0] - pts[0][0]) * 0.33)
    y0 = interpolate(pts, h0)
    y1 = interpolate(pts, h1)
    if y0 <= 0 or y1 <= y0:
        return None
    n_doubles = math.log(y1 / y0) / math.log(2)
    if n_doubles <= 0:
        return None
    return (h1 - h0) / n_doubles


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    sim_csv, ref_csv, schema = sys.argv[1], sys.argv[2], sys.argv[3]

    sim = read_sim(sim_csv)
    ref = read_ref(ref_csv, schema)
    if not sim or not ref:
        print(f"EMPTY DATA  sim={len(sim)}  ref={len(ref)}")
        return

    # Normalise so both are comparable:
    # - ctc_cellcount / growth_curve / single_cell: absolute counts
    # - viability: relative fold-change (sim/init == ref)
    if schema == "viability":
        sim_init = sim[0][1] if sim[0][1] > 0 else 1
        sim = [(h, c / sim_init) for h, c in sim]

    print(f"{'hours':>8}  {'sim':>10}  {'ref':>10}  {'abs_err':>10}  {'rel_err_pct':>12}")
    abs_errs, rel_errs = [], []
    for h, ry in ref:
        sy = interpolate(sim, h)
        abs_err = sy - ry
        rel_err = (abs_err / ry * 100.0) if ry else 0.0
        abs_errs.append(abs(abs_err))
        rel_errs.append(abs(rel_err))
        print(f"{h:>8.1f}  {sy:>10.3f}  {ry:>10.3f}  {abs_err:>+10.3f}  {rel_err:>+12.2f}")

    mean_abs = sum(abs_errs) / len(abs_errs)
    mean_rel = sum(rel_errs) / len(rel_errs)
    max_rel = max(rel_errs)
    dt_sim = doubling_time(sim)
    dt_ref = doubling_time(ref)
    print("-" * 60)
    print(f"mean |abs_err|     = {mean_abs:.3f}")
    print(f"mean |rel_err|     = {mean_rel:.2f} %")
    print(f"max  |rel_err|     = {max_rel:.2f} %")
    if dt_sim and dt_ref:
        print(f"doubling_time sim  = {dt_sim:.2f} h")
        print(f"doubling_time ref  = {dt_ref:.2f} h")
        print(f"doubling_time Δ    = {dt_sim - dt_ref:+.2f} h "
              f"({(dt_sim - dt_ref) / dt_ref * 100:+.1f} %)")


if __name__ == "__main__":
    main()
