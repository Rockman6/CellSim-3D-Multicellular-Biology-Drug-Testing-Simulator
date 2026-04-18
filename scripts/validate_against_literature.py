#!/usr/bin/env python3
"""
CellSim validation against published HeLa cell data.

Reads:
  logs/export/cellsim_timeseries.csv              (simulation output)
  data/reference/hela_reference.csv               (literature ranges)
  data/reference/growth_curves/*.csv              (experimental datasets)

Prints pass/fail for each parameter.
Writes a comparison summary to logs/export/validation_report.txt
"""
import csv, os, sys, math

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)
SIM_CSV  = os.path.join(ROOT, "logs/export/cellsim_timeseries.csv")
REF_CSV  = os.path.join(ROOT, "data/reference/hela_reference.csv")
GC_DIR   = os.path.join(ROOT, "data/reference/growth_curves")
REPORT   = os.path.join(ROOT, "logs/export/validation_report.txt")

GREEN = "\033[92m"
RED   = "\033[91m"
YELLOW= "\033[93m"
RESET = "\033[0m"

def load_reference():
    ref = {}
    with open(REF_CSV) as f:
        lines = [l for l in f if not l.strip().startswith('#') and l.strip()]
    reader = csv.DictReader(lines)
    for row in reader:
        name = row.get('parameter','').strip()
        if not name or name == 'parameter': continue
        try:
            ref[name] = {
                'value': float(row['value']) if row.get('value') else None,
                'min': float(row['min']) if row.get('min') else None,
                'max': float(row['max']) if row.get('max') else None,
                'unit': row.get('unit',''), 'source': row.get('source',''),
            }
        except (ValueError, KeyError): pass
    return ref

def load_sim():
    rows = []
    with open(SIM_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader: rows.append(row)
    if not rows:
        print("ERROR: No simulation data. Run CellSim first."); sys.exit(1)
    return rows

def load_growth_ref():
    """Load the single-cell-origin expected growth curve."""
    path = os.path.join(GC_DIR, "hela_single_cell_origin.csv")
    pts = []
    with open(path) as f:
        lines = [l for l in f if not l.strip().startswith('#') and l.strip()]
    reader = csv.DictReader(lines)
    for row in reader:
        pts.append({
            'hours': float(row['hours']),
            'expected': float(row['expected_cells']),
            'min': float(row['min_cells']),
            'max': float(row['max_cells']),
        })
    return pts

def load_metabolomics():
    path = os.path.join(GC_DIR, "hela_metabolomics.csv")
    met = {}
    with open(path) as f:
        lines = [l for l in f if not l.strip().startswith('#') and l.strip()]
    reader = csv.DictReader(lines)
    for row in reader:
        name = row.get('metabolite','').strip()
        if not name: continue
        try:
            met[name] = {
                'conc': float(row['concentration_mM']),
                'std': float(row.get('std_mM',0)),
                'source': row.get('source',''),
            }
        except (ValueError, KeyError): pass
    return met

def check(name, sim_val, lo, hi, unit, source, out_lines):
    ok = lo <= sim_val <= hi
    tag = f"{GREEN}PASS{RESET}" if ok else f"{RED}FAIL{RESET}"
    line = f"  {tag}  {name:35s} = {sim_val:12.4f}  [{lo:.4f} - {hi:.4f}] {unit:6s}  {source}"
    print(line)
    plain = f"  {'PASS' if ok else 'FAIL'}  {name:35s} = {sim_val:12.4f}  [{lo:.4f} - {hi:.4f}] {unit:6s}  {source}"
    out_lines.append(plain)
    return ok

def main():
    ref = load_reference()
    sim_rows = load_sim()
    last = sim_rows[-1]
    gc_ref = load_growth_ref()
    met_ref = load_metabolomics()

    bio_h = float(last['bio_hours'])
    n_cells = int(last['cell_count'])
    out = []
    n_pass, n_fail, n_skip = 0, 0, 0

    header = (f"\n{'='*80}\n  CellSim Validation Report\n"
              f"  Bio-time: {bio_h:.1f}h  |  Wall: {float(last['wall_sec']):.0f}s  |  Cells: {n_cells}\n"
              f"{'='*80}")
    print(header); out.append(header)

    # ── Metabolite concentrations vs Park 2016 ──
    sec = "\n--- Metabolite Concentrations (vs Park 2016 LC-MS) ---"
    print(sec); out.append(sec)
    met_map = {
        'ATP': 'ATP_mM', 'ADP': 'ADP_mM', 'AMP': 'AMP_mM',
        'pyruvate': 'pyruvate_mM', 'lactate': 'lactate_mM', 'citrate': 'citrate_mM',
        'glucose-6-phosphate': None, 'NAD+': 'NAD_plus_mM', 'NADH': 'NADH_mM',
    }
    for met_name, m in met_ref.items():
        sim_col = met_map.get(met_name)
        if sim_col and sim_col in last:
            lo = m['conc'] - 2*m['std']
            hi = m['conc'] + 2*m['std']
            ok = check(met_name, float(last[sim_col]), max(0,lo), hi, 'mM', m['source'], out)
            if ok: n_pass += 1
            else: n_fail += 1

    # ── Derived ratios ──
    sec = "\n--- Derived Ratios ---"
    print(sec); out.append(sec)
    r = ref
    for name, col in [('energy_charge','energy_charge'), ('NAD_ratio','NAD_ratio')]:
        if name in r and col in last:
            ok = check(name, float(last[col]), r[name]['min'], r[name]['max'], r[name]['unit'], r[name]['source'], out)
            if ok: n_pass += 1
            else: n_fail += 1

    # ── Membrane ──
    sec = "\n--- Membrane Potential ---"
    print(sec); out.append(sec)
    if 'membrane_potential' in r and 'membrane_potential_mV' in last:
        ok = check('membrane_potential', float(last['membrane_potential_mV']),
                    r['membrane_potential']['min'], r['membrane_potential']['max'],
                    'mV', r['membrane_potential']['source'], out)
        if ok: n_pass += 1
        else: n_fail += 1

    # ── Growth curve (single cell origin) ──
    sec = "\n--- Growth Curve (vs single-cell-origin reference) ---"
    print(sec); out.append(sec)
    for pt in gc_ref:
        th = pt['hours']
        # Find closest sim row
        best = None; best_diff = 1e9
        for row in sim_rows:
            h = float(row['bio_hours'])
            if abs(h - th) < best_diff:
                best_diff = abs(h - th)
                best = row
        if best and best_diff < 5:
            sim_n = float(best['cell_count'])
            ok = check(f"cells_at_{th:.0f}h", sim_n, pt['min'], pt['max'],
                       'cells', 'Chao2019/ATCC', out)
            if ok: n_pass += 1
            else: n_fail += 1
        elif th <= bio_h + 5:
            line = f"  {YELLOW}SKIP{RESET}  cells_at_{th:.0f}h (no matching sim timepoint)"
            print(line); out.append(line); n_skip += 1

    # ── Telomere ──
    sec = "\n--- Telomere ---"
    print(sec); out.append(sec)
    if 'telomere_bp' in last:
        gen = int(float(last.get('generation', 0)))
        telo = float(last['telomere_bp'])
        expected_telo = 10000 - gen * 100
        line = f"  INFO  telomere={telo:.0f} bp  gen={gen}  expected={expected_telo:.0f} bp (100 bp/div)"
        print(line); out.append(line)
        if gen > 0:
            loss_per_div = (10000 - telo) / gen
            ok = check('telomere_loss_per_div', loss_per_div, 50, 200, 'bp', 'Harley 1990', out)
            if ok: n_pass += 1
            else: n_fail += 1

    # ── Phase distribution ──
    sec = "\n--- Cell Cycle Phase Distribution ---"
    print(sec); out.append(sec)
    if n_cells >= 4:  # need enough cells for meaningful fractions
        fG1 = int(last['phase_G1']) / n_cells
        fS  = int(last['phase_S'])  / n_cells
        fG2 = int(last['phase_G2']) / n_cells
        fM  = int(last['phase_M'])  / n_cells
        for name, val in [('G1_fraction',fG1),('S_fraction',fS),('G2_fraction',fG2),('M_fraction',fM)]:
            if name in r:
                ok = check(name, val, r[name]['min'], r[name]['max'], 'frac', r[name]['source'], out)
                if ok: n_pass += 1
                else: n_fail += 1
    else:
        line = f"  {YELLOW}SKIP{RESET}  phase distribution (need >=4 cells, have {n_cells})"
        print(line); out.append(line); n_skip += 1

    # ── Fluxes ──
    sec = "\n--- Metabolic Fluxes ---"
    print(sec); out.append(sec)
    for name, col in [('glycolysis_flux','glycolysis_flux'),('TCA_flux','TCA_flux'),('OxPhos_ATP_rate','OxPhos_flux')]:
        if name in r and col in last:
            ok = check(name, float(last[col]), r[name]['min'], r[name]['max'], 'mM/s', r[name]['source'], out)
            if ok: n_pass += 1
            else: n_fail += 1

    # ── Summary ──
    total = n_pass + n_fail
    pct = 100 * n_pass / max(total, 1)
    summary = (f"\n{'='*80}\n"
               f"  RESULTS: {n_pass} PASS / {n_fail} FAIL / {n_skip} SKIP  ({pct:.0f}% pass rate)\n"
               f"  Doubling time target: 22h  |  Expected cells at {bio_h:.0f}h: ~{2**(bio_h/22):.0f}  |  Actual: {n_cells}\n"
               f"{'='*80}\n")
    print(summary); out.append(summary)

    # Write report
    with open(REPORT, 'w') as f:
        f.write('\n'.join(out))
    print(f"  Report saved to: {REPORT}")

if __name__ == "__main__":
    main()
