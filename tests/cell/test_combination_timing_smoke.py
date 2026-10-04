#!/usr/bin/env python3
"""Phase-2 gates: two drugs at once, and what the engine does and does not show.

Within-line, so the per-line prediction ceiling in docs/VALIDATION.md
does not apply: every arm uses the same two drugs at the same doses in
the same cell line, and only the timing differs.

What is gated:

* **Antagonism**, which is robust. Cisplatin damages DNA, driving
  p53 -> p21 and arresting the cycle; paclitaxel kills only cells that
  attempt mitosis while it is present. A cytostatic removes the cells a
  phase-specific partner needs, so the pair must kill LESS than
  independent action predicts.
* **The timing confound**, also robust and worth pinning. The apparent
  advantage of one order over the other shrinks toward 1 as the readout
  lengthens, because cisplatin acts slowly and placing it second
  truncates its effect inside a fixed window.

What is deliberately NOT gated: a biological sequence effect. With both
drugs at their IC50 the engine produces none, and above that the one it
produces mostly fades with a longer readout; the gates pin exactly that.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, calibrate_cycle_scale, simulate  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

LINE = get_line("A549")
P = Params()
K = calibrate_cycle_scale(LINE, P)
N_CELLS = 64
HALF, FULL = 36.0, 72.0
CIS, PAC = 15.0, 0.04
# Held at the potency the doses were designed for (the fit before the
# October 2026 refit); both drugs enter the engine only as potency x
# concentration, so this is exactly a rescaling of CIS and PAC.
CISPLATIN = dataclasses.replace(get_drug("cisplatin"), k_damage_per_uM_h=0.000951)
PACLITAXEL = dataclasses.replace(get_drug("paclitaxel"), partition=0.153)
PAIR = [CISPLATIN, PACLITAXEL]


def _surviving(drugs, dose_fn, total=FULL):
    treated = simulate(LINE, drugs, 0.0, t_end_h=total, n_cells=N_CELLS, k_cyc=K, p=P,
                       dose_fn=dose_fn, record_every_h=total)
    control = simulate(LINE, None, 0.0, t_end_h=total, n_cells=N_CELLS, k_cyc=K, p=P,
                       record_every_h=total)
    return float(treated.alive_weight[-1, 0] / max(control.alive_weight[-1, 0], 1e-12))


def test_engine_carries_two_drugs_independently():
    """Each drug needs its own intracellular concentration; a combination
    must not collapse into a single shared pool."""
    res = simulate(LINE, PAIR, [CIS, PAC], t_end_h=12.0, n_cells=8, k_cyc=K, p=P)
    from cellsim.cell.engine import N_BASE
    assert res.final.Y.shape[1] == N_BASE + 2, res.final.Y.shape
    cis_in, pac_in = res.final.Y[:, N_BASE], res.final.Y[:, N_BASE + 1]
    assert cis_in.max() > 0 and pac_in.max() > 0, "both drugs should enter cells"
    assert cis_in.max() > 100 * pac_in.max(), (
        "the two pools should track their own very different doses, not be shared")


def test_a_cytostatic_antagonises_a_phase_specific_partner():
    """The pair must kill LESS than independent action predicts, because
    the platinum arrests the cells the taxane needs in mitosis."""
    solo_cis = _surviving([CISPLATIN], lambda t: (CIS if t < HALF else 0.0,))
    solo_pac = _surviving([PACLITAXEL], lambda t: (PAC if t < HALF else 0.0,))
    together = _surviving(PAIR, lambda t: (CIS, PAC) if t < HALF else (0.0, 0.0))
    independent = solo_cis * solo_pac
    assert together > independent * 1.1, (
        f"expected antagonism: combination {together:.3f} vs independent "
        f"expectation {independent:.3f}")
    assert together < min(solo_cis, solo_pac), (
        "the combination should still beat either drug alone")


# Engine 72 h IC50 at the design potency above (cellsim.cell.engine.ic50).
IC50_CIS, IC50_PAC = 13.24, 0.0418
SEEDS = (1, 2, 3)


def _ordering_ratio(cis_uM: float, pac_uM: float) -> tuple[float, float]:
    """Mean over seeds of cisplatin-first / paclitaxel-first survival at a
    72 h and a 120 h readout. One 128-cell run carries +-0.02-0.03 of
    noise on this ratio — the size of the effect itself at matched
    doses — so it is only ever read as an average over seeds."""
    def cis_first(t):
        return (cis_uM, 0.0) if t < HALF else ((0.0, pac_uM) if t < FULL else (0.0, 0.0))

    def pac_first(t):
        return (0.0, pac_uM) if t < HALF else ((cis_uM, 0.0) if t < FULL else (0.0, 0.0))
    r72, r120 = [], []
    for seed in SEEDS:
        a = simulate(LINE, PAIR, 0.0, t_end_h=120.0, n_cells=128, k_cyc=K, p=P,
                     dose_fn=cis_first, record_every_h=24.0, seed=seed).alive_weight[:, 0]
        b = simulate(LINE, PAIR, 0.0, t_end_h=120.0, n_cells=128, k_cyc=K, p=P,
                     dose_fn=pac_first, record_every_h=24.0, seed=seed).alive_weight[:, 0]
        r72.append(a[3] / b[3])
        r120.append(a[5] / b[5])
    return float(sum(r72) / len(r72)), float(sum(r120) / len(r120))


def test_no_ordering_effect_with_both_drugs_at_their_ic50():
    """The engine must not invent a sequence effect: with each drug at its
    own IC50 the two orders leave the same colony at any readout."""
    r72, r120 = _ordering_ratio(IC50_CIS, IC50_PAC)
    assert abs(r72 - 1.0) < 0.05 and abs(r120 - 1.0) < 0.05, (
        f"ordering ratio {r72:.3f} at 72 h, {r120:.3f} at 120 h; expected ~1")


def test_an_ordering_effect_above_ic50_fades_with_a_longer_readout():
    """With cisplatin above its IC50 the platinum-first order looks better
    at 72 h, but the gap shrinks once the slow drug given second has had
    time to act — so a sequence claim read off one fixed endpoint would
    mostly be measuring when the assay stopped."""
    r72, r120 = _ordering_ratio(1.25 * IC50_CIS, IC50_PAC)
    assert r72 < 0.95, f"expected an apparent ordering effect at 72 h, got {r72:.3f}"
    assert r120 - r72 > 0.03, (
        f"the effect should fade with a longer readout: {r72:.3f} -> {r120:.3f}")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  ok   {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL {fn.__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} combination gates pass")
    raise SystemExit(1 if failed else 0)
