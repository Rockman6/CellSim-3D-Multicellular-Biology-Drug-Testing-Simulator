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

What is deliberately NOT gated: a biological sequence effect. The engine
does not produce one once the readout is long enough, and a test
asserting otherwise would be pinning an artefact.
"""
from __future__ import annotations

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
CISPLATIN, PACLITAXEL = get_drug("cisplatin"), get_drug("paclitaxel")
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


def test_the_apparent_sequence_effect_is_a_readout_artefact():
    """Giving the slow drug second truncates its effect inside a fixed
    window. Lengthening the readout must shrink the ordering difference
    toward 1 — so a sequence claim read off a single endpoint would be
    an artefact of when the assay stopped."""
    def cis_first(t):
        return (CIS, 0.0) if t < HALF else ((0.0, PAC) if t < FULL else (0.0, 0.0))

    def pac_first(t):
        return (0.0, PAC) if t < HALF else ((CIS, 0.0) if t < FULL else (0.0, 0.0))

    short = _surviving(PAIR, cis_first, 72.0) / _surviving(PAIR, pac_first, 72.0)
    long = _surviving(PAIR, cis_first, 120.0) / _surviving(PAIR, pac_first, 120.0)
    assert short < 0.98, f"expected an apparent ordering effect at 72 h, got {short:.2f}"
    assert abs(long - 1.0) < abs(short - 1.0), (
        f"the effect should shrink with a longer readout: {short:.2f} -> {long:.2f}")


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
