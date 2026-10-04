#!/usr/bin/env python3
"""Phase-3 gates: an MDM2 inhibitor, where p53 status IS the mechanism.

Nutlin-3a occupies MDM2's p53-binding pocket, so p53 accumulates without
any DNA damage (Vassilev et al. 2004 Science 303:844). Whether that kills
a cell depends entirely on whether functional p53 is there to accumulate,
which makes this the one drug in the panel whose line-to-line differences
the engine can derive from mechanism rather than from a fitted marker.

Three cases, and the third is the one worth having:

* **TP53 wild-type** — p53 rises, PUMA follows, the cell dies.
* **TP53 mutant** — transactivation is gone, so a stabilised p53 does
  nothing. Blocking MDM2 harder cannot help.
* **HeLa** — wild-type TP53 *by sequence*, but HPV18 E6 hands p53 to the
  E6AP ligase instead of MDM2, and in cervical cancer cells that switch
  is complete (Hengstermann et al. 2001 PNAS 98:1218). So an MDM2
  inhibitor should fail there too, and it does in GDSC. A model that read
  only the genotype would get this line wrong.

Nothing here is fitted to the p53 split: `nutlin-3a`'s one fitted
constant is its partition coefficient, fitted on two wild-type lines
only. The mutant and HeLa predictions are consequences.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import Params, dose_response  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402

NUTLIN = get_drug("nutlin-3a")
# GDSC screened nutlin-3a up to 10 uM; 30 uM is beyond any tested dose,
# so "still alive at 30 uM" means resistant by any reasonable reading.
CONC_UM = np.array([1.0, 3.0, 10.0, 30.0])
N_CELLS = 32


def _viability(line_name: str) -> np.ndarray:
    viab, _ = dose_response(get_line(line_name), NUTLIN, CONC_UM, t_end_h=72.0,
                            n_cells_per_conc=N_CELLS, seed=1, p=Params())
    return viab[1:]


def test_mdm2_inhibition_kills_p53_wild_type_lines():
    for name in ("A549", "HCT116"):
        v = _viability(name)
        assert v[-1] < 0.5, f"{name} (TP53 wild-type) should die to nutlin: {np.round(v, 3)}"
        assert v[0] >= v[-1], f"{name}: viability should fall with dose, got {np.round(v, 3)}"


def test_mdm2_inhibition_does_nothing_without_functional_p53():
    """A stabilised mutant p53 cannot transactivate PUMA, so no dose helps."""
    v = _viability("HT-29")
    assert v[-1] > 0.8, f"TP53-mutant HT-29 should resist nutlin: {np.round(v, 3)}"


def test_hela_resists_although_its_TP53_is_wild_type():
    """The prediction that genotype alone would get wrong: HPV E6 degrades
    p53 independently of MDM2, so blocking MDM2 cannot rescue it."""
    hela = get_line("HeLa")
    assert hela.p53_functional, "HeLa is TP53 wild-type by sequence"
    assert hela.p53_mdm2_independent_deg_per_h > 0, "the E6 route must be modelled"
    v = _viability("HeLa")
    assert v[-1] > 0.8, f"HeLa should resist nutlin despite wild-type TP53: {np.round(v, 3)}"


def test_removing_the_e6_route_restores_sensitivity():
    """Specificity: the resistance must come from the E6 route and not
    from something else about the line."""
    import dataclasses
    hela = get_line("HeLa")
    no_e6 = dataclasses.replace(hela, p53_mdm2_independent_deg_per_h=0.0)
    viab, _ = dose_response(no_e6, NUTLIN, CONC_UM, t_end_h=72.0, n_cells_per_conc=N_CELLS,
                            seed=1, p=Params())
    assert viab[-1] < 0.5, (
        f"with the E6 route switched off, HeLa should respond: {np.round(viab[1:], 3)}")


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
    print(f"{len(fns) - failed}/{len(fns)} MDM2/p53 gates pass")
    raise SystemExit(1 if failed else 0)
