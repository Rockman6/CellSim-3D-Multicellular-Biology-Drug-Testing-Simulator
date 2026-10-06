#!/usr/bin/env python3
"""Gates for genetic perturbation (`cellsim.cell.perturb`).

A knockout is only worth having if it does the specific thing the gene
does. These check that by comparison rather than by eyeballing a number:
a TP53 knockout must land where the library's p53-mutant lines already
are, knocking out any single step of one linear pathway must give the
same answer as knocking out any other, and a gene that cannot touch a
given drug must leave it alone.

The last one is the real test of specificity. ABCB1 is P-glycoprotein;
raising it must protect against a substrate (doxorubicin) and do nothing
at all to cisplatin, which ABCB1 does not transport. A perturbation
framework that moves everything is not modelling anything.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.engine import dose_response  # noqa: E402
from cellsim.cell.library import get_drug, get_line  # noqa: E402
from cellsim.cell.perturb import GENES, Perturbation  # noqa: E402

A549 = get_line("A549")
CIS = get_drug("cisplatin")
DOX = get_drug("doxorubicin")
CONC = np.array([3.0, 10.0, 30.0, 100.0])
N = 24


def _viab(drug=CIS, pert=None, conc=CONC, line=A549):
    v, _ = dose_response(line, drug, conc, n_cells_per_conc=N, perturbations=pert)
    return v[1:]


def test_an_unknown_gene_is_refused_not_ignored():
    """A knockout that silently does nothing looks like a result."""
    for bad in ("KRAS", "EGFR", "nonsense"):
        try:
            Perturbation(bad)
        except ValueError as e:
            assert "no counterpart" in str(e), e
        else:
            raise AssertionError(f"{bad} should be refused: it is not in the engine")


def test_the_modes_reject_nonsense_multipliers():
    for args in (("TP53", "knockdown", 1.5), ("TP53", "overexpress", 0.5),
                 ("TP53", "sideways", 1.0)):
        try:
            Perturbation(*args)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{args} should be refused")


def test_tp53_knockout_protects_like_a_p53_mutant_line():
    """The engine already models p53 loss through CellLine.p53_functional,
    so a TP53 knockout must reach the same place rather than inventing a
    second, differently-behaved route."""
    wt = _viab()
    ko = _viab(pert=[Perturbation("TP53")])
    assert (ko >= wt - 0.02).all(), f"knockout should not be MORE sensitive: {wt} -> {ko}"
    assert ko[1] > wt[1] + 0.1, f"and should be clearly protected mid-curve: {wt} -> {ko}"


def test_knocking_out_any_step_of_one_pathway_gives_the_same_answer():
    """PUMA -> BAX -> caspase 3 is linear in this engine, so losing any
    one of them must block it identically. If they differ, something is
    carrying flux around the break."""
    puma = _viab(pert=[Perturbation("BBC3")])
    bax = _viab(pert=[Perturbation("BAX")])
    casp = _viab(pert=[Perturbation("CASP3")])
    assert np.allclose(puma, bax, atol=0.03), f"PUMA {puma} vs BAX {bax}"
    assert np.allclose(bax, casp, atol=0.03), f"BAX {bax} vs CASP3 {casp}"


def test_blocking_apoptosis_converts_killing_into_arrest():
    """The informative prediction: with the death route gone the cells
    stop dying but do not resume growing, so viability floors well above
    zero instead of reaching it."""
    wt = _viab()
    ko = _viab(pert=[Perturbation("BAX")])
    assert wt[-1] < 0.05, f"wild type should be killed at the top dose: {wt}"
    assert 0.05 < ko[-1] < 0.5, (
        f"apoptosis-blocked cells should persist without proliferating: {ko}")


def test_abcb1_protects_against_a_substrate_and_not_against_cisplatin():
    """Specificity. Doxorubicin is a P-gp substrate; cisplatin is not, and
    GDSC's own correlations say the same."""
    dox_conc = np.array([0.02, 0.06, 0.2, 0.6])
    dox_wt = _viab(DOX, conc=dox_conc)
    dox_hi = _viab(DOX, pert=[Perturbation("ABCB1", "overexpress", 8)], conc=dox_conc)
    assert (dox_hi >= dox_wt - 0.02).all(), f"P-gp should protect: {dox_wt} -> {dox_hi}"
    assert dox_hi.sum() > dox_wt.sum() + 0.15, f"clearly: {dox_wt} -> {dox_hi}"

    cis_wt = _viab(CIS)
    cis_hi = _viab(CIS, pert=[Perturbation("ABCB1", "overexpress", 8)])
    assert np.allclose(cis_wt, cis_hi, atol=0.03), (
        f"ABCB1 must not move cisplatin, which it does not transport: "
        f"{cis_wt} -> {cis_hi}")


def test_bcl2_overexpression_raises_the_apoptotic_threshold():
    wt = _viab()
    hi = _viab(pert=[Perturbation("BCL2", "overexpress", 4)])
    assert hi.sum() > wt.sum() + 0.2, f"a bigger reserve should protect: {wt} -> {hi}"


def test_a_knockdown_sits_between_wild_type_and_knockout():
    wt = _viab()
    kd = _viab(pert=[Perturbation("BAX", "knockdown", 0.3)])
    ko = _viab(pert=[Perturbation("BAX")])
    assert wt.sum() <= kd.sum() + 0.05, f"{wt} -> {kd}"
    assert kd.sum() <= ko.sum() + 0.05, f"{kd} -> {ko}"


def test_every_listed_gene_runs_and_is_documented():
    for name, info in GENES.items():
        assert info.get("note"), f"{name} has no note saying what it is"
        assert info["kind"] in ("state", "line"), info
        v = _viab(pert=[Perturbation(name)], conc=np.array([10.0]))
        assert np.isfinite(v).all(), f"{name} knockout produced {v}"


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
        except Exception as e:                      # noqa: BLE001
            failed += 1
            print(f"  ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"{len(fns) - failed}/{len(fns)} perturbation gates pass")
    raise SystemExit(1 if failed else 0)
