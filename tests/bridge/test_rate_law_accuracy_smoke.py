"""The molecule→cell bridge must carry ACCURACY, not just precision.

`binding_to_hill` turns a binding ΔG into a K_d rate-law prior that a
Campaign-2 cell model consumes. The danger (identical to the one fixed
one layer down in docking): the input σ is REPRODUCIBILITY, and a tight
reproducibility bar on a wrong ΔG would hand the cell model an
over-confident K_d. On biotin/streptavidin a ±0.29 kcal/mol seed bar hid
an ~11.8 kcal/mol true error.

These tests pin that:
  * a docking-sourced prior on an untrustworthy target class WIDENS its
    CI to the measured accuracy and flags `trust`, even when the input σ
    is tiny — the biotin trap cannot propagate up;
  * an FEP-sourced prior is `uncalibrated` (docking's reliability table
    does not apply to FEP, and CellSim's FEP is not yet GPU-validated) —
    accuracy is flagged UNKNOWN, never faked;
  * a well-behaved (trypsin) target reports a trustworthy verdict.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.bridge import binding_to_hill  # noqa: E402

RT = 0.0019872041 * 298.15


def _ci_from_sigma(Kd_M, sigma):
    lo = Kd_M * math.exp(-1.96 * sigma / RT)
    hi = Kd_M * math.exp(+1.96 * sigma / RT)
    return lo, hi


def test_tight_seed_bar_on_untrustworthy_target_does_not_give_tight_ci():
    """The biotin case: tiny input σ, but streptavidin is
    do_not_trust_absolute (mean err ~5 kcal/mol). The CI must reflect the
    accuracy, not the reproducibility."""
    prior = binding_to_hill(
        -7.3, uncertainty_kcalmol=0.29,   # the misleadingly tight seed bar
        receptor="benchmarks/dock/1stp.pdb", method="vina")

    assert prior.trust == "do_not_trust_absolute"
    assert prior.accuracy_known
    assert prior.accuracy_kcalmol > 3.0

    # CI must match the ACCURACY σ (~4.98), not the input σ (0.29).
    Kd_M = prior.parameters["Kd_M"]
    lo, hi = prior.parameter_ci95["Kd_M"]
    exp_lo, exp_hi = _ci_from_sigma(Kd_M, prior.accuracy_kcalmol)
    assert abs(lo - exp_lo) / exp_lo < 1e-9
    assert abs(hi - exp_hi) / exp_hi < 1e-9

    # And it must be VASTLY wider than the tight-seed CI would have been.
    tight_lo, tight_hi = _ci_from_sigma(Kd_M, 0.29)
    assert (hi / lo) > 100 * (tight_hi / tight_lo)


def test_fep_source_is_uncalibrated_not_faked():
    """FEP uses a different method than docking, and CellSim's absolute
    FEP is not yet GPU-validated → accuracy UNKNOWN, never borrowed."""
    prior = binding_to_hill(-12.0, uncertainty_kcalmol=0.4,
                            method="FEP-DDM-MBAR")
    assert prior.trust == "uncalibrated"
    assert prior.accuracy_kcalmol is None
    assert not prior.accuracy_known
    assert prior.accuracy_basis == "fep-unvalidated"
    # CI still exists from precision, but the note warns it is not accuracy.
    assert prior.parameter_ci95  # precision CI present
    assert "rank-order" in prior.notes


def test_trustworthy_target_reports_trustworthy_verdict():
    """Trypsin/benzamidine is the well-behaved class — absolute usable."""
    prior = binding_to_hill(-6.0, uncertainty_kcalmol=0.2,
                            receptor="benchmarks/dock/3ptb.pdb",
                            method="vina")
    assert prior.trust == "trustworthy_absolute"
    assert prior.accuracy_known
    assert prior.accuracy_kcalmol < 1.5


def test_unknown_receptor_stays_uncalibrated():
    """NEVER GUESS: a receptor not in the table is uncalibrated."""
    prior = binding_to_hill(-8.0, uncertainty_kcalmol=0.5,
                            receptor="benchmarks/dock/9xyz.pdb",
                            method="vina")
    assert prior.trust == "uncalibrated"
    assert prior.accuracy_kcalmol is None
    assert prior.accuracy_basis == "docking-uncalibrated-target"
