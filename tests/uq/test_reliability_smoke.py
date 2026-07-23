"""Predictions must carry an ACCURACY estimate, not just precision.

The uncertainty CellSim showed users came from Monte-Carlo over Vina
seeds — that is reproducibility ("does Vina agree with itself?"), not
accuracy ("is Vina right?"). Measured on biotin/streptavidin the seed
bar was ±0.29 kcal/mol while the true error vs experiment was 11.8:
the displayed bar understated the real error ~41×.

These tests pin the reliability lookup that fixes it, and — most
importantly — pin that an UNCALIBRATED receptor reports UNKNOWN rather
than borrowing a neighbouring target class's number.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.uq import reliability_for  # noqa: E402


def test_calibrated_targets_report_measured_error():
    """Each calibrated receptor returns its measured class error."""
    rel = reliability_for("benchmarks/dock/3ptb.pdb")
    assert rel.calibrated
    assert rel.target_class == "trypsin_like"
    assert rel.verdict == "trustworthy_absolute"
    assert 0.0 < rel.mean_abs_err_kcalmol < 1.5

    rel = reliability_for("benchmarks/dock/1stp.pdb")
    assert rel.calibrated
    assert rel.verdict == "do_not_trust_absolute", (
        "streptavidin is the ultra-tight-binder class where Vina "
        "saturates; absolute ΔG must not be advertised as usable")
    assert rel.worst_abs_err_kcalmol > 10.0


def test_uncalibrated_target_reports_unknown_not_a_guess():
    """The load-bearing safety property: never invent an error bar."""
    rel = reliability_for("benchmarks/dock/9zzz.pdb")
    assert rel.calibrated is False
    assert rel.verdict == "uncalibrated"
    assert rel.mean_abs_err_kcalmol is None, (
        "an uncalibrated target must NOT inherit another class's "
        "error bar — unknown must stay unknown")
    assert "UNKNOWN" in rel.summary()


def test_receptor_reference_forms_all_resolve():
    """Path, bare id, and upper-case filename resolve identically."""
    a = reliability_for("benchmarks/dock/1m17.pdb")
    b = reliability_for("1m17")
    c = reliability_for("1M17.PDB")
    assert a.target_class == b.target_class == c.target_class
    assert a.target_class == "kinase_atp_site"


def test_accuracy_and_precision_are_not_confused():
    """Streptavidin's measured accuracy must dwarf a typical seed CI.

    Pins the actual bug: seed scatter there is ~0.3 kcal/mol while the
    measured accuracy error is ~5 (worst 11.6). If these ever converge,
    someone has wired precision into the accuracy field.
    """
    rel = reliability_for("benchmarks/dock/1stp.pdb")
    assert rel.mean_abs_err_kcalmol > 2.0, (
        "measured accuracy for the saturating class should be large; "
        "a small value suggests seed-scatter leaked into this field")


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
