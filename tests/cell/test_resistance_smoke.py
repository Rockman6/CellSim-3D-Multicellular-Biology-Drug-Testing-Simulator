"""Resistance evolution — treatment selects for what was already there.

Anchors:
  * the logistic resistant-fraction formula agrees with the raw
    two-exponential ratio (independent derivation);
  * treatment inverts selection: s < 0 without drug (fitness cost makes
    resistance lose), s > 0 under drug — the drug-holiday rationale;
  * relapse happens when the resistant clone regrows the population;
  * a resistant clone that is still killed gives a durable response.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.cell import (  # noqa: E402
    ResistantClone,
    resistance_outcome,
    CellFateParams,
)

_DAY = 86400.0


def test_logistic_matches_raw_two_exponential_ratio():
    """f(t) from the logistic form must equal N_r/(N_s+N_r) computed raw."""
    out = resistance_outcome(0.9, ResistantClone(initial_fraction=1e-3))
    f0 = out.initial_resistant_fraction
    for t in (0.0, 5 * _DAY, 30 * _DAY, 120 * _DAY):
        n_s = (1 - f0) * math.exp(out.k_sensitive_per_s * t)
        n_r = f0 * math.exp(out.k_resistant_per_s * t)
        raw = n_r / (n_s + n_r)
        assert abs(out.resistant_fraction_at(t) - raw) < 1e-12


def test_drug_inverts_the_direction_of_selection():
    """Fitness cost makes resistance LOSE without drug and WIN with it."""
    clone = ResistantClone(occupancy_scale=0.05, fitness_cost=0.2)
    no_drug = resistance_outcome(0.0, clone)
    on_drug = resistance_outcome(0.9, clone)
    assert no_drug.selection_coefficient_per_s < 0     # resistance costs
    assert on_drug.selection_coefficient_per_s > 0     # drug selects it
    # And the fraction moves accordingly.
    assert no_drug.resistant_fraction_at(60 * _DAY) < \
        no_drug.initial_resistant_fraction
    assert on_drug.resistant_fraction_at(60 * _DAY) > \
        on_drug.initial_resistant_fraction


def test_relapse_occurs_when_resistant_clone_survives():
    """Sensitive cells die, resistant ones grow → population returns."""
    out = resistance_outcome(0.95, ResistantClone(
        occupancy_scale=0.02, fitness_cost=0.1, initial_fraction=1e-4))
    assert out.fate_sensitive == "dying"
    assert out.fate_resistant == "proliferating"
    t_rel = out.time_to_relapse_s()
    assert t_rel is not None
    assert 0 < t_rel < 365 * _DAY
    # Population dips below baseline before relapsing.
    assert out.total_population_at(t_rel * 0.3) < 1.0
    # ...and the resistant clone dominates by then.
    assert out.resistant_fraction_at(t_rel) > 0.5


def test_durable_response_when_resistant_clone_also_dies():
    """If the resistance mechanism is weak, both clones die → no relapse."""
    out = resistance_outcome(0.99, ResistantClone(
        occupancy_scale=0.95, fitness_cost=0.1, initial_fraction=1e-4))
    assert out.fate_resistant == "dying"
    assert out.time_to_relapse_s() is None


def test_no_response_reports_no_relapse():
    """Untreated population never dips, so there is no relapse to report."""
    out = resistance_outcome(0.0, ResistantClone())
    assert out.time_to_relapse_s() is None


def test_rejects_bad_clone_parameters():
    for bad in (dict(occupancy_scale=0.0), dict(occupancy_scale=1.5),
                dict(fitness_cost=1.0), dict(initial_fraction=1.0)):
        try:
            ResistantClone(**bad)
            assert False, bad
        except ValueError:
            pass
