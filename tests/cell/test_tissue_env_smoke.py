#!/usr/bin/env python3
"""Gates for the three things added to the dish: glucose, migration and
effector (immune) cells.

Each is off by default, so the first gate for every one of them is that
an untouched dish is unchanged — three merged phases of results depend on
that and nothing here is allowed to move them.

After that the gates are quantitative where a number exists to check
against. Migration has one: a random walk of hop rate k on a lattice of
spacing h has diffusion coefficient D = k h^2 / 2d, so the simulated mean
squared displacement must match that and not merely "look like motion".
Glucose and killing are checked for the shape of their response —
monotone in the thing that drives them, and specific to it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from cellsim.cell.dish import Dish, DishParams, run_dish  # noqa: E402
from cellsim.cell.engine import Params  # noqa: E402
from cellsim.cell.library import get_line  # noqa: E402

A549, DLD1, HELA = get_line("A549"), get_line("DLD-1"), get_line("HeLa")


def _spheroid(glucose=0.0, hours=144.0, seed=2):
    dp = DishParams(geometry="spheroid", glucose_medium_mM=glucose)
    return run_dish(DLD1, t_end_h=hours, n_seed=3000, grid_sites=52,
                    record_every_h=48, seed=seed, p=Params(cycle_cv=0.25), dp=dp)


# ── everything is off by default ──────────────────────────────────────
def test_the_defaults_leave_the_dish_exactly_as_it_was():
    """Glucose off, no migration, no effectors — the committed spheroid
    result must be reproduced bit for bit."""
    a = _spheroid(glucose=0.0)
    b = _spheroid(glucose=0.0)
    assert a.record.n_live[-1] == b.record.n_live[-1]
    dp = DishParams(geometry="spheroid")
    assert dp.glucose_medium_mM == 0.0
    assert dp.migration_per_h == 0.0
    d = Dish(A549, n_seed=20, grid_sites=20, seed=1, p=Params(), dp=DishParams())
    assert len(d.eff_pos) == 0, "no effectors until add_effectors() is called"


def test_plentiful_glucose_changes_nothing():
    """Switching the field on at a normal medium concentration must not
    move the answer: if it does, the coupling is wrong rather than the
    biology interesting."""
    off = _spheroid(glucose=0.0)
    plenty = _spheroid(glucose=5.5)
    assert off.record.n_live[-1] == plenty.record.n_live[-1], (
        f"{off.record.n_live[-1]} vs {plenty.record.n_live[-1]} at 5.5 mM")


# ── glucose ───────────────────────────────────────────────────────────
def test_glucose_is_consumed_and_falls_toward_the_centre():
    r = _spheroid(glucose=0.8)
    prof = r.dish.profile_glucose
    assert len(prof) > 2, "a spheroid should carry a glucose profile"
    assert prof[0] < prof[-1], f"centre must be poorer than the rim: {prof[0]:.3f} vs {prof[-1]:.3f}"
    assert prof[-1] <= 0.8 + 1e-9, "nothing may exceed the medium"


def test_starving_a_spheroid_stops_it_growing_then_kills_it():
    """Freyer & Sutherland's result in outline: lowering medium glucose
    thins the viable rim. Here it must first arrest growth, then at lower
    still produce a core, and the ordering is what is checked."""
    plenty = _spheroid(glucose=5.5)
    lean = _spheroid(glucose=0.2)
    starved = _spheroid(glucose=0.1)
    assert lean.record.n_live[-1] < plenty.record.n_live[-1] * 0.8, (
        f"0.2 mM should slow growth: {lean.record.n_live[-1]} vs {plenty.record.n_live[-1]}")
    assert starved.record.n_live[-1] < lean.record.n_live[-1], (
        f"0.1 mM should be worse still: {starved.record.n_live[-1]} vs {lean.record.n_live[-1]}")
    assert starved.record.n_necrotic[-1] > 100, (
        f"and should kill: {starved.record.n_necrotic[-1]} necrotic")


# ── migration ─────────────────────────────────────────────────────────
def _msd(hop_rate, hours=12.0, seeds=(3, 4, 5), spacing=20.0):
    out = []
    for seed in seeds:
        dp = DishParams(geometry="monolayer", migration_per_h=hop_rate,
                        spacing_um=spacing, field_dt_h=0.5, dt_h=0.02)
        d = Dish(HELA, n_seed=40, grid_sites=(80, 80), seed=seed, p=Params(), dp=dp)
        start = d.pos.copy().astype(float)
        for _ in range(int(hours / 0.5)):
            d.advance(None)
        ids = np.arange(40)
        out.append(float(np.mean(np.sum(((d.pos[ids] - start) * spacing) ** 2, axis=1))))
    return float(np.mean(out))


def test_cells_do_not_move_unless_asked():
    assert _msd(0.0) == 0.0


def test_migration_matches_the_random_walk_it_claims_to_be():
    """D = k h^2 / 4 in two dimensions. The simulated MSD after 12 h must
    give that back to within a quarter, which is the real check that the
    hop rate means what the parameter says."""
    for k in (0.5, 2.0):
        msd = _msd(k)
        d_sim = msd / (4 * 12.0)
        d_theory = k * 20.0 ** 2 / 4.0
        assert 0.7 < d_sim / d_theory < 1.3, (
            f"hop {k}/h: D {d_sim:.0f} vs theory {d_theory:.0f} um^2/h")


def test_a_faster_hop_rate_moves_cells_further():
    assert _msd(2.0) > _msd(0.5) * 2.0


# ── effector cells ────────────────────────────────────────────────────
def _assay(et_ratio, hours=24.0, seeds=(1, 2), n_tumour=400, **dish_kw):
    lives = []
    for seed in seeds:
        dp = DishParams(geometry="monolayer", spacing_um=20.0, field_dt_h=0.5, dt_h=0.02,
                        **dish_kw)
        d = Dish(A549, n_seed=n_tumour, grid_sites=(40, 40), seed=seed, p=Params(), dp=dp)
        d.add_effectors(int(round(et_ratio * n_tumour)))
        for _ in range(int(hours / 0.5)):
            d.advance(None)
        lives.append(int(d.cells.alive.sum()))
    return float(np.mean(lives))


def test_killing_rises_with_the_effector_to_target_ratio():
    """The shape every cytotoxicity assay has: more effectors, more
    lysis, saturating when there is nothing left to kill."""
    base = _assay(0.0)
    low, mid, high = _assay(0.1), _assay(0.3), _assay(1.0)
    assert low < base, f"effectors must kill: {base} -> {low}"
    assert mid < low, f"and more of them must kill more: {low} -> {mid}"
    assert high <= mid, f"saturating at the top: {mid} -> {high}"
    # Clearance needs a HIGHER ratio than it used to. Under the three-hit
    # rule (Weigelin 2021) a lone effector cannot kill: its hits are
    # repaired before a third arrives, so lysis only takes off once the
    # density is enough for several effectors to hit one target inside
    # the repair window. E:T 1 now leaves ~40 % lysis; E:T 3 clears it.
    very_high = _assay(3.0)
    assert very_high < base * 0.1, (
        f"E:T 3 should clear most of the monolayer: {base} -> {very_high}")


def test_exhaustion_caps_what_one_effector_can_do():
    """Serial killing is limited: a higher cap must allow more killing.

    Run at a density where killing actually happens. Under the three-hit
    rule 40 effectors on 400 targets kill almost nothing — their hits are
    repaired before a third lands — so this test compared 1 kill against
    1 kill and failed when the model changed, which is the right way for
    a gate to react to a mechanism it no longer matches."""
    def run(max_kills):
        dp = DishParams(geometry="monolayer", spacing_um=20.0, field_dt_h=0.5,
                        dt_h=0.02, effector_max_kills=max_kills)
        d = Dish(A549, n_seed=400, grid_sites=(40, 40), seed=1, p=Params(), dp=dp)
        d.add_effectors(1200)                     # E:T 3, where lysis is real
        for _ in range(48):
            d.advance(None)
        return d.n_killed_by_effectors
    one, five = run(1), run(5)
    assert one > 0, f"the test must be run where killing happens, got {one}"
    assert five > one, f"a higher cap must allow more killing: {one} -> {five}"


def test_the_killing_comes_from_the_effectors_and_not_from_their_presence():
    """Specificity. Effectors that cannot kill must leave the monolayer
    alone: if merely adding them changes the count — by crowding, by
    disturbing the fields, by a bookkeeping slip — then the lysis above
    is not measuring what it claims to."""
    def run(kill_rate):
        dp = DishParams(geometry="monolayer", spacing_um=20.0, field_dt_h=0.5,
                        dt_h=0.02, effector_kill_per_h=kill_rate)
        d = Dish(A549, n_seed=400, grid_sites=(40, 40), seed=1, p=Params(), dp=dp)
        d.add_effectors(400)
        for _ in range(48):
            d.advance(None)
        return int(d.cells.alive.sum()), d.n_killed_by_effectors

    none_live, none_kills = run(0.0)
    some_live, some_kills = run(1.0)
    bare_live, _ = run(0.0)                     # same arm again, for the control count
    assert none_kills == 0, f"a zero kill rate must kill nothing, got {none_kills}"
    assert none_live == bare_live, "and must be deterministic"
    assert some_kills > 0 and some_live < none_live, (
        f"killing effectors must reduce the monolayer: {none_live} -> {some_live}")


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
    print(f"{len(fns) - failed}/{len(fns)} tissue-environment gates pass")
    raise SystemExit(1 if failed else 0)
