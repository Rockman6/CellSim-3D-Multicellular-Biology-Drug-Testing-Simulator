"""Agent-based cells on a lattice — individual cells living, dividing, dying.

Everything upstream is a population average: a rate, a fraction, a
concentration. This module drops to individuals. Each cell occupies a
lattice site, senses the drug concentration AT ITS OWN POSITION, and each
step either divides into a free neighbouring site, dies, or does nothing —
sampled from the rates the `fate` layer computes.

Three things emerge here that a rate equation cannot express:

  * **Contact inhibition.** A cell can only divide if a neighbouring site
    is free, so a dense colony grows at its SURFACE, not exponentially.
    The switch from exponential to surface-limited growth is emergent, not
    imposed.
  * **Spatial structure.** With a drug gradient (vessel on one edge), the
    kill zone and the surviving population are actual regions on the grid,
    and survivors regrow into the space the dead cells vacated.
  * **Stochastic extinction.** A small population can die out by chance
    even when the average rate says it should grow — the reason minimal
    residual disease is a probability, not a threshold.

Rate → per-step probability uses the exact Poisson relation
p = 1 − exp(−k·Δt), not the p ≈ k·Δt approximation, so the simulation is
correct at any step size (a test anchors this against the analytic curve).

Deterministic given a seed (numpy default_rng), matching the project's
reproducibility discipline.

Non-AI: stochastic simulation of closed-form rates. No learned model.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from src.cell.fate import CellFateParams, fate_from_occupancy
from src.cell.occupancy import hill_occupancy

# 8-neighbour (Moore) connectivity for division into free space.
_NEIGHBOURS = [(-1, -1), (-1, 0), (-1, 1), (0, -1),
               (0, 1), (1, -1), (1, 0), (1, 1)]


@dataclass
class Clone:
    """A cell lineage with its own drug sensitivity and fitness."""

    name: str = "sensitive"
    occupancy_scale: float = 1.0      # engagement achieved vs the drug field
    fitness_cost: float = 0.0         # penalty on drug-free proliferation

    def __post_init__(self):
        if not (0.0 < self.occupancy_scale <= 1.0):
            raise ValueError("occupancy_scale must be in (0,1]")
        if not (0.0 <= self.fitness_cost < 1.0):
            raise ValueError("fitness_cost must be in [0,1)")


@dataclass
class ColonyTrajectory:
    """Population history and the final spatial state of a colony."""

    t_h: List[float]
    n_total: List[int]
    n_by_clone: Dict[str, List[int]]
    n_divisions: int
    n_deaths: int
    grid_size: Tuple[int, int]
    occupied: List[Tuple[int, int, str]]      # final (x, y, clone)
    seed: int
    extinct: bool

    @property
    def final_count(self) -> int:
        return self.n_total[-1] if self.n_total else 0

    def clone_fraction(self, name: str) -> float:
        series = self.n_by_clone.get(name)
        if not series or not self.final_count:
            return 0.0
        return series[-1] / self.final_count

    def occupancy_grid(self):
        """Dense 2-D array of clone indices (−1 = empty), for rendering."""
        import numpy as np
        names = sorted(self.n_by_clone.keys())
        idx = {n: i for i, n in enumerate(names)}
        g = -np.ones(self.grid_size, dtype=int)
        for x, y, cl in self.occupied:
            g[x, y] = idx[cl]
        return g, names

    def summary(self) -> str:
        state = "EXTINCT" if self.extinct else f"{self.final_count} cells"
        clones = ", ".join(
            f"{n} {100*self.clone_fraction(n):.0f}%"
            for n in sorted(self.n_by_clone) if self.clone_fraction(n) > 0)
        return (f"{state} after {self.t_h[-1]:.0f} h "
                f"({self.n_divisions} divisions, {self.n_deaths} deaths)"
                + (f"; {clones}" if clones else ""))


def simulate_colony(
    *,
    grid_size: Tuple[int, int] = (60, 60),
    n_seed_cells: int = 1,
    clones: Optional[List[Clone]] = None,
    seed_fractions: Optional[List[float]] = None,
    prior=None,
    drug_field: Optional[Callable[[int, int], float]] = None,
    params: Optional[CellFateParams] = None,
    dt_h: float = 1.0,
    n_steps: int = 200,
    rng_seed: int = 0,
    record_every: int = 1,
    seed_scattered: bool = False,
) -> ColonyTrajectory:
    """Simulate individual cells dividing and dying on a lattice.

    Args:
        grid_size: (nx, ny) lattice.
        n_seed_cells: how many cells to place initially (centred cluster).
        clones: lineages present (default: one sensitive clone).
        seed_fractions: initial composition (default: all first clone).
        prior: Hill K_d prior — required if `drug_field` is given.
        drug_field: (x, y) → drug concentration (M) at that site. None
            means no drug anywhere.
        params: cell-fate pharmacodynamics.
        dt_h: step size in hours. n_steps: number of steps.
        rng_seed: deterministic seed.
        seed_scattered: place the seed cells at random free sites
            (a DILUTE culture, where every cell has room to divide) rather
            than as one packed cluster (a COLONY, where interior cells are
            contact-inhibited from the start). The two give genuinely
            different growth curves — exponential vs surface-limited.

    Returns a `ColonyTrajectory`.
    """
    import numpy as np

    nx, ny = grid_size
    if nx < 1 or ny < 1:
        raise ValueError("grid_size must be positive")
    if dt_h <= 0 or n_steps < 1:
        raise ValueError("dt_h > 0 and n_steps ≥ 1 required")
    if drug_field is not None and prior is None:
        raise ValueError("a drug_field needs a K_d prior to compute occupancy")

    cl_list = clones or [Clone()]
    fracs = seed_fractions or ([1.0] + [0.0] * (len(cl_list) - 1))
    if len(fracs) != len(cl_list):
        raise ValueError("seed_fractions must match clones")

    p_base = params or CellFateParams()
    rng = np.random.default_rng(rng_seed)

    Kd = float(prior.parameters["Kd_M"]) if prior is not None else None
    n_hill = float(prior.parameters.get("n_hill", 1.0)) if prior else 1.0

    # Per-clone fate parameters (fitness cost on proliferation).
    clone_params = {
        c.name: CellFateParams(
            k_prolif_per_s=p_base.k_prolif_per_s * (1.0 - c.fitness_cost),
            k_death_basal_per_s=p_base.k_death_basal_per_s,
            k_maxkill_per_s=p_base.k_maxkill_per_s,
            cytostatic_fraction=p_base.cytostatic_fraction,
            transduction_tau=p_base.transduction_tau,
            source=p_base.source, trust=p_base.trust)
        for c in cl_list}
    clone_by_name = {c.name: c for c in cl_list}

    # Cache rates per (clone, site) — the drug field is static.
    rate_cache: Dict[Tuple[str, int, int], Tuple[float, float]] = {}

    def rates_at(clone_name: str, x: int, y: int) -> Tuple[float, float]:
        """(division rate, death rate) per hour at this site."""
        key = (clone_name, x, y)
        hit = rate_cache.get(key)
        if hit is not None:
            return hit
        p = clone_params[clone_name]
        occ = 0.0
        if drug_field is not None:
            conc = max(0.0, float(drug_field(x, y)))
            occ = hill_occupancy(conc, Kd, n_hill) \
                * clone_by_name[clone_name].occupancy_scale
        f = fate_from_occupancy(min(1.0, occ), p)
        from src.cell.fate import effect_from_occupancy
        E = effect_from_occupancy(min(1.0, occ), p.transduction_tau)
        k_div = max(0.0, p.k_prolif_per_s * (1.0 - p.cytostatic_fraction * E))
        k_die = max(0.0, p.k_death_basal_per_s + p.k_maxkill_per_s * E)
        out = (k_div * 3600.0, k_die * 3600.0)     # per hour
        rate_cache[key] = out
        return out

    # ---- seed the colony ------------------------------------------------
    occupied: Dict[Tuple[int, int], str] = {}
    cx, cy = nx // 2, ny // 2
    weights = np.array(fracs, dtype=float)
    weights = weights / weights.sum()

    if seed_scattered:
        n_place = min(n_seed_cells, nx * ny)
        flat = rng.choice(nx * ny, size=n_place, replace=False)
        for f in flat:
            x, y = int(f // ny), int(f % ny)
            occupied[(x, y)] = cl_list[
                int(rng.choice(len(cl_list), p=weights))].name

    # Spiral outward from the centre until enough sites are filled.
    ring = 0
    placed = len(occupied)
    while not seed_scattered and placed < n_seed_cells \
            and ring < max(nx, ny):
        for dx in range(-ring, ring + 1):
            for dy in range(-ring, ring + 1):
                if placed >= n_seed_cells:
                    break
                if max(abs(dx), abs(dy)) != ring:
                    continue
                x, y = cx + dx, cy + dy
                if 0 <= x < nx and 0 <= y < ny and (x, y) not in occupied:
                    name = cl_list[int(rng.choice(len(cl_list), p=weights))].name
                    occupied[(x, y)] = name
                    placed += 1
        ring += 1

    t_h: List[float] = [0.0]
    n_total: List[int] = [len(occupied)]
    n_by_clone: Dict[str, List[int]] = {
        c.name: [sum(1 for v in occupied.values() if v == c.name)]
        for c in cl_list}
    n_div = n_death = 0

    # ---- step ------------------------------------------------------------
    for step in range(1, n_steps + 1):
        if not occupied:
            break
        sites = list(occupied.keys())
        # Shuffle so update order carries no spatial bias.
        order = rng.permutation(len(sites))
        for i in order:
            site = sites[i]
            name = occupied.get(site)
            if name is None:          # died earlier this step
                continue
            x, y = site
            k_div, k_die = rates_at(name, x, y)
            p_die = 1.0 - math.exp(-k_die * dt_h)
            p_div = 1.0 - math.exp(-k_div * dt_h)
            u = rng.random()
            if u < p_die:
                del occupied[site]
                n_death += 1
                continue
            if rng.random() < p_div:
                free = [(x + a, y + b) for a, b in _NEIGHBOURS
                        if 0 <= x + a < nx and 0 <= y + b < ny
                        and (x + a, y + b) not in occupied]
                if free:              # contact inhibition when none free
                    tgt = free[int(rng.integers(len(free)))]
                    occupied[tgt] = name
                    n_div += 1

        if step % record_every == 0 or step == n_steps:
            t_h.append(step * dt_h)
            n_total.append(len(occupied))
            for c in cl_list:
                n_by_clone[c.name].append(
                    sum(1 for v in occupied.values() if v == c.name))

    return ColonyTrajectory(
        t_h=t_h, n_total=n_total, n_by_clone=n_by_clone,
        n_divisions=n_div, n_deaths=n_death, grid_size=(nx, ny),
        occupied=[(x, y, nm) for (x, y), nm in occupied.items()],
        seed=rng_seed, extinct=not occupied)
