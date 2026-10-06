"""cellsim.cell.dish — the engine's cells, in space.

The well-mixed engine answers "how many cells survive"; the dish answers
"which ones, and where". Each lattice site holds at most one cell and
each cell carries a full engine state row, so the biology validated in
`engine.py` (cell cycle, p53, apoptosis, drug uptake and efflux) runs per
cell, unchanged. Space adds what a stirred population cannot have:

* **Contact inhibition.** A cell divides by pushing its neighbours toward
  the nearest free site, and only if that site is within `push_sites`.
  Without room its G1 growth signal drops to zero and it waits in G1.
* **Fields.** Oxygen and drug diffuse in from the medium and are consumed
  by the cells they pass, so cells deep in an aggregate see less of both.
* **Hypoxic quiescence and necrosis.** G1 progression slows with oxygen
  and stops below a threshold; near-anoxia kills by necrosis within
  hours, which leaves debris rather than a cleared apoptotic body.
* **Inheritance, and optionally evolution.** A daughter takes its
  mother's heterogeneity (reserve, uptake, cycle speed) and lineage, so
  clones form contiguous patches. With `DishParams.mutation_rate` a
  daughter may instead differ heritably from its mother, which is what
  lets resistance arise during treatment rather than only be selected
  from the tail the population began with.

Geometries
----------
``monolayer``: a 2-D sheet under well-stirred medium. Oxygen and drug are
uniform; only contact inhibition is spatial.
``spheroid``: a 3-D aggregate in medium. Oxygen and drug are solved on a
radial grid about the aggregate's centroid, which is exact for a
spherically symmetric aggregate and is how the spheroid literature
models it (Grimes et al. 2014 J R Soc Interface 11:20131124).

Every environmental constant is an input with a source (DishParams); the
only per-line inputs are the engine's. Non-AI: ODEs, finite volumes and
a stochastic lattice. Deterministic for a given seed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable, Optional, Sequence, Union

import numpy as np
from scipy import ndimage
from scipy.linalg import solve_banded

from cellsim.cell.engine import (IX, N_BASE, Params, Population, _as_drug_tuple, _cell_events,
                                 _phase, _rk4, calibrate_cycle_scale, make_population)
from cellsim.cell.library import CellLine, Drug

EMPTY, NECROTIC, APOPTOTIC = -1, -2, -3


@dataclass(frozen=True)
class DishParams:
    """Spatial and environmental constants. Cell biology stays in Params."""
    geometry: str = "monolayer"          # 'monolayer' or 'spheroid'
    spacing_um: float = 15.0             # one site = one cell; DLD-1 suspended radius ~7.5 um (Grimes 2014)
    # How far a dividing cell can push its neighbours to reach free space,
    # in sites. Small values give the surface-limited growth of strong
    # contact inhibition; inf removes any mechanical limit, leaving
    # nutrients to decide (the classical spheroid assumption). This is the
    # dish's one CALIBRATED constant: DLD-1 spheroids grow ~15 um/day in
    # radius once their core is anoxic (Grimes 2014, days 6-17), and a
    # reach of 1 gives 15.5 um/day where 2, 3 and inf give 28, 35 and 59
    # (scripts/validate_dish.py). Solid stress is known to limit spheroid
    # proliferation this way (Helmlinger et al. 1997 Nat Biotechnol 15:778).
    push_sites: float = 1.0
    field_dt_h: float = 0.24             # field update and lattice housekeeping interval
    dt_h: float = 0.02                   # engine step, as in engine.simulate
    # Oxygen, from DLD-1 spheroids (Grimes et al. 2014 J R Soc Interface
    # 11:20131124): D = 2e-9 m2/s, medium 100 mmHg, consumption
    # a*Omega = 7.29e-7 m3 kg-1 s-1 x 3.0318e7 mmHg kg m-3 = 22.1 mmHg/s in
    # packed viable tissue, which gives their measured 233 um diffusion limit.
    o2_medium_mmHg: float = 100.0
    o2_D_um2_per_s: float = 2000.0
    o2_uptake_mmHg_per_s: float = 22.1
    o2_Km_mmHg: float = 1.0              # half-saturation of mitochondrial uptake (~1 uM)
    # Oxygen-dependent proliferation, PhysiCell's documented defaults
    # (Ghaffarizadeh et al. 2018 PLoS Comput Biol 14:e1005991): G1 runs at
    # full speed at or above 38 mmHg (5 % O2), slowing linearly to a halt
    # at 5 mmHg. Necrosis at near-anoxia, because Grimes et al. define the
    # DLD-1 necrotic core as the anoxic one and the oxygen numbers above
    # are theirs; PhysiCell's generic 5 / 2.5 mmHg put the core ~30 um
    # outside the anoxic boundary (scripts/validate_dish.py). Full rate
    # 1/(6 h): cells survive about 6 h without oxygen (PhysiCell).
    o2_g1_full_mmHg: float = 38.0
    o2_g1_zero_mmHg: float = 5.0
    o2_necrosis_onset_mmHg: float = 1.0
    o2_necrosis_full_mmHg: float = 0.5
    necrosis_rate_per_h: float = 1.0 / 6.0
    # Apoptotic bodies occupy their site until cleared: PhysiCell's 516 min
    # apoptosis duration (Macklin et al. 2012 J Theor Biol 301:122).
    apoptotic_clearance_h: float = 8.6
    boundary_layer_um: float = 0.0       # unstirred medium between aggregate and bulk
    extracellular_fraction: float = 0.4  # interstitial share of a packed site's volume
    # Small-molecule interstitial diffusivity, 1e-6 cm2/s (Nugent & Jain
    # 1984 Cancer Res 44:238), the same default as cellsim.cell.tissue.
    drug_D_um2_per_s: float = 100.0
    # Heritable change at division. With probability `mutation_rate` a
    # daughter's drug accumulation is multiplied by a log-normal factor of
    # log-sd `mutation_effect_sd` (its reserve likewise), and the change is
    # inherited by that daughter's own descendants. This is what makes
    # resistance EVOLVE rather than merely be selected from the tail the
    # population started with. Both default to 0, which leaves inheritance
    # exactly clonal.
    #
    # The rate is a modelling choice, not a measurement: resistance here
    # arises from many routes at once (expression changes, copy number,
    # point mutations), whose combined rate per division is not a published
    # constant. scripts/experiment_evolution.py reports what a given rate
    # produces, against the 2-8x fold-resistance that derived lines
    # actually show (McDermott et al. 2014 Front Oncol 4:40).
    mutation_rate: float = 0.0
    mutation_effect_sd: float = 0.0
    # ── glucose ────────────────────────────────────────────────────────
    # The second thing a spheroid runs out of, and often the first: in
    # EMT6/Ro spheroids the viable rim thickness tracks the medium's
    # glucose as much as its oxygen (Freyer & Sutherland 1986 Cancer Res
    # 46:3504; 1986 Cancer Res 46:3513), and lowering glucose from 16.5 to
    # 0.8 mM thins the rim and raises respiration.
    #
    # Standard DMEM is 25 mM, RPMI 11 mM; 5.5 mM is physiological. The
    # diffusivity is the small-molecule value and consumption is per unit
    # volume of packed cells. `glucose_medium_mM = 0` switches the whole
    # field off, which is the default, so every earlier result stands.
    glucose_medium_mM: float = 0.0
    glucose_D_um2_per_s: float = 500.0        # ~5e-6 cm2/s, Casciari 1988
    glucose_uptake_mM_per_s: float = 0.012    # packed tissue, Casciari 1992
    glucose_Km_mM: float = 0.04               # GLUT/hexokinase half-saturation
    # Below this a cell cannot hold its ATP and dies; between this and
    # `glucose_quiescent_mM` it survives but does not cycle.
    glucose_death_mM: float = 0.06
    glucose_quiescent_mM: float = 0.30
    # ── migration ──────────────────────────────────────────────────────
    # A random walk on the lattice: each live cell that has an empty
    # neighbouring site hops into one at this rate. On a lattice of
    # spacing h, a hop rate k gives a diffusion coefficient
    # D = k h^2 / (2 d) in d dimensions, which is how a measured cell
    # motility converts to this number. HeLa in the Cell Tracking
    # Challenge movie moves of order 10 um/h; 0 (the default) keeps every
    # earlier result, where cells move only when pushed by a division.
    migration_per_h: float = 0.0
    # Chemotaxis toward oxygen: the fraction of hops biased up the local
    # gradient rather than taken at random (0 = unbiased, 1 = always up).
    migration_o2_bias: float = 0.0
    # ── effector (immune) cells ───────────────────────────────────────
    # Cytotoxic cells added on top of the tumour layer, as T or NK cells
    # sit on a monolayer in a killing assay. They are much smaller than
    # tumour cells and crawl over them, so they live on their own layer
    # and do not compete for tumour sites.
    #
    # Each effector random-walks, biased toward tumour cells; on reaching
    # one it forms a conjugate and kills it at `effector_kill_per_h`, is
    # then busy for `effector_conjugate_h` before it can engage again, and
    # stops killing after `effector_max_kills` (exhaustion). CTLs take
    # about an hour from contact to target death and kill a few targets a
    # day each (Halle et al. 2016 Immunity 44:233), with serial killing
    # capped by exhaustion. A killed target dies by apoptosis, exactly as
    # a drug-killed one does.
    effector_speed_per_h: float = 3.0          # hops/h on the effector layer
    effector_tumour_bias: float = 0.7          # fraction of hops toward the nearest target
    effector_kill_per_h: float = 1.0           # kill rate while conjugated
    effector_conjugate_h: float = 1.0          # busy after each kill
    effector_max_kills: int = 5
    effector_lifetime_h: float = 96.0          # effectors die off in culture

    @property
    def site_volume_um3(self) -> float:
        return self.spacing_um ** 3


# ── one-dimensional field solvers (radial shells or planar slabs) ─────
def _geometry(edges: np.ndarray, kind: str):
    """Face areas, shell volumes and centre-to-centre distances."""
    r = edges
    if kind == "sphere":
        area = 4.0 * math.pi * r ** 2
        vol = 4.0 / 3.0 * math.pi * (r[1:] ** 3 - r[:-1] ** 3)
    elif kind == "plane":
        area = np.ones_like(r)
        vol = r[1:] - r[:-1]
    else:
        raise ValueError(f"unknown field geometry {kind!r}")
    centres = 0.5 * (r[1:] + r[:-1])
    gaps = np.empty(len(r))
    gaps[0] = np.inf                       # symmetry / no-flux face
    gaps[1:-1] = centres[1:] - centres[:-1]
    gaps[-1] = r[-1] - centres[-1]         # outer face to the Dirichlet boundary
    return area, vol, gaps


def solve_steady(edges: np.ndarray, kind: str, D: float, boundary: float, uptake: np.ndarray,
                 Km: float, guess: Optional[np.ndarray] = None, iters: int = 50) -> np.ndarray:
    """Steady reaction-diffusion with saturable consumption,

        D * div(grad c) = uptake(r) * c / (Km + c),

    no flux at r = 0 and c = `boundary` at r = edges[-1]. `uptake` is the
    maximal consumption per unit volume in each shell (concentration per
    second). Finite volumes on shells, Newton iterations on a tridiagonal
    Jacobian; returns the shell-centre values."""
    area, vol, gaps = _geometry(edges, kind)
    n = len(vol)
    tin = D * area[:-1] / gaps[:-1]        # conductance of each shell's inner face
    tin[0] = 0.0
    tout = D * area[1:] / gaps[1:]         # ... and outer face (last one to the boundary)
    c = np.full(n, boundary, float) if guess is None or len(guess) != n else guess.copy()
    for _ in range(iters):
        c = np.maximum(c, 0.0)
        flux_in = np.zeros(n)
        flux_in[1:] += tin[1:] * (c[:-1] - c[1:])
        flux_out = np.zeros(n)
        flux_out[:-1] = tout[:-1] * (c[1:] - c[:-1])
        flux_out[-1] = tout[-1] * (boundary - c[-1])
        sink = vol * uptake * c / (Km + c)
        F = flux_in + flux_out - sink
        diag = -tin - tout - vol * uptake * Km / (Km + c) ** 2
        ab = np.zeros((3, n))
        ab[0, 1:] = tout[:-1]              # super-diagonal
        ab[1] = diag
        ab[2, :-1] = tin[1:]               # sub-diagonal
        dc = solve_banded((1, 1), ab, -F)
        c = c + dc
        if np.max(np.abs(dc)) < 1e-9 * max(boundary, 1e-12):
            break
    return np.maximum(c, 0.0)


def step_transient(edges: np.ndarray, kind: str, D: float, boundary: float, c: np.ndarray,
                   dt_s: float, a: np.ndarray, b: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """One backward-Euler step of eps * dc/dt = D * div(grad c) - a*c + b,
    with the same boundary conditions. `a` (1/s) and `b` (conc/s) are the
    linearised cellular uptake per unit volume; `eps` the extracellular
    fraction in each shell. Unconditionally stable."""
    area, vol, gaps = _geometry(edges, kind)
    n = len(vol)
    tin = D * area[:-1] / gaps[:-1]
    tin[0] = 0.0
    tout = D * area[1:] / gaps[1:]
    m = eps * vol / dt_s
    diag = m + tin + tout + vol * a
    rhs = m * c + vol * b
    rhs[-1] += tout[-1] * boundary
    ab = np.zeros((3, n))
    ab[0, 1:] = -tout[:-1]
    ab[1] = diag
    ab[2, :-1] = -tin[1:]
    return np.maximum(solve_banded((1, 1), ab, rhs), 0.0)


# ── the dish ───────────────────────────────────────────────────────────
DoseFn = Callable[[float], Union[float, Sequence[float]]]


@dataclass
class DishRecord:
    """Population-level readouts over time."""
    t_h: list = field(default_factory=list)
    n_live: list = field(default_factory=list)
    n_apoptotic: list = field(default_factory=list)
    n_necrotic: list = field(default_factory=list)
    phase_frac: list = field(default_factory=list)      # live cells, G1/S/G2/M
    quiescent_frac: list = field(default_factory=list)  # live G1 cells with no growth signal
    radius_um: list = field(default_factory=list)       # equivalent sphere / disc of all occupied sites
    necrotic_radius_um: list = field(default_factory=list)
    hypoxic_frac: list = field(default_factory=list)    # live cells below o2_g1_full
    mean_uptake: list = field(default_factory=list)     # live cells, 1 = the starting population
    mean_reserve: list = field(default_factory=list)
    n_mutations: list = field(default_factory=list)
    min_o2_mmHg: list = field(default_factory=list)
    confluence: list = field(default_factory=list)      # occupied share of the lattice
    blocked_divisions: list = field(default_factory=list)


class Dish:
    """A lattice of cells, each integrated by the engine. Use `run_dish`
    unless you need to step it yourself."""

    def __init__(self, line: CellLine, drug=None, *, n_seed: int = 100,
                 grid_sites=None, seed: int = 1, p: Params = Params(),
                 dp: DishParams = DishParams(), k_cyc: Optional[float] = None,
                 layout: str = "auto"):
        if dp.geometry not in ("monolayer", "spheroid"):
            raise ValueError(f"geometry must be 'monolayer' or 'spheroid', not {dp.geometry!r}")
        if abs(dp.field_dt_h / dp.dt_h - round(dp.field_dt_h / dp.dt_h)) > 1e-9:
            raise ValueError("field_dt_h must be a whole number of engine steps")
        self.line, self.p, self.dp = line, p, dp
        self.drugs = _as_drug_tuple(drug)
        self.n_drugs = max(1, len(self.drugs))
        self.k_cyc = calibrate_cycle_scale(line, p) if k_cyc is None else k_cyc
        self.rng = np.random.default_rng(seed)
        self.t_h = 0.0
        self.blocked = 0
        self.n_mutations = 0

        if isinstance(grid_sites, (tuple, list)):          # explicit lattice dimensions
            dims = tuple(int(v) for v in grid_sites)
            self.shape = (dims[0], dims[1], 1) if dp.geometry == "monolayer" else (
                dims + (dims[-1],) * (3 - len(dims)))[:3]
        elif dp.geometry == "monolayer":
            side = grid_sites or 100
            self.shape = (side, side, 1)
        else:
            side = grid_sites or 64
            self.shape = (side, side, side)
        self.occ = np.full(self.shape, EMPTY, dtype=np.int64)
        self.body_t = np.zeros(self.shape)              # time of death, for apoptotic bodies
        n_sites = int(np.prod(self.shape))
        if n_seed > n_sites:
            raise ValueError(f"{n_seed} cells do not fit on a {self.shape} lattice")

        # Initial cells: the engine's own asynchronous population.
        from cellsim.cell.engine import _quiescent_fraction
        pop = make_population(n_seed, 0.0, self.rng, het_sigma=p.het_sigma,
                              n_drugs=self.n_drugs, cycle_cv=p.cycle_cv,
                              quiescent_fraction=_quiescent_fraction(line, p))
        self.cells = pop
        self.lineage = np.arange(n_seed)
        # Birth time of each cell; -inf for the founders, whose cycle began
        # before the run. Every division by a cell born in the run logs one
        # complete cycle (birth -> own division), as a time-lapse would.
        self.birth_t = np.full(n_seed, -np.inf)
        self.cycles: list[tuple[float, float]] = []      # (birth_h, division_h)
        self.pos = self._seed_positions(n_seed, layout)
        self.occ[tuple(self.pos.T)] = np.arange(n_seed)

        # Fields: per-cell local values and, for a spheroid, radial profiles.
        self.o2 = np.full(n_seed, dp.o2_medium_mmHg)
        self.glucose = np.full(n_seed, dp.glucose_medium_mM)
        self.gs = np.ones(n_seed)
        self.profile_r_um = np.zeros(0)
        self.profile_o2 = np.zeros(0)
        self.profile_glucose = np.zeros(0)
        self.profile_drug = [np.zeros(0) for _ in range(self.n_drugs)]
        self._near = None                                # nearest-free-site map
        self._pending: dict = {}                         # daughters placed this step
        self._offsets = self._local_offsets(2)
        # Effector cells: none until add_effectors() is called.
        self.eff_pos = np.zeros((0, 3), dtype=np.int64)
        self.eff_alive = np.zeros(0, bool)
        self.eff_kills = np.zeros(0, dtype=np.int64)
        self.eff_busy_until = np.zeros(0)
        self.eff_born = np.zeros(0)
        self.n_killed_by_effectors = 0
        self.record = DishRecord()

    # ── set-up helpers ──
    def _seed_positions(self, n: int, layout: str) -> np.ndarray:
        nx, ny, nz = self.shape
        if layout == "auto":
            layout = "scattered" if self.dp.geometry == "monolayer" else "ball"
        if layout == "scattered":
            flat = self.rng.choice(nx * ny * nz, size=n, replace=False)
            return np.column_stack(np.unravel_index(flat, self.shape)).astype(np.int64)
        if layout in ("ball", "cluster"):
            grid = np.indices(self.shape).reshape(3, -1).T
            centre = (np.array(self.shape) - 1) / 2.0
            d = np.linalg.norm(grid - centre, axis=1)
            d = d + 1e-6 * self.rng.random(len(d))       # break ties reproducibly
            return grid[np.argsort(d)[:n]].astype(np.int64)
        raise ValueError(f"unknown layout {layout!r}")

    def _local_offsets(self, radius: int) -> np.ndarray:
        r = np.arange(-radius, radius + 1)
        dz = r if self.shape[2] > 1 else np.array([0])
        off = np.array(np.meshgrid(r, r, dz, indexing="ij")).reshape(3, -1).T
        off = off[np.any(off != 0, axis=1)]
        return off[np.argsort(np.linalg.norm(off, axis=1), kind="stable")]

    # ── geometry of the aggregate ──
    def centroid(self) -> np.ndarray:
        occupied = np.argwhere(self.occ != EMPTY)
        return occupied.mean(axis=0) if len(occupied) else (np.array(self.shape) - 1) / 2.0

    def equivalent_radius_um(self, mask: np.ndarray) -> float:
        n = int(mask.sum())
        if self.dp.geometry == "monolayer":
            return math.sqrt(n * self.dp.spacing_um ** 2 / math.pi)
        return (3.0 * n * self.dp.site_volume_um3 / (4.0 * math.pi)) ** (1.0 / 3.0)

    # ── fields ──
    def _radial_bins(self):
        """Shell index of every lattice site and of every cell."""
        h = self.dp.spacing_um
        c = self.centroid()
        grid = np.indices(self.shape).reshape(3, -1).T
        d_sites = np.linalg.norm(grid - c, axis=1) * h
        occupied = (self.occ != EMPTY).reshape(-1)
        r_out = (d_sites[occupied].max() if occupied.any() else 0.0) + 0.5 * h
        r_out += self.dp.boundary_layer_um
        dr = 0.5 * h
        nb = max(4, int(math.ceil(r_out / dr)))
        edges = np.linspace(0.0, nb * dr, nb + 1)
        site_bin = np.minimum((d_sites / dr).astype(int), nb)    # nb = beyond the boundary
        cell_d = np.linalg.norm(self.pos - c, axis=1) * h
        cell_bin = np.minimum((cell_d / dr).astype(int), nb - 1)
        return edges, site_bin, cell_bin

    def update_fields(self, dose: Sequence[float], dt_h: float) -> None:
        dp = self.dp
        n = len(self.pos)
        if dp.geometry == "monolayer":
            # A thin sheet under stirred medium: no gradients.
            self.o2 = np.full(n, dp.o2_medium_mmHg)
            self.glucose = np.full(n, dp.glucose_medium_mM)
            for i in range(self.n_drugs):
                self._set_c_out(i, np.full(n, float(dose[i])))
            return
        edges, site_bin, cell_bin = self._radial_bins()
        nb = len(edges) - 1
        sites_per_bin = np.bincount(site_bin, minlength=nb + 1)[:nb].astype(float)
        occ_flat = self.occ.reshape(-1)
        occupied_per_bin = np.bincount(site_bin[occ_flat != EMPTY], minlength=nb + 1)[:nb]
        live_per_bin = np.bincount(cell_bin, minlength=nb)
        with np.errstate(invalid="ignore", divide="ignore"):
            phi_live = np.where(sites_per_bin > 0, live_per_bin / sites_per_bin, 0.0)
            phi_occ = np.where(sites_per_bin > 0, occupied_per_bin / sites_per_bin, 0.0)
        # Oxygen: quasi-steady (it equilibrates in seconds).
        uptake = dp.o2_uptake_mmHg_per_s * np.clip(phi_live, 0.0, 1.0)
        guess = self.profile_o2 if len(self.profile_o2) == nb else None
        o2 = solve_steady(edges, "sphere", dp.o2_D_um2_per_s, dp.o2_medium_mmHg, uptake,
                          dp.o2_Km_mmHg, guess=guess)
        self.profile_r_um = 0.5 * (edges[1:] + edges[:-1])
        self.profile_o2 = o2
        self.o2 = o2[cell_bin]
        # Glucose, the same quasi-steady reaction-diffusion problem with
        # its own diffusivity and uptake. Solved only when the medium
        # carries any, so the default costs nothing.
        if dp.glucose_medium_mM > 0:
            g_uptake = dp.glucose_uptake_mM_per_s * np.clip(phi_live, 0.0, 1.0)
            guess_g = (self.profile_glucose if len(self.profile_glucose) == nb else None)
            glc = solve_steady(edges, "sphere", dp.glucose_D_um2_per_s,
                               dp.glucose_medium_mM, g_uptake, dp.glucose_Km_mM,
                               guess=guess_g)
            self.profile_glucose = glc
            self.glucose = glc[cell_bin]
        else:
            self.glucose = np.full(len(cell_bin), dp.glucose_medium_mM)
        # Drug: transient, with cellular uptake as a linearised sink.
        eps = 1.0 - np.clip(phi_occ, 0.0, 1.0) * (1.0 - dp.extracellular_fraction)
        shell_vol = _geometry(edges, "sphere")[1]
        v_cell = (1.0 - dp.extracellular_fraction) * dp.site_volume_um3
        for i, dg in enumerate(self.drugs):
            prof = self.profile_drug[i]
            if len(prof) != nb:                       # aggregate grew: carry the profile over
                prof = (np.interp(self.profile_r_um, self._last_r, prof)
                        if len(prof) and hasattr(self, "_last_r") else np.zeros(nb))
            cin = self.cells.Y[:, N_BASE + i]
            het = self.cells.het_uptake
            B = max(dg.accumulation_ratio, 1e-9)
            alpha = v_cell * B * het / dg.tau_uptake_h / 3600.0
            beta = v_cell * (B / dg.partition) * cin / dg.tau_uptake_h / 3600.0
            if self.p.efflux_vmax_uM_per_h > 0 and getattr(dg, "pgp_substrate", False):
                vmax = self.p.efflux_vmax_uM_per_h * self.line.efflux_level * het
                beta = beta - v_cell * (B / dg.partition) * vmax * cin / (self.p.efflux_km_uM + cin) / 3600.0
            a = np.bincount(cell_bin, weights=alpha, minlength=nb) / shell_vol
            b = np.bincount(cell_bin, weights=beta, minlength=nb) / shell_vol
            if dt_h > 0:
                prof = step_transient(edges, "sphere", dp.drug_D_um2_per_s, float(dose[i]), prof,
                                      dt_h * 3600.0, a, b, eps)
            self.profile_drug[i] = prof
            self._set_c_out(i, prof[cell_bin])
        self._last_r = self.profile_r_um.copy()

    def _set_c_out(self, i: int, values: np.ndarray) -> None:
        if self.cells.C_out.ndim == 1:
            self.cells.C_out[:] = values
        else:
            self.cells.C_out[:, i] = values

    # ── behaviour set by the environment ──
    def update_growth_signal(self) -> None:
        dp = self.dp
        span = max(dp.o2_g1_full_mmHg - dp.o2_g1_zero_mmHg, 1e-9)
        g_o2 = np.clip((self.o2 - dp.o2_g1_zero_mmHg) / span, 0.0, 1.0)
        # Glucose gates the cycle the same way, and the scarcer of the two
        # governs: a cell with oxygen but no sugar does not divide either.
        if dp.glucose_medium_mM > 0:
            gspan = max(dp.glucose_quiescent_mM - dp.glucose_death_mM, 1e-9)
            g_glc = np.clip((self.glucose - dp.glucose_death_mM) / gspan, 0.0, 1.0)
            g_o2 = np.minimum(g_o2, g_glc)
        blocked = self.occ != EMPTY
        if blocked.all():                            # nowhere left to go
            self._near = None
            self.gs = np.zeros(len(self.pos))
            return
        dist, idx = ndimage.distance_transform_edt(blocked, return_indices=True)
        self._near = idx
        reach = dist[tuple(self.pos.T)] if len(self.pos) else np.zeros(0)
        room = reach <= dp.push_sites + 1e-9
        self.gs = g_o2 * room

    def _aux(self) -> dict:
        c = self.cells
        aux = {"C_out": c.C_out, "arrest_h": c.arrest_h, "forced_D": c.forced_D,
               "het_bcl2": c.het_bcl2, "het_uptake": c.het_uptake, "gs": self.gs}
        if c.cycle_rate is not None:
            aux["cycle_rate"] = c.cycle_rate
        if c.quiescent is not None:
            aux["quiescent"] = c.quiescent
        return aux

    # ── division, death, clearance ──
    def _free_site_near(self, m: np.ndarray) -> Optional[np.ndarray]:
        """Nearest free site to `m`, from the map built this interval and
        refreshed locally where earlier placements have used it up."""
        if self._near is None:
            return None
        f = self._near[:, m[0], m[1], m[2]]
        if self.occ[tuple(f)] == EMPTY and not np.array_equal(f, m):
            return f
        cand = f + self._offsets
        ok = np.all((cand >= 0) & (cand < np.array(self.shape)), axis=1)
        cand = cand[ok]
        free = self.occ[tuple(cand.T)] == EMPTY
        if free.any():
            cand = cand[free]
            return cand[np.argmin(np.linalg.norm(cand - m, axis=1))]
        blocked = self.occ != EMPTY
        if blocked.all():
            return None
        _, idx = ndimage.distance_transform_edt(blocked, return_indices=True)
        self._near = idx
        f = idx[:, m[0], m[1], m[2]]
        return f if self.occ[tuple(f)] == EMPTY else None

    def _place(self, j: int, new_index: int) -> Optional[np.ndarray]:
        """Push a chain of neighbours one site toward the nearest free site
        and return the site freed next to mother `j`, now taken by her
        daughter; None if there is no room within reach."""
        m = self.pos[j].copy()
        f = self._free_site_near(m)
        if f is None or np.linalg.norm(f - m) > self.dp.push_sites + 1e-9:
            return None
        L = int(np.max(np.abs(f - m)))
        path = [m + np.rint(k * (f - m) / L).astype(np.int64) for k in range(L + 1)]
        end = L
        for k in range(1, L):                       # stop at any free site on the way
            if self.occ[tuple(path[k])] == EMPTY:
                end = k
                break
        for k in range(end, 1, -1):                 # shift the chain outward
            src, dst = tuple(path[k - 1]), tuple(path[k])
            who = self.occ[src]
            self.occ[dst] = who
            if who in self._pending:                # a daughter placed earlier this step
                self._pending[who] = path[k]
            elif who >= 0:
                self.pos[who] = path[k]
            elif who == APOPTOTIC:
                self.body_t[dst] = self.body_t[src]
        site = path[1]
        self.occ[tuple(site)] = new_index
        self._pending[new_index] = site
        return site

    def _divide(self, divide: np.ndarray) -> None:
        mothers = np.flatnonzero(divide)
        n0 = len(self.pos)
        for j in mothers:                       # the mother's cycle is complete either way
            if np.isfinite(self.birth_t[j]):
                self.cycles.append((float(self.birth_t[j]), self.t_h))
            self.birth_t[j] = self.t_h
        placed = []
        self._pending = {}
        for j in mothers:
            site = self._place(int(j), n0 + len(placed))
            if site is None:
                self.blocked += 1
                continue
            placed.append(int(j))
        if not placed:
            return
        sites = [self._pending[n0 + q] for q in range(len(placed))]
        self._pending = {}
        idx = np.array(placed)
        c = self.cells
        c.Y = np.vstack([c.Y, c.Y[idx]])
        c.weight = np.concatenate([c.weight, np.ones(len(idx))])
        c.alive = np.concatenate([c.alive, np.ones(len(idx), bool)])
        c.in_M = np.concatenate([c.in_M, np.zeros(len(idx), bool)])
        c.t_in_M = np.concatenate([c.t_in_M, np.zeros(len(idx))])
        c.arrest_h = np.concatenate([c.arrest_h, np.zeros(len(idx))])
        c.generation = np.concatenate([c.generation, c.generation[idx]])
        c.C_out = np.concatenate([c.C_out, c.C_out[idx]])
        c.forced_D = np.concatenate([c.forced_D, c.forced_D[idx]])
        new_bcl2, new_uptake = c.het_bcl2[idx].copy(), c.het_uptake[idx].copy()
        if self.dp.mutation_rate > 0 and self.dp.mutation_effect_sd > 0:
            # One daughter may differ from its mother; the other keeps her
            # values, so a division is at most one heritable change.
            sd = self.dp.mutation_effect_sd
            for arr in (new_bcl2, new_uptake):
                hit = self.rng.random(len(arr)) < self.dp.mutation_rate
                if hit.any():
                    arr[hit] *= np.exp(self.rng.normal(0.0, sd, int(hit.sum())))
                    self.n_mutations += int(hit.sum())
        c.het_bcl2 = np.concatenate([c.het_bcl2, new_bcl2])
        c.het_uptake = np.concatenate([c.het_uptake, new_uptake])
        if c.cycle_rate is not None:
            c.cycle_rate = np.concatenate([c.cycle_rate, c.cycle_rate[idx]])
        if c.quiescent is not None:
            # Each daughter decides independently, so the sister of a
            # quiescent cell is not automatically quiescent — which is the
            # observation the bifurcation was described from.
            from cellsim.cell.engine import _quiescent_fraction
            draw = self.rng.random(len(idx)) < _quiescent_fraction(self.line, self.p)
            c.quiescent = np.concatenate([c.quiescent, draw])
        self.lineage = np.concatenate([self.lineage, self.lineage[idx]])
        self.birth_t = np.concatenate([self.birth_t, np.full(len(idx), self.t_h)])
        self.pos = np.vstack([self.pos, np.array(sites)])
        self.o2 = np.concatenate([self.o2, self.o2[idx]])
        self.glucose = np.concatenate([self.glucose, self.glucose[idx]])
        self.gs = np.concatenate([self.gs, self.gs[idx]])

    def _bury(self, die: np.ndarray, kind: int) -> None:
        sites = tuple(self.pos[die].T)
        self.occ[sites] = kind
        if kind == APOPTOTIC:
            self.body_t[sites] = self.t_h

    def _compact(self) -> None:
        keep = self.cells.alive
        if keep.all():
            return
        c = self.cells
        for name in ("Y", "weight", "alive", "in_M", "t_in_M", "arrest_h", "generation",
                     "C_out", "forced_D", "het_bcl2", "het_uptake"):
            setattr(c, name, getattr(c, name)[keep])
        if c.cycle_rate is not None:
            c.cycle_rate = c.cycle_rate[keep]
        if c.quiescent is not None:
            c.quiescent = c.quiescent[keep]
        self.lineage, self.pos = self.lineage[keep], self.pos[keep]
        self.birth_t = self.birth_t[keep]
        self.o2, self.gs = self.o2[keep], self.gs[keep]
        self.glucose = self.glucose[keep]
        live = self.occ >= 0
        self.occ[live] = EMPTY
        self.occ[tuple(self.pos.T)] = np.arange(len(self.pos))

    def _necrosis(self, dt_h: float) -> None:
        dp = self.dp
        span = max(dp.o2_necrosis_onset_mmHg - dp.o2_necrosis_full_mmHg, 1e-9)
        frac = np.clip((dp.o2_necrosis_onset_mmHg - self.o2) / span, 0.0, 1.0)
        if dp.glucose_medium_mM > 0:
            # Starvation kills on the same clock as anoxia, and whichever
            # is scarcer decides. Freyer & Sutherland's central finding is
            # that a spheroid's necrotic core can be set by glucose rather
            # than oxygen, so neither can be the only route.
            gfrac = np.clip((dp.glucose_death_mM - self.glucose) /
                            max(dp.glucose_death_mM, 1e-9), 0.0, 1.0)
            frac = np.maximum(frac, gfrac)
        rate = dp.necrosis_rate_per_h * frac
        die = self.cells.alive & (self.rng.random(len(rate)) < 1.0 - np.exp(-rate * dt_h))
        if die.any():
            self.cells.alive[die] = False
            self._bury(die, NECROTIC)

    def _clear_bodies(self) -> None:
        old = (self.occ == APOPTOTIC) & (self.t_h - self.body_t >= self.dp.apoptotic_clearance_h)
        self.occ[old] = EMPTY

    # ── time stepping ──
    def advance(self, dose_fn: Optional[DoseFn]) -> None:
        """One field interval: fields, growth signal, engine steps with
        divisions and deaths, necrosis, clearance, compaction."""
        dp = self.dp
        dose = _dose_vector(dose_fn, self.t_h, self.n_drugs)
        self.update_fields(dose, dp.field_dt_h)
        self.update_growth_signal()
        n_sub = int(round(dp.field_dt_h / dp.dt_h))
        for _ in range(n_sub):
            _rk4(self.cells.Y, self._aux(), self.line, self.drugs or None, self.p, self.k_cyc, dp.dt_h)
            self.t_h += dp.dt_h
            divide, die = _cell_events(self.cells, self.drugs or None, self.p, dp.dt_h,
                                       rng=self.rng, line=self.line)
            if die.any():
                self._bury(die, APOPTOTIC)
            if divide.any():
                self._divide(divide & self.cells.alive)
        self._necrosis(dp.field_dt_h)
        self._clear_bodies()
        if len(self.eff_pos):
            self._effectors(dp.field_dt_h)
        if dp.migration_per_h > 0:
            # At most one hop per round, so a round must be short enough
            # that the chance of wanting two is small. One round per field
            # interval capped the effective rate at 1/interval: asked for
            # 2 hops/h it delivered 1.26, and the measured diffusion
            # coefficient came out 126 um^2/h against a theoretical 200.
            n_rounds = max(1, int(np.ceil(dp.migration_per_h * dp.field_dt_h / 0.1)))
            for _ in range(n_rounds):
                self._migrate(dp.field_dt_h / n_rounds)
        self._compact()

    def add_effectors(self, n: int) -> None:
        """Scatter `n` effector cells over the occupied part of the dish,
        as when T cells are pipetted onto a culture."""
        occupied = np.argwhere(self.occ != EMPTY)
        if len(occupied) == 0 or n <= 0:
            return
        pick = occupied[self.rng.integers(len(occupied), size=n)]
        self.eff_pos = np.concatenate([self.eff_pos, pick.astype(np.int64)])
        self.eff_alive = np.concatenate([self.eff_alive, np.ones(n, bool)])
        self.eff_kills = np.concatenate([self.eff_kills, np.zeros(n, dtype=np.int64)])
        self.eff_busy_until = np.concatenate([self.eff_busy_until, np.full(n, self.t_h)])
        self.eff_born = np.concatenate([self.eff_born, np.full(n, self.t_h)])

    def _effectors(self, dt_h: float) -> None:
        """Move, engage and kill, in sub-rounds short enough that a hop
        and a kill in the same round are rare."""
        dp = self.dp
        n_rounds = max(1, int(np.ceil(max(dp.effector_speed_per_h, dp.effector_kill_per_h)
                                      * dt_h / 0.2)))
        h = dt_h / n_rounds
        shape = np.array(self.shape)
        near = self._offsets[np.linalg.norm(self._offsets, axis=1) <= 1.0 + 1e-9]
        p_hop = 1.0 - np.exp(-dp.effector_speed_per_h * h)
        p_kill = 1.0 - np.exp(-dp.effector_kill_per_h * h)
        for r in range(n_rounds):
            t = self.t_h - dt_h + (r + 1) * h
            # culture lifetime
            self.eff_alive &= (t - self.eff_born) < dp.effector_lifetime_h
            active = self.eff_alive & (self.eff_kills < dp.effector_max_kills) & (
                self.eff_busy_until <= t)
            for e in np.flatnonzero(active):
                here = self.eff_pos[e]
                who = self.occ[tuple(here)]
                # a live tumour cell at this site: engage
                if who >= 0 and self.cells.alive[who]:
                    if self.rng.random() < p_kill:
                        die = np.zeros(len(self.cells.alive), bool)
                        die[who] = True
                        self.cells.alive[who] = False
                        self._bury(die, APOPTOTIC)
                        self.eff_kills[e] += 1
                        self.eff_busy_until[e] = t + dp.effector_conjugate_h
                        self.n_killed_by_effectors += 1
                    continue                      # stays conjugated, does not move
                if self.rng.random() >= p_hop:
                    continue
                cand = here + near
                cand = cand[np.all((cand >= 0) & (cand < shape), axis=1)]
                if len(cand) == 0:
                    continue
                if self.rng.random() < dp.effector_tumour_bias and len(self.cells.alive):
                    # toward a neighbouring live tumour cell if there is one.
                    # The length check matters: once the effectors have
                    # killed everything, `alive` is empty and indexing it
                    # even with a clamped index raises.
                    ids = self.occ[tuple(cand.T)]
                    live_there = (ids >= 0) & self.cells.alive[np.clip(ids, 0, len(self.cells.alive) - 1)]
                    targ = cand[live_there]
                    if len(targ):
                        self.eff_pos[e] = targ[self.rng.integers(len(targ))]
                        continue
                self.eff_pos[e] = cand[self.rng.integers(len(cand))]

    def _migrate(self, dt_h: float) -> None:
        """Each live cell with a free neighbouring site hops into one, with
        probability 1 - exp(-k dt). Moves are applied one cell at a time in
        a random order so two cells never claim the same site."""
        dp = self.dp
        p_hop = 1.0 - np.exp(-dp.migration_per_h * dt_h)
        live = np.flatnonzero(self.cells.alive)
        if len(live) == 0:
            return
        movers = live[self.rng.random(len(live)) < p_hop]
        if len(movers) == 0:
            return
        self.rng.shuffle(movers)
        near = self._offsets[np.linalg.norm(self._offsets, axis=1) <= 1.0 + 1e-9]
        shape = np.array(self.shape)
        for j in movers:
            here = self.pos[j]
            cand = here + near
            ok = np.all((cand >= 0) & (cand < shape), axis=1)
            cand = cand[ok]
            if len(cand) == 0:
                continue
            free = cand[self.occ[tuple(cand.T)] == EMPTY]
            if len(free) == 0:
                continue                      # boxed in: contact inhibition of motion
            if dp.migration_o2_bias > 0 and len(self.profile_o2) and self.rng.random() < dp.migration_o2_bias:
                # up the oxygen gradient = outward from the centre
                c = (shape - 1) / 2.0
                out = np.linalg.norm(free - c, axis=1)
                dst = free[int(np.argmax(out))]
            else:
                dst = free[self.rng.integers(len(free))]
            self.occ[tuple(here)] = EMPTY
            self.occ[tuple(dst)] = j
            self.pos[j] = dst

    def observe(self) -> None:
        dp, r = self.dp, self.record
        c = self.cells
        n = len(self.pos)
        r.t_h.append(round(self.t_h, 6))
        r.n_live.append(n)
        r.n_apoptotic.append(int((self.occ == APOPTOTIC).sum()))
        r.n_necrotic.append(int((self.occ == NECROTIC).sum()))
        ph = np.where(c.in_M, 3, _phase(c.Y)) if n else np.zeros(0, int)
        r.phase_frac.append([float(np.mean(ph == q)) if n else 0.0 for q in range(4)])
        r.quiescent_frac.append(float(np.mean((ph == 0) & (self.gs < 0.05))) if n else 0.0)
        r.radius_um.append(self.equivalent_radius_um(self.occ != EMPTY))
        r.necrotic_radius_um.append(self.equivalent_radius_um(self.occ == NECROTIC))
        r.hypoxic_frac.append(float(np.mean(self.o2 < dp.o2_g1_full_mmHg)) if n else 0.0)
        r.min_o2_mmHg.append(float(self.profile_o2.min()) if len(self.profile_o2)
                             else dp.o2_medium_mmHg)
        r.confluence.append(float(np.mean(self.occ != EMPTY)))
        r.blocked_divisions.append(self.blocked)
        r.mean_uptake.append(float(np.mean(c.het_uptake)) if n else float("nan"))
        r.mean_reserve.append(float(np.mean(c.het_bcl2)) if n else float("nan"))
        r.n_mutations.append(self.n_mutations)

    def snapshot(self) -> dict:
        """Plain-JSON state for a renderer (schema cellsim.dish.stream/v1)."""
        h = self.dp.spacing_um
        c = self.cells
        ph = np.where(c.in_M, 3, _phase(c.Y)) if len(self.pos) else np.zeros(0, int)
        bodies = np.argwhere(self.occ == APOPTOTIC)
        debris = np.argwhere(self.occ == NECROTIC)
        snap = {
            "t_h": round(self.t_h, 4),
            "cells": {
                "x_um": np.round(self.pos[:, 0] * h, 2).tolist(),
                "y_um": np.round(self.pos[:, 1] * h, 2).tolist(),
                "z_um": np.round(self.pos[:, 2] * h, 2).tolist(),
                "phase": ph.astype(int).tolist(),
                "p53": np.round(c.Y[:, IX["p53"]], 4).tolist(),
                "caspase3": np.round(c.Y[:, IX["C3"]], 4).tolist(),
                "damage": np.round(c.Y[:, IX["D"]], 4).tolist(),
                "o2_mmHg": np.round(self.o2, 2).tolist(),
                "growth_signal": np.round(self.gs, 3).tolist(),
                "lineage": self.lineage.astype(int).tolist(),
            },
            "apoptotic_um": np.round(bodies * h, 2).tolist(),
            "necrotic_um": np.round(debris * h, 2).tolist(),
        }
        if len(self.profile_r_um):
            snap["profile"] = {"r_um": np.round(self.profile_r_um, 2).tolist(),
                               "o2_mmHg": np.round(self.profile_o2, 3).tolist(),
                               "drug_uM": [np.round(pr, 6).tolist() for pr in self.profile_drug]}
        return snap


def _dose_vector(dose_fn: Optional[DoseFn], t_h: float, n: int) -> list:
    if dose_fn is None:
        return [0.0] * n
    v = dose_fn(t_h)
    if np.ndim(v) == 0:
        return [float(v)] * n if n == 1 else [float(v)] + [0.0] * (n - 1)
    v = [float(x) for x in v]
    if len(v) != n:
        raise ValueError(f"dose_fn returned {len(v)} values for {n} drugs")
    return v


@dataclass
class DishResult:
    record: DishRecord
    dish: Dish


def run_dish(line: CellLine, drug=None, dose_fn: Optional[DoseFn] = None, *,
             t_end_h: float = 72.0, n_seed: int = 100, grid_sites=None,
             record_every_h: float = 1.0, seed: int = 1, p: Params = Params(),
             dp: DishParams = DishParams(), layout: str = "auto",
             on_record: Optional[Callable[[Dish], None]] = None) -> DishResult:
    """Grow a dish (monolayer or spheroid) for `t_end_h` under `dose_fn`
    (µM per drug, as a function of time in hours; omit for no drug)."""
    dish = Dish(line, drug, n_seed=n_seed, grid_sites=grid_sites, seed=seed, p=p, dp=dp,
                layout=layout)
    stride = max(1, int(round(record_every_h / dp.field_dt_h)))
    n_steps = int(round(t_end_h / dp.field_dt_h))
    dish.update_fields(_dose_vector(dose_fn, 0.0, dish.n_drugs), 0.0)
    dish.update_growth_signal()
    dish.observe()
    if on_record is not None:
        on_record(dish)
    for k in range(1, n_steps + 1):
        dish.advance(dose_fn)
        if k % stride == 0 or k == n_steps:
            dish.observe()
            if on_record is not None:
                on_record(dish)
    return DishResult(dish.record, dish)


# ── state stream for a renderer ─────────────────────────────────────────
DISH_SCHEMA = "cellsim.dish.stream/v1"
DISH_CELL_FIELDS = ("x_um", "y_um", "z_um", "phase", "p53", "caspase3", "damage",
                    "o2_mmHg", "growth_signal", "lineage")


def run_dish_stream(line: CellLine, drug=None, schedule=None, *, geometry: str = "monolayer",
                    t_end_h: float = 72.0, n_seed: int = 100, grid_sites: Optional[int] = None,
                    record_every_h: float = 2.0, seed: int = 1, p: Params = Params(),
                    dp: Optional[DishParams] = None) -> list[dict]:
    """Run a dish and return [header, *snapshots] (schema cellsim.dish.stream/v1).

    Each snapshot carries every live cell's position and state, the
    apoptotic bodies and necrotic debris, the population readouts and,
    for a spheroid, the radial oxygen and drug profiles. Phases:
    0 = G1, 1 = S, 2 = G2, 3 = M. A renderer needs nothing else."""
    import dataclasses
    dp = dp or DishParams(geometry=geometry)
    dp = replace(dp, geometry=geometry)
    drugs = _as_drug_tuple(drug)
    dose_fn = schedule if schedule is not None else None
    snaps: list[dict] = []

    def on_record(d: Dish) -> None:
        r = d.record
        snap = d.snapshot()
        snap["drug_uM"] = _dose_vector(dose_fn, d.t_h, d.n_drugs)
        snap["population"] = {
            "n_live": r.n_live[-1], "n_apoptotic": r.n_apoptotic[-1],
            "n_necrotic": r.n_necrotic[-1],
            "fold_change": r.n_live[-1] / max(r.n_live[0], 1),
            "phase_frac": dict(zip(("G1", "S", "G2", "M"), r.phase_frac[-1])),
            "quiescent_frac": r.quiescent_frac[-1], "hypoxic_frac": r.hypoxic_frac[-1],
            "radius_um": r.radius_um[-1], "necrotic_radius_um": r.necrotic_radius_um[-1],
            "confluence": r.confluence[-1]}
        snaps.append(snap)

    res = run_dish(line, drug, dose_fn, t_end_h=t_end_h, n_seed=n_seed, grid_sites=grid_sites,
                   record_every_h=record_every_h, seed=seed, p=p, dp=dp, on_record=on_record)
    segs = getattr(schedule, "segments", None)
    header = {
        "schema": DISH_SCHEMA, "geometry": geometry, "line": line.name,
        "drug": [dg.name for dg in drugs] or None,
        "schedule": [list(s) for s in segs] if segs is not None else None,
        "t_end_h": t_end_h, "record_every_h": record_every_h, "n_seed": n_seed, "seed": seed,
        "grid_shape": list(res.dish.shape), "spacing_um": dp.spacing_um,
        "params": dataclasses.asdict(p), "dish_params": dataclasses.asdict(dp),
        "doubling_time_h": line.doubling_time_h, "p53_functional": line.p53_functional,
        "fields": list(DISH_CELL_FIELDS),
    }
    return [header, *snaps]


def main(argv: Optional[list[str]] = None) -> int:
    import argparse
    from pathlib import Path

    from cellsim.cell.library import get_drug, get_line
    from cellsim.cell.stream import Schedule, write_jsonl

    ap = argparse.ArgumentParser(prog="cellsim dish", description=(
        "Grow cells in space (a monolayer or a spheroid) under a dosing schedule and "
        "write the per-cell state stream (JSON Lines, schema " + DISH_SCHEMA + ")."))
    ap.add_argument("--line", default="A549")
    ap.add_argument("--geometry", choices=("monolayer", "spheroid"), default="monolayer")
    ap.add_argument("--drug", default="none", help="drug name, or 'none' for untreated")
    ap.add_argument("--schedule", default="0:0",
                    help="'start_h:conc_uM,...' piecewise constant, e.g. '24:10,48:0'")
    ap.add_argument("--hours", type=float, default=72.0)
    ap.add_argument("--cells", type=int, default=200, help="cells at the start")
    ap.add_argument("--grid", type=int, default=None, help="lattice side in sites")
    ap.add_argument("--every", type=float, default=2.0, help="record interval, hours")
    ap.add_argument("--cycle-cv", type=float, default=0.25,
                    help="cell-to-cell CV of cycle time (CTC HeLa: 0.22-0.33)")
    ap.add_argument("--push", type=float, default=None,
                    help="mechanical reach of a dividing cell, sites (default per geometry)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True, help="output .jsonl")
    a = ap.parse_args(argv)
    drug = None if a.drug.lower() == "none" else get_drug(a.drug)
    dp = DishParams(geometry=a.geometry)
    if a.push is not None:
        dp = replace(dp, push_sites=a.push)
    recs = run_dish_stream(get_line(a.line), drug, Schedule.parse(a.schedule),
                           geometry=a.geometry, t_end_h=a.hours, n_seed=a.cells,
                           grid_sites=a.grid, record_every_h=a.every, seed=a.seed,
                           p=Params(cycle_cv=a.cycle_cv), dp=dp)
    write_jsonl(recs, a.out)
    last = recs[-1]["population"]
    print(f"wrote {len(recs) - 1} snapshots to {a.out}: {last['n_live']} live cells "
          f"(x{last['fold_change']:.2f}), {last['n_apoptotic']} apoptotic, "
          f"{last['n_necrotic']} necrotic, radius {last['radius_um']:.0f} um")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
