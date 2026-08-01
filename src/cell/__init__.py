"""src/cell — Campaign-2 seed: cell-level readouts from Campaign-1 priors.

This is the FIRST consumer of the `src/bridge` rate-law priors — the
point where a molecular free energy becomes a cell-level quantity. It is
deliberately minimal (a single-compartment, single-site occupancy model)
so the end-to-end path exists and is honest before it is elaborated:

    binding ΔG  →  K_d prior (+ accuracy, + trust)  →  target occupancy
    (src/fep, src/dock)   (src/bridge)                 (here)

Everything carried through: the uncertainty on ΔG becomes a CI on
occupancy, and the target-class `trust` verdict rides along so a readout
computed from an untrustworthy K_d is never mistaken for a decision-grade
number.

Non-AI: closed-form Hill occupancy. No learned surrogate.

NOT yet here (later Campaign-2 work; see docs/campaign1_closeout.md):
  - multiple compartments (extracellular / cytoplasm) + membrane
    permeability, so occupancy is computed at the LOCAL concentration;
  - competitive / multi-drug binding at a shared site;
  - dynamics (ODE integration over time).
"""
from __future__ import annotations

from src.cell.occupancy import (
    OccupancyResult,
    occupancy_from_prior,
    hill_occupancy,
)
from src.cell.compartments import (
    Permeability,
    CellGeometry,
    CellOccupancyResult,
    spherical_cell_geometry,
    equilibration_tau_s,
    intracellular_concentration,
    occupancy_in_cell,
)
from src.cell.competition import (
    LigandOccupancy,
    CompetitiveResult,
    competitive_occupancy,
)
from src.cell.complexation import (
    ComplexationResult,
    complex_equilibrium,
)
from src.cell.dynamics import (
    RateConstants,
    TimeCourse,
    rate_constants_from_prior,
    simulate_binding_in_cell,
)
from src.cell.montecarlo import (
    MCResult,
    monte_carlo_propagate,
)
from src.cell.partitioning import (
    PartitionResult,
    COMPARTMENT_PH,
    neutral_fraction,
    neutral_fraction_polyprotic,
    accumulation_ratio,
    accumulation_ratio_polyprotic,
    partition_across_membrane,
)
from src.cell.transport import (
    EffluxPump,
    EffluxSteadyState,
    efflux_steady_state,
    pump_from_michaelis_prior,
)
from src.cell.binding_sink import (
    BindingSink,
    SinkResult,
    total_from_free,
    free_from_total,
)
from src.cell.steady_state import (
    IntracellularSteadyState,
    solve_steady_state,
)
from src.cell.disposition import (
    DispositionTimeCourse,
    simulate_disposition,
)
from src.cell.tissue import (
    TissueProfile,
    penetration_depth_um,
    penetration_profile_first_order,
    penetration_profile_saturable,
)
from src.cell.agents import (
    Clone,
    ColonyTrajectory,
    simulate_colony,
)
from src.cell.clearance import (
    HepaticClearance,
    ExposureProfile,
    RegimenOutcome,
    single_dose_exposure,
    repeat_dose_exposure,
    accumulation_ratio as dosing_accumulation_ratio,
    evaluate_regimen,
)
from src.cell.resistance import (
    ResistantClone,
    ResistanceOutcome,
    resistance_outcome,
)
from src.cell.fate import (
    CellFateParams,
    CellFate,
    TissueFateProfile,
    effect_from_occupancy,
    fate_from_occupancy,
    critical_occupancy,
    tissue_fate_profile,
)

__all__ = [
    "OccupancyResult", "occupancy_from_prior", "hill_occupancy",
    "Permeability", "CellGeometry", "CellOccupancyResult",
    "spherical_cell_geometry", "equilibration_tau_s",
    "intracellular_concentration", "occupancy_in_cell",
    "LigandOccupancy", "CompetitiveResult", "competitive_occupancy",
    "ComplexationResult", "complex_equilibrium",
    "RateConstants", "TimeCourse", "rate_constants_from_prior",
    "simulate_binding_in_cell",
    "MCResult", "monte_carlo_propagate",
    "PartitionResult", "COMPARTMENT_PH", "neutral_fraction",
    "neutral_fraction_polyprotic", "accumulation_ratio",
    "accumulation_ratio_polyprotic", "partition_across_membrane",
    "EffluxPump", "EffluxSteadyState", "efflux_steady_state",
    "pump_from_michaelis_prior",
    "BindingSink", "SinkResult", "total_from_free", "free_from_total",
    "IntracellularSteadyState", "solve_steady_state",
    "DispositionTimeCourse", "simulate_disposition",
    "TissueProfile", "penetration_depth_um",
    "penetration_profile_first_order", "penetration_profile_saturable",
    "CellFateParams", "CellFate", "TissueFateProfile",
    "effect_from_occupancy", "fate_from_occupancy", "critical_occupancy",
    "tissue_fate_profile",
    "ResistantClone", "ResistanceOutcome", "resistance_outcome",
    "HepaticClearance", "ExposureProfile", "RegimenOutcome",
    "single_dose_exposure", "repeat_dose_exposure",
    "dosing_accumulation_ratio", "evaluate_regimen",
    "Clone", "ColonyTrajectory", "simulate_colony",
]
