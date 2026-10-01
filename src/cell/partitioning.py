"""pH partitioning / ion trapping — why C_in ≠ C_out at steady state.

The passive permeation model (`compartments.py`) equilibrates the drug to
C_in = C_out because it treats the drug as a single neutral species. Real
ionisable drugs don't: only the NEUTRAL form crosses the membrane freely,
so at steady state it is the neutral concentration that equalises across
the membrane — while the TOTAL concentration differs between compartments
whenever their pH differs. A weak base is more protonated (charged,
trapped) in an acidic compartment, so it ACCUMULATES there; a weak acid
accumulates in a basic compartment. This is the mechanism behind
lysosomotropic drugs (chloroquine, many antihistamines).

Physics (exact Henderson–Hasselbalch, monoprotic):

    base  BH⁺ ⇌ B + H⁺ :  f_neutral(pH) = 1 / (1 + 10^(pKa − pH))
    acid  HA  ⇌ A⁻ + H⁺:  f_neutral(pH) = 1 / (1 + 10^(pH − pKa))

At steady state [B]_in = [B]_out (neutral equalises), so

    accumulation R = total_in / total_out = f_neutral(pH_out) / f_neutral(pH_in)

R > 1 means the drug concentrates inside. R is exact given (pKa, pH_in,
pH_out, ion type); it is passive (no active transport) and monoprotic
(diprotic bases such as chloroquine trap even more strongly — deferred).

Provenance (never a magic number): pKa is a molecular property (input;
would come from PROPKA / DimorphiteDL, see `src/chem`). Compartment pH
values are cited physiological constants below.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict

# Cited physiological pH values. Sources:
#  * blood/plasma 7.4, cytosol ~7.2 — Roos & Boron 1981, Physiol Rev 61:296
#    (cytosol also BioNumbers BNID 105980).
#  * lysosome ~4.7 — Ohkuma & Poole 1978, PNAS 75:3327.
#  * mitochondrial matrix ~8.0 — Llopis et al 1998, PNAS 95:6803.
#  * early endosome ~6.3 — Maxfield & Yamashiro 1987.
COMPARTMENT_PH: Dict[str, float] = {
    "blood": 7.4,
    "cytosol": 7.2,
    "lysosome": 4.7,
    "early_endosome": 6.3,
    "mitochondrion": 8.0,
}

_VALID_IONS = ("base", "acid", "neutral")


def neutral_fraction(pH: float, pKa: float, ion_type: str) -> float:
    """Fraction of a monoprotic drug in the membrane-permeant neutral form.

    At pH = pKa the drug is exactly half ionised (f_neutral = 0.5).
    A neutral drug is always fully permeant (f_neutral = 1).
    """
    if ion_type not in _VALID_IONS:
        raise ValueError(f"ion_type must be one of {_VALID_IONS}")
    if ion_type == "neutral":
        return 1.0
    if ion_type == "base":
        exponent = pKa - pH        # more acidic (low pH) → more protonated → less neutral
    else:                          # acid
        exponent = pH - pKa        # more basic (high pH) → more deprotonated → less neutral
    return 1.0 / (1.0 + 10.0 ** exponent)


def neutral_fraction_polyprotic(pH: float, pKas, ion_type: str) -> float:
    """Neutral fraction for a POLYPROTIC drug with several ionisable sites.

    Chloroquine-class diprotic bases trap far harder than the monoprotic
    model allows, because each protonation site multiplies the number of
    charged states available in an acidic compartment.

    For a base with pKa₁…pKa_n the successive-protonation partition
    function (relative to the neutral form) is

        Z = 1 + Σ_k Π_{i≤k} 10^(pKa_i − pH)      →   f_neutral = 1/Z

    and for an acid the exponents flip sign (successive deprotonations).
    Reduces exactly to `neutral_fraction` for a single pKa.

    HONESTY LIMIT — this is an UPPER BOUND on passive trapping. For
    chloroquine (pKa 10.2, 8.4) it predicts ~2.3e5× lysosomal
    accumulation, while measured values are ~1e3-1e4×. The gap is real
    physics deliberately omitted here: an accumulating base RAISES
    lysosomal pH, shrinking the gradient that drives the trapping
    (self-limiting), and finite drug / lysosomal volume / binding cap it
    further. Treat large polyprotic ratios as "traps strongly, bounded
    above by this number", not as a calibrated prediction.
    """
    if ion_type not in _VALID_IONS:
        raise ValueError(f"ion_type must be one of {_VALID_IONS}")
    if ion_type == "neutral":
        return 1.0
    ks = [float(p) for p in pKas]
    if not ks:
        raise ValueError("need at least one pKa")
    # Protonate strongest site first (descending) for a base; for an acid
    # deprotonate the strongest acid (lowest pKa) first (ascending).
    ks.sort(reverse=(ion_type == "base"))
    Z = 1.0
    term = 1.0
    for pKa in ks:
        exponent = (pKa - pH) if ion_type == "base" else (pH - pKa)
        term *= 10.0 ** exponent
        Z += term
    return 1.0 / Z


def accumulation_ratio_polyprotic(pKas, pH_in: float, pH_out: float,
                                  ion_type: str) -> float:
    """Steady-state total_in/total_out for a polyprotic drug."""
    fn_in = neutral_fraction_polyprotic(pH_in, pKas, ion_type)
    fn_out = neutral_fraction_polyprotic(pH_out, pKas, ion_type)
    if fn_in <= 0.0:
        return math.inf
    return fn_out / fn_in


@dataclass
class PartitionResult:
    """Steady-state accumulation of an ionisable drug across a membrane."""

    accumulation_ratio: float          # total_in / total_out
    total_in_M: float
    total_out_M: float
    f_neutral_in: float
    f_neutral_out: float
    pKa: float
    pH_in: float
    pH_out: float
    ion_type: str
    pKa_source: str = "input"

    def summary(self) -> str:
        direction = ("accumulates INSIDE" if self.accumulation_ratio > 1.0
                     else "excluded from inside" if self.accumulation_ratio < 1.0
                     else "no net gradient")
        return (f"{self.ion_type} pKa={self.pKa:g}: "
                f"C_in/C_out = {self.accumulation_ratio:.2f}× "
                f"(pH {self.pH_out}→{self.pH_in}) — {direction} "
                f"[total_in={self.total_in_M:.2e} M]")


def accumulation_ratio(pKa: float, pH_in: float, pH_out: float,
                       ion_type: str) -> float:
    """Steady-state total_in/total_out from Henderson-Hasselbalch.

    R = f_neutral(pH_out) / f_neutral(pH_in): the neutral form equalises
    across the membrane, so the compartment where the drug is MORE ionised
    holds more total drug.
    """
    fn_in = neutral_fraction(pH_in, pKa, ion_type)
    fn_out = neutral_fraction(pH_out, pKa, ion_type)
    if fn_in <= 0.0:
        return math.inf
    return fn_out / fn_in


def partition_across_membrane(
    C_out_M: float, pKa: float, *,
    pH_in: float, pH_out: float = 7.4, ion_type: str = "base",
    pKa_source: str = "input",
) -> PartitionResult:
    """Total intracellular concentration of an ionisable drug at steady state.

    Args:
        C_out_M: extracellular TOTAL concentration.
        pKa: the drug's (monoprotic) pKa.
        pH_in: intracellular / target-compartment pH (use COMPARTMENT_PH).
        pH_out: extracellular pH (default blood 7.4).
        ion_type: 'base' | 'acid' | 'neutral'.

    The resulting `total_in_M` is what an intracellular target actually
    sees — feed it to the occupancy model in place of a naive C_out.
    """
    if C_out_M < 0:
        raise ValueError("C_out_M must be ≥ 0")
    R = accumulation_ratio(pKa, pH_in, pH_out, ion_type)
    return PartitionResult(
        accumulation_ratio=R,
        total_in_M=R * C_out_M,
        total_out_M=C_out_M,
        f_neutral_in=neutral_fraction(pH_in, pKa, ion_type),
        f_neutral_out=neutral_fraction(pH_out, pKa, ion_type),
        pKa=pKa, pH_in=pH_in, pH_out=pH_out, ion_type=ion_type,
        pKa_source=pKa_source,
    )
