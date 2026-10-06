"""cellsim.cell.perturb — knock a gene out, or push one up.

The engine already carries the pathway as state variables: p53, MDM2,
p21, PUMA, BAX, the caspases, the cyclins. A genetic perturbation is
therefore not a new mechanism but a constraint on one that is already
there — a knockout holds a node at zero, a knockdown scales its
production, an overexpression adds a constant source.

    Perturbation("TP53", mode="knockout")          # p53 never accumulates
    Perturbation("BCL2", mode="overexpress", x=4)  # a bigger reserve
    Perturbation("CDKN1A", mode="knockdown", x=0.2)  # 20 % of normal p21

That framing is what makes these checkable rather than decorative. A
TP53 knockout must reproduce what the library's `p53_functional=False`
lines already do; a CDKN1A knockout must remove the arrest that p21
causes and nothing else. `tests/cell/test_perturb_smoke.py` holds both.

GENE NAMES map to engine nodes, not the other way round: `TP53` is the
`p53` variable, `CDKN1A` is `p21`, `BCL2` is the line's anti-apoptotic
reserve. Where a gene has no counterpart in the engine the perturbation
is REFUSED rather than silently ignored, because a knockout that does
nothing is worse than an error — it looks like a result.

What this is NOT: a model of transcription, of knockout efficiency, or
of the compensation a real cell mounts over days. It is "this node is
absent/reduced/raised, what does the rest of the network do", which is
the question a pathway perturbation usually asks.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Union

import numpy as np

__all__ = ["Perturbation", "GENES", "apply_to_params", "describe"]

# Gene symbol -> what it is in the engine. Each entry says which node it
# controls and how a perturbation reaches it, because they are not
# uniform: some are state variables held at a level, others are rate
# constants in Params, and one is a property of the cell line.
GENES: dict[str, dict] = {
    "TP53":   {"node": "p53",  "kind": "state",
               "note": "p53 itself; a knockout also removes transactivation of "
                       "MDM2, p21 and PUMA, which is what a hotspot mutant does"},
    "MDM2":   {"node": "MDM2", "kind": "state",
               "note": "the ligase that degrades p53; losing it stabilises p53"},
    "CDKN1A": {"node": "p21",  "kind": "state",
               "note": "p21, the CDK inhibitor that enforces p53's arrest"},
    "BBC3":   {"node": "Puma", "kind": "state",
               "note": "PUMA, the BH3-only protein p53 induces to start apoptosis"},
    "BAX":    {"node": "Bax",  "kind": "state",
               "note": "the pore-former; without it MOMP cannot happen"},
    "CASP3":  {"node": "C3",   "kind": "state",
               "note": "executioner caspase; the engine's commitment point"},
    "CASP9":  {"node": "C9",   "kind": "state", "note": "apical caspase of the intrinsic route"},
    "ATM":    {"node": "ATM",  "kind": "state", "note": "the damage sensor that activates p53"},
    "BCL2":   {"node": "bcl2_level", "kind": "line",
               "note": "the anti-apoptotic reserve BAX must exceed, a per-line property"},
    "ABCB1":  {"node": "efflux_level", "kind": "line",
               "note": "P-glycoprotein; raising it pumps substrate drugs back out"},
    "CCNE1":  {"node": "CycE", "kind": "state", "note": "cyclin E, the G1/S switch"},
    "CCNB1":  {"node": "CycB", "kind": "state", "note": "cyclin B; without it cells cannot enter mitosis"},
    "RB1":    {"node": "Rb",   "kind": "state", "note": "the restriction-point brake"},
}

_MODES = ("knockout", "knockdown", "overexpress")


@dataclass(frozen=True)
class Perturbation:
    """One gene, perturbed one way.

    `x` is the multiplier for knockdown (0-1) and overexpression (>1); it
    is ignored for a knockout, which is a knockdown to zero.
    """
    gene: str
    mode: str = "knockout"
    x: float = 1.0

    def __post_init__(self):
        g = self.gene.upper()
        if g not in GENES:
            raise ValueError(
                f"{self.gene!r} has no counterpart in this engine. Known genes: "
                f"{', '.join(sorted(GENES))}. A perturbation of a gene the model "
                f"does not represent would do nothing, which looks like a result "
                f"and is worse than this error.")
        if self.mode not in _MODES:
            raise ValueError(f"mode must be one of {_MODES}, not {self.mode!r}")
        if self.mode == "knockdown" and not (0.0 <= self.x < 1.0):
            raise ValueError(f"a knockdown needs 0 <= x < 1, got {self.x}")
        if self.mode == "overexpress" and self.x <= 1.0:
            raise ValueError(f"an overexpression needs x > 1, got {self.x}")
        object.__setattr__(self, "gene", g)

    @property
    def factor(self) -> float:
        """What the node is multiplied by."""
        return {"knockout": 0.0, "knockdown": self.x, "overexpress": self.x}[self.mode]

    @property
    def info(self) -> dict:
        return GENES[self.gene]

    def __str__(self) -> str:
        if self.mode == "knockout":
            return f"{self.gene} knockout"
        return f"{self.gene} {self.mode} x{self.x:g}"


def describe(perturbations: Optional[Sequence[Perturbation]]) -> str:
    if not perturbations:
        return "unperturbed"
    return "; ".join(str(p) for p in perturbations)


def apply_to_params(line, p, perturbations: Optional[Sequence[Perturbation]]):
    """Return (line, params, state_factors) with the perturbations applied.

    Line-level genes (BCL2, ABCB1) change the `CellLine`; node-level genes
    become a multiplier applied to that state variable every step, which
    the engine's `rhs` enforces via `aux["state_factor"]`.
    """
    import dataclasses

    if not perturbations:
        return line, p, None
    from cellsim.cell.engine import IX, NV

    factors = np.ones(NV)
    line_changes: dict = {}
    for pert in perturbations:
        info = pert.info
        if info["kind"] == "line":
            field = info["node"]
            line_changes[field] = getattr(line, field) * pert.factor
        else:
            factors[IX[info["node"]]] *= pert.factor
    if line_changes:
        line = dataclasses.replace(line, **line_changes)
    # TP53 knocked out means no transactivation either, exactly as a
    # loss-of-function mutant: the engine already models that through
    # CellLine.p53_functional, so reuse it rather than inventing a second
    # route to the same place.
    if any(q.gene == "TP53" and q.factor == 0.0 for q in perturbations):
        line = dataclasses.replace(line, p53_functional=False)
    return line, p, (factors if not np.all(factors == 1.0) else None)
