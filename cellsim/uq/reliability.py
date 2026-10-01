"""Target-class reliability lookup — attach an ACCURACY estimate to a
prediction, not just a reproducibility estimate.

The uncertainty CellSim reported to users came from Monte-Carlo over
Vina seeds. That is PRECISION ("does Vina agree with itself?"), not
ACCURACY ("is Vina right?"). The two diverge badly: on
biotin/streptavidin the seed-scatter bar was ±0.29 kcal/mol while the
true error vs experiment was 11.82 — the displayed bar understated the
real error ~41×. Vina reproducibly returns the same wrong number, and a
tight bar on a wrong number invites exactly the trust it doesn't merit.

This module reads `benchmarks/dock/reliability_table.yaml` (measured
absolute error per target class) and answers, for a given receptor:

    "how far off is CellSim's dG likely to be here, and should the
     absolute number be trusted at all?"

Design rule — NEVER GUESS. A receptor that is not in the table is
reported as `calibrated=False` with `mean_abs_err_kcalmol=None` and a
verdict of `uncalibrated`. We do not borrow a neighbouring class's
number: a new kinase is not guaranteed to behave like EGFR, and
inventing an error bar is worse than admitting we lack one.

Non-AI: this is a lookup of measured experimental residuals. No model.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional
import logging
import re

logger = logging.getLogger(__name__)

_DEFAULT_TABLE = (Path(__file__).resolve().parents[2]
                  / "benchmarks" / "dock" / "reliability_table.yaml")

# Human-readable one-liners for each verdict.
_VERDICT_TEXT = {
    "trustworthy_absolute":
        "absolute ΔG usable for decisions",
    "rank_order_only":
        "use ordering only — absolute ΔG unreliable",
    "do_not_trust_absolute":
        "absolute ΔG NOT meaningful for this target class",
    "uncalibrated":
        "target not calibrated — accuracy UNKNOWN; treat ΔG as "
        "rank-order only",
}


@dataclass
class Reliability:
    """Accuracy estimate for a prediction on a given receptor."""

    calibrated: bool
    target_class: Optional[str] = None
    verdict: str = "uncalibrated"
    mean_abs_err_kcalmol: Optional[float] = None
    worst_abs_err_kcalmol: Optional[float] = None
    n: int = 0
    note: str = ""
    receptor_key: Optional[str] = None

    @property
    def verdict_text(self) -> str:
        return _VERDICT_TEXT.get(self.verdict, self.verdict)

    def as_dict(self) -> dict:
        d = asdict(self)
        d["verdict_text"] = self.verdict_text
        return d

    def summary(self) -> str:
        """One line a biologist can act on."""
        if not self.calibrated:
            return ("accuracy: UNKNOWN (uncalibrated target) — "
                    "treat ΔG as rank-order only")
        return (f"accuracy: ±{self.mean_abs_err_kcalmol:.1f} kcal/mol "
                f"typical (worst {self.worst_abs_err_kcalmol:.1f}, "
                f"n={self.n}, class={self.target_class}) — "
                f"{self.verdict_text}")


def _receptor_key(receptor: str | Path) -> str:
    """Normalise a receptor reference to a bare 4-char PDB-ish id.

    Accepts a path ('benchmarks/dock/1stp.pdb'), a filename ('1STP.pdb')
    or a bare id ('1stp'). Returns lowercase stem.
    """
    stem = Path(str(receptor)).name
    stem = re.sub(r"\.(pdb|pdbqt|cif|ent)(\.gz)?$", "", stem,
                  flags=re.IGNORECASE)
    return stem.strip().lower()


def load_table(path: str | Path | None = None) -> dict:
    """Load the reliability table YAML. Returns {} if unavailable."""
    p = Path(path) if path else _DEFAULT_TABLE
    if not p.exists():
        logger.warning("reliability table not found at %s", p)
        return {}
    try:
        import yaml
        return yaml.safe_load(p.read_text()) or {}
    except Exception as e:  # noqa: BLE001
        logger.warning("could not read reliability table: %s", e)
        return {}


def reliability_for(receptor: str | Path,
                    table_path: str | Path | None = None) -> Reliability:
    """Measured accuracy for `receptor`, or an explicit 'uncalibrated'.

    Matching is by PDB id against `receptor_pdb_ids` in the table. No
    fuzzy/nearest-class fallback by design — see module docstring.
    """
    key = _receptor_key(receptor)
    data = load_table(table_path)
    for entry in (data.get("classes") or []):
        ids = [str(i).lower() for i in (entry.get("receptor_pdb_ids") or [])]
        if key in ids:
            return Reliability(
                calibrated=True,
                target_class=entry.get("target_class"),
                verdict=entry.get("verdict", "uncalibrated"),
                mean_abs_err_kcalmol=entry.get("mean_abs_err_kcalmol"),
                worst_abs_err_kcalmol=entry.get("worst_abs_err_kcalmol"),
                n=int(entry.get("n") or 0),
                note=(entry.get("note") or "").strip(),
                receptor_key=key,
            )
    return Reliability(calibrated=False, receptor_key=key)


def annotate_result(result, receptor: str | Path | None = None,
                    table_path: str | Path | None = None):
    """Attach `.reliability` to a DockingResult (or any object).

    Uses `result.receptor_source` when `receptor` isn't given. Returns
    the same object for chaining. Never raises — a missing table must
    not break a docking run.
    """
    try:
        ref = receptor if receptor is not None else getattr(
            result, "receptor_source", None)
        if ref is None:
            return result
        result.reliability = reliability_for(ref, table_path)
    except Exception as e:  # noqa: BLE001
        logger.debug("reliability annotation failed: %s", e)
    return result


__all__ = ["Reliability", "reliability_for", "annotate_result", "load_table"]
