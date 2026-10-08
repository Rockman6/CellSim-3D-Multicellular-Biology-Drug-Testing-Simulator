"""cellsim.cell.compound — a drug the library does not have.

The library holds ten drugs, each with a mechanism, transport constants
from the literature and one potency constant fitted to GDSC. A drug typed
into the Lab by name, or drawn as a structure nobody has made, has none of
that. This module builds the best `Drug` the engine can run from what IS
known, and records where every number came from, so an invented molecule
can never pass for a validated prediction.

What can be known, and from where:

* **What it hits** — from the curated mechanism of action (ChEMBL,
  `classify`). Only mechanisms the engine's network expresses are
  accepted; anything else is refused with the target named, rather than
  run as something it is not.
* **How strongly** — NOT from structure. For a known class the potency is
  *borrowed* from the class's library drug, and the result is labelled
  uncalibrated. Within one class that can be off by orders of magnitude
  (gemcitabine and 5-fluorouracil are both antimetabolites and their
  fitted constants differ ~400-fold), which is why the Lab offers a
  strength multiplier and, given one measured IC50,
  `cellsim.calibrate.calibrate_potency` fits it properly.
* **How it gets in and out** — from structure, coarsely
  (`transport_from_properties`): whether it is likely a P-glycoprotein
  substrate (the "rule of fours", Didziapetris et al. 2003 J Drug Target
  11:391), and whether it is too polar to cross a membrane unaided
  (Veber et al. 2002 J Med Chem 45:2615). Uptake timing barely moves a
  72-hour outcome; P-gp status can, so only that one changes the run.
"""
from __future__ import annotations

import dataclasses
import re
from dataclasses import dataclass, field
from typing import Optional, Sequence

from cellsim.cell.library import DRUGS, Drug

__all__ = ["Classification", "classify", "transport_from_properties", "Transport",
           "from_mechanism", "MECHANISM_RULES", "designed", "ACTIONS", "DESIGNABLE"]

# Ordered: the first rule any mechanism text matches wins, so the specific
# ones come first (doxorubicin is curated as both "DNA inhibitor" and
# "DNA topoisomerase II alpha inhibitor"; it is a TopII poison).
# (pattern, engine mechanism, library drug the potency is borrowed from, words)
MECHANISM_RULES: tuple[tuple[str, str, str, str], ...] = (
    (r"topoisomerase\s+II\b|topoisomerase\s+2\b", "topo2", "etoposide",
     "poisons topoisomerase II, breaking DNA mostly while it is copied"),
    (r"topoisomerase\s+I\b|topoisomerase\s+1\b", "s_phase", "sn-38",
     "traps topoisomerase I, so replication forks break DNA"),
    (r"tubulin|microtubule", "tubulin", "paclitaxel",
     "binds tubulin and freezes the mitotic spindle"),
    (r"thymidylate synthase|ribonucleo\w*.?\s*reductase|DNA polymerase|"
     r"dihydrofolate reductase|DNA primase", "s_phase", "gemcitabine",
     "starves or corrupts DNA synthesis (an antimetabolite)"),
    (r"\bmdm2\b|p53.{0,20}mdm2|mdm2.{0,20}p53", "mdm2", "nutlin-3a",
     "blocks MDM2, so p53 is no longer destroyed"),
    (r"cyclin.dependent kinase [46]\b|\bCDK[46]\b", "cdk46", "palbociclib",
     "blocks cyclin D-CDK4/6, so the cell cannot pass the G1 restriction point"),
    (r"\bDNA (inhibitor|alkylating|cross-?linking)|alkylating agent|cross-?linking agent",
     "dna_adduct", "cisplatin", "binds DNA and damages it directly"),
)


@dataclass
class Classification:
    mechanism: str              # engine mechanism
    reference: str              # library drug whose potency is borrowed
    words: str                  # what it does, for a non-specialist
    matched: str                # the curated mechanism text that decided it


def classify(mechanism_texts: Sequence[str]) -> Optional[Classification]:
    """The engine mechanism for a set of curated mechanism-of-action texts,
    or None when none of them is something the engine models."""
    for pattern, mech, ref, words in MECHANISM_RULES:
        rx = re.compile(pattern, re.I)
        for text in mechanism_texts:
            if rx.search(text or ""):
                return Classification(mech, ref, words, text)
    return None


@dataclass
class Transport:
    pgp_substrate: Optional[bool]          # None = the rule cannot say
    permeable: Optional[bool]
    notes: list = field(default_factory=list)


def transport_from_properties(*, mw: float, tpsa: Optional[float] = None,
                              logp: Optional[float] = None,
                              n_plus_o: Optional[int] = None) -> Transport:
    """Coarse transport calls from structure, each with its reason.

    P-glycoprotein, rule of fours (Didziapetris et al. 2003): N+O >= 8 and
    MW > 400 marks a likely substrate; N+O <= 4 and MW < 400 a likely
    non-substrate; anything between is not called. On the library's ten
    it makes eight calls and declines SN-38 and gemcitabine. Seven agree
    with the library. The eighth shows what the rule cannot see: it calls
    nutlin-3a a substrate, and P-gp does transport it, but nutlin also
    blocks the pump at its active dose (Michaelis et al. 2009 Cancer Res
    69:416), so its net efflux is small and the library does not flag it.
    The rule predicts transport, not net efflux; a drug that is also a
    P-gp inhibitor will be over-called.

    Permeability: polar surface above 140 A^2 predicts poor passive
    permeation (Veber 2002). It is reported, not acted on — several
    anticancer drugs above it (paclitaxel 221, doxorubicin 206) still
    reach their targets within hours, and a 72-h outcome is insensitive
    to uptake time.
    """
    t = Transport(pgp_substrate=None, permeable=None)
    if n_plus_o is not None:
        if n_plus_o >= 8 and mw > 400:
            t.pgp_substrate = True
            t.notes.append(f"likely pumped out by P-glycoprotein (N+O = {n_plus_o}, "
                           f"MW {mw:.0f}: rule of fours)")
        elif n_plus_o <= 4 and mw < 400:
            t.pgp_substrate = False
            t.notes.append(f"unlikely to be pumped out by P-glycoprotein (N+O = {n_plus_o}, "
                           f"MW {mw:.0f})")
        else:
            t.notes.append("P-glycoprotein: the rule of fours cannot call this one")
    if tpsa is not None:
        t.permeable = tpsa <= 140
        if tpsa > 140:
            t.notes.append(f"very polar (surface {tpsa:.0f} A^2): would cross membranes "
                           f"slowly on its own; the engine assumes it gets in like its "
                           f"reference drug")
    if logp is not None and logp < -1:
        t.notes.append(f"very water-loving (logP {logp:.1f}): likely needs a transporter "
                       f"to enter cells, as gemcitabine does")
    return t


def from_mechanism(name: str, classification: Classification, *,
                   strength: float = 1.0, transport: Optional[Transport] = None) -> Drug:
    """A runnable Drug: the reference drug's constants, renamed, with its
    fitted potency multiplied by `strength` and its P-gp status replaced
    when the structure gives one. The `source` says all of that."""
    ref = DRUGS[classification.reference]
    changes = {"name": name}
    if ref.fit_target and strength != 1.0:
        changes[ref.fit_target] = getattr(ref, ref.fit_target) * strength
    if transport is not None and transport.pgp_substrate is not None:
        changes["pgp_substrate"] = transport.pgp_substrate
    changes["source"] = (
        f"UNCALIBRATED. Mechanism from the curated record ({classification.matched!r}). "
        f"Potency borrowed from {ref.name}"
        + (f" x{strength:g}" if strength != 1.0 else "")
        + ", whose constant was fitted to GDSC for that drug, not this one."
        + (" P-gp status from structure (rule of fours)." if "pgp_substrate" in changes else ""))
    changes["fit_target"] = ref.fit_target
    return dataclasses.replace(ref, **changes)


# ── design mode ───────────────────────────────────────────────────────
# What a designed drug can do. The first five reuse a library drug's
# whole mechanism (and its fitted potency, times `strength`); 'block' and
# 'boost' act on any one species through the engine's 'target' mechanism.
ACTIONS = {
    "damage_dna":   ("cisplatin",  "damages DNA directly, like cisplatin"),
    "damage_s":     ("sn-38",      "damages DNA only while it is copied, like SN-38"),
    "poison_top2":  ("etoposide",  "poisons topoisomerase II, like etoposide"),
    "freeze":       ("paclitaxel", "freezes the mitotic spindle, like paclitaxel"),
    "block_mdm2":   ("nutlin-3a",  "blocks MDM2 so p53 survives, like nutlin-3a"),
    "block_cdk46":  ("palbociclib", "blocks CDK4/6 so the cell stops in G1, like palbociclib"),
}
# Species a designed drug may block or boost: every one on the map except
# the drug slot and the damage index (damage has its own actions above).
DESIGNABLE = ("CycD", "Rb", "E2F", "CycE", "CycA", "CycB", "p21", "ATM", "p53",
              "MDM2_mRNA", "MDM2", "Puma", "Bax", "MOMP", "CytC", "Smac", "C9", "C3", "Exec")


def designed(name: str, action: str, *, species: str = "", kd_uM: float = 0.1,
             effect: float = 1.0, strength: float = 1.0,
             transport: Optional[Transport] = None) -> Drug:
    """A drug that does what the user says it does. Nothing about it is
    measured: the `source` says so, and the Lab shows it as a hypothesis.

    `action` is a key of ACTIONS, or 'block' / 'boost' with `species` one
    of DESIGNABLE, `kd_uM` its binding constant and `effect` the fraction
    blocked (block, 0-1) or the fold reached (boost, > 1) at full binding.
    """
    if action in ACTIONS:
        ref, words = ACTIONS[action]
        c = Classification(DRUGS[ref].mechanism, ref, words, f"designed: {action}")
        d = from_mechanism(name, c, strength=strength, transport=transport)
        return dataclasses.replace(d, source=f"DESIGNED, not measured: {words}; "
                                             f"strength {ref} x{strength:g}.")
    if action not in ("block", "boost"):
        raise ValueError(f"unknown action {action!r}; one of {sorted(ACTIONS) + ['block', 'boost']}")
    if species not in DESIGNABLE:
        raise ValueError(f"cannot {action} {species!r}; one of {DESIGNABLE}")
    if action == "block" and not 0.0 <= effect <= 1.0:
        raise ValueError("a blocker removes between 0 and all of its target (effect 0-1)")
    if action == "boost" and effect < 1.0:
        raise ValueError("a booster multiplies its target (effect >= 1)")
    if kd_uM <= 0:
        raise ValueError("kd_uM must be positive")
    pgp = bool(transport and transport.pgp_substrate)
    return Drug(name, "target", Kd_target_uM=kd_uM, target_species=species,
                target_mode=action, target_effect=effect, tau_uptake_h=0.3,
                partition=1.0, pgp_substrate=pgp,
                source=f"DESIGNED, not measured: {action}s {species} "
                       f"({effect:g}{'' if action == 'block' else '-fold'} at full binding, "
                       f"Kd {kd_uM:g} uM).")
