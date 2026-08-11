# Provenance audit — Campaign 1 decision constants

**Audit date:** 2026-07
**Scope:** every hardcoded numeric constant that influences a *decision*
or a *predicted value* in the Campaign-1 scoring / triage / strain /
ADMET / metabolism / uncertainty logic.
**Standard (MISSION.md):** *"A grep for the value should land on a
literature citation, a physics derivation, or a labelled project
convention — never a naked magic number."*

**Method:** grep the decision modules for numeric literals, trace each
to its source, classify.

## Verdict

**Mostly compliant, with a corrected false claim and a small number of
labelled conventions.** The physics and cheminformatics constants are
properly cited. The exceptions are (a) the scope doc's claim that
provenance is *enforced* by a `src/uq/Prediction` class — that class
does not exist; corrected — and (b) two med-chem triage cutoffs that
are project policy, now explicitly labelled as such rather than implied
to be standards.

## What carries provenance today

There is **no unified `Prediction` envelope** (the scope doc claimed one;
it was never built). Instead each result dataclass carries provenance
*data*:

| Envelope | Provenance carried |
|---|---|
| `DockingResult` (`src/dock/vina.py`) | `tool_versions` (vina + meeko), seed, exhaustiveness, box, receptor+ligand hashes, `cache_key()` |
| `SoMResult` (`src/quantum/metabolism.py`) | `tool_versions` (xtb) |
| `BindingDGResult` (`src/fep/binding.py`) | force-field IDs, seed, window schedule, restraint params |

So provenance is *attached* but not *enforced by a single type*. A
unified `Prediction` type is future work (tracked; not claimed as
shipped).

## Decision-constant register

### Cited — literature or physics (compliant)

| Constant | Value | Where | Source |
|---|---|---|---|
| logP (Crippen) | — | `chem/admet.py` | Wildman & Crippen 1999, JCICS 39:868 |
| TPSA | — | `chem/admet.py` | Ertl, Rohde, Selzer 2000, JMC 43:3714 |
| Rule-of-5 (MW≤500, logP≤5, HBA≤10, HBD≤5) | as listed | `chem/admet.py` | Lipinski et al 1997, Adv Drug Deliv Rev 23:3 |
| QED | — | `chem/admet.py` | Bickerton et al 2012, Nat Chem 4:90 |
| ESOL solubility coefficients | — | `chem/admet.py` | Delaney 2004, JCICS 44:1000 |
| BBB rule-of-three | — | `chem/admet.py` | Pardridge 2003, J Neurochem 87:1 |
| Strain bands | 1.5 / 3.0 / 7.0 | `dock/strain.py` | Perola & Charifson 2004, JMC 47:2499; Buttenschoen 2024, Chem Sci 15:3130 |
| Kd from ΔG | Kd = exp(ΔG/RT) | `bridge/__init__.py` | thermodynamics (physics derivation) |
| FEP standard-state correction | +kT·ln(V_std/V_harm) | `fep/binding.py` | Gilson/Boresch; validated vs openmmtools SSC |
| Reliability trust bar | 1.5 / 3.0 kcal/mol | `benchmarks/dock/reliability_table.yaml` | RT·ln(10)=1.36 kcal/mol = 1 log Kd; FEP success target, Wang et al 2015, JACS 137:2695 *(added + cited this audit)* |

### Labelled project conventions (compliant — transparent, not implied-standard)

| Constant | Value | Where | Status |
|---|---|---|---|
| `follow_up` ΔG cutoff | −7.3 kcal/mol (~IC50 4 µM) | `dock/batch.py` | project triage policy; rationale in-comment; **labelled a convention this audit** |
| `drop` ΔG floor | −6.0 kcal/mol (~IC50 40 µM) | `dock/batch.py` | project triage policy; labelled |
| Kd human-readable band edges | 1 nM / 1 µM / … | `dock/batch.py` | display labels, not decisions |

### Corrected this audit

- Scope-doc criterion 6 claimed provenance *enforced by* `src/uq/
  Prediction` — **no such class exists**. Corrected to describe the
  real per-envelope provenance data + note the unified type as future
  work.
- The reliability thresholds added earlier this session (1.5 / 3.0
  kcal/mol) were uncited — **now cited** (RT·ln(10) + Wang 2015).

## Remaining gaps (honest)

1. **No enforcement.** Nothing *prevents* a future contributor from
   adding a naked magic number. A `Prediction` envelope type (or a CI
   lint that flags un-commented numeric literals in decision modules)
   would make the discipline enforced rather than convention. Not built.
2. **The triage cutoffs are policy, not physics.** They are honest and
   documented, but a different program would set them differently. They
   should be surfaced as tunable parameters, not baked constants — minor
   refactor, not done.
