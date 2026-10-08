#!/usr/bin/env python3
"""Gates for drugs the library does not have (`cellsim.cell.compound`).

The mechanism texts below are ChEMBL's curated records, copied verbatim
in October 2026, so these run offline. Three things are checked:

* every library drug's own curated mechanism maps back to the engine
  mechanism the library gives it — the rules are tested against the ten
  cases where the right answer is already known;
* a target the engine does not model (EGFR, PARP) is refused, not
  quietly run as something else;
* a borrowed-potency drug runs, says it is uncalibrated, and behaves like
  its reference drug when the multiplier is 1.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _tiers import run  # noqa: E402
from cellsim.cell.compound import classify, from_mechanism, transport_from_properties  # noqa: E402
from cellsim.cell.engine import simulate  # noqa: E402
from cellsim.cell.library import DRUGS, get_line  # noqa: E402

# ChEMBL mechanism_of_action, by parent molecule (October 2026).
CHEMBL = {
    "cisplatin": ["DNA inhibitor"],
    "doxorubicin": ["DNA inhibitor", "DNA topoisomerase II alpha inhibitor"],
    "etoposide": ["DNA topoisomerase II inhibitor"],
    "sn-38": ["DNA topoisomerase I inhibitor"],
    "gemcitabine": ["DNA polymerase (alpha/delta/epsilon) inhibitor",
                    "Ribonucleoside-diphosphate reductase RR1 inhibitor"],
    "5-fluorouracil": ["DNA inhibitor", "RNA inhibitor", "Thymidylate synthase inhibitor"],
    "paclitaxel": ["Tubulin inhibitor"],
    "docetaxel": ["Tubulin inhibitor", "Tubulin stabiliser"],
    "vinorelbine": ["Tubulin inhibitor"],
    "idasanutlin": ["Tumour suppressor p53/oncoprotein Mdm2 inhibitor"],
    "palbociclib": ["CDK6/cyclin D1 inhibitor", "Cyclin-dependent kinase 4/cyclin D1 inhibitor"],
}


def test_library_drugs_classify_to_their_own_mechanism():
    wrong = []
    for name, texts in CHEMBL.items():
        c = classify(texts)
        want = DRUGS["nutlin-3a" if name == "idasanutlin" else name].mechanism
        if c is None or c.mechanism != want:
            wrong.append(f"{name}: got {c and c.mechanism}, library says {want}")
    # 5-FU is curated as "DNA inhibitor" AND "Thymidylate synthase inhibitor";
    # the library models it as S-phase damage, so the specific rule must win.
    assert not wrong, "\n".join(wrong)


def test_unmodelled_targets_are_refused():
    for texts in (["Epidermal growth factor receptor erbB1 inhibitor"],
                  ["PARP 1, 2 and 3 inhibitor"],
                  ["Apoptosis regulator Bcl-2 inhibitor"],
                  []):
        assert classify(texts) is None, f"{texts} should not map onto the engine"


def test_rule_of_fours_on_known_substrates():
    # formulas from PubChem; N+O counted from them
    assert transport_from_properties(mw=853.9, n_plus_o=15).pgp_substrate is True    # paclitaxel
    assert transport_from_properties(mw=543.5, n_plus_o=12).pgp_substrate is True    # doxorubicin
    assert transport_from_properties(mw=300.1, n_plus_o=2).pgp_substrate is False    # cisplatin
    assert transport_from_properties(mw=263.2, n_plus_o=7).pgp_substrate is None     # gemcitabine
    polar = transport_from_properties(mw=263.2, tpsa=150.0, logp=-1.5, n_plus_o=7)
    assert polar.permeable is False and len(polar.notes) == 3


def test_borrowed_drug_runs_and_says_it_is_uncalibrated():
    c = classify(["DNA topoisomerase I inhibitor"])
    d = from_mechanism("topotecan", c)
    assert d.name == "topotecan" and d.mechanism == "s_phase"
    assert d.source.startswith("UNCALIBRATED") and "sn-38" in d.source
    line = get_line("A549")
    same = simulate(line, d, 0.01, t_end_h=48, n_cells=32, seed=3)
    ref = simulate(line, DRUGS["sn-38"], 0.01, t_end_h=48, n_cells=32, seed=3)
    assert abs(same.alive_weight[-1, 0] - ref.alive_weight[-1, 0]) < 1e-9
    # the multiplier must act on the potency, whatever the reference's value
    weak, strong = (simulate(line, from_mechanism("x", c, strength=s), 1.0,
                             t_end_h=48, n_cells=32, seed=3) for s in (0.01, 100.0))
    assert strong.alive_weight[-1, 0] < weak.alive_weight[-1, 0], \
        (strong.alive_weight[-1, 0], weak.alive_weight[-1, 0])


if __name__ == "__main__":
    raise SystemExit(run(globals(), "compound"))
