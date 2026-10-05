---
title: 'CellSim: a mechanistic simulator of cultured cells under drug, with its errors measured'
tags:
  - Python
  - systems biology
  - pharmacology
  - cancer
  - agent-based modelling
authors:
  - name: Henry
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 4 October 2026
bibliography: paper.bib
---

# Summary

`CellSim` simulates a population of cultured cancer cells responding to a
drug, from the molecular events inside each cell to the colony that grows
or dies in the dish. Each cell carries a cell-cycle oscillator, the
ATM/p53/MDM2 damage response, and the intrinsic apoptotic cascade, driven
by the intracellular concentration of one or more drugs. The same cells
can be run well-mixed, where a whole dose-response curve is a single
simulation, or on a lattice, where oxygen and drug diffuse in from the
medium, cells compete for space, and a spheroid develops a hypoxic rim
and a necrotic core.

The package answers questions *within* a cell line, which is where
mechanism carries the explanatory weight: whether to pulse a drug or hold
it, how long after wash-out a culture recovers, whether two drugs
antagonise, what a schedule selects for, and how deep a drug penetrates
an aggregate. Every result is a time course of individual cells, not a
fitted curve, and every simulation writes a JSON Lines stream that the
bundled browser viewer renders without running any biology of its own.

The cell simulator depends only on NumPy and SciPy and installs with
`pip`; a dose-response takes seconds on a laptop and a simulated week of
spheroid growth takes minutes.

# Statement of need

Tools for this question tend to sit at one of two extremes. Agent-based
frameworks such as PhysiCell [@Ghaffarizadeh:2018] and Chaste
[@Cooper:2020] model space and population dynamics in depth, but treat
the inside of a cell phenomenologically, so a drug enters as a death rate
rather than as a mechanism. Systems-biology models of the p53 network
[@GevaZatorsky:2006; @Purvis:2012] or of apoptosis [@Albeck:2008] capture
the molecular detail, but describe one cell in a well-stirred bath and
are not connected to a measurable colony outcome. A researcher who wants
to ask "should I pulse this drug?" has to bridge that gap themselves.

`CellSim` connects the two, and reports how far the connection can be
trusted. Its validation is the part we would ask others to copy:

* **The fitted quantities are named and counted.** One potency constant
  per drug, fitted on two cell lines against GDSC [@Yang:2013]; one
  spatial constant, calibrated against published spheroid growth
  [@Grimes:2014]. Everything else is a literature input cited at the
  point of use.
* **Predictions are compared with public measurements, and the misses are
  published next to the passes.** Re-running the Cell Tracking Challenge
  HeLa time-lapse [@Ulman:2017] reproduces colony size at 46 h but not
  the spread of single-cell cycle times; the spheroid's necrotic core is
  larger than the measured anoxic radius implies. Both are recorded.
* **A negative result shaped the scope.** Predicting which cell line is
  more sensitive turned out to be unreachable from canonical markers: an
  optimal linear fit of eleven mechanism-chosen markers to IC50 across
  156–683 GDSC lines explains 9–20 % of the variance. Per-line ranking is
  therefore stated as out of scope rather than attempted.
* **Claims are gated over seeds and doses.** Two results published during
  development were withdrawn when re-measured — one rested on a single
  surviving lineage, the other on one random seed — and the test suite
  now pins the corrected versions.

The result is a simulator whose answers come with the conditions under
which they were checked, which is what makes a model usable by someone
who did not build it.

# Validated behaviour

Exposure-time dependence is reproduced from mechanism rather than
fitted: the concentration needed for half kill falls as $1/T$ for
cisplatin and doxorubicin, the concentration-time law measured for
cell-cycle phase-non-specific agents [@Ozawa:1989], while paclitaxel
cannot halve an asynchronous population below 24 h of exposure at any
concentration, matching the plateau and the exposure dependence reported
for taxanes [@Liebmann:1993; @Georgiadis:1997]. A cytostatic paired with
a mitosis-specific partner antagonises it, and the apparent advantage of
one drug order over the other is shown to shrink as the readout
lengthens — a timing confound of fixed-endpoint assays rather than a
biological sequence effect.

In the lattice, oxygen solved against a packed aggregate reproduces the
233 µm diffusion limit measured for DLD-1 spheroids [@Grimes:2014], and
heritable change at division lets resistance arise during treatment:
starting from an identical population, where selection has nothing to act
on, continuous exposure leaves survivors with roughly twice the parental
IC50, while the same total exposure given as brief high-dose pulses
eradicates the culture before resistance appears.

# Acknowledgements

`CellSim` builds on public data from GDSC, DepMap, the Cell Tracking
Challenge and the PDB, and on the open-source scientific Python stack.

# References
