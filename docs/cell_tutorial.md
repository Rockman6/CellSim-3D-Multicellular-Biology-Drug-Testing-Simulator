# Tutorial: should I pulse this drug, or leave it on?

*For someone who has not used CellSim before. Twenty minutes, a laptop,
no conda. Every command here is run and its real output pasted in.*

This is the question CellSim is built for. It is a question about one
cell line compared with itself under different schedules, which is where
mechanism does the work — and it is deliberately **not** the question
"which of my cell lines is most sensitive", which CellSim refuses to
answer because that turned out to be unreachable from the available
markers ([why](VALIDATION.md)).

## Install

The cell simulator needs only NumPy and SciPy:

```bash
pip install cellsim
```

The molecular layers (docking, MD, quantum) need the conda environment in
`environment.yml`; nothing in this tutorial does.

## 1. Where does this drug act at all?

Start with a dose-response curve. `curve()` returns a tidy table — a
`pandas.DataFrame` if pandas is installed, a list of dicts otherwise.

```python
from cellsim.api import curve

curve("A549", "paclitaxel")
```

```
 line       drug  conc_uM  viability  ic50_uM  exposure_h  readout_h
 A549 paclitaxel 0.000271   0.951613 0.032681        72.0       72.0
 A549 paclitaxel 0.001256   1.064516 0.032681        72.0       72.0
 A549 paclitaxel 0.012558   0.935484 0.032681        72.0       72.0
 A549 paclitaxel 0.027055   0.661290 0.032681        72.0       72.0
 A549 paclitaxel 0.058288   0.006048 0.032681        72.0       72.0
 A549 paclitaxel 0.270551   0.002016 0.032681        72.0       72.0
 A549 paclitaxel 2.705512   0.000000 0.032681        72.0       72.0
```

(Abridged: the call returns all thirteen concentrations.) The IC50 is
about 0.033 µM. Notice how steep it is — 0.66 viability at 0.027 µM and
0.006 at 0.058. Once tubulin occupancy crosses the arrest threshold,
every dividing cell stalls, so the curve is much steeper than a Hill
slope of 1. Concentrations default to four decades either
side of it, so you do not have to guess the range.

> **How much is this number worth?** The engine's IC50s carry a
> calibrated interval: ±2.6× covers 64 % of held-out lines, ±4.8× covers
> 85 %. That is wide, and honestly so — GDSC's own repeat screens of the
> same line and drug disagree by a median 5.3×. Treat a predicted IC50 as
> an order of magnitude, not a measurement.

## 2. Ask the schedule question properly

The tempting experiment is to fix a concentration and vary the hours. The
informative one is the **iso-effect curve**: for each exposure time, what
concentration is needed to halve the colony?

```python
from cellsim.api import exposure

exposure("A549", "paclitaxel", hours=[1, 3, 6, 12, 24, 48, 72])
```

```
 exposure_h  readout_h   c50_uM  conc_x_hours_uM_h  reaches_half_kill
        1.0       72.0      inf                inf              False
        3.0       72.0      inf                inf              False
        6.0       72.0      inf                inf              False
       12.0       72.0      inf                inf              False
       24.0       72.0 0.084238           2.021707               True
       48.0       72.0 0.028141           1.350790               True
       72.0       72.0 0.029696           2.138079               True
```

Read the first column: **below 24 hours there is no concentration that
halves the colony.** Not a high one, not any one. Paclitaxel only kills
cells that attempt mitosis while it is present, so a short pulse misses
most of an asynchronous population however concentrated it is.

Now run the same thing for a platinum:

```python
exposure("A549", "cisplatin", hours=[1, 3, 6, 12, 24, 48, 72])
```

```
 exposure_h     c50_uM  conc_x_hours_uM_h  reaches_half_kill
        1.0 318.253177         318.253177               True
        3.0 108.958853         326.876560               True
        6.0  56.296001         337.776003               True
       12.0  27.748253         332.979031               True
       24.0  16.775299         402.607180               True
       48.0  12.296506         590.232301               True
       72.0  10.000000         720.000000               True
```

Every exposure works, and `conc_x_hours` stays within 318–403 µM·h from
1 to 24 h — the concentration-times-time law measured for cisplatin
(Ozawa et al. 1989). Halve the time, double the dose, same kill. Beyond
24 h it drifts upward because repair starts to keep pace with the slow
accumulation of adducts.

**That contrast is the answer to the question.** A platinum can be
pulsed: what matters is total exposure. A taxane cannot: it needs time
above a threshold, and compressing the exposure throws the drug away.

## 3. Check what happens after you wash it out

```python
from cellsim.api import washout

washout("A549", "cisplatin", conc_uM=40.0, exposure_h=24.0, hours=168.0)
```

```
   t_h  drug_present  fold_change  control_fold_change  relative_to_control
   0.0          True     1.000000             1.000000             1.000000
  24.0         False     0.687500             2.166667             0.317308
  48.0         False     0.208333             4.666667             0.044643
  96.0         False     0.437500            21.166667             0.020669
 168.0         False     5.250000           190.666667             0.027535
```

The colony is cut to a fifth of its starting size by 48 h — a day *after*
the drug came off, because committed cells keep dying — and then regrows,
reaching 5.25× its starting size by day seven. Relative to the untreated
control it still looks like 3 % survival, because the control ran away
to 191×.

Both facts matter and they disagree about what happened: the culture was
not eradicated, it was set back about four days. If your assay reports
only the ratio to control, you cannot tell those apart. Try 60 µM for the
same 24 h and the colony ends at 0.63× — below where it started.

## 4. Two drugs: check the readout, not just the result

```python
from cellsim.api import combination

combination("A549", "cisplatin", "paclitaxel", 10.0, 0.03, readout_h=72.0)
```

```
                      arm  surviving_fraction  surviving_sd  ratio_to_independent  ratio_sd
          cisplatin alone            0.769236      0.022913                   NaN       NaN
         paclitaxel alone            0.465439      0.022066                   NaN       NaN
             simultaneous            0.422007      0.023502              1.179152  0.040188
cisplatin then paclitaxel            0.436547      0.015645              1.220167  0.021519
paclitaxel then cisplatin            0.447851      0.019187              1.251612  0.028100
```

`ratio_to_independent` above 1 is antagonism. All three arrangements kill
*less* than the two drugs acting independently would, because cisplatin
arrests the cycle through p53 → p21 and removes the mitotic cells
paclitaxel needs.

Note `surviving_sd`: these are averaged over seeds, and the spread is
±0.02. The gap between the two orderings (0.437 against 0.448) is about
half of that, so **there is no sequence effect here** — and if you ran
this once and saw those two numbers you might report one. Lengthen the readout
and even the apparent difference shrinks, because the slow drug placed
second simply had less time to act.

## 5. Into three dimensions

Monolayers see the drug you add. A spheroid does not.

```python
from cellsim.api import spheroid

spheroid("DLD-1", days=10, n_seed=4000, grid_sites=56, record_every_h=48.0)
```

```
 day  n_live  n_necrotic  radius_um  necrotic_radius_um  hypoxic_fraction  min_o2_mmHg
 0.0    4000           0 147.711753            0.000000          0.000000    57.219622
 2.0    5737           0 166.580185            0.000000          0.000000    40.556088
 4.0    7928           0 185.545148            0.000000          0.143416    20.961440
 6.0   10572           0 204.227447            0.000000          0.340049     3.525840
 8.0   13564         203 223.018810           54.688213          0.518947     0.778873
10.0   16808         990 242.951055           92.741359          0.625535     0.794362
```

Read it as a sequence. Oxygen at the centre falls steadily; a hypoxic
fraction appears on day 4, while every cell is still alive; the centre
goes anoxic on day 6; and only on day 8 does a necrotic core appear. By
day 10 nearly two thirds of the cells are hypoxic and the core is 93 µm
across. Hypoxia precedes death by days — which is why a drug that needs
dividing cells loses potency in an aggregate long before anything in it
has died. The oxygen numbers are DLD-1's own measured ones
(Grimes et al. 2014); the one spatial constant is calibrated against that
study's growth rate.

To see it, write a stream and open it in the viewer:

```bash
cellsim dish --line DLD-1 --geometry spheroid --hours 240 \
    --cells 4000 --grid 56 --every 8 --out spheroid.jsonl
```

Open `web/viewer/index.html` and press **Open .jsonl**. Tick *cut open*
to slice the spheroid and colour by *oxygen* to see the gradient.

## What to be careful about

- **Two known misses**, both in [VALIDATION.md](VALIDATION.md): a 24 h
  paclitaxel exposure is far more effective in the engine than in
  published data, and simulated necrotic cores are larger than the
  measured anoxic radius implies. If your question sits on either, read
  that section first.
- **No drug retention after wash-out.** The engine empties a cell as fast
  as it fills it, so very short exposures are modelled as doing nothing.
- **Absolute IC50s are order-of-magnitude.** The schedule *comparisons*
  above are the trustworthy part, because both arms share whatever the
  constant gets wrong.
- **Do not use it to rank your cell lines.** Measured as unreachable;
  there is deliberately no function for it.

## Where to go next

- [`VALIDATION.md`](VALIDATION.md) — every claim, its reproducer, and its
  caveat, including two results withdrawn when re-measured.
- `scripts/experiment_*.py` — the full sweeps behind the summary numbers.
- [`PLAN.md`](PLAN.md) — what is coming, most relevantly calibrating the
  engine to *your* plate-reader data rather than to ours.
