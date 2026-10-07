/* CellSim Lab — the real engine, running in the browser.
 *
 * Pyodide loads CPython with numpy and scipy, then micropip installs the
 * published `cellsim` wheel from PyPI. Nothing is simulated twice and
 * nothing is faked: the dose-response below is the same `dose_response`
 * that `pip install cellsim` gives you at a terminal.
 *
 * Nothing is uploaded. The page is static files; the computation happens
 * on the viewer's own machine.
 */
const $ = (id) => document.getElementById(id);
const state = { py: null, lines: [], drugs: [], genes: [], busy: false };

function setBusy(on, label) {
  state.busy = on;
  $("run").disabled = on || !state.py;
  $("bar").classList.toggle("indet", on);
  $("bar").querySelector("i").style.width = on ? "" : "0";
  if (label) $("engine-state").textContent = label;
}

/* ── boot ──────────────────────────────────────────────────────── */
async function boot() {
  setBusy(true, "downloading Python…");
  state.py = await loadPyodide({
    indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.2/full/",
  });
  $("engine-state").textContent = "loading numpy and scipy…";
  await state.py.loadPackage(["numpy", "scipy", "micropip"]);
  $("engine-state").textContent = "installing cellsim…";
  await state.py.runPythonAsync(`
import micropip
await micropip.install("cellsim")
`);
  const meta = JSON.parse(await state.py.runPythonAsync(`
import json
from cellsim.cell.library import CELL_LINES, DRUGS
from cellsim.cell.perturb import GENES
json.dumps({
  "lines": [{"name": n, "tissue": c.tissue, "p53": c.p53_functional,
             "doubling": c.doubling_time_h} for n, c in CELL_LINES.items()],
  "drugs": [{"name": n, "mechanism": d.mechanism} for n, d in DRUGS.items()],
  "genes": [{"name": g, "note": v["note"]} for g, v in GENES.items()],
})
`));
  Object.assign(state, meta);
  fillMenus();
  setBusy(false, "engine ready");
  $("run").disabled = false;
}

const MECHANISM = {
  dna_adduct: "crosslinks DNA; damage accumulates whatever the cell is doing",
  topo2: "poisons topoisomerase II; worst during DNA replication",
  tubulin: "freezes the spindle, so it only bites cells entering mitosis",
  s_phase: "only damages cells that are actively copying their DNA",
  mdm2: "frees p53 — works only where p53 is intact",
};

function fillMenus() {
  const ls = $("line");
  state.lines.forEach((l) => ls.add(new Option(`${l.name} — ${l.tissue}`, l.name)));
  const ds = $("drug");
  state.drugs.forEach((d) => ds.add(new Option(d.name, d.name)));
  const gs = $("gene");
  state.genes.forEach((g) => gs.add(new Option(g.name, g.name)));
  ls.value = "A549"; ds.value = "cisplatin";
  describeLine(); describeDrug(); describeGene();
}

function describeLine() {
  const l = state.lines.find((x) => x.name === $("line").value);
  $("line-hint").innerHTML = l
    ? `${l.tissue}. p53 is <b>${l.p53 ? "intact" : "mutated"}</b>, and one cell
       takes about ${l.doubling} h to divide.`
    : "";
}
function describeDrug() {
  const d = state.drugs.find((x) => x.name === $("drug").value);
  $("drug-hint").textContent = d ? (MECHANISM[d.mechanism] || d.mechanism) : "";
}
function describeGene() {
  const g = state.genes.find((x) => x.name === $("gene").value);
  $("gene-hint").textContent = g ? g.note : "Leave this alone for normal cells.";
}

/* ── run ───────────────────────────────────────────────────────── */
async function run() {
  if (state.busy) return;
  setBusy(true, "simulating…");
  $("placeholder").hidden = true;
  const cfg = {
    line: $("line").value,
    drug: $("drug").value,
    gene: $("gene").value || null,
    hours: Number($("hours").value),
    exposure: $("exposure").value ? Number($("exposure").value) : null,
  };
  try {
    const out = JSON.parse(await state.py.runPythonAsync(`
import json, numpy as np
from cellsim.cell.engine import Params, calibrate_cycle_scale, dose_response, ic50
from cellsim.cell.library import get_drug, get_line
cfg = ${JSON.stringify(JSON.stringify(cfg))}
cfg = json.loads(cfg)
line, drug = get_line(cfg["line"]), get_drug(cfg["drug"])
pert = None
if cfg["gene"]:
    from cellsim.cell.perturb import Perturbation
    pert = [Perturbation(cfg["gene"])]
p = Params()
k = calibrate_cycle_scale(line, p)
guess = {"paclitaxel":0.03,"docetaxel":0.03,"vinorelbine":0.05,"doxorubicin":0.05,
         "etoposide":1.0,"sn-38":0.03,"gemcitabine":0.05,"5-fluorouracil":5.0,
         "nutlin-3a":5.0}.get(drug.name, 10.0)
centre = ic50(line, drug, guess_uM=guess, n_cells_per_conc=16, k_cyc=k, p=p,
              n_seeds=1, exposure_h=cfg["exposure"], perturbations=pert)
if not np.isfinite(centre):
    centre = guess * 50
conc = centre * np.logspace(-1.5, 1.5, 8)
viab, _ = dose_response(line, drug, conc, t_end_h=cfg["hours"], n_cells_per_conc=24,
                        k_cyc=k, p=p, exposure_h=cfg["exposure"], perturbations=pert)
# the band: repeat under different seeds, because the engine is stochastic
runs = []
for sd in (1, 2, 3):
    v, _ = dose_response(line, drug, conc, t_end_h=cfg["hours"], n_cells_per_conc=24,
                         seed=sd, k_cyc=k, p=p, exposure_h=cfg["exposure"],
                         perturbations=pert)
    runs.append(v[1:])
runs = np.array(runs)
ics = [ic50(line, drug, guess_uM=guess, n_cells_per_conc=16, k_cyc=k, p=p, seed=sd,
            n_seeds=1, exposure_h=cfg["exposure"], perturbations=pert) for sd in (1,2,3)]
fin = [x for x in ics if np.isfinite(x)]
json.dumps({
  "conc": conc.tolist(),
  "mid": np.median(runs, axis=0).tolist(),
  "lo": runs.min(axis=0).tolist(),
  "hi": runs.max(axis=0).tolist(),
  "ic50": float(np.median(fin)) if fin else None,
  "ic50_lo": float(min(fin)) if fin else None,
  "ic50_hi": float(max(fin)) if fin else None,
  "reached": len(fin) > 1,
})
`));
    draw(out, cfg);
    setBusy(false, "engine ready");
  } catch (e) {
    setBusy(false, "engine ready");
    $("result").innerHTML =
      `<div class="note warn"><b>That run did not complete.</b><br>${String(e).slice(0, 400)}</div>`;
  }
}

/* ── results ───────────────────────────────────────────────────── */
function draw(o, cfg) {
  const band = {
    x: o.conc.concat([...o.conc].reverse()),
    y: o.hi.concat([...o.lo].reverse()),
    fill: "toself", fillcolor: "rgba(53,179,162,.16)",
    line: { width: 0 }, hoverinfo: "skip", showlegend: false, type: "scatter",
  };
  const mid = {
    x: o.conc, y: o.mid, mode: "lines+markers", type: "scatter",
    line: { color: "#35b3a2", width: 2 },
    marker: { size: 7, color: "#35b3a2" },
    name: "surviving fraction",
    hovertemplate: "%{x:.3g} µM → %{y:.0%}<extra></extra>",
  };
  Plotly.newPlot("stage", [band, mid], {
    paper_bgcolor: "#090d10", plot_bgcolor: "#090d10",
    font: { color: "#8d9cab", family: "ui-monospace, Menlo, monospace", size: 11 },
    margin: { l: 58, r: 22, t: 46, b: 52 },
    title: {
      text: `${cfg.line}${cfg.gene ? " · " + cfg.gene + " knocked out" : ""} · ${cfg.drug}`,
      font: { color: "#e4ebf2", size: 13.5 }, x: 0.02,
    },
    xaxis: {
      title: "concentration (µM)", type: "log",
      gridcolor: "#1d262e", zerolinecolor: "#2a3540",
    },
    yaxis: {
      title: "fraction surviving", range: [0, 1.08], tickformat: ".0%",
      gridcolor: "#1d262e", zerolinecolor: "#2a3540",
    },
    showlegend: false,
  }, { displayModeBar: false, responsive: true });

  const pct = (v) => (v * 100).toFixed(0) + "%";
  const ic = o.ic50;
  $("result").innerHTML = `
    <div class="metric">
      <div class="k">Concentration that halves the colony</div>
      <div class="v">${ic ? ic.toPrecision(3) + " µM" : "not reached"}</div>
      ${ic ? `<div class="band">range over repeats
         ${o.ic50_lo.toPrecision(3)} – ${o.ic50_hi.toPrecision(3)} µM</div>` : ""}
    </div>
    ${ic ? "" : `<div class="note warn">No concentration tested halved the
       colony. For a spindle poison that is the expected answer below about
       24 h of exposure — the drug can only act on cells that reach mitosis
       while it is present.</div>`}
    <div class="metric">
      <div class="k">At the highest dose tested</div>
      <div class="v">${pct(o.mid[o.mid.length - 1])}</div>
      <div class="band">still alive</div>
    </div>
    <div class="note">
      The shaded band is three repeats of the same experiment. The engine
      draws a finite sample of cell lineages and survival is set by the
      resistant few, so repeats move by 10–20 %. <b>Two numbers closer
      than the band are not different.</b>
    </div>
    <table class="data">
      <tr><th>µM</th><th>alive</th><th>range</th></tr>
      ${o.conc.map((c, i) => `<tr>
        <td>${c.toPrecision(3)}</td><td>${pct(o.mid[i])}</td>
        <td>${pct(o.lo[i])}–${pct(o.hi[i])}</td></tr>`).join("")}
    </table>`;
  $("overlay").hidden = false;
  $("overlay").innerHTML =
    `<b>${cfg.line}</b> · ${cfg.drug} · ${cfg.hours} h` +
    (cfg.exposure ? ` · drug present ${cfg.exposure} h` : " · drug present throughout");
}

/* ── presets ───────────────────────────────────────────────────── */
const PRESETS = {
  pulse:     { line: "A549", drug: "paclitaxel", exposure: "12", hours: 72, gene: "" },
  p53:       { line: "A549", drug: "cisplatin", exposure: "", hours: 72, gene: "TP53" },
  resistant: { line: "DLD-1", drug: "doxorubicin", exposure: "", hours: 72, gene: "" },
  immune:    { line: "A549", drug: "cisplatin", exposure: "", hours: 72, gene: "BCL2" },
};

document.addEventListener("DOMContentLoaded", () => {
  $("hours").oninput = () => ($("hours-val").textContent = $("hours").value + " h");
  $("line").onchange = describeLine;
  $("drug").onchange = describeDrug;
  $("gene").onchange = describeGene;
  $("run").onclick = run;
  document.querySelectorAll(".preset").forEach((b) => {
    b.onclick = () => {
      const p = PRESETS[b.dataset.preset];
      if (!p) return;
      $("line").value = p.line; $("drug").value = p.drug;
      $("gene").value = p.gene; $("exposure").value = p.exposure;
      $("hours").value = p.hours; $("hours-val").textContent = p.hours + " h";
      describeLine(); describeDrug(); describeGene();
      run();
    };
  });
  boot().catch((e) => {
    setBusy(false, "engine failed to start");
    $("placeholder").innerHTML =
      `<h2>The engine could not start</h2><div>${String(e).slice(0, 300)}</div>`;
  });
});
