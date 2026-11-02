// The lab journal: tasks, what you observed on each slide, the growth curves
// and dose–response the CellSim engine produced, and your own results.
import { SPECIMENS, STAGES } from "./specimens.js";
import { WELLS } from "./props.js";
import { fmtConc } from "./items.js";

const KEY = "cellsim-lab:journal:v1";
const FACTS = {
  onion: { kind: "plant", shape: "Long bricks in neat rows", wall: "Yes", chloro: "None", vacuole: "Yes, large" },
  elodea: { kind: "plant", shape: "Small bricks, full of green dots", wall: "Yes", chloro: "Yes, many", vacuole: "Yes, large" },
  cheek: { kind: "animal", shape: "Flat and irregular", wall: "No", chloro: "None", vacuole: "No" },
  redonion: { kind: "plant", shape: "Long bricks, red sap", wall: "Yes", chloro: "None", vacuole: "Yes, red" },
  root: { kind: "plant", shape: "Small boxes, big nuclei", wall: "Yes", chloro: "None", vacuole: "Small" },
};

export class Journal {
  constructor() {
    try { this.data = JSON.parse(localStorage.getItem(KEY)) || {}; } catch (e) { this.data = {}; }
    this.data.obs = this.data.obs || {};
    this.data.results = this.data.results || {};
  }
  save() { try { localStorage.setItem(KEY, JSON.stringify(this.data)); } catch (e) { /* private window */ } }
  observe(spec, patch) {
    const d = (this.data.obs[spec] = this.data.obs[spec] || {});
    if (patch.mag) d.mag = Math.max(d.mag || 0, patch.mag);
    for (const k of ["stain", "cover", "size"]) if (patch[k] != null) d[k] = patch[k];
    this.save();
  }
  result(key, value) { this.data.results[key] = { ...value, at: Date.now() }; this.save(); }
  render(tab, tasks, dishes, plate) {
    if (tab === "tasks") {
      const v = tasks.view();
      if (!v) return `<p class="note">Free roam: no task list. Pick an experiment from the title screen to get one.</p>`;
      return `<h3>${v.title}</h3><ol class="checklist">${v.steps.map((s) => `<li class="${s.done ? "done" : ""}">${s.text}</li>`).join("")}</ol>${v.finished ? `<p class="note">Experiment complete.</p>` : ""}`;
    }
    if (tab === "obs") {
      const order = ["onion", "redonion", "elodea", "root", "cheek"], o = this.data.obs;
      const seen = (s) => (o[s]?.mag || 0) >= 100;
      const rows = [
        ["Stain", (s) => o[s]?.stain || "—"],
        ["Highest power in focus", (s) => (o[s]?.mag ? `${o[s].mag}×` : "—")],
        ["Shape", (s) => (seen(s) ? FACTS[s].shape : "—")],
        ["Cell wall", (s) => (seen(s) ? FACTS[s].wall : "—")],
        ["Chloroplasts", (s) => (seen(s) ? FACTS[s].chloro : "—")],
        ["Large vacuole", (s) => (seen(s) ? FACTS[s].vacuole : "—")],
        ["Your size estimate", (s) => (o[s]?.size ? `about ${o[s].size} µm` : "—")],
        ["Real size", (s) => SPECIMENS[s].sizeText],
      ];
      return `<div class="table-wrap"><table class="jt"><thead><tr><th></th>${order.map((s) => `<th>${SPECIMENS[s].name}<br><span class="muted">${FACTS[s].kind} cell</span></th>`).join("")}</tr></thead>
        <tbody>${rows.map(([l, f]) => `<tr><th>${l}</th>${order.map((s) => `<td>${f(s)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>
        ${this.mitosisHTML()}${this.osmosisHTML()}
        <p class="note">Rows fill in once you've seen a sample in focus at 100× or more. These slides are drawn from measured cell sizes and stage proportions; they are not simulated by CellSim.</p>`;
    }
    // Culture: dose–response, cell counts, and live cells over time, all from the engine.
    const parts = [this.doseHTML(plate), this.countHTML()];
    const withData = dishes.filter((d) => d.state.series.length);
    if (withData.length) {
      parts.push(`<h3>Growth</h3>` + chart(withData) + withData.map((d) => {
        const last = d.state.series[d.state.series.length - 1];
        const drugs = d.state.events.length ? d.state.events.map((e) => `${e.t} h: ${Object.entries(e.conc).map(([k, v]) => `${k} ${fmtConc(v)}`).join(", ")}`).join(" · ") : "no drug";
        return `<p><b>${d.name}</b>: ${last.n_live} live cells at ${d.state.t} h (started with ${d.state.series[0].n_live}) · ${last.n_apoptotic} apoptotic · treatment: ${drugs}</p>`;
      }).join(""));
    }
    const body = parts.filter(Boolean).join("");
    if (!body) return `<p class="note">No culture data yet. Look at a dish under the microscope, incubate one, or read the 6-well plate, and the CellSim engine's numbers appear here.</p>`;
    return body + `<p class="note">Every number here comes from cellsim.cell.dish, running in your browser.</p>`;
  }

  mitosisHTML() {
    const m = this.data.results.mitosis;
    if (!m || !m.total) return "";
    const dividing = m.total - m.counts[0];
    return `<h3>Mitosis in the root tip</h3><p>You classified <b>${m.total}</b> cells: ${STAGES.map((s, i) => `${s.toLowerCase()} ${m.counts[i]}`).join(", ")}.
      Mitotic index = dividing ÷ all = ${dividing} ÷ ${m.total} = <b>${(dividing / m.total * 100).toFixed(1)} %</b>. You got ${m.right} of ${m.total} right first time.</p>
      <p class="note">In a growing onion root tip about 10–15 % of cells near the tip are dividing. Prophase is the longest stage, so it is the one you see most.</p>`;
  }
  osmosisHTML() {
    const o = this.data.results.osmosis;
    if (!o || !o.log?.length) return "";
    return `<h3>Osmosis in red onion cells</h3><ul>${o.log.map((e) => `<li>${e}</li>`).join("")}</ul>
      <p class="note">Water moves out of the cells into a more concentrated solution, so the protoplast shrinks from the wall (plasmolysis); in pure water it moves back in. The wall stops the cell bursting.</p>`;
  }
  countHTML() {
    const c = this.data.results.count;
    if (!c) return "";
    return `<h3>Cell count</h3><p>${c.text}</p>`;
  }
  doseHTML(plate) {
    if (!plate || !plate.state.read) return "";
    const st = plate.state, r = st.read;
    const rows = r.wells.map((w, i) => `<tr><th>${WELLS[i]}</th><td>${w.dose > 0 ? `${w.drug} ${fmtConc(w.dose)}` : "none (control)"}</td><td>${w.live}</td><td>${w.viability.toFixed(0)} %</td></tr>`).join("");
    return `<h3>Dose–response · ${plate.name} at ${r.t} h</h3>${doseChart(r)}
      <p>${r.ic50 ? `IC50 ≈ <b>${fmtConc(r.ic50)}</b> of ${r.drug}: the dose that leaves half as many live cells as the control.` : r.note}</p>
      <div class="table-wrap"><table class="jt"><thead><tr><th>Well</th><th>Drug</th><th>Live cells</th><th>Viability</th></tr></thead><tbody>${rows}</tbody></table></div>`;
  }
}

function chart(dishes) {
  const W = 640, H = 240, L = 46, B = 30, T = 12, Rr = 14;
  const tmax = Math.max(24, ...dishes.map((d) => d.state.t));
  const nmax = Math.max(10, ...dishes.flatMap((d) => d.state.series.map((s) => s.n_live)));
  const x = (t) => L + (t / tmax) * (W - L - Rr), y = (n) => H - B - (n / nmax) * (H - B - T);
  const colors = ["#5b8cff", "#f0b44c", "#58c46e", "#ff7466"];
  const ticksT = [0, tmax / 2, tmax].map((t) => `<text x="${x(t)}" y="${H - 10}" text-anchor="middle">${Math.round(t)} h</text>`).join("");
  const ticksN = [0, nmax / 2, nmax].map((n) => `<text x="${L - 6}" y="${y(n) + 4}" text-anchor="end">${Math.round(n)}</text><line x1="${L}" x2="${W - Rr}" y1="${y(n)}" y2="${y(n)}" stroke="rgba(255,255,255,0.08)"/>`).join("");
  const lines = dishes.map((d, i) => {
    const pts = d.state.series.map((s) => `${x(s.t_h).toFixed(1)},${y(s.n_live).toFixed(1)}`).join(" ");
    const last = d.state.series[d.state.series.length - 1];
    return `<polyline points="${pts}" fill="none" stroke="${colors[i % 4]}" stroke-width="2.5"/><circle cx="${x(last.t_h)}" cy="${y(last.n_live)}" r="4" fill="${colors[i % 4]}"/><text x="${x(last.t_h) - 6}" y="${y(last.n_live) - 8}" text-anchor="end" style="fill:${colors[i % 4]}">${d.name}</text>`;
  }).join("");
  return `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="Live cells over time">${ticksN}${ticksT}${lines}<text x="${L}" y="${T}" >live cells</text></svg>`;
}

// Viability against dose on a log scale, with the 50 % line and the IC50.
function doseChart(r) {
  const pts = r.wells.filter((w) => w.dose > 0).sort((a, b) => a.dose - b.dose);
  if (!pts.length) return "";
  const W = 640, H = 260, L = 50, B = 34, T = 14, Rr = 40;
  const lo = Math.floor(Math.log10(pts[0].dose)), hi = Math.ceil(Math.log10(pts[pts.length - 1].dose) + 0.01);
  const x = (d) => L + ((Math.log10(d) - lo) / Math.max(1, hi - lo)) * (W - L - Rr), y = (v) => H - B - (Math.min(v, 120) / 120) * (H - B - T);
  let g = "";
  for (let e = lo; e <= hi; e++) g += `<line x1="${x(10 ** e)}" x2="${x(10 ** e)}" y1="${T}" y2="${H - B}" stroke="rgba(255,255,255,0.08)"/><text x="${x(10 ** e)}" y="${H - 12}" text-anchor="middle">${fmtConc(10 ** e)}</text>`;
  for (const v of [0, 50, 100]) g += `<text x="${L - 6}" y="${y(v) + 4}" text-anchor="end">${v}%</text>`;
  g += `<line x1="${L}" x2="${W - Rr}" y1="${y(50)}" y2="${y(50)}" stroke="rgba(240,180,76,0.6)" stroke-dasharray="5 4"/>`;
  g += `<line x1="${L}" x2="${W - Rr}" y1="${y(100)}" y2="${y(100)}" stroke="rgba(255,255,255,0.15)"/>`;
  const line = pts.map((p) => `${x(p.dose).toFixed(1)},${y(p.viability).toFixed(1)}`).join(" ");
  g += `<polyline points="${line}" fill="none" stroke="#5b8cff" stroke-width="2.5"/>`;
  g += pts.map((p) => `<circle cx="${x(p.dose)}" cy="${y(p.viability)}" r="4.5" fill="#5b8cff"/>`).join("");
  if (r.ic50) g += `<line x1="${x(r.ic50)}" x2="${x(r.ic50)}" y1="${y(50)}" y2="${H - B}" stroke="#f0b44c" stroke-width="1.5"/><text x="${x(r.ic50) + 6}" y="${y(50) - 6}" style="fill:#f0b44c">IC50</text>`;
  g += `<text x="${L}" y="${T - 2}">live cells, % of control</text>`;
  return `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img" aria-label="Dose–response curve">${g}</svg>`;
}
