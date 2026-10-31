// The lab journal: tasks, what you observed on each slide, and the growth
// curves the Cell Sim engine produced for each culture dish.
import { SPECIMENS } from "./specimens.js";

const KEY = "cellsim-lab:journal:v1";
const FACTS = {
  onion: { kind: "plant", shape: "Long bricks in neat rows", wall: "Yes", chloro: "None", vacuole: "Yes, large" },
  elodea: { kind: "plant", shape: "Small bricks, full of green dots", wall: "Yes", chloro: "Yes, many", vacuole: "Yes, large" },
  cheek: { kind: "animal", shape: "Flat and irregular", wall: "No", chloro: "None", vacuole: "No" },
};

export class Journal {
  constructor() {
    try { this.data = JSON.parse(localStorage.getItem(KEY)) || {}; } catch (e) { this.data = {}; }
    this.data.obs = this.data.obs || {};
  }
  save() { try { localStorage.setItem(KEY, JSON.stringify(this.data)); } catch (e) { /* private window */ } }
  observe(spec, patch) {
    const d = (this.data.obs[spec] = this.data.obs[spec] || {});
    if (patch.mag) d.mag = Math.max(d.mag || 0, patch.mag);
    for (const k of ["stain", "cover", "size"]) if (patch[k] != null) d[k] = patch[k];
    this.save();
  }
  render(tab, tasks, dishes) {
    if (tab === "tasks") {
      const v = tasks.view();
      if (!v) return `<p class="note">Free roam: no task list. Pick an experiment from the title screen to get one.</p>`;
      return `<h3>${v.title}</h3><ol class="checklist">${v.steps.map((s) => `<li class="${s.done ? "done" : ""}">${s.text}</li>`).join("")}</ol>${v.finished ? `<p class="note">Experiment complete.</p>` : ""}`;
    }
    if (tab === "obs") {
      const order = ["onion", "elodea", "cheek"], o = this.data.obs;
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
        <p class="note">Rows fill in once you've seen a sample in focus at 100× or more. These slides are drawn from measured cell sizes; they are not simulated.</p>`;
    }
    // Culture: live cell count over time for each dish, straight from the engine.
    const withData = dishes.filter((d) => d.state.series.length);
    if (!withData.length) return `<p class="note">No culture data yet. Look at a dish under the microscope, or incubate one, and the Cell Sim engine's counts appear here.</p>`;
    return chart(withData) + withData.map((d) => {
      const last = d.state.series[d.state.series.length - 1];
      const drugs = d.state.events.length ? d.state.events.map((e) => `${e.t} h: ${Object.entries(e.conc).map(([k, v]) => `${k} ${v.toFixed(2)} µM`).join(", ")}`).join(" · ") : "no drug";
      return `<p><b>${d.name}</b>: ${last.n_live} live cells at ${d.state.t} h (started with ${d.state.series[0].n_live}) · ${last.n_apoptotic} apoptotic · treatment: ${drugs}</p>`;
    }).join("") + `<p class="note">Every number here comes from cellsim.cell.dish, running in your browser.</p>`;
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
