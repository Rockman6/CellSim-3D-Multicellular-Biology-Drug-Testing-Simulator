/* Inside a cell — every molecule the engine tracks, animated from the
 * engine's own values.
 *
 * The page runs cellsim.cell.trace in the browser (Pyodide, same wheel as
 * the dose-response page), which records about twenty species for a few
 * cells every 15 minutes. Each circle's brightness is that species' level
 * at the current frame, scaled to the largest level it reaches in this
 * run so that a change is visible whatever its absolute size. Nothing is
 * drawn that the engine did not compute.
 */
const $ = (id) => document.getElementById(id);
const S = { py: null, trace: null, cell: 0, frame: 0, playing: false, timer: null, scale: {}, rest: {} };
const NS = "http://www.w3.org/2000/svg";

/* ── the map: where each species sits, which zone it belongs to ───── */
const ZONES = {
  drug:   { colour: "#a78bfa", label: "drug" },
  stress: { colour: "#e0a640", label: "damage response" },
  cycle:  { colour: "#4ea8de", label: "cell cycle" },
  death:  { colour: "#e0675a", label: "death" },
};
// `lab` places the label: below by default, or left/right where vertical
// arrows would otherwise run through it (the whole death column) and
// where a zone label sits underneath (the drug in the medium).
const NODES = {
  Cout: { x: 95,  y: 58,  z: "drug",   name: "Drug",         sub: "in the medium", lab: "right" },
  Cin:  { x: 95,  y: 195, z: "drug",   name: "Drug",         sub: "inside the cell" },
  D:    { x: 235, y: 195, z: "stress", name: "DNA damage",   sub: "" },
  ATM:  { x: 375, y: 195, z: "stress", name: "ATM",          sub: "damage sensor" },
  p53:  { x: 515, y: 195, z: "stress", name: "p53",          sub: "the alarm" },
  MDM2: { x: 600, y: 305, z: "stress", name: "MDM2",         sub: "turns p53 off" },
  p21:  { x: 375, y: 310, z: "stress", name: "p21",          sub: "brake on division", lab: "left" },
  CycD: { x: 95,  y: 470, z: "cycle",  name: "Cyclin D",     sub: "start" },
  E2F:  { x: 225, y: 470, z: "cycle",  name: "E2F",          sub: "opens the gate" },
  CycE: { x: 355, y: 470, z: "cycle",  name: "Cyclin E",     sub: "enter S" },
  CycA: { x: 485, y: 470, z: "cycle",  name: "Cyclin A",     sub: "copy DNA" },
  CycB: { x: 615, y: 470, z: "cycle",  name: "Cyclin B",     sub: "divide" },
  Puma: { x: 840, y: 80,  z: "death",  name: "PUMA",         sub: "death signal",        lab: "left" },
  Bax:  { x: 840, y: 165, z: "death",  name: "BAX",          sub: "pierces mitochondria", lab: "left" },
  MOMP: { x: 840, y: 250, z: "death",  name: "Mitochondria", sub: "opened",              lab: "left" },
  CytC: { x: 775, y: 345, z: "death",  name: "Cytochrome c", sub: "released" },
  Smac: { x: 905, y: 345, z: "death",  name: "Smac",         sub: "released" },
  C9:   { x: 840, y: 435, z: "death",  name: "Caspase-9",    sub: "",                    lab: "left" },
  C3:   { x: 840, y: 520, z: "death",  name: "Caspase-3",    sub: "the executioner",     lab: "left" },
};
// [from, to, kind] — "act" pushes the target up, "inh" holds it down
const EDGES = [
  ["Cout", "Cin", "act"], ["D", "ATM", "act"], ["ATM", "p53", "act"],
  ["p53", "MDM2", "act"], ["MDM2", "p53", "inh"], ["p53", "p21", "act"], ["p21", "CycE", "inh"],
  ["CycD", "E2F", "act"], ["E2F", "CycE", "act"], ["CycE", "CycA", "act"], ["CycA", "CycB", "act"],
  ["p53", "Puma", "act"], ["Puma", "Bax", "act"], ["Bax", "MOMP", "act"],
  ["MOMP", "CytC", "act"], ["MOMP", "Smac", "act"], ["CytC", "C9", "act"],
  ["C9", "C3", "act"], ["Smac", "C3", "act"],
];
const R = 24;

function el(tag, attrs = {}, parent) {
  const e = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

function drawEdge(layer, a, b, kind, bendBy = 0) {
  const A = NODES[a], B = NODES[b];
  const dx = B.x - A.x, dy = B.y - A.y, L = Math.hypot(dx, dy);
  const ux = dx / L, uy = dy / L;
  // two curved strokes between p53 and MDM2 so the loop reads as a loop
  const bend = bendBy || ((a === "p53" && b === "MDM2") || (a === "MDM2" && b === "p53") ? 16 : 0);
  const x1 = A.x + ux * R, y1 = A.y + uy * R, x2 = B.x - ux * (R + 3), y2 = B.y - uy * (R + 3);
  const mx = (x1 + x2) / 2 - uy * bend, my = (y1 + y2) / 2 + ux * bend;
  const p = el("path", { class: `edge ${kind === "inh" ? "inhib" : ""}`,
    d: `M${x1},${y1} Q${mx},${my} ${x2},${y2}`,
    "marker-end": `url(#${kind === "inh" ? "i" : "a"})` }, layer);
  p.dataset.from = a; p.dataset.kind = kind;
  return p;
}

// Where the drug's arrow goes, by mechanism. It used to point at DNA damage
// for every drug, which told the viewer that nutlin and paclitaxel damage
// DNA — neither does. A spindle poison acts in mitosis, so its arrow goes
// to cyclin B, the division node.
function setDrugEdge(mechanism, target, mode) {
  if (!S.drugLayer) return;
  S.drugLayer.innerHTML = "";
  const to = { dna_adduct: ["D", "act"], topo2: ["D", "act"], s_phase: ["D", "act"],
               mdm2: ["MDM2", "inh"], tubulin: ["CycB", "inh"] }[mechanism]
    || (mechanism === "target" && NODES[target] ? [target, mode === "boost" ? "act" : "inh"] : null);
  if (!to) return;
  if (to[0] === "D") { drawEdge(S.drugLayer, "Cin", "D", to[1]); return; }
  // Anything else is routed along the empty lanes between the zones, so the
  // arrow never runs through another molecule or its label.
  const C = NODES.Cin, T = NODES[to[0]], pts = [[C.x + R, C.y - 6]];
  if (NODES[to[0]].z === "death") {
    pts.push([140, 110], [672, 110], [672, T.y - 40], [T.x, T.y - 40], [T.x, T.y - R - 3]);
  } else if (NODES[to[0]].z === "cycle") {
    pts.push([140, C.y - 6], [140, 388], [T.x, 388], [T.x, T.y - R - 3]);
  } else if (T.y < 250) {                       // ATM, p53: from above
    pts.push([140, 150], [T.x, 150], [T.x, T.y - R - 3]);
  } else {                                      // p21, MDM2: the lane under the row
    pts.push([140, C.y - 6], [140, 268], [T.x, 268], [T.x, T.y - R - 3]);
  }
  const p = el("path", { class: `edge drug ${to[1] === "inh" ? "inhib" : ""}`,
    d: "M" + pts.map((q) => q.join(",")).join(" L"),
    "marker-end": `url(#${to[1] === "inh" ? "i" : "a"})` }, S.drugLayer);
  p.dataset.from = "Cin"; p.dataset.kind = to[1];
}

function drawMap() {
  const svg = $("map");
  svg.innerHTML = "";
  const defs = el("defs", {}, svg);
  for (const [id, colour] of [["a", "#3a4856"], ["aon", "#6fd3c4"]]) {
    const m = el("marker", { id, viewBox: "0 0 10 10", refX: 9, refY: 5,
      markerWidth: 6, markerHeight: 6, orient: "auto-start-reverse" }, defs);
    el("path", { d: "M0,0 L10,5 L0,10 z", fill: colour }, m);
  }
  for (const [id, colour] of [["i", "#3a4856"], ["ion", "#e0675a"]]) {
    const m = el("marker", { id, viewBox: "0 0 10 10", refX: 5, refY: 5,
      markerWidth: 7, markerHeight: 7, orient: "auto" }, defs);
    el("path", { d: "M5,0 L5,10", stroke: colour, "stroke-width": 3 }, m);
  }
  // zone frames, so the three programmes read as three programmes
  // the cell-cycle heading sits left: the drug's arrow enters that zone on the right
  const frames = [["damage response", 30, 125, 625, 245, "end"], ["cell cycle", 30, 405, 640, 125, "start"],
                  ["death", 690, 30, 270, 545, "end"]];
  for (const [lab, x, y, w, h, side] of frames) {
    el("rect", { class: "zone", x, y, width: w, height: h, rx: 10 }, svg);
    el("text", { class: "zone-label", x: side === "end" ? x + w - 12 : x + 12, y: y + 18,
                 "text-anchor": side }, svg).textContent = lab;
  }
  const edgeLayer = el("g", {}, svg);
  for (const [a, b, kind] of EDGES) drawEdge(edgeLayer, a, b, kind);
  S.drugLayer = el("g", {}, svg);
  setDrugEdge("dna_adduct");
  for (const [key, n] of Object.entries(NODES)) {
    const g = el("g", { class: "node", transform: `translate(${n.x},${n.y})` }, svg);
    g.dataset.key = key;
    el("circle", { class: "halo", r: R + 10, fill: ZONES[n.z].colour, opacity: 0 }, g);
    el("circle", { class: "core", r: R, fill: "#1d262e", stroke: ZONES[n.z].colour }, g);
    const side = n.lab || "below";
    const [lx, ly, anchor] = side === "left" ? [-R - 9, -1, "end"]
                           : side === "right" ? [R + 9, -1, "start"] : [0, R + 15, "middle"];
    el("text", { x: lx, y: ly, "text-anchor": anchor }, g).textContent = n.name;
    if (n.sub) el("text", { class: "sub", x: lx, y: ly + 13, "text-anchor": anchor }, g).textContent = n.sub;
    el("text", { class: "val", y: 4, "text-anchor": "middle" }, g);
  }
}

/* ── colour by level ───────────────────────────────────────────── */
function mix(hex, t) {
  const c = parseInt(hex.slice(1), 16), base = [0x1d, 0x26, 0x2e];
  const rgb = [(c >> 16) & 255, (c >> 8) & 255, c & 255];
  return `rgb(${rgb.map((v, i) => Math.round(base[i] + (v - base[i]) * t)).join(",")})`;
}
// How bright each species is drawn. FIXED, so that "lit" means the same
// thing in every scenario: the first version scaled each species to its own
// maximum in the run, which made a healthy cell's flat, resting p53 (0.096)
// glow exactly as brightly as nutlin's (2.0) — and told the viewer "p53 has
// risen" about a cell with no damage at all.
//
// Scales come from the engine's measured ranges across untreated,
// cisplatin, paclitaxel, nutlin and doxorubicin runs. Signalling proteins
// that sit at a non-zero resting level are drawn as their RISE above that
// level (taken from the run's first frame, before any drug has acted), so a
// resting cell is dark and a responding one lights up.
const DISPLAY = {
  // signalling: [scale, measured from rest]
  ATM: [1.0, true], p53: [0.30, true], MDM2: [0.60, true], MDM2_mRNA: [0.60, true],
  p21: [0.50, true],
  // the cycle oscillates in every healthy cell; show its natural swing
  CycD: [1.5], Rb: [0.9], E2F: [1.0], CycE: [0.75], CycA: [0.6], CycB: [0.9],
  // damage and the death cascade start at zero
  D: [0.30], Puma: [0.50], Bax: [0.30], MOMP: [0.15], CytC: [0.50], Smac: [0.50],
  C9: [0.40], C3: [0.40], Exec: [0.70],
};
function level(key, f) {
  const cell = S.trace.cells[S.cell];
  if (key === "Cout") {
    const t = S.trace.t_h[f], ex = S.trace.exposure_h;
    return S.trace.conc_uM > 0 && (ex == null || t < ex) ? 1 : 0;
  }
  const v = cell.series[key][f];
  // the drug's own level depends on the dose, so it alone is scaled per run
  if (key === "Cin") return Math.max(0, Math.min(1, v / (S.scale.Cin || 1)));
  const [scale, fromRest] = DISPLAY[key] || [1.0];
  const rest = fromRest ? (S.rest[key] || 0) : 0;
  return Math.max(0, Math.min(1, (v - rest) / Math.max(scale - rest, 1e-6)));
}

function render() {
  if (!S.trace) return;
  const f = S.frame, cell = S.trace.cells[S.cell], t = S.trace.t_h[f];
  const alive = cell.alive[f];
  for (const g of $("map").querySelectorAll(".node")) {
    const key = g.dataset.key, lv = level(key, f), z = ZONES[NODES[key].z].colour;
    g.querySelector(".core").setAttribute("fill", mix(z, alive || key === "Cout" ? lv : lv * 0.35));
    g.querySelector(".halo").setAttribute("opacity", (alive ? lv : 0) * 0.32);
    const raw = key === "Cout" ? S.trace.conc_uM * lv : cell.series[key][f];
    g.querySelector(".val").textContent =
      key === "Cout" ? (lv ? `${S.trace.conc_uM}` : "0")
      // two significant figures, so nanomolar paclitaxel (~0.02) is not
      // rounded to "0.0" and read as no drug at all
      : key === "Cin" ? (raw >= 10 ? raw.toFixed(0) : raw === 0 ? "0" : raw.toPrecision(2))
      : raw.toFixed(2);
  }
  for (const p of $("map").querySelectorAll(".edge")) {
    const on = alive && level(p.dataset.from, f) > 0.45;
    p.classList.toggle("on", on);
    const inh = p.dataset.kind === "inh";
    p.setAttribute("marker-end", `url(#${inh ? (on ? "ion" : "i") : (on ? "aon" : "a")})`);
    if (inh && on) p.style.stroke = "#e0675a"; else p.style.stroke = "";
  }
  $("clock").textContent = `${t.toFixed(1)} h`;
  $("scrub").value = f;
  drawCell(f);
  explain(f);
  drawLog(t);
  for (const c of document.querySelectorAll(".track .cursor")) {
    const x = (f / (S.trace.t_h.length - 1)) * 100;
    c.setAttribute("x1", `${x}%`); c.setAttribute("x2", `${x}%`);
  }
}

/* ── the cell itself ───────────────────────────────────────────── */
const PHASE = { G1: ["#4ea8de", "G1 — growing, getting ready"],
                S:  ["#35b3a2", "S — copying its DNA"],
                G2: ["#7c9cf0", "G2 — checking the copy"],
                M:  ["#f0c050", "M — dividing"],
                dead: ["#596572", "dead"] };
function drawCell(f) {
  const cell = S.trace.cells[S.cell], ph = cell.phase[f], svg = $("cellview");
  svg.innerHTML = "";
  const [col, words] = PHASE[ph] || PHASE.G1;
  const arrested = cell.arrested[f];
  const justDivided = cell.events.some(e => e.kind === "divided" &&
    Math.abs(e.t_h - S.trace.t_h[f]) < 0.6);
  if (ph === "dead") {
    // apoptotic bodies: the cell breaks into small blebs
    const pts = [[90, 100, 22], [125, 92, 16], [108, 128, 18], [140, 125, 12], [80, 132, 11]];
    for (const [x, y, r] of pts) el("circle", { cx: x, cy: y, r, fill: col, opacity: 0.55 }, svg);
  } else if (ph === "M" || justDivided) {
    const gap = justDivided ? 34 : arrested ? 10 : 18;
    for (const dx of [-gap, gap])
      el("circle", { cx: 110 + dx, cy: 110, r: 50, fill: col, opacity: 0.22,
                     stroke: col, "stroke-width": 2 }, svg);
    if (!justDivided) {   // chromosomes on the spindle
      el("line", { x1: 110, y1: 70, x2: 110, y2: 150, stroke: "#fff", "stroke-opacity": 0.5,
                   "stroke-width": arrested ? 4 : 2, "stroke-dasharray": "3 4" }, svg);
    }
  } else {
    el("circle", { cx: 110, cy: 110, r: 68, fill: col, opacity: 0.18, stroke: col, "stroke-width": 2 }, svg);
    const nr = ph === "G2" ? 32 : ph === "S" ? 28 : 24;    // nucleus grows as DNA is copied
    el("circle", { cx: 110, cy: 110, r: nr, fill: col, opacity: 0.45 }, svg);
    // damage shown as specks in the nucleus
    const dmg = level("D", f);
    for (let i = 0; i < Math.round(dmg * 9); i++) {
      const a = i * 2.4, rr = nr * 0.6;
      el("circle", { cx: 110 + Math.cos(a) * rr * (0.4 + (i % 3) * 0.25),
                     cy: 110 + Math.sin(a) * rr * (0.4 + (i % 3) * 0.25),
                     r: 2.6, fill: "#e0a640" }, svg);
    }
  }
  if (arrested && ph !== "dead")
    el("circle", { cx: 110, cy: 110, r: 92, fill: "none", stroke: "#e0675a",
                   "stroke-width": 2, "stroke-dasharray": "6 5" }, svg);
  $("phase").textContent = arrested && ph !== "dead" ? "stuck in mitosis" : words;
  $("phase").style.color = arrested && ph !== "dead" ? "#e0675a" : col;
}

/* ── plain words for what is happening ─────────────────────────── */
// Known to be wrong, measured, and shown to the viewer rather than hidden.
// Watching single cells is how these were found: every earlier check looked
// only at 72 h, where the endpoints are calibrated but the speeds never were.
// What the measurements say about a drug, shown beside the explanation the
// viewer is already reading. "known" is a fault the engine still has;
// "measured" is what real cells do, for comparison with what is on screen.
// Both nutlin and cisplatin carried "known" faults until October 2026 —
// nutlin's was a placeholder constant 77x too potent plus a p53 that killed
// without damage (docs/VALIDATION.md); cisplatin's turned out to be inside
// the measured range.
const NOTES = {
  "nutlin-3a": ["measured", "In A549 and HCT116, nutlin mostly stops division rather " +
    "than killing: Tovar et al. (2006) found 7.6 % and 8.8 % of cells dying after 72 hours. " +
    "The engine is fitted to that. Lines with extra copies of the MDM2 gene, such as SJSA-1 " +
    "(88 % dying), respond very differently and are not in this library."],
  "cisplatin": ["measured", "p53 rises about two-fold here. Across twelve cell lines with " +
    "normal p53, DNA breaks raise it between 1.25- and 5-fold (Stewart-Ornstein & Lahav " +
    "2017), so this is inside the measured range, at its low end. In single HCT116 cells " +
    "at a similar dose about half die, spread over three days (Paek et al. 2016)."],
};

function explain(f) {
  const cell = S.trace.cells[S.cell], lv = (k) => level(k, f), ph = cell.phase[f];
  let msg;
  if (ph === "dead") msg = "This cell has died. Caspase-3 dismantled it from the inside — " +
    "apoptosis, the orderly kind of cell death — and it broke into small fragments.";
  else if (lv("C3") > 0.5) msg = "Caspase-3 is active. Once it switches on there is no way back: " +
    "it is cutting up the cell's own proteins.";
  else if (lv("MOMP") > 0.5) msg = "The mitochondria have been pierced and are leaking cytochrome c " +
    "and Smac. Those switch on the caspases.";
  else if (cell.arrested[f]) msg = "Stuck in mitosis. The drug has frozen the spindle, so the " +
    "chromosomes cannot be pulled apart. The longer it waits, the more likely it is to die — " +
    "or to give up and slip back out without dividing.";
  else if (lv("Puma") > 0.5) msg = "p53 has switched on PUMA, the signal that starts the death " +
    "programme. Whether the cell dies now depends on how much protection it has.";
  else if (lv("p53") > 0.6 && lv("D") < 0.1) msg = "p53 has risen without any DNA damage, " +
    "because the drug is blocking MDM2, the protein that normally destroys it. p53 raised " +
    "this way stops division through p21, but it rarely starts the death programme: it lacks " +
    "the chemical marks that damage puts on it, and those are what point it at PUMA.";
  else if (lv("p53") > 0.6) msg = "p53 has risen — the cell's alarm. It turns on p21 to stop " +
    "division and, if the damage is bad enough, PUMA to start the death programme.";
  else if (lv("D") > 0.4) msg = "The drug is damaging DNA. ATM senses the breaks and starts " +
    "raising p53.";
  else if (lv("Cin") > 0.4) msg = "The drug is inside the cell.";
  else msg = { G1: "Growing and preparing. Cyclin D is pushing the cell toward the gate " +
                   "(Rb/E2F) where it commits to dividing.",
               S: "Past the gate and copying its DNA. Cyclin A is driving replication.",
               G2: "DNA copied. Cyclin B is building up for the division itself.",
               M: "Dividing — the chromosomes are being pulled into two new cells." }[ph] || "";
  const note = NOTES[S.trace.drug] || S.trace.customNote;
  $("explain").innerHTML = "";
  $("explain").append(msg);
  if (note) {
    const [kind, text] = note;
    const w = document.createElement("div");
    w.className = kind === "measured" ? "measured" : "known";
    const b = document.createElement("b");
    b.textContent = { known: "Known limitation. ", measured: "What real cells do. ",
                      design: "Your hypothesis. ", uncal: "Not calibrated. " }[kind];
    w.append(b, text);
    $("explain").append(w);
  }
}

function drawLog(t) {
  const evs = S.trace.cells[S.cell].events;
  $("log").innerHTML = evs.length
    ? evs.map(e => `<li class="${e.t_h <= t + 1e-6 ? "past" : "future"}"><span class="t">${e.t_h.toFixed(1)} h</span>${e.text}</li>`).join("")
    : `<li class="hint">No division, arrest or death in ${S.trace.t_h.at(-1).toFixed(0)} h.</li>`;
}

/* ── small multiples ───────────────────────────────────────────── */
const TRACKS = [["Cin", "drug inside", "#a78bfa"], ["D", "DNA damage", "#e0a640"],
                ["p53", "p53", "#e0a640"], ["C3", "caspase-3", "#e0675a"]];
function drawTracks() {
  const cell = S.trace.cells[S.cell], n = S.trace.t_h.length;
  $("tracks").innerHTML = TRACKS.map(([k, lab]) =>
    `<div class="track"><div class="tl">${lab}</div><svg viewBox="0 0 100 34" preserveAspectRatio="none" data-k="${k}"></svg></div>`).join("");
  for (const [k, , colour] of TRACKS) {
    const svg = $("tracks").querySelector(`svg[data-k="${k}"]`);
    // the same rule as the map's brightness, so the chart and the glow
    // never disagree about how large a change was
    const pts = cell.series[k].map((_, i) =>
      `${(i / (n - 1)) * 100},${32 - level(k, i) * 30}`).join(" ");
    el("polyline", { points: pts, fill: "none", stroke: colour, "stroke-width": 1.4,
                     "vector-effect": "non-scaling-stroke" }, svg);
    el("line", { class: "cursor", x1: "0%", x2: "0%", y1: 0, y2: 34,
                 "vector-effect": "non-scaling-stroke" }, svg);
  }
}

/* ── cell picker ───────────────────────────────────────────────── */
function fateOf(cell) {
  const d = cell.events.find(e => e.kind === "died");
  if (d) return `died ${d.t_h.toFixed(0)} h`;
  const n = cell.events.filter(e => e.kind === "divided").length;
  if (cell.arrested.some(Boolean)) return "arrested";
  return n ? `${n} division${n > 1 ? "s" : ""}` : "no change";
}
function drawPicker() {
  $("cellpick").innerHTML = S.trace.cells.map((c, i) =>
    `<button data-i="${i}" class="${i === S.cell ? "on" : ""}">Cell ${i + 1}<span class="fate">${fateOf(c)}</span></button>`).join("");
  for (const b of $("cellpick").querySelectorAll("button")) b.onclick = () => {
    S.cell = Number(b.dataset.i); drawPicker(); drawTracks(); render();
  };
}

/* ── playback ──────────────────────────────────────────────────── */
function play(on) {
  S.playing = on;
  $("play").textContent = on ? "❚❚" : "▶";
  clearInterval(S.timer);
  if (!on) return;
  S.timer = setInterval(() => {
    const step = Number($("speed").value);
    if (S.frame >= S.trace.t_h.length - 1) { play(false); return; }
    S.frame = Math.min(S.frame + step, S.trace.t_h.length - 1);
    render();
  }, 80);
}

/* ── any drug: look-up ─────────────────────────────────────────── */
const esc = (s) => String(s).replace(/[&<>"']/g, (c) =>
  ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const ALIASES = { "fluorouracil": "5-fluorouracil", "5-fu": "5-fluorouracil",
                  "sn38": "sn-38", "nutlin-3": "nutlin-3a", "nutlin": "nutlin-3a" };

// A drug the library already has is used from the library: that one is
// calibrated, a looked-up copy of it would not be.
function libraryName(...names) {
  const lib = [...$("drug").options].map((o) => o.value).filter((v) => v && v !== "__custom__");
  for (const n of names) {
    if (!n) continue;
    const k = n.trim().toLowerCase();
    const hit = lib.find((x) => x.toLowerCase() === k) || ALIASES[k];
    if (hit && lib.includes(hit)) return hit;
  }
  return null;
}

async function doLookup() {
  const out = $("lookup-out");
  out.hidden = false;
  out.innerHTML = `<span class="hint">Looking it up…</span>`;
  $("lookup").disabled = true;
  try {
    const r = await Lookup.lookup($("q").value);
    if (r.kind === "empty") { out.hidden = true; return; }
    if (r.kind === "invalid") {
      out.innerHTML = `<div class="verdict bad"><b>This cannot exist as written.</b></div>
        <div>${esc(r.error)}</div>`;
      return;
    }
    const lib = r.kind === "known" ? libraryName(r.query, r.name) : null;
    const texts = (r.chembl && r.chembl.texts) || [];
    const p = r.props;
    const py = JSON.parse(await S.py.runPythonAsync(`
import json
from cellsim.cell.compound import classify, transport_from_properties
a = json.loads(${JSON.stringify(JSON.stringify({ texts, p }))})
c = classify(a["texts"])
t = transport_from_properties(mw=a["p"]["mw"], tpsa=a["p"]["tpsa"], logp=a["p"]["logp"],
                              n_plus_o=a["p"]["n_plus_o"])
json.dumps({"c": c and c.__dict__, "notes": t.notes, "pgp": t.pgp_substrate})`));
    const num = (v, d = 1) => (v == null || isNaN(v) ? "–" : Number(v).toFixed(d));
    let html = `<span class="mol">${r.svg}</span>
      <h5>${esc(r.name)}</h5>
      <div class="props">${p.formula ? esc(p.formula) + " · " : ""}MW ${num(p.mw, 0)} · logP ${num(p.logp)}
        · polar surface ${num(p.tpsa, 0)} Å²</div>
      <div class="from">${esc(r.propsFrom)}${r.chembl ? ` · ChEMBL ${esc(r.chembl.id)}` : ""}</div>`;
    S.found = null;
    if (lib) {
      html += `<div class="verdict">In the library, with its strength fitted to measured
        GDSC data.</div><button class="use primary" data-lib="${esc(lib)}">Use ${esc(lib)}</button>`;
    } else if (py.c) {
      S.found = { name: r.name, texts, p, ref: py.c.reference };
      html += `<div class="verdict"><b>What it does:</b> ${esc(py.c.words)}
          <div class="from">ChEMBL: “${esc(py.c.matched)}”</div></div>
        <div class="strength">Strength
          <select id="strength">${[0.01, 0.1, 1, 10, 100].map((s) =>
            `<option value="${s}"${s === 1 ? " selected" : ""}>×${s} of ${esc(py.c.reference)}</option>`).join("")}
          </select></div>
        <div class="from">Not calibrated: borrowed from ${esc(py.c.reference)}, which works the
          same way. Drugs of one class can differ a hundred-fold in strength.</div>
        <button class="use primary" data-custom="1">Use ${esc(r.name)}</button>`;
    } else if (r.chemblError) {
      html += `<div class="verdict no"><b>ChEMBL did not answer</b> (${esc(r.chemblError)}),
        so what this drug targets is unknown for now. Try the look-up again in a moment.</div>`;
    } else if (texts.length) {
      html += `<div class="verdict no"><b>Not something the engine models yet.</b> ChEMBL says:
        ${texts.map((t) => "“" + esc(t) + "”").join(", ")}. Running it as anything else
        would show a mechanism this drug does not have.</div>`;
    } else if (r.kind === "known") {
      html += `<div class="verdict no"><b>No recorded target.</b> PubChem knows this compound,
        but no curated record says what it binds, so there is nothing to wire it to.</div>`;
    } else {
      html += `<div class="verdict no"><b>A molecule nobody has recorded.</b> Its properties
        above are computed from your structure. What it would bind cannot be computed
        reliably from structure, so it cannot be run as a drug.</div>`;
    }
    if (py.notes.length) html += `<ul>${py.notes.map((n) => `<li>${esc(n)}</li>`).join("")}</ul>`;
    if (!lib && !py.c) html += `<button class="to-design" data-design="1">Design what it does →</button>`;
    out.innerHTML = html;
    const des = out.querySelector("button[data-design]");
    if (des) des.onclick = () => designFromLookup(r.kind === "new" ? "your molecule" : r.name, py.pgp);
    const use = out.querySelector("button.use");
    if (use) use.onclick = () => useFound(use.dataset.lib || null);
  } catch (e) {
    out.innerHTML = `<div class="verdict bad">The look-up did not complete:
      ${esc(String(e).slice(0, 200))}</div>`;
  } finally {
    $("lookup").disabled = false;
  }
}

function useFound(lib) {
  const sel = $("drug");
  if (lib) { sel.value = lib; return; }
  S.custom = { ...S.found, strength: Number($("strength").value) };
  let opt = [...sel.options].find((o) => o.value === "__custom__");
  if (!opt) { opt = new Option("", "__custom__"); sel.add(opt); }
  opt.text = `${S.custom.name} (looked up, uncalibrated)`;
  sel.value = "__custom__";
}

/* ── design mode ───────────────────────────────────────────────── */
// The molecules a design may block or boost (cellsim.cell.compound.DESIGNABLE):
// everything on the map except the drug itself and the damage index.
const DESIGNABLE = ["CycD", "Rb", "E2F", "CycE", "CycA", "CycB", "p21", "ATM", "p53",
  "MDM2_mRNA", "MDM2", "Puma", "Bax", "MOMP", "CytC", "Smac", "C9", "C3", "Exec"];

function designUI() {
  const a = $("d-action").value, onMolecule = a === "block" || a === "boost";
  $("d-target").hidden = !onMolecule;
  $("d-strength-row").hidden = onMolecule;
  const was = $("d-effect").dataset.kind;
  if (onMolecule && was !== a) {
    $("d-effect").innerHTML = (a === "block"
      ? [[1, "all of it"], [0.8, "most of it — 80 %"], [0.5, "half of it"]]
      : [[2, "doubles it"], [5, "five-fold"], [10, "ten-fold"]])
      .map(([v, t]) => `<option value="${v}">${t}</option>`).join("");
    $("d-effect").dataset.kind = a;
  }
}

function describeDesign(g) {
  const sp = $("d-species").selectedOptions[0];
  if (g.action === "block" || g.action === "boost") {
    const what = sp ? sp.textContent : g.species;
    return g.action === "block"
      ? `blocks ${what} — ${Math.round(g.effect * 100)} % of it when fully bound, binding with Kd ${g.kd} µM`
      : `boosts ${what} ${g.effect}-fold when fully bound, binding with Kd ${g.kd} µM`;
  }
  return $("d-action").selectedOptions[0].textContent + ` at ×${g.strength} its strength`;
}

function useDesign() {
  const g = { action: $("d-action").value, species: $("d-species").value,
              kd: Number($("d-kd").value), effect: Number($("d-effect").value),
              strength: Number($("d-strength").value) };
  const name = $("d-name").value.trim() || "my drug";
  S.custom = { design: g, name, pgp: S.designPgp ?? null, words: describeDesign(g) };
  const sel = $("drug");
  let opt = [...sel.options].find((o) => o.value === "__custom__");
  if (!opt) { opt = new Option("", "__custom__"); sel.add(opt); }
  opt.text = `${name} (your design)`;
  sel.value = "__custom__";
  $("lookup-out").hidden = true;
}

// From the look-up: a molecule with no usable target can still be given one.
function designFromLookup(name, pgp) {
  $("design-step").open = true;
  $("d-name").value = name;
  S.designPgp = pgp;
  $("d-from").hidden = false;
  $("d-from").textContent = `Designing what ${name} does. Its structure ` +
    (pgp === true ? "suggests P-glycoprotein will pump it out of cells, and the engine will."
     : pgp === false ? "suggests P-glycoprotein will not pump it out."
     : "does not say whether P-glycoprotein pumps it out; the engine assumes not.");
  $("design-step").scrollIntoView({ behavior: "smooth", block: "start" });
}

/* ── engine ────────────────────────────────────────────────────── */
async function boot() {
  $("engine-state").textContent = "downloading Python…";
  S.py = await loadPyodide({ indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.2/full/" });
  $("engine-state").textContent = "loading numpy and scipy…";
  await S.py.loadPackage(["numpy", "scipy", "micropip"]);
  $("engine-state").textContent = "installing cellsim…";
  // Same wheel as the dose-response page: built from the deployed commit.
  const wheel = (await (await fetch("wheel/LATEST")).text()).trim();
  if (!wheel.endsWith(".whl")) throw new Error(`no wheel found (LATEST = "${wheel}")`);
  const url = new URL(`wheel/${wheel}`, location.href).href;
  await S.py.runPythonAsync(`import micropip\nawait micropip.install(${JSON.stringify(url)})`);
  const meta = JSON.parse(await S.py.runPythonAsync(`
import json
from cellsim.cell.library import CELL_LINES, DRUGS
from cellsim.cell.trace import SPECIES
json.dumps({"lines": list(CELL_LINES), "drugs": list(DRUGS), "species": SPECIES})`));
  for (const k of DESIGNABLE) $("d-species").add(new Option(meta.species[k] || k, k));
  $("d-species").value = "Bax";
  for (const n of meta.lines) $("line").add(new Option(n, n));
  for (const n of meta.drugs) $("drug").add(new Option(n, n));
  $("line").value = "A549"; $("drug").value = "cisplatin";
  $("engine-state").textContent = "engine ready";
  $("run").disabled = false;
  $("lookup").disabled = false;
  $("d-use").disabled = false;
}

async function simulate() {
  $("run").disabled = true; play(false);
  $("engine-state").textContent = "simulating…";
  $("bar").classList.add("indet");
  const custom = $("drug").value === "__custom__" ? S.custom : null;
  const cfg = { line: $("line").value, drug: $("drug").value || null,
                conc: Number($("conc").value) || 0, custom,
                exposure: $("exposure").value ? Number($("exposure").value) : null };
  try {
    const out = JSON.parse(await S.py.runPythonAsync(`
import json
from cellsim.cell.trace import trace_cells
from cellsim.cell.compound import classify, from_mechanism, transport_from_properties
cfg = json.loads(${JSON.stringify(JSON.stringify(cfg))})
drug = cfg["drug"]
if cfg["custom"] and cfg["custom"].get("design"):
    from cellsim.cell.compound import Transport, designed
    c = cfg["custom"]; g = c["design"]
    t = None if c.get("pgp") is None else Transport(pgp_substrate=c["pgp"], permeable=None)
    drug = designed(c["name"], g["action"], species=g["species"], kd_uM=g["kd"],
                    effect=g["effect"], strength=g["strength"], transport=t)
elif cfg["custom"]:
    c = cfg["custom"]; p = c["p"]
    drug = from_mechanism(c["name"], classify(c["texts"]), strength=c["strength"],
                          transport=transport_from_properties(mw=p["mw"], tpsa=p["tpsa"],
                                                              logp=p["logp"], n_plus_o=p["n_plus_o"]))
tr = trace_cells(cfg["line"], drug, cfg["conc"] if drug else 0.0,
                 hours=72.0, n_cells=6, every_h=0.25, exposure_h=cfg["exposure"])
d = tr.to_dict(); d["exposure_h"] = cfg["exposure"]
if drug:
    from cellsim.cell.library import get_drug
    dg = get_drug(drug) if isinstance(drug, str) else drug
    d["mechanism"], d["target"], d["target_mode"] = dg.mechanism, dg.target_species, dg.target_mode
json.dumps(d)`));
    out.customNote = custom && custom.design ? ["design", `${custom.name} ${custom.words}. ` +
      `Nothing about it is measured: this is the engine playing out your idea, at the dose ` +
      `set in step 1, in a cell that is otherwise the Lab's calibrated one.` +
      (custom.design.action === "boost" ? " A booster multiplies what the cell is already " +
        "making — if the cell makes none of it at the moment (PUMA, say, while p53 is quiet), " +
        "there is nothing to multiply." : "")]
      : custom ? ["uncal", `${custom.name}'s strength is ` +
      `borrowed from ${custom.ref} (×${custom.strength}), which acts the same way; drugs of one ` +
      `class can differ a hundred-fold. Read this as what the mechanism does, not as what ` +
      `${custom.name} does at this dose.`] : null;
    S.trace = out;
    setDrugEdge(out.mechanism || null, out.target, out.target_mode);
    // scale each species to the largest level it reaches in this run, with
    // a floor, so a flat-zero species is not inflated into noise
    // the drug inside is dose-dependent, so it alone is scaled to this run
    let m = 0;
    for (const c of out.cells) for (const v of c.series.Cin) if (v > m) m = v;
    S.scale.Cin = Math.max(m, 1e-6);
    // resting levels: the median across cells of the first frame, before
    // any drug has had time to act
    S.rest = {};
    for (const k of Object.keys(out.species)) {
      const first = out.cells.map(c => c.series[k][0]).sort((a, b) => a - b);
      S.rest[k] = first[Math.floor(first.length / 2)];
    }
    S.cell = 0; S.frame = 0;
    $("placeholder").hidden = true;
    $("scrub").max = out.t_h.length - 1;
    for (const id of ["scrub", "play"]) $(id).disabled = false;
    drawPicker(); drawTracks(); render(); play(true);
  } catch (e) {
    $("explain").textContent = `The simulation did not complete: ${String(e).slice(0, 300)}`;
  } finally {
    $("run").disabled = false; $("bar").classList.remove("indet");
    $("engine-state").textContent = "engine ready";
  }
}

const PRESETS = {
  none: { drug: "", conc: 0, exposure: "" },
  cis:  { drug: "cisplatin", conc: 20, exposure: "" },
  pac:  { drug: "paclitaxel", conc: 0.1, exposure: "" },
  nut:  { drug: "nutlin-3a", conc: 10, exposure: "" },
};

document.addEventListener("DOMContentLoaded", () => {
  drawMap();
  $("run").onclick = simulate;
  $("lookup").onclick = doLookup;
  $("d-action").onchange = designUI;
  $("d-use").onclick = useDesign;
  designUI();
  $("q").onkeydown = (e) => { if (e.key === "Enter" && !$("lookup").disabled) doLookup(); };
  $("play").onclick = () => play(!S.playing);
  $("scrub").oninput = () => { play(false); S.frame = Number($("scrub").value); render(); };
  for (const b of document.querySelectorAll(".preset")) b.onclick = () => {
    const p = PRESETS[b.dataset.preset];
    $("line").value = "A549"; $("drug").value = p.drug; $("conc").value = p.conc;
    $("exposure").value = p.exposure;
    if (!$("run").disabled) simulate();
  };
  boot().catch(e => {
    $("engine-state").textContent = "engine failed to start";
    $("explain").textContent = String(e).slice(0, 300);
  });
});
