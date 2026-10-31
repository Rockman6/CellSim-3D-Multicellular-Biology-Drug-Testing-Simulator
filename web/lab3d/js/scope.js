// Looking through the eyepiece: objective lenses, coarse and fine focus,
// light and diaphragm, moving the slide, measuring, inspecting a living
// cell's engine state, naming the stage of dividing cells, and counting
// cells in a hemocytometer.
import { makeField, makeHemoField, SPECIMENS, FOV_UM, STAGES, CORNERS } from "./specimens.js";
import { makeCultureField } from "./culture.js";

const $ = (id) => document.getElementById(id);
const FOCUS_AT = { 4: 60, 10: 61.5, 40: 62.4 };  // nearly parfocal, like a real microscope
const DEPTH = { 4: 5, 10: 2.6, 40: 1.1 };        // depth of field shrinks at high power
const STEP = { coarse: 2.5, fine: 0.25 };
const STAGE_NOTES = [
  "Interphase: one round nucleus with a dark nucleolus; the DNA is spread out, not yet visible as chromosomes.",
  "Prophase: the chromosomes condense into tangled threads inside the nucleus.",
  "Metaphase: the chromosomes line up across the middle of the cell.",
  "Anaphase: the two sets of chromosomes are pulled apart, towards the ends of the cell.",
  "Telophase: two new nuclei form at the ends, and a cell plate grows between them.",
];

export function createScope({ onFocus, onExit, toast, observe, onEvent = () => {}, resolution = 720 }) {
  const canvas = $("field"), ep = $("eyepiece"), msg = $("epMsg");
  const W = resolution;
  canvas.width = W; canvas.height = W;
  const ctx = canvas.getContext("2d");
  const s = { spec: null, field: null, obj: 4, focus: 28, light: 0.7, diaphragm: 0.35, cx: 0, cy: 0, t: 0, dirty: true, done: {}, told: {}, sel: null, osmoT: 0 };

  const offBy = () => Math.abs(s.focus - FOCUS_AT[s.obj]) / DEPTH[s.obj];
  const inFocus = () => offBy() < 0.8;
  const fov = () => FOV_UM[s.obj];
  const onSample = () => s.field && s.field.inSample(s.cx, s.cy);

  function open(spec) {
    if (spec.kind === "wait") {
      s.spec = null; s.field = null;
      $("scopeTitle").textContent = spec.title;
      $("cultureTime").textContent = "";
      ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.fillStyle = "#cfc9c0"; ctx.fillRect(0, 0, W, W);
      msg.textContent = spec.text; msg.hidden = false;
      for (const id of ["cellCard", "measureBtn", "classifyPanel", "countPanel", "osmoReadout", "measurePanel"]) $(id).hidden = true;
      $("wellPick").hidden = !spec.wells;
      return;
    }
    msg.textContent = "Nothing here. Drag to move the slide.";
    const keep = s.spec && spec.keepView && s.spec.kind === spec.kind;
    s.spec = spec;
    s.field = spec.kind === "culture" ? makeCultureField(spec.header, spec.snap, { tint: spec.tint }) : spec.kind === "hemo" ? makeHemoField(spec.hemo) : makeField(spec.prep);
    if (!keep) Object.assign(s, { obj: 4, focus: 28, cx: 0, cy: 0, done: {}, told: {} });
    Object.assign(s, { t: 0, dirty: true, sel: null });
    $("scopeTitle").textContent = spec.title;
    $("cultureTime").textContent = spec.kind === "culture" ? ` · ${spec.time} h in culture` : "";
    $("measurePanel").hidden = true;
    $("measureLine").hidden = true;
    $("cellCard").hidden = true;
    $("mCount").value = ""; $("mResult").textContent = "";
    $("measureBtn").hidden = spec.kind !== "slide" || !spec.prep.specimen;
    $("classifyPanel").hidden = !s.field.classify;
    $("countPanel").hidden = !s.field.hemo;
    $("osmoReadout").hidden = !s.field.osmosis;
    $("wellPick").hidden = !spec.wells;
    if (spec.wells) renderWells();
    if (s.field.classify) { $("clsResult").textContent = ""; renderTally(); }
    // Which cells were counted or named is kept by index: the field is drawn fresh each time, from the same seed.
    if (s.field.hemo) { s.field.cells.forEach((c, i) => (c.counted = spec.hemo.counted.has(i))); $("countResult").textContent = ""; renderCount(); }
    if (s.field.classify) s.field.cells.forEach((c, i) => { if (spec.tally.named.has(i)) c.named = true; });
    if (s.field.classify && s.field.clumped) toast(s.field.soft ? "The cells are piled in layers. Squash the slide under a paper towel to spread them out." : "The cells are stuck together in clumps: the root tip wasn't softened in warm acid.", "warn");
    sync(); apply();
  }

  function apply() {
    const blur = Math.min(14, Math.pow(offBy(), 1.15) * 1.6);
    const bright = 0.45 + s.light * 0.75 - s.diaphragm * 0.3;
    const contrast = 0.8 + s.diaphragm * 0.8 - Math.max(0, s.light - 0.7) * 0.8;
    canvas.style.filter = `blur(${blur.toFixed(2)}px) brightness(${bright.toFixed(2)}) contrast(${contrast.toFixed(2)})`;
    $("magTotal").textContent = `${s.obj * 10}×`;
    const um = fov();
    $("fovText").textContent = um >= 1000 ? `${(um / 1000).toFixed(1)} mm` : `${um} µm`;
    $("mFov").textContent = `${um} µm`;
    check();
  }

  function sync() {
    for (const b of document.querySelectorAll("[data-obj]")) b.setAttribute("aria-pressed", String(+b.dataset.obj === s.obj));
    $("light").value = Math.round(s.light * 100);
    $("diaphragm").value = Math.round(s.diaphragm * 100);
  }

  function check() {
    if (!s.spec) return;
    const good = inFocus() && onSample();
    if (good && !s.done[s.obj]) {
      s.done[s.obj] = true;
      onFocus(s.obj, s.spec);
      if (s.spec.kind === "slide" && s.spec.prep.specimen) observe(s.spec.prep.specimen, { mag: s.obj * 10, stain: s.spec.stainName, cover: s.spec.prep.cover });
    }
    if (good && s.obj >= 10 && s.spec.kind === "slide") {
      const p = s.spec.prep;
      if (p.cover === "flat" && !s.told.bubble && bubbleInView()) { s.told.bubble = true; toast("The circles with thick dark rims are air bubbles, from dropping the cover slip flat.", "warn"); }
      if (p.specimen && p.stain !== SPECIMENS[p.specimen].best && !s.told.stain) {
        s.told.stain = true;
        toast({ cheek: "Hard to see? Cheek cells are almost clear. Methylene blue would stain them.", onion: "The nuclei are faint. Iodine stains onion cells best.",
          root: "The chromosomes are faint. Toluidine blue stains DNA dark purple.", redonion: "No stain needed: the red sap already shows each cell.",
          elodea: "A stain darkens the chloroplasts. Plain water shows the green best." }[p.specimen], "warn");
      }
    }
    msg.hidden = onSample() || !s.field;
    renderTasks();
  }

  function bubbleInView() {
    const h = fov() / 2;
    return (s.field.bubbles || []).some((b) => Math.abs(b.x - s.cx) < h + b.rad && Math.abs(b.y - s.cy) < h + b.rad);
  }

  function renderTasks() {
    const T = s.spec && s.spec.kind === "hemo"
      ? [[4, "Find the ruled grid with the 4× lens"], [10, "Switch to 10×: one corner square fills the view"]]
      : [[4, "Focus with the 4× lens (coarse knob)"], [10, "Switch to 10× and sharpen"], [40, "Switch to 40×, fine knob only"]];
    const extra = !s.spec ? "" : s.spec.kind === "culture" ? `<li>Click a cell to see what the engine says about it</li>`
      : s.field?.classify ? `<li>Click a cell, then pick its stage below</li>` : s.field?.hemo ? `<li>Click each cell in the four corner squares</li>` : "";
    $("scopeTasks").innerHTML = T.map(([k, t]) => `<li class="${s.done[k] ? "done" : ""}">${t}</li>`).join("") + extra;
  }

  function focusStep(kind, dir) {
    if (!s.spec) return;
    if (kind === "coarse" && s.obj === 40) { toast("At 40× use only the fine knob. The coarse knob moves the lens fast and could crack the slide.", "warn"); return false; }
    s.focus = Math.max(0, Math.min(100, s.focus + dir * STEP[kind]));
    apply();
    return true;
  }
  for (const b of document.querySelectorAll("[data-focus]")) {
    const [kind, d] = b.dataset.focus.split(":");
    let timer = null, rep = null;
    const stop = () => { clearTimeout(timer); clearInterval(rep); };
    b.addEventListener("pointerdown", (e) => { e.preventDefault(); if (focusStep(kind, +d) === false) return; timer = setTimeout(() => (rep = setInterval(() => focusStep(kind, +d), 55)), 320); });
    for (const ev of ["pointerup", "pointerleave", "pointercancel"]) b.addEventListener(ev, stop);
    b.addEventListener("click", (e) => { if (e.detail === 0) focusStep(kind, +d); });
  }
  for (const b of document.querySelectorAll("[data-obj]")) b.addEventListener("click", () => {
    const next = +b.dataset.obj;
    if (next === s.obj) return;
    if (next === 40 && s.obj === 4) toast("Tip: go up one lens at a time, 4× then 10× then 40×.");
    s.obj = next; s.dirty = true; sync(); apply();
  });
  $("light").addEventListener("input", (e) => { s.light = e.target.value / 100; apply(); });
  $("diaphragm").addEventListener("input", (e) => { s.diaphragm = e.target.value / 100; apply(); });
  ep.addEventListener("wheel", (e) => { e.preventDefault(); focusStep(e.shiftKey ? "coarse" : "fine", e.deltaY < 0 ? 1 : -1); }, { passive: false });

  // Drag to move the slide; a click without dragging picks a cell.
  let drag = null;
  ep.addEventListener("pointerdown", (e) => { drag = { x: e.clientX, y: e.clientY, moved: 0 }; ep.setPointerCapture(e.pointerId); });
  ep.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const k = fov() / ep.getBoundingClientRect().width;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag.moved += Math.abs(dx) + Math.abs(dy);
    s.cx = Math.max(-4200, Math.min(4200, s.cx - dx * k));
    s.cy = Math.max(-3200, Math.min(3200, s.cy - dy * k));
    drag.x = e.clientX; drag.y = e.clientY; s.dirty = true;
  });
  ep.addEventListener("pointerup", (e) => {
    if (drag && drag.moved < 5 && s.field && s.field.pickCell) pick(e);
    drag = null; check();
  });

  function pick(e) {
    const r = ep.getBoundingClientRect();
    const k = fov() / r.width;
    const x = s.cx + (e.clientX - r.left - r.width / 2) * k, y = s.cy + (e.clientY - r.top - r.height / 2) * k;
    const c = s.field.pickCell(x, y);
    if (s.field.classify) return selectCell(c);
    if (s.field.hemo) return countCell(c);
    inspect(c);
  }

  function inspect(c) {
    const card = $("cellCard");
    if (!c) { card.hidden = true; return; }
    const fate = c.c3 > 0.35 ? "dying (apoptosis switched on)" : c.phase === 3 ? "dividing" : c.gs < 0.05 ? "waiting: no room to divide" : "growing";
    card.innerHTML = `<b>One A549 cell · ${fate}</b>
      <dl><dt>Cell-cycle phase</dt><dd>${s.field.phaseName(c.phase)}</dd><dt>p53</dt><dd>${c.p53.toFixed(3)}</dd>
      <dt>Caspase-3</dt><dd>${c.c3.toFixed(3)}</dd><dt>DNA damage</dt><dd>${c.dmg.toFixed(3)}</dd>
      <dt>Oxygen</dt><dd>${c.o2.toFixed(0)} mmHg</dd><dt>Growth signal</dt><dd>${c.gs.toFixed(2)}</dd><dt>Lineage</dt><dd>#${c.lin}</dd></dl>
      <p class="note">Values from the CellSim engine for this cell at ${s.spec.time} h. <a href="../lab/cell.html" target="_blank" rel="noopener">Watch inside one cell →</a></p>`;
    card.hidden = false;
  }

  /* ---------- naming the stage of dividing cells ---------- */
  $("clsButtons").innerHTML = STAGES.map((n, i) => `<button class="mbtn small" data-stage="${i}">${n}</button>`).join("");
  function selectCell(c) {
    if (!c) { s.sel = null; s.dirty = true; return; }
    if (s.obj < 40) toast("Use the 40× lens to see the chromosomes well enough.", "warn");
    s.sel = c; s.dirty = true;
    $("clsResult").textContent = c.named != null ? `Already named: ${STAGES[c.stage]}.` : "Which stage is this cell in?";
  }
  $("clsButtons").addEventListener("click", (e) => {
    const b = e.target.closest("[data-stage]");
    if (!b || !s.field?.classify) return;
    const c = s.sel;
    if (!c) return toast("Click a cell in the eyepiece first.");
    if (c.named != null) return toast("You've already named that one. Pick another cell.");
    const guess = +b.dataset.stage, right = guess === c.stage;
    c.named = guess;
    const t = s.spec.tally;
    t.named.add(s.field.cells.indexOf(c));
    t.counts[c.stage]++; t.total++; if (right) t.right++;
    $("clsResult").innerHTML = `${right ? "<b>Right.</b> " : `<b>Not quite:</b> that one is in ${STAGES[c.stage].toLowerCase()}. `}${STAGE_NOTES[c.stage]}`;
    renderTally();
    onEvent("classify", t);
    s.sel = null; s.dirty = true;
  });
  function renderTally() {
    const t = s.spec.tally, dividing = t.total - t.counts[0];
    $("clsTally").innerHTML = t.total
      ? `${STAGES.map((n, i) => `${n} ${t.counts[i]}`).join(" · ")}<br>Mitotic index: ${dividing} ÷ ${t.total} = <b>${(dividing / t.total * 100).toFixed(1)} %</b>`
      : "None named yet.";
  }

  /* ---------- counting cells in the hemocytometer ---------- */
  function countCell(c) {
    if (!c) return;
    if (c.sq < 0) return toast("That cell is outside the four corner squares. Count only those.", "warn");
    c.counted = !c.counted;
    const i = s.field.cells.indexOf(c);
    if (c.counted) s.spec.hemo.counted.add(i); else s.spec.hemo.counted.delete(i);
    s.dirty = true;
    renderCount();
    const got = CORNERS.map((_, i) => s.field.cells.filter((x) => x.sq === i && x.counted).length);
    const truth = s.field.truth.map((q) => q.live + q.dead);
    if (got.every((n, i) => n >= truth[i] - 1 && truth[i] > 0)) onEvent("count:4");
  }
  function counts() {
    return CORNERS.map((_, i) => ({ live: s.field.cells.filter((x) => x.sq === i && x.counted && !(x.dead && s.spec.hemo.trypan)).length, dead: s.field.cells.filter((x) => x.sq === i && x.counted && x.dead && s.spec.hemo.trypan).length }));
  }
  function renderCount() {
    const c = counts();
    $("countTable").innerHTML = `<tr><th></th><th>Live</th><th>Dead (blue)</th></tr>` + CORNERS.map((q, i) => `<tr><th>${q.name}</th><td>${c[i].live}</td><td>${c[i].dead}</td></tr>`).join("");
  }
  $("countCalc").addEventListener("click", () => {
    if (!s.field?.hemo) return;
    const c = counts(), h = s.spec.hemo;
    const live = c.reduce((a, q) => a + q.live, 0), dead = c.reduce((a, q) => a + q.dead, 0);
    if (!live && !dead) return ($("countResult").textContent = "Click the cells in the corner squares first.");
    const empty = c.filter((q) => !q.live && !q.dead).length;
    const perSq = live / 4, perMl = perSq * 1e4 * h.dilution;
    const viab = live + dead ? (live / (live + dead)) * 100 : 0;
    const truth = s.field.truth, tl = truth.reduce((a, q) => a + q.live, 0), td = truth.reduce((a, q) => a + q.dead, 0);
    const html = `Live cells: ${live} in 4 squares = ${perSq.toFixed(1)} per square.<br>× 10<sup>4</sup> (each square holds 0.1 µL) × ${h.dilution} (${h.trypan ? "diluted 1:1 in trypan blue" : "dilution"}) = <b>${fmt(perMl)} live cells per mL</b>.<br>
      ${h.trypan ? `Viability: ${live} ÷ ${live + dead} = <b>${viab.toFixed(0)} %</b>.` : "No trypan blue: you can't tell live from dead."}<br>
      <span class="muted">In the squares there were ${tl} live and ${td} dead cells. The suspension really holds ${fmt(h.livePerMl * h.dilution)} live cells per mL, from the CellSim engine's dish at ${h.t} h.
      Cells land in the chamber at random, so a count of ${live} is uncertain by about ±${Math.round(Math.sqrt(live))} (the square root of the count); counting more squares narrows it.</span>
      ${empty ? `<br>Count all four squares for a fair average (${empty} still empty).` : ""}`;
    $("countResult").innerHTML = html;
    onEvent("count", { live, dead, perMl, viab, truePerMl: h.livePerMl * h.dilution, t: h.t, dish: h.dish });
  });
  function fmt(n) { if (n <= 0) return "0"; const e = Math.floor(Math.log10(n)); return `${(n / 10 ** e).toFixed(1)} × 10<sup>${e}</sup>`; }

  /* ---------- wells of a plate ---------- */
  function renderWells() {
    const sp = s.spec;
    $("wellPick").innerHTML = sp.wells.map((w, i) => `<button data-well="${i}" aria-pressed="${i === sp.well}">${w}</button>`).join("");
  }
  $("wellPick").addEventListener("click", (e) => {
    const b = e.target.closest("[data-well]");
    if (b) onEvent("well", +b.dataset.well);
  });

  $("measureBtn").addEventListener("click", () => {
    const show = $("measurePanel").hidden;
    $("measurePanel").hidden = !show; $("measureLine").hidden = !show;
    if (show && s.spec.prep.specimen) $("mAsk").textContent = SPECIMENS[s.spec.prep.specimen].measureWhat === "long" ? "About how many cells fit end to end along it?" : "About how many cells would fit side by side along it?";
  });
  $("measureForm").addEventListener("submit", (e) => {
    e.preventDefault();
    const n = parseFloat($("mCount").value), spc = s.spec.prep.specimen;
    if (!(n > 0)) { $("mResult").textContent = "Type how many cells you counted, for example 2 or 7.5."; return; }
    const sp = SPECIMENS[spc], size = Math.round(fov() / n), [lo, hi] = sp.sizeRange;
    const close = size >= lo * 0.7 && size <= hi * 1.3;
    $("mResult").innerHTML = `${fov()} µm ÷ ${n} = <b>about ${size} µm</b>. Real ${sp.name.toLowerCase()} cells are ${sp.sizeText}. ${close ? "Close!" : "Quite far off. Count again, or use a higher lens."}`;
    observe(spc, { size });
  });
  $("scopeExit").addEventListener("click", () => onExit());

  function update(dt) {
    if (!s.field) return;
    const view = { cx: s.cx, cy: s.cy, fov: fov(), t: s.t };
    if (s.field.animated(view)) { s.t += dt; view.t = s.t; s.dirty = true; }
    if (!s.dirty) return;
    s.dirty = false;
    s.field.draw(ctx, W, view);
    // A ring round the cell you picked, drawn in the field's own µm space.
    if (s.sel) {
      const k = W / view.fov;
      ctx.lineWidth = Math.max(1.5, 2 / k); ctx.strokeStyle = "#f0b44c";
      ctx.beginPath(); ctx.arc(s.sel.x, s.sel.y, Math.max(14, s.sel.w || 10), 0, Math.PI * 2); ctx.stroke();
    }
    if (s.field.osmosis && s.field.lastShrunk) {
      const { shrunk, of } = s.field.lastShrunk;
      $("osmoReadout").innerHTML = `In view: <b>${shrunk}</b> of ${of} cells plasmolysed · under the cover slip: ${s.spec.medium}`;
      const now = performance.now();
      if (of >= 4 && inFocus() && s.obj >= 10 && now - s.osmoT > 400) { s.osmoT = now; onEvent("osmo", shrunk / of); }
    }
  }
  // state is for automated checks, like window.lab.
  return { open, update, get spec() { return s.spec; }, get state() { return s; } };
}
