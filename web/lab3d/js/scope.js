// Looking through the eyepiece: objective lenses, coarse and fine focus,
// light and diaphragm, moving the slide, measuring, and inspecting a
// living cell's engine state.
import { makeField, SPECIMENS, FOV_UM } from "./specimens.js";
import { makeCultureField } from "./culture.js";

const $ = (id) => document.getElementById(id);
const FOCUS_AT = { 4: 60, 10: 61.5, 40: 62.4 };  // nearly parfocal, like a real microscope
const DEPTH = { 4: 5, 10: 2.6, 40: 1.1 };        // depth of field shrinks at high power
const STEP = { coarse: 2.5, fine: 0.25 };

export function createScope({ onFocus, onExit, toast, observe, resolution = 720 }) {
  const canvas = $("field"), ep = $("eyepiece"), msg = $("epMsg");
  const W = resolution;
  canvas.width = W; canvas.height = W;
  const ctx = canvas.getContext("2d");
  const s = { spec: null, field: null, obj: 4, focus: 28, light: 0.7, diaphragm: 0.35, cx: 0, cy: 0, t: 0, dirty: true, done: {}, told: {} };

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
      $("cellCard").hidden = true; $("measureBtn").hidden = true;
      return;
    }
    msg.textContent = "Nothing here. Drag to move the slide.";
    s.spec = spec;
    s.field = spec.kind === "culture" ? makeCultureField(spec.header, spec.snap, { tint: spec.tint }) : makeField(spec.prep);
    Object.assign(s, { obj: 4, focus: 28, cx: 0, cy: 0, t: 0, dirty: true, done: {}, told: {} });
    $("scopeTitle").textContent = spec.title;
    $("cultureTime").textContent = spec.kind === "culture" ? ` · ${spec.time} h in culture` : "";
    $("measurePanel").hidden = true;
    $("measureLine").hidden = true;
    $("cellCard").hidden = true;
    $("mCount").value = ""; $("mResult").textContent = "";
    $("measureBtn").hidden = spec.kind === "culture" || !spec.prep.specimen;
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
        toast(p.specimen === "cheek" ? "Hard to see? Cheek cells are almost clear. Methylene blue would stain them."
          : p.specimen === "onion" ? "The nuclei are faint. Iodine stains onion cells best." : "A stain darkens the chloroplasts. Plain water shows the green best.", "warn");
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
    const T = [[4, "Focus with the 4× lens (coarse knob)"], [10, "Switch to 10× and sharpen"], [40, "Switch to 40×, fine knob only"]];
    $("scopeTasks").innerHTML = T.map(([k, t]) => `<li class="${s.done[k] ? "done" : ""}">${t}</li>`).join("") +
      (s.spec && s.spec.kind === "culture" ? `<li>Click a cell to see what the engine says about it</li>` : "");
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

  // Drag to move the slide; a click without dragging inspects a cell.
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
    if (drag && drag.moved < 5 && s.field && s.field.pickCell) inspect(e);
    drag = null; check();
  });

  function inspect(e) {
    const r = ep.getBoundingClientRect();
    const k = fov() / r.width;
    const x = s.cx + (e.clientX - r.left - r.width / 2) * k, y = s.cy + (e.clientY - r.top - r.height / 2) * k;
    const c = s.field.pickCell(x, y);
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
  }
  return { open, update };
}
