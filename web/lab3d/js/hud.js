// Heads-up display: hotbar, prompts, task tracker, toasts.
const $ = (id) => document.getElementById(id);

export const HUD = {
  hotbar(inv, selected, iconFor) {
    $("hotbar").innerHTML = inv.slots.slice(0, 9).map((it, i) => slotHTML(it, i, i === selected, iconFor)).join("");
  },
  inventoryGrids(inv, pick, iconFor) {
    $("invGrid").innerHTML = inv.slots.slice(9).map((it, i) => slotHTML(it, i + 9, pick === i + 9, iconFor, true)).join("");
    $("invHotbar").innerHTML = inv.slots.slice(0, 9).map((it, i) => slotHTML(it, i, pick === i, iconFor, true)).join("");
  },
  prompt(actions, holding) {
    const el = $("prompt");
    if (!actions || (!actions.left && !actions.right)) { el.innerHTML = ""; return; }
    let h = actions.target ? `<span class="target">${esc(actions.target)}</span>` : "";
    if (actions.right) h += `<span class="p${actions.right.blocked ? " blocked" : ""}"><span class="mouse">Right</span>${esc(actions.right.label)}</span>`;
    if (actions.left) h += `<span class="p"><span class="mouse">Left</span>${esc(actions.left.label)}</span>`;
    el.innerHTML = h;
  },
  held(item) {
    const el = $("heldInfo");
    if (!item) { el.innerHTML = ""; return; }
    const info = item.info();
    el.innerHTML = `<b>${esc(item.name)}</b>${info ? esc(info) : ""}`;
  },
  itemName(name) {
    const el = $("itemName");
    el.textContent = name || "";
    el.classList.toggle("show", !!name);
    clearTimeout(this._nameT);
    this._nameT = setTimeout(() => el.classList.remove("show"), 1600);
  },
  toast(text, kind = "") {
    const box = $("toasts");
    const t = document.createElement("div");
    t.className = `toast ${kind}`;
    t.textContent = text;
    box.prepend(t);
    while (box.children.length > 4) box.lastChild.remove();
    setTimeout(() => t.remove(), 4200);
  },
  tasks(view) {
    const el = $("tasks");
    if (!view) { el.innerHTML = ""; return; }
    const lis = view.steps.map((s, i) => {
      const cls = s.done ? "done" : i === view.now ? "now" : i > view.now + 2 ? "later" : "";
      return `<li class="${cls}">${esc(s.text)}</li>`;
    }).join("");
    const doneN = view.steps.filter((s) => s.done).length;
    el.innerHTML = `<p class="sub">Experiment · ${doneN}/${view.steps.length}</p><h3>${esc(view.title)}</h3><ol>${lis}</ol>${view.finished ? `<p class="sub" style="color:var(--good)">Complete. Well done.</p>` : ""}`;
  },
  engine(state, text) {
    const el = $("engineChip");
    el.hidden = !text;
    el.className = `engine-chip ${state}`;
    el.textContent = text;
  },
  hold(p) {
    const r = $("holdRing");
    r.hidden = p == null;
    if (p != null) r.style.setProperty("--p", p.toFixed(3));
  },
};

function slotHTML(it, i, sel, iconFor, clickable) {
  const lvl = it && it.container ? `<span class="lvl"><i style="width:${(it.container.fraction * 100).toFixed(0)}%;--c:${levelColor(it)}"></i></span>` : "";
  return `<div class="slot${sel ? " sel" : ""}${clickable && sel ? " pick" : ""}" data-slot="${i}" role="listitem" title="${it ? esc(it.name) : "Empty"}">
    ${i < 9 ? `<span class="n">${i + 1}</span>` : ""}${it ? `<img src="${iconFor(it)}" alt="${esc(it.name)}">` : ""}${lvl}</div>`;
}

function levelColor(it) {
  const { color } = it.container.color();
  return `#${color.getHexString()}`;
}

function esc(s) { return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]); }
