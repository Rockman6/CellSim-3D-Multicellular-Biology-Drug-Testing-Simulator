// Every object you can take or use, and what it does. getActions() decides
// what left- and right-click do for what you hold and what you look at.
import * as THREE from "three";
import * as P from "./props.js";
import { Container, SOLUTIONS, fmtVol } from "./liquids.js";
import { fit, enableShadows, link } from "./assets.js";

const DROP = 0.05; // mL in one drop from a dropper or pipette
const HIT_MAT = new THREE.MeshBasicMaterial({ visible: false });

let uid = 0;
export class Item {
  constructor(kind, name, obj, extra = {}) {
    this.id = `${kind}-${++uid}`;
    this.kind = kind;
    this.name = name;
    this.obj = obj;
    this.state = {};
    this.pickable = true;
    Object.assign(this, extra);
    link(obj, "item", this);
    // An invisible box a little wider than the item, so thin things like
    // forceps or a cover slip are easy to aim at.
    const b = new THREE.Box3().setFromObject(obj);
    if (!b.isEmpty()) {
      const size = b.getSize(new THREE.Vector3()), c = b.getCenter(new THREE.Vector3());
      const hit = new THREE.Mesh(new THREE.BoxGeometry(Math.max(size.x + 0.012, 0.03), Math.max(size.y + 0.003, 0.012), Math.max(size.z + 0.012, 0.03)), HIT_MAT);
      hit.position.set(c.x, b.min.y + Math.max(size.y + 0.003, 0.012) / 2, c.z);
      hit.name = "hitbox";
      obj.add(hit);
    }
    obj.traverse((m) => { if (m.isMesh) link(m, "itemRoot", obj); });
  }
  refresh() {
    if (this.liquid) this.liquid.update(this.container);
    if (this.kind === "slide") P.updateSlideLook(this);
  }
  info() {
    if (this.kind === "slide") return slideInfo(this);
    if (this.kind === "dish") return `${this.container.describe()} · ${this.state.t} h in culture${drugText(this.container)}`;
    if (this.container) return this.container.describe();
    if (this.kind === "forceps" && this.state.carry) return `Holding ${this.state.carry === "onion" ? "a thin layer of onion skin" : "one pondweed leaf"}`;
    if (this.kind === "swab") return this.state.used ? "Used: drop it in the disinfectant" : this.state.cells ? "Has cheek cells on it" : "Clean";
    if (this.kind === "slidebox") return `${this.state.left} clean slides left`;
    if (this.kind === "coverbox") return `${this.state.left} cover slips left`;
    return "";
  }
}

// Lower-case a name for a sentence, but keep names like "Dish A" as they are.
function lower(it) { return it.kind === "dish" ? it.name : it.name.toLowerCase(); }

function drugText(c) {
  const d = c.drugs();
  const parts = Object.entries(d).filter(([, v]) => v > 1e-4).map(([k, v]) => `${k} ${v < 0.1 ? v.toFixed(3) : v.toFixed(2)} µM`);
  return parts.length ? ` · ${parts.join(", ")}` : "";
}

function slideInfo(s) {
  const st = s.state, bits = [];
  bits.push(st.sample ? { onion: "Onion skin", cheek: "Cheek cells", elodea: "Pondweed leaf" }[st.sample] : "No sample");
  if (!s.container.empty) bits.push(s.container.describe());
  bits.push(st.cover ? (st.cover === "angle" ? "cover slip, no bubbles" : "cover slip, air bubbles") : "no cover slip");
  return bits.join(" · ");
}

/* ---------- making items ---------- */
function liquidItem(kind, name, made, extra) {
  return new Item(kind, name, made.obj, { container: made.container, liquid: made.liquid, ...extra });
}

export function newSlide() {
  const it = new Item("slide", "Glass slide", P.makeSlide(), { container: new Container(0.4) });
  it.state = { sample: null, cover: null, label: "" };
  return it;
}
export function newCoverSlip() { return new Item("coverslip", "Cover slip", P.makeCoverSlip()); }
export function newSwab() { const it = new Item("swab", "Cotton swab", P.makeSwab()); it.state = { cells: false, used: false }; return it; }
export function newPencil() { return new Item("pencil", "Pencil", P.makePencil()); }

export function createItems(game, assets, room) {
  const Y = room.TABLE_Y, B = room.BENCH.h;
  const items = [];
  const put = (it, x, y, z, ry = 0) => { it.obj.position.set(x, y, z); it.obj.rotation.y = ry; enableShadows(it.obj); game.scene.add(it.obj); it.refresh(); items.push(it); return it; };

  // Slide preparation area, along the near edge of the big bench.
  put(new Item("tile", "White tile", P.makeTile()), -2.3, Y, -1.66);
  const sbox = put(new Item("slidebox", "Box of slides", P.makeBox({ w: 0.085, h: 0.03, d: 0.035, color: 0x2a5fb0, labelLines: ["Microscope slides", "76 × 26 mm"] })), -2.75, Y, -1.62);
  sbox.state.left = 10;
  const cbox = put(new Item("coverbox", "Cover slips", P.makeBox({ w: 0.032, h: 0.014, d: 0.032, color: 0xdfe6ea, labelLines: ["Cover slips"] })), -2.6, Y, -1.9);
  cbox.state.left = 20;
  put(new Item("towels", "Paper towels", P.makePaperTowels()), -1.95, Y, -1.98);
  const stains = [["water", 25], ["saline", 25], ["iodine", 22], ["mb", 22]];
  stains.forEach(([s, v], i) => put(liquidItem("dropper", `${P.shortName(s)} dropper bottle`, P.makeDropperBottle(s, v), { solution: s }), -1.82 + i * 0.075, Y, -1.66));
  const onionObj = enableShadows(fit(assets.gltf.yellow_onion, { height: 0.085 }));
  put(new Item("onion", "Onion", onionObj), -3.25, Y, -1.7).state.cut = false;
  put(new Item("scalpel", "Scalpel", P.makeScalpel()), -3.02, Y, -1.52, 0.3);
  put(new Item("forceps", "Forceps", P.makeForceps()), -2.95, Y, -1.45, 0.2);
  put(new Item("swabcup", "Cup of cotton swabs", P.makeSwabCup()), -3.55, Y, -1.66);
  put(liquidItem("beaker", "Disinfectant beaker", P.makeBeaker({ disinfectant: 150 })), -3.78, Y, -1.95);
  put(liquidItem("beaker", "Empty beaker", P.makeBeaker({})), -1.25, Y, -1.92);
  put(new Item("tank", "Pondweed tank", P.makeTank(), { pickable: false }), -4.5, Y, -0.55, Math.PI / 2);
  const lens = enableShadows(fit(assets.gltf.magnifying_glass_01, { maxDim: 0.27, rotX: -Math.PI / 2 }));
  put(new Item("lens", "Hand lens", lens), -1.15, Y, -1.6, 0.5);
  const cup = enableShadows(fit(assets.gltf.stationery_supplies, { maxDim: 0.3 }));
  put(new Item("pencilcup", "Pens and pencils", cup, { pickable: false }), -0.95, Y, -1.9, 0.4);
  put(new Item("sanitizer", "Hand sanitizer", P.makeSanitizer(), { pickable: false }), 0.05, Y, -2.35);

  // Stock bottles on the steel shelves.
  const shelfLevels = findShelfLevels(game, room.shelves);
  const stocks = [["water", 480], ["saline", 450], ["iodine", 400], ["mb", 400], ["ethanol", 450]];
  stocks.forEach(([s, v], i) => {
    const lvl = shelfLevels[Math.min(shelfLevels.length - 1, 1 + Math.floor(i / 3))];
    put(liquidItem("stock", `${P.shortName(s)} stock bottle`, P.makeStockBottle(s, v), { solution: s }), -3.95 + (i % 3) * 0.32, lvl, 3.66);
  });

  // The microscope bench: drug stocks and the goggles.
  [["paclitaxel", -2.1], ["doxorubicin", -1.98], ["cisplatin", -1.86]].forEach(([d, x]) =>
    put(liquidItem("drug", `${P.shortName(d)} vial`, P.makeDrugVial(d, 8), { solution: d }), x, B, 0.68));
  put(new Item("goggles", "Safety goggles", P.makeGoggles()), -0.42, B, 1.15, -0.4);

  // On the teacher's desk: the journal and the method sheet.
  const journal = enableShadows(fit(assets.gltf.binder_notebook, { maxDim: 0.32 }));
  journal.traverse((n) => { if (n.name === "binder_notebook") n.visible = false; });
  put(new Item("journal", "Lab journal", journal, { pickable: false }), 3.85, room.deskTop, -2.75, 0.2);
  const clip = enableShadows(fit(assets.gltf.clipboard, { maxDim: 0.34, rotX: -Math.PI / 2 }));
  put(new Item("clipboard", "Method sheet", clip, { pickable: false }), 3.8, room.deskTop, -1.95, -0.3);

  // Two culture dishes waiting in the incubator: A to treat, B as the control.
  const dishes = ["Dish A", "Dish B"].map((name, i) => {
    const it = liquidItem("dish", name, P.makeCultureDish(name), {});
    it.state = { label: name, seed: 11, t: 0, events: [], snap: null, header: null, series: [], line: "A549" };
    enableShadows(it.obj);
    it.refresh();
    items.push(it);
    return it;
  });
  return { items, dishes };
}

// Shelf heights, found by dropping a ray down through the shelf unit.
function findShelfLevels(game, shelves) {
  const ray = new THREE.Raycaster();
  const b = new THREE.Box3().setFromObject(shelves);
  const x = (b.min.x + b.max.x) / 2 - 0.2, z = (b.min.z + b.max.z) / 2;
  const levels = [];
  let y = b.max.y + 0.05;
  for (let k = 0; k < 6 && y > 0.1; k++) {
    ray.set(new THREE.Vector3(x, y, z), new THREE.Vector3(0, -1, 0));
    const h = ray.intersectObject(shelves, true)[0];
    if (!h) break;
    levels.push(h.point.y);
    y = h.point.y - 0.06;
  }
  return levels.sort((a, b) => a - b).filter((v, i, a) => i === 0 || v - a[i - 1] > 0.15);
}

/* ---------- what clicking does ---------- */
export function getActions(game, held, hit) {
  let t = hit && hit.item, f = hit && hit.furniture;
  const out = { left: null, right: null, target: t ? t.name : f ? f.userData.target.name : "" };
  if (t && t.pickable) out.left = { label: `Take ${lower(t)}`, run: () => game.take(t) };
  // A slide or dish on the stage counts as the microscope for right-click.
  if (t && t === game.stageItem) {
    out.target = `Microscope · ${t.name}`;
    out.right = held ? null : furnitureAction(game, game.room.scope);
    out.left = { label: `Take the ${lower(t)} off the stage`, run: () => game.takeFromStage() };
    return out;
  }
  if (f && f.userData.target.kind === "microscope" && game.stageItem) out.left = { label: `Take the ${lower(game.stageItem)} off the stage`, run: () => game.takeFromStage() };
  if (f && f.userData.target.kind === "drawer" && f.userData.open && f.userData.spare) out.left = { label: `Take ${f.userData.spare.label}`, run: () => game.takeSpare(f) };
  out.right = rightAction(game, held, t, f, hit);
  // "the Dish A" reads badly: names don't take "the".
  for (const a of [out.left, out.right]) if (a) a.label = a.label.replace(/\bthe (Dish [AB])\b/g, "$1");
  return out;
}

function rightAction(game, held, t, f, hit) {
  const k = held && held.kind;
  // Using what you hold on an item.
  if (held && t) {
    const r = useOnItem(game, held, t);
    if (r) return r;
  }
  if (held && f) {
    const r = useOnFurniture(game, held, f);
    if (r) return r;
  }
  if (f && !held) return furnitureAction(game, f);
  if (t && !held) {
    const r = itemAction(game, t);
    if (r) return r;
  }
  // Boxes, cups and fixed equipment work whatever you're holding.
  if (t && held && ["slidebox", "coverbox", "swabcup", "pencilcup", "sanitizer", "journal", "clipboard", "goggles"].includes(t.kind)) {
    const r = itemAction(game, t);
    if (r) return r;
  }
  if (f && held && !["microscope", "incubator", "glassbin", "sink"].includes(f.userData.target.kind)) return furnitureAction(game, f);
  // Something in your hand that works on its own.
  if (held) {
    if (k === "goggles") return { label: "Put on the safety goggles", run: () => game.wearGoggles(held) };
    if (k === "swab" && !held.state.cells && !held.state.used) return { label: "Rub the inside of your cheek", run: () => { held.state.cells = true; game.emit("swabbed"); game.toast("Cheek cells collected on the swab"); } };
    if (k === "lens") return { label: "Look closely (hold)", whileHeld: true, start: () => (game.player.zoom = 3.2), end: () => (game.player.zoom = 1) };
  }
  // Putting it down.
  if (held && hit && hit.surface) return { label: `Put the ${lower(held)} down`, run: () => game.place(held, hit.point) };
  if (f && held) return furnitureAction(game, f);
  if (t && held) { const r = itemAction(game, t); if (r) return r; }
  return null;
}

function useOnItem(game, held, t) {
  const k = held.kind;
  const slide = t.kind === "slide";
  if (k === "slide" && t.kind === "tile") return { label: "Put the slide on the white tile", run: () => { game.place(held, tileTop(t)); game.emit("slide-on-tile"); } };
  if (slide && k === "forceps" && held.state.carry) return { label: `Lay the ${held.state.carry === "onion" ? "onion skin" : "leaf"} flat on the slide`, run: () => {
    if (t.state.sample) return game.toast("This slide already has a sample. Rinse it at the sink first.", "warn");
    t.state.sample = held.state.carry; held.state.carry = null; t.refresh(); game.emit(`sample:${t.state.sample}`);
    if (t.state.sample === "elodea" && t.container.empty) game.toast("Tip: a leaf goes into a drop of water, or it dries out.", "warn");
  } };
  if (t.kind === "onion" && k === "scalpel") return { label: "Cut a small piece from the onion", run: () => { t.state.cut = true; game.emit("onion-cut"); game.toast("Onion cut. Now peel the thin inner skin with the forceps."); } };
  if (t.kind === "onion" && k === "forceps") return { label: "Peel a thin layer of skin", run: () => {
    if (!t.state.cut) return game.toast("Cut the onion with the scalpel first.", "warn");
    if (held.state.carry) return game.toast("The forceps are already holding something.", "warn");
    held.state.carry = "onion"; game.toast("A see-through layer of onion skin, one cell thick");
  } };
  if (t.kind === "tank" && k === "forceps") return { label: "Pull off one small leaf", run: () => {
    if (held.state.carry) return game.toast("The forceps are already holding something.", "warn");
    held.state.carry = "elodea"; game.toast("One young pondweed leaf, two cells thick");
  } };
  if (slide && k === "swab" && held.state.cells) return { label: "Smear the swab across the slide", run: () => {
    if (t.state.sample) return game.toast("This slide already has a sample.", "warn");
    t.state.sample = "cheek"; held.state.cells = false; held.state.used = true; t.refresh(); game.emit("sample:cheek");
    game.toast("Now drop the used swab in the disinfectant.");
  } };
  if (t.kind === "beaker" && k === "swab" && t.container.c.disinfectant) return { label: "Drop the swab in the disinfectant", run: () => { game.dispose(held); game.emit("dispose:swab"); } };
  if (slide && (k === "dropper" || k === "drug")) {
    if (k === "drug") return { label: "Drugs go in the culture dish", run: () => game.toast("Cytotoxic drugs are only for the living cells in a culture dish.", "warn") };
    return { label: t.state.cover ? `Add a drop of ${P.shortName(held.solution).toLowerCase()} at the edge of the cover slip` : `Add one drop of ${P.shortName(held.solution).toLowerCase()}`, run: () => {
      if (held.container.empty) return game.toast("This bottle is empty. Refill it from the stock bottle on the shelf.", "warn");
      const before = t.container.volume;
      held.container.pour(t.container, DROP);
      if (t.container.volume - before < 1e-6) return game.toast("The slide can't hold any more. Blot it with a paper towel.", "warn");
      held.refresh(); t.refresh();
      game.emit(`drop:${held.solution}`);
      if (t.state.cover) game.toast("The stain runs under the cover slip");
    } };
  }
  if (t.kind === "dish" && k === "drug") return { label: `Add one drop of ${P.shortName(held.solution).toLowerCase()} to ${t.name}`, run: () => {
    if (held.container.empty) return game.toast("The vial is empty.", "warn");
    held.container.pour(t.container, DROP);
    held.refresh(); t.refresh();
    game.dishDosed(t);
  } };
  if (t.kind === "dish" && k === "dropper") return { label: "Stains are for slides", run: () => game.toast("A stain would kill the living cells. Use it on slides only.", "warn") };
  if (slide && k === "coverslip") return { label: "Lower the cover slip: hold for slow, tap to drop", hold: 1.1, run: (slow) => {
    if (t.state.cover) return game.toast("This slide already has a cover slip.", "warn");
    t.state.cover = slow ? "angle" : "flat"; game.dispose(held, { quiet: true }); t.refresh();
    game.emit(`cover:${t.state.cover}`);
    game.toast(slow ? "Lowered slowly at an angle: no air bubbles" : "Dropped flat: that trapped air bubbles", slow ? "good" : "warn");
  } };
  if (slide && k === "pencil") return { label: "Label the frosted end", run: () => {
    t.state.label = t.state.sample ? { onion: "Onion", cheek: "Cheek", elodea: "Leaf" }[t.state.sample] : "Slide"; t.refresh(); game.emit("label"); } };
  if (slide && k === "towels") return { label: "Blot away the extra liquid", run: () => {
    const extra = t.container.volume - DROP;
    if (extra <= 0.001) return game.toast("There's no extra liquid to blot.");
    t.container.take(extra); t.refresh(); game.toast("Blotted. The cover slip sits flat now.");
  } };
  if ((k === "stock" && t.kind === "dropper") || (k === "dropper" && t.kind === "stock")) {
    const stock = k === "stock" ? held : t, bottle = k === "stock" ? t : held;
    if (stock.solution !== bottle.solution) return { label: "Different solutions", run: () => game.toast("Only refill a bottle from the stock of the same solution.", "warn") };
    return { label: "Refill the dropper bottle from the stock", run: () => {
      const v = stock.container.pour(bottle.container, bottle.container.capacity - bottle.container.volume);
      stock.refresh(); bottle.refresh(); game.toast(`Refilled with ${fmtVol(v)}`);
    } };
  }
  if ((k === "stock" && t.kind === "beaker") || (k === "beaker" && t.kind === "stock")) {
    const from = k === "stock" ? held : t, to = k === "stock" ? t : held;
    return { label: `Pour 20 mL of ${P.shortName(from.solution).toLowerCase()} into the beaker`, run: () => {
      const v = from.container.pour(to.container, 20); from.refresh(); to.refresh();
      game.toast(v > 0 ? `Poured ${fmtVol(v)}. The beaker now holds ${to.container.describe()}` : "The beaker is full.");
    } };
  }
  if (k === "dropper" && t.kind === "beaker") return { label: `Add a drop to the beaker`, run: () => { held.container.pour(t.container, DROP); held.refresh(); t.refresh(); } };
  return null;
}

function tileTop(tile) { return tile.obj.position.clone().add(new THREE.Vector3(0, 0.006, 0)); }

function useOnFurniture(game, held, f) {
  const kind = f.userData.target.kind, k = held.kind;
  if (kind === "microscope" && (k === "slide" || k === "dish") && !game.stageItem) return { label: `Clip the ${lower(held)} onto the stage`, run: () => game.putOnStage(held) };
  if (kind === "incubator" && k === "dish") return { label: `Put ${held.name} in the incubator`, run: () => game.putInIncubator(held) };
  if (kind === "glassbin" && (k === "slide" || k === "coverslip")) return { label: `Throw the ${lower(held)} into the glass waste`, run: () => { game.dispose(held); game.emit(`dispose:${k}`); } };
  if (kind === "sink" && k === "slide") return { label: "Rinse the slide clean", run: () => {
    held.state = { sample: null, cover: null, label: "" }; held.container = new Container(0.4); held.refresh(); game.runTap(); game.toast("Rinsed and dried. A clean slide again.");
  } };
  if (kind === "sink" && k === "beaker") return { label: "Fill the beaker with tap water", run: () => { held.container.add("tap", 200); held.refresh(); game.runTap(); } };
  return null;
}

function furnitureAction(game, f) {
  const t = f.userData.target;
  switch (t.kind) {
    case "microscope": return game.stageItem
      ? { label: "Look through the eyepiece", run: () => game.enterScope() }
      : { label: "Look through the eyepiece", run: () => game.toast("Clip a slide or a culture dish onto the stage first.", "warn") };
    case "incubator": return { label: "Open the incubator", run: () => game.openIncubator() };
    case "sink": return { label: "Wash your hands", run: () => { game.runTap(); game.washHands(); } };
    case "door": return { label: "Leave the lab", run: () => game.leave() };
    case "switch": return { label: game.lightsOn ? "Turn the lights off" : "Turn the lights on", run: () => game.toggleLights() };
    case "blinds": return { label: f.userData.down ? "Open the blinds" : "Close the blinds", run: () => game.toggleBlinds(f) };
    case "stool": return game.player.sitting === f ? null : { label: "Sit down", run: () => game.player.sit(f) };
    case "clock": return { label: "Check the time", run: () => game.toast(`It's ${new Date().toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" })}`) };
    case "whiteboard": return { label: "Read today's method", run: () => game.openReader() };
    case "drawer": return { label: f.userData.open ? "Close the drawer" : "Open the drawer", run: () => game.toggleDrawer(f) };
    case "glassbin": return { label: "Glass waste: for used slides and cover slips", run: () => game.toast("Hold a used slide or cover slip, then right-click the bin.") };
  }
  return null;
}

function itemAction(game, t) {
  switch (t.kind) {
    case "slidebox": return { label: "Take a clean slide", run: () => {
      if (t.state.left <= 0) return game.toast("The box is empty. There are spares in the drawer cabinet.", "warn");
      t.state.left--; game.give(newSlide()); } };
    case "coverbox": return { label: "Take a cover slip", run: () => {
      if (t.state.left <= 0) return game.toast("No cover slips left. There are spares in the drawer cabinet.", "warn");
      t.state.left--; game.give(newCoverSlip()); } };
    case "swabcup": return { label: "Take a cotton swab", run: () => game.give(newSwab()) };
    case "pencilcup": return { label: "Take a pencil", run: () => game.give(newPencil()) };
    case "goggles": return { label: "Put on the safety goggles", run: () => game.wearGoggles(t) };
    case "sanitizer": return { label: "Clean your hands with sanitizer", run: () => game.washHands() };
    case "journal": return { label: "Open your lab journal", run: () => game.openJournal() };
    case "clipboard": return { label: "Read the method sheet", run: () => game.openReader() };
    case "tank": return { label: "Pondweed: take a leaf with the forceps", run: () => game.toast("Hold the forceps, then right-click the tank.") };
    case "onion": return { label: "Onion: cut it with the scalpel", run: () => game.toast(t.state.cut ? "Peel the thin inner skin with the forceps." : "Hold the scalpel, then right-click the onion.") };
    case "slide": return { label: "Slide: add a sample, a drop, a cover slip", run: () => game.toast("Hold forceps, a swab, a dropper bottle or a cover slip, then right-click the slide.") };
    case "dish": return { label: `${t.name}: ${t.state.t} h in culture`, run: () => game.toast("Clip the dish onto the microscope stage to look at the cells.") };
  }
  return null;
}
