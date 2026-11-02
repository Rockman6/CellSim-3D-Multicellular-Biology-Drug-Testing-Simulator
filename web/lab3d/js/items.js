// Every object you can take or use, and what it does. getActions() decides
// what left- and right-click do for what you hold and what you look at.
import * as THREE from "three";
import * as P from "./props.js";
import { Container, SOLUTIONS, fmtVol } from "./liquids.js";
import { fit, enableShadows, link } from "./assets.js";

const DROP = 0.05; // mL in one drop from a dropper bottle
const HIT_MAT = new THREE.MeshBasicMaterial({ visible: false });
export const SOFTEN_S = 20;   // seconds in warm acid, standing in for 5 minutes

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
    // forceps or a cover slip are easy to aim at, and a soft shadow under it.
    // Big open things (the water bath) get no box: it would hide what's inside them.
    const b = new THREE.Box3().setFromObject(obj);
    if (!b.isEmpty()) {
      const size = b.getSize(new THREE.Vector3()), c = b.getCenter(new THREE.Vector3());
      const hit = extra.hitbox === false ? new THREE.Object3D() : new THREE.Mesh(new THREE.BoxGeometry(Math.max(size.x + 0.012, 0.03), Math.max(size.y + 0.003, 0.012), Math.max(size.z + 0.012, 0.03)), HIT_MAT);
      hit.position.set(c.x, b.min.y + Math.max(size.y + 0.003, 0.012) / 2, c.z);
      hit.name = "hitbox";
      obj.add(hit);
      const sh = P.contactShadow(size.x, size.z);
      sh.position.set(c.x, b.min.y + 0.0004, c.z);
      obj.add(sh);
      link(obj, "blob", sh);
      this.foot = Math.max(0.012, Math.hypot(size.x, size.z) / 2 * 0.8);
    }
    obj.traverse((m) => { if (m.isMesh) link(m, "itemRoot", obj); });
  }
  refresh() {
    if (this.liquid) this.liquid.update(this.container);
    if (this.kind === "slide") P.updateSlideLook(this);
    if (this.kind === "pipette") P.updatePipetteLook(this);
    if (this.kind === "plate") for (const w of this.wells) w.liquid.update(w.container);
    if (this.kind === "hemocytometer") this.obj.userData.wet.visible = !this.container.empty;
    if (this.kind === "tipbox" && this.state.left != null) this.obj.userData.tips.count = this.state.left;
  }
  info() {
    if (this.kind === "slide") return slideInfo(this);
    if (this.kind === "dish") return `${this.state.detached ? `Cells in suspension: ${fmtVol(this.container.volume)}` : this.container.describe()} · ${this.state.t} h in culture${drugText(this.container)}`;
    if (this.kind === "plate") { const n = this.state.wells.filter((w) => w.events.length).length; return `${this.state.t} h in culture · ${n} of 6 wells treated`; }
    if (this.kind === "pipette") return pipetteInfo(this);
    if (this.kind === "hemocytometer") return this.container.empty ? "Empty: load 10 µL of stained cells" : `Loaded: ${this.container.describe()}`;
    if (this.container) return this.container.describe();
    if (this.kind === "forceps" && this.state.carry) return `Holding ${CARRY[this.state.carry]}`;
    if (this.kind === "swab") return this.state.used ? "Used: drop it in the disinfectant" : this.state.cells ? "Has cheek cells on it" : "Clean";
    if (this.kind === "slidebox") return `${this.state.left} clean slides left`;
    if (this.kind === "coverbox") return `${this.state.left} cover slips left`;
    if (this.kind === "tipbox") return `${this.state.left} sterile tips left`;
    return "";
  }
}

const CARRY = { onion: "a thin layer of onion skin", elodea: "one pondweed leaf", redonion: "a thin layer of red onion skin", root: "a 5 mm root tip", "root-soft": "a softened root tip" };

// Lower-case a name for a sentence, but keep names like "Dish A" as they
// are, and units and acronyms (1 M, EDTA, A549) in capitals.
export function sentence(name) { return name.replace(/\b([A-Z])([a-z]+)\b/g, (w, a, b) => a.toLowerCase() + b); }
function lower(it) { return it.kind === "dish" || it.kind === "plate" ? it.name : sentence(it.name); }
const lowName = (solution) => sentence(P.shortName(solution));

function drugText(c) {
  const d = c.drugs();
  const parts = Object.entries(d).filter(([, v]) => v > 1e-5).map(([k, v]) => `${k} ${fmtConc(v)}`);
  return parts.length ? ` · ${parts.join(", ")}` : "";
}
export function fmtConc(uM) { return uM >= 0.1 ? `${uM.toFixed(2)} µM` : uM >= 0.001 ? `${(uM * 1000).toFixed(uM >= 0.01 ? 0 : 1)} nM` : `${(uM * 1000).toFixed(2)} nM`; }

function slideInfo(s) {
  const st = s.state, bits = [];
  bits.push(st.sample ? { onion: "Onion skin", cheek: "Cheek cells", elodea: "Pondweed leaf", redonion: "Red onion skin", root: st.softened ? "Softened root tip" : "Root tip" }[st.sample] : "No sample");
  if (!s.container.empty) bits.push(s.container.describe());
  bits.push(st.cover ? (st.cover === "angle" ? "cover slip, no bubbles" : "cover slip, air bubbles") : "no cover slip");
  if (st.squashed) bits.push("squashed");
  if (st.edge) bits.push(`a drop of ${lowName(st.edge)} waiting at the edge`);
  return bits.join(" · ");
}

function pipetteInfo(p) {
  const st = p.state;
  const tip = !st.tip ? "no tip: fit one from the tip box" : st.load && !st.load.empty ? `holding ${fmtVol(st.load.volume)} of ${sentence(st.load.describe().split(" · ")[1])}` : st.used ? `used tip (it touched ${st.usedName})` : "clean tip";
  return `Set to ${Math.round(st.vol)} µL · ${tip} · R: set the volume`;
}

/* ---------- making items ---------- */
function liquidItem(kind, name, made, extra) {
  return new Item(kind, name, made.obj, { container: made.container, liquid: made.liquid, ...extra });
}

export function newSlide() {
  const it = new Item("slide", "Glass slide", P.makeSlide(), { container: new Container(0.4) });
  it.state = { sample: null, cover: null, label: "", osmo: [] };
  return it;
}
export function newCoverSlip() { return new Item("coverslip", "Cover slip", P.makeCoverSlip()); }
export function newSwab() { const it = new Item("swab", "Cotton swab", P.makeSwab()); it.state = { cells: false, used: false }; return it; }
export function newPencil() { return new Item("pencil", "Pencil", P.makePencil()); }
function dropper(s, v) { return liquidItem("dropper", `${P.shortName(s)} dropper bottle`, P.makeDropperBottle(s, v), { solution: s }); }

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
  [["water", 25], ["saline", 25], ["iodine", 22], ["mb", 22], ["salt", 25]].forEach(([s, v], i) => put(dropper(s, v), -1.82 + i * 0.075, Y, -1.66));
  const onionObj = enableShadows(fit(assets.gltf.yellow_onion, { height: 0.085 }));
  put(new Item("onion", "Onion", onionObj), -3.25, Y, -1.7).state.cut = false;
  const red = enableShadows(fit(P.makeRedOnion(assets.gltf.yellow_onion.clone()), { height: 0.08 }));
  put(new Item("redonion", "Red onion", red), -3.42, Y, -1.95, 1.2).state.cut = false;
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

  // Root tip station, on the bench arm by the windows: onion roots, acid, a warm bath.
  const rooted = P.makeRootedOnion(enableShadows(fit(assets.gltf.yellow_onion.clone(), { height: 0.07 })));
  const ro = put(liquidItem("rootonion", "Onion with roots", rooted, { pickable: false }), -4.45, Y, 0.05);
  ro.state = { cut: 0 };
  const bath = put(new Item("waterbath", "Water bath (60 °C)", P.makeWaterBath(), { pickable: false, hitbox: false }), -4.5, Y, 0.62, Math.PI / 2);
  put(liquidItem("watchglass", "Watch glass", P.makeWatchGlass()), -4.5, Y + bath.obj.userData.shelfY, 0.62).state = {};
  put(dropper("hcl", 25), -4.18, Y, 0.36);
  put(dropper("toluidine", 22), -4.18, Y, 0.47);

  // Stock bottles on the steel shelves.
  const shelfLevels = findShelfLevels(game, room.shelves);
  const stocks = [["water", 480], ["saline", 450], ["iodine", 400], ["mb", 400], ["ethanol", 450], ["salt", 450]];
  stocks.forEach(([s, v], i) => {
    const lvl = shelfLevels[Math.min(shelfLevels.length - 1, 1 + Math.floor(i / 3))];
    put(liquidItem("stock", `${P.shortName(s)} stock bottle`, P.makeStockBottle(s, v), { solution: s }), -3.95 + (i % 3) * 0.32, lvl, 3.66);
  });

  // The microscope bench: drug stocks and goggles, and the cell culture kit:
  // micropipette and tips, trypsin, trypan blue, a hemocytometer, a plate reader.
  [["paclitaxel", -2.1], ["doxorubicin", -1.98], ["cisplatin", -1.86]].forEach(([d, x]) =>
    put(liquidItem("drug", `${P.shortName(d)} vial`, P.makeDrugVial(d, 8), { solution: d }), x, B, 0.68));
  put(dropper("trypan", 20), -1.74, B, 0.68);
  put(new Item("goggles", "Safety goggles", P.makeGoggles()), -0.42, B, 1.15, -0.4);
  const pip = put(new Item("pipette", "Micropipette", P.makeMicropipette()), -2.02, B, 1.13, 0.1);
  pip.state = { vol: 20, tip: false, load: null, used: null };
  pip.refresh();
  const tipbox = new Item("tipbox", "Box of pipette tips", P.makeTipBox(), { pickable: false });
  tipbox.state.left = 96;
  put(tipbox, -2.2, B, 0.92);
  put(new Item("tipwaste", "Tip waste", P.makeTipWaste(), { pickable: false }), -1.76, B, 1.19);
  put(liquidItem("trypsin", "Trypsin-EDTA bottle", P.makeStockBottle("trypsin", 80, { capacity: 100, scale: 0.6 }), { solution: "trypsin" }), -2.27, B, 1.2);
  const rack = put(new Item("tuberack", "Tube rack", P.makeTubeRack(), { pickable: false }), -1.95, B, 0.93);
  const tubes = [["TB", "Trypan blue tube", "trypan", 0.01], ["1", "Empty tube", null, 0]];
  tubes.forEach(([lab, name, sol, v], i) => {
    const t = put(liquidItem("tube", name, P.makeMicrotube(lab, sol, v)), 0, 0, 0);
    seatInRack(rack, t, i);
  });
  put(liquidItem("hemocytometer", "Hemocytometer", { obj: P.makeHemocytometer(), container: new Container(0.02), liquid: null }), -1.25, B, 1.18, 0.2);
  put(new Item("platereader", "Plate reader", P.makePlateReader(), { pickable: false }), -1.45, B, 0.78);

  // On the teacher's desk: the journal and the method sheet.
  const journal = enableShadows(fit(assets.gltf.binder_notebook, { maxDim: 0.32 }));
  journal.traverse((n) => { if (n.name === "binder_notebook") n.visible = false; });
  put(new Item("journal", "Lab journal", journal, { pickable: false }), 3.85, room.deskTop, -2.75, 0.2);
  const clip = enableShadows(fit(assets.gltf.clipboard, { maxDim: 0.34, rotX: -Math.PI / 2 }));
  put(new Item("clipboard", "Method sheet", clip, { pickable: false }), 3.8, room.deskTop, -1.95, -0.3);

  // In the incubator: Dish A to treat, Dish B as the control, and a 6-well plate.
  const dishes = ["Dish A", "Dish B"].map((name) => {
    const it = liquidItem("dish", name, P.makeCultureDish(name), {});
    it.state = { label: name, seed: 11, t: 0, events: [], snap: null, header: null, series: [], line: "A549" };
    enableShadows(it.obj);
    it.refresh();
    items.push(it);
    return it;
  });
  const made = P.makeWellPlate("Plate 1");
  const plate = new Item("plate", "Plate 1", made.obj, { wells: made.wells });
  plate.state = { t: 0, seed: 21, line: "A549", wells: made.wells.map(() => ({ events: [], result: null })) };
  enableShadows(plate.obj);
  plate.refresh();
  items.push(plate);
  return { items, dishes, plate };
}

// Put a microtube upright in hole i of the rack.
export function seatInRack(rack, tube, i) {
  const h = P.RACK_HOLES[i];
  rack.obj.updateMatrixWorld(true);
  const p = rack.obj.localToWorld(h.clone().setY(0.012));
  tube.obj.position.copy(p);
  tube.obj.rotation.set(0, 0, 0);
  tube.state.rackHole = i;
  if (tube.obj.userData.blob) tube.obj.userData.blob.visible = false;
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
  if (t && t.kind === "plate" && hit.point) out.target = `${t.name} · well ${P.WELLS[P.wellAt(t.obj, t.obj.worldToLocal(hit.point.clone()))]}`;
  if (t && t.pickable) out.left = { label: `Take ${lower(t)}`, run: () => game.take(t) };
  // A slide or dish on the stage counts as the microscope for right-click.
  if (t && t === game.stageItem) {
    out.target = `Microscope · ${t.name}`;
    out.right = held ? null : furnitureAction(game, game.room.scope);
    out.left = { label: `Take the ${lower(t)} off the stage`, run: () => game.takeFromStage() };
    return tidy(out);
  }
  if (f && f.userData.target.kind === "microscope" && game.stageItem) out.left = { label: `Take the ${lower(game.stageItem)} off the stage`, run: () => game.takeFromStage() };
  if (f && f.userData.target.kind === "drawer" && f.userData.open && f.userData.spare) out.left = { label: `Take ${f.userData.spare.label}`, run: () => game.takeSpare(f) };
  out.right = rightAction(game, held, t, f, hit);
  return tidy(out);
}

// "the Dish A" reads badly: names don't take "the".
function tidy(out) {
  for (const a of [out.left, out.right]) if (a) a.label = a.label.replace(/\bthe (Dish [AB]|Plate 1)\b/g, "$1");
  return out;
}

function rightAction(game, held, t, f, hit) {
  const k = held && held.kind;
  // Using what you hold on an item.
  if (held && t) {
    const r = useOnItem(game, held, t, hit);
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
  // Putting it down: only on clear bench, or a slide onto the tile.
  if (held && hit && hit.surface && (!hit.item || hit.item.kind === "tile")) {
    const block = game.placeBlocked(held, hit.point, hit.item);
    if (block) return { label: `No room here: ${lower(block)} is in the way`, blocked: true, run: () => game.toast("Aim at a clear spot on the bench.", "warn") };
    return { label: `Put the ${lower(held)} down`, place: true, run: () => game.place(held, hit.point) };
  }
  if (k === "pipette" && !t && !f) return { label: `Set the volume (now ${Math.round(held.state.vol)} µL)`, run: () => game.openPipette(held) };
  if (f && held) return furnitureAction(game, f);
  if (t && held) { const r = itemAction(game, t); if (r) return r; }
  return null;
}

/* ---------- the micropipette ---------- */
// Containers a pipette can draw from or dispense into, with what to call them.
function pipetteTarget(t, hit) {
  if (t.kind === "plate") {
    const i = P.wellAt(t.obj, t.obj.worldToLocal(hit.point.clone()));
    return { container: t.wells[i].container, name: `well ${P.WELLS[i]}`, well: i };
  }
  if (t.kind === "hemocytometer") return { container: t.container, name: "the hemocytometer", chamber: true };
  if (t.container && ["drug", "dropper", "stock", "tube", "dish", "beaker", "trypsin", "watchglass"].includes(t.kind)) return { container: t.container, name: t.kind === "dish" ? t.name : `the ${lower(t)}` };
  return null;
}

function pipetteOn(game, p, t, hit) {
  const st = p.state, vol = st.vol / 1000;
  if (t.kind === "tipbox") {
    if (st.tip) return { label: "There's already a tip on: eject it into the tip waste", run: () => game.toast("Aim at the tip waste and right-click to eject the old tip.", "warn") };
    return { label: "Press the pipette into a fresh tip", run: () => {
      if (t.state.left <= 0) return game.toast("The tip box is empty.", "warn");
      t.state.left--; t.refresh();
      Object.assign(st, { tip: true, load: null, used: null, usedName: "" }); p.refresh(); game.sfx("click"); game.emit("tip");
    } };
  }
  if (t.kind === "tipwaste") {
    if (!st.tip) return { label: "Tip waste: used tips go here", run: () => game.toast("Fit a tip from the tip box first.") };
    return { label: "Eject the tip into the tip waste", run: () => { Object.assign(st, { tip: false, load: null, used: null }); p.refresh(); game.sfx("click"); game.emit("tip-eject"); } };
  }
  const tg = pipetteTarget(t, hit);
  if (!tg) return null;
  if (!st.tip) return { label: "Fit a tip first (from the tip box)", run: () => game.toast("Never put a bare pipette into a liquid: press it into a tip from the tip box first.", "warn") };
  const loaded = st.load && !st.load.empty;
  if (!loaded) {
    const what = tg.container.dominant();
    if (!what || tg.container.empty) return { label: `${cap(tg.name)} is empty`, run: () => game.toast("There's nothing to draw up.", "warn") };
    return { label: `Draw up ${Math.round(st.vol)} µL from ${tg.name}`, run: () => {
      // One tip per solution: a used tip would carry the last liquid into this one.
      if (st.used && st.used !== tg.container) return game.toast(`Eject this tip and fit a fresh one first: it has touched ${st.usedName}, and would carry it into ${tg.name}.`, "warn");
      if (tg.container.volume < vol - 1e-9) return game.toast(`There's only ${fmtVol(tg.container.volume)} left in ${tg.name}.`, "warn");
      st.load = new Container(0.2);
      tg.container.pour(st.load, vol);
      st.used = tg.container; st.usedName = tg.name;
      if (t.kind === "dish" && !t.state.detached && st.load.c.cells == null) game.toast("That's only medium: the cells are stuck to the dish. Trypsin lifts them off.", "warn");
      p.refresh(); t.refresh(); game.sfx("click");
      game.emit(`draw:${what}`);
    } };
  }
  return { label: `Dispense ${fmtVol(st.load.volume)} into ${tg.name}`, run: () => {
    if (tg.chamber && t.container.volume >= 0.0095) return game.toast("The chamber is already full: one 10 µL load is enough.", "warn");
    const v = st.load.volume;
    st.load.pour(tg.container, v);
    // After dispensing into something, the tip has touched it too.
    st.used = tg.container; st.usedName = tg.name;
    p.refresh(); t.refresh(); game.sfx("drip");
    if (tg.well != null) game.wellDosed(t, tg.well);
    else if (t.kind === "dish") game.dishDosed(t);
    else if (tg.chamber) { game.emit("hemo-loaded"); game.toast("The liquid runs under the cover glass and fills the chamber."); }
    else if (t.kind === "tube" && t.container.c.trypan && t.container.c.cells) { game.emit("trypan-mixed"); game.toast("Cells and trypan blue, mixed 1:1. Dead cells take up the blue; live ones keep it out."); }
  } };
}
function cap(s) { return s.charAt(0).toUpperCase() + s.slice(1); }

function useOnItem(game, held, t, hit) {
  const k = held.kind;
  const slide = t.kind === "slide";
  if (k === "pipette") { const r = pipetteOn(game, held, t, hit); if (r) return r; }
  if (k === "slide" && t.kind === "tile") return { label: "Put the slide on the white tile", run: () => { game.place(held, tileTop(t)); game.emit("slide-on-tile"); } };
  if (k === "tube" && t.kind === "tuberack") return { label: "Put the tube in the rack", run: () => {
    const used = new Set(game.items.filter((x) => x.kind === "tube" && x.obj.parent && x !== held).map((x) => x.state.rackHole));
    const free = P.RACK_HOLES.findIndex((_, i) => !used.has(i));
    if (free < 0) return game.toast("The rack is full.", "warn");
    game.place(held, t.obj.position.clone()); seatInRack(t, held, free); game.afterAction();
  } };
  if (slide && k === "forceps" && held.state.carry) return { label: `Lay the ${CARRY[held.state.carry].replace(/^(a|one) /, "")} flat on the slide`, run: () => {
    if (t.state.sample) return game.toast("This slide already has a sample. Rinse it at the sink first.", "warn");
    const c = held.state.carry;
    t.state.sample = c === "root-soft" ? "root" : c; t.state.softened = c === "root-soft"; held.state.carry = null; t.refresh();
    game.emit(`sample:${t.state.sample}`);
    if (t.state.sample === "elodea" && t.container.empty) game.toast("Tip: a leaf goes into a drop of water, or it dries out.", "warn");
    if (c === "root") game.toast("This root tip hasn't been softened: its cells will stay stuck together.", "warn");
  } };
  if ((t.kind === "onion" || t.kind === "redonion") && k === "scalpel") return { label: `Cut a small piece from the ${lower(t)}`, run: () => { t.state.cut = true; game.emit(`${t.kind}-cut`); game.sfx("place"); game.toast("Cut. Now peel the thin skin with the forceps."); } };
  if ((t.kind === "onion" || t.kind === "redonion") && k === "forceps") return { label: t.kind === "redonion" ? "Peel the thin red skin from the outside" : "Peel a thin layer of skin", run: () => {
    if (!t.state.cut) return game.toast("Cut it with the scalpel first.", "warn");
    if (held.state.carry) return game.toast("The forceps are already holding something.", "warn");
    held.state.carry = t.kind; game.toast(t.kind === "redonion" ? "A see-through layer of red cells, one cell thick" : "A see-through layer of onion skin, one cell thick");
  } };
  if (t.kind === "tank" && k === "forceps") return { label: "Pull off one small leaf", run: () => {
    if (held.state.carry) return game.toast("The forceps are already holding something.", "warn");
    held.state.carry = "elodea"; game.toast("One young pondweed leaf, two cells thick");
  } };
  if (t.kind === "rootonion" && k === "scalpel") return { label: "Cut off the last 5 mm of a root", run: () => {
    const roots = t.obj.userData.roots;
    if (t.state.cut >= roots.length) return game.toast("No roots left to cut.", "warn");
    roots[t.state.cut].scale.y = 0.72; roots[t.state.cut].position.y = (0.075 - 0.008) * 0.28;   // shorter from the tip, still attached at the top
    t.state.cut++; t.state.loose = (t.state.loose || 0) + 1;
    game.emit("root-cut"); game.toast("Cut. The 5 mm tip is where the cells divide. Pick it up with the forceps.");
  } };
  if (t.kind === "rootonion" && k === "forceps") return { label: "Pick up the cut root tip", run: () => {
    if (!t.state.loose) return game.toast("Cut a root tip with the scalpel first.", "warn");
    if (held.state.carry) return game.toast("The forceps are already holding something.", "warn");
    t.state.loose--; held.state.carry = "root"; game.toast("Holding a 5 mm root tip. Soften it in the warm acid.");
  } };
  if (t.kind === "watchglass") {
    if (k === "dropper" && held.solution === "hcl") return { label: "Add a drop of 1 M hydrochloric acid", run: () => {
      held.container.pour(t.container, DROP); held.refresh(); t.refresh(); game.sfx("drip");
      if (t.container.volume >= 0.149) game.emit("acid");
    } };
    if (k === "forceps" && held.state.carry === "root") return { label: "Put the root tip in the warm acid", run: () => {
      if (!t.container.c.hcl || t.container.volume < 0.099) return game.toast("Add at least 2 drops of hydrochloric acid to the watch glass first.", "warn");
      if (!game.inWaterBath(t)) return game.toast("Stand the watch glass on the shelf in the water bath first: the acid works at 60 °C.", "warn");
      held.state.carry = null; t.state.root = performance.now(); game.emit("root-acid");
      game.toast(`The acid breaks down what glues the cells together. Leave it ${SOFTEN_S} s (5 minutes, sped up).`);
    } };
    if (k === "forceps" && !held.state.carry && t.state.root) return { label: "Take the root tip out of the acid", run: () => {
      const left = SOFTEN_S - (performance.now() - t.state.root) / 1000;
      if (left > 0) return game.toast(`Not soft yet: ${Math.ceil(left)} s more.`, "warn");
      t.state.root = null; held.state.carry = "root-soft"; game.emit("root-soft");
      game.toast("Soft and pale. Lay it on a slide.");
    } };
  }
  if (slide && k === "swab" && held.state.cells) return { label: "Smear the swab across the slide", run: () => {
    if (t.state.sample) return game.toast("This slide already has a sample.", "warn");
    t.state.sample = "cheek"; held.state.cells = false; held.state.used = true; t.refresh(); game.emit("sample:cheek");
    game.toast("Now drop the used swab in the disinfectant.");
  } };
  if (t.kind === "beaker" && k === "swab" && t.container.c.disinfectant) return { label: "Drop the swab in the disinfectant", run: () => { game.dispose(held); game.emit("dispose:swab"); } };
  if (slide && (k === "dropper" || k === "drug")) {
    if (k === "drug") return { label: "Drugs go in the culture dish", run: () => game.toast("Cytotoxic drugs are only for the living cells in a culture dish.", "warn") };
    if (held.solution === "hcl") return { label: "Acid goes in the watch glass", run: () => game.toast("The acid softens a root tip in the watch glass, in the water bath.", "warn") };
    return { label: t.state.cover ? `Add a drop of ${lowName(held.solution)} at the edge of the cover slip` : `Add one drop of ${lowName(held.solution)}`, run: () => {
      if (held.container.empty) return game.toast("This bottle is empty. Refill it from the stock bottle on the shelf.", "warn");
      const before = t.container.volume;
      held.container.pour(t.container, DROP);
      if (t.container.volume - before < 1e-6) return game.toast("The slide can't hold any more. Blot it with a paper towel.", "warn");
      held.refresh(); t.refresh(); game.sfx("drip");
      game.emit(`drop:${held.solution}`);
      if (t.state.cover) {
        t.state.edge = held.solution;
        game.toast(t.state.sample === "redonion" ? "Now touch a paper towel to the opposite edge to draw it under the cover slip." : "The stain runs under the cover slip");
      }
    } };
  }
  if (t.kind === "dish" && k === "drug") return { label: `Add one drop of ${lowName(held.solution)} to ${t.name}`, run: () => {
    if (held.container.empty) return game.toast("The vial is empty.", "warn");
    held.container.pour(t.container, DROP);
    held.refresh(); t.refresh(); game.sfx("drip");
    game.dishDosed(t);
  } };
  if (t.kind === "dish" && k === "dropper") return { label: "Stains are for slides", run: () => game.toast("A stain would kill the living cells. Use it on slides only.", "warn") };
  if (t.kind === "dish" && k === "trypsin") return { label: `Add trypsin to lift the cells off ${t.name}`, run: () => {
    if (t.state.detached) return game.toast("The cells are already in suspension.");
    // Pour off the medium, trypsin for 5 minutes, then fresh medium stops it.
    held.container.take(0.5); held.refresh();
    t.container = new Container(5.2, { cells: 1 });
    t.container.meta = { dish: t.id };
    t.liquid.update(t.container);
    t.state.detached = true;
    game.emit("trypsin"); game.sfx("pour");
    game.toast("Medium off, 0.5 mL trypsin, 5 minutes at 37 °C: the cells round up and float. Fresh medium stops the trypsin: 1 mL of cell suspension.");
  } };
  if (slide && k === "coverslip") return { label: "Lower the cover slip: hold for slow, tap to drop", hold: 1.1, run: (slow) => {
    if (t.state.cover) return game.toast("This slide already has a cover slip.", "warn");
    t.state.cover = slow ? "angle" : "flat"; game.dispose(held, { quiet: true }); t.refresh(); game.sfx("glass");
    t.state.osmo.push({ t: performance.now(), osm: t.container.osmolarity() });
    t.state.mediumName = t.container.dominant();
    game.emit(`cover:${t.state.cover}`);
    game.toast(slow ? "Lowered slowly at an angle: no air bubbles" : "Dropped flat: that trapped air bubbles", slow ? "good" : "warn");
  } };
  if (slide && k === "pencil") return { label: "Label the frosted end", run: () => {
    t.state.label = t.state.sample ? { onion: "Onion", cheek: "Cheek", elodea: "Leaf", redonion: "Red onion", root: "Root tip" }[t.state.sample] : "Slide"; t.refresh(); game.emit("label"); } };
  if (slide && k === "towels") {
    if (t.state.sample === "root" && t.state.cover && !t.state.squashed) return { label: "Squash: press down hard through the folded towel (hold)", hold: 1.0, run: (full) => {
      if (!full) return game.toast("Press harder: hold right-click until the ring fills.", "warn");
      t.state.squashed = true; t.refresh(); game.emit("squash");
      game.toast(t.state.softened ? "Squashed flat: the cells spread out into a single layer." : "Squashed, but the unsoftened tip stays in a clump.", t.state.softened ? "good" : "warn");
    } };
    if (t.state.cover && t.state.edge) return { label: "Draw the drop under the cover slip with the towel", run: () => {
      // The towel soaks up the old liquid as the new drop runs in from the other side.
      const old = t.container.take(t.container.volume);
      const edge = t.state.edge;
      t.container = new Container(0.4, { [edge]: DROP * 0.85 });
      for (const kk in old) t.container.add(kk, (old[kk] / Math.max(1e-9, Object.values(old).reduce((a, b) => a + b, 0))) * DROP * 0.15);
      t.state.edge = null;
      t.state.osmo.push({ t: performance.now(), osm: t.container.osmolarity() });
      t.state.mediumName = edge;
      t.refresh(); game.emit(`irrigate:${edge}`);
      game.toast(`The towel draws the ${lowName(edge)} under the cover slip. Watch the cells.`);
    } };
    return { label: "Blot away the extra liquid", run: () => {
      const extra = t.container.volume - DROP;
      if (extra <= 0.001) return game.toast("There's no extra liquid to blot.");
      t.container.take(extra); t.refresh(); game.toast("Blotted. The cover slip sits flat now.");
    } };
  }
  if ((k === "stock" && t.kind === "dropper") || (k === "dropper" && t.kind === "stock")) {
    const stock = k === "stock" ? held : t, bottle = k === "stock" ? t : held;
    if (stock.solution !== bottle.solution) return { label: "Different solutions", run: () => game.toast("Only refill a bottle from the stock of the same solution.", "warn") };
    return { label: "Refill the dropper bottle from the stock", run: () => {
      const v = stock.container.pour(bottle.container, bottle.container.capacity - bottle.container.volume);
      stock.refresh(); bottle.refresh(); game.sfx("pour"); game.toast(`Refilled with ${fmtVol(v)}`);
    } };
  }
  if ((k === "stock" && t.kind === "beaker") || (k === "beaker" && t.kind === "stock")) {
    const from = k === "stock" ? held : t, to = k === "stock" ? t : held;
    return { label: `Pour 20 mL of ${lowName(from.solution)} into the beaker`, run: () => {
      const v = from.container.pour(to.container, 20); from.refresh(); to.refresh(); game.sfx("pour");
      game.toast(v > 0 ? `Poured ${fmtVol(v)}. The beaker now holds ${to.container.describe()}` : "The beaker is full.");
    } };
  }
  if (k === "dropper" && t.kind === "beaker") return { label: `Add a drop to the beaker`, run: () => { held.container.pour(t.container, DROP); held.refresh(); t.refresh(); game.sfx("drip"); } };
  if (k === "plate" && t.kind === "platereader") return { label: "Slide the plate into the reader", run: () => game.readPlate(held, t) };
  if (k === "watchglass" && t.kind === "waterbath") return { label: "Stand the watch glass on the shelf in the warm water", run: () => {
    game.place(held, t.obj.position.clone().setY(t.obj.position.y + t.obj.userData.shelfY));
  } };
  return null;
}

function tileTop(tile) { return tile.obj.position.clone().add(new THREE.Vector3(0, 0.006, 0)); }

function useOnFurniture(game, held, f) {
  const kind = f.userData.target.kind, k = held.kind;
  if (kind === "microscope" && ["slide", "dish", "plate", "hemocytometer"].includes(k) && !game.stageItem) return { label: `Clip the ${lower(held)} onto the stage`, run: () => game.putOnStage(held) };
  if (kind === "incubator" && (k === "dish" || k === "plate")) return { label: `Put ${held.name} in the incubator`, run: () => game.putInIncubator(held) };
  if (kind === "glassbin" && (k === "slide" || k === "coverslip")) return { label: `Throw the ${lower(held)} into the glass waste`, run: () => { game.dispose(held); game.emit(`dispose:${k}`); } };
  if (kind === "glassbin" && k === "pipette" && held.state.tip) return { label: "Tips go in the tip waste", run: () => game.toast("Plastic tips go in the tip waste on the microscope bench, not the glass bin.", "warn") };
  if (kind === "sink" && k === "slide") return { label: "Rinse the slide clean", run: () => {
    held.state = { sample: null, cover: null, label: "", osmo: [] }; held.container = new Container(0.4); held.refresh(); game.runTap(); game.toast("Rinsed and dried. A clean slide again.");
  } };
  if (kind === "sink" && k === "beaker") return { label: "Fill the beaker with tap water", run: () => { held.container.add("tap", 200); held.refresh(); game.runTap(); } };
  return null;
}

function furnitureAction(game, f) {
  const t = f.userData.target;
  switch (t.kind) {
    case "microscope": return game.stageItem
      ? { label: "Look through the eyepiece", run: () => game.enterScope() }
      : { label: "Look through the eyepiece", run: () => game.toast("Clip a slide, a dish, a plate or the hemocytometer onto the stage first.", "warn") };
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
    case "onion": case "redonion": return { label: `${t.name}: cut it with the scalpel`, run: () => game.toast(t.state.cut ? "Peel the thin skin with the forceps." : "Hold the scalpel, then right-click the onion.") };
    case "rootonion": return { label: "Roots: cut a tip with the scalpel", run: () => game.toast("Hold the scalpel and right-click the onion to cut a root tip, then pick it up with the forceps.") };
    case "waterbath": return { label: "Water bath: 60 °C", run: () => game.toast("Warm acid softens a root tip in 5 minutes. The watch glass sits on the shelf inside.") };
    case "watchglass": return { label: t.state.root ? "Root tip softening in warm acid" : "Watch glass", run: () => {
      if (!t.state.root) return game.toast("Add 1 M hydrochloric acid, then put a root tip in it with the forceps.");
      const left = SOFTEN_S - (performance.now() - t.state.root) / 1000;
      game.toast(left > 0 ? `Softening: ${Math.ceil(left)} s more.` : "Soft now. Take it out with the forceps.");
    } };
    case "slide": return { label: "Slide: add a sample, a drop, a cover slip", run: () => game.toast("Hold forceps, a swab, a dropper bottle or a cover slip, then right-click the slide.") };
    case "dish": return { label: `${t.name}: ${t.state.t} h in culture`, run: () => game.toast("Clip the dish onto the microscope stage to look at the cells.") };
    case "plate": return { label: `${t.name}: ${t.state.t} h in culture`, run: () => game.toast("Hold the micropipette with a tip and right-click a well to add a drug to it.") };
    case "tipbox": return { label: `Pipette tips: ${t.state.left} left`, run: () => game.toast("Hold the micropipette, then right-click the box to fit a tip.") };
    case "tipwaste": return { label: "Tip waste", run: () => game.toast("Hold the micropipette with a used tip and right-click here to eject it.") };
    case "tuberack": return { label: "Tube rack", run: () => game.toast("Hold a microtube and right-click the rack to stand it in a hole.") };
    case "platereader": return { label: "Plate reader: hold a plate to read it", run: () => game.toast("Bring the 6-well plate from the incubator, then right-click the reader.") };
    case "hemocytometer": return { label: "Hemocytometer", run: () => game.toast("Load 10 µL of cells mixed with trypan blue with the micropipette, then clip it onto the microscope.") };
  }
  return null;
}
