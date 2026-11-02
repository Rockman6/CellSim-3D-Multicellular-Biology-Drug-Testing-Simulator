// Find anything in the lab by name, typos and all. Every item and every piece
// of furniture is scored against what you typed: whole words, the start of a
// word, and spelling mistakes (one wrong, missing, extra or swapped letter in
// every four still matches). Each result says where the thing is now.
import * as THREE from "three";

// Other words people use for things, by item or furniture kind.
const ALIASES = {
  goggles: ["safety glasses", "glasses", "eye protection"],
  sink: ["tap", "faucet", "wash hands", "water"],
  microscope: ["scope", "eyepiece", "magnify"],
  incubator: ["cells", "culture", "warm", "37"],
  glassbin: ["bin", "trash", "waste", "rubbish", "garbage", "sharps"],
  dropper: ["stain", "bottle", "drops"],
  stock: ["refill", "reagent", "big bottle"],
  drug: ["medicine", "chemo", "cytotoxic", "treatment"],
  dish: ["cells", "culture", "petri dish", "cancer cells", "a549"],
  plate: ["6 well plate", "well plate", "wells", "multiwell"],
  pipette: ["pipette", "pipettor", "p200"],
  tipbox: ["pipette tips", "tips"],
  tipwaste: ["tip bin", "waste", "used tips"],
  platereader: ["reader", "plate reader", "viability", "results"],
  hemocytometer: ["haemocytometer", "counting chamber", "neubauer", "cell counter"],
  tube: ["microtube", "eppendorf", "tube"],
  trypsin: ["detach", "trypsin edta"],
  tank: ["elodea", "pondweed", "plant", "leaf", "aquarium"],
  onion: ["bulb", "vegetable"],
  redonion: ["red onion", "purple onion", "osmosis"],
  rootonion: ["roots", "root tip", "mitosis", "growing onion"],
  waterbath: ["water bath", "heater", "60"],
  watchglass: ["watch glass", "acid dish"],
  slidebox: ["slides", "glass slide", "microscope slide"],
  slide: ["glass slide", "microscope slide"],
  coverbox: ["coverslip", "cover slips", "cover glass"],
  coverslip: ["coverslip", "cover glass"],
  towels: ["tissue", "paper", "blot", "kitchen roll"],
  swabcup: ["cotton bud", "q tip", "swab", "cotton swab"],
  swab: ["cotton bud", "q tip"],
  forceps: ["tweezers", "pincers", "tongs"],
  scalpel: ["knife", "blade", "cutter"],
  lens: ["magnifying glass", "magnifier"],
  pencilcup: ["pen", "pencil", "marker"],
  pencil: ["pen", "marker"],
  sanitizer: ["hand gel", "alcohol gel", "clean hands"],
  journal: ["notebook", "notes", "results", "data"],
  clipboard: ["method", "instructions", "worksheet"],
  beaker: ["glass", "container", "cup"],
  tile: ["white tile", "work surface"],
  drawer: ["spares", "storage", "cabinet"],
  stool: ["chair", "seat"],
  door: ["exit", "leave", "way out"],
  switch: ["lights", "lamp", "light switch"],
  clock: ["time"],
  whiteboard: ["board", "method", "instructions"],
  blinds: ["window", "curtain", "shade", "sunlight"],
};

export function norm(s) {
  return String(s).toLowerCase().normalize("NFKD").replace(/[̀-ͯ]/g, "").replace(/[^a-z0-9µ ]+/g, " ").replace(/\s+/g, " ").trim();
}

// Optimal string alignment distance: Levenshtein plus swapped neighbours.
function osa(a, b) {
  const m = a.length, n = b.length;
  if (!m) return n;
  if (!n) return m;
  let p2 = null, p1 = Array.from({ length: n + 1 }, (_, j) => j), cur;
  for (let i = 1; i <= m; i++) {
    cur = [i];
    for (let j = 1; j <= n; j++) {
      const cost = a[i - 1] === b[j - 1] ? 0 : 1;
      let v = Math.min(p1[j] + 1, cur[j - 1] + 1, p1[j - 1] + cost);
      if (p2 && i > 1 && j > 1 && a[i - 1] === b[j - 2] && a[i - 2] === b[j - 1]) v = Math.min(v, p2[j - 2] + 1);
      cur[j] = v;
    }
    p2 = p1; p1 = cur;
  }
  return p1[n];
}

// How well one typed word matches one word of a name: 1 is exact.
function wordScore(q, w) {
  if (q === w) return 1;
  if (w.startsWith(q)) return q.length >= 3 ? 0.96 : 0.7 + 0.08 * q.length;
  const full = 1 - osa(q, w) / Math.max(q.length, w.length);
  let pre = 0;
  if (q.length >= 3) for (const len of [q.length - 1, q.length, q.length + 1]) {
    if (len < 2 || len > w.length) continue;
    pre = Math.max(pre, 1 - osa(q, w.slice(0, len)) / Math.max(q.length, len));
  }
  return Math.max(full, pre * 0.93);
}

export function similarity(query, terms) {
  const q = norm(query);
  if (!q) return 0;
  const qWords = q.split(" "), qJoined = q.replace(/ /g, "");
  let best = 0;
  for (const t of terms) {
    if (t === q) return 1;
    let s = t.includes(q) ? 0.94 : 0;
    const tJoined = t.replace(/ /g, "");
    s = Math.max(s, 1 - osa(qJoined, tJoined) / Math.max(qJoined.length, tJoined.length));
    const words = t.split(" ");
    let sum = 0;
    for (const a of qWords) { let m = 0; for (const w of words) m = Math.max(m, wordScore(a, w)); sum += m; }
    s = Math.max(s, (sum / qWords.length) * 0.97);
    best = Math.max(best, s);
  }
  return best;
}

/* ---------- where things are ---------- */
const PLACES = {
  drawer: "against the back wall: spare slides, cover slips and stains", sink: "set into the big bench",
  door: "in the corner by the whiteboard", switch: "on the wall by the door", clock: "on the wall above the big bench",
  whiteboard: "on the wall by the teacher's desk", blinds: "on the windows", glassbin: "on the floor by the bench arm", stool: "by the benches",
};

function zone(game, p) {
  const room = game.room;
  if (p.y < 0.15) return "on the floor";
  if (Math.abs(p.y - room.TABLE_Y) < 0.12 && (p.z < -1.25 || p.x < -3.95)) {
    if (p.x > -1.1 && p.z < -1.6) return "on the big bench, by the sink";
    if (p.x < -3.95) return "on the big bench, by the windows";
    return "on the big bench";
  }
  const B = room.BENCH;
  if (Math.abs(p.x - B.x) < B.w / 2 + 0.1 && Math.abs(p.z - B.z) < B.d / 2 + 0.1) return "on the microscope bench";
  if (p.x > 3.2 && p.z < -1.2) return "on the teacher's desk";
  if (p.x < -2.9 && p.z > 3.2) return "on the steel shelves";
  return "in the lab";
}

function entry(game, { name, kind, object, item, where }) {
  const extra = ALIASES[kind] || [];
  const sol = item && item.solution ? [item.solution] : [];
  const terms = [...new Set([name, ...extra, ...sol].map(norm).filter(Boolean))];
  return { name, kind, object, item, where, terms };
}

// Everything you could be looking for, with where it is right now.
export function buildIndex(game) {
  const out = [];
  const room = game.room;
  for (const it of game.items) {
    const slot = game.inv.slots.indexOf(it);
    if (slot >= 0) { out.push(entry(game, { name: it.name, kind: it.kind, item: it, object: null, where: slot < 9 ? `in your hotbar, slot ${slot + 1}` : "in your inventory (E)" })); continue; }
    if (it === game.stageItem) { out.push(entry(game, { name: it.name, kind: it.kind, item: it, object: room.scope, where: "on the microscope stage" })); continue; }
    if (game.incubator.dishes.includes(it)) { out.push(entry(game, { name: it.name, kind: it.kind, item: it, object: room.incubator, where: "inside the incubator" })); continue; }
    if (!it.obj.parent) continue;   // used up or thrown away
    out.push(entry(game, { name: it.name, kind: it.kind, item: it, object: it.obj, where: zone(game, it.obj.getWorldPosition(new THREE.Vector3())) }));
  }
  const seen = new Map();
  for (const t of room.targets) {
    const d = t.userData.target;
    if (!d) continue;
    const key = d.kind === "drawer" ? "drawer" : d.name;
    const p = t.getWorldPosition(new THREE.Vector3());
    const dist = p.distanceTo(game.player.pos);
    // Four stools or three blinds: keep the nearest of each.
    if (seen.has(key) && seen.get(key).dist <= dist) continue;
    const e = entry(game, { name: d.kind === "drawer" ? "Drawer cabinet" : d.name, kind: d.kind, object: d.kind === "drawer" ? room.cabinet : t, where: PLACES[d.kind] || zone(game, p) });
    e.dist = dist;
    seen.set(key, e);
  }
  out.push(...seen.values());
  return out;
}

export function search(game, query, max = 7) {
  const list = buildIndex(game);
  const pos = game.player.pos;
  return list
    .map((e) => ({ ...e, score: similarity(query, e.terms), dist: e.object ? e.object.getWorldPosition(new THREE.Vector3()).distanceTo(pos) : 0 }))
    .filter((e) => e.score >= (norm(query).length <= 3 ? 0.7 : 0.62))
    .sort((a, b) => b.score - a.score || a.dist - b.dist)
    .slice(0, max);
}

/* ---------- following a marker to it ---------- */
const $ = (id) => document.getElementById(id);

export class Tracker {
  constructor(game) {
    this.game = game;
    this.target = null;
    this.box = new THREE.Box3();
    this.v = new THREE.Vector3();
  }
  start(e) {
    this.target = e;
    this.game.outliner.set("track", e.object);
    $("trackName").textContent = e.name;
    $("trackMarker").hidden = false;
    this.game.dirty = true;
  }
  stop() {
    if (!this.target) return;
    this.target = null;
    this.game.outliner.set("track", null);
    $("trackMarker").hidden = true;
    this.game.dirty = true;
  }
  // Is the thing you aim at the thing you were looking for?
  reached(aim) {
    const t = this.target;
    if (!t || !aim) return false;
    if (t.item && aim.item === t.item) return true;
    if (aim.furniture && aim.furniture === t.object) return true;
    if (aim.furniture && t.object && aim.furniture.parent === t.object) return true;
    return false;
  }
  update(camera) {
    const t = this.target;
    if (!t) return;
    const c = this.box.setFromObject(t.object).getCenter(this.v);
    const dist = c.distanceTo(camera.position);
    $("trackDist").textContent = `${dist < 10 ? dist.toFixed(1) : dist.toFixed(0)} m`;
    const p = c.clone().project(camera);
    const behind = c.clone().applyMatrix4(camera.matrixWorldInverse).z > 0;
    const W = innerWidth, H = innerHeight, m = 48;
    let x = (p.x * 0.5 + 0.5) * W, y = (-p.y * 0.5 + 0.5) * H;
    const on = !behind && x > m && x < W - m && y > m && y < H - m;
    const el = $("trackMarker");
    el.classList.toggle("edge", !on);
    if (!on) {
      // Pin it to the screen edge, pointing the way to turn.
      let dx = x - W / 2, dy = y - H / 2;
      if (behind) { dx = -dx; dy = -dy; if (Math.abs(dx) < 1 && Math.abs(dy) < 1) dy = H; }
      const k = Math.min((W / 2 - m) / Math.max(1e-3, Math.abs(dx)), (H / 2 - m) / Math.max(1e-3, Math.abs(dy)));
      x = W / 2 + dx * k; y = H / 2 + dy * k;
      el.style.setProperty("--a", `${Math.atan2(dy, dx)}rad`);
    }
    el.style.transform = `translate(${x.toFixed(0)}px, ${y.toFixed(0)}px)`;
  }
}
