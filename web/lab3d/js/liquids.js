// Liquids: every container holds real volumes of named solutions. Colour
// mixes by volume weighted by how strongly each solution colours the mix,
// and drug concentrations dilute by volume, which is what the engine gets.
import * as THREE from "three";

export const SOLUTIONS = {
  water:        { name: "Distilled water", color: "#d7ebf2", strength: 0.02, alpha: 0.22 },
  tap:          { name: "Tap water", color: "#d2e6ee", strength: 0.02, alpha: 0.22 },
  saline:       { name: "Saline 0.9%", color: "#dcecf2", strength: 0.02, alpha: 0.22 },
  iodine:       { name: "Iodine solution", color: "#8a4a12", strength: 1, alpha: 0.9, stain: "iodine" },
  mb:           { name: "Methylene blue", color: "#1b48b0", strength: 1.3, alpha: 0.92, stain: "mb" },
  ethanol:      { name: "Ethanol 70%", color: "#e9f1f3", strength: 0.01, alpha: 0.18 },
  disinfectant: { name: "Disinfectant 1%", color: "#bfe3a6", strength: 0.3, alpha: 0.5 },
  medium:       { name: "Culture medium (phenol red)", color: "#dc4a63", strength: 0.6, alpha: 0.7 },
  // Drug stocks. µM is the stock concentration; the dish dilutes it.
  paclitaxel:   { name: "Paclitaxel 10 µM", color: "#f2f1ea", strength: 0.01, alpha: 0.25, drug: "paclitaxel", uM: 10 },
  doxorubicin:  { name: "Doxorubicin 20 µM", color: "#cf2a1e", strength: 0.9, alpha: 0.8, drug: "doxorubicin", uM: 20 },
  cisplatin:    { name: "Cisplatin 200 µM", color: "#f4f1d6", strength: 0.02, alpha: 0.25, drug: "cisplatin", uM: 200 },
};

export class Container {
  constructor(capacity, contents = {}) {
    this.capacity = capacity; // mL
    this.c = { ...contents };
  }
  get volume() { let v = 0; for (const k in this.c) v += this.c[k]; return v; }
  get fraction() { return Math.min(1, this.volume / this.capacity); }
  get empty() { return this.volume < 1e-6; }
  add(id, mL) {
    const room = this.capacity - this.volume;
    const v = Math.max(0, Math.min(mL, room));
    this.c[id] = (this.c[id] || 0) + v;
    return v;
  }
  // Remove `mL` of the mixture and return what was removed, in proportion.
  take(mL) {
    const total = this.volume;
    if (total <= 0) return {};
    const f = Math.min(1, mL / total), out = {};
    for (const k in this.c) {
      out[k] = this.c[k] * f;
      this.c[k] -= out[k];
      if (this.c[k] < 1e-9) delete this.c[k];
    }
    return out;
  }
  pour(other, mL) {
    const room = other.capacity - other.volume;
    const got = this.take(Math.min(mL, room));
    let moved = 0;
    for (const k in got) moved += other.add(k, got[k]);
    return moved;
  }
  dominant() {
    let best = null, w = 0;
    for (const k in this.c) { const s = this.c[k] * (SOLUTIONS[k]?.strength ?? 0.05); if (s > w) { w = s; best = k; } }
    return best;
  }
  // The stain that will show under the microscope, if any.
  stain() {
    let best = null, w = 0;
    for (const k in this.c) if (SOLUTIONS[k]?.stain && this.c[k] > w) { w = this.c[k]; best = SOLUTIONS[k].stain; }
    return best;
  }
  // µM of each drug in this container.
  drugs() {
    const v = this.volume, out = {};
    if (v <= 0) return out;
    for (const k in this.c) { const s = SOLUTIONS[k]; if (s?.drug) out[s.drug] = (out[s.drug] || 0) + (s.uM * this.c[k]) / v; }
    return out;
  }
  color(target = new THREE.Color()) {
    let wsum = 0, r = 0, g = 0, b = 0, a = 0;
    const tmp = new THREE.Color();
    for (const k in this.c) {
      const s = SOLUTIONS[k] || SOLUTIONS.water;
      const w = this.c[k] * s.strength;
      tmp.set(s.color);
      r += tmp.r * w; g += tmp.g * w; b += tmp.b * w; a += s.alpha * w; wsum += w;
    }
    if (wsum <= 0) return { color: target.set(SOLUTIONS.water.color), alpha: 0.2 };
    target.setRGB(r / wsum, g / wsum, b / wsum);
    return { color: target, alpha: a / wsum };
  }
  describe() {
    const v = this.volume;
    if (v <= 0) return "Empty";
    const parts = Object.entries(this.c).sort((a, b) => b[1] - a[1]).map(([k, mL]) => `${SOLUTIONS[k]?.name || k}`);
    return `${fmtVol(v)} · ${parts.slice(0, 2).join(" + ")}`;
  }
}

export function fmtVol(mL) {
  if (mL >= 10) return `${mL.toFixed(0)} mL`;
  if (mL >= 1) return `${mL.toFixed(1)} mL`;
  return `${(mL * 1000).toFixed(0)} µL`;
}

/* A liquid column inside a cylindrical vessel: its height follows the fill
   level and its colour follows the mix. */
export function liquidColumn({ radius, bottom, depth, radiusTop = radius, segments = 28 }) {
  const mat = new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.08, transparent: true, opacity: 0.8, depthWrite: false });
  const geo = new THREE.CylinderGeometry(1, 1, 1, segments, 1, false);
  geo.translate(0, 0.5, 0);
  const mesh = new THREE.Mesh(geo, mat);
  mesh.position.y = bottom;
  mesh.renderOrder = 3;
  return {
    mesh,
    update(container) {
      const f = container.fraction;
      mesh.visible = f > 0.002;
      const h = Math.max(0.0005, depth * f);
      const r = radius + (radiusTop - radius) * f;
      mesh.scale.set((radius + r) / 2, h, (radius + r) / 2);
      const { color, alpha } = container.color();
      mat.color.copy(color);
      mat.opacity = Math.min(0.95, 0.25 + alpha);
    },
  };
}
