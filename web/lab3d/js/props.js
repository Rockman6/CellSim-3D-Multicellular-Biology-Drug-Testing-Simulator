// Small lab equipment built in code at real size (metres), so each piece can
// hold a liquid that drains and changes colour. Every prop's origin is the
// middle of its base, so it can be set down on any surface.
import * as THREE from "three";
import { Container, SOLUTIONS, liquidColumn } from "./liquids.js";
import { link } from "./assets.js";

const V2 = THREE.Vector2;
export const MAT = {
  glass: () => new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.04, metalness: 0, transparent: true, opacity: 0.22, ior: 1.5, clearcoat: 1, clearcoatRoughness: 0.05, depthWrite: false, side: THREE.DoubleSide }),
  amber: () => new THREE.MeshPhysicalMaterial({ color: 0x6b3510, roughness: 0.08, transparent: true, opacity: 0.72, clearcoat: 1, depthWrite: false }),
  blueGlass: () => new THREE.MeshPhysicalMaterial({ color: 0x1d3f8c, roughness: 0.08, transparent: true, opacity: 0.7, clearcoat: 1, depthWrite: false }),
  steel: () => new THREE.MeshStandardMaterial({ color: 0xc9ced3, metalness: 1, roughness: 0.28 }),
  rubber: () => new THREE.MeshStandardMaterial({ color: 0x202326, roughness: 0.75 }),
  plastic: (c) => new THREE.MeshStandardMaterial({ color: c, roughness: 0.45 }),
  paper: (c = 0xf6f4ee) => new THREE.MeshStandardMaterial({ color: c, roughness: 0.9 }),
};

function mesh(geo, mat, cast = true) {
  const m = new THREE.Mesh(geo, mat);
  m.castShadow = cast;
  m.receiveShadow = true;
  return m;
}

// A label printed on paper, as a canvas texture.
export function labelTexture(lines, { w = 256, h = 128, bg = "#f7f5ef", fg = "#1d2326", band = null, size = 30 } = {}) {
  const c = document.createElement("canvas");
  c.width = w; c.height = h;
  const g = c.getContext("2d");
  g.fillStyle = bg; g.fillRect(0, 0, w, h);
  if (band) { g.fillStyle = band; g.fillRect(0, 0, w, h * 0.22); }
  g.fillStyle = fg;
  g.textAlign = "center"; g.textBaseline = "middle";
  const n = lines.length;
  lines.forEach((t, i) => {
    const s = i === 0 ? size : size * 0.62;
    g.font = `${i === 0 ? 700 : 500} ${s}px "IBM Plex Sans", Helvetica, sans-serif`;
    g.fillText(t, w / 2, h * (band ? 0.6 : 0.5) + (i - (n - 1) / 2) * size * 0.95, w - 12);
  });
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 4;
  return t;
}

function lathe(points, segments = 32) {
  return new THREE.LatheGeometry(points.map(([r, y]) => new V2(r, y)), segments);
}

/* ---------- slide, cover slip and their boxes ---------- */
export function makeSlide() {
  const g = new THREE.Group();
  const L = 0.076, W = 0.026, T = 0.001;
  const glass = mesh(new THREE.BoxGeometry(L, T, W), MAT.glass());
  glass.position.y = T / 2;
  const frost = mesh(new THREE.BoxGeometry(0.02, T * 1.05, W), MAT.paper(0xeef0ef));
  frost.position.set(-L / 2 + 0.01, T / 2, 0);
  g.add(glass, frost);
  // The sample, the drop of liquid and the cover slip, shown as they are added.
  const sample = mesh(new THREE.CircleGeometry(0.0065, 24), new THREE.MeshStandardMaterial({ color: 0xe9dfc8, roughness: 0.6, transparent: true, opacity: 0.75 }), false);
  sample.rotation.x = -Math.PI / 2; sample.position.set(0.008, T + 0.0002, 0); sample.visible = false;
  const drop = mesh(new THREE.SphereGeometry(1, 20, 10, 0, Math.PI * 2, 0, Math.PI / 2), new THREE.MeshPhysicalMaterial({ color: 0xd7ebf2, roughness: 0.05, transparent: true, opacity: 0.5, clearcoat: 1, depthWrite: false }), false);
  drop.position.set(0.008, T, 0); drop.visible = false; drop.renderOrder = 3; drop.scale.setScalar(0.003);
  const cover = mesh(new THREE.BoxGeometry(0.022, 0.0003, 0.022), MAT.glass(), false);
  cover.position.set(0.008, T + 0.0007, 0); cover.visible = false;
  const label = mesh(new THREE.PlaneGeometry(0.018, 0.022), new THREE.MeshStandardMaterial({ transparent: true, roughness: 0.9 }), false);
  label.rotation.x = -Math.PI / 2; label.position.set(-L / 2 + 0.01, T * 1.1, 0); label.visible = false;
  g.add(sample, drop, cover, label);
  link(g, "parts", { sample, drop, cover, label });
  return g;
}

export function updateSlideLook(slide) {
  const { sample, drop, cover, label } = slide.obj.userData.parts;
  const s = slide.state;
  sample.visible = !!s.sample;
  if (s.sample) sample.material.color.set(s.sample === "elodea" ? 0x4a9a3e : s.sample === "cheek" ? 0xf0dcd4 : 0xece3cf);
  const v = slide.container.volume;
  drop.visible = v > 0.001;
  if (drop.visible) {
    const { color, alpha } = slide.container.color();
    drop.material.color.copy(color);
    drop.material.opacity = Math.min(0.9, 0.3 + alpha * 0.6);
    const r = 0.0035 + Math.sqrt(v) * 0.012;               // a 50 µL drop is about 6 mm across
    const flat = s.cover ? 0.04 : 1;
    drop.scale.set(s.cover ? 0.0105 : r, s.cover ? 0.0003 : r * 0.35 * flat, s.cover ? 0.0105 : r);
  }
  cover.visible = !!s.cover;
  cover.rotation.z = 0;
  if (s.label && !label.material.map) {
    label.material.map = labelTexture([s.label], { w: 128, h: 160, bg: "#eef0ef", size: 26 });
    label.material.needsUpdate = true;
  }
  label.visible = !!s.label;
}

export function makeBox({ w, h, d, color, labelLines, band }) {
  const g = new THREE.Group();
  const body = mesh(new THREE.BoxGeometry(w, h, d), MAT.plastic(color));
  body.position.y = h / 2;
  g.add(body);
  if (labelLines) {
    const lab = mesh(new THREE.PlaneGeometry(w * 0.8, h * 0.7), new THREE.MeshStandardMaterial({ map: labelTexture(labelLines, { band }), roughness: 0.85 }), false);
    lab.position.set(0, h / 2, d / 2 + 0.0006);
    g.add(lab);
  }
  return g;
}

export function makeCoverSlip() {
  const g = new THREE.Group();
  const s = mesh(new THREE.BoxGeometry(0.022, 0.0004, 0.022), MAT.glass(), false);
  s.position.y = 0.0002;
  g.add(s);
  return g;
}

/* ---------- bottles and beakers ---------- */
// A 30 mL dropper bottle with a rubber bulb; its liquid follows the volume.
export function makeDropperBottle(solution, volume = 25) {
  const s = SOLUTIONS[solution];
  const g = new THREE.Group();
  const glassMat = solution === "iodine" ? MAT.amber() : solution === "mb" ? MAT.blueGlass() : MAT.glass();
  const r = 0.0175, h = 0.062;
  g.add(mesh(lathe([[0.0001, 0], [r, 0], [r, h], [r * 0.55, h + 0.009], [0.0075, h + 0.012], [0.0075, h + 0.016]]), glassMat));
  const cap = mesh(new THREE.CylinderGeometry(0.0095, 0.0095, 0.012, 24), MAT.rubber());
  cap.position.y = h + 0.022;
  const bulb = mesh(new THREE.SphereGeometry(0.0085, 20, 14), MAT.rubber());
  bulb.scale.y = 1.35; bulb.position.y = h + 0.037;
  g.add(cap, bulb);
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0004, r + 0.0004, 0.03, 32, 1, true, -0.9, 1.8), new THREE.MeshStandardMaterial({ map: labelTexture([shortName(solution)], { size: 34, band: s.color }), roughness: 0.85 }), false);
  lab.position.y = 0.03;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0015, bottom: 0.002, depth: h - 0.004 });
  g.add(liq.mesh);
  const container = new Container(30, { [solution]: volume });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

// A 500 mL screw-cap reagent bottle with a printed label.
export function makeStockBottle(solution, volume = 450) {
  const s = SOLUTIONS[solution];
  const g = new THREE.Group();
  const r = 0.042, h = 0.15;
  const glassMat = solution === "iodine" ? MAT.amber() : MAT.glass();
  g.add(mesh(lathe([[0.0001, 0], [r, 0.002], [r, h], [0.03, h + 0.018], [0.016, h + 0.025], [0.016, h + 0.035]], 40), glassMat));
  const cap = mesh(new THREE.CylinderGeometry(0.019, 0.019, 0.022, 28), MAT.plastic(0x2a5fb0));
  cap.position.y = h + 0.046;
  g.add(cap);
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0006, r + 0.0006, 0.07, 40, 1, true, -1.1, 2.2), new THREE.MeshStandardMaterial({ map: labelTexture([shortName(solution), "Stock · 500 mL"], { w: 320, h: 160, band: s.color, size: 36 }), roughness: 0.85 }), false);
  lab.position.y = 0.065;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.002, bottom: 0.003, depth: h - 0.006 });
  g.add(liq.mesh);
  const container = new Container(500, { [solution]: volume });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

// A 10 mL amber vial of a cytotoxic drug stock, with a dropper cap.
export function makeDrugVial(drug, volume = 8) {
  const s = SOLUTIONS[drug];
  const g = new THREE.Group();
  const r = 0.012, h = 0.045;
  g.add(mesh(lathe([[0.0001, 0], [r, 0], [r, h], [0.007, h + 0.006], [0.007, h + 0.01]]), MAT.amber()));
  const cap = mesh(new THREE.CylinderGeometry(0.0085, 0.0085, 0.01, 20), MAT.plastic(0xd8d8d2));
  cap.position.y = h + 0.015;
  const tip = mesh(new THREE.ConeGeometry(0.004, 0.012, 16), MAT.plastic(0xd8d8d2));
  tip.position.y = h + 0.026;
  g.add(cap, tip);
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0004, r + 0.0004, 0.028, 32, 1, true, -1, 2), new THREE.MeshStandardMaterial({ map: labelTexture([shortName(drug), `${s.uM} µM · CYTOTOXIC`], { w: 256, h: 128, band: "#d8342a", size: 30 }), roughness: 0.85 }), false);
  lab.position.y = 0.024;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0012, bottom: 0.0015, depth: h - 0.003 });
  g.add(liq.mesh);
  const container = new Container(10, { [drug]: volume });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

// A 250 mL low-form beaker with a pouring lip and printed scale.
export function makeBeaker(contents = {}, capacity = 250) {
  const g = new THREE.Group();
  const r = 0.035, h = 0.095;
  g.add(mesh(lathe([[0.0001, 0], [r, 0], [r, h], [r + 0.003, h + 0.003]], 40), MAT.glass()));
  const grad = document.createElement("canvas");
  grad.width = 64; grad.height = 256;
  const gg = grad.getContext("2d");
  gg.fillStyle = "rgba(255,255,255,0.9)";
  for (let i = 1; i <= 5; i++) { const y = 256 - i * 44; gg.fillRect(0, y, i % 2 ? 22 : 34, 3); gg.font = "18px sans-serif"; if (i % 2 === 0) gg.fillText(`${i * 50}`, 30, y + 6); }
  const gt = new THREE.CanvasTexture(grad);
  const marks = mesh(new THREE.CylinderGeometry(r + 0.0003, r + 0.0003, h, 24, 1, true, -0.35, 0.7), new THREE.MeshBasicMaterial({ map: gt, transparent: true, depthWrite: false }), false);
  marks.position.y = h / 2;
  g.add(marks);
  const liq = liquidColumn({ radius: r - 0.0015, bottom: 0.0015, depth: h - 0.004 });
  g.add(liq.mesh);
  const container = new Container(capacity, contents);
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

/* ---------- culture dish ---------- */
export function makeCultureDish(name) {
  const g = new THREE.Group();
  const r = 0.03, h = 0.013;
  const base = mesh(lathe([[0.0001, 0], [r, 0], [r, h * 0.8]], 48), MAT.glass(), false);
  const lid = mesh(lathe([[0.0001, h + 0.001], [r + 0.0015, h + 0.001], [r + 0.0015, h * 0.45]], 48), MAT.glass(), false);
  g.add(base, lid);
  const lab = mesh(new THREE.PlaneGeometry(0.03, 0.012), new THREE.MeshStandardMaterial({ map: labelTexture([name, "A549 lung cells"], { w: 256, h: 100, size: 30 }), roughness: 0.9, transparent: true }), false);
  lab.rotation.x = -Math.PI / 2; lab.position.y = h + 0.0013;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0012, bottom: 0.0008, depth: 0.004 });
  g.add(liq.mesh);
  const container = new Container(5.2, { medium: 5 });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

/* ---------- pondweed tank ---------- */
export function makeTank() {
  const g = new THREE.Group();
  const w = 0.36, h = 0.24, d = 0.22;
  const glass = mesh(new THREE.BoxGeometry(w, h, d), MAT.glass(), false);
  glass.position.y = h / 2;
  const water = mesh(new THREE.BoxGeometry(w - 0.006, h * 0.82, d - 0.006), new THREE.MeshPhysicalMaterial({ color: 0x9fc9c4, transparent: true, opacity: 0.3, roughness: 0.05, depthWrite: false }), false);
  water.position.y = h * 0.41 + 0.003;
  water.renderOrder = 2;
  const gravel = mesh(new THREE.BoxGeometry(w - 0.008, 0.02, d - 0.008), new THREE.MeshStandardMaterial({ color: 0x8b8173, roughness: 1 }));
  gravel.position.y = 0.012;
  g.add(gravel, glass, water);
  // Elodea: stems with whorls of three small leaves.
  const stemMat = new THREE.MeshStandardMaterial({ color: 0x3f7f2e, roughness: 0.7 });
  const leafMat = new THREE.MeshStandardMaterial({ color: 0x4f9a35, roughness: 0.6, side: THREE.DoubleSide });
  const leafGeo = new THREE.PlaneGeometry(0.004, 0.012);
  leafGeo.translate(0, 0.006, 0);
  let seed = 3;
  const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  const leaves = new THREE.InstancedMesh(leafGeo, leafMat, 9 * 30 * 3);
  let li = 0;
  const m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), e = new THREE.Euler();
  for (let s = 0; s < 9; s++) {
    const x0 = (rnd() - 0.5) * (w - 0.06), z0 = (rnd() - 0.5) * (d - 0.06);
    const pts = [];
    for (let i = 0; i <= 6; i++) pts.push(new THREE.Vector3(x0 + Math.sin(i * 0.9 + s) * 0.012, 0.02 + i * 0.03, z0 + Math.cos(i * 0.7 + s) * 0.01));
    const curve = new THREE.CatmullRomCurve3(pts);
    g.add(mesh(new THREE.TubeGeometry(curve, 24, 0.0012, 5), stemMat, false));
    for (let k = 0; k < 30; k++) {
      const p = curve.getPoint(k / 30);
      for (let j = 0; j < 3; j++) {
        e.set(0.9 + rnd() * 0.3, (j / 3) * Math.PI * 2 + k * 0.7, 0);
        q.setFromEuler(e);
        m4.compose(p, q, new THREE.Vector3(1, 1, 1));
        leaves.setMatrixAt(li++, m4);
      }
    }
  }
  leaves.count = li;
  g.add(leaves);
  const lab = mesh(new THREE.PlaneGeometry(0.12, 0.05), new THREE.MeshStandardMaterial({ map: labelTexture(["Pondweed", "Elodea"], { w: 256, h: 110, size: 34 }), roughness: 0.9 }), false);
  lab.position.set(0, h - 0.04, d / 2 + 0.001);
  g.add(lab);
  return g;
}

/* ---------- tools ---------- */
export function makeTile() {
  const g = new THREE.Group();
  const t = mesh(new THREE.BoxGeometry(0.15, 0.006, 0.15), new THREE.MeshPhysicalMaterial({ color: 0xf7f7f4, roughness: 0.15, clearcoat: 0.8 }));
  t.position.y = 0.003;
  g.add(t);
  return g;
}

export function makeSwab({ used = false } = {}) {
  const g = new THREE.Group();
  const stick = mesh(new THREE.CylinderGeometry(0.0012, 0.0012, 0.14, 8), MAT.paper(0xf3efe4));
  stick.rotation.z = Math.PI / 2; stick.position.y = 0.003;
  const tip = mesh(new THREE.SphereGeometry(0.0035, 12, 8), MAT.paper(used ? 0xf4e7e0 : 0xffffff));
  tip.scale.x = 2; tip.position.set(0.072, 0.003, 0);
  g.add(stick, tip);
  link(g, "tip", tip);
  return g;
}

export function makeSwabCup() {
  const g = new THREE.Group();
  g.add(mesh(lathe([[0.0001, 0], [0.026, 0], [0.033, 0.09], [0.031, 0.09], [0.024, 0.004], [0.0001, 0.004]], 32), MAT.plastic(0xf2f4f6)));
  for (let i = 0; i < 6; i++) {
    const s = makeSwab();
    s.rotation.z = Math.PI / 2 - 0.08 * (i - 2.5);
    s.rotation.y = i * 1.05;
    s.position.set(Math.cos(i) * 0.008, 0.07, Math.sin(i) * 0.008);
    g.add(s);
  }
  const lab = mesh(new THREE.PlaneGeometry(0.05, 0.03), new THREE.MeshStandardMaterial({ map: labelTexture(["Sterile swabs"], { w: 256, h: 140, size: 30 }), roughness: 0.9 }), false);
  lab.position.set(0, 0.045, 0.0305); lab.rotation.x = -0.08;
  g.add(lab);
  return g;
}

export function makeForceps() {
  const g = new THREE.Group();
  for (const s of [-1, 1]) {
    const arm = mesh(new THREE.BoxGeometry(0.11, 0.0015, 0.004), MAT.steel());
    arm.position.set(0, 0.002, s * 0.003);
    arm.rotation.y = s * 0.04;
    g.add(arm);
  }
  return g;
}

export function makeScalpel() {
  const g = new THREE.Group();
  const handle = mesh(new THREE.BoxGeometry(0.11, 0.003, 0.009), MAT.steel());
  handle.position.y = 0.0015;
  const blade = mesh(new THREE.BoxGeometry(0.028, 0.0006, 0.007), new THREE.MeshStandardMaterial({ color: 0xe2e6ea, metalness: 1, roughness: 0.15 }));
  blade.position.set(0.068, 0.0015, 0.0005);
  g.add(handle, blade);
  return g;
}

export function makePencil() {
  const g = new THREE.Group();
  const body = mesh(new THREE.CylinderGeometry(0.0036, 0.0036, 0.15, 6), MAT.plastic(0xf2c230));
  body.rotation.z = Math.PI / 2; body.position.y = 0.0036;
  const tip = mesh(new THREE.ConeGeometry(0.0036, 0.016, 6), MAT.plastic(0xe7c49a));
  tip.rotation.z = -Math.PI / 2; tip.position.set(0.083, 0.0036, 0);
  g.add(body, tip);
  return g;
}

export function makePaperTowels() {
  const g = new THREE.Group();
  const stack = mesh(new THREE.BoxGeometry(0.21, 0.03, 0.11), MAT.paper(0xf2efe6));
  stack.position.y = 0.015;
  g.add(stack);
  return g;
}

export function makeGlassBin() {
  const g = new THREE.Group();
  const body = mesh(new THREE.BoxGeometry(0.3, 0.36, 0.3), MAT.plastic(0xe9eaea));
  body.position.y = 0.18;
  const lid = mesh(new THREE.BoxGeometry(0.31, 0.02, 0.31), MAT.plastic(0xf0c419));
  lid.position.y = 0.37;
  const lab = mesh(new THREE.PlaneGeometry(0.22, 0.12), new THREE.MeshStandardMaterial({ map: labelTexture(["Glass waste", "Slides and cover slips only"], { w: 320, h: 170, band: "#f0c419", size: 36 }), roughness: 0.9 }), false);
  lab.position.set(0, 0.22, 0.1506);
  g.add(body, lid, lab);
  return g;
}

// Wraparound safety goggles, as they rest on the bench with the strap behind.
export function makeGoggles() {
  const g = new THREE.Group();
  const r = 0.085, h = 0.05, th = 1.15;   // lens radius, height, half the wrap angle
  const lens = mesh(new THREE.CylinderGeometry(r, r, h, 40, 1, true, -th, th * 2), new THREE.MeshPhysicalMaterial({ color: 0xcfe9f2, roughness: 0.05, transparent: true, opacity: 0.45, clearcoat: 1, depthWrite: false, side: THREE.DoubleSide }), false);
  const frameMat = new THREE.MeshStandardMaterial({ color: 0x2f6fb8, roughness: 0.45, side: THREE.DoubleSide });
  // The soft frame: a flat curved strip along the top and bottom of the lens.
  for (const y of [-h / 2, h / 2]) {
    const rim = mesh(new THREE.RingGeometry(r - 0.016, r + 0.004, 40, 1, -Math.PI / 2 - th, th * 2), frameMat);
    rim.rotation.x = -Math.PI / 2;
    rim.position.y = y;
    g.add(rim);
    const edge = mesh(new THREE.CylinderGeometry(r + 0.004, r + 0.004, 0.006, 40, 1, true, -th, th * 2), frameMat);
    edge.position.y = y;
    g.add(edge);
  }
  for (const s of [-1, 1]) {
    const side = mesh(new THREE.BoxGeometry(0.008, h + 0.008, 0.012), frameMat);
    side.position.set(s * r * Math.sin(th), 0, r * Math.cos(th));
    side.rotation.y = s * th;
    g.add(side);
  }
  for (let i = -2; i <= 2; i++) {
    const vent = mesh(new THREE.BoxGeometry(0.008, 0.004, 0.006), frameMat);
    vent.position.set(Math.sin(i * 0.35) * (r + 0.003), h / 2 + 0.005, Math.cos(i * 0.35) * (r + 0.003));
    vent.rotation.y = i * 0.35;
    g.add(vent);
  }
  const strap = mesh(new THREE.CylinderGeometry(r - 0.002, r - 0.002, 0.016, 40, 1, true, th, Math.PI * 2 - th * 2), new THREE.MeshStandardMaterial({ color: 0x2b3550, roughness: 0.9, side: THREE.DoubleSide }));
  strap.scale.z = 0.7;
  g.add(lens, strap);
  // Centre the goggles over their footprint, resting on the lower rim.
  for (const c of g.children) c.position.add(new THREE.Vector3(0, h / 2 + 0.004, -r * 0.2));
  return g;
}

// A white pump bottle of alcohol gel.
export function makeSanitizer() {
  const g = new THREE.Group();
  const white = new THREE.MeshPhysicalMaterial({ color: 0xf3f5f7, roughness: 0.35, clearcoat: 0.6 });
  const r = 0.032, h = 0.11;
  g.add(mesh(lathe([[0.0001, 0], [r - 0.003, 0], [r, 0.004], [r, h], [r * 0.8, h + 0.015], [0.012, h + 0.02], [0.012, h + 0.028]], 40), white));
  const collar = mesh(new THREE.CylinderGeometry(0.0145, 0.0145, 0.014, 28), MAT.plastic(0xe6e8ea));
  collar.position.y = h + 0.035;
  const stem = mesh(new THREE.CylinderGeometry(0.004, 0.004, 0.02, 12), MAT.plastic(0xe6e8ea));
  stem.position.y = h + 0.05;
  const head = mesh(new THREE.BoxGeometry(0.018, 0.014, 0.022), MAT.plastic(0xf6f7f8));
  head.position.y = h + 0.066;
  const nozzle = mesh(new THREE.BoxGeometry(0.03, 0.007, 0.009), MAT.plastic(0xf6f7f8));
  nozzle.position.set(0.02, h + 0.069, 0);
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0005, r + 0.0005, 0.06, 40, 1, true, -1.1, 2.2), new THREE.MeshStandardMaterial({ map: labelTexture(["Hand sanitizer", "70% alcohol gel"], { w: 320, h: 170, band: "#2f9e6b", size: 36 }), roughness: 0.85 }), false);
  lab.position.y = 0.055;
  g.add(collar, stem, head, nozzle, lab);
  return g;
}

/* Venetian blinds `width` wide (along z) and `height` tall, origin at the
   bottom middle. The slats hang from the head rail in a group, linked as
   "stack", that is scaled in y to raise them, so the rail keeps its shape. */
export function makeBlinds({ width, height }) {
  const g = new THREE.Group();
  const mat = new THREE.MeshStandardMaterial({ color: 0xe9e8e2, roughness: 0.5, metalness: 0.1 });
  const railH = 0.04, pitch = 0.024;
  const rail = mesh(new THREE.BoxGeometry(0.05, railH, width), mat);
  rail.position.y = height - railH / 2;
  const stack = new THREE.Group();
  stack.position.y = height - railH;
  const n = Math.floor((height - railH - 0.02) / pitch);
  const slats = new THREE.InstancedMesh(new THREE.BoxGeometry(0.025, 0.0015, width - 0.01), mat, n);
  slats.castShadow = true; slats.receiveShadow = true;
  const m4 = new THREE.Matrix4(), q = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(0, 0, 1), 0.5), one = new THREE.Vector3(1, 1, 1);
  for (let i = 0; i < n; i++) slats.setMatrixAt(i, m4.compose(new THREE.Vector3(0, -(i + 1) * pitch, 0), q, one));
  const bottom = mesh(new THREE.BoxGeometry(0.03, 0.012, width - 0.01), mat);
  bottom.position.y = -(n + 1) * pitch;
  const cordMat = MAT.paper(0xdcdad2);
  for (const z of [-width * 0.35, width * 0.35]) {
    const cord = mesh(new THREE.CylinderGeometry(0.0012, 0.0012, (n + 1) * pitch, 6), cordMat, false);
    cord.position.set(0, -(n + 1) * pitch / 2, z);
    stack.add(cord);
  }
  stack.add(slats, bottom);
  const wand = mesh(new THREE.CylinderGeometry(0.004, 0.004, 0.6, 10), new THREE.MeshPhysicalMaterial({ color: 0xf2f6f8, roughness: 0.1, transparent: true, opacity: 0.6 }));
  wand.position.set(0.03, height - railH - 0.3, width / 2 - 0.05);
  g.add(rail, stack, wand);
  link(g, "stack", stack);
  return g;
}

export function makeLightSwitch() {
  const g = new THREE.Group();
  const plate = mesh(new THREE.BoxGeometry(0.086, 0.086, 0.008), MAT.plastic(0xf4f4f1));
  const rocker = mesh(new THREE.BoxGeometry(0.03, 0.05, 0.01), MAT.plastic(0xe9e9e4));
  rocker.position.z = 0.006;
  g.add(plate, rocker);
  link(g, "rocker", rocker);
  return g;
}

// A whiteboard with today's method written on it.
export function makeWhiteboard(lines) {
  const g = new THREE.Group();
  const W = 1.8, H = 1.0;
  const c = document.createElement("canvas");
  c.width = 1800; c.height = 1000;
  const x = c.getContext("2d");
  x.fillStyle = "#fbfbf8"; x.fillRect(0, 0, 1800, 1000);
  x.fillStyle = "#1f3c8c";
  x.font = "700 64px 'Chakra Petch', 'Marker Felt', sans-serif";
  x.fillText(lines[0], 70, 120);
  x.font = "44px 'Marker Felt', 'Comic Sans MS', sans-serif";
  lines.slice(1).forEach((t, i) => { x.fillStyle = i % 3 === 2 ? "#b4302a" : "#26292c"; x.fillText(t, 80, 210 + i * 66); });
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 8;
  const board = mesh(new THREE.PlaneGeometry(W, H), new THREE.MeshStandardMaterial({ map: t, roughness: 0.25 }), false);
  const frame = mesh(new THREE.BoxGeometry(W + 0.04, H + 0.04, 0.02), MAT.steel());
  frame.position.z = -0.011;
  const tray = mesh(new THREE.BoxGeometry(W * 0.9, 0.02, 0.06), MAT.steel());
  tray.position.set(0, -H / 2 - 0.02, 0.03);
  g.add(frame, board, tray);
  return g;
}

export function shortName(solution) {
  return { water: "Distilled water", tap: "Tap water", saline: "Saline 0.9%", iodine: "Iodine", mb: "Methylene blue", ethanol: "Ethanol 70%",
    disinfectant: "Disinfectant", medium: "Medium", paclitaxel: "Paclitaxel", doxorubicin: "Doxorubicin", cisplatin: "Cisplatin" }[solution] || solution;
}
