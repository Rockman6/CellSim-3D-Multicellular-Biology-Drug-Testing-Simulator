// Small lab equipment built in code at real size (metres), so each piece can
// hold a liquid that drains and changes colour. Every prop's origin is the
// middle of its base, so it can be set down on any surface. Glassware has
// real wall thickness and a rounded lip, which is what makes glass read as
// glass: thick bright edges, thin clear walls.
import * as THREE from "three";
import { Container, SOLUTIONS, liquidColumn } from "./liquids.js";
import { link } from "./assets.js";
import { mergeGeometries } from "three/addons/utils/BufferGeometryUtils.js";

const V2 = THREE.Vector2;
export const MAT = {
  // forceSinglePass: both faces in one draw, not two (three.js draws
  // see-through two-sided things twice otherwise).
  glass: () => new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.03, metalness: 0, transparent: true, opacity: 0.17, ior: 1.5, clearcoat: 1, clearcoatRoughness: 0.04, envMapIntensity: 1.4, depthWrite: false, side: THREE.DoubleSide, forceSinglePass: true }),
  plasticClear: () => new THREE.MeshPhysicalMaterial({ color: 0xf4f8fa, roughness: 0.12, transparent: true, opacity: 0.22, clearcoat: 0.6, depthWrite: false, side: THREE.DoubleSide, forceSinglePass: true }),
  amber: () => new THREE.MeshPhysicalMaterial({ color: 0x6b3510, roughness: 0.06, transparent: true, opacity: 0.74, clearcoat: 1, envMapIntensity: 1.3, depthWrite: false }),
  blueGlass: () => new THREE.MeshPhysicalMaterial({ color: 0x1d3f8c, roughness: 0.06, transparent: true, opacity: 0.72, clearcoat: 1, envMapIntensity: 1.3, depthWrite: false }),
  steel: () => new THREE.MeshStandardMaterial({ color: 0xc9ced3, metalness: 1, roughness: 0.28 }),
  brushed: () => new THREE.MeshStandardMaterial({ color: 0xb8bdc2, metalness: 0.9, roughness: 0.42 }),
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

/* Merge a group's plain child meshes that share a material into one mesh
   each. It looks the same and costs one draw call instead of many, which
   matters on older laptops. Meshes with their own material are left alone. */
function compact(g) {
  const byMat = new Map();
  for (const m of g.children) {
    if (!m.isMesh || m.isInstancedMesh || m.children.length || !m.visible) continue;
    if (!byMat.has(m.material)) byMat.set(m.material, []);
    byMat.get(m.material).push(m);
  }
  for (const [mat, list] of byMat) {
    if (list.length < 2) continue;
    const geos = list.map((m) => {
      m.updateMatrix();
      const geo = (m.geometry.index ? m.geometry.toNonIndexed() : m.geometry.clone()).applyMatrix4(m.matrix);
      for (const k of Object.keys(geo.attributes)) if (!["position", "normal", "uv"].includes(k)) geo.deleteAttribute(k);
      if (!geo.attributes.uv) geo.setAttribute("uv", new THREE.Float32BufferAttribute(new Float32Array(geo.attributes.position.count * 2), 2));
      return geo;
    });
    const merged = mergeGeometries(geos, false);
    if (!merged) continue;
    for (const m of list) g.remove(m);
    g.add(mesh(merged, mat, list.some((m) => m.castShadow)));
  }
  return g;
}

// A label printed on paper, as a canvas texture.
export function labelTexture(lines, { w = 256, h = 128, bg = "#f7f5ef", fg = "#1d2326", band = null, size = 30, hazard = null } = {}) {
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
  // A small hazard diamond in the corner: red rim, black symbol.
  if (hazard) {
    const d = h * 0.26, x = w - d * 0.85, y = h - d * 0.85;
    g.save(); g.translate(x, y); g.rotate(Math.PI / 4);
    g.fillStyle = "#fff"; g.fillRect(-d / 2, -d / 2, d, d);
    g.lineWidth = d * 0.12; g.strokeStyle = "#d1242a"; g.strokeRect(-d / 2, -d / 2, d, d);
    g.restore();
    g.fillStyle = "#111"; g.font = `700 ${d * 0.6}px sans-serif`; g.fillText(hazard, x, y + d * 0.04);
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.anisotropy = 4;
  return t;
}

function lathe(points, segments = 32) {
  return new THREE.LatheGeometry(points.map(([r, y]) => new V2(r, y)), segments);
}

/* A vessel with real walls: `outer` runs from the middle of the base out and
   up to the rim; the inside follows it `t` in, and a rounded lip joins them. */
function vessel(outer, t, segments = 40) {
  const pts = outer.map(([r, y]) => [r, y]);
  const [rt, yt] = outer[outer.length - 1];
  pts.push([rt - t * 0.15, yt + t * 0.45], [rt - t * 0.5, yt + t * 0.55], [rt - t * 0.85, yt + t * 0.45]);
  for (let i = outer.length - 1; i >= 1; i--) {
    const [r, y] = outer[i];
    pts.push([Math.max(0.0002, r - t), Math.max(y, t)]);
  }
  pts.push([0.0001, t]);
  return lathe(pts, segments);
}

// Push the rim out into a pouring spout, centred on +z.
function addSpout(geo, yFrom, yTo, amount) {
  const p = geo.attributes.position;
  for (let i = 0; i < p.count; i++) {
    const y = p.getY(i);
    if (y < yFrom) continue;
    const x = p.getX(i), z = p.getZ(i), r = Math.hypot(x, z);
    if (r < 1e-5) continue;
    const w = Math.max(0, 1 - Math.abs(Math.atan2(x, z)) / 0.5);
    const k = amount * w * w * Math.min(1, (y - yFrom) / (yTo - yFrom));
    p.setXYZ(i, (x / r) * (r + k), y + k * 0.4, (z / r) * (r + k));
  }
  geo.computeVertexNormals();
}

// A screw cap with grip ribs round its side.
function ribbedCap(r, h, color, ribs = 30) {
  const geo = new THREE.CylinderGeometry(r, r, h, ribs * 4, 1);
  const p = geo.attributes.position;
  for (let i = 0; i < p.count; i++) {
    const x = p.getX(i), z = p.getZ(i), rr = Math.hypot(x, z);
    if (rr < r * 0.98) continue;   // leave the flat top and bottom alone
    const k = 1 + 0.045 * Math.max(0, Math.cos(Math.atan2(x, z) * ribs));
    p.setX(i, x * k); p.setZ(i, z * k);
  }
  geo.computeVertexNormals();
  return mesh(geo, MAT.plastic(color));
}

// White enamel scale marks printed on glass, every `step` up to `max`.
function scaleMarks(r, h, max, step, { unit = "mL", y0 = 0.1, span = 0.8, arc = 0.75 } = {}) {
  const c = document.createElement("canvas");
  c.width = 96; c.height = 384;
  const g = c.getContext("2d");
  g.fillStyle = "rgba(255,255,255,0.95)";
  g.font = "600 22px 'IBM Plex Sans', sans-serif";
  const n = Math.round(max / step);
  for (let i = 1; i <= n; i++) {
    const y = 384 * (1 - y0) - (i / n) * 384 * span;
    const major = i % 2 === 0 || n <= 5;
    g.fillRect(0, y - 1.5, major ? 38 : 22, 3);
    if (major) g.fillText(`${i * step}`, 44, y + 7);
  }
  g.fillText(unit, 44, 384 * (1 - y0) - 384 * span - 26);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  const m = mesh(new THREE.CylinderGeometry(r + 0.0002, r + 0.0002, h, 24, 1, true, -arc / 2, arc), new THREE.MeshBasicMaterial({ map: t, transparent: true, depthWrite: false }), false);
  m.position.y = h / 2;
  return m;
}

// A rounded rectangle outline, for extruded trays, tiles and boxes.
function roundedRect(w, d, r) {
  const s = new THREE.Shape(), x = -w / 2, y = -d / 2;
  r = Math.min(r, w / 2, d / 2);
  s.moveTo(x + r, y); s.lineTo(x + w - r, y); s.quadraticCurveTo(x + w, y, x + w, y + r);
  s.lineTo(x + w, y + d - r); s.quadraticCurveTo(x + w, y + d, x + w - r, y + d);
  s.lineTo(x + r, y + d); s.quadraticCurveTo(x, y + d, x, y + d - r);
  s.lineTo(x, y + r); s.quadraticCurveTo(x, y, x + r, y);
  return s;
}

// A box with softly rounded edges, origin at the middle of its base.
export function roundedBoxGeometry(w, h, d, r = Math.min(w, d) * 0.08, bevel = Math.min(h * 0.25, r * 0.6, 0.002)) {
  const geo = new THREE.ExtrudeGeometry(roundedRect(w - 2 * bevel, d - 2 * bevel, Math.max(0.0002, r - bevel)),
    { depth: Math.max(0.0001, h - 2 * bevel), bevelEnabled: bevel > 0, bevelThickness: bevel, bevelSize: bevel, bevelSegments: 2, curveSegments: 6 });
  geo.rotateX(-Math.PI / 2);
  geo.translate(0, bevel, 0);
  return geo;
}

/* A soft dark patch under a thing standing on a surface: on every graphics
   setting, it is what makes objects sit on the bench instead of floating. */
let shadowTex = null, shadowMat = null;
export function contactShadow(w, d) {
  if (!shadowTex) {
    const c = document.createElement("canvas");
    c.width = c.height = 64;
    const g = c.getContext("2d");
    const grd = g.createRadialGradient(32, 32, 4, 32, 32, 32);
    grd.addColorStop(0, "rgba(0,0,0,0.85)"); grd.addColorStop(0.55, "rgba(0,0,0,0.35)"); grd.addColorStop(1, "rgba(0,0,0,0)");
    g.fillStyle = grd; g.fillRect(0, 0, 64, 64);
    shadowTex = new THREE.CanvasTexture(c);
    shadowMat = new THREE.MeshBasicMaterial({ map: shadowTex, color: 0x000000, transparent: true, opacity: 0.5, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2 });
  }
  const m = new THREE.Mesh(new THREE.PlaneGeometry(1, 1), shadowMat);
  m.rotation.x = -Math.PI / 2;
  m.scale.set(w * 1.45 + 0.01, d * 1.45 + 0.01, 1);
  m.position.y = 0.0004;
  m.renderOrder = 1;
  m.name = "shadow";
  m.raycast = () => {};
  return m;
}

/* ---------- slide, cover slip and their boxes ---------- */
export function makeSlide() {
  const g = new THREE.Group();
  const L = 0.076, W = 0.026, T = 0.001;
  const glass = mesh(roundedBoxGeometry(L, T, W, 0.0008, 0.0003), MAT.glass());
  const frost = mesh(new THREE.BoxGeometry(0.02, T * 1.05, W - 0.0004), MAT.paper(0xeef0ef));
  frost.position.set(-L / 2 + 0.0105, T / 2, 0);
  g.add(glass, frost);
  // The sample, the drop of liquid and the cover slip, shown as they are added.
  const sample = mesh(new THREE.CircleGeometry(0.0065, 24), new THREE.MeshStandardMaterial({ color: 0xe9dfc8, roughness: 0.6, transparent: true, opacity: 0.75 }), false);
  sample.rotation.x = -Math.PI / 2; sample.position.set(0.008, T + 0.0002, 0); sample.visible = false;
  const drop = mesh(new THREE.SphereGeometry(1, 20, 10, 0, Math.PI * 2, 0, Math.PI / 2), new THREE.MeshPhysicalMaterial({ color: 0xd7ebf2, roughness: 0.05, transparent: true, opacity: 0.5, clearcoat: 1, depthWrite: false }), false);
  drop.position.set(0.008, T, 0); drop.visible = false; drop.renderOrder = 3; drop.scale.setScalar(0.003);
  const cover = mesh(roundedBoxGeometry(0.022, 0.00025, 0.022, 0.0004, 0.0001), MAT.glass(), false);
  cover.position.set(0.008, T + 0.0005, 0); cover.visible = false;
  const label = mesh(new THREE.PlaneGeometry(0.018, 0.022), new THREE.MeshStandardMaterial({ transparent: true, roughness: 0.9 }), false);
  label.rotation.x = -Math.PI / 2; label.position.set(-L / 2 + 0.0105, T * 1.1, 0); label.visible = false;
  g.add(sample, drop, cover, label);
  link(g, "parts", { sample, drop, cover, label });
  return g;
}

const SAMPLE_LOOK = {
  elodea: { color: 0x4a9a3e, sx: 1.2, sz: 0.55 }, cheek: { color: 0xf0dcd4, sx: 1, sz: 1 }, onion: { color: 0xece3cf, sx: 1, sz: 1 },
  redonion: { color: 0xb3416d, sx: 1, sz: 0.9 }, root: { color: 0xf2efe6, sx: 0.55, sz: 0.2 },
};
export function updateSlideLook(slide) {
  const { sample, drop, cover, label } = slide.obj.userData.parts;
  const s = slide.state;
  sample.visible = !!s.sample;
  if (s.sample) {
    const look = SAMPLE_LOOK[s.sample] || SAMPLE_LOOK.onion;
    sample.material.color.set(look.color);
    // A squashed root tip spreads into a thin pale smear.
    const k = s.squashed ? 1.7 : 1;
    sample.scale.set(look.sx * k, look.sz * k * (s.squashed ? 2.4 : 1), 1);
    if (s.sample === "root" && slide.container.stain() === "toluidine") sample.material.color.set(0x7d5aa6);
  }
  const v = slide.container.volume;
  drop.visible = v > 0.001;
  if (drop.visible) {
    const { color, alpha } = slide.container.color();
    drop.material.color.copy(color);
    drop.material.opacity = Math.min(0.9, 0.3 + alpha * 0.6);
    const r = 0.0035 + Math.sqrt(v) * 0.012;               // a 50 µL drop is about 6 mm across
    drop.scale.set(s.cover ? 0.0105 : r, s.cover ? 0.0003 : r * 0.35, s.cover ? 0.0105 : r);
  }
  cover.visible = !!s.cover;
  if (s.label && !label.material.map) {
    label.material.map = labelTexture([s.label], { w: 128, h: 160, bg: "#eef0ef", size: 26 });
    label.material.needsUpdate = true;
  }
  label.visible = !!s.label;
}

// A cardboard or plastic box with a lid seam and a printed label on the front and top.
export function makeBox({ w, h, d, color, labelLines, band }) {
  const g = new THREE.Group();
  const mat = MAT.plastic(color);
  mat.roughness = 0.7;
  const body = mesh(roundedBoxGeometry(w, h * 0.72, d, 0.0015, 0.0006), mat);
  const lid = mesh(roundedBoxGeometry(w + 0.0012, h * 0.34, d + 0.0012, 0.0018, 0.0007), mat);
  lid.position.y = h * 0.66;
  g.add(body, lid);
  if (labelLines) {
    const tex = labelTexture(labelLines, { band });
    const front = mesh(new THREE.PlaneGeometry(w * 0.82, h * 0.5), new THREE.MeshStandardMaterial({ map: tex, roughness: 0.85 }), false);
    front.position.set(0, h * 0.36, d / 2 + 0.0004);
    const top = mesh(new THREE.PlaneGeometry(w * 0.82, d * 0.7), new THREE.MeshStandardMaterial({ map: tex, roughness: 0.85 }), false);
    top.rotation.x = -Math.PI / 2; top.position.y = h + 0.0002;
    g.add(front, top);
  }
  return compact(g);
}

export function makeCoverSlip() {
  const g = new THREE.Group();
  const s = mesh(roundedBoxGeometry(0.022, 0.0004, 0.022, 0.0004, 0.00012), MAT.glass(), false);
  g.add(s);
  return g;
}

/* ---------- bottles and beakers ---------- */
function bottleGlass(solution) {
  const kind = SOLUTIONS[solution]?.glass;
  return kind === "amber" ? MAT.amber() : kind === "blue" ? MAT.blueGlass() : MAT.glass();
}

// A 30 mL dropper bottle: glass with a ribbed collar, a rubber teat, and the
// glass dropper tube visible inside. Its liquid follows the volume.
export function makeDropperBottle(solution, volume = 25) {
  const s = SOLUTIONS[solution];
  const g = new THREE.Group();
  const r = 0.0175, h = 0.062;
  g.add(mesh(vessel([[0.0001, 0], [r - 0.0022, 0], [r, 0.0022], [r, h], [r * 0.86, h + 0.006], [r * 0.52, h + 0.0105], [0.0079, h + 0.0125], [0.0079, h + 0.019]], 0.0013), bottleGlass(solution)));
  const collar = ribbedCap(0.0097, 0.012, 0x24282c, 24);
  collar.position.y = h + 0.024;
  const teat = mesh(lathe([[0.0001, 0], [0.0062, 0], [0.0068, 0.004], [0.0083, 0.01], [0.0085, 0.016], [0.0069, 0.022], [0.0036, 0.0262], [0.0001, 0.027]], 24), MAT.rubber());
  teat.position.y = h + 0.029;
  const tube = mesh(new THREE.CylinderGeometry(0.0017, 0.0011, h + 0.02, 12, 1, true), MAT.glass(), false);
  tube.position.y = (h + 0.02) / 2 + 0.009;
  g.add(collar, teat, tube);
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0004, r + 0.0004, 0.032, 32, 1, true, -0.95, 1.9), new THREE.MeshStandardMaterial({ map: labelTexture([shortName(solution), s.note || "30 mL"], { size: 34, band: s.color, hazard: s.hazard }), roughness: 0.85 }), false);
  lab.position.y = 0.031;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0016, bottom: 0.0014, depth: h - 0.003 });
  g.add(liq.mesh);
  const container = new Container(30, { [solution]: volume });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

// A 500 mL screw-cap reagent bottle with a pour ring, scale and a printed label.
export function makeStockBottle(solution, volume = 450, { capacity = 500, scale = 1 } = {}) {
  const s = SOLUTIONS[solution];
  const g = new THREE.Group();
  const r = 0.042 * scale, h = 0.15 * scale;
  g.add(mesh(vessel([[0.0001, 0], [r - 0.004, 0], [r, 0.004], [r, h], [r * 0.9, h + 0.012 * scale], [r * 0.66, h + 0.022 * scale], [0.02 * scale, h + 0.028 * scale], [0.02 * scale, h + 0.04 * scale]], 0.0022 * scale, 48), bottleGlass(solution)));
  const ring = mesh(new THREE.CylinderGeometry(0.0235 * scale, 0.0235 * scale, 0.006 * scale, 32), MAT.plastic(0xf0f0ea));
  ring.position.y = h + 0.03 * scale;
  const cap = ribbedCap(0.0245 * scale, 0.024 * scale, s.cap || 0x2a5fb0, 36);
  cap.position.y = h + 0.045 * scale;
  g.add(ring, cap);
  g.add(scaleMarks(r, h, capacity, capacity / 5, { y0: 0.12, span: 0.78, arc: 0.55 }));
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0006, r + 0.0006, 0.07 * scale, 40, 1, true, Math.PI - 1.1, 2.2), new THREE.MeshStandardMaterial({ map: labelTexture([shortName(solution), s.note || `Stock · ${capacity} mL`], { w: 320, h: 160, band: s.color, size: 36, hazard: s.hazard }), roughness: 0.85 }), false);
  lab.position.y = 0.065 * scale;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0026 * scale, bottom: 0.0024 * scale, depth: h - 0.006 * scale });
  g.add(liq.mesh);
  const container = new Container(capacity, { [solution]: volume });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

// A 10 mL amber vial of a cytotoxic drug stock, with a dropper cap.
export function makeDrugVial(drug, volume = 8) {
  const s = SOLUTIONS[drug];
  const g = new THREE.Group();
  const r = 0.012, h = 0.045;
  g.add(mesh(vessel([[0.0001, 0], [r - 0.0015, 0], [r, 0.0015], [r, h], [r * 0.86, h + 0.004], [0.0068, h + 0.007], [0.0068, h + 0.012]], 0.001), MAT.amber()));
  const cap = ribbedCap(0.0085, 0.01, 0xe9e9e3, 20);
  cap.position.y = h + 0.017;
  const tip = mesh(lathe([[0.0001, 0], [0.0042, 0], [0.0024, 0.008], [0.0011, 0.013], [0.0001, 0.0135]], 16), MAT.plastic(0xe9e9e3));
  tip.position.y = h + 0.022;
  g.add(cap, tip);
  const lab = mesh(new THREE.CylinderGeometry(r + 0.0004, r + 0.0004, 0.028, 32, 1, true, -1, 2), new THREE.MeshStandardMaterial({ map: labelTexture([shortName(drug), `${s.uM} µM · CYTOTOXIC`], { w: 256, h: 128, band: "#d8342a", size: 30, hazard: "!" }), roughness: 0.85 }), false);
  lab.position.y = 0.024;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0012, bottom: 0.0012, depth: h - 0.003 });
  g.add(liq.mesh);
  const container = new Container(10, { [drug]: volume });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

// A 250 mL low-form beaker with a pouring spout and an enamel scale.
export function makeBeaker(contents = {}, capacity = 250) {
  const g = new THREE.Group();
  const r = 0.035, h = 0.095;
  const geo = vessel([[0.0001, 0], [r - 0.003, 0], [r, 0.003], [r, h - 0.002], [r + 0.0022, h + 0.0015]], 0.0016, 48);
  addSpout(geo, h - 0.012, h + 0.003, 0.007);
  g.add(mesh(geo, MAT.glass()));
  const marks = scaleMarks(r, h, capacity, 50, { y0: 0.04, span: capacity / 300 * 0.92, arc: 0.7 });
  marks.rotation.y = Math.PI;
  g.add(marks);
  const liq = liquidColumn({ radius: r - 0.0018, bottom: 0.0017, depth: h - 0.006 });
  g.add(liq.mesh);
  const container = new Container(capacity, contents);
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

/* ---------- culture dish ---------- */
// A 60 mm polystyrene culture dish with its lid on and a handwritten label.
export function makeCultureDish(name, line = "A549 lung cells") {
  const g = new THREE.Group();
  const r = 0.03, h = 0.0125;
  const base = mesh(vessel([[0.0001, 0], [r - 0.0008, 0], [r, 0.0008], [r, h]], 0.0008, 56), MAT.plasticClear(), false);
  const lid = mesh(lathe([[0.0001, h + 0.0016], [r + 0.0016, h + 0.0016], [r + 0.0019, h + 0.0012], [r + 0.0019, h * 0.45], [r + 0.0012, h * 0.45], [r + 0.0012, h + 0.0009], [0.0001, h + 0.0009]], 56), MAT.plasticClear(), false);
  g.add(base, lid);
  const lab = mesh(new THREE.PlaneGeometry(0.03, 0.012), new THREE.MeshStandardMaterial({ map: labelTexture([name, line], { w: 256, h: 100, size: 30 }), roughness: 0.9, transparent: true }), false);
  lab.rotation.x = -Math.PI / 2; lab.position.y = h + 0.0019;
  g.add(lab);
  const liq = liquidColumn({ radius: r - 0.0012, bottom: 0.0009, depth: 0.003 });
  g.add(liq.mesh);
  const container = new Container(5.2, { medium: 5 });
  liq.update(container);
  return { obj: g, container, liquid: liq };
}

/* ---------- 6-well plate ---------- */
// A standard 6-well plate (127.8 × 85.5 mm), each well with 2 mL of medium.
// Wells are A1–A3 along the far row and B1–B3 along the near one.
export const WELLS = ["A1", "A2", "A3", "B1", "B2", "B3"];
export function makeWellPlate(name) {
  const g = new THREE.Group();
  const W = 0.1278, D = 0.0855, H = 0.0205, R = 0.0178, pitch = 0.039;
  const clear = MAT.plasticClear();
  const wall = new THREE.Shape(); wall.copy(roundedRect(W, D, 0.005));
  wall.holes.push(roundedRect(W - 0.0024, D - 0.0024, 0.004));
  const skirt = new THREE.ExtrudeGeometry(wall, { depth: H, bevelEnabled: false, curveSegments: 6 });
  skirt.rotateX(-Math.PI / 2);
  g.add(mesh(skirt, clear, false));
  const deckShape = roundedRect(W - 0.0024, D - 0.0024, 0.004);
  const centres = WELLS.map((_, i) => new THREE.Vector3((i % 3 - 1) * pitch, 0, (i < 3 ? -0.5 : 0.5) * pitch));
  for (const c of centres) { const hole = new THREE.Path(); hole.absarc(c.x, -c.z, R, 0, Math.PI * 2, true); deckShape.holes.push(hole); }
  const deck = new THREE.ExtrudeGeometry(deckShape, { depth: 0.0012, bevelEnabled: false, curveSegments: 24 });
  deck.rotateX(-Math.PI / 2); deck.translate(0, H - 0.0012, 0);
  g.add(mesh(deck, clear, false));
  const wells = centres.map((c, i) => {
    const wallM = mesh(new THREE.CylinderGeometry(R, R, H - 0.0015, 40, 1, true), clear, false);
    wallM.position.set(c.x, 0.0015 + (H - 0.0015) / 2, c.z);
    const floor = mesh(new THREE.CircleGeometry(R, 40), clear, false);
    floor.rotation.x = -Math.PI / 2; floor.position.set(c.x, 0.0015, c.z);
    g.add(wallM, floor);
    const liq = liquidColumn({ radius: R - 0.0004, bottom: 0.0016, depth: 0.0035, segments: 32 });
    liq.mesh.position.x = c.x; liq.mesh.position.z = c.z;
    g.add(liq.mesh);
    const container = new Container(3, { medium: 2 });
    liq.update(container);
    return { name: WELLS[i], centre: c, container, liquid: liq };
  });
  // Moulded letters and numbers along the edges, and a label on the lid.
  const c = document.createElement("canvas");
  c.width = 512; c.height = 342;
  const x = c.getContext("2d");
  x.fillStyle = "rgba(255,255,255,0.85)"; x.font = "600 26px 'IBM Plex Sans', sans-serif"; x.textAlign = "center"; x.textBaseline = "middle";
  [0, 1, 2].forEach((k) => x.fillText(`${k + 1}`, 256 + (k - 1) * 156, 16));
  ["A", "B"].forEach((k, j) => x.fillText(k, 14, 171 + (j - 0.5) * 156));
  x.fillStyle = "#f7f5ef"; x.fillRect(300, 300, 200, 38);
  x.fillStyle = "#1d2326"; x.font = "700 24px 'IBM Plex Sans', sans-serif"; x.fillText(name, 400, 320);
  const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
  const print = mesh(new THREE.PlaneGeometry(W - 0.003, D - 0.003), new THREE.MeshBasicMaterial({ map: tex, transparent: true, depthWrite: false }), false);
  print.rotation.x = -Math.PI / 2; print.position.y = H + 0.0026;
  const lidTop = mesh(roundedBoxGeometry(W + 0.002, 0.0011, D + 0.002, 0.006, 0.0003), clear, false);
  lidTop.position.y = H + 0.0014;
  g.add(lidTop, print);
  compact(g);
  link(g, "wells", wells);
  return { obj: g, wells };
}

// The well nearest a point given in the plate's own space.
export function wellAt(plateObj, localPoint) {
  let best = 0, d2 = Infinity;
  plateObj.userData.wells.forEach((w, i) => { const d = (w.centre.x - localPoint.x) ** 2 + (w.centre.z - localPoint.z) ** 2; if (d < d2) { d2 = d; best = i; } });
  return best;
}

/* ---------- micropipette, tips and tip waste ---------- */
// A 2–200 µL air-displacement micropipette. Its display shows the set volume;
// a yellow tip appears when one is fitted, with any liquid drawn up in it.
export function makeMicropipette() {
  const g = new THREE.Group();
  const body = new THREE.Group();
  g.add(body);
  const grey = MAT.plastic(0xdfe3e7), blue = MAT.plastic(0x41607f), white = MAT.plastic(0xf3f4f2);
  grey.roughness = 0.35; blue.roughness = 0.4;
  body.add(mesh(lathe([[0.0001, 0], [0.0023, 0], [0.0031, 0.02], [0.0046, 0.045], [0.0001, 0.045]], 24), grey));
  body.add(mesh(lathe([[0.0001, 0.045], [0.0058, 0.045], [0.0064, 0.05], [0.0064, 0.074], [0.0001, 0.074]], 24), white));
  body.add(mesh(lathe([[0.0001, 0.074], [0.0098, 0.076], [0.0112, 0.1], [0.0126, 0.15], [0.0122, 0.185], [0.0106, 0.196], [0.0001, 0.197]], 32), blue));
  const hook = mesh(roundedBoxGeometry(0.008, 0.012, 0.016, 0.003, 0.0015), blue);
  hook.position.set(0, 0.168, 0.016);
  const ejector = mesh(roundedBoxGeometry(0.006, 0.01, 0.01, 0.002, 0.001), white);
  ejector.position.set(0, 0.188, -0.013);
  const rod = mesh(new THREE.CylinderGeometry(0.0028, 0.0028, 0.016, 12), grey);
  rod.position.y = 0.205;
  const button = mesh(lathe([[0.0001, 0], [0.0069, 0], [0.0071, 0.006], [0.0062, 0.01], [0.0001, 0.0105]], 24), MAT.plastic(0xf2c230));
  button.position.y = 0.212;
  body.add(hook, ejector, rod, button);
  // The volume counter, behind a small window on the handle.
  const dc = document.createElement("canvas");
  dc.width = 128; dc.height = 48;
  const dtex = new THREE.CanvasTexture(dc); dtex.colorSpace = THREE.SRGBColorSpace;
  const display = mesh(new THREE.PlaneGeometry(0.012, 0.0045), new THREE.MeshBasicMaterial({ map: dtex }), false);
  display.position.set(0, 0.128, 0.0122);
  body.add(display);
  const tip = mesh(lathe([[0.0004, -0.052], [0.0009, -0.05], [0.0024, -0.012], [0.0029, 0.002], [0.0029, 0.009]], 20), new THREE.MeshPhysicalMaterial({ color: 0xf6d34a, roughness: 0.2, transparent: true, opacity: 0.55, depthWrite: false, side: THREE.DoubleSide, forceSinglePass: true }), false);
  tip.visible = false;
  const load = mesh(new THREE.ConeGeometry(0.0016, 0.012, 14), new THREE.MeshPhysicalMaterial({ color: 0xffffff, transparent: true, opacity: 0.85, depthWrite: false }), false);
  load.rotation.x = Math.PI; load.position.y = -0.044; load.visible = false;
  body.add(tip, load);
  compact(body);
  link(g, "parts", { body, tip, load, dc, dtex });
  link(g, "pose", (p) => {
    // Lying on the bench, or upright in the hand.
    if (p === "hand") { body.rotation.set(0, 0, 0); body.position.set(0, 0, 0); }
    else { body.rotation.set(0, 0, -Math.PI / 2); body.position.set(-0.1, 0.0126, 0); }
  });
  g.userData.pose("rest");
  return g;
}

export function updatePipetteLook(it) {
  const { tip, load, dc, dtex } = it.obj.userData.parts;
  const st = it.state;
  tip.visible = !!st.tip;
  load.visible = !!(st.tip && st.load && !st.load.empty);
  if (load.visible) { const { color } = st.load.color(); load.material.color.copy(color); }
  if (st.shown !== st.vol) {
    st.shown = st.vol;
    const g = dc.getContext("2d");
    g.fillStyle = "#e9ecd8"; g.fillRect(0, 0, 128, 48);
    g.fillStyle = "#1a1d1f"; g.font = "600 38px 'IBM Plex Mono', monospace"; g.textAlign = "center"; g.textBaseline = "middle";
    g.fillText(String(Math.round(st.vol)).padStart(3, "0"), 64, 26);
    dtex.needsUpdate = true;
  }
}

// A rack of 96 yellow 200 µL tips with its lid open.
export function makeTipBox() {
  const g = new THREE.Group();
  const W = 0.08, D = 0.122, H = 0.055;
  g.add(mesh(roundedBoxGeometry(W, H, D, 0.004, 0.0015), new THREE.MeshPhysicalMaterial({ color: 0x2f74c8, roughness: 0.3, transparent: true, opacity: 0.85 })));
  const rack = mesh(roundedBoxGeometry(W - 0.003, 0.003, D - 0.003, 0.003, 0.0008), MAT.plastic(0xf1f2ef));
  rack.position.y = H;
  g.add(rack);
  const tipGeo = lathe([[0.0004, -0.04], [0.0009, -0.038], [0.0024, -0.005], [0.0029, 0.006], [0.0029, 0.012], [0.0025, 0.012]], 10);
  const tips = new THREE.InstancedMesh(tipGeo, new THREE.MeshStandardMaterial({ color: 0xf6d34a, roughness: 0.35, transparent: true, opacity: 0.9 }), 96);
  const m4 = new THREE.Matrix4();
  let k = 0;
  for (let i = 0; i < 12; i++) for (let j = 0; j < 8; j++) tips.setMatrixAt(k++, m4.makeTranslation((j - 3.5) * 0.009, H + 0.003, (i - 5.5) * 0.009));
  tips.castShadow = true;
  g.add(tips);
  const lid = mesh(roundedBoxGeometry(W + 0.002, 0.016, D + 0.002, 0.004, 0.0012), MAT.plasticClear(), false);
  const hinge = new THREE.Group();
  hinge.position.set(0, H + 0.003, -D / 2);
  lid.position.z = D / 2 + 0.001;
  lid.rotation.x = Math.PI;
  hinge.add(lid);
  hinge.rotation.x = -1.95;
  g.add(hinge);
  link(g, "tips", tips);
  return g;
}

export function makeTipWaste() {
  const g = new THREE.Group();
  g.add(mesh(vessel([[0.0001, 0], [0.034, 0], [0.037, 0.003], [0.04, 0.09]], 0.002, 32), MAT.plastic(0xf2f3f1)));
  const lab = mesh(new THREE.CylinderGeometry(0.0395, 0.038, 0.035, 32, 1, true, -0.9, 1.8), new THREE.MeshStandardMaterial({ map: labelTexture(["Used tips", "Tip waste"], { band: "#d8342a" }), roughness: 0.85 }), false);
  lab.position.y = 0.05;
  g.add(lab);
  // A few used tips already inside.
  const yellow = MAT.plastic(0xf6d34a);
  for (let i = 0; i < 7; i++) {
    const t = mesh(new THREE.ConeGeometry(0.0028, 0.05, 10), yellow, false);
    t.position.set(Math.cos(i * 2.3) * 0.015, 0.025, Math.sin(i * 2.3) * 0.015);
    t.rotation.set(Math.PI + Math.sin(i) * 0.6, i, Math.cos(i * 1.7) * 0.6);
    g.add(t);
  }
  return compact(g);
}

/* ---------- plate reader ---------- */
export function makePlateReader() {
  const g = new THREE.Group();
  const W = 0.36, H = 0.16, D = 0.3;
  g.add(mesh(roundedBoxGeometry(W, H, D, 0.02, 0.006), MAT.plastic(0xe8eaec)));
  const slot = mesh(roundedBoxGeometry(0.16, 0.026, 0.004, 0.004, 0.001), MAT.plastic(0x23272b), false);
  slot.position.set(-0.06, 0.035, D / 2);
  const tray = mesh(roundedBoxGeometry(0.15, 0.006, 0.05, 0.003, 0.001), MAT.plastic(0x30353a));
  tray.position.set(-0.06, 0.04, D / 2 + 0.02);
  g.add(slot, tray);
  const c = document.createElement("canvas");
  c.width = 256; c.height = 128;
  const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
  const screen = mesh(new THREE.PlaneGeometry(0.09, 0.045), new THREE.MeshBasicMaterial({ map: tex, toneMapped: false }), false);
  screen.position.set(0.1, H * 0.62, D / 2 + 0.0008);
  g.add(screen);
  const go = MAT.plastic(0x58c46e), grey = MAT.plastic(0x8a9096);
  for (let i = 0; i < 3; i++) {
    const b = mesh(new THREE.CylinderGeometry(0.006, 0.006, 0.004, 16), i === 0 ? go : grey);
    b.rotation.x = Math.PI / 2; b.position.set(0.075 + i * 0.022, H * 0.3, D / 2 + 0.002);
    g.add(b);
  }
  const lab = mesh(new THREE.PlaneGeometry(0.14, 0.022), new THREE.MeshStandardMaterial({ map: labelTexture(["Plate reader", "Cell viability"], { w: 320, h: 64, size: 26, bg: "#e8eaec" }), roughness: 0.6 }), false);
  lab.position.set(-0.06, H * 0.75, D / 2 + 0.0008);
  g.add(lab);
  const show = (lines, color = "#7dffb0") => {
    const x = c.getContext("2d");
    x.fillStyle = "#11171a"; x.fillRect(0, 0, 256, 128);
    x.fillStyle = color; x.font = "600 30px 'IBM Plex Mono', monospace";
    lines.forEach((t, i) => x.fillText(t, 14, 44 + i * 40));
    tex.needsUpdate = true;
  };
  show(["READY", "Insert a plate"]);
  link(g, "screen", show);
  return compact(g);
}

/* ---------- hemocytometer and microtubes ---------- */
// A Neubauer counting chamber: a thick glass slide with two mirrored, ruled
// counting areas under a thick cover glass.
export function makeHemocytometer() {
  const g = new THREE.Group();
  const L = 0.075, W = 0.032, T = 0.004;
  g.add(mesh(roundedBoxGeometry(L, T, W, 0.001, 0.0004), MAT.glass()));
  const c = document.createElement("canvas");
  c.width = 128; c.height = 128;
  const x = c.getContext("2d");
  x.fillStyle = "#cfd6da"; x.fillRect(0, 0, 128, 128);
  x.strokeStyle = "rgba(40,50,58,0.75)"; x.lineWidth = 1;
  for (let i = 0; i <= 12; i++) { x.beginPath(); x.moveTo(16 + i * 8, 16); x.lineTo(16 + i * 8, 112); x.stroke(); x.beginPath(); x.moveTo(16, 16 + i * 8); x.lineTo(112, 16 + i * 8); x.stroke(); }
  const tex = new THREE.CanvasTexture(c);
  const mirror = new THREE.MeshStandardMaterial({ map: tex, metalness: 0.9, roughness: 0.15 }), slot = MAT.plastic(0x5c666c);
  for (const s of [-1, 1]) {
    const pad = mesh(new THREE.PlaneGeometry(0.012, 0.009), mirror, false);
    pad.rotation.x = -Math.PI / 2; pad.position.set(0.004, T + 0.0001, s * 0.0065);
    const groove = mesh(new THREE.PlaneGeometry(0.022, 0.0012), slot, false);
    groove.rotation.x = -Math.PI / 2; groove.position.set(0.004, T + 0.00012, s * 0.0118);
    g.add(pad, groove);
  }
  const mid = mesh(new THREE.PlaneGeometry(0.022, 0.0012), slot, false);
  mid.rotation.x = -Math.PI / 2; mid.position.set(0.004, T + 0.00012, 0);
  g.add(mid);
  const cover = mesh(roundedBoxGeometry(0.022, 0.0006, 0.026, 0.0005, 0.0002), MAT.glass(), false);
  cover.position.set(0.004, T + 0.0002, 0);
  g.add(cover);
  const wet = mesh(new THREE.PlaneGeometry(0.0115, 0.0085), new THREE.MeshPhysicalMaterial({ color: 0x8fb4e8, transparent: true, opacity: 0.5, roughness: 0.05, depthWrite: false }), false);
  wet.rotation.x = -Math.PI / 2; wet.position.set(0.004, T + 0.00015, -0.0065); wet.visible = false;
  g.add(wet);
  const lab = mesh(new THREE.PlaneGeometry(0.02, 0.024), new THREE.MeshStandardMaterial({ map: labelTexture(["Neubauer", "0.1 mm"], { w: 128, h: 150, size: 24, bg: "#eef0ef" }), roughness: 0.9 }), false);
  lab.rotation.x = -Math.PI / 2; lab.position.set(-0.025, T + 0.0001, 0);
  g.add(lab);
  compact(g);
  link(g, "wet", wet);
  return g;
}

// A 1.5 mL microtube with its cap shut; small volumes sit in its cone.
export function makeMicrotube(label, solution = null, volume = 0) {
  const g = new THREE.Group();
  g.add(mesh(vessel([[0.0001, 0], [0.0009, 0.0006], [0.0034, 0.012], [0.0046, 0.019], [0.0047, 0.037], [0.0053, 0.038]], 0.0005, 24), MAT.plasticClear(), false));
  const cap = mesh(lathe([[0.0001, 0.0405], [0.0057, 0.0405], [0.0058, 0.0385], [0.0047, 0.0385], [0.0047, 0.034], [0.0001, 0.034]], 24), new THREE.MeshPhysicalMaterial({ color: 0xf6f8f9, roughness: 0.2, transparent: true, opacity: 0.6 }), false);
  const hinge = mesh(new THREE.BoxGeometry(0.004, 0.0012, 0.004), MAT.plasticClear(), false);
  hinge.position.set(0, 0.0385, 0.0065);
  g.add(cap, hinge);
  const top = mesh(new THREE.CircleGeometry(0.0045, 20), new THREE.MeshBasicMaterial({ map: labelTexture([label], { w: 96, h: 96, size: 40, bg: "#fafafa" }) }), false);
  top.rotation.x = -Math.PI / 2; top.position.y = 0.0407;
  g.add(top);
  // Liquid in the cone: a cone's volume goes with height cubed.
  const mat = new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.1, transparent: true, opacity: 0.85, depthWrite: false });
  const cone = new THREE.Mesh(new THREE.ConeGeometry(1, 1, 18, 1), mat);
  cone.rotation.x = Math.PI;
  cone.renderOrder = 3;
  g.add(cone);
  const liquid = {
    mesh: cone,
    update(container) {
      const v = container.volume;
      cone.visible = v > 1e-5;
      const h = Math.min(0.03, 0.012 * Math.cbrt(v / 0.03));
      const r = Math.min(0.0044, 0.0034 * (h / 0.012));
      cone.scale.set(r, h, r);
      cone.position.y = 0.0006 + h / 2;
      const { color, alpha } = container.color();
      mat.color.copy(color); mat.opacity = Math.min(0.95, 0.35 + alpha);
    },
  };
  const container = new Container(1.5, solution ? { [solution]: volume } : {});
  liquid.update(container);
  return { obj: g, container, liquid };
}

// A white rack for microtubes, two rows of five holes.
export const RACK_HOLES = Array.from({ length: 10 }, (_, i) => new THREE.Vector3(((i % 5) - 2) * 0.016, 0.03, (i < 5 ? -0.5 : 0.5) * 0.016));
export function makeTubeRack() {
  const g = new THREE.Group();
  g.add(mesh(roundedBoxGeometry(0.092, 0.03, 0.04, 0.003, 0.0012), MAT.plastic(0xf2f3ef)));
  const c = document.createElement("canvas");
  c.width = 368; c.height = 160;
  const x = c.getContext("2d");
  x.fillStyle = "#f2f3ef"; x.fillRect(0, 0, 368, 160);
  for (const h of RACK_HOLES) { x.fillStyle = "#2a2d30"; x.beginPath(); x.arc(184 + h.x / 0.092 * 368, 80 + h.z / 0.04 * 160, 22, 0, Math.PI * 2); x.fill(); }
  const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace;
  const top = mesh(new THREE.PlaneGeometry(0.088, 0.038), new THREE.MeshStandardMaterial({ map: t, roughness: 0.6 }), false);
  top.rotation.x = -Math.PI / 2; top.position.y = 0.0302;
  g.add(top);
  return g;
}

/* ---------- water bath and watch glass ---------- */
// A small stainless water bath held at 60 °C, with a perforated shelf just
// under the water for a watch glass.
export function makeWaterBath() {
  const g = new THREE.Group();
  const W = 0.34, H = 0.17, D = 0.24, t = 0.006;
  const steel = MAT.brushed();
  g.add(mesh(roundedBoxGeometry(W, 0.012, D, 0.012, 0.003), steel));
  for (const [w, d, x, z] of [[W, t, 0, -D / 2 + t / 2], [W, t, 0, D / 2 - t / 2], [t, D, -W / 2 + t / 2, 0], [t, D, W / 2 - t / 2, 0]]) {
    const wall = mesh(new THREE.BoxGeometry(w, H, d), steel);
    wall.position.set(x, H / 2, z);
    g.add(wall);
  }
  const water = mesh(new THREE.PlaneGeometry(W - 2 * t, D - 2 * t), new THREE.MeshPhysicalMaterial({ color: 0xa9cdd2, roughness: 0.04, transparent: true, opacity: 0.55, clearcoat: 1, depthWrite: false }), false);
  water.rotation.x = -Math.PI / 2; water.position.y = H - 0.035;
  water.renderOrder = 2;
  water.raycast = () => {};   // aim straight through the water at what's on the shelf
  const shelf = mesh(new THREE.BoxGeometry(W * 0.6, 0.003, D * 0.6), steel, false);
  shelf.position.y = H - 0.04;
  g.add(water, shelf);
  const c = document.createElement("canvas");
  c.width = 256; c.height = 96;
  const x = c.getContext("2d");
  x.fillStyle = "#12181b"; x.fillRect(0, 0, 256, 96);
  x.fillStyle = "#ff9d5c"; x.font = "600 44px 'IBM Plex Mono', monospace"; x.fillText("60.0°C", 16, 52);
  x.fillStyle = "#9fe0ff"; x.font = "500 22px 'IBM Plex Mono', monospace"; x.fillText("SET 60.0", 18, 84);
  const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
  const panel = mesh(new THREE.PlaneGeometry(0.08, 0.03), new THREE.MeshBasicMaterial({ map: tex, toneMapped: false }), false);
  panel.position.set(-0.08, H * 0.5, D / 2 + 0.0008);
  g.add(panel);
  compact(g);
  link(g, "shelfY", H - 0.0385);
  return g;
}

export function makeWatchGlass() {
  const g = new THREE.Group();
  const pts = [];
  for (let i = 0; i <= 10; i++) { const r = 0.0001 + (i / 10) * 0.038; pts.push([r, 0.0012 + 6.5 * r * r]); }
  for (let i = 10; i >= 0; i--) { const r = 0.0001 + (i / 10) * 0.038 - (i ? 0.0004 : 0); pts.push([Math.max(0.0001, r), 6.5 * r * r]); }
  g.add(mesh(lathe(pts, 40), MAT.glass()));
  const mat = new THREE.MeshPhysicalMaterial({ color: 0xffffff, roughness: 0.05, transparent: true, opacity: 0.6, depthWrite: false });
  const pool = new THREE.Mesh(new THREE.CylinderGeometry(1, 0.6, 1, 28), mat);
  pool.renderOrder = 3;
  g.add(pool);
  const liquid = {
    mesh: pool,
    update(container) {
      const v = container.volume;
      pool.visible = v > 1e-4;
      const r = Math.min(0.03, 0.012 * Math.sqrt(v));
      pool.scale.set(r, Math.max(0.0005, 0.0025 * Math.sqrt(v)), r);
      pool.position.y = 0.0014 + pool.scale.y / 2;
      const { color, alpha } = container.color();
      mat.color.copy(color); mat.opacity = Math.min(0.9, 0.3 + alpha);
    },
  };
  const container = new Container(5, {});
  liquid.update(container);
  return { obj: g, container, liquid };
}

/* ---------- onions ---------- */
// An onion sitting on a jar of water, with white roots growing down into it.
export function makeRootedOnion(onionModel) {
  const g = new THREE.Group();
  const r = 0.03, h = 0.075;
  g.add(mesh(vessel([[0.0001, 0], [r - 0.003, 0], [r, 0.003], [r, h - 0.006], [r - 0.004, h - 0.002], [r - 0.0035, h + 0.004]], 0.0018, 40), MAT.glass()));
  const liq = liquidColumn({ radius: r - 0.002, bottom: 0.0019, depth: h - 0.012 });
  g.add(liq.mesh);
  const container = new Container(150, { tap: 120 });
  liq.update(container);
  onionModel.position.y = h - 0.012;
  g.add(onionModel);
  const rootMat = new THREE.MeshStandardMaterial({ color: 0xf4f1e6, roughness: 0.55 });
  const roots = [];
  for (let i = 0; i < 14; i++) {
    const a = i * 2.4, rr = 0.004 + (i % 4) * 0.002;
    const top = new THREE.Vector3(Math.cos(a) * rr, h - 0.008, Math.sin(a) * rr);
    const len = 0.035 + ((i * 37) % 11) * 0.0025;
    const pts = [top];
    for (let k = 1; k <= 4; k++) pts.push(new THREE.Vector3(top.x + Math.cos(a) * 0.003 * k + Math.sin(i + k) * 0.0015, top.y - (len * k) / 4, top.z + Math.sin(a) * 0.003 * k + Math.cos(i * 2 + k) * 0.0015));
    const tipGeo = new THREE.SphereGeometry(0.001, 8, 6).translate(pts[4].x, pts[4].y, pts[4].z);
    const root = mesh(mergeGeometries([new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), 10, 0.0007, 5), tipGeo]), rootMat, false);
    g.add(root);
    roots.push(root);
  }
  link(g, "roots", roots);
  return { obj: g, container, liquid: liq };
}

// The same onion model, tinted to the purple-red of a red onion.
export function makeRedOnion(onionModel) {
  onionModel.traverse((m) => {
    if (!m.isMesh) return;
    m.material = m.material.clone();
    m.material.color = new THREE.Color(0xb1446e);
  });
  return onionModel;
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
  return compact(g);
}

/* ---------- tools ---------- */
// A glazed white ceramic tile with rounded edges.
export function makeTile() {
  const g = new THREE.Group();
  g.add(mesh(roundedBoxGeometry(0.15, 0.006, 0.15, 0.004, 0.0018), new THREE.MeshPhysicalMaterial({ color: 0xf7f7f4, roughness: 0.12, clearcoat: 0.9, clearcoatRoughness: 0.08 })));
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
  const stickMat = MAT.paper(0xf3efe4), tipMat = MAT.paper(0xffffff);
  for (let i = 0; i < 6; i++) {
    const s = makeSwab();
    s.rotation.z = Math.PI / 2 - 0.08 * (i - 2.5);
    s.rotation.y = i * 1.05;
    s.position.set(Math.cos(i) * 0.008, 0.07, Math.sin(i) * 0.008);
    s.updateMatrix();
    for (const m of [...s.children]) { m.applyMatrix4(s.matrix); m.material = m === s.userData.tip ? tipMat : stickMat; g.add(m); }
  }
  const lab = mesh(new THREE.PlaneGeometry(0.05, 0.03), new THREE.MeshStandardMaterial({ map: labelTexture(["Sterile swabs"], { w: 256, h: 140, size: 30 }), roughness: 0.9 }), false);
  lab.position.set(0, 0.045, 0.0305); lab.rotation.x = -0.08;
  g.add(lab);
  return compact(g);
}

// Fine-pointed steel forceps: two tapered arms joined at the back, with
// ridged grips, their tips just touching.
export function makeForceps() {
  const g = new THREE.Group();
  const L = 0.115;
  const arm = new THREE.Shape();
  arm.moveTo(0, -0.0024); arm.lineTo(L * 0.62, -0.0021); arm.quadraticCurveTo(L * 0.9, -0.0012, L, -0.00035);
  arm.lineTo(L, 0.00035); arm.quadraticCurveTo(L * 0.9, 0.0012, L * 0.62, 0.0021); arm.lineTo(0, 0.0024); arm.lineTo(0, -0.0024);
  const geo = new THREE.ExtrudeGeometry(arm, { depth: 0.0009, bevelEnabled: true, bevelThickness: 0.0002, bevelSize: 0.0002, bevelSegments: 1, curveSegments: 8 });
  geo.rotateX(-Math.PI / 2);
  geo.translate(-L / 2, 0, 0);
  const steel = MAT.steel();
  for (const s of [-1, 1]) {
    const a = mesh(geo, steel);
    a.position.set(0, s > 0 ? 0.0034 : 0.0011, 0);
    a.rotation.z = s * 0.028;
    g.add(a);
    for (let i = 0; i < 6; i++) {
      const ridge = mesh(new THREE.BoxGeometry(0.0009, 0.00025, 0.0042), steel, false);
      ridge.position.set(-0.012 + i * 0.0026, (s > 0 ? 0.0034 : 0.0011) + (s > 0 ? 0.0012 : -0.0001), 0);
      g.add(ridge);
    }
  }
  const back = mesh(new THREE.BoxGeometry(0.006, 0.0035, 0.0048), steel);
  back.position.set(-L / 2 + 0.003, 0.0026, 0);
  g.add(back);
  return compact(g);
}

// A scalpel: flat steel handle with a gripped section and a curved blade.
export function makeScalpel() {
  const g = new THREE.Group();
  const steel = MAT.steel();
  const handle = mesh(roundedBoxGeometry(0.105, 0.0028, 0.0085, 0.002, 0.0008), steel);
  g.add(handle);
  const grip = MAT.brushed();
  for (let i = 0; i < 9; i++) {
    const ridge = mesh(new THREE.BoxGeometry(0.0012, 0.0004, 0.0074), grip, false);
    ridge.position.set(-0.025 + i * 0.0035, 0.003, 0);
    g.add(ridge);
  }
  const b = new THREE.Shape();
  b.moveTo(0, -0.0028); b.lineTo(0.018, -0.0028); b.quadraticCurveTo(0.03, -0.0026, 0.036, 0.0032);
  b.quadraticCurveTo(0.024, 0.0038, 0.012, 0.0034); b.lineTo(0, 0.0034); b.lineTo(0, -0.0028);
  const bg = new THREE.ExtrudeGeometry(b, { depth: 0.0004, bevelEnabled: false, curveSegments: 10 });
  bg.rotateX(-Math.PI / 2);
  const blade = mesh(bg, new THREE.MeshStandardMaterial({ color: 0xe6eaee, metalness: 1, roughness: 0.12 }));
  blade.position.set(0.05, 0.0012, 0);
  g.add(blade);
  return compact(g);
}

// A sharpened yellow pencil with a ferrule and an eraser.
export function makePencil() {
  const g = new THREE.Group();
  const body = mesh(new THREE.CylinderGeometry(0.0036, 0.0036, 0.15, 6), MAT.plastic(0xf2c230));
  body.rotation.z = Math.PI / 2; body.position.y = 0.0036;
  const wood = mesh(new THREE.ConeGeometry(0.0036, 0.016, 6), MAT.plastic(0xe7c49a));
  wood.rotation.z = -Math.PI / 2; wood.position.set(0.083, 0.0036, 0);
  const lead = mesh(new THREE.ConeGeometry(0.0011, 0.005, 8), MAT.plastic(0x2a2a2c));
  lead.rotation.z = -Math.PI / 2; lead.position.set(0.0889, 0.0036, 0);
  const ferrule = mesh(new THREE.CylinderGeometry(0.0038, 0.0038, 0.008, 16), MAT.brushed());
  ferrule.rotation.z = Math.PI / 2; ferrule.position.set(-0.079, 0.0036, 0);
  const eraser = mesh(new THREE.CylinderGeometry(0.0035, 0.0035, 0.008, 16), MAT.plastic(0xe58a8a));
  eraser.rotation.z = Math.PI / 2; eraser.position.set(-0.087, 0.0036, 0);
  g.add(body, wood, lead, ferrule, eraser);
  return g;
}

// A stack of folded paper towels, each sheet a little out of line.
export function makePaperTowels() {
  const g = new THREE.Group();
  const paper = MAT.paper(0xf2efe6);
  for (let i = 0; i < 8; i++) {
    const s = mesh(roundedBoxGeometry(0.21, 0.0034, 0.11, 0.004, 0.0012), paper, i === 7);
    s.position.set(Math.sin(i * 2.1) * 0.0025, i * 0.0036, Math.cos(i * 1.7) * 0.002);
    s.rotation.y = Math.sin(i * 3.3) * 0.02;
    g.add(s);
  }
  const fold = mesh(new THREE.BoxGeometry(0.0015, 0.0005, 0.108), MAT.paper(0xe4e0d6), false);
  fold.position.set(0.02, 0.0292, 0);
  g.add(fold);
  return compact(g);
}

export function makeGlassBin() {
  const g = new THREE.Group();
  const body = mesh(roundedBoxGeometry(0.3, 0.36, 0.3, 0.02, 0.006), MAT.plastic(0xe9eaea));
  const lid = mesh(roundedBoxGeometry(0.31, 0.022, 0.31, 0.022, 0.006), MAT.plastic(0xf0c419));
  lid.position.y = 0.355;
  const lab = mesh(new THREE.PlaneGeometry(0.22, 0.12), new THREE.MeshStandardMaterial({ map: labelTexture(["Glass waste", "Slides and cover slips only"], { w: 320, h: 170, band: "#f0c419", size: 36 }), roughness: 0.9 }), false);
  lab.position.set(0, 0.22, 0.1506);
  g.add(body, lid, lab);
  return g;
}

// Wraparound safety goggles, as they rest on the bench with the strap behind.
export function makeGoggles() {
  const g = new THREE.Group();
  const r = 0.085, h = 0.05, th = 1.15;   // lens radius, height, half the wrap angle
  const lens = mesh(new THREE.CylinderGeometry(r, r, h, 40, 1, true, -th, th * 2), new THREE.MeshPhysicalMaterial({ color: 0xcfe9f2, roughness: 0.05, transparent: true, opacity: 0.45, clearcoat: 1, depthWrite: false, side: THREE.DoubleSide, forceSinglePass: true }), false);
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
  return compact(g);
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
    disinfectant: "Disinfectant", medium: "Medium", paclitaxel: "Paclitaxel", doxorubicin: "Doxorubicin", cisplatin: "Cisplatin",
    salt: "Salt solution 1 M", hcl: "Hydrochloric acid 1 M", toluidine: "Toluidine blue", trypsin: "Trypsin-EDTA", trypan: "Trypan blue 0.4%",
    cells: "Cell suspension" }[solution] || solution;
}
