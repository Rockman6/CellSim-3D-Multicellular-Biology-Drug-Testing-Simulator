// The lab room: walls with real window openings, floor, ceiling, lights and
// furniture. Furniture that does something is tagged with userData.target.
import * as THREE from "three";
import { fit, setMaterial, enableShadows, link } from "./assets.js";
import { makeWhiteboard, makeLightSwitch, makeGlassBin, makeBlinds, MAT } from "./props.js";

export const ROOM = { x0: -5, x1: 5, z0: -4, z1: 4, h: 3.0 };
const WINDOWS = [-2.7, -0.3, 2.2];          // window centres along the west wall (z)
const WIN = { w: 1.5, y0: 1.25, y1: 2.55 };

function canvasTex(w, h, draw, repeat) {
  const c = document.createElement("canvas");
  c.width = w; c.height = h;
  draw(c.getContext("2d"), w, h);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  if (repeat) { t.wrapS = t.wrapT = THREE.RepeatWrapping; t.repeat.set(...repeat); }
  return t;
}

// Freshly painted plaster: the scan's surface bumps, in a clean lab off-white.
function wallMaterial(assets, w, h) {
  const nor = assets.tex.wallNor.clone();
  nor.wrapS = nor.wrapT = THREE.RepeatWrapping; nor.repeat.set(w / 2.5, h / 2.5); nor.needsUpdate = true;
  return new THREE.MeshStandardMaterial({ normalMap: nor, normalScale: new THREE.Vector2(0.6, 0.6), color: 0xe3e6df, roughness: 0.9 });
}

export function buildRoom(scene, assets, preset) {
  const R = ROOM;
  const colliders = [];      // {x0,x1,z0,z1}
  const surfaces = [];       // meshes you can put things on
  const targets = [];        // furniture you can use
  const lamps = [];
  const add = (o) => { scene.add(o); return o; };
  const collide = (o, pad = 0.02) => {
    const b = new THREE.Box3().setFromObject(o);
    colliders.push({ x0: b.min.x - pad, x1: b.max.x + pad, z0: b.min.z - pad, z1: b.max.z + pad, top: b.max.y });
    return b;
  };
  const surface = (o) => { o.traverse((m) => { if (m.isMesh) surfaces.push(m); }); };
  const target = (o, t) => { o.userData.target = t; o.traverse((m) => { if (m.isMesh) link(m, "targetRoot", o); }); targets.push(o); return o; };

  /* ---------- floor and ceiling ---------- */
  const fmap = assets.tex.floorDiff, fnor = assets.tex.floorNor, frough = assets.tex.floorRough;
  for (const t of [fmap, fnor, frough]) { t.wrapS = t.wrapT = THREE.RepeatWrapping; t.repeat.set(5, 4); t.anisotropy = preset.maxAniso; }
  const floor = new THREE.Mesh(new THREE.PlaneGeometry(10, 8), new THREE.MeshStandardMaterial({ map: fmap, normalMap: fnor, roughnessMap: frough, color: 0xdedbd4 }));
  floor.rotation.x = -Math.PI / 2;
  floor.receiveShadow = true;
  floor.name = "floor";
  add(floor); surfaces.push(floor);

  const ceilTex = canvasTex(512, 512, (g, w, h) => {
    g.fillStyle = "#eeeeea"; g.fillRect(0, 0, w, h);
    for (let i = 0; i < 9000; i++) { g.fillStyle = `rgba(120,120,110,${Math.random() * 0.12})`; g.fillRect(Math.random() * w, Math.random() * h, 1.5, 1.5); }
    g.fillStyle = "#c9cac4"; g.fillRect(0, 0, w, 8); g.fillRect(0, 0, 8, h);
  }, [10 / 0.6, 8 / 0.6]);
  const ceiling = new THREE.Mesh(new THREE.PlaneGeometry(10, 8), new THREE.MeshStandardMaterial({ map: ceilTex, roughness: 0.95 }));
  ceiling.rotation.x = Math.PI / 2;
  ceiling.position.y = R.h;
  add(ceiling);

  /* ---------- walls ---------- */
  const T = 0.12;
  const box = (w, h, d, mat, x, y, z) => { const m = new THREE.Mesh(new THREE.BoxGeometry(w, h, d), mat); m.position.set(x, y, z); m.castShadow = true; m.receiveShadow = true; return add(m); };
  box(10 + 2 * T, R.h, T, wallMaterial(assets, 10, R.h), 0, R.h / 2, R.z0 - T / 2);           // north
  box(10 + 2 * T, R.h, T, wallMaterial(assets, 10, R.h), 0, R.h / 2, R.z1 + T / 2);           // south (door is in front of it)
  box(T, R.h, 8, wallMaterial(assets, 8, R.h), R.x1 + T / 2, R.h / 2, 0);                     // east
  // West wall in pieces around the window openings, so daylight comes in.
  const wm = wallMaterial(assets, 8, R.h);
  box(T, WIN.y0, 8, wm, R.x0 - T / 2, WIN.y0 / 2, 0);
  box(T, R.h - WIN.y1, 8, wm, R.x0 - T / 2, (R.h + WIN.y1) / 2, 0);
  const edges = [R.z0, ...WINDOWS.flatMap((z) => [z - WIN.w / 2, z + WIN.w / 2]), R.z1];
  for (let i = 0; i < edges.length; i += 2) {
    const a = edges[i], b = edges[i + 1];
    if (b - a > 0.001) box(T, WIN.y1 - WIN.y0, b - a, wm, R.x0 - T / 2, (WIN.y0 + WIN.y1) / 2, (a + b) / 2);
  }
  const frameMat = new THREE.MeshStandardMaterial({ color: 0xf2f2ef, roughness: 0.4, metalness: 0.2 });
  const glassMat = new THREE.MeshPhysicalMaterial({ color: 0xdfeef2, transparent: true, opacity: 0.12, roughness: 0.02, clearcoat: 1, depthWrite: false });
  for (const z of WINDOWS) {
    const wy = (WIN.y0 + WIN.y1) / 2, wh = WIN.y1 - WIN.y0;
    box(0.06, 0.05, WIN.w, frameMat, R.x0 + 0.01, WIN.y0, z);
    box(0.06, 0.05, WIN.w, frameMat, R.x0 + 0.01, WIN.y1, z);
    box(0.06, wh, 0.05, frameMat, R.x0 + 0.01, wy, z - WIN.w / 2);
    box(0.06, wh, 0.05, frameMat, R.x0 + 0.01, wy, z + WIN.w / 2);
    box(0.05, wh, 0.04, frameMat, R.x0 + 0.01, wy, z);
    box(0.25, 0.03, WIN.w + 0.1, frameMat, R.x0 + 0.1, WIN.y0 - 0.01, z);      // sill
    const pane = new THREE.Mesh(new THREE.PlaneGeometry(WIN.w, wh), glassMat);
    pane.rotation.y = Math.PI / 2; pane.position.set(R.x0 - 0.02, wy, z);
    add(pane);
  }
  // What you see outside: sky over trees and a low building across the yard.
  const outside = canvasTex(2048, 512, (g, w, h) => {
    const sky = g.createLinearGradient(0, 0, 0, h);
    sky.addColorStop(0, "#8fb8e0"); sky.addColorStop(0.6, "#d7e7f2"); sky.addColorStop(1, "#eef3f2");
    g.fillStyle = sky; g.fillRect(0, 0, w, h);
    g.fillStyle = "rgba(255,255,255,0.55)";
    for (let i = 0; i < 26; i++) { g.beginPath(); g.ellipse(Math.random() * w, 40 + Math.random() * 150, 60 + Math.random() * 120, 14 + Math.random() * 18, 0, 0, Math.PI * 2); g.fill(); }
    g.fillStyle = "#c9c2b5"; g.fillRect(0, h * 0.62, w, h * 0.2);
    g.fillStyle = "#a9a294"; for (let x = 20; x < w; x += 90) g.fillRect(x, h * 0.66, 50, 34);
    for (let i = 0; i < 160; i++) { const x = Math.random() * w, r = 30 + Math.random() * 50; g.fillStyle = ["#4e7a3e", "#5f8b47", "#3f6a34"][i % 3]; g.beginPath(); g.arc(x, h * 0.78 + Math.random() * 30, r, 0, Math.PI * 2); g.fill(); }
    g.fillStyle = "#6f8f52"; g.fillRect(0, h * 0.9, w, h * 0.1);
  });
  const backdrop = new THREE.Mesh(new THREE.PlaneGeometry(30, 7.5), new THREE.MeshBasicMaterial({ map: outside, toneMapped: false }));
  backdrop.rotation.y = Math.PI / 2; backdrop.position.set(R.x0 - 7, 2.2, 0);
  add(backdrop);

  // Skirting boards.
  const skirt = new THREE.MeshStandardMaterial({ color: 0x3b3f44, roughness: 0.6 });
  box(10, 0.1, 0.02, skirt, 0, 0.05, R.z0 + 0.01); box(10, 0.1, 0.02, skirt, 0, 0.05, R.z1 - 0.01);
  box(0.02, 0.1, 8, skirt, R.x1 - 0.01, 0.05, 0); box(0.02, 0.1, 8, skirt, R.x0 + 0.01, 0.05, 0);

  /* ---------- the big lab bench (Sketchfab "Chemistry Lab Table") ---------- */
  const table = new THREE.Group();
  table.add(setMaterial(assets.obj.tableTop, new THREE.MeshStandardMaterial({ color: 0x1c1e21, roughness: 0.3 })));
  table.add(setMaterial(assets.obj.tableCabinets, new THREE.MeshStandardMaterial({ color: 0xe6e9e6, roughness: 0.55 })));
  table.add(setMaterial(assets.obj.tableFixtures, new THREE.MeshStandardMaterial({ color: 0xc6cbd0, metalness: 1, roughness: 0.3 })));
  table.add(setMaterial(assets.obj.tableBase, new THREE.MeshStandardMaterial({ color: 0x2c2f33, roughness: 0.7 })));
  table.position.set(-2.75, 0, -1.79);
  add(table);
  surface(assets.obj.tableTop);
  colliders.push({ x0: -5, x1: 1.16, z0: -4, z1: -1.29, top: 1.1 }, { x0: -5, x1: -3.98, z0: -1.3, z1: 1.33, top: 1.1 });
  const TABLE_Y = 1.098;
  // The sink and tap, at the bench's built-in basin.
  const sink = new THREE.Mesh(new THREE.BoxGeometry(0.5, 0.42, 0.42), new THREE.MeshBasicMaterial({ visible: false }));
  sink.position.set(-0.63, TABLE_Y + 0.1, -2.08);
  add(target(sink, { kind: "sink", name: "Sink" }));
  const water = new THREE.Mesh(new THREE.CylinderGeometry(0.006, 0.008, 0.3, 10), new THREE.MeshPhysicalMaterial({ color: 0xdff1f7, transparent: true, opacity: 0.55, roughness: 0.02 }));
  water.position.set(-0.63, TABLE_Y + 0.15, -2.16); water.visible = false;
  add(water);
  link(sink, "water", water);

  /* ---------- microscope bench, built to match ---------- */
  const BENCH = { x: -1.0, z: 0.9, w: 2.7, d: 0.85, h: 0.92 };
  const bench = new THREE.Group();
  const white = new THREE.MeshStandardMaterial({ color: 0xe6e9e6, roughness: 0.55 });
  const body = new THREE.Mesh(new THREE.BoxGeometry(BENCH.w - 0.06, BENCH.h - 0.13, BENCH.d - 0.1), white);
  body.position.y = (BENCH.h - 0.13) / 2 + 0.1;
  const top = new THREE.Mesh(new THREE.BoxGeometry(BENCH.w, 0.03, BENCH.d), new THREE.MeshStandardMaterial({ color: 0x1c1e21, roughness: 0.3 }));
  top.position.y = BENCH.h - 0.015;
  const kick = new THREE.Mesh(new THREE.BoxGeometry(BENCH.w - 0.1, 0.1, BENCH.d - 0.16), skirt);
  kick.position.y = 0.05;
  bench.add(body, top, kick);
  const handleMat = MAT.steel();
  for (const side of [-1, 1]) for (let i = 0; i < 4; i++) {
    const fx = -BENCH.w / 2 + 0.36 + i * 0.66;
    const door = new THREE.Mesh(new THREE.BoxGeometry(0.62, BENCH.h - 0.2, 0.015), white);
    door.position.set(fx, (BENCH.h - 0.2) / 2 + 0.12, side * (BENCH.d / 2 - 0.045));
    const handle = new THREE.Mesh(new THREE.BoxGeometry(0.12, 0.012, 0.018), handleMat);
    handle.position.set(fx, BENCH.h - 0.17, side * (BENCH.d / 2 - 0.03));
    bench.add(door, handle);
  }
  enableShadows(bench);
  bench.position.set(BENCH.x, 0, BENCH.z);
  add(bench);
  surface(top);
  collide(bench);

  /* ---------- microscope (Sketchfab "Compound Microscope") ---------- */
  const scopeObj = fit(setMaterial(assets.obj.microscope, new THREE.MeshStandardMaterial({ color: 0xc9ccc8, roughness: 0.5, metalness: 0.15 })), { height: 0.42, rotY: Math.PI });
  scopeObj.position.set(-0.85, BENCH.h, 0.82);
  add(scopeObj);
  target(scopeObj, { kind: "microscope", name: "Microscope" });
  // Where a slide or dish sits on the stage, in the microscope's own space.
  scopeObj.userData.stage = new THREE.Vector3(0, 0.155, 0.02);

  /* ---------- incubator (from the Sketchfab "lab bench") ---------- */
  const inc = fit(setMaterial(assets.obj.incubator, new THREE.MeshStandardMaterial({ color: 0xe9dfb7, roughness: 0.5 })), { height: 0.45 });
  inc.position.set(0.05, BENCH.h, 0.8);
  add(inc);
  // A front display and a door handle, so it reads as an incubator.
  const ib = new THREE.Box3().setFromObject(inc);
  const disp = document.createElement("canvas");
  disp.width = 256; disp.height = 96;
  const dg = disp.getContext("2d");
  dg.fillStyle = "#12181b"; dg.fillRect(0, 0, 256, 96);
  dg.fillStyle = "#7dffb0"; dg.font = "600 40px 'IBM Plex Mono', monospace"; dg.fillText("37.0°C", 18, 46);
  dg.font = "500 24px 'IBM Plex Mono', monospace"; dg.fillStyle = "#9fe0ff"; dg.fillText("CO2 5.0%", 20, 82);
  const dt = new THREE.CanvasTexture(disp); dt.colorSpace = THREE.SRGBColorSpace;
  const panel = new THREE.Mesh(new THREE.PlaneGeometry(0.11, 0.042), new THREE.MeshBasicMaterial({ map: dt, toneMapped: false }));
  panel.position.set((ib.min.x + ib.max.x) / 2 - 0.05, ib.max.y - 0.06, ib.max.z + 0.002);
  const ih = new THREE.Mesh(new THREE.BoxGeometry(0.014, 0.12, 0.02), MAT.steel());
  ih.position.set(ib.max.x - 0.04, (ib.min.y + ib.max.y) / 2, ib.max.z + 0.01);
  add(panel); add(ih);
  target(inc, { kind: "incubator", name: "Incubator" });

  /* ---------- teacher's desk, journal and whiteboard ---------- */
  const desk = enableShadows(fit(assets.gltf.metal_office_desk, { height: 0.79, rotY: -Math.PI / 2 }));
  desk.position.set(3.95, 0, -2.3);
  add(desk); surface(desk);
  const deskBox = collide(desk);
  const board = makeWhiteboard([
    "Today: cells under the microscope",
    "1. Goggles on, hands washed",
    "2. Clean slide on the white tile",
    "3. Sample: onion skin, cheek cells or a pondweed leaf",
    "4. One drop: iodine (onion), methylene blue (cheek), water (leaf)",
    "5. Cover slip at 45°, lower it slowly: no bubbles!",
    "6. Clip it on the stage. Start at 4×, then 10×, then 40×",
    "Living cells: A549 dishes are in the incubator (37 °C).",
    "Cytotoxic drugs: gloves, and only drops into the dish.",
  ]);
  board.rotation.y = -Math.PI / 2;
  board.position.set(R.x1 - 0.02, 1.6, 0.4);
  add(board);
  target(board, { kind: "whiteboard", name: "Whiteboard" });

  /* ---------- storage ---------- */
  const shelves = enableShadows(fit(assets.gltf.steel_frame_shelves_01, { height: 2.14 }));
  shelves.position.set(-3.6, 0, 3.7);
  add(shelves); surface(shelves); collide(shelves);
  const cabinet = enableShadows(fit(assets.gltf.drawer_cabinet, { height: 1.88, rotY: Math.PI }));
  cabinet.position.set(-1.4, 0, 3.7);
  add(cabinet); surface(cabinet); collide(cabinet);
  const drawers = [];
  cabinet.traverse((n) => { if (/drawer_0\d$/.test(n.name)) drawers.push(n); });
  drawers.forEach((d, i) => { d.userData.closedZ = d.position.z; d.userData.open = false; target(d, { kind: "drawer", name: `Drawer ${i + 1}`, index: i }); });

  /* ---------- stools ---------- */
  const stools = [[-1.15, 1.75], [0.0, 1.75], [-2.3, -0.9], [-3.4, -0.9]].map(([x, z], i) => {
    const s = enableShadows(fit(i === 0 ? assets.gltf.metal_stool_01 : assets.gltf.metal_stool_01.clone(), { height: 0.88 }));
    s.position.set(x, 0, z);
    add(s);
    target(s, { kind: "stool", name: "Lab stool" });
    return s;
  });

  /* ---------- the door, the light switch, the clock ---------- */
  const door = new THREE.Group();
  const leafMat = new THREE.MeshStandardMaterial({ map: assets.tex.doorLeafColor, roughnessMap: assets.tex.doorLeafRough, roughness: 1 });
  door.add(setMaterial(assets.obj.doorLeaves, leafMat));
  door.add(setMaterial(assets.obj.doorFrame, new THREE.MeshStandardMaterial({ map: assets.tex.doorFrameColor, roughness: 0.6 })));
  door.add(setMaterial(assets.obj.doorGlass, new THREE.MeshPhysicalMaterial({ color: 0xdfeef2, transparent: true, opacity: 0.25, roughness: 0.05 }), { shadows: false }));
  const doorObj = fit(door, { height: 2.2, rotY: Math.PI });
  doorObj.position.set(3.0, 0, R.z1 - 0.07);
  add(doorObj);
  target(doorObj, { kind: "door", name: "Door" });

  const sw = makeLightSwitch();
  sw.rotation.y = Math.PI;
  sw.position.set(1.7, 1.25, R.z1 - 0.005);
  add(sw);
  target(sw, { kind: "switch", name: "Light switch" });

  const clock = enableShadows(fit(assets.gltf.wall_clock, { height: 0.32 }));
  clock.position.set(-2.2, 2.25, R.z0 + 0.03);
  add(clock);
  target(clock, { kind: "clock", name: "Wall clock" });
  const hands = {};
  clock.traverse((n) => { if (/minute_hand/.test(n.name)) hands.m = n; else if (/second_hand/.test(n.name)) hands.s = n; else if (/hours_hand/.test(n.name)) hands.h = n; });

  /* ---------- blinds on each window ---------- */
  const blinds = WINDOWS.map((z) => {
    const b = makeBlinds({ width: WIN.w + 0.06, height: WIN.y1 - WIN.y0 + 0.1 });
    b.position.set(R.x0 + 0.08, WIN.y0 - 0.05, z);
    b.userData.down = true;
    add(b);
    target(b, { kind: "blinds", name: "Window blinds" });
    return b;
  });

  /* ---------- ceiling lights ---------- */
  const fixturePos = [[-3, -2.2], [0, -2.2], [3, -2.2], [-3, 1.6], [0, 1.6], [3, 1.6]];
  const fixtures = fixturePos.map(([x, z], i) => {
    const f = fit(i === 0 ? assets.gltf.mounted_fluorescent_lights : assets.gltf.mounted_fluorescent_lights.clone(), { maxDim: 1.2 });
    f.position.set(x, R.h - 0.08, z);
    add(f);
    return f;
  });
  const tubeMat = new THREE.MeshBasicMaterial({ color: 0xf4f7ff, toneMapped: false });
  for (const [x, z] of fixturePos) {
    const glow = new THREE.Mesh(new THREE.PlaneGeometry(1.1, 0.5), tubeMat);
    glow.rotation.x = Math.PI / 2; glow.position.set(x, R.h - 0.09, z);
    add(glow);
  }
  const lampOrder = [1, 4, 0, 5, 2, 3];
  // All six exist; the graphics preset decides how many are switched on.
  for (let i = 0; i < 6; i++) {
    const [x, z] = fixturePos[lampOrder[i]];
    const l = new THREE.PointLight(0xf3f6ff, 7, 0, 2);
    l.position.set(x, R.h - 0.35, z);
    add(l);
    lamps.push(l);
  }

  /* ---------- glass waste ---------- */
  const bin = makeGlassBin();
  bin.position.set(-3.6, 0, 1.75);
  add(bin); collide(bin);
  target(bin, { kind: "glassbin", name: "Glass waste bin" });

  /* ---------- daylight and fill ---------- */
  const hemi = new THREE.HemisphereLight(0xf2f5ff, 0x6a675f, preset.env ? 0.35 : 1.25);
  add(hemi);
  const sun = new THREE.DirectionalLight(0xfff1df, 3.2);
  sun.position.set(-9, 5.5, -1.5);
  sun.target.position.set(0, 0, 0.5);
  add(sun); add(sun.target);
  if (preset.shadows) {
    sun.castShadow = true;
    sun.shadow.mapSize.set(preset.shadowSize, preset.shadowSize);
    Object.assign(sun.shadow.camera, { left: -7, right: 7, top: 6, bottom: -6, near: 1, far: 22 });
    sun.shadow.bias = -0.0004;
    sun.shadow.normalBias = 0.02;
  }

  // Walls stop you walking out.
  colliders.push(
    { x0: -99, x1: R.x0 + 0.05, z0: -99, z1: 99 }, { x0: R.x1 - 0.05, x1: 99, z0: -99, z1: 99 },
    { x0: -99, x1: 99, z0: -99, z1: R.z0 + 0.05 }, { x0: -99, x1: 99, z0: R.z1 - 0.05, z1: 99 },
  );

  return {
    colliders, surfaces, targets, lamps, sun, hemi, sink, water, scope: scopeObj, incubator: inc, clockHands: hands,
    drawers, blinds, door: doorObj, lightSwitch: sw, stools, fixtures, tubeMat,
    TABLE_Y, BENCH, deskTop: deskBox.max.y, shelves, cabinet,
  };
}
