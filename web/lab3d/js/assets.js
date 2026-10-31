// Loads every model, texture and the lighting map, reporting progress.
import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js";
import { OBJLoader } from "three/addons/loaders/OBJLoader.js";
import { RGBELoader } from "three/addons/loaders/RGBELoader.js";

// Poly Haven models (CC0), glTF at 1k textures.
const GLTF = ["yellow_onion", "metal_stool_01", "steel_frame_shelves_01", "mounted_fluorescent_lights", "drawer_cabinet",
  "clipboard", "stationery_supplies", "wall_clock", "magnifying_glass_01", "binder_notebook", "metal_office_desk"];

// Sketchfab models (CC Attribution), OBJ without materials: colours are set in room.js.
const OBJ = {
  tableTop: "sketchfab/lab-table/countertop.obj",
  tableCabinets: "sketchfab/lab-table/cabinets.obj",
  tableFixtures: "sketchfab/lab-table/fixtures.obj",
  tableBase: "sketchfab/lab-table/base.obj",
  microscope: "sketchfab/microscope/microscope.obj",
  incubator: "sketchfab/incubator/incubator.obj",
  doorLeaves: "sketchfab/door/leaves.obj",
  doorFrame: "sketchfab/door/frame.obj",
  doorGlass: "sketchfab/door/glass.obj",
};

const TEX = {
  floorDiff: "textures/floor_tiles_06/floor_tiles_06_diff_1k.jpg",
  floorNor: "textures/floor_tiles_06/floor_tiles_06_nor_gl_1k.jpg",
  floorRough: "textures/floor_tiles_06/floor_tiles_06_rough_1k.jpg",
  wallDiff: "textures/painted_plaster_wall/painted_plaster_wall_diff_1k.jpg",
  wallNor: "textures/painted_plaster_wall/painted_plaster_wall_nor_gl_1k.jpg",
  doorLeafColor: "models/sketchfab/door/01_01_BaseColor.png",
  doorLeafRough: "models/sketchfab/door/01_01_Roughness.png",
  doorFrameColor: "models/sketchfab/door/01_02_BaseColor.png",
};

export async function loadAssets(onProgress, { env = true } = {}) {
  const manager = new THREE.LoadingManager();
  manager.onProgress = (url, done, total) => onProgress(done / total, url.split("/").pop());
  const base = "assets/";
  const gltf = new GLTFLoader(manager), obj = new OBJLoader(manager), tex = new THREE.TextureLoader(manager);
  const out = { gltf: {}, obj: {}, tex: {}, hdr: null };
  const jobs = [];
  for (const id of GLTF) jobs.push(gltf.loadAsync(`${base}models/${id}/${id}_1k.gltf`).then((g) => (out.gltf[id] = g.scene)));
  for (const [k, p] of Object.entries(OBJ)) jobs.push(obj.loadAsync(base + "models/" + p).then((o) => (out.obj[k] = o)));
  for (const [k, p] of Object.entries(TEX)) jobs.push(tex.loadAsync(base + p).then((t) => (out.tex[k] = t)));
  if (env) jobs.push(new RGBELoader(manager).loadAsync(base + "hdri/hospital_room_1k.hdr").then((h) => (out.hdr = h)));
  await Promise.all(jobs);
  for (const k of ["floorDiff", "wallDiff", "doorLeafColor", "doorFrameColor"]) out.tex[k].colorSpace = THREE.SRGBColorSpace;
  return out;
}

/* Wrap a model so its bottom sits at y = 0 and it is centred in x/z, scaled so
   its height (or largest side) matches `size` metres, after an optional turn. */
export function fit(object, { height, maxDim, rotX = 0, rotY = 0, rotZ = 0 } = {}) {
  const inner = new THREE.Group();
  inner.add(object);
  inner.rotation.set(rotX, rotY, rotZ);
  inner.updateMatrixWorld(true);
  let box = new THREE.Box3().setFromObject(inner);
  const size = box.getSize(new THREE.Vector3());
  const s = height ? height / size.y : maxDim ? maxDim / Math.max(size.x, size.y, size.z) : 1;
  inner.scale.setScalar(s);
  inner.updateMatrixWorld(true);
  box = new THREE.Box3().setFromObject(inner);
  const c = box.getCenter(new THREE.Vector3());
  inner.position.set(-c.x, -box.min.y, -c.z);
  const outer = new THREE.Group();
  outer.add(inner);
  return outer;
}

export function setMaterial(object, material, { shadows = true } = {}) {
  object.traverse((m) => {
    if (!m.isMesh) return;
    m.material = material;
    if (!m.geometry.attributes.normal) m.geometry.computeVertexNormals();
    m.castShadow = shadows;
    m.receiveShadow = true;
  });
  return object;
}

// Store a reference to another object in userData without it being copied:
// three.js clones userData through JSON, which can't follow these links.
export function link(object, key, value) {
  Object.defineProperty(object.userData, key, { value, writable: true, configurable: true, enumerable: false });
}

export function enableShadows(object) {
  object.traverse((m) => { if (m.isMesh) { m.castShadow = true; m.receiveShadow = true; } });
  return object;
}
