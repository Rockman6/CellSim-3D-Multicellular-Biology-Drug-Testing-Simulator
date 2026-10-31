// CellSim Lab: starts the renderer, loads the lab, and runs the game.
import * as THREE from "three";
import { loadAssets } from "./assets.js";
import { buildRoom } from "./room.js";
import { Player } from "./player.js";
import { createItems, getActions, Item } from "./items.js";
import { HUD } from "./hud.js";
import { Tasks, EXPERIMENTS } from "./tasks.js";
import { Journal } from "./journal.js";
import { createScope } from "./scope.js";
import { Engine } from "./engine.js";
import { settings, saveSettings, detectPreset, activePreset } from "./settings.js";
import * as P from "./props.js";

const $ = (id) => document.getElementById(id);
const REACH = 2.6;
const TIPS = [
  "Left-click takes things into your pocket. Right-click uses them.",
  "Lower a cover slip slowly by holding right-click: a flat drop traps air bubbles.",
  "Start every slide with the 4× lens, then go up one lens at a time.",
  "The living cells in the incubator are simulated by the Cell Sim engine, cell by cell.",
  "Dish B is your control: treat Dish A, keep B untreated, and compare.",
  "Low graphics saves battery: no shadows, 30 frames a second, and no drawing when nothing moves.",
];

class Inventory {
  constructor() { this.slots = Array(36).fill(null); }
  add(it, prefer) {
    let i = prefer != null && !this.slots[prefer] ? prefer : this.slots.findIndex((s, k) => !s && k < 9);
    if (i < 0) i = this.slots.indexOf(null);
    if (i < 0) return -1;
    this.slots[i] = it;
    return i;
  }
  remove(it) { const i = this.slots.indexOf(it); if (i >= 0) this.slots[i] = null; return i; }
}

class Game {
  async boot() {
    if (location.protocol === "file:") { $("loading").hidden = true; $("fileWarn").hidden = false; return; }
    const probe = document.createElement("canvas").getContext("webgl2") || document.createElement("canvas").getContext("webgl");
    if (!probe) { $("loading").hidden = true; $("noGL").hidden = false; return; }
    this.detected = detectPreset(probe);
    this.preset = activePreset(this.detected.pick);
    this.canvas = $("view");
    this.renderer = new THREE.WebGLRenderer({ canvas: this.canvas, antialias: this.preset.key !== "low", powerPreference: this.preset.key === "low" ? "low-power" : "high-performance" });
    const r = this.renderer;
    r.toneMapping = THREE.ACESFilmicToneMapping;
    r.toneMappingExposure = 1.0;
    r.shadowMap.type = THREE.PCFSoftShadowMap;
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x222629);
    this.camera = new THREE.PerspectiveCamera(settings.fov, 1, 0.03, 60);
    this.scene.add(this.camera);
    this.hand = new THREE.Group();
    this.camera.add(this.hand);
    this.resize();
    addEventListener("resize", () => this.resize());

    let tip = 0;
    const tipT = setInterval(() => ($("loadTip").textContent = TIPS[tip++ % TIPS.length]), 2600);
    $("loadTip").textContent = TIPS[tip++];
    const assets = await loadAssets((f, name) => { $("loadBar").style.width = `${(f * 100).toFixed(0)}%`; $("loadText").textContent = `Loading ${name}`; }, { env: true });
    clearInterval(tipT);
    this.assets = assets;
    $("loadText").textContent = "Setting up the lab…";
    const pm = new THREE.PMREMGenerator(r);
    this.envMap = pm.fromEquirectangular(assets.hdr).texture;
    assets.hdr.dispose(); pm.dispose();

    this.room = buildRoom(this.scene, assets, this.preset);
    this.player = new Player(this.camera);
    this.inv = new Inventory();
    this.selected = 0;
    this.tasks = new Tasks((key) => { HUD.tasks(this.tasks.view()); if (key) HUD.toast("Step done ✓", "good"); });
    this.journal = new Journal();
    this.engine = new Engine((state, text) => HUD.engine(state, text));
    const made = createItems(this, assets, this.room);
    this.items = made.items;
    this.incubator = { dishes: [...made.dishes] };
    this.dishes = made.dishes;
    this.stageItem = null;
    this.lightsOn = true;
    this.anims = [];
    this.iconCache = new Map();
    this.scope = createScope({
      resolution: this.preset.key === "low" ? 560 : this.preset.key === "ultra" ? 1024 : 768,
      onFocus: (obj, spec) => this.onFocus(obj, spec),
      onExit: () => this.exitScope(),
      toast: (t, k) => HUD.toast(t, k),
      observe: (s, p) => this.journal.observe(s, p),
    });
    for (const b of this.room.blinds) this.setBlinds(b, false, true);
    this.applyPreset();
    this.refreshPickRoots();
    this.bindUI();
    this.bindInput();
    this.screen = "title";
    $("loading").hidden = true;
    const exp = new URLSearchParams(location.search).get("exp");
    if (exp && (EXPERIMENTS[exp] || exp === "free")) this.startGame(exp === "free" ? null : exp, { lock: false });
    else this.showTitle();
    this.last = performance.now();
    requestAnimationFrame((t) => this.frame(t));
    window.lab = this;   // for automated checks and curious people with the console open
  }

  /* ---------- graphics ---------- */
  resize() {
    const w = innerWidth, h = innerHeight;
    this.renderer.setPixelRatio(Math.min(devicePixelRatio || 1, this.preset.pixelRatio));
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.dirty = true;
  }
  applyPreset() {
    const p = this.preset, r = this.renderer, room = this.room;
    r.setPixelRatio(Math.min(devicePixelRatio || 1, p.pixelRatio));
    r.shadowMap.enabled = p.shadows;
    room.sun.castShadow = p.shadows;
    if (p.shadows) {
      if (room.sun.shadow.mapSize.x !== p.shadowSize) { room.sun.shadow.mapSize.set(p.shadowSize, p.shadowSize); room.sun.shadow.map?.dispose(); room.sun.shadow.map = null; }
      Object.assign(room.sun.shadow.camera, { left: -7, right: 7, top: 6, bottom: -6, near: 1, far: 22 });
      room.sun.shadow.camera.updateProjectionMatrix();
      room.sun.shadow.bias = -0.0004; room.sun.shadow.normalBias = 0.02;
    }
    this.scene.environment = p.env ? this.envMap : null;
    room.hemi.intensity = p.env ? 0.35 : 1.25;
    room.lamps.forEach((l, i) => (l.visible = this.lightsOn && i < p.lamps));
    this.scene.traverse((o) => { if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) => (m.needsUpdate = true)); });
    this.resize();
    $("qualityNote").textContent = `Using ${p.label} graphics${settings.quality === "auto" ? ` (picked for ${this.detected.gpu || "this computer"})` : ""}. ${p.fps ? `Capped at ${p.fps} fps.` : "No frame cap."}`;
  }

  /* ---------- screens ---------- */
  showTitle() {
    this.screen = "title";
    $("hud").hidden = true;
    $("title").hidden = false;
    if (matchMedia("(pointer: coarse)").matches && !matchMedia("(any-pointer: fine)").matches && !$("touchNote"))
      $("mainMenu").insertAdjacentHTML("afterbegin", `<p id="touchNote" class="note">CellSim Lab needs a keyboard and mouse. On this device you can look around, but not play.</p>`);
    $("mainMenu").hidden = false;
    $("experimentMenu").hidden = true;
    document.exitPointerLock?.();
  }
  startGame(id, { lock = true } = {}) {
    if (this.started) { location.search = `?exp=${id || "free"}`; return; }   // a fresh lab for a new experiment
    this.started = true;
    this.tasks.start(id);
    this.screen = "play";
    $("title").hidden = true;
    $("hud").hidden = false;
    this.player.pos.set(3, 0, 2.9);
    this.player.yaw = 0.8;
    this.player.pitch = -0.12;
    HUD.tasks(this.tasks.view());
    this.updateHotbar();
    if (id === "culture") this.engine.start().catch(() => {});
    if (lock) this.lock();
    HUD.toast("WASD to walk, mouse to look. Left-click takes, right-click uses.");
  }
  lock() { if (this.screen === "play" && !this.ui) this.canvas.requestPointerLock?.()?.catch?.(() => {}); }
  openUI(id) {
    if (this.ui) $(this.ui).hidden = true;
    this.uiStack = this.uiStack || [];
    if (this.ui) this.uiStack.push(this.ui);
    this.ui = id;
    $(id).hidden = false;
    this.player.keys.clear();
    if (document.pointerLockElement) { this.unlockingForUI = true; document.exitPointerLock(); }
  }
  closeUI() {
    if (!this.ui) return;
    $(this.ui).hidden = true;
    this.ui = this.uiStack && this.uiStack.length ? this.uiStack.pop() : null;
    if (this.ui) $(this.ui).hidden = false;
    else if (this.screen === "title") $("title").hidden = false;
    this.dirty = true;
  }
  closeAllUI() { while (this.ui) this.closeUI(); }

  bindUI() {
    // Title menu
    $("mainMenu").addEventListener("click", (e) => {
      const b = e.target.closest("[data-menu]");
      if (!b) return;
      const m = b.dataset.menu;
      if (m === "experiments") { $("mainMenu").hidden = true; this.renderExperiments(); $("experimentMenu").hidden = false; }
      else if (m === "free") this.startGame(null);
      else this.openUI(m);
    });
    $("experimentMenu").addEventListener("click", (e) => {
      const b = e.target.closest("[data-exp]");
      if (!b) return;
      if (b.dataset.exp === "back") { $("experimentMenu").hidden = true; $("mainMenu").hidden = false; return; }
      this.startGame(b.dataset.exp);
    });
    $("pause").addEventListener("click", (e) => {
      const b = e.target.closest("[data-pause]");
      if (!b) return;
      const a = b.dataset.pause;
      if (a === "resume") { this.closeAllUI(); this.lock(); }
      else if (a === "journal") this.openJournal();
      else if (a === "title") { this.closeAllUI(); this.leave(); }
      else this.openUI(a);
    });
    for (const b of document.querySelectorAll("[data-back]")) b.addEventListener("click", () => { this.closeUI(); if (!this.ui && this.screen === "play") this.lock(); });
    // Settings
    const q = $("setQuality"), fps = $("setFps"), sens = $("setSens"), fov = $("setFov"), inv = $("setInvert"), bob = $("setBob");
    q.value = settings.quality; fps.value = settings.fps; sens.value = settings.sensitivity; fov.value = settings.fov; inv.checked = settings.invertY; bob.checked = settings.bob;
    const onSet = () => {
      Object.assign(settings, { quality: q.value, fps: fps.value, sensitivity: +sens.value, fov: +fov.value, invertY: inv.checked, bob: bob.checked });
      saveSettings();
      this.preset = activePreset(this.detected.pick);
      this.applyPreset();
    };
    for (const el of [q, fps, sens, fov, inv, bob]) el.addEventListener("change", onSet);
    // Credits
    $("creditsBody").innerHTML = CREDITS;
    // Journal tabs
    $("journal").addEventListener("click", (e) => {
      const b = e.target.closest("[data-tab]");
      if (!b) return;
      this.journalTab = b.dataset.tab;
      this.renderJournal();
    });
    // Inventory: click a slot, then another, to swap them.
    const invClick = (e) => {
      const s = e.target.closest("[data-slot]");
      if (!s) return;
      const i = +s.dataset.slot;
      if (this.invPick == null) { if (this.inv.slots[i]) this.invPick = i; }
      else { const a = this.inv.slots; [a[i], a[this.invPick]] = [a[this.invPick], a[i]]; this.invPick = null; this.updateHotbar(); }
      HUD.inventoryGrids(this.inv, this.invPick, (it) => this.iconFor(it));
    };
    $("invGrid").addEventListener("click", invClick);
    $("invHotbar").addEventListener("click", invClick);
    // Incubator buttons
    $("incubatorBody").addEventListener("click", (e) => {
      const b = e.target.closest("[data-inc]");
      if (!b) return;
      if (b.dataset.inc === "wait") this.incubate(+b.dataset.h);
      else if (b.dataset.inc === "take") {
        const d = this.incubator.dishes.find((x) => x.id === b.dataset.id);
        if (d && this.give(d)) this.emit("dish-taken");   // take() also removes it from the incubator
        this.renderIncubator();
      }
    });
  }

  renderExperiments() {
    $("experimentMenu").innerHTML = Object.entries(EXPERIMENTS).map(([id, e]) =>
      `<button class="exp-card" data-exp="${id}"><span class="tag">${e.tag}</span><b>${e.title}</b><span>${e.blurb}</span></button>`).join("") +
      `<button class="mbtn small" data-exp="back">Back</button>`;
  }

  bindInput() {
    const c = this.canvas;
    c.addEventListener("click", () => { if (this.screen === "play" && !document.pointerLockElement) this.lock(); });
    document.addEventListener("pointerlockchange", () => {
      this.locked = document.pointerLockElement === c;
      if (!this.locked && this.screen === "play" && !this.ui && !this.unlockingForUI) this.openUI("pause");
      this.unlockingForUI = false;
      this.dirty = true;
    });
    document.addEventListener("mousemove", (e) => {
      if (!this.locked) return;
      this.player.look(e.movementX, e.movementY);
      this.lastInput = performance.now();
    });
    c.addEventListener("contextmenu", (e) => e.preventDefault());
    document.addEventListener("mousedown", (e) => { if (this.locked) { e.preventDefault(); this.press(e.button === 2 ? "right" : e.button === 0 ? "left" : null); } });
    document.addEventListener("mouseup", (e) => { if (this.locked && e.button === 2) this.release(); });
    document.addEventListener("wheel", (e) => { if (!this.locked) return; this.select((this.selected + (e.deltaY > 0 ? 1 : 8)) % 9); }, { passive: true });
    document.addEventListener("keydown", (e) => {
      if (e.code === "F3") { e.preventDefault(); $("debug").hidden = !$("debug").hidden; return; }
      if (this.screen === "scope") { if (e.code === "Escape") this.exitScope(); return; }
      if (this.screen !== "play") return;
      if (this.ui) {
        if (e.code === "Escape" || (e.code === "KeyE" && this.ui === "inventory") || (e.code === "KeyJ" && this.ui === "journal")) { this.closeUI(); if (!this.ui) this.lock(); }
        return;
      }
      if (/^Digit[1-9]$/.test(e.code)) this.select(+e.code.slice(5) - 1);
      else if (e.code === "KeyQ") this.dropHeld();
      else if (e.code === "KeyE") this.openInventory();
      else if (e.code === "KeyJ") this.openJournal();
      else this.player.keys.add(e.code);
      this.lastInput = performance.now();
    });
    document.addEventListener("keyup", (e) => this.player.keys.delete(e.code));
    addEventListener("blur", () => this.player.keys.clear());
  }

  /* ---------- inventory ---------- */
  get held() { return this.inv.slots[this.selected]; }
  select(i) {
    if (i === this.selected) return;
    this.selected = i;
    this.updateHotbar();
    HUD.itemName(this.held ? this.held.name : "");
  }
  updateHotbar() {
    HUD.hotbar(this.inv, this.selected, (it) => this.iconFor(it));
    HUD.held(this.held);
    // Show the selected item in your hand.
    while (this.hand.children.length) this.hand.remove(this.hand.children[0]);
    const it = this.held;
    if (it) {
      if (!it.size) {
        it.obj.removeFromParent();
        it.obj.position.set(0, 0, 0); it.obj.rotation.set(0, 0, 0); it.obj.scale.setScalar(1);
        it.obj.updateMatrixWorld(true);
        it.size = new THREE.Box3().setFromObject(it.obj).getSize(new THREE.Vector3());
      }
      const size = it.size;
      const s = THREE.MathUtils.clamp(0.15 / Math.max(size.x, size.y, size.z, 0.001), 0.5, 2.4);
      it.obj.position.set(0.2, -0.2, -0.42);
      it.obj.rotation.set(0.15, -0.5, 0);
      it.obj.scale.setScalar(s);
      this.hand.add(it.obj);
    }
    this.dirty = true;
  }
  take(it) {
    if (it === this.stageItem) { this.stageItem = null; }
    const empty = !this.held;
    const i = this.inv.add(it, empty ? this.selected : null);
    if (i < 0) { HUD.toast("Your pockets are full. Put something down first (Q).", "warn"); return false; }
    it.obj.removeFromParent();
    it.obj.scale.setScalar(1);
    this.refreshPickRoots();
    this.updateHotbar();
    HUD.itemName(it.name);
    if (it.kind === "dish" && this.incubator.dishes.includes(it)) this.incubator.dishes.splice(this.incubator.dishes.indexOf(it), 1);
    return true;
  }
  give(it) { if (!this.items.includes(it)) this.items.push(it); return this.take(it); }
  place(it, point) {
    this.inv.remove(it);
    it.obj.removeFromParent();
    it.obj.scale.setScalar(1);
    it.obj.position.copy(point);
    it.obj.rotation.set(0, this.player.yaw + Math.PI / 2, 0);
    this.scene.add(it.obj);
    this.refreshPickRoots();
    this.updateHotbar();
  }
  dropHeld() {
    const it = this.held;
    if (!it) return;
    const hit = this.aim;
    if (hit && hit.surface) this.place(it, hit.point);
    else {
      const f = new THREE.Vector3(-Math.sin(this.player.yaw), 0, -Math.cos(this.player.yaw));
      this.place(it, this.player.pos.clone().addScaledVector(f, 0.6).setY(0));
    }
  }
  dispose(it, { quiet = false } = {}) {
    this.inv.remove(it);
    it.obj.removeFromParent();
    this.items.splice(this.items.indexOf(it), 1);
    this.refreshPickRoots();
    this.updateHotbar();
    if (!quiet) HUD.toast(`${it.name} thrown away`);
  }
  openInventory() { this.invPick = null; HUD.inventoryGrids(this.inv, null, (it) => this.iconFor(it)); $("wearing").textContent = this.player.goggles ? "Wearing: safety goggles" : "Not wearing goggles"; this.openUI("inventory"); }

  // A small picture of an item for its hotbar slot, rendered once.
  iconFor(it) {
    const key = `${it.kind}:${it.name}`;
    if (this.iconCache.has(key)) return this.iconCache.get(key);
    const N = 96, r = this.renderer;
    if (!this.iconRT) {
      this.iconRT = new THREE.WebGLRenderTarget(N, N, { colorSpace: THREE.SRGBColorSpace });
      this.iconScene = new THREE.Scene();
      this.iconScene.add(new THREE.HemisphereLight(0xffffff, 0x777777, 2.2));
      const d = new THREE.DirectionalLight(0xffffff, 2.5); d.position.set(1, 2, 2); this.iconScene.add(d);
      this.iconCam = new THREE.PerspectiveCamera(30, 1, 0.001, 10);
    }
    const o = it.obj.clone();
    o.position.set(0, 0, 0); o.rotation.set(0, 0.6, 0); o.scale.setScalar(1);
    this.iconScene.add(o);
    this.iconScene.environment = this.envMap;
    const b = new THREE.Box3().setFromObject(o), c = b.getCenter(new THREE.Vector3()), sz = b.getSize(new THREE.Vector3()).length();
    this.iconCam.position.copy(c).add(new THREE.Vector3(0.45, 0.55, 0.9).multiplyScalar(sz * 1.25));
    this.iconCam.lookAt(c);
    r.setRenderTarget(this.iconRT);
    r.setClearColor(0x000000, 0);
    r.clear();
    r.render(this.iconScene, this.iconCam);
    const px = new Uint8Array(N * N * 4);
    r.readRenderTargetPixels(this.iconRT, 0, 0, N, N, px);
    r.setRenderTarget(null);
    this.iconScene.remove(o);
    const cv = document.createElement("canvas"); cv.width = cv.height = N;
    const g = cv.getContext("2d"), img = g.createImageData(N, N);
    for (let y = 0; y < N; y++) img.data.set(px.subarray((N - 1 - y) * N * 4, (N - y) * N * 4), y * N * 4);
    g.putImageData(img, 0, 0);
    const url = cv.toDataURL();
    this.iconCache.set(key, url);
    return url;
  }

  /* ---------- aiming and clicking ---------- */
  refreshPickRoots() {
    this.pickRoots = [...this.items.filter((it) => it.obj.parent && it.obj.parent !== this.hand && !this.inv.slots.includes(it)).map((it) => it.obj), ...this.room.targets, ...this.room.surfaces];
  }
  computeAim() {
    const ray = (this.ray = this.ray || new THREE.Raycaster());
    ray.setFromCamera({ x: 0, y: 0 }, this.camera);
    ray.far = REACH;
    const hits = ray.intersectObjects(this.pickRoots, true);
    for (const h of hits) {
      const o = h.object;
      if (!shown(o)) continue;   // rays hit hidden meshes too; skip them
      const up = h.face && h.face.normal.clone().transformDirection(o.matrixWorld).y > 0.7;
      const itemRoot = o.userData.itemRoot;
      // Items count as surfaces too, so you can set one thing on another.
      if (itemRoot && itemRoot.userData.item) return { item: itemRoot.userData.item, point: h.point, dist: h.distance, surface: up };
      let t = o;
      while (t && !t.userData.target) t = t.parent;
      if (t) return { furniture: t, point: h.point, dist: h.distance };
      if (up) return { surface: true, point: h.point, dist: h.distance };
      return null;
    }
    return null;
  }
  press(button) {
    if (!button) return;
    const a = this.actions;
    if (!a) return;
    const act = button === "left" ? a.left : a.right;
    if (!act) return;
    this.lastInput = performance.now();
    if (act.hold) { this.holding = { act, t: 0 }; HUD.hold(0); return; }
    if (act.whileHeld) { this.holdingLens = act; act.start(); return; }
    act.run();
    this.afterAction();
  }
  release() {
    if (this.holding) { const h = this.holding; this.holding = null; HUD.hold(null); h.act.run(h.t >= h.act.hold); this.afterAction(); }
    if (this.holdingLens) { this.holdingLens.end(); this.holdingLens = null; }
  }
  afterAction() { this.updateHotbar(); this.refreshPickRoots(); this.dirty = true; }
  emit(ev) { this.tasks.event(ev); }
  toast(text, kind) { HUD.toast(text, kind); }

  /* ---------- lab equipment ---------- */
  wearGoggles(it) {
    this.player.goggles = true;
    $("goggleFrame").hidden = false;
    this.dispose(it, { quiet: true });
    this.emit("goggles");
    HUD.toast("Goggles on", "good");
  }
  washHands() { this.player.cleanHands = true; this.emit("hands"); HUD.toast("Hands clean", "good"); }
  runTap() {
    const w = this.room.water;
    w.visible = true;
    clearTimeout(this.tapT);
    this.tapT = setTimeout(() => { w.visible = false; this.dirty = true; }, 1600);
  }
  toggleLights() {
    this.lightsOn = !this.lightsOn;
    this.room.lamps.forEach((l, i) => (l.visible = this.lightsOn && i < this.preset.lamps));
    this.room.tubeMat.color.set(this.lightsOn ? 0xf4f7ff : 0x3a3d40);
    this.room.lightSwitch.userData.rocker.rotation.x = this.lightsOn ? 0 : 0.3;
    this.room.hemi.intensity = (this.preset.env ? 0.35 : 1.25) * (this.lightsOn ? 1 : 0.4);
  }
  setBlinds(b, down, instant) {
    const stack = b.userData.stack;
    b.userData.down = down;
    const sy = down ? 1 : 0.12;
    if (instant) stack.scale.y = sy;
    else this.tween(0.6, stack.scale.y, sy, (v) => { stack.scale.y = v; });
    const open = this.room.blinds.filter((x) => !x.userData.down).length;
    this.room.sun.intensity = 0.6 + (open / this.room.blinds.length) * 2.6;
  }
  toggleBlinds(b) { this.setBlinds(b, !b.userData.down); }
  toggleDrawer(d) {
    d.userData.open = !d.userData.open;
    if (!d.userData.spare) d.userData.spare = { ...SPARES[d.userData.target.index % SPARES.length], left: 3 };
    const from = d.position.z, to = d.userData.closedZ + (d.userData.open ? 0.3 : 0);
    this.tween(0.35, from, to, (v) => (d.position.z = v));
    if (d.userData.open && d.userData.spare) HUD.toast(`Inside: ${d.userData.spare.label} (left-click to take)`);
  }
  takeSpare(d) {
    const s = d.userData.spare;
    if (!s || s.left <= 0) return HUD.toast("This drawer is empty.");
    s.left--;
    this.give(s.make());
    if (s.left <= 0) d.userData.spare = { ...s, label: "nothing left", left: 0 };
  }
  tween(dur, a, b, fn) { this.anims.push({ t: 0, dur, a, b, fn }); }

  putOnStage(it) {
    this.inv.remove(it);
    const scope = this.room.scope;
    it.obj.removeFromParent();
    it.obj.scale.setScalar(1);
    it.obj.position.copy(scope.userData.stage);
    it.obj.rotation.set(0, Math.PI / 2, 0);
    scope.add(it.obj);
    this.stageItem = it;
    this.emit(it.kind === "dish" ? "dish-stage" : "stage");
    HUD.toast(`${it.name} clipped onto the stage. Right-click the microscope to look.`);
    this.afterAction();
  }
  takeFromStage() { if (this.stageItem) this.take(this.stageItem); }

  async enterScope() {
    const it = this.stageItem;
    if (!it) return;
    this.screen = "scope";
    document.exitPointerLock?.();
    this.unlockingForUI = true;
    $("hud").hidden = true;
    $("scope").hidden = false;
    if (it.kind === "slide") {
      const st = it.state, stain = it.container.stain() || "water";
      const names = { onion: "Onion skin", cheek: "Cheek cells", elodea: "Pondweed leaf" };
      const stainName = stain === "iodine" ? "Iodine" : stain === "mb" ? "Methylene blue" : it.container.empty ? "Dry" : "Water";
      this.scope.open({ kind: "slide", title: `${st.sample ? names[st.sample] : "Empty slide"} · ${stainName.toLowerCase()}`, stainName,
        prep: { specimen: st.sample, stain, cover: st.cover || "none" } });
      return;
    }
    // A living culture: show what the engine says the dish looks like now.
    this.scope.open({ kind: "wait", title: it.name, text: this.engine.version ? "Simulating the cells…" : "Starting the Cell Sim engine (the first time downloads Python, about 30 s)…" });
    try {
      await this.simulate(it);
      if (this.screen === "scope" && this.stageItem === it) this.scope.open(this.cultureSpec(it));
    } catch (e) {
      this.scope.open({ kind: "wait", title: it.name, text: `The engine could not run: ${e.message}` });
    }
  }
  cultureSpec(d) {
    const { color } = d.container.color();
    const tint = `#${color.clone().lerp(new THREE.Color(0xe8e4dc), 0.82).getHexString()}`;
    return { kind: "culture", title: `${d.name} · A549 · ${d.state.events.length ? "treated" : "untreated"}`, header: d.state.header, snap: d.state.snap, time: d.state.t, tint };
  }
  exitScope() {
    this.screen = "play";
    $("scope").hidden = true;
    $("hud").hidden = false;
    this.lock();
    this.dirty = true;
  }
  onFocus(obj, spec) {
    if (spec.kind === "slide") this.emit(`focus:${obj}`);
    else {
      const d = this.stageItem;
      if (obj >= 10) {
        this.emit("observe:culture");
        if (d && d.state.t > 0 && d.state.events.length) this.emit("observe:after");
        if (d && d.name === "Dish B" && d.state.t > 0) this.emit("compare");
      }
    }
  }

  // Run the engine for a dish up to its current time.
  async simulate(d) {
    const st = d.state;
    if (st.snap && !st.needs) return;
    const drugs = [...new Set(st.events.flatMap((e) => Object.keys(e.conc)))];
    const out = await this.engine.dish({ label: d.name, line: st.line, drugs, events: st.events, t: st.t, n_seed: 80, grid: 60, seed: st.seed, every: 2 });
    st.header = out.header; st.snap = out.snap; st.series = out.series; st.needs = false;
  }
  dishDosed(d) {
    const conc = d.container.drugs(), st = d.state;
    const ev = st.events.find((e) => e.t === st.t);
    if (ev) ev.conc = conc; else st.events.push({ t: st.t, conc });
    st.needs = true;
    this.emit("drug-drop");
    HUD.toast(`${d.name} now has ${Object.entries(conc).map(([k, v]) => `${k} ${v.toFixed(2)} µM`).join(", ")}`);
  }
  putInIncubator(d) {
    this.inv.remove(d);
    d.obj.removeFromParent();
    this.incubator.dishes.push(d);
    this.afterAction();
    HUD.toast(`${d.name} is in the incubator. Right-click it to let time pass.`);
  }
  openIncubator() { this.renderIncubator(); this.openUI("incubatorPanel"); }
  renderIncubator() {
    const ds = this.incubator.dishes;
    $("incubatorBody").innerHTML = (ds.length ? ds.map((d) => `<div class="dish-row"><b>${d.name}</b><span class="note">${d.state.t} h in culture · ${d.state.events.length ? "treated" : "no drug"}</span>
        <div class="btns"><button class="mbtn small" data-inc="take" data-id="${d.id}">Take it out</button></div></div>`).join("") : `<p class="note">The incubator is empty.</p>`) +
      (ds.length ? `<p class="glabel">Let time pass for every dish inside</p><div class="btns" style="display:flex;gap:6px">${[6, 12, 24].map((h) => `<button class="mbtn small primary" data-inc="wait" data-h="${h}">Wait ${h} h</button>`).join("")}</div>` : "") +
      `<p class="note">${this.engine.version ? `Cell Sim ${this.engine.version} simulates each dish.` : "The Cell Sim engine starts the first time you incubate or look at a dish."}</p>`;
  }
  async incubate(h) {
    const ds = [...this.incubator.dishes];
    for (const d of ds) { d.state.t += h; d.state.needs = true; }
    this.renderIncubator();
    if (ds.some((d) => d.state.events.length)) this.emit("incubate");
    try {
      for (const d of ds) {
        await this.simulate(d);
        const last = d.state.series[d.state.series.length - 1];
        HUD.toast(`${d.name} at ${d.state.t} h: ${last.n_live} live cells, ${last.n_apoptotic} dying`, "good");
      }
    } catch (e) { HUD.toast(`The engine could not run: ${e.message}`, "warn"); }
    this.renderIncubator();
  }
  openJournal() { this.journalTab = this.journalTab || "tasks"; this.renderJournal(); this.openUI("journal"); }
  renderJournal() {
    for (const b of document.querySelectorAll("[data-tab]")) b.setAttribute("aria-selected", String(b.dataset.tab === this.journalTab));
    $("journalBody").innerHTML = this.journal.render(this.journalTab, this.tasks, this.dishes);
  }
  openReader() {
    const e = this.tasks.exp;
    $("readerBody").innerHTML = e
      ? `<p class="kicker">Method</p><h2>${e.title}</h2><ol>${e.steps.map(([, t]) => `<li>${t}</li>`).join("")}</ol>`
      : `<p class="kicker">Method</p><h2>Cells under the microscope</h2><p>Pick an experiment from the title screen for a step-by-step method. In free roam, everything in the lab works.</p>`;
    this.openUI("reader");
  }
  leave() {
    this.journal.save();
    this.player.keys.clear();
    this.showTitle();
    $("mainMenu").querySelector('[data-menu="experiments"]').textContent = "Start a new experiment";
  }

  /* ---------- the loop ---------- */
  frame(now) {
    requestAnimationFrame((t) => this.frame(t));
    const p = this.preset;
    if (p.fps && now - (this.lastDraw || 0) < 1000 / p.fps - 1.5) return;
    const dt = Math.min(0.05, (now - this.last) / 1000);
    this.last = now;
    let animating = false;
    for (let i = this.anims.length - 1; i >= 0; i--) {
      const a = this.anims[i];
      a.t = Math.min(1, a.t + dt / a.dur);
      const k = a.t < 0.5 ? 2 * a.t * a.t : 1 - Math.pow(-2 * a.t + 2, 2) / 2;
      a.fn(a.a + (a.b - a.a) * k, k);
      if (a.t >= 1) this.anims.splice(i, 1);
      animating = true;
    }
    if (this.screen === "scope") { this.scope.update(dt); this.lastDraw = now; return; }
    if (this.screen === "title") {
      const t = now / 1000 * 0.06;
      this.camera.position.set(-0.6 + Math.cos(t) * 3.4, 1.85, 0.2 + Math.sin(t) * 2.6);
      this.camera.lookAt(-1.2, 1.0, -0.6);
      this.render(now);
      return;
    }
    if (this.screen !== "play") return;
    if (!this.ui) this.player.update(dt, this.room.colliders);
    this.updateClock();
    // Aim and prompts.
    if (!this.ui) {
      // Frames are skipped when nothing moves, so update positions here too:
      // an item just put down must be hit where it now is.
      this.scene.updateMatrixWorld();
      this.aim = this.computeAim();
      this.actions = getActions(this, this.held, this.aim);
      this.tick = (this.tick || 0) + 1;
      const live = this.locked || this.forceLock;
      const sig = live ? `${this.actions.left?.label}|${this.actions.right?.label}|${this.actions.target}` : "unlocked";
      if (sig !== this.promptSig) { this.promptSig = sig; HUD.prompt(live ? this.actions : { target: "Click to look around" }); }
    }
    if (this.holding) {
      this.holding.t += dt;
      HUD.hold(Math.min(1, this.holding.t / this.holding.act.hold));
      animating = true;
    }
    if (this.held) HUD.held(this.held);
    // Only draw when something changed: saves battery when you stand still.
    const recent = now - (this.lastInput || 0) < 300;
    const tick = Math.floor(now / 1000) !== Math.floor((this.lastDraw || 0) / 1000);
    if (this.player.moved || recent || animating || this.dirty || tick || this.room.water.visible) this.render(now);
    this.debug(now, dt);
  }
  render(now) {
    this.renderer.render(this.scene, this.camera);
    this.lastDraw = now;
    this.dirty = false;
    this.frames = (this.frames || 0) + 1;
  }
  updateClock() {
    const h = this.room.clockHands, d = new Date();
    if (!h.s) return;
    const sec = d.getSeconds(), min = d.getMinutes() + sec / 60, hr = (d.getHours() % 12) + min / 60;
    h.s.rotation.z = -sec / 60 * Math.PI * 2;
    if (h.m) h.m.rotation.z = -min / 60 * Math.PI * 2;
    if (h.h) h.h.rotation.z = -hr / 12 * Math.PI * 2;
  }
  debug(now) {
    if ($("debug").hidden) return;
    if (now - (this.dbgT || 0) < 500) return;
    const fps = ((this.frames || 0) - (this.dbgF || 0)) / ((now - (this.dbgT || now - 500)) / 1000);
    this.dbgT = now; this.dbgF = this.frames;
    const i = this.renderer.info, mem = performance.memory ? `${(performance.memory.usedJSHeapSize / 1048576).toFixed(0)} MB JS heap` : "";
    $("debug").textContent = `${fps.toFixed(0)} fps drawn · ${this.preset.label} graphics\n${i.render.calls} draw calls · ${(i.render.triangles / 1000).toFixed(0)}k triangles\n${i.memory.geometries} geometries · ${i.memory.textures} textures ${mem}\nGPU: ${this.detected.gpu || "unknown"}\nEngine: ${this.engine.version || "not loaded"}`;
  }
}

function shown(o) {
  for (let x = o; x; x = x.parent) if (!x.visible) return false;
  return true;
}

const SPARES = [
  { label: "a new box of slides", make: () => { const it = new Item("slidebox", "Box of slides", P.makeBox({ w: 0.085, h: 0.03, d: 0.035, color: 0x2a5fb0, labelLines: ["Microscope slides", "76 × 26 mm"] })); it.state.left = 10; return it; } },
  { label: "a new pack of cover slips", make: () => { const it = new Item("coverbox", "Cover slips", P.makeBox({ w: 0.032, h: 0.014, d: 0.032, color: 0xdfe6ea, labelLines: ["Cover slips"] })); it.state.left = 20; return it; } },
  { label: "a fresh 30 mL bottle of iodine", make: () => { const m = P.makeDropperBottle("iodine", 30); return new Item("dropper", "Iodine dropper bottle", m.obj, { container: m.container, liquid: m.liquid, solution: "iodine" }); } },
  { label: "a fresh 30 mL bottle of methylene blue", make: () => { const m = P.makeDropperBottle("mb", 30); return new Item("dropper", "Methylene blue dropper bottle", m.obj, { container: m.container, liquid: m.liquid, solution: "mb" }); } },
];

const CREDITS = `
<p><b>Cell engine:</b> <a href="https://github.com/Rockman6/CellSim-3D-Multicellular-Biology-Drug-Testing-Simulator" target="_blank" rel="noopener">Cell Sim</a> (MIT), running in your browser with <a href="https://pyodide.org" target="_blank" rel="noopener">Pyodide</a>.</p>
<p><b>3D engine:</b> <a href="https://threejs.org" target="_blank" rel="noopener">three.js</a> r160 (MIT).</p>
<p><b>Models, textures and lighting from Poly Haven (CC0):</b> yellow onion, metal stool, steel shelves, fluorescent lights, drawer cabinet, clipboard, stationery, wall clock, magnifying glass, binder notebook, metal office desk, floor tiles, painted plaster, hospital room HDRI.</p>
<p><b>Models from Sketchfab (CC Attribution 4.0):</b>
"Chemistry Lab Table" by Jawahar Yokesh · "lab bench" by Konstantin_Keller (the incubator) · "Compound Microscope" by kg007 · "Door Classroom" by Mehdi Shahsavan.</p>
<p><b>Built in code:</b> the room, slides, cover slips, bottles, beakers, culture dishes, the pondweed tank, the lab tools, the safety goggles, the hand sanitizer and the window blinds.</p>`;

new Game().boot().catch((e) => {
  console.error(e);
  $("loadText").textContent = `Something went wrong while loading: ${e.message}`;
});
