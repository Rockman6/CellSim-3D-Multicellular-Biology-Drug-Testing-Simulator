// CellSim: starts the renderer, loads the lab, and runs the game.
import * as THREE from "three";
import { loadAssets, enableShadows } from "./assets.js";
import { buildRoom } from "./room.js";
import { Player } from "./player.js";
import { createItems, getActions, Item, fmtConc, sentence } from "./items.js";
import { HUD } from "./hud.js";
import { Tasks, EXPERIMENTS } from "./tasks.js";
import { Journal } from "./journal.js";
import { createScope } from "./scope.js";
import { Engine } from "./engine.js";
import { Outliner } from "./outline.js";
import { search, buildIndex, norm, Tracker } from "./finder.js";
import { sfx, unlockAudio } from "./audio.js";
import { settings, saveSettings, detectPreset, activePreset } from "./settings.js";
import * as P from "./props.js";

const $ = (id) => document.getElementById(id);
const REACH = 2.4;
const DISH_AREA_UM2 = 2.15e9;   // growth area of a 60 mm dish, 21.5 cm²
const GLASSY = new Set(["slide", "coverslip", "beaker", "dropper", "stock", "drug", "dish", "tube", "hemocytometer", "watchglass", "trypsin", "plate"]);
const TIPS = [
  "Left-click takes things into your pocket. Right-click uses them.",
  "Can't find something? Press F and type its name. Typos are fine.",
  "Lower a cover slip slowly by holding right-click: a flat drop traps air bubbles.",
  "Start every slide with the 4× lens, then go up one lens at a time.",
  "The living cells in the incubator are simulated by the CellSim engine, cell by cell.",
  "A fresh pipette tip for every solution, or you carry one into the next.",
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
    // The lab barely moves, so the sun's shadows are redrawn only when something does.
    r.shadowMap.autoUpdate = false;
    this.shadowsDirty = true;
    r.info.autoReset = false;   // F3 counts the whole frame, outline pass included
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x222629);
    this.camera = new THREE.PerspectiveCamera(settings.fov, 1, 0.03, 60);
    this.scene.add(this.camera);
    this.hand = new THREE.Group();
    this.camera.add(this.hand);
    this.outliner = new Outliner(r, { scale: this.preset.key === "low" ? 0.4 : 0.6 });
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
    this.player.onStep = (loud) => sfx.step(loud);
    this.inv = new Inventory();
    this.selected = 0;
    this.tasks = new Tasks((key) => { HUD.tasks(this.tasks.view()); if (key) { HUD.toast("Step done ✓", "good"); sfx.ok(); } });
    this.journal = new Journal();
    this.engine = new Engine((state, text) => HUD.engine(state, text));
    const made = createItems(this, assets, this.room);
    this.items = made.items;
    this.plate = made.plate;
    this.incubator = { dishes: [...made.dishes, made.plate] };
    this.dishes = made.dishes;
    this.stageItem = null;
    this.lightsOn = true;
    this.anims = [];
    this.iconCache = new Map();
    this.tracker = new Tracker(this);
    this.ghost = new THREE.Mesh(new THREE.RingGeometry(0.86, 1, 40), new THREE.MeshBasicMaterial({ color: 0x58c46e, transparent: true, opacity: 0.75, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -4 }));
    this.ghost.rotation.x = -Math.PI / 2;
    this.ghost.visible = false;
    this.ghost.raycast = () => {};
    this.scene.add(this.ghost);
    this.scope = createScope({
      resolution: this.preset.key === "low" ? 560 : this.preset.key === "ultra" ? 1024 : 768,
      onFocus: (obj, spec) => this.onFocus(obj, spec),
      onExit: () => this.exitScope(),
      toast: (t, k) => this.toast(t, k),
      observe: (s, p) => this.journal.observe(s, p),
      onEvent: (name, data) => this.onScopeEvent(name, data),
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
    const w = innerWidth, h = innerHeight, pr = Math.min(devicePixelRatio || 1, this.preset.pixelRatio);
    this.renderer.setPixelRatio(pr);
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.outliner.resize(w, h, pr);
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
    this.shadowsDirty = true;
    this.resize();
    $("qualityNote").textContent = `Using ${p.label} graphics${settings.quality === "auto" ? ` (picked for ${this.detected.gpu || "this computer"})` : ""}. ${p.fps ? `Capped at ${p.fps} fps.` : "No frame cap."}`;
  }

  /* ---------- screens ---------- */
  showTitle() {
    this.screen = "title";
    $("hud").hidden = true;
    $("title").hidden = false;
    if (matchMedia("(pointer: coarse)").matches && !matchMedia("(any-pointer: fine)").matches && !$("touchNote"))
      $("mainMenu").insertAdjacentHTML("afterbegin", `<p id="touchNote" class="note">CellSim needs a keyboard and mouse. On this device you can look around, but not play.</p>`);
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
    if (["culture", "dose", "count"].includes(id)) this.engine.start().catch(() => {});
    if (lock) this.lock();
    HUD.toast("WASD to walk, mouse to look. Left-click takes, right-click uses. F finds anything.");
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
      else if (a === "finder") { this.closeAllUI(); this.openFinder(); }
      else if (a === "title") { this.closeAllUI(); this.leave(); }
      else this.openUI(a);
    });
    for (const b of document.querySelectorAll("[data-back]")) b.addEventListener("click", () => { this.closeUI(); if (!this.ui && this.screen === "play") this.lock(); });
    // Settings
    const q = $("setQuality"), fps = $("setFps"), sens = $("setSens"), fov = $("setFov"), inv = $("setInvert"), bob = $("setBob"), snd = $("setSound");
    q.value = settings.quality; fps.value = settings.fps; sens.value = settings.sensitivity; fov.value = settings.fov; inv.checked = settings.invertY; bob.checked = settings.bob; snd.checked = settings.sound;
    const onSet = () => {
      Object.assign(settings, { quality: q.value, fps: fps.value, sensitivity: +sens.value, fov: +fov.value, invertY: inv.checked, bob: bob.checked, sound: snd.checked });
      saveSettings();
      if (settings.sound) unlockAudio();
      this.preset = activePreset(this.detected.pick);
      this.applyPreset();
    };
    for (const el of [q, fps, sens, fov, inv, bob, snd]) el.addEventListener("change", onSet);
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
        if (d && this.give(d)) this.emit(d.kind === "plate" ? "plate-taken" : "dish-taken");   // take() also removes it from the incubator
        this.renderIncubator();
      }
    });
    // Find anything
    const inp = $("findInput");
    inp.addEventListener("input", () => { this.findSel = 0; this.renderFind(); });
    inp.addEventListener("keydown", (e) => {
      const n = (this.findList || []).length;
      if (e.key === "ArrowDown") { e.preventDefault(); this.findSel = Math.min(n - 1, (this.findSel || 0) + 1); this.renderFind(false); }
      else if (e.key === "ArrowUp") { e.preventDefault(); this.findSel = Math.max(0, (this.findSel || 0) - 1); this.renderFind(false); }
      else if (e.key === "Enter") { e.preventDefault(); if (n) this.chooseFind(this.findSel || 0); }
    });
    $("findResults").addEventListener("click", (e) => { const li = e.target.closest("[data-i]"); if (li) this.chooseFind(+li.dataset.i); });
    $("finder").addEventListener("click", (e) => { if (e.target.id === "finder") { this.closeUI(); this.lock(); } });
    // Micropipette volume
    const setVol = (v) => {
      const p = this.pipItem;
      if (!p) return;
      p.state.vol = THREE.MathUtils.clamp(Math.round(v), 2, 200);
      p.refresh();
      $("pipVol").textContent = p.state.vol; $("pipRange").value = p.state.vol;
      HUD.held(this.held);
      sfx.click();
      if (p.state.vol === 10) this.emit("vol:10");
    };
    $("pipRange").addEventListener("input", (e) => setVol(+e.target.value));
    $("pipettePanel").addEventListener("click", (e) => {
      const b = e.target.closest("[data-pip],[data-pipset]");
      if (!b || !this.pipItem) return;
      setVol(b.dataset.pipset ? +b.dataset.pipset : this.pipItem.state.vol + +b.dataset.pip);
    });
  }

  renderExperiments() {
    $("experimentMenu").innerHTML = Object.entries(EXPERIMENTS).map(([id, e]) =>
      `<button class="exp-card" data-exp="${id}"><span class="tag">${e.tag}</span><b>${e.title}</b><span>${e.blurb}</span></button>`).join("") +
      `<button class="mbtn small" data-exp="back">Back</button>`;
  }

  bindInput() {
    const c = this.canvas;
    c.addEventListener("click", () => { unlockAudio(); if (this.screen === "play" && !document.pointerLockElement) this.lock(); });
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
      unlockAudio();
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
      else if (e.code === "KeyF") { e.preventDefault(); this.tracker.target ? this.tracker.stop() : this.openFinder(); }
      else if (e.code === "KeyR" && this.held?.kind === "pipette") this.openPipette(this.held);
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
    // Show the selected item in your hand. It casts no shadow there: the
    // shadows are only redrawn when something in the lab moves.
    while (this.hand.children.length) this.hand.remove(this.hand.children[0]);
    const it = this.held;
    if (it) {
      it.obj.userData.pose?.("hand");
      if (!it.size) {
        it.obj.removeFromParent();
        it.obj.position.set(0, 0, 0); it.obj.rotation.set(0, 0, 0); it.obj.scale.setScalar(1);
        it.obj.updateMatrixWorld(true);
        it.size = new THREE.Box3().setFromObject(it.obj).getSize(new THREE.Vector3());
      }
      const size = it.size;
      const s = THREE.MathUtils.clamp(0.15 / Math.max(size.x, size.y, size.z, 0.001), 0.5, 2.4);
      it.obj.position.set(0.2, -0.2, -0.42);
      it.obj.rotation.set(it.kind === "pipette" ? 0.05 : 0.15, -0.5, it.kind === "pipette" ? 0.25 : 0);
      it.obj.scale.setScalar(it.kind === "pipette" ? 1.1 : s);
      it.obj.traverse((m) => { if (m.isMesh) m.castShadow = false; });
      if (it.obj.userData.blob) it.obj.userData.blob.visible = false;
      this.hand.add(it.obj);
    }
    this.dirty = true;
  }
  take(it) {
    if (it === this.stageItem) { this.stageItem = null; }
    const empty = !this.held;
    const i = this.inv.add(it, empty ? this.selected : null);
    if (i < 0) { this.toast("Your pockets are full. Put something down first (Q).", "warn"); return false; }
    it.obj.removeFromParent();
    it.obj.scale.setScalar(1);
    it.size = null;
    if (it.kind === "tube") it.state.rackHole = null;
    this.refreshPickRoots();
    this.updateHotbar();
    HUD.itemName(it.name);
    if (this.incubator.dishes.includes(it)) this.incubator.dishes.splice(this.incubator.dishes.indexOf(it), 1);
    this.shadowsDirty = true;
    sfx.take();
    return true;
  }
  give(it) { if (!this.items.includes(it)) this.items.push(it); return this.take(it); }
  place(it, point) {
    this.inv.remove(it);
    it.obj.removeFromParent();
    it.obj.userData.pose?.("rest");
    it.obj.scale.setScalar(1);
    it.obj.position.copy(point);
    it.obj.rotation.set(0, this.player.yaw + Math.PI / 2, 0);
    enableShadows(it.obj);
    if (it.obj.userData.blob) it.obj.userData.blob.visible = true;
    this.scene.add(it.obj);
    this.refreshPickRoots();
    this.updateHotbar();
    this.shadowsDirty = true;
    GLASSY.has(it.kind) ? sfx.glass() : sfx.place();
  }
  dropHeld() {
    const it = this.held;
    if (!it) return;
    const hit = this.aim;
    if (hit && hit.surface && (!hit.item || hit.item.kind === "tile") && !this.placeBlocked(it, hit.point, hit.item)) this.place(it, hit.point);
    else {
      const f = new THREE.Vector3(-Math.sin(this.player.yaw), 0, -Math.cos(this.player.yaw));
      this.place(it, this.player.pos.clone().addScaledVector(f, 0.6).setY(0));
    }
  }
  // Something already standing where you're about to put an item down.
  placeBlocked(held, p, onItem) {
    const r = held.foot || 0.03;
    for (const it of this.items) {
      if (it === held || it === onItem || it.obj.parent !== this.scene) continue;
      if (it.kind === "tube" && it.state.rackHole != null) continue;
      const q = it.obj.position;
      if (q.y > p.y + 0.02 || q.y < p.y - 0.25) continue;   // on another surface
      if (Math.hypot(q.x - p.x, q.z - p.z) < (r + (it.foot || 0.03)) * 0.8) return it;
    }
    return null;
  }
  dispose(it, { quiet = false } = {}) {
    this.inv.remove(it);
    it.obj.removeFromParent();
    this.items.splice(this.items.indexOf(it), 1);
    this.refreshPickRoots();
    this.updateHotbar();
    this.shadowsDirty = true;
    if (!quiet) this.toast(`${it.name} thrown away`);
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
    o.traverse((m) => { if (m.name === "shadow") m.visible = false; });
    this.iconScene.add(o);
    this.iconScene.environment = this.envMap;
    const b = new THREE.Box3().setFromObject(o), c = b.getCenter(new THREE.Vector3()), sz = b.getSize(new THREE.Vector3()).length();
    this.iconCam.position.copy(c).add(new THREE.Vector3(0.45, 0.55, 0.9).multiplyScalar(sz * 1.25));
    this.iconCam.lookAt(c);
    const auto = r.shadowMap.autoUpdate;
    r.setRenderTarget(this.iconRT);
    r.setClearColor(0x000000, 0);
    r.clear();
    r.render(this.iconScene, this.iconCam);
    const px = new Uint8Array(N * N * 4);
    r.readRenderTargetPixels(this.iconRT, 0, 0, N, N, px);
    r.setRenderTarget(null);
    r.shadowMap.autoUpdate = auto;
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
      // Items count as surfaces too, so you can set a slide on the tile.
      if (itemRoot && itemRoot.userData.item) return { item: itemRoot.userData.item, point: h.point, dist: h.distance, surface: up };
      let t = o;
      while (t && !t.userData.target) t = t.parent;
      if (t) return { furniture: t, point: h.point, dist: h.distance };
      if (up) return { surface: true, point: h.point, dist: h.distance };
      return null;
    }
    return null;
  }
  // The white outline, the crosshair, and the ring where a held item would go.
  updateAimCues() {
    const a = this.actions, aim = this.aim;
    const useful = !!(a && (a.left || a.right));
    const obj = !useful || !aim ? null : aim.item ? (aim.item === this.stageItem ? this.room.scope : aim.item.obj) : aim.furniture || null;
    if (this.outliner.set("hover", obj)) this.dirty = true;
    const ch = (this.crosshair = this.crosshair || document.querySelector(".crosshair"));
    ch.classList.toggle("on", useful && !a.right?.blocked);
    ch.classList.toggle("blocked", !!a?.right?.blocked);
    const ring = a && (a.right?.place || a.right?.blocked) && aim?.point;
    if (ring) {
      this.ghost.position.copy(aim.point).y += 0.0015;
      this.ghost.scale.setScalar(this.held?.foot || 0.03);
      this.ghost.material.color.set(a.right.blocked ? 0xff7466 : 0x58c46e);
    }
    if (!!ring !== this.ghost.visible || ring) this.dirty = true;
    this.ghost.visible = !!ring;
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
  afterAction() { this.updateHotbar(); this.refreshPickRoots(); this.dirty = true; this.shadowsDirty = true; }
  emit(ev) { this.tasks.event(ev); }
  toast(text, kind) { HUD.toast(text, kind); if (kind === "warn") sfx.warn(); }
  sfx(name) { sfx[name]?.(); }

  /* ---------- finding things ---------- */
  openFinder() {
    if (this.screen !== "play") return;
    this.openUI("finder");
    const inp = $("findInput");
    inp.value = "";
    this.findSel = 0;
    this.renderFind();
    setTimeout(() => inp.focus(), 0);
  }
  renderFind(rebuild = true) {
    const q = $("findInput").value;
    let head = "";
    if (rebuild) {
      if (q.trim()) this.findList = search(this, q);
      else {
        // Nothing typed yet: what the current step needs.
        const v = this.tasks.view();
        const step = v && !v.finished ? v.steps[v.now].text : "";
        const words = new Set(norm(step).split(" ").map((w) => w.replace(/s$/, "")));
        const pos = this.player.pos;
        this.findList = step ? buildIndex(this)
          .filter((e) => e.terms.some((t) => t.length > 3 && t.split(" ").every((w) => words.has(w.replace(/s$/, "")))))
          .map((e) => ({ ...e, dist: e.object ? e.object.getWorldPosition(new THREE.Vector3()).distanceTo(pos) : 0 }))
          .sort((a, b) => a.dist - b.dist).slice(0, 5) : [];
        this.findHead = this.findList.length ? "For this step" : "";
      }
    }
    head = q.trim() ? "" : this.findHead;
    const list = this.findList || [];
    const sel = Math.min(this.findSel || 0, Math.max(0, list.length - 1));
    $("findResults").innerHTML = list.length
      ? (head ? `<li class="sect">${head}</li>` : "") + list.map((e, i) => `<li role="option" data-i="${i}" aria-selected="${i === sel}"><b>${esc(e.name)}</b><span class="where">${esc(capital(e.where || "in the lab"))}</span><span class="dist">${e.object ? `${e.dist.toFixed(1)} m` : "with you"}</span></li>`).join("")
      : `<li class="empty">${q.trim() ? "Nothing like that in the lab. Try another word." : "Type what you're looking for."}</li>`;
  }
  chooseFind(i) {
    const e = (this.findList || [])[i];
    if (!e) return;
    this.closeUI();
    if (!e.object) {
      const slot = this.inv.slots.indexOf(e.item);
      if (slot >= 0 && slot < 9) { this.select(slot); this.toast(`${e.name} is in your hotbar, slot ${slot + 1}: it's in your hand now.`); }
      else this.toast(`${e.name} is in your inventory. Press E to move it to your hotbar.`);
    } else {
      this.tracker.start(e);
      this.toast(`Follow the blue marker to the ${sentence(e.name)} (${e.where || "in the lab"}). Press F to stop.`);
    }
    this.lock();
  }

  /* ---------- the micropipette ---------- */
  openPipette(p) {
    this.pipItem = p;
    $("pipVol").textContent = Math.round(p.state.vol);
    $("pipRange").value = p.state.vol;
    this.openUI("pipettePanel");
  }

  /* ---------- lab equipment ---------- */
  wearGoggles(it) {
    this.player.goggles = true;
    $("goggleFrame").hidden = false;
    this.dispose(it, { quiet: true });
    this.emit("goggles");
    this.toast("Goggles on", "good");
  }
  washHands() { this.player.cleanHands = true; this.emit("hands"); this.toast("Hands clean", "good"); }
  runTap() {
    const w = this.room.water;
    w.visible = true;
    sfx.pour();
    clearTimeout(this.tapT);
    this.tapT = setTimeout(() => { w.visible = false; this.dirty = true; }, 1600);
  }
  toggleLights() {
    this.lightsOn = !this.lightsOn;
    this.room.lamps.forEach((l, i) => (l.visible = this.lightsOn && i < this.preset.lamps));
    this.room.tubeMat.color.set(this.lightsOn ? 0xf4f7ff : 0x3a3d40);
    this.room.lightSwitch.userData.rocker.rotation.x = this.lightsOn ? 0 : 0.3;
    this.room.hemi.intensity = (this.preset.env ? 0.35 : 1.25) * (this.lightsOn ? 1 : 0.4);
    sfx.click();
  }
  setBlinds(b, down, instant) {
    const stack = b.userData.stack;
    b.userData.down = down;
    const sy = down ? 1 : 0.12;
    if (instant) stack.scale.y = sy;
    else this.tween(0.6, stack.scale.y, sy, (v) => { stack.scale.y = v; });
    const open = this.room.blinds.filter((x) => !x.userData.down).length;
    this.room.sun.intensity = 0.6 + (open / this.room.blinds.length) * 2.6;
    this.shadowsDirty = true;
  }
  toggleBlinds(b) { this.setBlinds(b, !b.userData.down); }
  toggleDrawer(d) {
    d.userData.open = !d.userData.open;
    if (!d.userData.spare) d.userData.spare = { ...SPARES[d.userData.target.index % SPARES.length], left: 3 };
    const from = d.position.z, to = d.userData.closedZ + (d.userData.open ? 0.3 : 0);
    this.tween(0.35, from, to, (v) => (d.position.z = v));
    sfx.place();
    if (d.userData.open && d.userData.spare) this.toast(`Inside: ${d.userData.spare.label} (left-click to take)`);
  }
  takeSpare(d) {
    const s = d.userData.spare;
    if (!s || s.left <= 0) return this.toast("This drawer is empty.");
    s.left--;
    this.give(s.make());
    if (s.left <= 0) d.userData.spare = { ...s, label: "nothing left", left: 0 };
  }
  tween(dur, a, b, fn) { this.anims.push({ t: 0, dur, a, b, fn }); }
  inWaterBath(wg) {
    const bath = this.items.find((i) => i.kind === "waterbath");
    if (!bath || wg.obj.parent !== this.scene) return false;
    const a = wg.obj.position, b = bath.obj.position;
    return Math.hypot(a.x - b.x, a.z - b.z) < 0.12 && a.y > b.y + 0.05;
  }

  putOnStage(it) {
    this.inv.remove(it);
    const scope = this.room.scope;
    it.obj.removeFromParent();
    it.obj.scale.setScalar(1);
    it.obj.position.copy(scope.userData.stage);
    it.obj.rotation.set(0, Math.PI / 2, 0);
    if (it.obj.userData.blob) it.obj.userData.blob.visible = false;
    scope.add(it.obj);
    this.stageItem = it;
    this.emit({ dish: "dish-stage", plate: "plate-stage", hemocytometer: "hemo-stage" }[it.kind] || "stage");
    this.toast(`${it.name} clipped onto the stage. Right-click the microscope to look.`);
    sfx.glass();
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
      const names = { onion: "Onion skin", cheek: "Cheek cells", elodea: "Pondweed leaf", redonion: "Red onion skin", root: "Root tip squash" };
      const stainName = { iodine: "Iodine", mb: "Methylene blue", toluidine: "Toluidine blue" }[stain] || (it.container.empty ? "Dry" : P.shortName(it.container.dominant()));
      st.tally = st.tally || { counts: [0, 0, 0, 0, 0], total: 0, right: 0, named: new Set() };
      const last = st.osmo[st.osmo.length - 1];
      this.scope.open({ kind: "slide", title: `${st.sample ? names[st.sample] : "Empty slide"} · ${sentence(stainName)}`, stainName, tally: st.tally,
        medium: last ? `${sentence(P.shortName(st.mediumName || "water"))} (${last.osm.toFixed(2)} osmol/L)` : "nothing",
        prep: { specimen: st.sample, stain, cover: st.cover || "none", softened: st.softened, squashed: st.squashed, osmo: st.osmo } });
      return;
    }
    if (it.kind === "hemocytometer") return this.enterHemo(it);
    if (it.kind === "plate") return this.enterWell(it, it.state.view || 0);
    // A living culture: show what the engine says the dish looks like now.
    if (it.state.detached) return this.scope.open({ kind: "wait", title: it.name, text: "The cells are floating in suspension now. Count them in the hemocytometer." });
    this.scope.open({ kind: "wait", title: it.name, text: this.engine.version ? "Simulating the cells…" : "Starting the CellSim engine (the first time downloads Python, about 30 s)…" });
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
  // One well of the 6-well plate, simulated on its own.
  async enterWell(plate, i) {
    plate.state.view = i;
    const st = plate.state, ws = st.wells[i];
    const conc = ws.events.length ? ws.events[ws.events.length - 1].conc : {};
    const label = Object.keys(conc).length ? Object.entries(conc).map(([k, v]) => `${k} ${fmtConc(v)}`).join(", ") : "no drug";
    const base = { title: `${plate.name} · well ${P.WELLS[i]} · ${label}`, wells: P.WELLS, well: i };
    this.scope.open({ kind: "wait", ...base, text: this.engine.version ? "Simulating the cells in this well…" : "Starting the CellSim engine…" });
    try {
      if (!ws.result || ws.result.t !== st.t) {
        const out = await this.runWell(plate, i);
        ws.result = { header: out.header, snap: out.snap, t: st.t };
      }
      if (this.screen !== "scope" || this.stageItem !== plate || st.view !== i) return;
      this.scope.open({ kind: "culture", ...base, header: ws.result.header, snap: ws.result.snap, time: st.t, tint: "#e6d6d9", keepView: true });
    } catch (e) {
      this.scope.open({ kind: "wait", ...base, text: `The engine could not run: ${e.message}` });
    }
  }
  runWell(plate, i) {
    const st = plate.state, ws = st.wells[i];
    const drugs = [...new Set(ws.events.flatMap((e) => Object.keys(e.conc)))];
    return this.engine.dish({ label: `well ${P.WELLS[i]}`, line: st.line, drugs, events: ws.events, t: st.t, n_seed: 40, grid: 40, seed: st.seed, every: 12 });
  }
  // The hemocytometer: cells from a dish's suspension, at the engine's count.
  async enterHemo(it) {
    const c = it.container;
    const title = "Hemocytometer";
    const dish = this.items.find((d) => d.id === c.meta?.dish);
    if (c.empty || !dish || !c.c.cells) return this.scope.open({ kind: "wait", title, text: c.empty ? "The chamber is empty. Load 10 µL of cells mixed with trypan blue first." : "There are no cells in this liquid. Lift cells off a dish with trypsin first." });
    this.scope.open({ kind: "wait", title, text: this.engine.version ? "Getting the dish's cells from the engine…" : "Starting the CellSim engine (the first time downloads Python, about 30 s)…" });
    try {
      await this.simulate(dish);
    } catch (e) { return this.scope.open({ kind: "wait", title, text: `The engine could not run: ${e.message}` }); }
    const last = dish.state.series[dish.state.series.length - 1];
    const h = dish.state.header, patch = (h.grid_shape[0] * h.spacing_um) * (h.grid_shape[1] * h.spacing_um);
    const scale = DISH_AREA_UM2 / patch;   // the engine simulates one patch of the dish
    const live = last.n_live * scale, dead = (last.n_apoptotic + last.n_necrotic) * scale + live * 0.04;   // trypsin kills a few
    const frac = c.c.cells / c.volume;   // 1 mL of suspension, diluted in the chamber liquid
    it.state.counted = it.state.counted || new Set();
    const hemo = { livePerMl: live * frac, deadPerMl: dead * frac, trypan: !!c.c.trypan, dilution: 1 / frac, t: dish.state.t, dish: dish.name, seed: 7 + dish.state.t, counted: it.state.counted };
    if (this.screen === "scope" && this.stageItem === it) this.scope.open({ kind: "hemo", title: `Hemocytometer · ${dish.name} cells at ${dish.state.t} h`, hemo });
  }
  exitScope() {
    this.screen = "play";
    $("scope").hidden = true;
    $("hud").hidden = false;
    this.lock();
    this.dirty = true;
  }
  onFocus(obj, spec) {
    if (spec.kind === "slide" || spec.kind === "hemo") this.emit(`focus:${obj}`);
    else {
      const d = this.stageItem;
      if (obj >= 10 && d && d.kind === "dish") {
        this.emit("observe:culture");
        if (d.state.t > 0 && d.state.events.length) this.emit("observe:after");
        if (d.name === "Dish B" && d.state.t > 0) this.emit("compare");
      }
    }
  }
  onScopeEvent(name, data) {
    const it = this.stageItem;
    if (name === "well" && it?.kind === "plate") this.enterWell(it, data);
    else if (name === "classify") {
      this.journal.result("mitosis", { counts: data.counts, total: data.total, right: data.right });
      if (data.total >= 25) this.emit("classify:25");
      if (data.counts.slice(1).every((n) => n > 0)) this.emit("stages:all");
    } else if (name === "count:4") this.emit("count:4");
    else if (name === "count") {
      const e = (n) => `${(n / 10 ** Math.floor(Math.log10(n || 1))).toFixed(1)} × 10^${Math.floor(Math.log10(n || 1))}`;
      this.journal.result("count", { text: `${data.dish} at ${data.t} h: you counted ${data.live} live and ${data.dead} dead cells in 4 squares, so <b>${e(data.perMl)} live cells per mL</b> and ${data.viab.toFixed(0)} % viability. The engine's dish really held ${e(data.truePerMl)} live cells per mL.` });
      this.emit("count:done");
    } else if (name === "osmo" && it?.kind === "slide") {
      const st = it.state, last = st.osmo[st.osmo.length - 1];
      const log = (this.journal.data.results.osmosis?.log || []).slice(-8);
      const where = `${P.shortName(st.mediumName || "water")} (${last ? last.osm.toFixed(2) : "?"} osmol/L)`;
      if (data >= 0.6 && !st.plas) { st.plas = true; log.push(`In ${where}: ${(data * 100).toFixed(0)} % of the cells in view plasmolysed.`); this.journal.result("osmosis", { log }); this.emit("osmo:plasmolysed"); }
      if (st.plas && data <= 0.1 && last && last.osm < 0.3) { st.plas = false; log.push(`Back in ${where}: the cells swelled against their walls again.`); this.journal.result("osmosis", { log }); this.emit("osmo:recovered"); }
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
    this.toast(`${d.name} now has ${Object.entries(conc).map(([k, v]) => `${k} ${fmtConc(v)}`).join(", ")}`);
  }
  wellDosed(plate, i) {
    const st = plate.state, ws = st.wells[i], conc = plate.wells[i].container.drugs();
    const ev = ws.events.find((e) => e.t === st.t);
    if (ev) ev.conc = conc; else ws.events.push({ t: st.t, conc });
    ws.result = null; st.read = null;
    const treated = st.wells.filter((w) => w.events.length);
    const doses = new Set(treated.map((w) => Object.values(w.events[w.events.length - 1].conc).reduce((a, b) => a + b, 0).toPrecision(2)));
    if (doses.size >= 3) this.emit("dose:3");
    if (doses.size >= 3 && treated.length < 6) this.emit("control");
    if (treated.length === 6) this.toast("Every well has a drug now: with no untreated control, there's nothing to compare against.", "warn");
    this.toast(`Well ${P.WELLS[i]}: ${Object.entries(conc).map(([k, v]) => `${k} ${fmtConc(v)}`).join(", ")}`);
  }
  putInIncubator(d) {
    this.inv.remove(d);
    d.obj.removeFromParent();
    this.incubator.dishes.push(d);
    this.afterAction();
    this.toast(`${d.name} is in the incubator. Right-click it to let time pass.`);
  }
  openIncubator() { this.renderIncubator(); this.openUI("incubatorPanel"); }
  renderIncubator() {
    const ds = this.incubator.dishes;
    const what = (d) => d.kind === "plate" ? `${d.state.wells.filter((w) => w.events.length).length} of 6 wells treated` : d.state.detached ? "cells lifted off" : d.state.events.length ? "treated" : "no drug";
    $("incubatorBody").innerHTML = (ds.length ? ds.map((d) => `<div class="dish-row"><b>${d.name}</b><span class="note">${d.state.t} h in culture · ${what(d)}</span>
        <div class="btns"><button class="mbtn small" data-inc="take" data-id="${d.id}">Take it out</button></div></div>`).join("") : `<p class="note">The incubator is empty.</p>`) +
      (ds.length ? `<p class="glabel">Let time pass for everything inside</p><div class="btns" style="display:flex;gap:6px;flex-wrap:wrap">${[6, 12, 24, 72].map((h) => `<button class="mbtn small primary" data-inc="wait" data-h="${h}">Wait ${h} h</button>`).join("")}</div>` : "") +
      `<p class="note">${this.engine.version ? `CellSim ${this.engine.version} simulates each dish and well.` : "The CellSim engine starts the first time you incubate or look at a dish."}</p>`;
  }
  async incubate(h) {
    const ds = [...this.incubator.dishes];
    for (const d of ds) { d.state.t += h; d.state.needs = true; }
    this.renderIncubator();
    if (ds.some((d) => d.kind === "dish" && d.state.events.length)) this.emit("incubate");
    if (ds.some((d) => d.kind === "plate" && d.state.t >= 72 && d.state.wells.some((w) => w.events.length))) this.emit("incubate:72");
    try {
      for (const d of ds) {
        if (d.kind !== "dish" || d.state.detached) continue;
        await this.simulate(d);
        const last = d.state.series[d.state.series.length - 1];
        this.toast(`${d.name} at ${d.state.t} h: ${last.n_live} live cells, ${last.n_apoptotic} dying`, "good");
      }
      for (const d of ds) if (d.kind === "plate") this.toast(`${d.name} at ${d.state.t} h. Read it in the plate reader to see how each well did.`, "good");
    } catch (e) { this.toast(`The engine could not run: ${e.message}`, "warn"); }
    this.renderIncubator();
  }
  // The plate reader: every well simulated to now, live cells relative to the untreated wells.
  async readPlate(plate, reader) {
    if (this.reading) return;
    const show = reader.obj.userData.screen, st = plate.state;
    if (st.t < 24) this.toast("The drug has barely had time to act: 72 h in the incubator gives a clear result.", "warn");
    this.reading = true;
    this.inv.remove(plate);
    plate.obj.removeFromParent();
    plate.obj.userData.pose?.("rest");
    plate.obj.position.copy(reader.obj.localToWorld(new THREE.Vector3(-0.06, 0.044, 0.19)));
    plate.obj.rotation.set(0, reader.obj.rotation.y, 0);
    enableShadows(plate.obj);
    if (plate.obj.userData.blob) plate.obj.userData.blob.visible = true;
    this.scene.add(plate.obj);
    this.afterAction();
    sfx.click();
    try {
      const wells = [];
      for (let i = 0; i < 6; i++) {
        show([`READING ${P.WELLS[i]}`, `${i + 1} of 6 wells`], "#f0b44c");
        this.dirty = true;
        const out = await this.runWell(plate, i);
        st.wells[i].result = { header: out.header, snap: out.snap, t: st.t };
        const last = out.series[out.series.length - 1], ev = st.wells[i].events;
        const conc = ev.length ? ev[ev.length - 1].conc : {};
        const dose = Object.values(conc).reduce((a, b) => a + b, 0);
        wells.push({ live: last.n_live, dose, drug: Object.keys(conc).filter((k) => conc[k] > 0).join(" + ") || "none" });
      }
      const ctrl = wells.filter((w) => w.dose <= 0);
      const base = ctrl.length ? ctrl.reduce((a, w) => a + w.live, 0) / ctrl.length : Math.max(...wells.map((w) => w.live));
      for (const w of wells) w.viability = (w.live / Math.max(1, base)) * 100;
      const treated = wells.filter((w) => w.dose > 0).sort((a, b) => a.dose - b.dose);
      const drugs = [...new Set(treated.map((w) => w.drug))];
      let ic50 = null, note = "";
      for (let k = 0; k + 1 < treated.length; k++) {
        const a = treated[k], b = treated[k + 1];
        if (a.viability >= 50 && b.viability < 50) {
          const la = Math.log10(a.dose), lb = Math.log10(b.dose);
          ic50 = 10 ** (la + ((a.viability - 50) * (lb - la)) / (a.viability - b.viability));
          break;
        }
      }
      if (!ic50) note = !treated.length ? "No well had a drug in it." : treated.every((w) => w.viability >= 50) ? "Every dose left more than half the cells alive: try higher doses (a bigger volume)." : treated[0].viability < 50 ? "Even the lowest dose killed more than half the cells: try smaller volumes." : "The curve doesn't cross 50 % cleanly. Try doses spread over a wider range.";
      if (!ctrl.length) note += " There was no untreated well, so viability is relative to the best well.";
      if (drugs.length > 1) note += " The wells had different drugs: use one drug for a clean curve.";
      st.read = { t: st.t, wells, ic50, drug: drugs.join(", ") || "none", note: note.trim() };
      show(["DONE", ic50 ? `IC50 ${fmtConc(ic50)}` : "See the journal"], "#7dffb0");
      this.emit("read");
      this.toast(`Plate read. Open the journal (J), Cell culture tab, for the dose–response curve.`, "good");
    } catch (e) {
      show(["ERROR", "Engine failed"], "#ff7466");
      this.toast(`The engine could not run: ${e.message}`, "warn");
    }
    this.reading = false;
    this.dirty = true;
  }
  openJournal() { this.journalTab = this.journalTab || "tasks"; this.renderJournal(); this.openUI("journal"); }
  renderJournal() {
    for (const b of document.querySelectorAll("[data-tab]")) b.setAttribute("aria-selected", String(b.dataset.tab === this.journalTab));
    $("journalBody").innerHTML = this.journal.render(this.journalTab, this.tasks, this.dishes, this.plate);
    if (this.journalTab === "culture" && this.plate.state.read) this.emit("dr-curve");
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
    this.tracker.stop();
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
      this.shadowsDirty = true;
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
      if (!live) this.actions = { target: "" };
      this.updateAimCues();
      const sig = live ? `${this.actions.left?.label}|${this.actions.right?.label}|${this.actions.target}` : "unlocked";
      if (sig !== this.promptSig) { this.promptSig = sig; HUD.prompt(live ? this.actions : { target: "Click to look around" }); }
      if (this.tracker.target) {
        const t = this.tracker.target;
        if (this.tracker.reached(this.aim)) { this.toast(`Found it: the ${sentence(t.name)}`, "good"); this.tracker.stop(); }
        else if (t.item && !t.item.obj.parent) this.tracker.stop();
      }
    }
    if (this.tracker.target) this.tracker.update(this.camera);
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
    if (this.shadowsDirty && this.preset.shadows) this.renderer.shadowMap.needsUpdate = true;
    this.shadowsDirty = false;
    this.renderer.info.reset();
    this.renderer.render(this.scene, this.camera);
    if (this.screen === "play") this.outliner.draw(this.scene, this.camera);
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
function esc(s) { return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]); }
function capital(s) { return s.charAt(0).toUpperCase() + s.slice(1); }

const SPARES = [
  { label: "a new box of slides", make: () => { const it = new Item("slidebox", "Box of slides", P.makeBox({ w: 0.085, h: 0.03, d: 0.035, color: 0x2a5fb0, labelLines: ["Microscope slides", "76 × 26 mm"] })); it.state.left = 10; return it; } },
  { label: "a new pack of cover slips", make: () => { const it = new Item("coverbox", "Cover slips", P.makeBox({ w: 0.032, h: 0.014, d: 0.032, color: 0xdfe6ea, labelLines: ["Cover slips"] })); it.state.left = 20; return it; } },
  { label: "a fresh 30 mL bottle of iodine", make: () => { const m = P.makeDropperBottle("iodine", 30); return new Item("dropper", "Iodine dropper bottle", m.obj, { container: m.container, liquid: m.liquid, solution: "iodine" }); } },
  { label: "a fresh 30 mL bottle of methylene blue", make: () => { const m = P.makeDropperBottle("mb", 30); return new Item("dropper", "Methylene blue dropper bottle", m.obj, { container: m.container, liquid: m.liquid, solution: "mb" }); } },
  { label: "a fresh bottle of toluidine blue", make: () => { const m = P.makeDropperBottle("toluidine", 30); return new Item("dropper", "Toluidine blue dropper bottle", m.obj, { container: m.container, liquid: m.liquid, solution: "toluidine" }); } },
];

const CREDITS = `
<p><b>Cell engine:</b> <a href="https://github.com/Rockman6/CellSim" target="_blank" rel="noopener">CellSim</a> (MIT), running in your browser with <a href="https://pyodide.org" target="_blank" rel="noopener">Pyodide</a>.</p>
<p><b>3D engine:</b> <a href="https://threejs.org" target="_blank" rel="noopener">three.js</a> r160 (MIT).</p>
<p><b>Models, textures and lighting from Poly Haven (CC0):</b> yellow onion, metal stool, steel shelves, fluorescent lights, drawer cabinet, clipboard, stationery, wall clock, magnifying glass, binder notebook, metal office desk, floor tiles, painted plaster, hospital room HDRI.</p>
<p><b>Models from Sketchfab (CC Attribution 4.0):</b>
"Chemistry Lab Table" by Jawahar Yokesh · "lab bench" by Konstantin_Keller (the incubator) · "Compound Microscope" by kg007 · "Door Classroom" by Mehdi Shahsavan.</p>
<p><b>Built in code:</b> the room, slides, cover slips, bottles, beakers, culture dishes, the 6-well plate, the micropipette and tips, microtubes, the hemocytometer, the plate reader, the water bath, the pondweed tank, the lab tools, the safety goggles, the hand sanitizer and the window blinds. The sound effects are made in code too.</p>
<p><b>Not simulated:</b> the onion, root tip, red onion, cheek and pondweed slides are drawn at measured sizes and stage proportions; osmosis follows a simple osmotic model. The living A549 cells are CellSim's.</p>`;

new Game().boot().catch((e) => {
  console.error(e);
  $("loadText").textContent = `Something went wrong while loading: ${e.message}`;
});
