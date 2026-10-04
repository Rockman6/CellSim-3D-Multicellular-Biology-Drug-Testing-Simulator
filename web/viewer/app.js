/* CellSim dish viewer: renders a cellsim stream; runs no biology. */
(function () {
  "use strict";
  const $ = (id) => document.getElementById(id);
  const S = window.CellSimStream;
  const Charts = window.CellSimCharts;

  const FIELDS = {
    phase: { label: "cycle phase", kind: "phase" },
    p53: { label: "p53 (normalised)", range: [0, 0.6] },
    caspase3: { label: "caspase-3 activity", range: [0, 1] },
    damage: { label: "DNA damage index", range: [0, 1] },
    o2_mmHg: { label: "oxygen, mmHg", range: [0, 100] },
    growth_signal: { label: "G1 growth signal", range: [0, 1] },
    lineage: { label: "lineage (clone)", kind: "lineage" },
  };
  // perceptually ordered ramp (viridis stops)
  const RAMP = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];

  const state = { model: null, frame: 0, playing: false, timer: null, shown: [] };
  let renderer, scene, camera, liveMesh, apoMesh, necMesh, light;
  const orbit = { theta: 0.6, phi: 1.0, dist: 800, target: null };
  const tmpM = typeof THREE !== "undefined" ? new THREE.Matrix4() : null;
  const tmpC = typeof THREE !== "undefined" ? new THREE.Color() : null;

  function css(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || "#888";
  }

  // ── loading ──
  async function loadBytes(buffer, name) {
    try {
      const text = await S.bytesToText(buffer);
      const model = S.normalize(S.parseJsonl(text));
      model.name = name;
      show(model);
    } catch (err) {
      $("meta").textContent = `Could not read ${name}: ${err.message}`;
    }
  }

  async function loadUrl(url) {
    $("meta").textContent = `Loading ${url} …`;
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      await loadBytes(await res.arrayBuffer(), url.split("/").pop());
    } catch (err) {
      $("meta").textContent = `Could not load ${url}: ${err.message}`;
    }
  }

  async function loadDemoList() {
    try {
      const res = await fetch("../demo/index.json");
      if (!res.ok) return;
      const list = await res.json();
      for (const d of list) {
        const o = document.createElement("option");
        o.value = "../demo/" + d.file;
        o.textContent = d.label;
        $("demo").appendChild(o);
      }
    } catch (_) { /* no demos alongside this copy */ }
  }

  // ── scene ──
  function ensureRenderer() {
    if (renderer || typeof THREE === "undefined") return;
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    renderer.setClearColor(new THREE.Color(css("--stage")), 1);
    $("stage").appendChild(renderer.domElement);
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(40, 1, 1, 100000);
    scene.add(new THREE.AmbientLight(0xffffff, 0.55));
    light = new THREE.DirectionalLight(0xffffff, 0.75);
    scene.add(light);
    new ResizeObserver(resize).observe($("stage-wrap"));
    resize();
    attachOrbit(renderer.domElement);
  }

  function resize() {
    if (!renderer) return;
    const r = $("stage-wrap").getBoundingClientRect();
    renderer.setSize(r.width, r.height, false);
    camera.aspect = r.width / Math.max(1, r.height);
    camera.updateProjectionMatrix();
    renderOnce();
  }

  function placeCamera() {
    const t = orbit.target || new THREE.Vector3();
    const sp = Math.sin(orbit.phi);
    camera.position.set(t.x + orbit.dist * sp * Math.cos(orbit.theta),
      t.y + orbit.dist * sp * Math.sin(orbit.theta), t.z + orbit.dist * Math.cos(orbit.phi));
    camera.up.set(0, 0, 1);
    camera.lookAt(t);
    light.position.copy(camera.position);
  }

  function attachOrbit(canvas) {
    let drag = null;
    canvas.addEventListener("pointerdown", (e) => { drag = { x: e.clientX, y: e.clientY, moved: false }; });
    window.addEventListener("pointerup", (e) => {
      if (drag && !drag.moved) pick(e);
      drag = null;
    });
    window.addEventListener("pointermove", (e) => {
      if (!drag) return;
      const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      if (Math.abs(dx) + Math.abs(dy) > 2) drag.moved = true;
      drag.x = e.clientX; drag.y = e.clientY;
      orbit.theta -= dx * 0.008;
      orbit.phi = Math.min(Math.PI - 0.05, Math.max(0.05, orbit.phi - dy * 0.008));
      renderOnce();
    });
    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      orbit.dist *= Math.exp(e.deltaY * 0.001);
      renderOnce();
    }, { passive: false });
  }

  function buildMeshes(model) {
    for (const m of [liveMesh, apoMesh, necMesh]) if (m) { scene.remove(m); m.geometry.dispose(); }
    liveMesh = apoMesh = necMesh = null;
    if (model.kind !== "dish") return;
    const h = model.header;
    const r = 0.45 * h.spacing_um;
    let maxLive = 1, maxApo = 1, maxNec = 1;
    for (const s of model.snaps) {
      maxLive = Math.max(maxLive, s.cells.x_um.length);
      maxApo = Math.max(maxApo, s.apoptotic_um.length);
      maxNec = Math.max(maxNec, s.necrotic_um.length);
    }
    const detail = maxLive > 20000 ? 0 : 1;
    const mat = new THREE.MeshLambertMaterial({ color: 0xffffff });
    liveMesh = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(r, detail), mat, maxLive);
    apoMesh = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(0.55 * r, 0),
      new THREE.MeshLambertMaterial({ color: new THREE.Color(css("--c-apo")) }), maxApo);
    necMesh = new THREE.InstancedMesh(new THREE.IcosahedronGeometry(0.8 * r, 0),
      new THREE.MeshLambertMaterial({ color: new THREE.Color(css("--c-nec")), transparent: true, opacity: 0.55 }), maxNec);
    scene.add(liveMesh, apoMesh, necMesh);
    // Frame the cells, not the lattice: centre on the largest snapshot's
    // cells and stand back far enough to see all of them.
    const big = model.snaps.reduce((a, b) => (b.cells.x_um.length > a.cells.x_um.length ? b : a));
    const c = big.cells, n = Math.max(1, c.x_um.length);
    const mean = (xs) => xs.reduce((a, b) => a + b, 0) / n;
    const cx = mean(c.x_um), cy = mean(c.y_um), cz = mean(c.z_um);
    let reach = h.spacing_um;
    for (let i = 0; i < c.x_um.length; i++) {
      reach = Math.max(reach, Math.hypot(c.x_um[i] - cx, c.y_um[i] - cy, c.z_um[i] - cz));
    }
    orbit.target = new THREE.Vector3(cx, cy, cz);
    orbit.dist = 3.2 * reach;
    orbit.phi = h.geometry === "monolayer" ? 0.55 : 1.1;
  }

  function colourFor(field, cells, i) {
    if (field === "phase") return state.phaseColours[cells.phase[i]];
    if (field === "lineage") {
      const hue = ((cells.lineage[i] * 0.618033988749895) % 1);
      return tmpC.setHSL(hue, 0.65, 0.55).clone();
    }
    const [lo, hi] = FIELDS[field].range;
    const v = Math.min(1, Math.max(0, (cells[field][i] - lo) / (hi - lo)));
    const k = v * (RAMP.length - 1), a = Math.floor(k), b = Math.min(RAMP.length - 1, a + 1), f = k - a;
    const c = RAMP[a].map((x, j) => (x + (RAMP[b][j] - x) * f) / 255);
    return new THREE.Color(c[0], c[1], c[2]);
  }

  function renderFrame(i) {
    const model = state.model;
    if (!model) return;
    state.frame = i;
    const snap = model.snaps[i];
    $("time").value = String(i);
    $("tlabel").textContent = `t = ${snap.t_h.toFixed(1)} h`;
    if (model.kind === "dish" && liveMesh) {
      const field = $("color").value;
      const cut = $("slice").checked;
      const zc = orbit.target ? orbit.target.z : 0;
      const cells = snap.cells;
      const shown = [];
      for (let k = 0; k < cells.x_um.length; k++) {
        if (cut && cells.z_um[k] > zc + 1e-6) continue;
        shown.push(k);
      }
      shown.forEach((k, j) => {
        tmpM.makeTranslation(cells.x_um[k], cells.y_um[k], cells.z_um[k]);
        liveMesh.setMatrixAt(j, tmpM);
        liveMesh.setColorAt(j, colourFor(field, cells, k));
      });
      liveMesh.count = shown.length;
      liveMesh.instanceMatrix.needsUpdate = true;
      if (liveMesh.instanceColor) liveMesh.instanceColor.needsUpdate = true;
      state.shown = shown;
      fillBodies(apoMesh, snap.apoptotic_um, cut, zc);
      fillBodies(necMesh, snap.necrotic_um, cut, zc);
      legend(field);
    }
    drawCharts();
    renderOnce();
  }

  function fillBodies(mesh, pts, cut, zc) {
    let j = 0;
    for (const p of pts) {
      if (cut && p[2] > zc + 1e-6) continue;
      tmpM.makeTranslation(p[0], p[1], p[2]);
      mesh.setMatrixAt(j++, tmpM);
    }
    mesh.count = j;
    mesh.instanceMatrix.needsUpdate = true;
  }

  function renderOnce() {
    if (!renderer || !scene) return;
    placeCamera();
    renderer.render(scene, camera);
  }

  function legend(field) {
    const box = $("legend");
    box.hidden = false;
    const f = FIELDS[field];
    if (f.kind === "phase") {
      box.innerHTML = S.PHASES.map((p, i) =>
        `<span class="swatch" style="background:${state.phaseCss[i]}"></span>${p}`).join("&nbsp;&nbsp;") +
        `&nbsp;&nbsp;<span class="swatch" style="background:${css("--c-apo")}"></span>apoptotic` +
        `&nbsp;&nbsp;<span class="swatch" style="background:${css("--c-nec")}"></span>necrotic`;
    } else if (f.kind === "lineage") {
      box.textContent = "Each colour is one founder cell and all its descendants.";
    } else {
      const stops = RAMP.map((c) => `rgb(${c.join(",")})`).join(",");
      box.innerHTML = `${f.label}<div class="ramp" style="background:linear-gradient(90deg,${stops})"></div>` +
        `<span class="mono">${f.range[0]}</span><span class="mono" style="float:right">${f.range[1]}</span>`;
    }
  }

  // ── picking ──
  function pick(e) {
    if (!state.model || state.model.kind !== "dish" || !liveMesh) return;
    const rect = renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1,
      -((e.clientY - rect.top) / rect.height) * 2 + 1);
    const ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, camera);
    const hit = ray.intersectObject(liveMesh)[0];
    const box = $("inspect");
    if (!hit || hit.instanceId === undefined) { box.hidden = true; return; }
    const k = state.shown[hit.instanceId];
    const c = state.model.snaps[state.frame].cells;
    const rows = [["cell", k], ["phase", S.PHASES[c.phase[k]]], ["lineage", c.lineage[k]],
      ["p53", c.p53[k]], ["caspase-3", c.caspase3[k]], ["DNA damage", c.damage[k]],
      ["oxygen (mmHg)", c.o2_mmHg[k]], ["growth signal", c.growth_signal[k]]];
    box.innerHTML = "<table>" + rows.map(([a, b]) =>
      `<tr><td class="muted">${a}</td><td class="mono">${typeof b === "number" ? +b.toFixed(3) : b}</td></tr>`).join("") +
      "</table>";
    box.hidden = false;
  }

  // ── charts ──
  function drawCharts() {
    const m = state.model;
    if (!m) return;
    const s = m.series, cur = s.t[state.frame];
    if (m.kind === "dish") {
      Charts.draw($("chart-counts"), { t: s.t, cursor: cur, series: [
        { values: s.live, color: css("--c-live"), label: "live" },
        { values: s.apoptotic, color: css("--c-apo"), label: "apoptotic" },
        { values: s.necrotic, color: css("--c-nec"), label: "necrotic" }] });
      if (m.header.geometry === "spheroid") {
        $("third-caption").textContent = "Radius and necrotic core (µm)";
        Charts.draw($("chart-third"), { t: s.t, cursor: cur, series: [
          { values: s.radius, color: css("--c-live"), label: "radius" },
          { values: s.necroticRadius, color: css("--c-nec"), label: "necrotic core" }] });
        const prof = m.snaps[state.frame].profile;
        $("fig-profile").hidden = !prof;
        if (prof) {
          const dose = Math.max(1e-12, (m.snaps[state.frame].drug_uM || [0])[0]);
          const o2Medium = (m.header.dish_params || {}).o2_medium_mmHg || 100;
          const series = [{ values: prof.o2_mmHg.map((v) => 100 * v / o2Medium),
            color: css("--c-o2"), label: "oxygen, % of medium" }];
          if (prof.drug_uM && prof.drug_uM[0] && prof.drug_uM[0].length && dose > 1e-12) {
            series.push({ values: prof.drug_uM[0].map((v) => 100 * v / dose), color: css("--c-drug"), label: "drug, % of medium" });
          }
          Charts.draw($("chart-profile"), { t: prof.r_um, series, yMin: 0, yMax: 100, xUnit: "µm" });
        }
      } else {
        $("third-caption").textContent = "Confluence and quiescent fraction";
        $("fig-profile").hidden = true;
        Charts.draw($("chart-third"), { t: s.t, cursor: cur, yMin: 0, yMax: 1, series: [
          { values: s.confluence, color: css("--c-live"), label: "confluence" },
          { values: s.quiescent, color: css("--c-nec"), label: "quiescent" }] });
      }
    } else {
      Charts.draw($("chart-counts"), { t: s.t, cursor: cur, series: [
        { values: s.fold, color: css("--c-live"), label: "fold change" }] });
      $("third-caption").textContent = "Mean state of live cells";
      $("fig-profile").hidden = true;
      Charts.draw($("chart-third"), { t: s.t, cursor: cur, series: [
        { values: s.p53, color: css("--c-g1"), label: "p53" },
        { values: s.c3, color: css("--c-m"), label: "caspase-3" },
        { values: s.damage, color: css("--c-g2"), label: "damage" }] });
    }
    Charts.draw($("chart-phase"), { t: s.t, cursor: cur, yMin: 0, yMax: 1, series: S.PHASES.map((p, i) =>
      ({ values: s.phase[i], color: state.phaseCss[i], label: p })) });
  }

  // ── wiring ──
  function show(model) {
    state.model = model;
    state.phaseCss = ["--c-g1", "--c-s", "--c-g2", "--c-m"].map(css);
    if (typeof THREE !== "undefined") state.phaseColours = state.phaseCss.map((c) => new THREE.Color(c));
    $("drop").hidden = true;
    const h = model.header;
    const drug = Array.isArray(h.drug) ? h.drug.join(" + ") : (h.drug || "no drug");
    const geo = model.kind === "dish" ? `${h.geometry}, ` : "well-mixed population, ";
    $("meta").innerHTML = `<strong>${h.line}</strong> · ${drug} · ${geo}${model.snaps.length} snapshots to ` +
      `${model.series.t[model.series.t.length - 1]} h <span class="muted">(${model.name})</span>`;
    for (const id of ["play", "time", "csv"]) $(id).disabled = false;
    $("color").disabled = model.kind !== "dish";
    $("slice").disabled = !(model.kind === "dish" && h.geometry === "spheroid");
    $("slice").checked = model.kind === "dish" && h.geometry === "spheroid";
    $("time").max = String(model.snaps.length - 1);
    $("legend").hidden = model.kind !== "dish";
    $("inspect").hidden = true;
    if (model.kind === "dish") { ensureRenderer(); buildMeshes(model); }
    renderFrame(0);
  }

  function togglePlay() {
    const m = state.model;
    if (!m) return;
    state.playing = !state.playing;
    $("play").textContent = state.playing ? "Pause" : "Play";
    clearInterval(state.timer);
    if (state.playing) {
      if (state.frame >= m.snaps.length - 1) state.frame = -1;
      state.timer = setInterval(() => {
        if (state.frame >= m.snaps.length - 1) { togglePlay(); return; }
        renderFrame(state.frame + 1);
      }, 220);
    }
  }

  function downloadCsv() {
    const m = state.model;
    if (!m) return;
    const blob = new Blob([S.toCSV(m)], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = (m.name || "cellsim").replace(/\.jsonl(\.gz)?$/, "") + "-timeseries.csv";
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function init() {
    $("file").addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (f) await loadBytes(await f.arrayBuffer(), f.name);
    });
    $("demo").addEventListener("change", (e) => { if (e.target.value) loadUrl(e.target.value); });
    $("play").addEventListener("click", togglePlay);
    $("time").addEventListener("input", (e) => renderFrame(+e.target.value));
    $("color").addEventListener("change", () => renderFrame(state.frame));
    $("slice").addEventListener("change", () => renderFrame(state.frame));
    $("csv").addEventListener("click", downloadCsv);
    const wrap = $("stage-wrap");
    wrap.addEventListener("dragover", (e) => { e.preventDefault(); $("drop").classList.add("hover"); });
    wrap.addEventListener("dragleave", () => $("drop").classList.remove("hover"));
    wrap.addEventListener("drop", async (e) => {
      e.preventDefault();
      $("drop").classList.remove("hover");
      const f = e.dataTransfer.files[0];
      if (f) await loadBytes(await f.arrayBuffer(), f.name);
    });
    window.addEventListener("resize", () => state.model && drawCharts());
    loadDemoList();
    const src = new URLSearchParams(location.search).get("src");
    if (src) loadUrl(src);
  }

  init();
})();
