/* Runs the Cell Sim engine off the main thread, so the lab keeps moving.
   Pyodide loads CPython with NumPy and SciPy, then micropip installs the
   cellsim wheel built from this repository (the same one ../lab/ uses). */
const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.26.2/full/";
let py = null;

const say = (text) => postMessage({ type: "status", text });

async function boot() {
  say("Downloading Python…");
  importScripts(PYODIDE + "pyodide.js");
  py = await loadPyodide({ indexURL: PYODIDE });
  say("Loading NumPy and SciPy…");
  await py.loadPackage(["numpy", "scipy", "micropip"]);
  const base = new URL("../../lab/wheel/", self.location.href);
  const wheel = (await (await fetch(new URL("LATEST", base))).text()).trim();
  if (!wheel.endsWith(".whl")) throw new Error(`no engine wheel found (LATEST = "${wheel}")`);
  say(`Installing ${wheel}…`);
  await py.runPythonAsync(`import micropip\nawait micropip.install(${JSON.stringify(new URL(wheel, base).href)})`);
  const bridge = await (await fetch(new URL("lab_engine.py", self.location.href))).text();
  await py.runPythonAsync(bridge);
  const version = py.runPython("import importlib.metadata as m; m.version('cellsim')");
  postMessage({ type: "ready", version });
}

const booting = boot().catch((e) => postMessage({ type: "fatal", message: String(e && e.message || e) }));

onmessage = async (ev) => {
  const { id, cmd, cfg } = ev.data;
  await booting;
  if (!py) return postMessage({ type: "error", id, message: "The engine did not start." });
  try {
    if (cmd === "dish") {
      py.globals.set("_cfg", JSON.stringify(cfg));
      const out = py.runPython("lab_run_dish(_cfg)");
      postMessage({ type: "result", id, data: out });
    }
  } catch (e) {
    postMessage({ type: "error", id, message: String(e && e.message || e).split("\n").slice(-3).join(" ") });
  }
};
