/* Runs the CellSim engine off the main thread, so the lab keeps moving.
   Pyodide loads CPython with NumPy and SciPy, then micropip installs the
   cellsim wheel built from this repository (the same one ../lab/ uses). */
const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.26.2/full/";
let py = null;

const say = (text) => postMessage({ type: "status", text });

// A missing file is an error, never its 404 page read as the file. The site
// moved when the repository was renamed, so an old link lands here.
async function getText(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`the engine file ${url.pathname.split("/").slice(-2).join("/")} is missing at this address (HTTP ${r.status}). If you opened an old link, the lab is now at rockman6.github.io/CellSim/lab3d/`);
  return r.text();
}

// One readable line: the last line of a Python traceback says what failed.
const brief = (e) => String((e && e.message) || e).trim().split("\n").filter(Boolean).pop()?.slice(0, 300) || "unknown error";

async function boot() {
  say("Downloading Python…");
  importScripts(PYODIDE + "pyodide.js");
  py = await loadPyodide({ indexURL: PYODIDE });
  say("Loading NumPy and SciPy…");
  await py.loadPackage(["numpy", "scipy", "micropip"]);
  const base = new URL("../../lab/wheel/", self.location.href);
  const wheel = (await getText(new URL("LATEST", base))).trim();
  if (!/^cellsim-[\w.]+-py3-none-any\.whl$/.test(wheel)) throw new Error(`no engine wheel is named in lab/wheel/LATEST (it reads "${wheel.slice(0, 60)}")`);
  say(`Installing ${wheel}…`);
  await py.runPythonAsync(`import micropip\nawait micropip.install(${JSON.stringify(new URL(wheel, base).href)})`);
  const bridge = await getText(new URL("lab_engine.py", self.location.href));
  await py.runPythonAsync(bridge);
  const version = py.runPython("import importlib.metadata as m; m.version('cellsim')");
  postMessage({ type: "ready", version });
}

// Python can be up while the engine is not, so readiness is its own flag.
let ready = false;
const booting = boot().then(() => (ready = true), (e) => postMessage({ type: "fatal", message: brief(e) }));

onmessage = async (ev) => {
  const { id, cmd, cfg } = ev.data;
  await booting;
  if (!ready) return postMessage({ type: "error", id, message: "The engine did not start." });
  try {
    if (cmd === "dish") {
      py.globals.set("_cfg", JSON.stringify(cfg));
      const out = py.runPython("lab_run_dish(_cfg)");
      postMessage({ type: "result", id, data: out });
    }
  } catch (e) {
    postMessage({ type: "error", id, message: brief(e) });
  }
};
