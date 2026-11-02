// Player settings and graphics presets. Low is meant for older Macs and
// laptops on battery; Ultra for a fast graphics card.
const KEY = "cellsim-lab:settings:v1";

export const PRESETS = {
  // Reflections stay on even at Low: they cost one texture read, and without
  // them every metal surface renders black.
  low:    { label: "Low",    pixelRatio: 0.75, shadows: false, shadowSize: 0,    env: true,  lamps: 1, fps: 30, maxAniso: 1 },
  medium: { label: "Medium", pixelRatio: 1,    shadows: true,  shadowSize: 1024, env: true,  lamps: 2, fps: 60, maxAniso: 4 },
  high:   { label: "High",   pixelRatio: 1.5,  shadows: true,  shadowSize: 2048, env: true,  lamps: 4, fps: 60, maxAniso: 8 },
  ultra:  { label: "Ultra",  pixelRatio: 2,    shadows: true,  shadowSize: 4096, env: true,  lamps: 6, fps: 0,  maxAniso: 16 },
};

const DEFAULTS = { quality: "auto", fps: "preset", sensitivity: 1, fov: 72, invertY: false, bob: true, sound: true };

export const settings = (() => {
  try { return { ...DEFAULTS, ...(JSON.parse(localStorage.getItem(KEY)) || {}) }; } catch (e) { return { ...DEFAULTS }; }
})();

export function saveSettings() {
  try { localStorage.setItem(KEY, JSON.stringify(settings)); } catch (e) { /* private window: keep for this visit */ }
}

// Pick a preset from what the browser says about the graphics chip and memory.
export function detectPreset(gl) {
  let gpu = "";
  try {
    const ext = gl.getExtension("WEBGL_debug_renderer_info");
    gpu = ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
  } catch (e) { /* some browsers hide it */ }
  const mem = navigator.deviceMemory || 8;
  let pick = "medium";
  if (/NVIDIA|GeForce|RTX|Radeon RX|AMD Radeon Pro/i.test(gpu)) pick = "high";
  else if (/Apple M\d (Pro|Max|Ultra)/i.test(gpu)) pick = "high";
  else if (/Apple M\d|Apple GPU/i.test(gpu)) pick = "medium";
  else if (/Intel|Mali|Adreno|PowerVR|SwiftShader|llvmpipe/i.test(gpu)) pick = "low";
  if (mem <= 4) pick = "low";
  return { pick, gpu };
}

export function activePreset(detected) {
  const key = settings.quality === "auto" ? detected : settings.quality;
  const p = { ...PRESETS[key] || PRESETS.medium, key };
  if (settings.fps !== "preset") p.fps = +settings.fps;
  return p;
}
