// Sound effects, made in code (no files to download): footsteps on a tiled
// floor, picking things up and putting them down, glass, a drip, a pipette
// click. The audio context starts on your first click, as browsers require.
import { settings } from "./settings.js";

let ctx = null, master = null, noise = null;

export function unlockAudio() {
  if (ctx || !settings.sound) return;
  try {
    ctx = new (window.AudioContext || window.webkitAudioContext)();
    master = ctx.createGain();
    master.gain.value = 0.55;
    master.connect(ctx.destination);
    // One second of white noise, reused by every footstep.
    noise = ctx.createBuffer(1, ctx.sampleRate, ctx.sampleRate);
    const d = noise.getChannelData(0);
    for (let i = 0; i < d.length; i++) d[i] = Math.random() * 2 - 1;
  } catch (e) { ctx = null; }
}

function ready() {
  if (!settings.sound || !ctx) return false;
  if (ctx.state === "suspended") ctx.resume();
  return true;
}

function env(gain, t, attack, decay, peak) {
  gain.gain.setValueAtTime(0.0001, t);
  gain.gain.exponentialRampToValueAtTime(peak, t + attack);
  gain.gain.exponentialRampToValueAtTime(0.0001, t + attack + decay);
}

function burst({ f = 1200, q = 1.2, peak = 0.2, decay = 0.07, offset = Math.random() * 0.8 } = {}) {
  const t = ctx.currentTime;
  const src = ctx.createBufferSource(), bp = ctx.createBiquadFilter(), g = ctx.createGain();
  src.buffer = noise;
  bp.type = "bandpass"; bp.frequency.value = f; bp.Q.value = q;
  env(g, t, 0.004, decay, peak);
  src.connect(bp).connect(g).connect(master);
  src.start(t, offset, decay + 0.02);
}

function tone({ f = 800, to = f, type = "sine", peak = 0.12, decay = 0.08, attack = 0.003 } = {}) {
  const t = ctx.currentTime;
  const o = ctx.createOscillator(), g = ctx.createGain();
  o.type = type;
  o.frequency.setValueAtTime(f, t);
  if (to !== f) o.frequency.exponentialRampToValueAtTime(to, t + attack + decay);
  env(g, t, attack, decay, peak);
  o.connect(g).connect(master);
  o.start(t); o.stop(t + attack + decay + 0.03);
}

export const sfx = {
  // A shoe on hard floor tiles: a low thud and a short, bright scuff.
  step(loud = 0.6) {
    if (!ready()) return;
    const v = 0.06 + loud * 0.12;
    burst({ f: 180 + Math.random() * 60, q: 0.9, peak: v, decay: 0.05 });
    burst({ f: 2200 + Math.random() * 900, q: 1.6, peak: v * 0.35, decay: 0.035 });
  },
  take() { if (ready()) { burst({ f: 900, q: 2, peak: 0.08, decay: 0.04 }); tone({ f: 520, to: 700, peak: 0.04, decay: 0.05 }); } },
  place() { if (ready()) burst({ f: 420, q: 1.4, peak: 0.16, decay: 0.06 }); },
  glass() { if (ready()) { tone({ f: 3100 + Math.random() * 300, type: "sine", peak: 0.05, decay: 0.18 }); tone({ f: 4700, peak: 0.02, decay: 0.12 }); } },
  drip() { if (ready()) tone({ f: 1500, to: 520, peak: 0.08, decay: 0.07 }); },
  click() { if (ready()) { burst({ f: 3200, q: 3, peak: 0.07, decay: 0.02 }); tone({ f: 1900, peak: 0.03, decay: 0.03 }); } },
  pour() { if (ready()) burst({ f: 700, q: 0.6, peak: 0.07, decay: 0.45 }); },
  ok() { if (ready()) { tone({ f: 880, peak: 0.05, decay: 0.09 }); setTimeout(() => ready() && tone({ f: 1320, peak: 0.045, decay: 0.12 }), 90); } },
  warn() { if (ready()) tone({ f: 330, to: 260, type: "triangle", peak: 0.06, decay: 0.14 }); },
};
