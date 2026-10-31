// Living cells under the microscope, drawn from a Cell Sim dish snapshot
// (schema cellsim.dish.stream/v1). Positions, phases and states are the
// engine's; only the cell outlines are drawn here, in a phase-contrast style.
const PHASES = ["G1", "S", "G2", "M"];

function hash(i, j) {
  let h = Math.imul(i | 0, 374761393) ^ Math.imul(j | 0, 668265263);
  h = Math.imul(h ^ (h >>> 13), 1274126177);
  return ((h ^ (h >>> 16)) >>> 0) / 4294967296;
}

export function makeCultureField(header, snap, { tint = "#e9dfe0" } = {}) {
  const sp = header.spacing_um || 15;
  const side = (header.grid_shape ? header.grid_shape[0] : 60) * sp;
  const off = side / 2;
  const C = snap.cells;
  const n = C.x_um.length;
  const cells = [];
  for (let i = 0; i < n; i++) {
    const x = C.x_um[i] - off, y = C.y_um[i] - off;
    const gx = Math.round(C.x_um[i] / sp), gy = Math.round(C.y_um[i] / sp);
    cells.push({ i, x, y, phase: C.phase[i], p53: C.p53[i], c3: C.caspase3[i], dmg: C.damage[i], o2: C.o2_mmHg[i], gs: C.growth_signal[i], lin: C.lineage[i],
      a: hash(gx, gy) * Math.PI, e: 0.65 + hash(gy, gx) * 0.3, k: hash(gx * 7, gy * 3) });
  }
  const bodies = (snap.apoptotic_um || []).map(([x, y]) => ({ x: x - off, y: y - off, k: hash(x | 0, y | 0) }));
  const debris = (snap.necrotic_um || []).map(([x, y]) => ({ x: x - off, y: y - off, k: hash(y | 0, x | 0) }));

  function draw(ctx, W, v) {
    const s = W / v.fov;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    // Phase contrast: grey field, tinted by the medium.
    const g = ctx.createRadialGradient(W / 2, W / 2, W * 0.1, W / 2, W / 2, W * 0.7);
    g.addColorStop(0, tint); g.addColorStop(1, "#b9b3ac");
    ctx.fillStyle = g; ctx.fillRect(0, 0, W, W);
    ctx.setTransform(s, 0, 0, s, W / 2 - v.cx * s, W / 2 - v.cy * s);
    const half = v.fov / 2 + 20;
    const vis = (o) => Math.abs(o.x - v.cx) < half && Math.abs(o.y - v.cy) < half;
    const lw = Math.max(0.8, 0.9 / s);
    // Spread interphase cells: soft grey body, bright halo, darker nucleus.
    for (const c of cells) {
      if (!vis(c) || c.phase === 3) continue;
      const dying = c.c3 > 0.35;
      const r = (dying ? 7 : 9.5) * (0.9 + c.k * 0.25);
      ctx.save();
      ctx.translate(c.x, c.y); ctx.rotate(c.a);
      ctx.beginPath(); ctx.ellipse(0, 0, r, r * c.e, 0, 0, Math.PI * 2);
      ctx.fillStyle = dying ? "rgba(70,66,62,0.55)" : "rgba(112,108,104,0.42)";
      ctx.fill();
      ctx.lineWidth = lw * 1.4; ctx.strokeStyle = "rgba(255,255,250,0.75)"; ctx.stroke();
      if (s > 0.12) {
        ctx.beginPath(); ctx.ellipse(r * 0.08, 0, r * 0.42, r * 0.36 * c.e + 0.6, 0.3, 0, Math.PI * 2);
        ctx.fillStyle = dying ? "rgba(40,36,34,0.8)" : "rgba(78,74,70,0.6)"; ctx.fill();
        if (s > 0.5 && !dying) { ctx.fillStyle = "rgba(40,38,36,0.8)"; ctx.beginPath(); ctx.arc(r * 0.15, -0.5, 1.1, 0, Math.PI * 2); ctx.fill(); }
      }
      ctx.restore();
    }
    // Mitotic cells round up and look bright, with their chromosomes in a line.
    for (const c of cells) {
      if (!vis(c) || c.phase !== 3) continue;
      ctx.beginPath(); ctx.arc(c.x, c.y, 6.2, 0, Math.PI * 2);
      ctx.fillStyle = "rgba(240,238,230,0.85)"; ctx.fill();
      ctx.lineWidth = lw * 2.2; ctx.strokeStyle = "rgba(60,58,55,0.75)"; ctx.stroke();
      if (s > 0.3) { ctx.save(); ctx.translate(c.x, c.y); ctx.rotate(c.a); ctx.fillStyle = "rgba(50,48,46,0.75)"; ctx.fillRect(-3.6, -0.9, 7.2, 1.8); ctx.restore(); }
    }
    // Apoptotic bodies: little bright blebs.
    for (const b of bodies) {
      if (!vis(b)) continue;
      for (let j = 0; j < 4; j++) {
        const ang = b.k * 6.28 + j * 1.7, d = 2.2 + j * 0.6;
        ctx.beginPath(); ctx.arc(b.x + Math.cos(ang) * d, b.y + Math.sin(ang) * d, 1.6 + ((j * 0.37 + b.k) % 1) * 1.4, 0, Math.PI * 2);
        ctx.fillStyle = "rgba(245,243,236,0.9)"; ctx.fill();
        ctx.lineWidth = lw; ctx.strokeStyle = "rgba(55,52,50,0.7)"; ctx.stroke();
      }
    }
    // Necrotic debris: dark granular remains.
    ctx.fillStyle = "rgba(38,34,32,0.7)";
    for (const d of debris) {
      if (!vis(d)) continue;
      for (let j = 0; j < 7; j++) ctx.fillRect(d.x + Math.cos(d.k * 9 + j) * 4, d.y + Math.sin(d.k * 7 + j * 2) * 4, 2, 2);
    }
  }

  return {
    culture: true,
    bubbles: [],
    inSample: (x, y) => Math.abs(x) < off + 30 && Math.abs(y) < off + 30,
    animated: () => false,
    draw,
    // The engine's numbers for the cell under the cursor.
    pickCell(x, y) {
      let best = null, d2 = 400;   // within 20 µm of the click
      for (const c of cells) { const d = (c.x - x) ** 2 + (c.y - y) ** 2; if (d < d2) { d2 = d; best = c; } }
      return best;
    },
    counts() {
      const ph = [0, 0, 0, 0];
      for (const c of cells) ph[c.phase]++;
      return { live: n, apoptotic: bodies.length, necrotic: debris.length, phases: ph };
    },
    phaseName: (p) => PHASES[p],
  };
}
