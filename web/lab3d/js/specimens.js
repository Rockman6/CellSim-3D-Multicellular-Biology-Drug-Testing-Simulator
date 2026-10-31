/* Specimens and what they look like under the microscope (ported from the voxel lab).
   Each field is laid out in micrometres (µm) at real cell sizes, from a fixed
   seed, so the same slide always looks the same. */
// Seeded random numbers (mulberry32), so every slide looks the same each visit.
function rng(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

  export const STAINS = {
    water: { name: "Water", swatch: "#9fd3e6" },
    iodine: { name: "Iodine", swatch: "#b5651d" },
    mb: { name: "Methylene blue", swatch: "#2f62c9" },
  };

  export const SPECIMENS = {
    onion: {
      name: "Onion skin",
      kind: "plant",
      best: "iodine",
      seed: 11,
      placeText: "Peel a thin, see-through layer from the inside of an onion scale and lay it flat on the slide.",
      stainTip: "Tip: iodine stains onion cells yellow-brown, so the walls and nuclei show up.",
      sizeText: "about 200–400 µm long",
      sizeRange: [200, 400],
      measureWhat: "long",
    },
    cheek: {
      name: "Cheek cells",
      kind: "animal",
      best: "mb",
      seed: 23,
      placeText: "Gently rub a clean swab on the inside of your cheek, then wipe it across the middle of the slide.",
      stainTip: "Tip: methylene blue stains the nucleus dark blue. Without a stain, cheek cells are almost invisible.",
      sizeText: "about 50–70 µm across",
      sizeRange: [45, 75],
      measureWhat: "across",
    },
    elodea: {
      name: "Pondweed leaf",
      kind: "plant",
      best: "water",
      seed: 37,
      placeText: "Pull off one small, young leaf near the tip of the pondweed and lay it flat on the slide.",
      stainTip: "Tip: water is enough. The green chloroplasts are already easy to see.",
      sizeText: "about 50–100 µm long",
      sizeRange: [50, 110],
      measureWhat: "long",
    },
  };

  // Field of view of a standard school microscope with a 10× eyepiece.
  export const FOV_UM = { 4: 4500, 10: 1800, 40: 450 };

  const BG = "#f4f1e6";

  // ---------- geometry helpers ----------
  function noisyOutline(r, rx, ry, n, rough) {
    const ph = [r() * 6.28, r() * 6.28, r() * 6.28];
    const pts = [];
    for (let i = 0; i < n; i++) {
      const a = (i / n) * Math.PI * 2;
      const k = 1 + rough * (0.5 * Math.sin(a * 3 + ph[0]) + 0.3 * Math.sin(a * 5 + ph[1]) + 0.2 * Math.sin(a * 9 + ph[2]));
      pts.push([Math.cos(a) * rx * k, Math.sin(a) * ry * k]);
    }
    return pts;
  }

  function inPoly(x, y, pts) {
    let inside = false;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const [xi, yi] = pts[i], [xj, yj] = pts[j];
      if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  }

  function tracePoly(ctx, pts) {
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    ctx.closePath();
  }

  /* Rows of brick-like cells with shared walls (onion skin, pondweed leaf).
     Row edges wave gently, and end walls lean, so the walls never leave gaps. */
  function brickCells(r, o) {
    const rows = [];
    for (let y = o.y0; y < o.y1 + 200; y += o.h0 + r() * o.hVar) rows.push({ y, amp: 2 + r() * 4, ph: r() * 6.28 });
    const edge = (i, x) => rows[i].y + rows[i].amp * Math.sin(x * 0.0035 + rows[i].ph);
    const cells = [];
    for (let i = 0; i < rows.length - 1; i++) {
      const midY = (rows[i].y + rows[i + 1].y) / 2;
      const lenScale = o.lenAt ? o.lenAt(midY) : 1;
      let x = o.x0 - r() * o.L0;
      let s = (r() - 0.5) * o.slant;
      while (x < o.x1) {
        const L = (o.L0 + r() * o.LVar) * lenScale;
        const s2 = (r() - 0.5) * o.slant;
        const x2 = x + L;
        const cx = x + L / 2;
        if (o.keep(cx, midY)) {
          const poly = [
            [x + s, edge(i, x + s)], [x2 + s2, edge(i, x2 + s2)],
            [x2 - s2, edge(i + 1, x2 - s2)], [x - s, edge(i + 1, x - s)],
          ];
          cells.push({
            poly,
            x0: Math.min(x + s, x - s), x1: Math.max(x2 + s2, x2 - s2),
            y0: rows[i].y - 8, y1: rows[i + 1].y + 8,
            top: rows[i].y, bot: rows[i + 1].y, left: x, right: x2,
            tone: r(),
          });
        }
        x = x2; s = s2;
      }
    }
    return cells;
  }

  function visible(c, v) {
    return c.x1 > v.l && c.x0 < v.r && c.y1 > v.t && c.y0 < v.b;
  }

  // Air bubbles, if the cover slip went down flat.
  function makeBubbles(r, cover, outline) {
    const out = [];
    const n = cover === "flat" ? 9 : r() < 0.5 ? 1 : 0;
    let tries = 0;
    while (out.length < n && tries++ < 200) {
      const x = (r() - 0.5) * 5200, y = (r() - 0.5) * 3600;
      if (!inPoly(x, y, outline)) continue;
      const rad = cover === "flat" ? 90 + r() * r() * 650 : 35 + r() * 25;
      out.push({ x, y, rad });
    }
    return out;
  }

  function drawBubbles(ctx, bubbles, v, s) {
    for (const b of bubbles) {
      if (b.x + b.rad < v.l || b.x - b.rad > v.r || b.y + b.rad < v.t || b.y - b.rad > v.b) continue;
      ctx.beginPath();
      ctx.arc(b.x, b.y, b.rad, 0, Math.PI * 2);
      ctx.fillStyle = "rgba(252,251,246,0.92)";
      ctx.fill();
      ctx.lineWidth = Math.max(b.rad * 0.11, 2 / s);
      ctx.strokeStyle = "rgba(38,38,36,0.78)";
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(b.x - b.rad * 0.25, b.y - b.rad * 0.25, b.rad * 0.35, Math.PI, Math.PI * 1.6);
      ctx.lineWidth = Math.max(b.rad * 0.05, 1 / s);
      ctx.strokeStyle = "rgba(255,255,255,0.9)";
      ctx.stroke();
    }
  }

  function makeDust(r) {
    const out = [];
    for (let i = 0; i < 60; i++) out.push({ x: (r() - 0.5) * 7000, y: (r() - 0.5) * 5000, rad: 2 + r() * r() * 12 });
    return out;
  }
  function drawDust(ctx, dust, v, s) {
    ctx.fillStyle = "rgba(60,55,45,0.55)";
    for (const d of dust) {
      if (d.x < v.l || d.x > v.r || d.y < v.t || d.y > v.b) continue;
      const w = Math.max(d.rad, 1.2 / s);
      ctx.fillRect(d.x - w / 2, d.y - w / 2, w, w * 0.7);
    }
  }

  function viewRect(v) {
    const h = v.fov / 2;
    return { l: v.cx - h, r: v.cx + h, t: v.cy - h, b: v.cy + h };
  }

  function begin(ctx, W, v) {
    const s = W / v.fov;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, W, W);
    ctx.setTransform(s, 0, 0, s, W / 2 - v.cx * s, W / 2 - v.cy * s);
    return s;
  }

  // ---------- onion skin ----------
  const ONION_STYLE = {
    water: { cyto: "rgba(232,226,204,0.16)", wall: "rgba(120,110,88,0.30)", nuc: "rgba(140,130,110,0.20)", nucEdge: "rgba(105,96,78,0.30)" },
    iodine: { cyto: "rgba(226,188,108,0.34)", wall: "rgba(146,92,28,0.78)", nuc: "rgba(132,74,18,0.78)", nucEdge: "rgba(96,52,10,0.85)" },
    mb: { cyto: "rgba(150,180,226,0.20)", wall: "rgba(66,98,160,0.48)", nuc: "rgba(28,58,150,0.82)", nucEdge: "rgba(18,40,110,0.9)" },
  };

  function onionField(prep) {
    const r = rng(SPECIMENS.onion.seed);
    const outline = noisyOutline(r, 3300, 2300, 80, 0.13);
    const cells = brickCells(r, {
      x0: -3700, x1: 3700, y0: -2700, y1: 2700, h0: 54, hVar: 18, L0: 200, LVar: 160, slant: 26,
      keep: (x, y) => inPoly(x, y, outline),
    });
    for (const c of cells) {
      c.nx = c.left + (c.right - c.left) * (0.2 + r() * 0.6);
      c.ny = c.top + (c.bot - c.top) * (0.3 + r() * 0.4);
      c.nr = 8 + r() * 3;
      c.nucleoli = 1 + Math.floor(r() * 2);
    }
    const bubbles = makeBubbles(r, prep.cover, outline);
    const dust = makeDust(r);
    const st = ONION_STYLE[prep.stain] || ONION_STYLE.water;

    return {
      bubbles,
      inSample: (x, y) => inPoly(x, y, outline),
      animated: () => false,
      draw(ctx, W, v) {
        const s = begin(ctx, W, v);
        const vr = viewRect(v);
        const vis = cells.filter((c) => visible(c, vr));
        ctx.fillStyle = st.cyto;
        ctx.beginPath();
        for (const c of vis) tracePoly(ctx, c.poly);
        ctx.fill();
        ctx.lineWidth = Math.max(2.6, 0.9 / s);
        ctx.strokeStyle = st.wall;
        ctx.lineJoin = "round";
        ctx.beginPath();
        for (const c of vis) tracePoly(ctx, c.poly);
        ctx.stroke();
        // nuclei
        ctx.fillStyle = st.nuc;
        ctx.strokeStyle = st.nucEdge;
        ctx.lineWidth = Math.max(1.4, 0.7 / s);
        for (const c of vis) {
          const rx = Math.max(c.nr * 1.25, 1.2 / s), ry = Math.max(c.nr, 1.2 / s);
          ctx.beginPath();
          ctx.ellipse(c.nx, c.ny, rx, ry, 0, 0, Math.PI * 2);
          ctx.fill();
          if (s > 0.3) ctx.stroke();
          if (s > 0.6) {
            ctx.fillStyle = st.nucEdge;
            for (let k = 0; k < c.nucleoli; k++) {
              ctx.beginPath();
              ctx.arc(c.nx + (k ? 3 : -2), c.ny + (k ? 1.5 : -1), 1.8, 0, Math.PI * 2);
              ctx.fill();
            }
            ctx.fillStyle = st.nuc;
          }
        }
        drawDust(ctx, dust, vr, s);
        drawBubbles(ctx, bubbles, vr, s);
      },
    };
  }

  // ---------- pondweed (Elodea) leaf ----------
  const ELODEA_STYLE = {
    water: { cyto: "rgba(196,224,166,0.38)", wall: "rgba(88,128,66,0.62)", chl: "#3d8c37", chlEdge: "#2a6a27", leaf: "rgba(110,170,80,0.55)" },
    iodine: { cyto: "rgba(214,196,120,0.40)", wall: "rgba(130,96,40,0.70)", chl: "#2c3a2b", chlEdge: "#1a2419", leaf: "rgba(120,120,60,0.6)" },
    mb: { cyto: "rgba(150,200,190,0.36)", wall: "rgba(60,100,140,0.60)", chl: "#2f6f5a", chlEdge: "#1f4f40", leaf: "rgba(80,140,130,0.55)" },
  };

  function elodeaField(prep) {
    const r = rng(SPECIMENS.elodea.seed);
    const LEN = 2900, WID = 820;
    const halfWidth = (x) => {
      const t = x / LEN;
      if (Math.abs(t) >= 1) return 0;
      const base = WID * Math.pow(1 - t * t, 0.55) * (t > 0 ? 1 - 0.35 * t : 1);
      const tooth = 14 * ((((x / 95) % 1) + 1) % 1); // small saw teeth along the edge
      return base + tooth;
    };
    const outline = [];
    for (let i = 0; i <= 120; i++) { const x = -LEN + (2 * LEN * i) / 120; outline.push([x, -halfWidth(x)]); }
    for (let i = 120; i >= 0; i--) { const x = -LEN + (2 * LEN * i) / 120; outline.push([x, halfWidth(x)]); }
    const inLeaf = (x, y) => Math.abs(y) < halfWidth(x) - 4;

    const cells = brickCells(r, {
      x0: -LEN, x1: LEN, y0: -WID - 60, y1: WID + 60, h0: 21, hVar: 7, L0: 60, LVar: 45, slant: 10,
      lenAt: (y) => (Math.abs(y) < 45 ? 2.2 : 1), // long narrow cells along the midrib
      keep: inLeaf,
    });
    for (const c of cells) {
      const w = c.right - c.left, h = c.bot - c.top;
      const n = Math.max(6, Math.round((w * h) / 70));
      c.chl = [];
      for (let k = 0; k < n; k++) {
        const edge = r() < 0.62;
        c.chl.push({ edge, p: r(), d: 3 + r() * 3, fx: r(), fy: r(), a: r() * 3.14 });
      }
      c.dir = r() < 0.5 ? -1 : 1;
      c.speed = 0.012 + r() * 0.025; // fraction of the cell's edge per second
    }
    const bubbles = makeBubbles(r, prep.cover, outline);
    const dust = makeDust(r);
    const st = ELODEA_STYLE[prep.stain] || ELODEA_STYLE.water;

    // Position on the inside edge of a cell, p in [0,1), d µm in from the wall.
    function along(c, p, d) {
      const x0 = c.left + d, x1 = c.right - d, y0 = c.top + d, y1 = c.bot - d;
      const w = x1 - x0, h = y1 - y0, P = 2 * (w + h);
      let q = (((p % 1) + 1) % 1) * P;
      if (q < w) return [x0 + q, y0];
      q -= w; if (q < h) return [x1, y0 + q];
      q -= h; if (q < w) return [x1 - q, y1];
      q -= w; return [x0, y1 - q];
    }

    return {
      bubbles,
      inSample: inLeaf,
      animated: (v) => v.fov <= 450,
      draw(ctx, W, v) {
        const s = begin(ctx, W, v);
        const vr = viewRect(v);
        if (s < 0.15) {
          // Low power: the leaf reads as a green shape with faint rows.
          ctx.fillStyle = st.leaf;
          ctx.beginPath(); tracePoly(ctx, outline); ctx.fill();
          ctx.fillStyle = "rgba(60,110,50,0.35)";
          ctx.fillRect(-LEN, -40, 2 * LEN * 0.95, 80);
        }
        const vis = cells.filter((c) => visible(c, vr));
        ctx.fillStyle = st.cyto;
        ctx.beginPath();
        for (const c of vis) tracePoly(ctx, c.poly);
        ctx.fill();
        if (s >= 0.15) {
          ctx.lineWidth = Math.max(1.8, 0.8 / s);
          ctx.strokeStyle = st.wall;
          ctx.beginPath();
          for (const c of vis) tracePoly(ctx, c.poly);
          ctx.stroke();
        }
        if (s >= 0.15) {
          const t = v.t || 0;
          const big = s >= 0.6;
          ctx.fillStyle = st.chl;
          ctx.strokeStyle = st.chlEdge;
          ctx.lineWidth = 0.6;
          const dot = Math.max(4, 1.1 / s);
          for (const c of vis) {
            for (const k of c.chl) {
              let x, y;
              if (k.edge) [x, y] = along(c, k.p + c.dir * c.speed * t, k.d);
              else { x = c.left + 5 + k.fx * (c.right - c.left - 10); y = c.top + 4 + k.fy * (c.bot - c.top - 8); }
              if (big) {
                ctx.beginPath();
                ctx.ellipse(x, y, 2.6, 2.0, k.a, 0, Math.PI * 2);
                ctx.fill();
                ctx.stroke();
              } else ctx.fillRect(x - dot / 2, y - dot / 2, dot, dot);
            }
          }
        }
        drawDust(ctx, dust, vr, s);
        drawBubbles(ctx, bubbles, vr, s);
      },
    };
  }

  // ---------- cheek cells ----------
  const CHEEK_STYLE = {
    water: { fill: "rgba(222,216,204,0.10)", edge: "rgba(120,114,102,0.22)", nuc: "rgba(140,134,122,0.16)", bact: "rgba(120,114,102,0.18)" },
    iodine: { fill: "rgba(226,196,122,0.22)", edge: "rgba(160,110,40,0.50)", nuc: "rgba(150,90,30,0.55)", bact: "rgba(150,100,40,0.4)" },
    mb: { fill: "rgba(122,160,222,0.32)", edge: "rgba(48,88,170,0.66)", nuc: "rgba(24,48,138,0.88)", bact: "rgba(30,60,150,0.85)" },
  };

  function cheekField(prep) {
    const r = rng(SPECIMENS.cheek.seed);
    const outline = noisyOutline(r, 2800, 2000, 72, 0.22);
    const cells = [];
    function addCell(x, y) {
      if (!inPoly(x, y, outline)) return;
      const R = 24 + r() * 10, rot = r() * 6.28, n = 9;
      const poly = [];
      for (let i = 0; i < n; i++) {
        const a = rot + (i / n) * Math.PI * 2;
        const rr = R * (0.82 + r() * 0.3);
        poly.push([x + Math.cos(a) * rr, y + Math.sin(a) * rr * 0.92]);
      }
      cells.push({
        poly, x0: x - R * 1.2, x1: x + R * 1.2, y0: y - R * 1.2, y1: y + R * 1.2,
        nx: x + (r() - 0.5) * 9, ny: y + (r() - 0.5) * 9, nr: 4 + r() * 1.2,
        fold: r() < 0.3 ? { a: r() * 6.28, R } : null, cx: x, cy: y,
      });
    }
    for (let i = 0; i < 1100; i++) addCell((r() - 0.5) * 6000, (r() - 0.5) * 4400);
    for (let k = 0; k < 40; k++) {
      const cx = (r() - 0.5) * 4800, cy = (r() - 0.5) * 3400, m = 6 + Math.floor(r() * 8);
      for (let i = 0; i < m; i++) addCell(cx + (r() - 0.5) * 140, cy + (r() - 0.5) * 140);
    }
    const bact = [];
    for (let k = 0; k < 45; k++) {
      const cx = (r() - 0.5) * 4800, cy = (r() - 0.5) * 3400;
      for (let i = 0; i < 16; i++) bact.push([cx + (r() - 0.5) * 50, cy + (r() - 0.5) * 50]);
    }
    const bubbles = makeBubbles(r, prep.cover, outline);
    const dust = makeDust(r);
    const st = CHEEK_STYLE[prep.stain] || CHEEK_STYLE.water;

    return {
      bubbles,
      inSample: (x, y) => inPoly(x, y, outline),
      animated: () => false,
      draw(ctx, W, v) {
        const s = begin(ctx, W, v);
        const vr = viewRect(v);
        const vis = cells.filter((c) => visible(c, vr));
        ctx.lineJoin = "round";
        ctx.lineWidth = Math.max(1.3, 0.8 / s);
        for (const c of vis) {
          ctx.beginPath(); tracePoly(ctx, c.poly);
          ctx.fillStyle = st.fill; ctx.fill();
          ctx.strokeStyle = st.edge; ctx.stroke();
          if (c.fold && s > 0.3) {
            ctx.beginPath();
            ctx.arc(c.cx, c.cy, c.fold.R * 0.7, c.fold.a, c.fold.a + 1.4);
            ctx.stroke();
          }
        }
        ctx.fillStyle = st.nuc;
        for (const c of vis) {
          ctx.beginPath();
          ctx.ellipse(c.nx, c.ny, Math.max(c.nr * 1.1, 1.2 / s), Math.max(c.nr, 1.2 / s), 0, 0, Math.PI * 2);
          ctx.fill();
        }
        if (s > 0.6) {
          ctx.fillStyle = st.bact;
          for (const [x, y] of bact) if (x > vr.l && x < vr.r && y > vr.t && y < vr.b) ctx.fillRect(x - 0.9, y - 0.9, 1.8, 1.8);
        }
        drawDust(ctx, dust, vr, s);
        drawBubbles(ctx, bubbles, vr, s);
      },
    };
  }

  // A slide with nothing on it: dust, and bubbles if the cover slip went down flat.
  function emptyField(prep) {
    const r = rng(5);
    const outline = noisyOutline(r, 3000, 2000, 48, 0.1);
    const bubbles = makeBubbles(r, prep.cover, outline);
    const dust = makeDust(r);
    return {
      bubbles,
      inSample: () => false,
      animated: () => false,
      draw(ctx, W, v) {
        const s = begin(ctx, W, v);
        const vr = viewRect(v);
        drawDust(ctx, dust, vr, s);
        drawBubbles(ctx, bubbles, vr, s);
      },
    };
  }

  export function makeField(prep) {
    if (prep.specimen === "onion") return onionField(prep);
    if (prep.specimen === "cheek") return cheekField(prep);
    if (prep.specimen === "elodea") return elodeaField(prep);
    return emptyField(prep);
  }
