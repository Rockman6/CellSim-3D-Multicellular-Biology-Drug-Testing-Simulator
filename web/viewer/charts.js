/* Minimal SVG line charts for the viewer: axes, a few series, a time cursor. */
(function (root) {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const PAD = { l: 38, r: 8, t: 8, b: 18 };

  function el(name, attrs, parent) {
    const e = document.createElementNS(NS, name);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function niceTicks(lo, hi, n) {
    if (!(hi > lo)) return [lo];
    const step0 = (hi - lo) / Math.max(1, n);
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= step0) || step0;
    const ticks = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * step; v += step) ticks.push(v);
    return ticks;
  }

  function label(v) {
    const a = Math.abs(v);
    if (a >= 1e4) return (v / 1e3).toFixed(0) + "k";
    if (a >= 1e3) return (v / 1e3).toFixed(1) + "k";
    if (a >= 10 || a === 0) return v.toFixed(0);
    if (a >= 1) return v.toFixed(1);
    return v.toPrecision(2);
  }

  /**
   * draw(svg, {t, series: [{values, color, label}], yMin, yMax, cursor, xUnit})
   * Colours are CSS colour strings (CSS variables allowed).
   */
  function draw(svg, opts) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    const box = svg.getBoundingClientRect();
    const W = Math.max(200, box.width || 340), H = Math.max(80, box.height || 130);
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const t = opts.t;
    if (!t || !t.length) return;
    const all = opts.series.flatMap((s) => s.values.filter((v) => Number.isFinite(v)));
    let yMin = opts.yMin !== undefined ? opts.yMin : Math.min(0, ...all);
    let yMax = opts.yMax !== undefined ? opts.yMax : Math.max(...all, 1e-9);
    if (yMax <= yMin) yMax = yMin + 1;
    const x0 = t[0], x1 = t[t.length - 1] > x0 ? t[t.length - 1] : x0 + 1;
    const X = (v) => PAD.l + (v - x0) / (x1 - x0) * (W - PAD.l - PAD.r);
    const Y = (v) => H - PAD.b - (v - yMin) / (yMax - yMin) * (H - PAD.t - PAD.b);

    for (const v of niceTicks(yMin, yMax, 3)) {
      el("line", { x1: PAD.l, x2: W - PAD.r, y1: Y(v), y2: Y(v), class: "grid" }, svg);
      const tx = el("text", { x: PAD.l - 4, y: Y(v) + 3, "text-anchor": "end" }, svg);
      tx.textContent = label(v);
    }
    const xt = niceTicks(x0, x1, 4);
    const unit = opts.xUnit === undefined ? "h" : opts.xUnit;
    xt.forEach((v, i) => {
      const tx = el("text", { x: X(v), y: H - 4, "text-anchor": "middle" }, svg);
      tx.textContent = label(v) + (i === xt.length - 1 && unit ? " " + unit : "");
    });
    let lx = PAD.l + 4;
    for (const s of opts.series) {
      const pts = [];
      s.values.forEach((v, i) => { if (Number.isFinite(v)) pts.push(`${X(t[i]).toFixed(1)},${Y(v).toFixed(1)}`); });
      if (pts.length) el("polyline", { points: pts.join(" "), class: "series", style: `stroke:${s.color}` }, svg);
      if (s.label) {
        el("rect", { x: lx, y: PAD.t, width: 8, height: 8, rx: 2, style: `fill:${s.color}` }, svg);
        const tx = el("text", { x: lx + 11, y: PAD.t + 7.5 }, svg);
        tx.textContent = s.label;
        lx += 16 + 6 * s.label.length;
      }
    }
    if (opts.cursor !== undefined) {
      el("line", { x1: X(opts.cursor), x2: X(opts.cursor), y1: PAD.t, y2: H - PAD.b, class: "cursor" }, svg);
    }
  }

  root.CellSimCharts = { draw };
})(typeof self !== "undefined" ? self : this);
