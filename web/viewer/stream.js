/* cellsim stream parsing — no DOM, no three.js, so it runs under Node too.
 *
 * Two schemas are understood:
 *   cellsim.dish.stream/v1  cells in space (cellsim dish)
 *   cellsim.cell.stream/v1  a well-mixed population (cellsim cell-stream)
 * Both are JSON Lines: a header record, then one snapshot per record tick.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.CellSimStream = api;
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const SCHEMAS = {
    "cellsim.dish.stream/v1": "dish",
    "cellsim.cell.stream/v1": "population",
  };
  const PHASES = ["G1", "S", "G2", "M"];

  /** Bytes -> text, gunzipping when the data starts with the gzip magic number. */
  async function bytesToText(buffer) {
    const bytes = new Uint8Array(buffer);
    if (bytes.length > 2 && bytes[0] === 0x1f && bytes[1] === 0x8b) {
      if (typeof DecompressionStream === "undefined") {
        throw new Error("This browser cannot decompress .gz files; open the .jsonl instead.");
      }
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      return await new Response(stream).text();
    }
    return new TextDecoder().decode(bytes);
  }

  function parseJsonl(text) {
    const records = [];
    const lines = text.split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i].trim();
      if (!line) continue;
      try {
        records.push(JSON.parse(line));
      } catch (e) {
        throw new Error(`line ${i + 1} is not valid JSON`);
      }
    }
    return records;
  }

  /** Records -> {kind, header, snaps, series}. Throws on an unknown schema. */
  function normalize(records) {
    if (!records.length) throw new Error("the file is empty");
    const header = records[0];
    const kind = SCHEMAS[header.schema];
    if (!kind) {
      throw new Error(`unknown schema ${JSON.stringify(header.schema)}; expected one of ` +
        Object.keys(SCHEMAS).join(", "));
    }
    const snaps = records.slice(1);
    if (!snaps.length) throw new Error("the stream has a header but no snapshots");
    const t = snaps.map((s) => s.t_h);
    const phase = PHASES.map((p) => snaps.map((s) => s.population.phase_frac[p]));
    const series = { t, phase, drug: snaps.map((s) => firstDose(s.drug_uM)) };
    if (kind === "dish") {
      const pop = snaps.map((s) => s.population);
      series.live = pop.map((p) => p.n_live);
      series.apoptotic = pop.map((p) => p.n_apoptotic);
      series.necrotic = pop.map((p) => p.n_necrotic);
      series.radius = pop.map((p) => p.radius_um);
      series.necroticRadius = pop.map((p) => p.necrotic_radius_um);
      series.confluence = pop.map((p) => p.confluence);
      series.quiescent = pop.map((p) => p.quiescent_frac);
      series.hypoxic = pop.map((p) => p.hypoxic_frac);
    } else {
      const pop = snaps.map((s) => s.population);
      series.alive = pop.map((p) => p.alive_weight);
      series.dead = pop.map((p) => p.dead_weight);
      series.fold = pop.map((p) => p.fold_change);
      series.p53 = pop.map((p) => p.mean.p53);
      series.c3 = pop.map((p) => p.mean.C3);
      series.damage = pop.map((p) => p.mean.D);
    }
    return { kind, header, snaps, series };
  }

  function firstDose(d) {
    if (Array.isArray(d)) return d.length ? d[0] : 0;
    return typeof d === "number" ? d : 0;
  }

  /** Population time series as CSV text. */
  function toCSV(model) {
    const s = model.series;
    let cols;
    if (model.kind === "dish") {
      cols = [["t_h", s.t], ["drug_uM", s.drug], ["n_live", s.live], ["n_apoptotic", s.apoptotic],
        ["n_necrotic", s.necrotic], ["radius_um", s.radius], ["necrotic_radius_um", s.necroticRadius],
        ["confluence", s.confluence], ["quiescent_frac", s.quiescent], ["hypoxic_frac", s.hypoxic]];
    } else {
      cols = [["t_h", s.t], ["drug_uM", s.drug], ["alive_weight", s.alive], ["dead_weight", s.dead],
        ["fold_change", s.fold], ["mean_p53", s.p53], ["mean_caspase3", s.c3], ["mean_damage", s.damage]];
    }
    PHASES.forEach((p, i) => cols.push([`frac_${p}`, s.phase[i]]));
    const rows = [cols.map((c) => c[0]).join(",")];
    for (let i = 0; i < s.t.length; i++) {
      rows.push(cols.map((c) => fmt(c[1][i])).join(","));
    }
    return rows.join("\n") + "\n";
  }

  function fmt(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return "";
    return typeof v === "number" ? String(Math.round(v * 1e6) / 1e6) : String(v);
  }

  return { SCHEMAS, PHASES, bytesToText, parseJsonl, normalize, toCSV };
});
