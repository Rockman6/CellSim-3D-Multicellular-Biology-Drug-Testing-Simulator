// Node check of web/viewer/stream.js against real streams: plain, gzipped,
// and both schemas. Usage: node viewer_stream_check.js <file.jsonl[.gz]>...
"use strict";
const fs = require("fs");
const path = require("path");
const S = require(path.join(__dirname, "..", "..", "web", "viewer", "stream.js"));

(async () => {
  let failed = 0;
  const check = (cond, msg) => { if (!cond) { failed++; console.log("  FAIL " + msg); } };
  for (const file of process.argv.slice(2)) {
    const text = await S.bytesToText(fs.readFileSync(file));
    const model = S.normalize(S.parseJsonl(text));
    const s = model.series;
    check(s.t.length === model.snaps.length && s.t.length > 1, `${file}: one time point per snapshot`);
    check(s.t.every((t, i) => i === 0 || t > s.t[i - 1]), `${file}: time increases`);
    for (let i = 0; i < s.t.length; i++) {
      const sum = s.phase.reduce((a, p) => a + p[i], 0);
      check(Math.abs(sum - 1) < 1e-6 || sum === 0, `${file}: phase fractions sum to 1 at ${s.t[i]} h`);
    }
    if (model.kind === "dish") {
      for (const snap of model.snaps) {
        const n = snap.cells.x_um.length;
        for (const f of model.header.fields) check(snap.cells[f].length === n, `${file}: field ${f} has one value per cell`);
        check(snap.population.n_live === n, `${file}: n_live matches the cell list`);
      }
    }
    const csv = S.toCSV(model).trim().split("\n");
    check(csv.length === s.t.length + 1, `${file}: CSV has a header and one row per snapshot`);
    const width = csv[0].split(",").length;
    check(csv.every((r) => r.split(",").length === width), `${file}: CSV rows are rectangular`);
    console.log(`  ok   ${path.basename(file)}: ${model.kind}, ${s.t.length} snapshots, CSV ${width} columns`);
  }
  let threw = false;
  try { S.normalize(S.parseJsonl('{"schema":"something/else"}\n{}')); } catch (e) { threw = /unknown schema/.test(e.message); }
  check(threw, "an unknown schema is rejected with a clear message");
  console.log(failed ? `${failed} viewer stream checks FAILED` : "viewer stream checks pass");
  process.exit(failed ? 1 : 0);
})();
