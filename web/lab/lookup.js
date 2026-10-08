/* Drug look-up for the Lab: a name or a structure in, everything the engine
   can honestly use out. All code ours; data from PubChem and ChEMBL (both
   public, both answer browser requests directly), chemistry from RDKit's own
   browser build, loaded only when someone looks something up.

   What comes back always says where each piece came from, because three
   very different things can be typed into the same box:
     - a known drug           -> PubChem record, ChEMBL's curated mechanism
     - a molecule nobody made -> RDKit computes its properties; no target
     - something impossible   -> RDKit says which atom breaks which rule */
"use strict";

const Lookup = (() => {
  const RDKIT = "https://cdn.jsdelivr.net/npm/@rdkit/rdkit@2026.9.1/dist/";
  const PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
  const CHEMBL = "https://www.ebi.ac.uk/chembl/api/data";
  let rdkit = null;

  async function getRDKit() {
    if (rdkit) return rdkit;
    await new Promise((ok, bad) => {
      const s = document.createElement("script");
      s.src = RDKIT + "RDKit_minimal.js"; s.onload = ok;
      s.onerror = () => bad(new Error("could not load RDKit"));
      document.head.appendChild(s);
    });
    rdkit = await window.initRDKitModule({ locateFile: () => RDKIT + "RDKit_minimal.wasm" });
    return rdkit;
  }

  // Up to three tries: ChEMBL's first request for an item it has not
  // cached can fail with a server error, and an error response carries no
  // CORS header, so the browser reports it as a blocked request. Seen on
  // the live site with both a filter query and a direct look-up; the
  // second or third try succeeds once ChEMBL has the item cached.
  const WAITS = [1200, 3000];
  async function json(url, attempt = 0) {
    let r = null;
    try { r = await fetch(url); } catch (e) { r = null; }
    if (r && r.status === 404) return null;
    if (r && r.ok) return r.json();
    if (attempt < WAITS.length) {
      await new Promise((ok) => setTimeout(ok, WAITS[attempt]));
      return json(url, attempt + 1);
    }
    throw new Error(r ? `${new URL(url).host} answered ${r.status}` : `${new URL(url).host} did not answer`);
  }

  // RDKit's own wording, made readable. The atom number is kept: it is
  // what someone fixing the structure needs.
  function explainError(raw) {
    const m = /Explicit valence for atom # (\d+) (\w+), (\d+), is greater than permitted/.exec(raw);
    if (m) return `Atom ${Number(m[1]) + 1} (${m[2]}) would have ${m[3]} bonds — more than ` +
      `${m[2]} can form. No such molecule can exist.`;
    if (/Can't kekulize/.test(raw)) return "The aromatic ring cannot be given a valid " +
      "arrangement of single and double bonds — check the lowercase (aromatic) atoms.";
    if (/unclosed ring/i.test(raw)) return "A ring is opened but never closed (a ring number " +
      "appears only once).";
    if (/Parse Error/.test(raw)) return "This is not a name PubChem knows, and not a valid " +
      "structure (SMILES) either.";
    return raw || "Not a valid structure.";
  }

  // Parse with RDKit; returns {mol} or {error}. The caller deletes mol.
  // One capture for the page's lifetime: opening a new one per parse
  // left later parses with an empty buffer and a generic message.
  let capture;
  async function parse(smiles) {
    const RD = await getRDKit();
    if (capture === undefined) capture = RD.set_log_capture ? RD.set_log_capture("rdApp.*") : null;
    const cap = capture;
    if (cap && cap.clear_buffer) cap.clear_buffer();
    const mol = RD.get_mol(smiles);
    if (mol && mol.is_valid()) return { mol };
    if (mol) mol.delete();
    const raw = cap ? (cap.get_buffer() || "").trim().split("\n").pop()
                      .replace(/^\[[\d:]+\]\s*/, "") : "";
    return { error: explainError(raw), raw };
  }

  function drawing(mol) {
    return mol.get_svg_with_highlights
      ? mol.get_svg_with_highlights(JSON.stringify({ width: 260, height: 180,
          backgroundColour: [1, 1, 1, 0] }))
      : mol.get_svg(260, 180);
  }

  function countNO(formula) {
    let n = 0;
    for (const [, el, k] of formula.matchAll(/([A-Z][a-z]?)(\d*)/g))
      if (el === "N" || el === "O") n += Number(k || 1);
    return n;
  }

  async function pubchemByName(name) {
    const props = "SMILES,IsomericSMILES,MolecularFormula,MolecularWeight,XLogP,TPSA,InChIKey,Title";
    const r = await json(`${PUBCHEM}/compound/name/${encodeURIComponent(name)}/property/${props}/JSON`);
    return r && r.PropertyTable.Properties[0];
  }

  async function pubchemByStructure(smiles) {
    // The structure goes in the query string, not the path: a SMILES can
    // hold "/" and "\". An unknown molecule comes back as CID 0, not a 404.
    const r = await fetch(`${PUBCHEM}/compound/smiles/cids/TXT?smiles=${encodeURIComponent(smiles)}`);
    if (!r.ok) return null;
    const cid = Number((await r.text()).trim().split("\n")[0]);
    if (!cid) return null;                  // PubChem has never seen it
    const props = "MolecularFormula,MolecularWeight,XLogP,TPSA,InChIKey,Title";
    const p = await json(`${PUBCHEM}/compound/cid/${cid}/property/${props}/JSON`);
    return p && p.PropertyTable.Properties[0];
  }

  async function chemblMechanisms(inchikey, name) {
    // The direct lookup by key is cached and fast; the filter query that
    // would avoid its 404 for unknown keys failed on live requests.
    let m = inchikey ? await json(`${CHEMBL}/molecule/${inchikey}.json`) : null;
    if (!m && name) {
      // metal complexes (cisplatin) carry different InChIKeys in the two databases
      const r = await json(`${CHEMBL}/molecule.json?pref_name__iexact=${encodeURIComponent(name)}&limit=1`);
      m = r && r.molecules[0];
    }
    if (!m) return null;
    const id = (m.molecule_hierarchy && m.molecule_hierarchy.parent_chembl_id) || m.molecule_chembl_id;
    // mechanisms are filed against the approved form (often a salt), so ask by parent
    const r = await json(`${CHEMBL}/mechanism.json?parent_molecule_chembl_id=${id}&limit=20`);
    const texts = [...new Set((r ? r.mechanisms : []).map(x => x.mechanism_of_action))];
    return { id, name: m.pref_name, texts };
  }

  function looksLikeStructure(q) {
    return /[()=#\[\]@+\\/]|\d/.test(q) || /^[BCNOPSFIcnops]+[a-z]*$/.test(q) && q.length <= 3;
  }

  /* The one entry point. Never throws for a bad query: the result says
     what went wrong in words. */
  async function lookup(query) {
    const q = query.trim();
    if (!q) return { kind: "empty" };
    let pc = null, smiles = null, from = null;
    if (!looksLikeStructure(q)) {
      pc = await pubchemByName(q);
      if (pc) { smiles = pc.SMILES || pc.IsomericSMILES; from = "PubChem"; }
    }
    if (!smiles) smiles = q;
    const parsed = await parse(smiles);
    if (parsed.error) return { kind: "invalid", query: q, error: parsed.error, raw: parsed.raw };
    const mol = parsed.mol;
    try {
      const d = JSON.parse(mol.get_descriptors());
      const svg = drawing(mol);
      if (!pc) { pc = await pubchemByStructure(smiles); if (pc) from = "PubChem"; }
      const props = {
        mw: pc ? Number(pc.MolecularWeight) : d.amw,
        logp: pc && pc.XLogP != null ? Number(pc.XLogP) : d.CrippenClogP,
        tpsa: pc && pc.TPSA != null ? Number(pc.TPSA) : d.tpsa,
        n_plus_o: pc ? countNO(pc.MolecularFormula) : d.lipinskiHBA,
        formula: pc ? pc.MolecularFormula : null,
      };
      let chembl = null, chemblError = null;
      if (pc) {
        try { chembl = await chemblMechanisms(pc.InChIKey, pc.Title || q); }
        catch (e) { chemblError = String(e.message || e); }
      }
      return { kind: pc ? "known" : "new", query: q,
               name: pc ? (pc.Title || q) : "your molecule", smiles, svg, props,
               propsFrom: pc ? "PubChem" : "computed from your structure (RDKit)",
               chembl, chemblError, from };
    } finally {
      mol.delete();
    }
  }

  return { lookup, parse, getRDKit };
})();
