// jsdom test of dist/inafinder.html: load files, scan (single thread), check results, downloads, saved-report round trip.
// node tests/ui_test.js <html> <outdir> <genome.fna> [annotation] [more genomes / .hmm files...]
// env: MODEL=bacterial|fungal|custom, EXPLORE=1, ROLES=name:N,name:C (custom roles)
// jsdom does not implement navigation; the "Not implemented: navigation" messages caused by download clicks are filtered out.
const fs = require("fs"), path = require("path"), zlib = require("zlib");
const JSDOM_DIR = path.join(__dirname, "../../genome_map_viewer/tests/jsdom_env/node_modules/jsdom");
const { JSDOM, VirtualConsole } = require(JSDOM_DIR);
const [, , htmlPath, outDir, ...inputs] = process.argv; fs.mkdirSync(outDir, { recursive: true });
function open(html, url) {
  const errors = []; const vc = new VirtualConsole(); const keep = m => { if (!/Not implemented: navigation/i.test(m)) errors.push(m); };
  vc.on("jsdomError", e => keep(String(e.message || e))); vc.on("error", e => keep(String(e)));
  const downloads = [];
  const dom = new JSDOM(html, { runScripts: "dangerously", pretendToBeVisual: true, url, virtualConsole: vc, beforeParse(w) {
    const noop = () => {}; const ctx = new Proxy({}, { get: (t, k) => (k in t ? t[k] : noop), set: (t, k, v) => { t[k] = v; return true; } });
    w.HTMLCanvasElement.prototype.getContext = () => ctx; w.TextEncoder = TextEncoder; w.TextDecoder = TextDecoder;
    w.URL.createObjectURL = b => { downloads.push(b); return "blob:" + downloads.length; };
  } });
  return { dom, w: dom.window, errors, downloads };
}
const readBlob = (w, b) => new Promise(res => { const r = new w.FileReader(); r.onload = () => res(r.result); r.readAsText(b); });
(async () => {
  const html = fs.readFileSync(htmlPath, "utf8");
  const A = open(html, "http://localhost/inafinder.html#selftest");
  await new Promise(r => A.w.addEventListener("load", r));
  const out = { errors: A.errors, selftest: JSON.parse(A.w.document.getElementById("selftest").textContent || "null") };
  const files = inputs.map(p => { let b = fs.readFileSync(p); let name = path.basename(p); if (name.endsWith(".gz")) { b = zlib.gunzipSync(b); name = name.slice(0, -3); } return new A.w.File([b.toString("utf8")], name); });
  const app = A.w.InaFinderApp;
  await app.addFiles(files);
  out.files = app.st.files.map(f => ({ name: f.name, kind: f.kind, annot: f.annot ? app.st.files.find(x => x.id === f.annot).name : null }));
  if (process.env.MODEL) { const r = A.w.document.querySelector(`input[name=model][value=${process.env.MODEL}]`); r.checked = true; r.dispatchEvent(new A.w.Event("change")); }
  if (process.env.ROLES) { for (const kv of process.env.ROLES.split(",")) { const [nm, role] = kv.split(":"); const i = app.st.params.roles.findIndex(r => r.profile === nm); const sel = A.w.document.querySelector(`select[data-role="${i}"]`); sel.value = role; sel.dispatchEvent(new A.w.Event("change")); } }
  if (process.env.EXPLORE) A.w.document.getElementById("optExplore").checked = true;
  out.model = { model: app.st.model, roles: app.st.params.roles.map(r => [r.profile, r.role]), explore: A.w.document.getElementById("optExplore").checked };
  out.runDisabledBefore = A.w.document.getElementById("runBtn").disabled;
  const t0 = Date.now(); await app.run({ mainThread: true }); out.runSeconds = (Date.now() - t0) / 1000;
  out.runs = app.st.runs.map(r => r.error ? { error: r.error } : { genome: r.result.genome, nORF: r.result.nORF, loci: r.result.loci.map(l => [l.contig, l.strand, l.cds_start, l.cds_end, l.structural_class, l.protein_length_aa, l.annotation_state, l.annot_locus_tag, l.contig_edge, l.period4_ratio]) });
  out.locRows = A.w.document.querySelectorAll("#locTab tbody tr").length;
  if (out.locRows) { app.selectLocus(0, 0); out.detailHasSvg = !!A.w.document.querySelector("#detail svg"); fs.writeFileSync(path.join(outDir, "detail.svg"), A.w.document.querySelector("#detail svg").outerHTML); }
  // downloads
  for (const k of ["tsv", "faa", "fna", "gff", "regions", "gbk", "json"]) { const n0 = A.downloads.length; app.DL[k](); const b = A.downloads[n0]; const t = k === "gbk" ? null : await readBlob(A.w, b); out["dl_" + k] = b ? b.size : 0; if (t !== null) fs.writeFileSync(path.join(outDir, "out." + k), t);
    if (k === "gbk") { const ab = await new Promise(res => { const r = new A.w.FileReader(); r.onload = () => res(r.result); r.readAsArrayBuffer(b); }); fs.writeFileSync(path.join(outDir, "out_gbk.zip"), Buffer.from(ab)); } }
  const rep = app.saveReport(); fs.writeFileSync(path.join(outDir, "report.html"), rep); out.reportBytes = rep.length;
  // reopen the report
  const B = open(rep, "http://localhost/report.html#selftest"); await new Promise(r => B.w.addEventListener("load", r));
  out.report = { errors: B.errors, selftest: JSON.parse(B.w.document.getElementById("selftest").textContent || "null"), locRows: B.w.document.querySelectorAll("#locTab tbody tr").length, resVisible: B.w.document.getElementById("resSec").style.display !== "none" };
  fs.writeFileSync(path.join(outDir, "ui_result.json"), JSON.stringify(out, null, 1));
  console.log(JSON.stringify(out, null, 1).slice(0, 3000)); process.exit(0);
})().catch(e => { console.error(e); process.exit(1); });
