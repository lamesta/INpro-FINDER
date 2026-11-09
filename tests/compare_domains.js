// Domain-level comparison with hmmsearch: node compare_domains.js genome.fna.gz ref.domtbl ref.orfinfo out.json
const fs = require("fs"), zlib = require("zlib"), path = require("path"), C = require("../src/inafinder-core.js");
const rd = f => { const b = fs.readFileSync(f); return (f.endsWith(".gz") ? zlib.gunzipSync(b) : b).toString("utf8"); };
const [, , fna, domtbl, orfinfo, out] = process.argv;
const hmmText = fs.readFileSync(path.join(__dirname, "../models/INA_N150.hmm"), "utf8") + fs.readFileSync(path.join(__dirname, "../models/INA_C60.hmm"), "utf8");
const model = C.buildModel("bacterial", hmmText, { N: "INA_N150", C: "INA_C60" });
const contigs = C.parseFasta(rd(fna)); const t0 = Date.now();
const orfs = C.sixFrameORFs(contigs, 50); const t1 = Date.now();
const res = {}; const times = {};
for (const role of ["N", "C"]) { const s = Date.now(); res[role] = C.searchProfile(model.profiles[role], orfs); times[role] = Date.now() - s; }
const key = o => `${o.contig} ${C.orfCoords(o)}`;
const mine = [];
for (const role of ["N", "C"]) for (const h of res[role].hits) for (const d of h.domains) mine.push({ orf: key(orfs[h.orf]), len: orfs[h.orf].len, model: model.profiles[role].name, hmm_from: d.hmm_from, hmm_to: d.hmm_to, ali_from: d.ali_from, ali_to: d.ali_to, score: d.score, bias: d.bias, ievalue: d.ievalue, seqE: h.seqE });
const info = new Map(); for (const l of fs.readFileSync(orfinfo, "utf8").split("\n")) { if (!l) continue; const f = l.split(/\s+/); const d = Object.fromEntries(f.slice(1).filter(x => x.includes("=")).map(x => x.split("="))); info.set(f[0], `${d.source} ${d.coords}`); }
const ref = [];
for (const l of fs.readFileSync(domtbl, "utf8").split("\n")) { if (!l || l[0] === "#") continue; const f = l.trim().split(/\s+/); ref.push({ orf: info.get(f[0]), len: +f[2], model: f[3], seqE: +f[6], hmm_from: +f[15], hmm_to: +f[16], ali_from: +f[17], ali_to: +f[18], env_from: +f[19], env_to: +f[20], score: +f[13], bias: +f[14], ievalue: +f[12] }); }
fs.writeFileSync(out, JSON.stringify({ Z: res.N.Z, filters: { N: [res.N.nMSV, res.N.nVit, res.N.nFwd], C: [res.C.nMSV, res.C.nVit, res.C.nFwd] }, ms: { orfs: t1 - t0, ...times }, mine, ref }));
console.log(JSON.stringify({ Z: res.N.Z, filters: { N: [res.N.nMSV, res.N.nVit, res.N.nFwd, res.N.hits.length], C: [res.C.nMSV, res.C.nVit, res.C.nFwd, res.C.hits.length] }, ms: { orfs: t1 - t0, ...times }, nMine: mine.length, nRef: ref.length }));
