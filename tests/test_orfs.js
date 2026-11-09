// Compare six-frame ORFs with esl-translate output: node test_orfs.js genome.fna[.gz] esl_orfs.faa
const fs = require("fs"), zlib = require("zlib"), C = require("../src/inafinder-core.js");
const rd = f => { const b = fs.readFileSync(f); return (f.endsWith(".gz") ? zlib.gunzipSync(b) : b).toString("utf8"); };
const contigs = C.parseFasta(rd(process.argv[2]));
const t0 = Date.now(); const orfs = C.sixFrameORFs(contigs, 50); const ms = Date.now() - t0;
const mine = new Map(orfs.map(o => [`${o.contig} ${C.orfCoords(o)} ${o.frame}`, o.aa]));
const ref = new Map(); let id = null;
for (const rec of C.parseFasta(rd(process.argv[3]))) { const m = rec.desc.match(/source=(\S+) coords=(\S+) length=(\d+) frame=(\d+)/); ref.set(`${m[1]} ${m[2]} ${m[4]}`, rec.seq); }
let same = 0, diffSeq = 0, onlyMine = [], onlyRef = [];
for (const [k, v] of mine) { if (!ref.has(k)) onlyMine.push(k); else if (ref.get(k) === v) same++; else diffSeq++; }
for (const k of ref.keys()) if (!mine.has(k)) onlyRef.push(k);
console.log(JSON.stringify({ contigs: contigs.length, mine: mine.size, ref: ref.size, identical: same, diffSeq, onlyMine: onlyMine.slice(0, 5), nOnlyMine: onlyMine.length, onlyRef: onlyRef.slice(0, 5), nOnlyRef: onlyRef.length, ms }));
