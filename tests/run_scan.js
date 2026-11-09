// node run_scan.js <genome_id> <genome.fna[.gz]> <annotation.gbff[.gz]|-> <out.json> [model=bacterial|model.json] [explore]
const fs = require("fs"), zlib = require("zlib"), path = require("path"), C = require("../src/inafinder-core.js");
const rd = f => { const b = fs.readFileSync(f); return (f.endsWith(".gz") ? zlib.gunzipSync(b) : b).toString("utf8"); };
const [, , gid, fna, ann, out, which, explore] = process.argv;
const M = path.join(__dirname, "../models");
let model;
if (!which || which === "bacterial") model = C.buildModel("bacterial", fs.readFileSync(M + "/INA_N150.hmm", "utf8") + fs.readFileSync(M + "/INA_C60.hmm", "utf8"), { N: "INA_N150", C: "INA_C60" });
else { const spec = JSON.parse(fs.readFileSync(which, "utf8")); model = C.buildModel(spec.name, spec.hmm_files.map(f => fs.readFileSync(path.resolve(path.dirname(which), f), "utf8")).join(""), spec.roles, spec.rules); }
try {
  const res = C.scanGenome({ name: gid, contigs: C.parseFasta(rd(fna)) }, model, ann && ann !== "-" ? C.parseAnnotationCDS(rd(ann)) : null, { exploratory: explore === "explore" });
  for (const l of res.loci) { delete l.orfs; }
  fs.writeFileSync(out, JSON.stringify(res));
  console.log(gid, res.nORF, res.loci.length, res.seconds.toFixed(1));
} catch (e) { fs.writeFileSync(out + ".err", String(e.stack || e)); console.log(gid, "ERROR", e.message); process.exit(1); }
