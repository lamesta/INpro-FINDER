// node tests/debug_region.js <genome.fna.gz> <contig> <from> <to> <model.json>
const fs = require("fs"), zlib = require("zlib"), path = require("path"), C = require("../src/inafinder-core.js");
const rd = f => { const b = fs.readFileSync(f); return (f.endsWith(".gz") ? zlib.gunzipSync(b) : b).toString("utf8"); };
const [, , fna, ctg, from, to, which] = process.argv;
const spec = JSON.parse(fs.readFileSync(which, "utf8"));
const model = C.buildModel(spec.name, spec.hmm_files.map(f => fs.readFileSync(path.resolve(path.dirname(which), f), "utf8")).join(""), spec.roles, Object.assign({}, spec.rules));
const cs = C.parseFasta(rd(fna)).filter(c => c.id === ctg);
const a = Math.max(1, +from), b = Math.min(cs[0].seq.length, +to);
const sub = [{ id: ctg, seq: cs[0].seq.slice(a - 1, b) }];
const res = C.scanGenome({ name: "dbg", contigs: sub }, model, null, { search: {} });
console.log("contig_len", cs[0].seq.length, "slice", a, b, "offset", a - 1);
for (const l of res.loci) {
  console.log(l.structural_class, l.strand, l.cds_start + a - 1, l.cds_end + a - 1, l.protein_length_aa, l.repeat_identity);
  for (const o of l.orfs) { console.log("  orf", o.strand, o.start + a - 1, o.end + a - 1, o.len, "mpos", o.mpos, o.best ? "BEST" : "");
    for (const h of o.hits) console.log("     ", h.profile, h.hmm_from, h.hmm_to, h.env_from, h.env_to, h.score, h.ievalue.toExponential(1), h.canonical); }
}
