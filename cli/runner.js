#!/usr/bin/env node
/* runner.js: Node entry used by inafinder.py.
 *   node runner.js scan <config.json>      one genome -> result JSON (config: name, fasta, annotation, model, search, exploratory, flank, out)
 *   node runner.js finalize <config.json>  per-genome JSONs -> loci.tsv, proteins.faa, genes.fna, ina.gff3, loci_20kb.fna, genbank/*.gbk, run.json, report data
 * The core files are looked up next to this file first, then in ../src (repository layout). */
const fs = require("fs"), path = require("path"), zlib = require("zlib"), crypto = require("crypto");
const find = n => [path.join(__dirname, n), path.join(__dirname, "..", "src", n)].find(p => fs.existsSync(p));
const C = require(find("inafinder-core.js"));
const IO = require(find("inafinder-io.js"));
const rd = f => { const b = fs.readFileSync(f); return (/\.gz$/i.test(f) ? zlib.gunzipSync(b) : b).toString("utf8"); };
const sha = t => crypto.createHash("sha256").update(t).digest("hex");
const [, , cmd, cfgPath] = process.argv;
const cfg = JSON.parse(fs.readFileSync(cfgPath, "utf8"));
try {
  if (cmd === "scan") {
    const t0 = Date.now(); const text = rd(cfg.fasta); const contigs = C.parseFasta(text);
    if (!contigs.length) throw new Error("no FASTA records found");
    let annotation = null, annMeta = null; const warnings = [];
    if (cfg.annotation) {
      const at = rd(cfg.annotation); annotation = C.parseAnnotationCDS(at); const ids = Object.keys(annotation), have = new Set(contigs.map(c => c.id));
      const matched = ids.filter(i => have.has(i)).length;
      annMeta = { file: path.basename(cfg.annotation), bytes: at.length, sha256: sha(at), contigs: ids.length, matched_contigs: matched, cds: ids.reduce((a, k) => a + annotation[k].length, 0) };
      if (!ids.length) warnings.push(`annotation ${path.basename(cfg.annotation)}: no CDS features found; every locus is reported as not annotated (check that this is the right file)`);
      else if (!matched) { warnings.push(`annotation ${path.basename(cfg.annotation)}: none of its sequence ids match the genome contigs; annotation ignored`); annotation = null; }
    }
    const model = C.buildModel(cfg.model.name, cfg.model.hmmText, cfg.model.roles, cfg.model.rules);
    const result = C.scanGenome({ name: cfg.name, contigs }, model, annotation, { search: cfg.search, exploratory: cfg.exploratory });
    result.warnings = warnings;
    const seqs = IO.locusSequences(contigs, result, annotation, cfg.flank || 20000);
    fs.writeFileSync(cfg.out, JSON.stringify({ job: cfg.name, result, seqs, meta: { file: path.basename(cfg.fasta), bytes: text.length, sha256: sha(text), annotation: annMeta, wall_seconds: (Date.now() - t0) / 1000 } }));
    process.stdout.write(JSON.stringify({ ok: true, genome: cfg.name, loci: result.loci.length, orfs: result.nORF, seconds: result.seconds, warnings }) + "\n");
  } else if (cmd === "finalize") {
    const runs = cfg.inputs.map(p => JSON.parse(fs.readFileSync(p, "utf8")));
    const ok = runs.filter(r => r.result);
    const o = cfg.outdir; fs.mkdirSync(path.join(o, "genbank"), { recursive: true });
    fs.writeFileSync(path.join(o, "loci.tsv"), IO.lociTSV(ok));
    fs.writeFileSync(path.join(o, "proteins.faa"), IO.proteinFASTA(ok));
    fs.writeFileSync(path.join(o, "genes.fna"), IO.geneFASTA(ok));
    fs.writeFileSync(path.join(o, "ina.gff3"), IO.gff3(ok));
    fs.writeFileSync(path.join(o, "loci_20kb.fna"), IO.regionFASTA(ok));
    let ng = 0; ok.forEach(r => r.result.loci.forEach((l, i) => { fs.writeFileSync(path.join(o, "genbank", `${IO.locusId(r.result.genome, l).replace(/[^A-Za-z0-9_.-]/g, "_")}.gbk`), IO.genbankRegion(r.result.genome, l, r.seqs[i])); ng++; }));
    const record = { software: "INpro FINDER", version: C.VERSION, created: new Date().toISOString(), parameters: cfg.params,
      genomes: runs.map(r => r.result ? { genome: r.result.genome, file: r.meta.file, bytes: r.meta.bytes, sha256: r.meta.sha256, annotation: r.meta.annotation, contigs: r.result.nContigs, length_bp: r.result.genomeLength,
        orfs: r.result.nORF, search: r.result.search, loci: r.result.loci.length, seconds: r.result.seconds, warnings: r.result.warnings } : { job: r.job, error: r.error }) };
    fs.writeFileSync(path.join(o, "run.json"), JSON.stringify(record, null, 1));
    fs.writeFileSync(cfg.report_data, JSON.stringify({ record, runs: runs.map(r => r.result ? { job: r.job, result: r.result, seqs: r.seqs, meta: r.meta } : r) }));
    process.stdout.write(JSON.stringify({ ok: true, genomes: ok.length, loci: ok.reduce((a, r) => a + r.result.loci.length, 0), genbank: ng }) + "\n");
  } else throw new Error("unknown command " + cmd);
} catch (e) { process.stdout.write(JSON.stringify({ ok: false, error: e.message }) + "\n"); process.exit(2); }
