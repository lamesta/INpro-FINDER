/* worker.js: runs one genome scan. Loaded after inafinder-core.js and inafinder-io.js, either inside a Web Worker
 * (blob URL built by app.js) or, when workers are unavailable, in the page itself (window.InaFinderWorker). */
(function () {
  const C = self.InaFinderCore, IO = self.InaFinderIO;
  async function readText(file) {
    const gz = /\.gz$/i.test(file.name || "");
    if (gz) {
      if (typeof DecompressionStream === "undefined") throw new Error(`${file.name}: this browser cannot read .gz files; decompress the file first`);
      return await new Response(file.stream().pipeThrough(new DecompressionStream("gzip"))).text();
    }
    if (typeof file === "string") return file;
    if (file.text) return await file.text();
    if (typeof FileReaderSync !== "undefined") return new FileReaderSync().readAsText(file);
    return await new Promise((res, rej) => { const r = new FileReader(); r.onload = () => res(r.result); r.onerror = () => rej(r.error); r.readAsText(file); });
  }
  async function sha256(text) {
    try { const h = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text)); return [...new Uint8Array(h)].map(b => b.toString(16).padStart(2, "0")).join(""); }
    catch (e) { return ""; }
  }
  function nucleotideFraction(contigs) {
    let n = 0, ok = 0; for (const c of contigs) { const s = c.seq.slice(0, 20000); for (let i = 0; i < s.length; i++) { n++; if ("ACGTUN".includes(s[i])) ok++; } if (n > 50000) break; }
    return n ? ok / n : 0;
  }
  async function runJob(job, post) {
    const t0 = Date.now();
    post({ type: "progress", id: job.id, stage: "reading", f: 0 });
    const text = await readText(job.genome);
    const contigs = C.parseFasta(text);
    if (!contigs.length) throw new Error("no FASTA records found");
    if (nucleotideFraction(contigs) < 0.9) throw new Error("this looks like a protein FASTA; INpro FINDER needs nucleotide genome sequences");
    let annotation = null, annMeta = null; const warnings = [];
    if (job.annotation) {
      const at = await readText(job.annotation);
      annotation = C.parseAnnotationCDS(at);
      const ids = Object.keys(annotation), have = new Set(contigs.map(c => c.id));
      const matched = ids.filter(id => have.has(id)).length;
      annMeta = { file: job.annotation.name, bytes: at.length, sha256: await sha256(at), contigs: ids.length, matched_contigs: matched, cds: ids.reduce((a, k) => a + annotation[k].length, 0) };
      if (!ids.length) warnings.push(`annotation ${job.annotation.name}: no CDS features found; every locus is reported as not annotated (check that this is the right file)`);
      else if (!matched) { warnings.push(`annotation ${job.annotation.name}: none of its sequence ids match the genome contigs; annotation ignored`); annotation = null; }
    }
    const model = C.buildModel(job.model.name, job.model.hmmText, job.model.roles, job.model.rules);
    const nProf = model.plist.length; let k = -1, last = "";
    const res = C.scanGenome({ name: job.name, contigs }, model, annotation, {
      search: job.search, exploratory: job.exploratory,
      progress: p => {
        if (p.stage === "orfs") post({ type: "progress", id: job.id, stage: `${p.n.toLocaleString()} ORFs`, f: 0.05 });
        else if (p.stage === "search") { if (p.profile !== last) { k++; last = p.profile; } post({ type: "progress", id: job.id, stage: `searching ${p.profile}`, f: 0.05 + 0.9 * (k + p.f) / nProf }); }
      } });
    res.warnings = warnings;
    const seqs = IO.locusSequences(contigs, res, annotation, job.flank || 20000);
    return { result: res, seqs, meta: { file: job.genome.name, bytes: text.length, sha256: await sha256(text), annotation: annMeta, wall_seconds: (Date.now() - t0) / 1000 } };
  }
  self.InaFinderWorker = { runJob, readText, sha256 };
  if (typeof WorkerGlobalScope !== "undefined" && self instanceof WorkerGlobalScope) {
    self.onmessage = async ev => {
      const job = ev.data;
      try { const out = await runJob(job, m => self.postMessage(m)); self.postMessage({ type: "done", id: job.id, out }); }
      catch (e) { self.postMessage({ type: "error", id: job.id, error: String((e && e.message) || e) }); }
    };
  }
})();
