/* app.js: INpro FINDER browser interface. Must be the last script of the page (it captures the page for "save report"). */
(function () {
  "use strict";
  const C = window.InaFinderCore, IO = window.InaFinderIO, WK = window.InaFinderWorker;
  const $ = id => document.getElementById(id);
  const PRISTINE = "<!DOCTYPE html>\n" + document.documentElement.outerHTML;
  const MODELS = JSON.parse($("if-models").textContent);
  const esc = s => String(s === undefined || s === null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const fmtE = v => (v === null || v === undefined || v === "") ? "—" : (+v).toExponential(1);
  const EDGE_LABEL = { "5prime_open": "5′ end open", "3prime_open": "3′ end open", both_ends: "both ends open" };
  const START_LABEL = { Met_rule: "first Met near the first profile match (prediction)", annotation_same_stop: "annotated start (annotated gene ends at the same stop codon)",
    ORF_start_no_Met: "no Met near the first match: the whole ORF is shown (start unknown; possible intron, frameshift or contig edge)", first_Met_in_ORF: "first Met of the ORF (exploratory call)" };
  const fmtSize = b => b > 1e6 ? (b / 1e6).toFixed(1) + " MB" : (b / 1e3).toFixed(0) + " kB";
  const CLASS_LABEL = { intact_full: "intact, full length", intact_short: "intact, 500–799 aa", intact_noN_canonical_domain: "intact, no canonical N", intact_noC_canonical_domain: "intact, no canonical C",
    interrupted_N_and_C_in_separate_ORFs: "interrupted", fragment_N_only: "fragment, N only", fragment_C_only: "fragment, C only", N_and_C_low_repeat_or_short: "N + C, short or weak repeats",
    repeat_only: "repeats only", weak_repeat_hit: "weak repeat match", repeat_only_no_profile_hit: "repeats, no profile match (exploratory)" };
  const CLASS_ORDER = Object.keys(CLASS_LABEL);
  const ANNOT_LABEL = { annotated_intact_CDS: "annotated, intact", annotated_CDS_different_start_or_stop: "annotated, other start/stop", annotated_pseudo_or_partial: "annotated pseudo/partial", annotated_CDS_same_stop: "annotated, same stop", no_CDS_annotated: "not annotated", not_checked: "no annotation loaded" };
  const st = { files: [], model: "bacterial", params: null, customProfiles: [], runs: [], running: false, workers: [], sel: null, hidden: new Set(), sortKey: null, sortDir: 1, savedReport: false };
  $("verTxt").textContent = "v" + C.VERSION;
  for (const k of ["bacterial", "fungal"]) $("mdesc_" + k).textContent = MODELS[k].spec.description;
  const download = (name, data, type) => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([data], { type })); a.download = name; document.body.appendChild(a); a.click(); a.remove(); };
  window.InaFinderDownload = download;

  // ------------------------------------------------------------------ files
  const GEN_RE = /\.(fa|fna|fasta|fsa|fas|ffn|seq|fsta)(\.gz)?$/i, ANN_RE = /\.(gff3?|gbff|gbk|gb|genbank|embl)(\.gz)?$/i, HMM_RE = /\.hmm$/i;
  async function sniff(file) {
    try {
      let head;
      if (/\.gz$/i.test(file.name)) { if (typeof DecompressionStream === "undefined") return ""; const r = file.stream().pipeThrough(new DecompressionStream("gzip")).getReader(); const { value } = await r.read(); r.cancel().catch(() => { }); head = new TextDecoder().decode(value || new Uint8Array()); }
      else head = await WK.readText(file.slice ? file.slice(0, 4096) : file);
      return head.slice(0, 4096);
    } catch (e) { return ""; }
  }
  async function kindOf(file) {
    if (HMM_RE.test(file.name)) return "hmm";
    const h = await sniff(file);
    if (/^HMMER3/.test(h)) return "hmm";
    if (/^##gff-version/.test(h) || /^LOCUS\s/.test(h)) return "annotation";
    if (h.startsWith(">")) { const seq = h.split(/\r?\n/).slice(1).join("").toUpperCase(); const nt = seq.length ? [...seq].filter(c => "ACGTUN".includes(c)).length / seq.length : 1; return nt >= 0.9 ? "genome" : "protein"; }
    if (GEN_RE.test(file.name)) return "genome"; if (ANN_RE.test(file.name)) return "annotation";
    return "unknown";
  }
  const stem = n => n.replace(/\.gz$/i, "").replace(/\.[^.]+$/, "").replace(/_genomic$/i, "").replace(/\.(fna|fa|fasta|gff3?|gbff|gbk)$/i, "").toLowerCase();
  async function addFiles(list) {
    for (const f of list) {
      const kind = await kindOf(f);
      const rec = { id: "f" + Math.random().toString(36).slice(2, 9), file: f, name: f.name, size: f.size, kind, annot: null };
      if (kind === "hmm") { try { const t = await WK.readText(f); const hs = C.parseHMM(t); for (const h of hs) st.customProfiles.push({ name: h.name, M: h.M, nseq: h.nseq, text: h.text, file: f.name, role: /c(_|-)?cap|cterm|_c\d*$|C60/i.test(h.name) ? "C" : "N" }); rec.note = `${hs.length} profile(s): ${hs.map(h => h.name).join(", ")}`; } catch (e) { rec.note = "not a valid HMMER3 file: " + e.message; rec.kind = "unknown"; } }
      st.files.push(rec);
    }
    // auto-assign annotations by file-name stem (or the only annotation for the only genome)
    const gens = st.files.filter(f => f.kind === "genome"), anns = st.files.filter(f => f.kind === "annotation");
    for (const g of gens) if (!g.annot) { const m = anns.find(a => stem(a.name) === stem(g.name)); if (m) g.annot = m.id; }
    if (gens.length === 1 && anns.length === 1 && !gens[0].annot) gens[0].annot = anns[0].id;
    if (st.customProfiles.length && st.model !== "custom") { document.querySelector('input[name=model][value=custom]').checked = true; st.model = "custom"; st.params = null; }
    renderFiles(); renderParams();
  }
  function renderFiles() {
    const tb = $("fileTab").querySelector("tbody"); const anns = st.files.filter(f => f.kind === "annotation");
    $("fileTab").style.display = st.files.length ? "" : "none";
    tb.innerHTML = st.files.map(f => `<tr data-id="${f.id}"><td>${esc(f.name)}</td><td>${fmtSize(f.size)}</td><td>${f.kind === "protein" ? '<span class="warn">protein FASTA (not usable)</span>' : esc(f.kind)}${f.note ? `<div class="small">${esc(f.note)}</div>` : ""}</td>
      <td>${f.kind === "genome" ? `<select data-ann="${f.id}"><option value="">none</option>${anns.map(a => `<option value="${a.id}" ${f.annot === a.id ? "selected" : ""}>${esc(a.name)}</option>`).join("")}</select>` : ""}</td>
      <td><button data-rm="${f.id}" title="remove">✕</button></td></tr>`).join("");
    tb.querySelectorAll("select[data-ann]").forEach(s => s.onchange = () => { st.files.find(f => f.id === s.dataset.ann).annot = s.value || null; });
    tb.querySelectorAll("button[data-rm]").forEach(b => b.onclick = () => { const f = st.files.find(x => x.id === b.dataset.rm); st.files = st.files.filter(x => x !== f); if (f.kind === "hmm") { st.customProfiles = st.customProfiles.filter(p => p.file !== f.name); st.params = null; renderParams(); } st.files.forEach(x => { if (x.annot === f.id) x.annot = null; }); renderFiles(); });
    const n = st.files.filter(f => f.kind === "genome").length;
    $("fileMsg").textContent = n ? `${n} genome(s) ready.` : (st.files.length ? "No genome FASTA loaded yet." : "");
    $("runBtn").disabled = !n || st.running || (st.model === "custom" && !st.customProfiles.some(p => p.role !== "ignore"));
  }
  const drop = $("drop");
  drop.addEventListener("dragover", e => { e.preventDefault(); drop.classList.add("over"); });
  drop.addEventListener("dragleave", () => drop.classList.remove("over"));
  drop.addEventListener("drop", e => { e.preventDefault(); drop.classList.remove("over"); addFiles([...e.dataTransfer.files]); });
  $("fileIn").onchange = e => { addFiles([...e.target.files]); e.target.value = ""; };

  // ------------------------------------------------------------------ model and parameters
  const GEN_FIELDS = [["min_orf_aa", "minimum ORF length (aa)"], ["candidate_min_hits", "weak matches for a candidate (n)"], ["candidate_min_len", "…with ORF length ≥ (aa)"], ["repeat_identity_min", "16-aa repeat identity ≥"],
    ["intact_full_min_aa", "intact, full length ≥ (aa)"], ["intact_short_min_aa", "intact, short ≥ (aa)"], ["repeat_only_min_aa", "repeats only ≥ (aa)"], ["locus_gap_bp", "ORFs in one locus within (bp)"],
    ["met_upstream_aa", "start Met: search from first match − (aa)"], ["met_downstream_aa", "…up to first match + (aa)"], ["met_fallback", "no Met near first match: start at"],
    ["exploratory_min_aa", "exploratory: protein length ≥ (aa)"], ["exploratory_identity_min", "exploratory: 16-aa repeat identity ≥"], ["max_period4_ratio", "exploratory: 4-aa / 16-aa identity ratio ≤"]];
  const GEN_DEFAULTS = Object.fromEntries(Object.entries(C.BACTERIAL_RULES).filter(([k]) => k !== "N" && k !== "C"));
  const SEARCH_FIELDS = [["F1", "MSV filter P"], ["F2", "Viterbi filter P"], ["F3", "Forward filter P"], ["E", "report sequences E ≤"], ["domE", "report domains c-E ≤"]];
  function defaultParams(kind) {
    if (kind === "custom") {
      const roles = st.customProfiles.map(p => ({ role: p.role, profile: p.name, rule: p.role === "N" ? { hmm_from_max: 10, hmm_to_min: Math.round(0.66 * p.M), ievalue_max: 1e-10, env_from_max: 300 } : { hmm_from_max: 15, hmm_to_min: Math.round(0.8 * p.M), ievalue_max: 1e-10, last_aa: 200 } }));
      return { roles, rules: Object.assign({}, GEN_DEFAULTS, JSON.parse(JSON.stringify(MODELS.bacterial.spec.rules))), search: Object.assign({}, C.SEARCH_DEFAULTS), name: "custom" };
    }
    const sp = MODELS[kind].spec; return { roles: JSON.parse(JSON.stringify(sp.roles)), rules: Object.assign({}, GEN_DEFAULTS, JSON.parse(JSON.stringify(sp.rules))), search: Object.assign({}, C.SEARCH_DEFAULTS), name: sp.name };
  }
  function renderParams() {
    if (!st.params) st.params = defaultParams(st.model);
    const P = st.params; const custom = st.model === "custom";
    let h = `<table class="tab"><thead><tr><th>profile</th><th>role</th><th>alignment starts at profile position ≤</th><th>…and reaches position ≥</th><th>i-Evalue ≤</th><th>position in protein</th></tr></thead><tbody>`;
    P.roles.forEach((r, i) => {
      const meta = custom ? st.customProfiles.find(p => p.name === r.profile) : null;
      h += `<tr><td><b>${esc(r.profile)}</b>${meta ? `<div class="small">${meta.M} positions, ${meta.nseq} sequences, ${esc(meta.file)}</div>` : ""}</td>
        <td>${custom ? `<select data-role="${i}"><option value="N" ${r.role === "N" ? "selected" : ""}>N-terminal domain</option><option value="C" ${r.role === "C" ? "selected" : ""}>C-terminal cap</option><option value="ignore" ${r.role === "ignore" ? "selected" : ""}>ignore</option></select>` : (r.role === "N" ? "N-terminal domain" : "C-terminal cap")}</td>
        <td><input data-rule="${i}" data-k="hmm_from_max" value="${r.rule.hmm_from_max}" size="5"></td><td><input data-rule="${i}" data-k="hmm_to_min" value="${r.rule.hmm_to_min}" size="5"></td>
        <td><input data-rule="${i}" data-k="ievalue_max" value="${r.rule.ievalue_max}" size="7"></td>
        <td>${r.role === "N" ? `starts within first <input data-rule="${i}" data-k="env_from_max" value="${r.rule.env_from_max ?? 300}" size="5"> aa` : `ends within last <input data-rule="${i}" data-k="last_aa" value="${r.rule.last_aa ?? 200}" size="5"> aa`}</td></tr>`;
    });
    h += "</tbody></table>";
    if (custom && !P.roles.length) h = `<div class="small warn">Drop one or more .hmm files (HMMER3, built with hmmbuild) into the file area to use custom profiles.</div>`;
    $("profTabWrap").innerHTML = h;
    $("ruleGrid").innerHTML = GEN_FIELDS.map(([k, lab]) => k === "met_fallback"
      ? `<label>${lab}<select data-gen="${k}"><option value="first_M" ${P.rules[k] !== "orf_start" ? "selected" : ""}>first Met in the ORF</option><option value="orf_start" ${P.rules[k] === "orf_start" ? "selected" : ""}>ORF start (keep whole ORF)</option></select></label>`
      : `<label>${lab}<input data-gen="${k}" value="${P.rules[k]}"></label>`).join("") +
      SEARCH_FIELDS.map(([k, lab]) => `<label>${lab}<input data-search="${k}" value="${P.search[k]}"></label>`).join("");
    $("profTabWrap").querySelectorAll("input[data-rule]").forEach(inp => inp.onchange = () => { P.roles[+inp.dataset.rule].rule[inp.dataset.k] = +inp.value; });
    $("profTabWrap").querySelectorAll("select[data-role]").forEach(s => s.onchange = () => { const r = P.roles[+s.dataset.role]; r.role = s.value; const cp = st.customProfiles.find(p => p.name === r.profile); if (cp) cp.role = s.value; if (s.value === "N") r.rule = { hmm_from_max: 10, hmm_to_min: Math.round(0.66 * cp.M), ievalue_max: 1e-10, env_from_max: 300 }; else if (s.value === "C") r.rule = { hmm_from_max: 15, hmm_to_min: Math.round(0.8 * cp.M), ievalue_max: 1e-10, last_aa: 200 }; renderParams(); renderFiles(); });
    $("ruleGrid").querySelectorAll("input[data-gen]").forEach(inp => inp.onchange = () => { P.rules[inp.dataset.gen] = +inp.value; });
    $("ruleGrid").querySelectorAll("select[data-gen]").forEach(sel => sel.onchange = () => { P.rules[sel.dataset.gen] = sel.value; });
    $("ruleGrid").querySelectorAll("input[data-search]").forEach(inp => inp.onchange = () => { P.search[inp.dataset.search] = +inp.value; });
    renderFiles();
  }
  document.querySelectorAll("input[name=model]").forEach(r => r.onchange = () => { st.model = r.value; st.params = null; renderParams(); });
  $("resetParams").onclick = () => { st.params = null; renderParams(); };
  function modelForRun() {
    const P = st.params; const roles = P.roles.filter(r => r.role !== "ignore");
    let hmmText;
    if (st.model === "custom") hmmText = roles.map(r => st.customProfiles.find(p => p.name === r.profile).text).join("");
    else { const m = MODELS[st.model]; hmmText = m.spec.hmm_files.map(f => m.hmm[f]).join(""); }
    return { name: P.name, hmmText, roles, rules: P.rules };
  }

  // ------------------------------------------------------------------ run
  function workerSource() { return ["src-core", "src-io", "src-worker"].map(id => $(id).textContent).join("\n;\n"); }
  function makeWorker() { const url = URL.createObjectURL(new Blob([workerSource()], { type: "text/javascript" })); return new Worker(url); }
  const progRows = {};
  function setProg(id, f, text, cls) { const r = progRows[id]; if (!r) return; r.querySelector(".bar>span").style.width = Math.round(100 * Math.min(1, f)) + "%"; r.querySelector(".pt").textContent = text; if (cls) r.className = "prog " + cls; }
  async function run(opts = {}) {
    const gens = st.files.filter(f => f.kind === "genome"); if (!gens.length) return;
    const model = modelForRun();
    try { C.buildModel(model.name, model.hmmText, model.roles, model.rules); } catch (e) { $("runInfo").textContent = "model error: " + e.message; return; }
    st.running = true; st.runs = []; st.sel = null; $("runBtn").disabled = true; $("stopBtn").disabled = false; st.savedReport = false;
    const runParams = { model: st.model, model_name: model.name, roles: model.roles, rules: model.rules, search: st.params.search, exploratory: $("optExplore").checked, started: new Date().toISOString(), hmm_sha256: await WK.sha256(model.hmmText) };
    st.runParams = runParams;
    $("progress").innerHTML = gens.map(g => `<div class="prog" id="p_${g.id}"><span style="min-width:260px">${esc(g.name)}</span><span class="bar"><span style="width:0"></span></span><span class="pt small">waiting</span></div>`).join("");
    gens.forEach(g => progRows[g.id] = $("p_" + g.id));
    const jobs = gens.map(g => ({ id: g.id, name: g.name.replace(/\.gz$/i, "").replace(/\.[^.]+$/, ""), genome: g.file, annotation: g.annot ? st.files.find(f => f.id === g.annot).file : null,
      model, search: st.params.search, exploratory: runParams.exploratory, flank: 20000 }));
    const t0 = Date.now(); let done = 0;
    const finish = (job, out, err) => {
      done++;
      if (err) { setProg(job.id, 1, "error: " + err, "warn"); st.runs.push({ job: job.name, error: err }); }
      else { st.runs.push(Object.assign({ job: job.name }, out)); setProg(job.id, 1, `${out.result.loci.length} locus/loci · ${out.result.seconds.toFixed(1)} s${out.result.warnings && out.result.warnings.length ? " · " + out.result.warnings.join("; ") : ""}`); renderResults(); }
      $("runInfo").textContent = `${done} / ${jobs.length} genomes · ${((Date.now() - t0) / 1000).toFixed(0)} s`;
    };
    let useWorkers = !opts.mainThread && typeof Worker !== "undefined";
    if (useWorkers) { try { const w = makeWorker(); w.terminate(); } catch (e) { useWorkers = false; } }
    if (useWorkers) {
      const nW = Math.max(1, Math.min(jobs.length, Math.min(6, (navigator.hardwareConcurrency || 2) - 1)));
      const queue = jobs.slice();
      await new Promise(resolve => {
        let active = 0;
        const next = w => {
          const job = queue.shift(); if (!job || !st.running) { w.terminate(); if (--active === 0) resolve(); return; }
          w.onmessage = ev => { const m = ev.data; if (m.type === "progress") setProg(m.id, m.f, m.stage); else { finish(job, m.out, m.type === "error" ? m.error : null); next(w); } };
          w.onerror = ev => { finish(job, null, ev.message || "worker error"); next(w); };
          w.postMessage(job);
        };
        st.workers = []; for (let i = 0; i < nW; i++) { const w = makeWorker(); st.workers.push(w); active++; next(w); }
      });
    } else {
      for (const job of jobs) {
        if (!st.running) break;
        await new Promise(r => setTimeout(r, 0));
        try { const out = await WK.runJob(job, m => setProg(m.id, m.f, m.stage)); finish(job, out, null); } catch (e) { finish(job, null, e.message || String(e)); }
      }
    }
    st.running = false; $("stopBtn").disabled = true; renderFiles();
    $("runInfo").textContent = `finished: ${st.runs.filter(r => !r.error).length} of ${jobs.length} genomes in ${((Date.now() - t0) / 1000).toFixed(0)} s${useWorkers ? "" : " (single thread)"}`;
    renderResults();
    return st.runs;
  }
  $("runBtn").onclick = () => run();
  $("stopBtn").onclick = () => { st.running = false; st.workers.forEach(w => w.terminate()); st.workers = []; $("runInfo").textContent = "stopped"; $("stopBtn").disabled = true; renderFiles(); };

  // ------------------------------------------------------------------ results
  const okRuns = () => st.runs.filter(r => !r.error && r.result);
  function allLoci() { const out = []; okRuns().forEach((r, ri) => r.result.loci.forEach((l, li) => out.push({ r, ri, li, l, genome: r.result.genome }))); return out; }
  function renderResults() {
    const runs = okRuns(); if (!runs.length && !st.runs.length) return;
    $("resSec").style.display = "";
    const rp = st.runParams || {};
    $("resHead").innerHTML = `Model <b>${esc(rp.model_name || "")}</b> · ${runs.length} genome(s) · ${allLoci().length} locus/loci · INpro FINDER ${C.VERSION}${st.savedReport ? " · saved report: load new genomes above to start a new scan" : ""}`;
    $("genTab").querySelector("tbody").innerHTML = st.runs.map(r => r.error ? `<tr><td>${esc(r.job)}</td><td colspan="6" class="warn">${esc(r.error)}</td></tr>` : (() => {
      const res = r.result, cnt = {}; res.loci.forEach(l => cnt[l.structural_class] = (cnt[l.structural_class] || 0) + 1);
      return `<tr><td>${esc(res.genome)}${r.meta && r.meta.annotation ? `<div class="small">annotation: ${esc(r.meta.annotation.file)} (${r.meta.annotation.matched_contigs}/${r.meta.annotation.contigs} contigs matched)</div>` : ""}</td><td>${res.nContigs}</td><td>${(res.genomeLength / 1e6).toFixed(2)}</td><td>${res.nORF.toLocaleString()}</td><td>${res.loci.length}</td>
        <td>${Object.entries(cnt).map(([k, v]) => `<span class="badge c-${k}">${esc(CLASS_LABEL[k] || k)} ${v}</span>`).join(" ")}</td><td>${res.seconds.toFixed(1)}</td></tr>`; })()).join("");
    const L = allLoci(); const cnt = {}; L.forEach(x => cnt[x.l.structural_class] = (cnt[x.l.structural_class] || 0) + 1);
    $("filterRow").innerHTML = "<span class='small'>show:</span>" + CLASS_ORDER.filter(k => cnt[k]).map(k => `<label class="small"><input type="checkbox" data-cls="${k}" ${st.hidden.has(k) ? "" : "checked"}> <span class="badge c-${k}">${esc(CLASS_LABEL[k])}</span> ${cnt[k]}</label>`).join(" ");
    $("filterRow").querySelectorAll("input[data-cls]").forEach(c => c.onchange = () => { if (c.checked) st.hidden.delete(c.dataset.cls); else st.hidden.add(c.dataset.cls); renderLoci(); });
    renderLoci();
  }
  function renderLoci() {
    let L = allLoci().filter(x => !st.hidden.has(x.l.structural_class));
    if (st.sortKey) { const k = st.sortKey, d = st.sortDir; const val = x => k === "genome" ? x.genome : k === "pos" ? x.l.contig + ":" + String(x.l.cds_start).padStart(10, "0") : (x.l[k] === null || x.l[k] === undefined ? Infinity : x.l[k]);
      L.sort((a, b) => { const va = val(a), vb = val(b); return (va < vb ? -1 : va > vb ? 1 : 0) * d; }); }
    $("locTab").querySelector("tbody").innerHTML = L.map(x => { const l = x.l, q = x.r.seqs && x.r.seqs[x.li];
      return `<tr data-ri="${x.ri}" data-li="${x.li}" class="${st.sel && st.sel.ri === x.ri && st.sel.li === x.li ? "sel" : ""}"><td>${esc(x.genome)}</td><td>${l.locus_index}</td><td>${esc(l.contig)}:${(q ? q.cds_start : l.cds_start).toLocaleString()}–${(q ? q.cds_end : l.cds_end).toLocaleString()} (${l.strand})</td>
      <td>${l.protein_length_aa}</td><td><span class="badge c-${l.structural_class}">${esc(CLASS_LABEL[l.structural_class] || l.structural_class)}</span></td><td>${l.repeat_identity.toFixed(2)} <span class="small">(${l.repeat_blocks_16aa})</span></td><td class="small">${esc(EDGE_LABEL[l.contig_edge] || "")}</td>
      <td>${fmtE(l.N_ievalue)}</td><td>${fmtE(l.C_ievalue)}</td><td>${esc(ANNOT_LABEL[l.annotation_state] || l.annotation_state)}${l.annot_locus_tag ? `<div class="small">${esc(l.annot_locus_tag)} ${esc(l.annot_product)}</div>` : ""}</td></tr>`; }).join("");
    $("locTab").querySelectorAll("tbody tr").forEach(tr => tr.onclick = () => { st.sel = { ri: +tr.dataset.ri, li: +tr.dataset.li }; renderLoci(); showDetail(); });
  }
  $("locTab").querySelectorAll("th").forEach(th => th.onclick = () => { const k = th.dataset.k; if (st.sortKey === k) st.sortDir *= -1; else { st.sortKey = k; st.sortDir = 1; } $("locTab").querySelectorAll("th").forEach(t => t.classList.toggle("sorted", t === th)); renderLoci(); });
  function selected() { if (!st.sel) return null; const r = okRuns()[st.sel.ri]; return r ? { r, l: r.result.loci[st.sel.li], q: r.seqs && r.seqs[st.sel.li] } : null; }
  function showDetail() {
    const s = selected(); if (!s) { $("detail").innerHTML = ""; return; } const { r, l, q } = s;
    const title = `${r.result.genome} · locus ${l.locus_index}`;
    const svg = l.orfs ? IO.architectureSVG(l, { title, width: Math.min(1180, Math.max(700, $("resSec").clientWidth - 40 || 1000)) }) : "";
    const hits = (l.orfs || []).flatMap((o, oi) => (o.hits || []).map(h => ({ o: oi + 1, h }))).filter(x => x.h.canonical || x.h.score > 15);
    $("detail").innerHTML = `<div class="det">${svg}<div class="row"><span class="small">save figure:</span><button data-fig="png">PNG</button><button data-fig="svg">SVG</button><button data-fig="pdf">PDF</button><span id="figMsg" class="small"></span></div>
      <div class="small">Start codon: ${esc(START_LABEL[l.start_source] || l.start_source)}. Repeat consensus: <code>${esc(l.repeat_consensus)}</code> (${l.repeat_blocks_16aa} blocks, identity ${l.repeat_identity}; 4-aa / 16-aa identity ratio ${l.period4_ratio ?? "n/a"}, INA about 0.5).
      ${l.contig_edge ? ` The ORF runs into the contig end without a stop codon (${esc(EDGE_LABEL[l.contig_edge])}): the gene may continue beyond the assembled sequence.` : ""}
      ${l.exploratory ? ` <b>Exploratory call</b>: no INA profile match; 16-aa repeats only. Treat as a hypothesis.` : ""}
      ${q ? ` Neighbouring genes in the ±20 kb export: ${q.neighbours.length} (${esc(q.neighbour_source)}).` : ""}</div>
      <table class="tab"><thead><tr><th>ORF</th><th>profile</th><th>profile positions</th><th>protein positions</th><th>bits</th><th>i-Evalue</th><th>canonical</th></tr></thead><tbody>
      ${hits.map(x => `<tr><td>${x.o}</td><td>${esc(x.h.profile)}</td><td>${x.h.hmm_from}–${x.h.hmm_to}</td><td>${x.h.env_from}–${x.h.env_to}</td><td>${x.h.score}</td><td>${fmtE(x.h.ievalue)}</td><td>${x.h.canonical ? "yes" : ""}</td></tr>`).join("")}</tbody></table>
      <div class="small">Protein (${l.protein_length_aa} aa) <button id="cpProt">copy</button></div><div class="seq">${esc(l.protein)}</div></div>`;
    $("cpProt").onclick = () => { if (navigator.clipboard) navigator.clipboard.writeText(l.protein); };
    $("detail").querySelectorAll("button[data-fig]").forEach(b => b.onclick = () => saveFigure(svg, `${r.result.genome}_locus${l.locus_index}_architecture`, b.dataset.fig, title));
  }
  async function saveFigure(svg, base, fmt, title) {
    const m = svg.match(/width="(\d+)" height="(\d+)"/); const w = +m[1], h = +m[2];
    if (fmt === "svg") return download(base + ".svg", svg, "image/svg+xml");
    $("figMsg").textContent = "rendering…";
    const img = new Image(); const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" }));
    await new Promise((res, rej) => { img.onload = res; img.onerror = () => rej(new Error("render failed")); img.src = url; });
    const sc = 4, cv = document.createElement("canvas"); cv.width = w * sc; cv.height = h * sc; const g = cv.getContext("2d"); g.fillStyle = "#fff"; g.fillRect(0, 0, cv.width, cv.height); g.scale(sc, sc); g.drawImage(img, 0, 0, w, h);
    if (fmt === "png") { const b = await new Promise(r => cv.toBlob(r, "image/png")); download(base + ".png", b, "image/png"); }
    else {
      let data, filter;
      if (typeof CompressionStream !== "undefined") { const px = g.getImageData(0, 0, cv.width, cv.height).data, rgb = new Uint8Array(cv.width * cv.height * 3); for (let i = 0, j = 0; i < px.length; i += 4) { rgb[j++] = px[i]; rgb[j++] = px[i + 1]; rgb[j++] = px[i + 2]; }
        const cs = new CompressionStream("deflate"); const wr = cs.writable.getWriter(); wr.write(rgb); wr.close(); data = new Uint8Array(await new Response(cs.readable).arrayBuffer()); filter = "FlateDecode"; }
      else { const b = await new Promise(r => cv.toBlob(r, "image/jpeg", 0.95)); data = new Uint8Array(await b.arrayBuffer()); filter = "DCTDecode"; }
      download(base + ".pdf", IO.buildPdf({ data, filter, imgW: cv.width, imgH: cv.height, pageW: w * 0.75, pageH: h * 0.75, title }), "application/pdf");
    }
    $("figMsg").textContent = "";
  }

  // ------------------------------------------------------------------ downloads and report
  function runRecord() {
    return { software: "INpro FINDER", version: C.VERSION, created: new Date().toISOString(), parameters: st.runParams,
      genomes: st.runs.map(r => r.error ? { job: r.job, error: r.error } : { genome: r.result.genome, file: r.meta.file, bytes: r.meta.bytes, sha256: r.meta.sha256, annotation: r.meta.annotation,
        contigs: r.result.nContigs, length_bp: r.result.genomeLength, orfs: r.result.nORF, search: r.result.search, loci: r.result.loci.length, seconds: r.result.seconds, warnings: r.result.warnings }) };
  }
  const stamp = () => new Date().toISOString().slice(0, 10);
  function saveReport() {
    const doc = new DOMParser().parseFromString(PRISTINE, "text/html");
    const old = doc.getElementById("if-results"); if (old) old.remove();
    const sc = doc.createElement("script"); sc.type = "application/json"; sc.id = "if-results";
    sc.textContent = JSON.stringify({ record: runRecord(), runs: st.runs.map(r => r.error ? r : { job: r.job, result: r.result, seqs: r.seqs, meta: r.meta }) }).replace(/</g, "\\u003c");
    const first = doc.getElementById("src-core"); first.parentNode.insertBefore(sc, first);
    const self = doc.getElementById("selftest"); if (self) self.textContent = "";
    return "<!DOCTYPE html>\n" + doc.documentElement.outerHTML;
  }
  const DL = {
    tsv: () => download(`ina_finder_loci_${stamp()}.tsv`, IO.lociTSV(okRuns()), "text/tab-separated-values"),
    faa: () => download(`ina_finder_proteins_${stamp()}.faa`, IO.proteinFASTA(okRuns()), "text/plain"),
    fna: () => download(`ina_finder_genes_${stamp()}.fna`, IO.geneFASTA(okRuns()), "text/plain"),
    gff: () => download(`ina_finder_${stamp()}.gff3`, IO.gff3(okRuns()), "text/plain"),
    regions: () => download(`ina_finder_loci_20kb_${stamp()}.fna`, IO.regionFASTA(okRuns()), "text/plain"),
    gbk: () => { const files = []; okRuns().forEach(r => r.result.loci.forEach((l, i) => { if (r.seqs && r.seqs[i]) files.push({ name: `${IO.locusId(r.result.genome, l).replace(/[^A-Za-z0-9_.-]/g, "_")}.gbk`, data: IO.genbankRegion(r.result.genome, l, r.seqs[i]) }); }));
      download(`ina_finder_loci_genbank_${stamp()}.zip`, IO.zipStore(files), "application/zip"); },
    json: () => download(`ina_finder_run_${stamp()}.json`, JSON.stringify(runRecord(), null, 1), "application/json"),
    report: () => download(`ina_finder_report_${stamp()}.html`, saveReport(), "text/html"),
  };
  document.querySelectorAll("button[data-dl]").forEach(b => b.onclick = () => DL[b.dataset.dl]());

  // ------------------------------------------------------------------ saved report: show results directly
  const emb = $("if-results");
  if (emb && emb.textContent.trim()) {
    try { const d = JSON.parse(emb.textContent); st.runs = d.runs; st.runParams = d.record.parameters; st.savedReport = true; renderResults(); } catch (e) { console.error(e); }
  }
  renderParams();
  window.InaFinderApp = { st, addFiles, run, renderResults, saveReport, runRecord, DL, selectLocus: (ri, li) => { st.sel = { ri, li }; renderLoci(); showDetail(); } };
  if (location.hash === "#selftest") $("selftest").textContent = JSON.stringify({ ok: true, version: C.VERSION, models: Object.keys(MODELS), savedRuns: st.runs.length });
})();
