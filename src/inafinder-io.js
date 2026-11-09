/* inafinder-io.js: output formats for INpro FINDER results (shared by the browser tool and the command-line version).
 * Locus sequences and flanks, loci TSV, protein / gene FASTA, GFF3, GenBank regions for clinker, one-page PDF writer. */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory(require("./inafinder-core.js"));
  else root.InaFinderIO = factory(root.InaFinderCore);
}(typeof self !== "undefined" ? self : this, function (C) {
  "use strict";
  const CMP = { A: "T", C: "G", G: "C", T: "A", R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K", B: "V", V: "B", D: "H", H: "D", N: "N" };
  function revcomp(s) { let o = ""; for (let i = s.length - 1; i >= 0; i--) o += CMP[s[i]] || "N"; return o; }
  const CODE = (() => { const b = "TCAG", aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; const m = {}; let i = 0; for (const x of b) for (const y of b) for (const z of b) m[x + y + z] = aa[i++]; return m; })();
  function translate(nt) { let p = ""; for (let i = 0; i + 2 < nt.length; i += 3) p += CODE[nt.slice(i, i + 3)] || "X"; return p; }
  const STARTS = new Set(["ATG", "GTG", "TTG"]);

  /** Simple gene calls for a region without annotation: ORFs >= minCodons with ATG/GTG/TTG start (longest per stop),
   *  then greedy non-overlapping selection by length (overlaps <= maxOverlap bp allowed). Coordinates 1-based in the region. */
  function simpleGeneCalls(nt, minCodons = 100, maxOverlap = 60) {
    const cands = [];
    for (const strand of ["+", "-"]) {
      const s = strand === "+" ? nt : revcomp(nt), L = s.length;
      for (let f = 0; f < 3; f++) {
        let start = -1;
        for (let i = f; i + 2 < L; i += 3) {
          const cod = s.slice(i, i + 3);
          if (start < 0 && STARTS.has(cod)) start = i;
          if (cod === "TAA" || cod === "TAG" || cod === "TGA") {
            if (start >= 0 && (i - start) / 3 >= minCodons) {
              const a = start + 1, b = i + 3; // includes stop
              cands.push(strand === "+" ? { s: a, e: b, st: "+" } : { s: L - b + 1, e: L - a + 1, st: "-" });
            }
            start = -1;
          }
        }
      }
    }
    cands.sort((x, y) => (y.e - y.s) - (x.e - x.s));
    const chosen = [];
    for (const c of cands) if (chosen.every(o => Math.min(o.e, c.e) - Math.max(o.s, c.s) + 1 <= maxOverlap)) chosen.push(c);
    return chosen.sort((x, y) => x.s - y.s);
  }

  /** For every locus of a scan result: CDS nucleotide sequence (with stop codon when present), protein, and the
   *  +/- flank region with neighbouring genes (annotation CDS when available, otherwise simple gene calls). */
  function locusSequences(contigs, res, annotation, flank = 20000) {
    const byId = new Map(contigs.map(c => [c.id, c.seq]));
    return res.loci.map(l => {
      const seq = byId.get(l.contig); const L = seq.length;
      let cs = l.cds_start, ce = l.cds_end;
      // add the stop codon if it follows the ORF end
      if (l.strand === "+" && ce + 3 <= L && ["TAA", "TAG", "TGA"].includes(seq.slice(ce, ce + 3))) ce += 3;
      if (l.strand === "-" && cs - 3 >= 1 && ["TTA", "CTA", "TCA"].includes(seq.slice(cs - 4, cs - 1))) cs -= 3;
      const ntF = seq.slice(cs - 1, ce); const cds_nt = l.strand === "+" ? ntF : revcomp(ntF);
      const rs = Math.max(1, cs - flank), re = Math.min(L, ce + flank);
      const region = seq.slice(rs - 1, re);
      let neighbours, source;
      const feats = annotation && annotation[l.contig];
      if (feats) {
        source = "annotation";
        neighbours = feats.filter(f => f.e >= rs && f.s <= re && !(f.e - f.s > 50000)).map(f => {
          const s0 = Math.max(rs, f.s), e0 = Math.min(re, f.e); const nt = seq.slice(s0 - 1, e0);
          const prot = (f.st === "+" ? translate(nt) : translate(revcomp(nt))).replace(/\*$/, "");
          return { s: s0 - rs + 1, e: e0 - rs + 1, st: f.st, tag: f.tag, prod: f.prod, pseudo: f.pseudo, partial: f.partial || s0 !== f.s || e0 !== f.e, translation: prot };
        });
      } else {
        source = "simple ORF calls (ATG/GTG/TTG start, >= 100 codons)";
        neighbours = simpleGeneCalls(region).map(g => { const nt = region.slice(g.s - 1, g.e); return { s: g.s, e: g.e, st: g.st, tag: "", prod: "ORF", translation: (g.st === "+" ? translate(nt) : translate(revcomp(nt))).replace(/\*$/, "") }; });
      }
      // drop neighbours that are the INA gene itself (>= 50% overlap, same strand)
      const ia = cs - rs + 1, ib = ce - rs + 1;
      neighbours = neighbours.filter(g => !(g.st === l.strand && Math.min(g.e, ib) - Math.max(g.s, ia) + 1 >= 0.5 * (g.e - g.s + 1)));
      return { cds_start: cs, cds_end: ce, cds_nt, protein: l.protein, region: { contig: l.contig, start: rs, end: re, nt: region }, ina_in_region: { s: ia, e: ib, st: l.strand }, neighbours, neighbour_source: source };
    });
  }

  const TSV_COLS = ["genome", "locus_id", "contig", "strand", "cds_start", "cds_end", "protein_length_aa", "structural_class", "repeat_identity", "repeat_blocks_16aa", "repeat_consensus",
    "N_canonical", "N_profile", "N_ievalue", "C_canonical", "C_profile", "C_ievalue", "n_orfs", "locus_start", "locus_end", "start_source", "contig_edge", "period4_ratio", "annotation_state", "annot_locus_tag", "annot_protein_id", "annot_product", "model"];
  const fmtE = v => (v === null || v === undefined || v === "") ? "" : (+v).toExponential(2);
  function locusId(genome, l) { return `${genome}|L${l.locus_index}`; }
  function lociTSV(runs) {
    const lines = [TSV_COLS.join("\t")];
    for (const r of runs) r.result.loci.forEach((l, i) => {
      const q = r.seqs && r.seqs[i];
      const row = { genome: r.result.genome, locus_id: locusId(r.result.genome, l), model: r.result.model, ...l, N_ievalue: fmtE(l.N_ievalue), C_ievalue: fmtE(l.C_ievalue),
        cds_start: q ? q.cds_start : l.cds_start, cds_end: q ? q.cds_end : l.cds_end };
      lines.push(TSV_COLS.map(c => (row[c] === undefined || row[c] === null ? "" : String(row[c]).replace(/[\t\n]/g, " "))).join("\t"));
    });
    return lines.join("\n") + "\n";
  }
  function wrap(s, n = 60) { let o = ""; for (let i = 0; i < s.length; i += n) o += s.slice(i, i + n) + "\n"; return o; }
  function proteinFASTA(runs) { let o = ""; for (const r of runs) r.result.loci.forEach(l => { o += `>${locusId(r.result.genome, l)} ${l.contig}:${l.cds_start}-${l.cds_end}(${l.strand}) ${l.structural_class} len=${l.protein_length_aa} repeat_identity=${l.repeat_identity}\n` + wrap(l.protein); }); return o; }
  function geneFASTA(runs) { let o = ""; for (const r of runs) (r.seqs || []).forEach((q, i) => { const l = r.result.loci[i]; o += `>${locusId(r.result.genome, l)} ${l.contig}:${q.cds_start}-${q.cds_end}(${l.strand}) ${l.structural_class}\n` + wrap(q.cds_nt); }); return o; }
  function regionFASTA(runs) { let o = ""; for (const r of runs) (r.seqs || []).forEach((q, i) => { const l = r.result.loci[i]; o += `>${locusId(r.result.genome, l)}_region ${q.region.contig}:${q.region.start}-${q.region.end} INA ${l.structural_class}\n` + wrap(q.region.nt); }); return o; }
  /** Product text by confidence: profile-supported classes vs low-confidence repeat classes. */
  function productOf(l) {
    if (l.structural_class === "repeat_only_no_profile_hit") return "16-aa repeat protein without INA profile match (exploratory, low confidence)";
    if (l.structural_class === "weak_repeat_hit") return "weak INA profile and repeat match (low confidence)";
    if (l.structural_class === "repeat_only") return "INA-like repeat protein without canonical terminal domains";
    return "ice nucleation protein (predicted)";
  }
  function gff3(runs, version) {
    let o = "##gff-version 3\n";
    for (const r of runs) {
      o += `# genome ${r.result.genome}; model ${r.result.model}; INpro FINDER ${version || C.VERSION}\n`;
      r.result.loci.forEach((l, i) => {
        const q = r.seqs && r.seqs[i]; const s = q ? q.cds_start : l.cds_start, e = q ? q.cds_end : l.cds_end; const id = locusId(r.result.genome, l).replace(/\|/g, "_");
        const attr = [`ID=${id}`, `Name=ina_${l.locus_index}`, `product=${productOf(l)}`, `structural_class=${l.structural_class}`, `repeat_identity=${l.repeat_identity}`, `period4_ratio=${l.period4_ratio ?? ""}`, `contig_edge=${l.contig_edge || "none"}`,
          `N_ievalue=${fmtE(l.N_ievalue)}`, `C_ievalue=${fmtE(l.C_ievalue)}`, `start_source=${l.start_source}`, `annotation_state=${l.annotation_state}`].join(";");
        o += [l.contig, "INAfinder", "gene", s, e, ".", l.strand, ".", `ID=${id}.gene;Name=ina_${l.locus_index}`].join("\t") + "\n";
        o += [l.contig, "INAfinder", "CDS", s, e, ".", l.strand, "0", attr + `;Parent=${id}.gene`].join("\t") + "\n";
      });
    }
    return o;
  }
  function gbLoc(s, e, st) { return st === "-" ? `complement(${s}..${e})` : `${s}..${e}`; }
  /** GenBank qualifier wrapped to 58 characters; text is broken at spaces when possible (parsers join wrapped text lines with a space). */
  function gbQual(k, v) {
    const t = `/${k}="${String(v).replace(/"/g, "'")}"`; const out = []; let line = t;
    while (line.length > 58) { let cut = 58; if (k !== "translation") { const sp = line.lastIndexOf(" ", 58); if (sp > 20) cut = sp; } out.push(line.slice(0, cut)); line = line.slice(cut).replace(/^ /, ""); }
    out.push(line); return out.map(x => "                     " + x).join("\n") + "\n";
  }
  /** GenBank text of one locus region (for clinker): INA CDS plus neighbouring genes. */
  function genbankRegion(genome, l, q, organism) {
    const name = locusId(genome, l).replace(/[^A-Za-z0-9_.-]/g, "_").slice(0, 60); const nt = q.region.nt.toLowerCase(); const L = nt.length;
    const d = new Date(), mon = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"][d.getMonth()];
    let o = `LOCUS       ${name.padEnd(16)} ${String(L).padStart(11)} bp    DNA     linear   UNK ${String(d.getDate()).padStart(2, "0")}-${mon}-${d.getFullYear()}\n`;
    o += `DEFINITION  ${organism || genome} ${q.region.contig}:${q.region.start}-${q.region.end}, INA locus (${l.structural_class}).\nACCESSION   ${name}\nVERSION     ${name}\nSOURCE      ${organism || genome}\n  ORGANISM  ${organism || genome}\nFEATURES             Location/Qualifiers\n`;
    o += `     source          1..${L}\n` + gbQual("note", `${genome} ${q.region.contig}:${q.region.start}-${q.region.end}`);
    const ia = q.ina_in_region;
    // partial ends: no start codon (5') or ORF running into the contig end without a stop codon (3')
    const p5 = l.start_source === "ORF_start_no_Met", p3 = l.contig_edge === "3prime_open" || l.contig_edge === "both_ends";
    const lo = ia.st === "+" ? (p5 ? "<" : "") + ia.s : (p3 ? "<" : "") + ia.s, hi = ia.st === "+" ? (p3 ? ">" : "") + ia.e : (p5 ? ">" : "") + ia.e;
    o += `     CDS             ${ia.st === "-" ? `complement(${lo}..${hi})` : `${lo}..${hi}`}\n` + gbQual("locus_tag", `INA_${l.locus_index}`) + gbQual("gene", productOf(l).startsWith("ice") ? "ina" : "ina_candidate") + gbQual("product", productOf(l)) +
      gbQual("note", `INpro FINDER: ${l.structural_class}; repeat identity ${l.repeat_identity}; N i-Evalue ${fmtE(l.N_ievalue) || "none"}; C i-Evalue ${fmtE(l.C_ievalue) || "none"}; start ${l.start_source}${l.contig_edge ? "; contig edge " + l.contig_edge : ""}`) + (p5 ? "                     /codon_start=1\n" : "") + gbQual("translation", l.protein);
    q.neighbours.forEach((g, k) => {
      o += `     CDS             ${gbLoc(g.s, g.e, g.st)}\n` + gbQual("locus_tag", g.tag || `ORF_${k + 1}`) + gbQual("product", g.prod || "hypothetical protein") +
        (g.pseudo ? gbQual("note", "annotated as pseudo") : "") + (q.neighbour_source !== "annotation" ? gbQual("note", q.neighbour_source) : "") + (g.translation ? gbQual("translation", g.translation) : "");
    });
    o += "ORIGIN\n";
    for (let i = 0; i < L; i += 60) { let line = String(i + 1).padStart(9) + " "; for (let j = i; j < Math.min(i + 60, L); j += 10) line += nt.slice(j, j + 10) + " "; o += line.trimEnd() + "\n"; }
    return o + "//\n";
  }

  /** One-page PDF with one image (FlateDecode RGB or DCTDecode JPEG). */
  function buildPdf(o) {
    const enc = s => { const b = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) b[i] = s.charCodeAt(i) & 255; return b; };
    const hex16 = t => "<FEFF" + Array.from(String(t || "")).map(ch => { const c = ch.codePointAt(0); const u = c > 0xFFFF ? [0xD800 + ((c - 0x10000) >> 10), 0xDC00 + ((c - 0x10000) & 1023)] : [c]; return u.map(x => x.toString(16).padStart(4, "0")).join(""); }).join("").toUpperCase() + ">";
    const f = v => (+v).toFixed(2).replace(/\.?0+$/, "");
    const content = `q ${f(o.pageW)} 0 0 ${f(o.pageH)} 0 0 cm /Im0 Do Q\n`;
    const parts = [], offs = []; let pos = 0;
    const push = b => { parts.push(b); pos += b.length; };
    const obj = (n, body) => { offs[n] = pos; push(enc(`${n} 0 obj\n${body}\nendobj\n`)); };
    push(new Uint8Array([0x25, 0x50, 0x44, 0x46, 0x2D, 0x31, 0x2E, 0x34, 0x0A, 0x25, 0xE2, 0xE3, 0xCF, 0xD3, 0x0A]));
    obj(1, "<< /Type /Catalog /Pages 2 0 R >>");
    obj(2, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    obj(3, `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${f(o.pageW)} ${f(o.pageH)}] /Resources << /XObject << /Im0 5 0 R >> >> /Contents 4 0 R >>`);
    obj(4, `<< /Length ${content.length} >>\nstream\n${content}endstream`);
    offs[5] = pos; push(enc(`5 0 obj\n<< /Type /XObject /Subtype /Image /Width ${o.imgW} /Height ${o.imgH} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /${o.filter} /Length ${o.data.length} >>\nstream\n`));
    push(o.data); push(enc("\nendstream\nendobj\n"));
    obj(6, `<< /Title ${hex16(o.title)} /Producer (INpro FINDER ${C.VERSION}) >>`);
    const xref = pos; let x = "xref\n0 7\n0000000000 65535 f \n"; for (let n = 1; n <= 6; n++) x += String(offs[n]).padStart(10, "0") + " 00000 n \n";
    push(enc(x + `trailer\n<< /Size 7 /Root 1 0 R /Info 6 0 R >>\nstartxref\n${xref}\n%%EOF\n`));
    const out = new Uint8Array(pos); let k = 0; for (const b of parts) { out.set(b, k); k += b.length; } return out;
  }

  /** Locus architecture as SVG (protein coordinates): ORFs of the locus, terminal domains, repeat blocks, other hits. */
  function architectureSVG(l, opts = {}) {
    const W = opts.width || 900, ML = 70, MR = 20, rowH = 26;
    const orfs = (l.orfs || []).slice().sort((a, b) => (l.strand === "+" ? a.start - b.start : b.end - a.end));
    const g0 = l.strand === "+" ? Math.min(...orfs.map(o => o.start)) : Math.max(...orfs.map(o => o.end));
    const g1 = l.strand === "+" ? Math.max(...orfs.map(o => o.end)) : Math.min(...orfs.map(o => o.start));
    const span = Math.abs(g1 - g0) + 1, sc = (W - ML - MR) / span;
    const gx = p => ML + (l.strand === "+" ? p - g0 : g0 - p) * sc;          // genomic position -> x (INA forward)
    const aaPos = (o, i) => (o.strand === "+" ? o.start + 3 * (i - 1) : o.end - 3 * (i - 1)); // ORF residue -> genomic
    // legend laid out first (wrapped to the figure width) so that the height is known
    const key = [["#2F4682", 1, "N-terminal domain (canonical)"], ["#C86A3C", 1, "C-terminal cap (canonical)"], ["#8A8F9A", 0.6, "weaker profile matches (N blue, C orange)"], ["#C9D3EA", 1, "16-aa repeat region"], ["#E6E8EF", 1, "predicted CDS"]];
    const keyPos = []; { let kx = ML, row = 0; for (const k of key) { const w = 30 + k[2].length * 6.2; if (kx > ML && kx + w > W - MR) { kx = ML; row++; } keyPos.push([kx, row]); kx += w; } }
    const nKeyRows = Math.max(...keyPos.map(p => p[1])) + 1;
    const H = 46 + orfs.length * rowH + 28 + 15 * nKeyRows + 6;
    const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
    let s = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Open Sauce Sans, Helvetica Neue, Arial, sans-serif" font-size="11">`;
    s += `<rect width="${W}" height="${H}" fill="#fff"/>`;
    s += `<text x="${ML}" y="16" font-size="12.5" font-weight="700" fill="#22242A">${esc(opts.title || "")}</text>`;
    s += `<text x="${ML}" y="32" fill="#5B616C">${esc(`${l.contig}:${l.cds_start}-${l.cds_end} (${l.strand}) · ${l.structural_class} · ${l.protein_length_aa} aa · repeat identity ${l.repeat_identity}`)}</text>`;
    orfs.forEach((o, r) => {
      const y = 50 + r * rowH, x0 = gx(o.strand === "+" ? o.start : o.end), x1 = gx(o.strand === "+" ? o.end : o.start);
      const cx0 = gx(o.strand === "+" ? o.cds_start : o.cds_end);
      s += `<text x="${ML - 6}" y="${y + 9}" text-anchor="end" fill="#5B616C" font-size="10">ORF ${r + 1}${o.best ? " *" : ""}</text>`;
      s += `<rect x="${Math.min(x0, x1)}" y="${y + 3}" width="${Math.abs(x1 - x0)}" height="10" fill="#F2F2F7" stroke="#B8B8C0" stroke-width="0.6"/>`;
      s += `<rect x="${Math.min(cx0, x1)}" y="${y + 3}" width="${Math.abs(x1 - cx0)}" height="10" fill="#E6E8EF" stroke="#8A8F9A" stroke-width="0.6"/>`;
      if (o.best && l.repeat_blocks_16aa > 0) { // repeat block region: protein 151 .. len-50
        const p0 = (o.mpos || 0) + 151, p1 = (o.mpos || 0) + l.protein_length_aa - 50; const a = gx(aaPos(o, p0)), b = gx(aaPos(o, p1));
        s += `<rect x="${Math.min(a, b)}" y="${y + 5}" width="${Math.abs(b - a)}" height="6" fill="#C9D3EA"/>`;
      }
      for (const h of o.hits || []) { const a = gx(aaPos(o, h.env_from)), b = gx(aaPos(o, h.env_to)); const col = h.model === "N" ? "#2F4682" : "#C86A3C";
        s += `<rect x="${Math.min(a, b)}" y="${y + (h.canonical ? 0 : 6)}" width="${Math.max(1, Math.abs(b - a))}" height="${h.canonical ? 16 : 4}" fill="${col}" opacity="${h.canonical ? 1 : 0.45}"><title>${esc(`${h.profile} ${h.hmm_from}-${h.hmm_to}, aa ${h.env_from}-${h.env_to}, ${h.score} bits, i-Evalue ${(+h.ievalue).toExponential(1)}${h.canonical ? " (canonical)" : ""}`)}</title></rect>`; }
      s += `<path d="M${x1} ${y + 2} l7 6 l-7 6 z" fill="#8A8F9A"/>`;
    });
    const yl = 50 + orfs.length * rowH + 12;
    key.forEach(([c, op, t], i) => { const [kx, row] = keyPos[i], ky = yl + 15 * row; s += `<rect x="${kx}" y="${ky}" width="12" height="9" fill="${c}" opacity="${op}"/><text x="${kx + 16}" y="${ky + 8}" font-size="10" fill="#22242A">${t}</text>`; });
    s += `<text x="${W - MR}" y="${yl + 15 * nKeyRows + 9}" text-anchor="end" font-size="9" fill="#5B616C">INA drawn left to right, 5' to 3'; ${Math.round(span / 1000 * 10) / 10} kb; * = representative ORF</text>`;
    return s + "</svg>";
  }


  // ------------------------------------------------------------------ minimal ZIP writer (stored, no compression)
  const CRC_T = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
  function crc32(b) { let c = 0xFFFFFFFF; for (let i = 0; i < b.length; i++) c = CRC_T[(c ^ b[i]) & 255] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; }
  function utf8(s) { return typeof TextEncoder !== "undefined" ? new TextEncoder().encode(s) : Uint8Array.from(Buffer.from(s, "utf8")); }
  /** files: [{name, data: string | Uint8Array}] -> Uint8Array (.zip) */
  function zipStore(files) {
    const parts = [], central = []; let off = 0;
    const d = new Date(), dt = ((d.getHours() << 11) | (d.getMinutes() << 5) | (d.getSeconds() >> 1)) & 0xFFFF, dd = (((d.getFullYear() - 1980) << 9) | ((d.getMonth() + 1) << 5) | d.getDate()) & 0xFFFF;
    for (const f of files) {
      const name = utf8(f.name), data = typeof f.data === "string" ? utf8(f.data) : f.data, crc = crc32(data);
      const h = new DataView(new ArrayBuffer(30));
      h.setUint32(0, 0x04034b50, true); h.setUint16(4, 20, true); h.setUint16(6, 0x0800, true); h.setUint16(8, 0, true); h.setUint16(10, dt, true); h.setUint16(12, dd, true);
      h.setUint32(14, crc, true); h.setUint32(18, data.length, true); h.setUint32(22, data.length, true); h.setUint16(26, name.length, true); h.setUint16(28, 0, true);
      parts.push(new Uint8Array(h.buffer), name, data);
      const c = new DataView(new ArrayBuffer(46));
      c.setUint32(0, 0x02014b50, true); c.setUint16(4, 20, true); c.setUint16(6, 20, true); c.setUint16(8, 0x0800, true); c.setUint16(10, 0, true); c.setUint16(12, dt, true); c.setUint16(14, dd, true);
      c.setUint32(16, crc, true); c.setUint32(20, data.length, true); c.setUint32(24, data.length, true); c.setUint16(28, name.length, true); c.setUint32(42, off, true);
      central.push(new Uint8Array(c.buffer), name);
      off += 30 + name.length + data.length;
    }
    const cdSize = central.reduce((a, b) => a + b.length, 0);
    const e = new DataView(new ArrayBuffer(22)); e.setUint32(0, 0x06054b50, true); e.setUint16(8, files.length, true); e.setUint16(10, files.length, true); e.setUint32(12, cdSize, true); e.setUint32(16, off, true);
    const all = parts.concat(central, [new Uint8Array(e.buffer)]); const out = new Uint8Array(all.reduce((a, b) => a + b.length, 0)); let k = 0; for (const b of all) { out.set(b, k); k += b.length; } return out;
  }

  return { productOf, revcomp, translate, simpleGeneCalls, locusSequences, lociTSV, proteinFASTA, geneFASTA, regionFASTA, gff3, genbankRegion, buildPdf, architectureSVG, locusId, TSV_COLS, zipStore, crc32 };
}));
