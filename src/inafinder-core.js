/* inafinder-core.js: annotation-independent detection of ice-nucleation (INA) genes in genome sequences.
 *
 * Re-implements, in plain JavaScript (browser and Node.js), the pipeline of
 *   ina_locus_pipeline/analyses/four_genera_complete_ncbi_20260930/scripts/hmm6f_one.sh   (esl-translate -l 50 | hmmsearch)
 *   ina_locus_pipeline/analyses/four_genera_complete_ncbi_20260930/scripts/classify_hmm6f_v2.py
 *   ina_locus_pipeline/src/ina_locus/core.py :: repeat_profile
 *
 * HMMER-compatible parts (HMMER 3.4 / Easel 0.49 conventions):
 *   - profile configuration: local multihit, BLOSUM62 background, insert emissions scored 0 (p7_ProfileConfig);
 *   - acceleration pipeline: generic MSV (P <= F1), generic Viterbi (P <= F2), Forward (P <= F3);
 *   - domains: Viterbi-trace segments, each rescored by unihit Forward/Backward on the segment with null2
 *     bias correction (p7_GNull2_ByExpectation), domain bit score and i-Evalue/c-Evalue as in p7_pipeline.c;
 *   - E-values: Z = number of ORFs searched (as hmmsearch on the esl-translate output), domZ = reported targets.
 * Differences from HMMER: domain boundaries come from the Viterbi trace, not from posterior-decoding
 * domain definition; filters are the generic (non-SIMD) versions without the composition-bias filter.
 * Scores therefore agree closely but not exactly with hmmsearch (see TEST_REPORT).
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.InaFinderCore = factory();
}(typeof self !== "undefined" ? self : this, function () {
  "use strict";
  const VERSION = "1.1.0";
  const LN2 = Math.LN2;
  const NEG = -Infinity;

  // ------------------------------------------------------------------ alphabet
  const AA = "ACDEFGHIKLMNPQRSTVWY";           // HMMER amino order; index 20 = X (any)
  const K = 20, XIDX = 20;
  const AAIDX = new Uint8Array(256).fill(XIDX);
  for (let i = 0; i < K; i++) { AAIDX[AA.charCodeAt(i)] = i; AAIDX[AA.toLowerCase().charCodeAt(i)] = i; }
  // Easel esl_composition_BL62 background, identical to pyhmmer Background(Alphabet.amino())
  const BG = Float64Array.from([0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
    0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364, 0.0540687, 0.0673417, 0.0114135, 0.0304133]);

  // ------------------------------------------------------------------ FASTA
  /** Parse (multi-)FASTA text. Returns [{id, desc, seq}] with seq upper-case and whitespace removed. */
  function parseFasta(text) {
    const out = []; let cur = null, parts = [];
    const flush = () => { if (cur) { cur.seq = parts.join("").replace(/[^A-Za-z*-]/g, "").toUpperCase(); out.push(cur); } parts = []; };
    let pos = 0; const n = text.length;
    while (pos < n) {
      let nl = text.indexOf("\n", pos); if (nl < 0) nl = n;
      let line = text.slice(pos, nl); if (line.endsWith("\r")) line = line.slice(0, -1);
      if (line.startsWith(">")) { flush(); const h = line.slice(1).trim(); const sp = h.search(/\s/); cur = { id: sp < 0 ? h : h.slice(0, sp), desc: sp < 0 ? "" : h.slice(sp + 1).trim(), seq: "" }; }
      else if (cur && line.length) parts.push(line);
      pos = nl + 1;
    }
    flush();
    return out;
  }

  // ------------------------------------------------------------------ translation (NCBI table 1; same codon -> aa map as table 11)
  const NT = { A: 1, C: 2, G: 4, T: 8, U: 8, R: 5, Y: 10, S: 6, W: 9, K: 12, M: 3, B: 14, D: 13, H: 11, V: 7, N: 15 };
  const NTCODE = new Uint8Array(256).fill(15);
  for (const [c, v] of Object.entries(NT)) { NTCODE[c.charCodeAt(0)] = v; NTCODE[c.toLowerCase().charCodeAt(0)] = v; }
  const COMP = new Uint8Array(16);
  for (let b = 0; b < 16; b++) COMP[b] = ((b & 1) ? 8 : 0) | ((b & 2) ? 4 : 0) | ((b & 4) ? 2 : 0) | ((b & 8) ? 1 : 0); // A<->T, C<->G
  const BASES = [1, 2, 4, 8]; // A C G T
  const STD = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"; // codons in A,C,G,T order
  const TRANS = new Uint8Array(4096); // (b1<<8)|(b2<<4)|b3 over 4-bit IUPAC masks -> aa char code; degenerate -> X unless unambiguous
  (function buildTrans() {
    const bi = b => [0, 1, 2, 3].filter(i => b & BASES[i]);
    for (let a = 1; a < 16; a++) for (let b = 1; b < 16; b++) for (let c = 1; c < 16; c++) {
      const set = new Set();
      for (const x of bi(a)) for (const y of bi(b)) for (const z of bi(c)) set.add(STD[x * 16 + y * 4 + z]);
      TRANS[(a << 8) | (b << 4) | c] = (set.size === 1 ? [...set][0] : "X").charCodeAt(0);
    }
  })();
  const STOPC = "*".charCodeAt(0);
  function bigStr(buf) { let s = ""; for (let i = 0; i < buf.length; i += 32768) s += String.fromCharCode.apply(null, buf.slice(i, i + 32768)); return s; }

  /** Six-frame stop-to-stop ORFs of at least minLen aa, as `esl-translate -l minLen` (standard code, no start-codon requirement).
   *  Coordinates are 1-based and exclude the stop codon; start <= end for both strands (`orfCoords` gives the esl-translate form). */
  function sixFrameORFs(contigs, minLen = 50) {
    const orfs = [];
    contigs.forEach((c, ci) => {
      const s = c.seq, L = s.length; if (L < 3) return;
      const fw = new Uint8Array(L); for (let i = 0; i < L; i++) fw[i] = NTCODE[s.charCodeAt(i)];
      for (let f = 0; f < 3; f++) { // frames 1..3
        let buf = [], st = f, first = true;
        for (let i = f; i + 2 < L; i += 3) {
          const aa = TRANS[(fw[i] << 8) | (fw[i + 1] << 4) | fw[i + 2]];
          if (aa === STOPC) { if (buf.length >= minLen) orfs.push({ contig: c.id, ci, frame: f + 1, strand: "+", start: st + 1, end: i, len: buf.length, aa: bigStr(buf), open5: first, open3: false }); buf = []; st = i + 3; first = false; }
          else buf.push(aa);
        }
        if (buf.length >= minLen) orfs.push({ contig: c.id, ci, frame: f + 1, strand: "+", start: st + 1, end: st + 3 * buf.length, len: buf.length, aa: bigStr(buf), open5: first, open3: true });
      }
      for (let f = 0; f < 3; f++) { // frames 4..6, reading the reverse complement from the last nucleotide
        let buf = [], st = L - 1 - f, first = true;
        for (let i = L - 1 - f; i - 2 >= 0; i -= 3) {
          const aa = TRANS[(COMP[fw[i]] << 8) | (COMP[fw[i - 1]] << 4) | COMP[fw[i - 2]]];
          if (aa === STOPC) { if (buf.length >= minLen) orfs.push({ contig: c.id, ci, frame: f + 4, strand: "-", start: i + 2, end: st + 1, len: buf.length, aa: bigStr(buf), open5: first, open3: false }); buf = []; st = i - 3; first = false; }
          else buf.push(aa);
        }
        if (buf.length >= minLen) orfs.push({ contig: c.id, ci, frame: f + 4, strand: "-", start: st - 3 * buf.length + 2, end: st + 1, len: buf.length, aa: bigStr(buf), open5: first, open3: true });
      }
    });
    return orfs;
  }
  function orfCoords(o) { return o.strand === "+" ? `${o.start}..${o.end}` : `${o.end}..${o.start}`; }
  function digitize(aa) { const d = new Uint8Array(aa.length + 1); for (let i = 0; i < aa.length; i++) d[i + 1] = AAIDX[aa.charCodeAt(i)]; return d; } // 1-based

  // ------------------------------------------------------------------ HMMER3 profile files
  /** Parse a HMMER3 ASCII file (one or more models). Probabilities are kept as natural logs (the file stores -ln p). */
  function parseHMM(text) {
    const models = []; const lines = text.split(/\r?\n/); let i = 0;
    const nums = s => s.trim().split(/\s+/).map(v => (v === "*" ? NEG : -parseFloat(v)));
    while (i < lines.length) {
      if (!/^HMMER3/.test(lines[i])) { i++; continue; }
      const m = { name: "", acc: "", desc: "", M: 0, stats: {}, nseq: 0, effn: 0, cksum: "", date: "" };
      const start = i;
      for (; i < lines.length && !/^HMM\s/.test(lines[i]); i++) {
        const l = lines[i], tag = l.slice(0, 6).trim(), val = l.slice(6).trim();
        if (tag === "NAME") m.name = val; else if (tag === "ACC") m.acc = val; else if (tag === "DESC") m.desc = val;
        else if (tag === "LENG") m.M = parseInt(val, 10); else if (tag === "ALPH") m.alph = val;
        else if (tag === "NSEQ") m.nseq = +val; else if (tag === "EFFN") m.effn = +val; else if (tag === "CKSUM") m.cksum = val; else if (tag === "DATE") m.date = val;
        else if (tag === "STATS") { const f = val.split(/\s+/); m.stats[f[1].toLowerCase()] = [parseFloat(f[2]), parseFloat(f[3])]; }
      }
      if (m.alph && m.alph.toLowerCase() !== "amino") throw new Error(`profile ${m.name}: only amino-acid profiles are supported`);
      if (!m.stats.msv || !m.stats.viterbi || !m.stats.forward) throw new Error(`profile ${m.name}: STATS LOCAL lines missing (build the profile with hmmbuild)`);
      i += 2; // "HMM" header line + transition header line
      const M = m.M;
      m.mat = Array.from({ length: M + 1 }, () => new Float64Array(K));
      m.ins = Array.from({ length: M + 1 }, () => new Float64Array(K));
      m.t = Array.from({ length: M + 1 }, () => new Float64Array(7));
      if (/^\s*COMPO/.test(lines[i])) i++;
      m.ins[0].set(nums(lines[i++]).slice(0, K)); m.t[0].set(nums(lines[i++]).slice(0, 7));
      for (let k = 1; k <= M; k++) {
        const f = lines[i++].trim().split(/\s+/);
        if (parseInt(f[0], 10) !== k) throw new Error(`profile ${m.name}: node ${k} not found`);
        for (let x = 0; x < K; x++) m.mat[k][x] = f[1 + x] === "*" ? NEG : -parseFloat(f[1 + x]);
        m.ins[k].set(nums(lines[i++]).slice(0, K)); m.t[k].set(nums(lines[i++]).slice(0, 7));
      }
      while (i < lines.length && lines[i].trim() !== "//") i++;
      m.text = lines.slice(start, i + 1).join("\n") + "\n";
      i++;
      models.push(m);
    }
    if (!models.length) throw new Error("no HMMER3 profile found in file");
    return models;
  }

  /** Configure a search profile (local multihit), following p7_ProfileConfig. Scores in nats. */
  function configureProfile(hmm) {
    const M = hmm.M;
    const P = { name: hmm.name, M, hmm, stats: hmm.stats };
    P.rsc = Array.from({ length: 21 }, () => new Float64Array(M + 1).fill(NEG)); // residue-major match scores
    for (let k = 1; k <= M; k++) {
      let num = 0, den = 0;
      for (let x = 0; x < K; x++) { const sc = hmm.mat[k][x] - Math.log(BG[x]); P.rsc[x][k] = sc; num += sc * BG[x]; den += BG[x]; }
      P.rsc[XIDX][k] = num / den; // esl_abc_FExpectScVec for X
    }
    const T = idx => { const a = new Float64Array(M + 1).fill(NEG); for (let k = 1; k < M; k++) a[k] = hmm.t[k][idx]; return a; };
    P.tMM = T(0); P.tMI = T(1); P.tMD = T(2); P.tIM = T(3); P.tII = T(4); P.tDM = T(5); P.tDD = T(6);
    const tp = k => Array.from(hmm.t[k], Math.exp);
    const occ = new Float64Array(M + 1); const t0 = tp(0); occ[1] = t0[1] + t0[0];
    for (let k = 2; k <= M; k++) { const tk = tp(k - 1); occ[k] = occ[k - 1] * (tk[0] + tk[1]) + (1 - occ[k - 1]) * tk[5]; }
    let Z = 0; for (let k = 1; k <= M; k++) Z += occ[k] * (M - k + 1);
    P.tBM = new Float64Array(M + 1).fill(NEG); for (let k = 1; k <= M; k++) P.tBM[k] = Math.log(occ[k] / Z);
    P.tbmkMSV = Math.log(2 / (M * (M + 1)));
    return P;
  }
  function lengthParams(L, nj) { const pmove = (2 + nj) / (L + 2 + nj); return { loop: Math.log(1 - pmove), move: Math.log(pmove) }; }
  function nullOne(L) { const p1 = L / (L + 1); return L * Math.log(p1) + Math.log(1 - p1); }
  function gumbelSurv(x, mu, lambda) { const y = lambda * (x - mu); return y > 50 ? Math.exp(-y) : -Math.expm1(-Math.exp(-y)); }
  function expLogSurv(x, tau, lambda) { return x < tau ? 0 : -lambda * (x - tau); }

  // ------------------------------------------------------------------ MSV (p7_GMSV)
  function msvScore(P, dsq, L, mm) {
    const M = P.M, loop = Math.log(L / (L + 3)), move = Math.log(3 / (L + 3)), tbmk = P.tbmkMSV, half = -LN2;
    mm.fill(NEG);
    let xN = 0, xB = move, xJ = NEG, xC = NEG;
    for (let i = 1; i <= L; i++) {
      const rs = P.rsc[dsq[i]]; const bm = xB + tbmk; let xE = NEG; let prev = NEG;
      for (let k = 1; k <= M; k++) { const cur = mm[k]; const v = (prev > bm ? prev : bm) + rs[k]; mm[k] = v; if (v > xE) xE = v; prev = cur; }
      const j1 = xJ + loop, j2 = xE + half; xJ = j1 > j2 ? j1 : j2;
      const c1 = xC + loop, c2 = xE + half; xC = c1 > c2 ? c1 : c2;
      xN += loop;
      const b1 = xN + move, b2 = xJ + move; xB = b1 > b2 ? b1 : b2;
    }
    return xC + move;
  }

  // ------------------------------------------------------------------ Viterbi (p7_GViterbi), local multihit, score only
  function viterbiScore(P, dsq, L, W) {
    const M = P.M, lp = lengthParams(L, 1), loop = lp.loop, move = lp.move, half = -LN2;
    const tMM = P.tMM, tMI = P.tMI, tMD = P.tMD, tIM = P.tIM, tII = P.tII, tDM = P.tDM, tDD = P.tDD, tBM = P.tBM;
    let pM = W.a, pI = W.b, pD = W.c, cM = W.d, cI = W.e, cD = W.f;
    pM.fill(NEG); pI.fill(NEG); pD.fill(NEG); cM[0] = cI[0] = cD[0] = NEG;
    let xN = 0, xB = move, xJ = NEG, xC = NEG;
    for (let i = 1; i <= L; i++) {
      const rs = P.rsc[dsq[i]]; let xE = NEG;
      for (let k = 1; k < M; k++) {
        let sc = pM[k - 1] + tMM[k - 1]; let v = pI[k - 1] + tIM[k - 1]; if (v > sc) sc = v;
        v = pD[k - 1] + tDM[k - 1]; if (v > sc) sc = v; v = xB + tBM[k]; if (v > sc) sc = v;
        const m = sc + rs[k]; cM[k] = m; if (m > xE) xE = m;
        let s2 = pM[k] + tMI[k]; v = pI[k] + tII[k]; if (v > s2) s2 = v; cI[k] = s2;
        let s3 = cM[k - 1] + tMD[k - 1]; v = cD[k - 1] + tDD[k - 1]; if (v > s3) s3 = v; cD[k] = s3;
      }
      let sc = pM[M - 1] + tMM[M - 1]; let v = pI[M - 1] + tIM[M - 1]; if (v > sc) sc = v;
      v = pD[M - 1] + tDM[M - 1]; if (v > sc) sc = v; v = xB + tBM[M]; if (v > sc) sc = v;
      cM[M] = sc + rs[M]; cI[M] = NEG;
      let d = cM[M - 1] + tMD[M - 1]; v = cD[M - 1] + tDD[M - 1]; if (v > d) d = v; cD[M] = d;
      if (cM[M] > xE) xE = cM[M]; if (d > xE) xE = d;
      const j1 = xJ + loop, j2 = xE + half; xJ = j1 > j2 ? j1 : j2;
      const c1 = xC + loop, c2 = xE + half; xC = c1 > c2 ? c1 : c2;
      xN += loop;
      const b1 = xN + move, b2 = xJ + move; xB = b1 > b2 ? b1 : b2;
      let t = pM; pM = cM; cM = t; t = pI; pI = cI; cI = t; t = pD; pD = cD; cD = t;
    }
    return xC + move;
  }
  function workspace(M) { const f = () => new Float64Array(M + 1); return { a: f(), b: f(), c: f(), d: f(), e: f(), f: f() }; }

  /** Viterbi with traceback; returns {score, domains: [{i1, i2, k1, k2}]} = the B..E segments of the optimal path. */
  function viterbiTrace(P, dsq, L) {
    const M = P.M, W = M + 1, lp = lengthParams(L, 1), loop = lp.loop, move = lp.move, half = -LN2;
    const { tMM, tMI, tMD, tIM, tII, tDM, tDD, tBM } = P;
    const MX = new Float64Array((L + 1) * W).fill(NEG), IX = new Float64Array((L + 1) * W).fill(NEG), DX = new Float64Array((L + 1) * W).fill(NEG);
    const XN = new Float64Array(L + 1), XB = new Float64Array(L + 1), XE = new Float64Array(L + 1).fill(NEG), XJ = new Float64Array(L + 1).fill(NEG), XC = new Float64Array(L + 1).fill(NEG);
    XN[0] = 0; XB[0] = move;
    for (let i = 1; i <= L; i++) {
      const rs = P.rsc[dsq[i]], o = i * W, p = (i - 1) * W; let xE = NEG;
      for (let k = 1; k <= M; k++) {
        let sc = MX[p + k - 1] + tMM[k - 1]; let v = IX[p + k - 1] + tIM[k - 1]; if (v > sc) sc = v;
        v = DX[p + k - 1] + tDM[k - 1]; if (v > sc) sc = v; v = XB[i - 1] + tBM[k]; if (v > sc) sc = v;
        const m = sc + rs[k]; MX[o + k] = m; if (m > xE) xE = m;
        if (k < M) { let s2 = MX[p + k] + tMI[k]; v = IX[p + k] + tII[k]; if (v > s2) s2 = v; IX[o + k] = s2; }
        let s3 = MX[o + k - 1] + tMD[k - 1]; v = DX[o + k - 1] + tDD[k - 1]; if (v > s3) s3 = v; DX[o + k] = s3;
      }
      if (DX[o + M] > xE) xE = DX[o + M];
      XE[i] = xE;
      const j1 = XJ[i - 1] + loop, j2 = xE + half; XJ[i] = j1 > j2 ? j1 : j2;
      const c1 = XC[i - 1] + loop, c2 = xE + half; XC[i] = c1 > c2 ? c1 : c2;
      XN[i] = XN[i - 1] + loop;
      const b1 = XN[i] + move, b2 = XJ[i] + move; XB[i] = b1 > b2 ? b1 : b2;
    }
    const eq = (a, b) => Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(a));
    const doms = []; let i = L, state = "C", k = 0, cur = null, guard = 0;
    while (guard++ < 10 * (L + 2) * (M + 2)) {
      if (state === "C") { if (i > 0 && eq(XC[i], XC[i - 1] + loop)) { i--; continue; } state = "E"; }
      else if (state === "E") {
        cur = { i1: 0, i2: 0, k1: 0, k2: 0 }; const o = i * W; let found = false;
        for (let kk = 1; kk <= M; kk++) if (eq(XE[i], MX[o + kk])) { state = "M"; k = kk; found = true; break; }
        if (!found && eq(XE[i], DX[o + M])) { state = "D"; k = M; found = true; }
        if (!found) throw new Error("Viterbi traceback failed at E");
      } else if (state === "M") {
        if (!cur.k2) { cur.k2 = k; cur.i2 = i; }
        cur.k1 = k; cur.i1 = i;
        const o = (i - 1) * W, v = MX[i * W + k] - P.rsc[dsq[i]][k];
        if (eq(v, XB[i - 1] + tBM[k])) { state = "B"; i--; }
        else if (k > 1 && eq(v, MX[o + k - 1] + tMM[k - 1])) { k--; i--; }
        else if (k > 1 && eq(v, IX[o + k - 1] + tIM[k - 1])) { state = "I"; k--; i--; }
        else if (k > 1 && eq(v, DX[o + k - 1] + tDM[k - 1])) { state = "D"; k--; i--; }
        else throw new Error("Viterbi traceback failed at M");
      } else if (state === "I") {
        const o = (i - 1) * W, v = IX[i * W + k];
        if (eq(v, MX[o + k] + tMI[k])) { state = "M"; i--; }
        else if (eq(v, IX[o + k] + tII[k])) i--;
        else throw new Error("Viterbi traceback failed at I");
      } else if (state === "D") {
        const o = i * W, v = DX[o + k];
        if (eq(v, MX[o + k - 1] + tMD[k - 1])) { state = "M"; k--; }
        else if (eq(v, DX[o + k - 1] + tDD[k - 1])) k--;
        else throw new Error("Viterbi traceback failed at D");
      } else if (state === "B") { doms.push(cur); cur = null; state = eq(XB[i], XN[i] + move) ? "N" : "J"; }
      else if (state === "J") { if (i > 0 && eq(XJ[i], XJ[i - 1] + loop)) i--; else state = "E"; }
      else break; // N
    }
    return { score: XC[L] + move, domains: doms.reverse() };
  }

  // ------------------------------------------------------------------ Forward / Backward in log space (p7_GForward / p7_GBackward)
  function logsum(a, b) { if (a === NEG) return b; if (b === NEG) return a; return a > b ? a + Math.log1p(Math.exp(b - a)) : b + Math.log1p(Math.exp(a - b)); }
  /** Forward over dsq[from..to] (1-based inclusive); profile length-configured for Lcfg; multihit or unihit (local in both cases). */
  function forward(P, dsq, from, to, Lcfg, multihit, keep) {
    const M = P.M, W = M + 1, L = to - from + 1, lp = lengthParams(Lcfg, multihit ? 1 : 0), loop = lp.loop, move = lp.move;
    const eJ = multihit ? -LN2 : NEG, eC = multihit ? -LN2 : 0;
    const { tMM, tMI, tMD, tIM, tII, tDM, tDD, tBM } = P;
    const MX = new Float64Array((L + 1) * W).fill(NEG), IX = new Float64Array((L + 1) * W).fill(NEG), DX = new Float64Array((L + 1) * W).fill(NEG);
    const XN = new Float64Array(L + 1).fill(NEG), XB = new Float64Array(L + 1).fill(NEG), XE = new Float64Array(L + 1).fill(NEG), XJ = new Float64Array(L + 1).fill(NEG), XC = new Float64Array(L + 1).fill(NEG);
    XN[0] = 0; XB[0] = move;
    for (let i = 1; i <= L; i++) {
      const rs = P.rsc[dsq[from + i - 1]], o = i * W, p = (i - 1) * W; let xE = NEG;
      for (let k = 1; k <= M; k++) {
        const sc = logsum(logsum(MX[p + k - 1] + tMM[k - 1], IX[p + k - 1] + tIM[k - 1]), logsum(XB[i - 1] + tBM[k], DX[p + k - 1] + tDM[k - 1]));
        MX[o + k] = sc + rs[k];
        if (k < M) IX[o + k] = logsum(MX[p + k] + tMI[k], IX[p + k] + tII[k]);
        DX[o + k] = logsum(MX[o + k - 1] + tMD[k - 1], DX[o + k - 1] + tDD[k - 1]);
        xE = logsum(xE, logsum(MX[o + k], DX[o + k])); // local exit from M_k and D_k
      }
      XE[i] = xE;
      XJ[i] = logsum(XJ[i - 1] + loop, xE + eJ);
      XC[i] = logsum(XC[i - 1] + loop, xE + eC);
      XN[i] = XN[i - 1] + loop;
      XB[i] = logsum(XN[i] + move, XJ[i] + move);
    }
    const sc = XC[L] + move;
    return keep ? { sc, MX, IX, DX, XN, XB, XE, XJ, XC, L, loop, move, eJ, eC } : { sc };
  }
  function backward(P, dsq, from, to, F) {
    const M = P.M, W = M + 1, L = F.L, loop = F.loop, move = F.move, eJ = F.eJ, eC = F.eC;
    const { tMM, tMI, tMD, tIM, tII, tDM, tDD, tBM } = P;
    const MX = new Float64Array((L + 1) * W).fill(NEG), IX = new Float64Array((L + 1) * W).fill(NEG), DX = new Float64Array((L + 1) * W).fill(NEG);
    const XN = new Float64Array(L + 1).fill(NEG), XB = new Float64Array(L + 1).fill(NEG), XE = new Float64Array(L + 1).fill(NEG), XJ = new Float64Array(L + 1).fill(NEG), XC = new Float64Array(L + 1).fill(NEG);
    XC[L] = move; XE[L] = XC[L] + eC;
    { const o = L * W; MX[o + M] = XE[L]; DX[o + M] = XE[L];
      for (let k = M - 1; k >= 1; k--) { MX[o + k] = logsum(XE[L], DX[o + k + 1] + tMD[k]); DX[o + k] = logsum(XE[L], DX[o + k + 1] + tDD[k]); } }
    for (let i = L - 1; i >= 0; i--) {
      const rsn = P.rsc[dsq[from + i]], o = i * W, n = (i + 1) * W;
      let xB = NEG; for (let k = 1; k <= M; k++) xB = logsum(xB, MX[n + k] + tBM[k] + rsn[k]); XB[i] = xB;
      if (i === 0) { XN[0] = logsum(XN[1] + loop, XB[0] + move); break; }
      XJ[i] = logsum(XJ[i + 1] + loop, XB[i] + move);
      XC[i] = XC[i + 1] + loop;
      XE[i] = logsum(XJ[i] + eJ, XC[i] + eC);
      XN[i] = logsum(XN[i + 1] + loop, XB[i] + move);
      MX[o + M] = XE[i]; DX[o + M] = XE[i];
      for (let k = M - 1; k >= 1; k--) {
        MX[o + k] = logsum(logsum(MX[n + k + 1] + tMM[k] + rsn[k + 1], IX[n + k] + tMI[k]), logsum(XE[i], DX[o + k + 1] + tMD[k]));
        IX[o + k] = logsum(MX[n + k + 1] + tIM[k] + rsn[k + 1], IX[n + k] + tII[k]);
        DX[o + k] = logsum(MX[n + k + 1] + tDM[k] + rsn[k + 1], logsum(DX[o + k + 1] + tDD[k], XE[i]));
      }
    }
    return { sc: XN[0], MX, IX, DX, XN, XJ, XC };
  }
  /** null2 odds ratios over an envelope (p7_GDecoding + p7_GNull2_ByExpectation). */
  function null2ForEnvelope(P, F, B) {
    const M = P.M, W = M + 1, L = F.L, tot = F.sc;
    const mOcc = new Float64Array(W), iOcc = new Float64Array(W); let xOcc = 0;
    const pm = new Float64Array(W), pi = new Float64Array(W);
    for (let i = 1; i <= L; i++) {
      const o = i * W; let den = 0;
      for (let k = 1; k <= M; k++) { pm[k] = Math.exp(F.MX[o + k] + B.MX[o + k] - tot); den += pm[k]; pi[k] = k < M ? Math.exp(F.IX[o + k] + B.IX[o + k] - tot) : 0; den += pi[k]; }
      const xs = Math.exp(F.XN[i - 1] + B.XN[i] + F.loop - tot) + Math.exp(F.XJ[i - 1] + B.XJ[i] + F.loop - tot) + Math.exp(F.XC[i - 1] + B.XC[i] + F.loop - tot);
      den += xs; const inv = 1 / den;
      for (let k = 1; k <= M; k++) { mOcc[k] += pm[k] * inv; iOcc[k] += pi[k] * inv; }
      xOcc += xs * inv;
    }
    for (let k = 1; k <= M; k++) { mOcc[k] /= L; iOcc[k] /= L; } xOcc /= L;
    const n2 = new Float64Array(21);
    for (let x = 0; x < K; x++) { let s = 0; const r = P.rsc[x]; for (let k = 1; k <= M; k++) s += mOcc[k] * Math.exp(r[k]) + iOcc[k]; n2[x] = s + xOcc; }
    let a = 0; for (let x = 0; x < K; x++) a += n2[x]; n2[XIDX] = a / K; // esl_abc_FAvgScVec for X
    return n2;
  }

  // ------------------------------------------------------------------ hmmsearch-like search of one profile against a set of ORFs
  const SEARCH_DEFAULTS = { F1: 0.02, F2: 1e-3, F3: 1e-5, E: 1e-5, domE: 1e-3, omega: 1 / 256 };
  function scoreTarget(P, dsq, L, o, Z) {
    const st = P.stats; const nullsc = nullOne(L);
    const fw = forward(P, dsq, 1, L, L, true, false).sc;
    if (Math.exp(expLogSurv((fw - nullsc) / LN2, st.forward[0], st.forward[1])) > o.F3) return null;
    const tr = viterbiTrace(P, dsq, L);
    const doms = [];
    for (const d of tr.domains) {
      const Ld = d.i2 - d.i1 + 1; if (Ld < 1) continue;
      const F = forward(P, dsq, d.i1, d.i2, L, false, true);
      const B = backward(P, dsq, d.i1, d.i2, F);
      const n2 = null2ForEnvelope(P, F, B);
      let corr = 0; for (let i = d.i1; i <= d.i2; i++) corr += Math.log(n2[dsq[i]]);
      const dombias = logsum(0, Math.log(o.omega) + corr);
      const bits = (F.sc + (L - Ld) * Math.log(L / (L + 3)) - (nullsc + dombias)) / LN2;
      doms.push({ hmm_from: d.k1, hmm_to: d.k2, ali_from: d.i1, ali_to: d.i2, env_from: d.i1, env_to: d.i2, score: bits, bias: dombias / LN2, lnP: expLogSurv(bits, st.forward[0], st.forward[1]), envsc: F.sc, corr });
    }
    let n2sum = 0; for (const d of doms) n2sum += d.corr;
    let seqScore = (fw - (nullsc + logsum(0, Math.log(o.omega) + n2sum))) / LN2;
    let sum = 0, Ld = 0, sb = 0; for (const d of doms) if (d.envsc - d.corr > 0) { sum += d.envsc; Ld += d.env_to - d.env_from + 1; sb += d.corr; }
    sum += (L - Ld) * Math.log(L / (L + 3));
    const sumScore = (sum - (nullsc + logsum(0, Math.log(o.omega) + sb))) / LN2; if (Ld > 0 && sumScore > seqScore) seqScore = sumScore;
    const seqLnP = expLogSurv(seqScore, st.forward[0], st.forward[1]);
    return { seqScore, seqE: Math.exp(seqLnP) * Z, domains: doms };
  }
  /** Returns {Z, nMSV, nVit, nFwd, hits: [{orf, seqScore, seqE, domains: [{hmm_from, hmm_to, ali_from, ali_to, env_from, env_to, score, bias, ievalue, cevalue}]}]} */
  function searchProfile(P, orfs, opts = {}, progress) {
    const o = Object.assign({}, SEARCH_DEFAULTS, opts); const Z = opts.Z || orfs.length;
    const st = P.stats; const W = workspace(P.M); const mrow = new Float64Array(P.M + 1);
    const hits = []; let nMSV = 0, nVit = 0, nFwd = 0;
    for (let n = 0; n < orfs.length; n++) {
      if (progress && n % 5000 === 0) progress(n / orfs.length);
      const orf = orfs[n]; const L = orf.len; const dsq = orf.dsq || (orf.dsq = digitize(orf.aa));
      const nullsc = nullOne(L);
      if (gumbelSurv((msvScore(P, dsq, L, mrow) - nullsc) / LN2, st.msv[0], st.msv[1]) > o.F1) continue; nMSV++;
      if (gumbelSurv((viterbiScore(P, dsq, L, W) - nullsc) / LN2, st.viterbi[0], st.viterbi[1]) > o.F2) continue; nVit++;
      const r = scoreTarget(P, dsq, L, o, Z); if (!r) continue; nFwd++;
      if (r.seqE > o.E) continue;
      hits.push(Object.assign({ orf: n }, r));
    }
    const domZ = hits.length;
    for (const h of hits) h.domains = h.domains.map(d => Object.assign(d, { ievalue: Math.exp(d.lnP) * Z, cevalue: Math.exp(d.lnP) * domZ })).filter(d => d.cevalue <= o.domE);
    return { Z, nMSV, nVit, nFwd, hits: hits.filter(h => h.domains.length) };
  }

  // ------------------------------------------------------------------ 16-aa repeat periodicity (port of ina_locus.core.repeat_profile)
  function repeatProfile(protein) {
    if (protein.length < 240) return { repeat_period_aa: 16, repeat_count: 0, best_offset: 0, mean_consensus_identity: 0, repeat_consensus: "" };
    const central = protein.slice(150, protein.length - 50);
    let best = [0, 0, 0, ""];
    for (let offset = 0; offset < 16; offset++) {
      const blocks = []; for (let idx = offset; idx < central.length - 15; idx += 16) blocks.push(central.slice(idx, idx + 16));
      if (blocks.length < 5) continue;
      let cons = "";
      for (let p = 0; p < 16; p++) { const cnt = new Map(); for (const b of blocks) { const ch = b[p]; cnt.set(ch, (cnt.get(ch) || 0) + 1); } let bc = "", bn = -1; for (const [ch, c] of cnt) if (c > bn) { bn = c; bc = ch; } cons += bc; }
      let tot = 0; for (const b of blocks) { let s = 0; for (let p = 0; p < 16; p++) if (b[p] === cons[p]) s++; tot += s / 16; }
      const cand = [tot / blocks.length, blocks.length, offset, cons];
      if (cmpTuple(cand, best) > 0) best = cand;
    }
    return { repeat_period_aa: 16, repeat_count: best[1], best_offset: best[2], mean_consensus_identity: pyRound(best[0], 5), repeat_consensus: best[3] };
  }
  /** Mean block identity to the consensus at period P (same method as repeatProfile; unrounded). */
  function repeatIdentityAtPeriod(protein, P) {
    if (protein.length < 240) return 0;
    const central = protein.slice(150, protein.length - 50); let best = 0;
    for (let offset = 0; offset < P; offset++) {
      const blocks = []; for (let idx = offset; idx < central.length - (P - 1); idx += P) blocks.push(central.slice(idx, idx + P));
      if (blocks.length < 5) continue;
      let cons = "";
      for (let q = 0; q < P; q++) { const cnt = new Map(); for (const b of blocks) cnt.set(b[q], (cnt.get(b[q]) || 0) + 1); let bc = "", bn = -1; for (const [ch, c] of cnt) if (c > bn) { bn = c; bc = ch; } cons += bc; }
      let tot = 0; for (const b of blocks) { let sc = 0; for (let q = 0; q < P; q++) if (b[q] === cons[q]) sc++; tot += sc / P; }
      best = Math.max(best, tot / blocks.length);
    }
    return best;
  }
  /** Ratio of 4-aa to 16-aa repeat identity. INA 16-aa repeats give about 0.5 (491 bacterial INA: 0.50 to 0.59; 17 fungal lineage I: 0.50 to 0.57);
   *  simple 4-aa (or 2-aa, 1-aa) repeats give about 1. Returns null when the 16-aa identity is 0. */
  function period4Ratio(protein, id16) {
    if (!(id16 > 0)) return null;
    return Math.round(repeatIdentityAtPeriod(protein, 4) / id16 * 1000) / 1000;
  }
  /** Python round(x, nd): round half to even on the exact binary value. */
  function pyRound(x, nd) {
    if (!isFinite(x) || x === 0) return x;
    const neg = x < 0, s = Math.abs(x).toFixed(100); const dot = s.indexOf(".");
    const ip = s.slice(0, dot), fp = s.slice(dot + 1), keep = fp.slice(0, nd), rest = fp.slice(nd).replace(/0+$/, "");
    let digits = (ip + keep).replace(/^0+(?=\d)/, "");
    let up = false;
    if (rest.length) { if (rest[0] > "5" || (rest[0] === "5" && rest.length > 1)) up = true; else if (rest === "5") up = (+digits[digits.length - 1]) % 2 === 1; }
    let n = BigInt(digits) + (up ? 1n : 0n);
    const v = Number(n) / Math.pow(10, nd);
    return neg ? -v : v;
  }
  function cmpTuple(a, b) { for (let i = 0; i < a.length; i++) { if (a[i] > b[i]) return 1; if (a[i] < b[i]) return -1; } return 0; }

  // ------------------------------------------------------------------ annotation (optional): CDS with pseudo/partial flags
  /** Parse GenBank or GFF3 text into {contigId: [{s, e, st, pseudo, partial, tag, pid, prod}]} (CDS only). */
  function parseAnnotationCDS(text) {
    const out = {}; const add = (c, f) => { (out[c] = out[c] || []).push(f); };
    if (/^LOCUS\s/m.test(text.slice(0, 5000))) {
      for (const r of text.split(/^\/\/\s*$/m)) {
        const vm = r.match(/^VERSION\s+(\S+)/m), lm = r.match(/^LOCUS\s+(\S+)/m); if (!lm) continue;
        const cid = vm ? vm[1] : lm[1];
        const fi = r.search(/^FEATURES/m), oi = r.search(/^ORIGIN/m); if (fi < 0) continue;
        const ft = r.slice(fi, oi < 0 ? r.length : oi).split(/\r?\n/).slice(1);
        let cur = null; const feats = [];
        for (const line of ft) {
          const fm = line.match(/^ {4,6}(\S+)\s+(\S.*)$/);
          if (fm && !line.trim().startsWith("/") && /^ {4,6}\S/.test(line)) { cur = { key: fm[1], loc: fm[2].trim(), q: [] }; feats.push(cur); continue; }
          if (!cur) continue;
          const t = line.trim(); if (!t) continue;
          if (t.startsWith("/")) cur.q.push(t);
          else if (cur.q.length) { const last = cur.q.length - 1; cur.q[last] += (cur.q[last].startsWith("/translation") ? "" : " ") + t; }
          else cur.loc += t;
        }
        for (const f of feats) {
          if (f.key !== "CDS") continue;
          const qv = k => { const x = f.q.find(s => s.startsWith("/" + k + "=") || s === "/" + k); if (x === undefined) return null; return x.slice(k.length + 2).replace(/^"|"$/g, ""); };
          const nums = (f.loc.match(/\d+/g) || []).map(Number); if (!nums.length) continue;
          add(cid, { s: Math.min(...nums), e: Math.max(...nums), st: /complement/.test(f.loc) ? "-" : "+", pseudo: qv("pseudo") !== null || qv("pseudogene") !== null || qv("translation") === null,
            partial: /[<>]/.test(f.loc), tag: qv("locus_tag") || "", pid: qv("protein_id") || "", prod: qv("product") || "" });
        }
      }
    } else {
      const pseudoParents = new Set(); const byId = new Map();
      for (const line of text.split(/\r?\n/)) {
        if (line.startsWith("##FASTA")) break;
        if (!line || line[0] === "#") continue;
        const f = line.split("\t"); if (f.length < 9) continue;
        const a = {}; for (const kv of f[8].split(";")) { const j = kv.indexOf("="); if (j > 0) { const key = kv.slice(0, j).trim(); let v = kv.slice(j + 1).trim(); try { v = decodeURIComponent(v); } catch (e) { } a[key] = v; } }
        if (f[2] === "pseudogene" && a.ID) pseudoParents.add(a.ID);
        if (f[2] !== "CDS") continue;
        const rec = { s: +f[3], e: +f[4], st: f[6] === "-" ? "-" : "+", pseudo: /^true$/i.test(a.pseudo || "") || (a.Parent ? pseudoParents.has(a.Parent) : false),
          partial: /^true$/i.test(a.partial || "") || !!a.start_range || !!a.end_range, tag: a.locus_tag || "", pid: a.protein_id || "", prod: a.product || "" };
        const key = f[0] + "\t" + (a.ID || a.protein_id || a.locus_tag || `${f[3]}-${f[4]}`);
        const y = byId.get(key);
        if (y && y.st === rec.st) { y.s = Math.min(y.s, rec.s); y.e = Math.max(y.e, rec.e); y.pseudo = y.pseudo || rec.pseudo; y.partial = y.partial || rec.partial; }
        else { byId.set(key, rec); add(f[0], rec); }
      }
    }
    return out;
  }

  // ------------------------------------------------------------------ models and classification
  /** Default thresholds = classify_hmm6f_v2.py (bacterial model). */
  const BACTERIAL_RULES = {
    N: { hmm_from_max: 10, hmm_to_min: 100, ievalue_max: 1e-10, env_from_max: 300 },
    C: { hmm_from_max: 15, hmm_to_min: 50, ievalue_max: 1e-10, last_aa: 200 },
    candidate_min_hits: 5, candidate_min_len: 300,
    repeat_identity_min: 0.45, intact_full_min_aa: 800, intact_short_min_aa: 500, repeat_only_min_aa: 300,
    locus_gap_bp: 8000, met_upstream_aa: 60, met_downstream_aa: 30, min_orf_aa: 50,
    // Added in 1.1.0 (defaults keep the 1.0.0 / classify_hmm6f_v2.py behaviour):
    met_fallback: "first_M",      // no Met near the first hit: "first_M" = first Met anywhere in the ORF (catalogue rule); "orf_start" = keep the whole ORF
    max_period4_ratio: 0.75,      // exploratory mode only: reject repeats whose 16-aa identity is explained by a 4-aa (or shorter) period
    exploratory_min_aa: 500,      // exploratory mode only: minimum protein length (300 aa gave 3.5 calls per bacterial genome, 500 aa 0.002; see TEST_REPORT)
    exploratory_identity_min: 0.45, // exploratory mode only: minimum 16-aa repeat identity
  };

  /** Complete detection on one genome.
   *  genome = {name, contigs: [{id, seq}]}; model = buildModel(...); annotation = parseAnnotationCDS(...) or null. */
  function scanGenome(genome, model, annotation, opts = {}) {
    const R = model.rules; const t0 = Date.now(); const log = opts.progress || (() => { });
    const orfs = sixFrameORFs(genome.contigs, R.min_orf_aa);
    log({ stage: "orfs", n: orfs.length });
    const search = {}; const hitsByOrf = new Map();
    for (const pr of model.plist) {
      const res = searchProfile(pr.P, orfs, opts.search || {}, f => log({ stage: "search", profile: pr.P.name, f }));
      search[pr.P.name] = { role: pr.role, Z: res.Z, nMSV: res.nMSV, nVit: res.nVit, nFwd: res.nFwd, nHits: res.hits.length };
      for (const h of res.hits) {
        const a = hitsByOrf.get(h.orf) || []; const L = orfs[h.orf].len, ru = pr.rule;
        for (const d of h.domains) {
          const canon = d.hmm_from <= ru.hmm_from_max && d.hmm_to >= ru.hmm_to_min && d.ievalue <= ru.ievalue_max &&
            (pr.role === "N" ? d.env_from <= ru.env_from_max : d.env_to >= L - ru.last_aa);
          a.push(Object.assign({ model: pr.role, profile: pr.P.name, canonical: canon }, d));
        }
        hitsByOrf.set(h.orf, a);
      }
    }
    const cand = [];
    for (const [n, hs] of hitsByOrf) {
      const orf = orfs[n];
      hs.sort((a, b) => (a.model === b.model ? a.env_from - b.env_from : a.model === "N" ? -1 : 1));
      const N = hs.filter(h => h.model === "N" && h.canonical);
      const C = hs.filter(h => h.model === "C" && h.canonical);
      const od = { orf: n, contig: orf.contig, start: orf.start, end: orf.end, strand: orf.strand, frame: orf.frame, orf_len: orf.len, aa: orf.aa,
        N: N.length ? N.reduce((a, b) => (b.ievalue < a.ievalue ? b : a)) : null, C: C.length ? C.reduce((a, b) => (b.ievalue < a.ievalue ? b : a)) : null,
        n_domain_hits: hs.length, first_hit: Math.min(...hs.map(h => h.env_from)), hits: hs };
      if (od.N || od.C || (hs.length >= R.candidate_min_hits && orf.len >= R.candidate_min_len)) cand.push(od);
    }
    for (const c of cand) {
      const aa = c.aa; let mpos, best_f = null;
      const fs = annotation ? (annotation[c.contig] || []).filter(f => f.st === c.strand && !(f.e - f.s > 50000)) : [];
      const same = fs.filter(f => !f.pseudo && !f.partial && ((c.strand === "+" && Math.abs(f.e - c.end) <= 3) || (c.strand === "-" && Math.abs(f.s - c.start) <= 3)));
      if (annotation && same.length) {
        const f = same[0]; mpos = Math.max(0, c.strand === "+" ? Math.floor((f.s - c.start) / 3) : Math.floor((c.end - f.e) / 3));
        c.start_source = "annotation_same_stop"; c.annotation_state = "annotated_intact_CDS"; best_f = f;
      } else {
        const x = (c.N ? c.N.env_from : c.first_hit) - 1;
        mpos = aa.indexOf("M", Math.max(0, x - R.met_upstream_aa));
        if (mpos < 0 || mpos > x + R.met_downstream_aa) mpos = R.met_fallback === "orf_start" ? 0 : aa.indexOf("M");
        mpos = Math.max(0, mpos); c.start_source = aa[mpos] === "M" ? "Met_rule" : "ORF_start_no_Met";
        if (annotation) {
          const s1 = c.strand === "+" ? c.start + 3 * mpos : c.start, e1 = c.strand === "+" ? c.end : c.end - 3 * mpos;
          const ov = fs.map(f => [Math.min(e1, f.e) - Math.max(s1, f.s) + 1, f]).filter(([o, f]) => o > 0.3 * (f.e - f.s + 1) || o > 0.3 * (e1 - s1 + 1));
          const pse = ov.filter(([, f]) => f.pseudo || f.partial).map(x => x[1]), intact = ov.filter(([, f]) => !(f.pseudo || f.partial)).map(x => x[1]);
          if (pse.length) { c.annotation_state = "annotated_pseudo_or_partial"; best_f = pse[0]; }
          else if (intact.length) { c.annotation_state = "annotated_CDS_different_start_or_stop"; best_f = intact[0]; }
          else c.annotation_state = "no_CDS_annotated";
        } else c.annotation_state = "not_checked";
      }
      c.protein = aa.slice(mpos); c.mpos = mpos; c.annot = best_f;
      const p = repeatProfile(c.protein); c.rep_id = p.mean_consensus_identity; c.rep_blocks = p.repeat_count; c.rep_cons = p.repeat_consensus;
      c.p4 = period4Ratio(c.protein, c.rep_id); c.open5 = !!orfs[c.orf].open5; c.open3 = !!orfs[c.orf].open3;
      if (c.strand === "+") { c.cds_start = c.start + 3 * mpos; c.cds_end = c.end; } else { c.cds_start = c.start; c.cds_end = c.end - 3 * mpos; }
    }
    cand.sort((a, b) => (a.contig < b.contig ? -1 : a.contig > b.contig ? 1 : 0) || (a.strand < b.strand ? -1 : a.strand > b.strand ? 1 : 0) || a.start - b.start);
    const groups = [];
    for (const c of cand) { const g = groups[groups.length - 1]; if (g && g[g.length - 1].contig === c.contig && g[g.length - 1].strand === c.strand && c.start - g[g.length - 1].end <= R.locus_gap_bp) g.push(c); else groups.push([c]); }
    const loci = groups.map((grp, li) => {
      let best = grp[0]; const key = c => [(c.N ? 1 : 0) + (c.C ? 1 : 0), c.protein.length];
      for (const c of grp.slice(1)) if (cmpTuple(key(c), key(best)) > 0) best = c;
      const an = best.annot || {};
      return { locus_index: li + 1, contig: best.contig, strand: best.strand, locus_start: Math.min(...grp.map(c => c.start)), locus_end: Math.max(...grp.map(c => c.end)), n_orfs: grp.length,
        structural_class: classify(best, grp, R), cds_start: best.cds_start, cds_end: best.cds_end, protein_length_aa: best.protein.length, repeat_identity: best.rep_id, repeat_blocks_16aa: best.rep_blocks,
        repeat_consensus: best.rep_cons, N_canonical: !!best.N, C_canonical: !!best.C, N_profile: best.N ? best.N.profile : "", C_profile: best.C ? best.C.profile : "", N_ievalue: best.N ? best.N.ievalue : null, C_ievalue: best.C ? best.C.ievalue : null,
        period4_ratio: best.p4, contig_edge: edgeLabel(best),
        start_source: best.start_source, annotation_state: best.annotation_state, annot_locus_tag: an.tag || "", annot_protein_id: an.pid || "", annot_product: an.prod || "",
        protein: best.protein, mpos: best.mpos,
        orfs: grp.map(c => ({ start: c.start, end: c.end, strand: c.strand, frame: c.frame, len: c.orf_len, mpos: c.mpos, cds_start: c.cds_start, cds_end: c.cds_end, N: c.N ? pick(c.N) : null, C: c.C ? pick(c.C) : null, hits: c.hits.map(pick), best: c === best })) };
    });
    // Exploratory mode (off by default): long ORFs without any profile-based candidate whose 16-aa repeat identity passes the
    // INA threshold and is not explained by a 4-aa (or shorter) period. Low confidence; not part of the validated catalogue rules.
    // Repetitive DNA is periodic in every frame, so ORFs overlapping a profile-based locus (either strand) are skipped, and among
    // overlapping exploratory ORFs (any frame or strand) only the longest protein is kept.
    let nExplored = 0;
    if (opts.exploratory || R.exploratory) {
      const used = new Set(cand.map(c => c.orf));
      const spans = {}; for (const l of loci) (spans[l.contig] = spans[l.contig] || []).push([l.locus_start, l.locus_end]);
      const xc = [];
      for (let n = 0; n < orfs.length; n++) {
        const o = orfs[n]; if (used.has(n) || o.len < R.exploratory_min_aa) continue;
        if ((spans[o.contig] || []).some(([a, b]) => o.start <= b && o.end >= a)) continue;
        nExplored++;
        let mpos = o.aa.indexOf("M"); if (mpos < 0 || mpos > R.met_upstream_aa) mpos = 0;
        const prot = o.aa.slice(mpos); if (prot.length < R.exploratory_min_aa) continue;
        const p = repeatProfile(prot); if (p.mean_consensus_identity < R.exploratory_identity_min) continue;
        const p4 = period4Ratio(prot, p.mean_consensus_identity); if (p4 === null || p4 > R.max_period4_ratio) continue;
        xc.push({ n, mpos, prot, p, p4 });
      }
      xc.sort((x, y) => y.prot.length - x.prot.length || x.n - y.n);
      const kept = [];
      for (const x of xc) { const o = orfs[x.n]; if (!kept.some(k => { const q = orfs[k.n]; return q.contig === o.contig && o.start <= q.end && o.end >= q.start; })) kept.push(x); }
      kept.sort((x, y) => x.n - y.n);
      for (const { n, mpos, prot, p, p4 } of kept) {
        const o = orfs[n];
        const cds_start = o.strand === "+" ? o.start + 3 * mpos : o.start, cds_end = o.strand === "+" ? o.end : o.end - 3 * mpos;
        let state = "not_checked", f0 = null;
        if (annotation) {
          const fs = (annotation[o.contig] || []).filter(f => f.st === o.strand && !(f.e - f.s > 50000));
          const ov = fs.map(f => [Math.min(cds_end, f.e) - Math.max(cds_start, f.s) + 1, f]).filter(([v, f]) => v > 0.3 * (f.e - f.s + 1) || v > 0.3 * (cds_end - cds_start + 1));
          const pse = ov.filter(([, f]) => f.pseudo || f.partial).map(x => x[1]), intact = ov.filter(([, f]) => !(f.pseudo || f.partial)).map(x => x[1]);
          const same = intact.filter(f => (o.strand === "+" && Math.abs(f.e - o.end) <= 3) || (o.strand === "-" && Math.abs(f.s - o.start) <= 3));
          if (same.length) { state = "annotated_CDS_same_stop"; f0 = same[0]; } else if (pse.length) { state = "annotated_pseudo_or_partial"; f0 = pse[0]; }
          else if (intact.length) { state = "annotated_CDS_different_start_or_stop"; f0 = intact[0]; } else state = "no_CDS_annotated";
        }
        const an = f0 || {}, c = { open5: !!o.open5, open3: !!o.open3 };
        loci.push({ locus_index: loci.length + 1, contig: o.contig, strand: o.strand, locus_start: o.start, locus_end: o.end, n_orfs: 1,
          structural_class: "repeat_only_no_profile_hit", cds_start, cds_end, protein_length_aa: prot.length, repeat_identity: p.mean_consensus_identity, repeat_blocks_16aa: p.repeat_count,
          repeat_consensus: p.repeat_consensus, N_canonical: false, C_canonical: false, N_profile: "", C_profile: "", N_ievalue: null, C_ievalue: null,
          period4_ratio: p4, contig_edge: edgeLabel(c), start_source: o.aa[mpos] === "M" ? "first_Met_in_ORF" : "ORF_start_no_Met", annotation_state: state,
          annot_locus_tag: an.tag || "", annot_protein_id: an.pid || "", annot_product: an.prod || "", protein: prot, mpos, exploratory: true,
          orfs: [{ start: o.start, end: o.end, strand: o.strand, frame: o.frame, len: o.len, mpos, cds_start, cds_end, N: null, C: null, hits: [], best: true }] });
      }
    }
    return { genome: genome.name, nContigs: genome.contigs.length, genomeLength: genome.contigs.reduce((a, c) => a + c.seq.length, 0), nORF: orfs.length, search, loci, exploratory: !!(opts.exploratory || R.exploratory), nExplored,
      seconds: (Date.now() - t0) / 1000, model: model.name, model_spec: model.spec, rules: R, version: VERSION };
  }
  function edgeLabel(c) { return c.open5 && c.open3 ? "both_ends" : c.open5 ? "5prime_open" : c.open3 ? "3prime_open" : ""; }
  function pick(h) { return { model: h.model, profile: h.profile, canonical: !!h.canonical, hmm_from: h.hmm_from, hmm_to: h.hmm_to, env_from: h.env_from, env_to: h.env_to, score: +h.score.toFixed(2), bias: +h.bias.toFixed(2), ievalue: h.ievalue }; }

  /** Structural class (cls_of in classify_hmm6f_v2.py). */
  function classify(best, grp, R) {
    const L = best.protein.length, rep = best.rep_id;
    const hasN = grp.some(c => c.N), hasC = grp.some(c => c.C);
    if (best.N && best.C) {
      if (rep >= R.repeat_identity_min && L >= R.intact_full_min_aa) return "intact_full";
      if (rep >= R.repeat_identity_min && L >= R.intact_short_min_aa) return "intact_short";
      return "N_and_C_low_repeat_or_short";
    }
    const gpos = (c, h) => (c.strand === "+" ? c.start + 3 * (h.env_from - 1) : c.end - 3 * (h.env_from - 1));
    const interrupted = () => {
      for (const cN of grp.filter(c => c.N)) for (const cC of grp.filter(c => c.C && c !== cN)) {
        const n = gpos(cN, cN.N), cc = gpos(cC, cC.C);
        if (cC.strand === "+") { if (n < cC.cds_start && cc > cN.cds_end) return true; }
        else if (n > cC.cds_end && cc < cN.cds_start) return true;
      }
      return false;
    };
    if (hasN && hasC && interrupted()) return "interrupted_N_and_C_in_separate_ORFs";
    if (best.C && rep >= R.repeat_identity_min && L >= R.intact_short_min_aa) return "intact_noN_canonical_domain";
    if (best.N && rep >= R.repeat_identity_min && L >= R.intact_short_min_aa) return "intact_noC_canonical_domain";
    if (hasN) return "fragment_N_only";
    if (hasC) return "fragment_C_only";
    return rep >= R.repeat_identity_min && L >= R.repeat_only_min_aa ? "repeat_only" : "weak_repeat_hit";
  }

  /** Model = profiles with a role (N-terminal domain or C-terminal cap) and a per-profile rule, plus general rules.
   *  `roles` is either {N: "profileName", C: "profileName"} (rules.N / rules.C apply) or a list
   *  [{role: "N"|"C", profile: "name", rule: {hmm_from_max, hmm_to_min, ievalue_max, env_from_max | last_aa}}].
   *  hmmText may hold several profiles (concatenated HMMER3 files). */
  function buildModel(name, hmmText, roles, rules, meta) {
    const R = Object.assign(JSON.parse(JSON.stringify(BACTERIAL_RULES)), JSON.parse(JSON.stringify(rules || {})));
    const hs = parseHMM(hmmText); const byName = Object.fromEntries(hs.map(h => [h.name, h]));
    const list = Array.isArray(roles) ? roles : Object.entries(roles).filter(([, nm]) => nm).map(([role, nm]) => ({ role, profile: nm, rule: R[role] }));
    const plist = list.map(x => {
      if (x.role !== "N" && x.role !== "C") throw new Error(`model ${name}: role must be N or C (got ${x.role})`);
      if (!byName[x.profile]) throw new Error(`model ${name}: profile ${x.profile} not found in the HMM file(s)`);
      const rule = Object.assign({}, x.role === "N" ? BACTERIAL_RULES.N : BACTERIAL_RULES.C, x.rule || {});
      return { role: x.role, P: configureProfile(byName[x.profile]), rule };
    });
    if (!plist.length) throw new Error(`model ${name}: no profiles`);
    const profiles = {}; for (const p of plist) if (!profiles[p.role]) profiles[p.role] = p.P;
    return { name, plist, profiles, rules: R, hmmText, meta: meta || {}, spec: list.map(x => ({ role: x.role, profile: x.profile, rule: plist.find(p => p.P.name === x.profile).rule })) };
  }

  return { VERSION, parseFasta, sixFrameORFs, orfCoords, digitize, parseHMM, configureProfile, msvScore, viterbiScore, viterbiTrace, forward, backward, null2ForEnvelope,
    scoreTarget, searchProfile, repeatProfile, repeatIdentityAtPeriod, period4Ratio, parseAnnotationCDS, scanGenome, classify, buildModel, BACTERIAL_RULES, SEARCH_DEFAULTS, BG, workspace, nullOne, gumbelSurv, expLogSurv };
}));
