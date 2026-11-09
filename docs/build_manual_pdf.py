"""Build docs/USER_GUIDE.pdf (INpro FINDER) from docs/USER_GUIDE.md.

Usage (from the package root):  python docs/build_manual_pdf.py
Requires: markdown, weasyprint, cairosvg; Node.js for the two example figures (made with the tool's own code).
The Markdown file stays the single source of text; this script adds a title page, a paged table of contents and figures.
"""
import json, re, subprocess
from pathlib import Path

import cairosvg
import markdown
from weasyprint import HTML

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"; IMG = DOCS / "img"; IMG.mkdir(exist_ok=True)
VERSION = "1.1.0"; DATE = "October 2026"
FONTS = ROOT.parent / "genome_map_viewer" / "src" / "fonts.css"
GENOMES = {
    "r1079": ("/Users/tommasolamesta/Desktop/PROJECTS/PSAIN/Genomes/P_syringae_R10.79/assembly.fasta", "P. syringae R10.79 · locus 1"),
    "xfrag": ("/Users/tommasolamesta/Desktop/PROJECTS/PSAIN/ina_locus_pipeline/analyses/four_genera_complete_ncbi_20260930/genomes/Xanthomonas/GCA_030582555.1/genomic.fna.gz",
              "X. fragariae GCA_030582555.1 · locus 1"),
}

# 1. example architecture figures, drawn by the tool itself (bacterial model, no annotation)
JS = r"""
const fs = require("fs"), zlib = require("zlib"); const [, , root, fna, title, out] = process.argv;
global.InaFinderCore = require(root + "/src/inafinder-core.js"); const C = global.InaFinderCore, IO = require(root + "/src/inafinder-io.js");
const rd = f => { const b = fs.readFileSync(f); return (f.endsWith(".gz") ? zlib.gunzipSync(b) : b).toString("utf8"); };
const m = root + "/models/";
const model = C.buildModel("bacterial", fs.readFileSync(m + "INA_N150.hmm", "utf8") + fs.readFileSync(m + "INA_C60.hmm", "utf8"), { N: "INA_N150", C: "INA_C60" });
const r = C.scanGenome({ name: "g", contigs: C.parseFasta(rd(fna)) }, model, null);
fs.writeFileSync(out, IO.architectureSVG(r.loci[0], { title, width: 800 }));
console.log(JSON.stringify({ n: r.loci.length, cls: r.loci[0].structural_class, len: r.loci[0].protein_length_aa }));
"""
(IMG / "_arch.js").write_text(JS)
facts = {}
for key, (fna, title) in GENOMES.items():
    svg = IMG / f"arch_{key}.svg"
    if not (IMG / f"arch_{key}.png").exists():
        res = subprocess.run(["node", str(IMG / "_arch.js"), str(ROOT), fna, title, str(svg)], capture_output=True, text=True, check=True)
        facts[key] = json.loads(res.stdout)
        cairosvg.svg2png(url=str(svg), write_to=str(IMG / f"arch_{key}.png"), output_width=2000)
(IMG / "_arch.js").unlink()

md_text = (DOCS / "USER_GUIDE.md").read_text(encoding="utf-8")
md = markdown.Markdown(extensions=["tables", "fenced_code", "toc", "sane_lists"], extension_configs={"toc": {"toc_depth": "2"}})
body = md.convert(md_text)
body = re.sub(r"<h1[^>]*>.*?</h1>\s*", "", body, count=1, flags=re.S)
body = re.sub(r"<p>Contents:</p>\s*<ol>.*?</ol>\s*", "<!--TOC-->", body, count=1, flags=re.S)
h2 = re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body)
toc = '<nav class="toc"><h2 class="toc-title">Contents</h2><ul>' + "".join(f'<li><a href="#{i}">{t}</a></li>' for i, t in h2) + "</ul></nav>"
body = body.replace("<!--TOC-->", toc)

figs = f"""
<figure class="fig">
  <img src="{(IMG / 'arch_r1079.png').as_uri()}" alt="INA locus of P. syringae R10.79">
  <figcaption><b>Figure 1.</b> Detail view of the INA gene of <i>P. syringae</i> R10.79 (intact_full, 1,342 aa; annotated by Bakta as
  "hypothetical protein", MGGFHP_01089). Dark blue: canonical N-terminal domain (INA_N150); orange: canonical C-terminal cap (INA_C60);
  light bars: weaker matches of single repeat units; pale blue band: the region used for the 16-aa repeat test.</figcaption>
</figure>
<figure class="fig">
  <img src="{(IMG / 'arch_xfrag.png').as_uri()}" alt="Unannotated INA gene of X. fragariae">
  <figcaption><b>Figure 2.</b> An INA gene without an annotated CDS: <i>Xanthomonas fragariae</i> GCA_030582555.1, CP114897.1, 1,525 aa,
  intact_full. The NCBI GenBank file of this assembly has no CDS feature at this position; the six-frame search finds the gene from the
  sequence alone.</figcaption>
</figure>
"""
body = re.sub(r'(<h2 id="7-output-files">)', figs + r"\1", body, count=1)
fig3 = f"""
<figure class="fig small">
  <img src="{(ROOT / 'validation' / 'fig_js_vs_hmmer.png').as_uri()}" alt="Domain scores, JavaScript engine vs HMMER">
  <figcaption><b>Figure 3.</b> Domain scores of the INpro FINDER search against <code>hmmsearch</code> (HMMER 3.4) for 460 matched domains in
  27 genomes (r = 0.998).</figcaption>
</figure>
"""
body = re.sub(r'(<h2 id="11-limits-and-good-practice">)', fig3 + r"\1", body, count=1)
assert body.count("<figure") == 3 and 'class="toc"' in body, "figure anchors not found"

css = (FONTS.read_text(encoding="utf-8") if FONTS.exists() else "") + """
@page { size: A4; margin: 22mm 20mm 20mm 20mm;
  @top-left { content: url(""" + (IMG / "INpro_FINDER_header.png").as_uri() + """); image-resolution: 850dpi; vertical-align: middle; }
  @top-right { content: "INpro FINDER """ + VERSION + """ · User guide"; font: 8pt 'Open Sauce Sans', Helvetica, sans-serif; color: #7a7f8a; vertical-align: middle; }
  @bottom-center { content: "page " counter(page) " of " counter(pages); font: 8pt 'Open Sauce Sans', Helvetica, sans-serif; color: #7a7f8a; } }
@page :first { @top-left { content: none; } @top-right { content: none; } @bottom-center { content: none; } }
html { font-family: 'Open Sauce Sans', 'Helvetica Neue', Arial, sans-serif; font-size: 9.6pt; line-height: 1.45; color: #1f2329; }
.title { page-break-after: always; padding-top: 48mm; }
.title img.logo { width: 118mm; display: block; margin: 0 0 10mm -2mm; }
.title h1 { font-size: 30pt; margin: 0 0 4mm 0; color: #2F4682; }
.title .sub { font-size: 13pt; color: #3a3f48; margin-bottom: 18mm; }
.title .meta { font-size: 9.5pt; color: #5b616c; line-height: 1.7; }
nav.toc { page-break-after: always; }
nav.toc ul { list-style: none; padding: 0; }
nav.toc li { margin: 2.2mm 0; font-size: 10.5pt; }
nav.toc a { color: #1f2329; text-decoration: none; }
nav.toc a::after { content: leader('.') target-counter(attr(href), page); }
h2 { font-size: 14pt; color: #2F4682; margin: 9mm 0 3mm 0; border-bottom: 0.6pt solid #c9d3ea; padding-bottom: 1.5mm; page-break-after: avoid; }
h3 { font-size: 11pt; margin: 5mm 0 2mm 0; page-break-after: avoid; }
p, li { orphans: 3; widows: 3; }
table { border-collapse: collapse; width: 100%; margin: 2.5mm 0 4mm 0; font-size: 8.4pt; page-break-inside: auto; }
th, td { border-bottom: 0.5pt solid #d5d8de; padding: 1.4mm 1.8mm; text-align: left; vertical-align: top; }
th { background: #eef1f7; font-weight: 600; }
tr { page-break-inside: avoid; }
code, pre { font-variant-ligatures: none; font-feature-settings: "liga" 0, "calt" 0; }
code { font-family: 'Courier New', Courier, monospace; font-size: 8.3pt; background: #f3f4f7; padding: 0 0.6mm; border-radius: 1.5pt; }
pre { background: #f3f4f7; padding: 2.5mm 3mm; border-radius: 2pt; font-size: 8pt; line-height: 1.4; white-space: pre-wrap; }
pre code { background: none; padding: 0; }
figure.fig { margin: 5mm 0 6mm 0; page-break-inside: avoid; }
figure.fig img { width: 100%; border: 0.5pt solid #d5d8de; }
figure.fig.small img { width: 62%; display: block; margin: 0 auto; }
figcaption { font-size: 8.4pt; color: #3a3f48; margin-top: 2mm; line-height: 1.4; }
"""
title = f"""<section class="title"><img class="logo" src="{(ROOT / 'branding' / 'INpro_FINDER_logo_transparent.svg').as_uri()}" alt="INpro FINDER">
<div class="sub">Annotation-independent search for ice-nucleation genes in genome sequences<br>User guide, version {VERSION}</div>
<div class="meta">Section for Microbiology, Aarhus University · PSAIN project<br>{DATE}<br>
Browser tool: <code>inafinder.html</code> · Command line: <code>inafinder.py</code></div></section>"""
html = f"<!doctype html><html><head><meta charset='utf-8'><title>INpro FINDER {VERSION} user guide</title><style>{css}</style></head><body>{title}{body}</body></html>"
out = DOCS / "USER_GUIDE.pdf"
HTML(string=html, base_url=str(DOCS)).write_pdf(out)
print(json.dumps({"pdf": str(out), "bytes": out.stat().st_size, "figures": facts}))
