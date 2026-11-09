#!/usr/bin/env python3
"""INpro FINDER: annotation-independent search for ice-nucleation (INA) genes in genome sequences.

  python inafinder.py viewer [-o inafinder.html]
      write the browser tool (one self-contained HTML file; open it and drag in your genomes)

  python inafinder.py scan genome1.fna[.gz] [genome2.fna ...] -o results/ [options]
      --model bacterial | fungal | my_model.json   (default bacterial)
      --annotation GENOME=ANNOTATION                optional GFF3/GenBank per genome (repeatable), or
      --annotation-dir DIR                          annotation files matched to genomes by file-name stem
      --jobs N                                      genomes scanned in parallel (default: CPUs - 1)
      --exploratory                                 also report long ORFs with 16-aa repeats but no profile match
      writes loci.tsv, proteins.faa, genes.fna, ina.gff3, loci_20kb.fna, genbank/*.gbk, run.json and report.html

  python inafinder.py extract DIR
      write the source files and model files carried by this script

Needs Python >= 3.8 (standard library only) and Node.js >= 18 for `scan`. Without Node.js use `viewer`:
the browser tool runs the same code and gives the same results.
"""
import argparse, json, os, re, shutil, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

VERSION = "1.1.0"
_EMBED = None  # filled in the single-file edition (dist/inafinder.py); otherwise files are read from the repository
_REPO = Path(__file__).resolve().parents[1]
_FILES = {
    "src/inafinder-core.js": "src/inafinder-core.js", "src/inafinder-io.js": "src/inafinder-io.js", "src/ui/worker.js": "src/ui/worker.js", "src/ui/app.js": "src/ui/app.js",
    "src/ui/body.html": "src/ui/body.html", "src/ui/style.css": "src/ui/style.css", "src/fonts.css": "../genome_map_viewer/src/fonts.css", "cli/runner.js": "cli/runner.js",
    "models/bacterial_model.json": "models/bacterial_model.json", "models/fungal_model.json": "models/fungal_model.json",
    "models/INA_N150.hmm": "models/INA_N150.hmm", "models/INA_C60.hmm": "models/INA_C60.hmm", "models/INAf_N.hmm": "models/INAf_N.hmm", "models/INAf_Ccap.hmm": "models/INAf_Ccap.hmm",
}


def get(name):
    if _EMBED is not None:
        import base64, zlib
        return zlib.decompress(base64.b64decode(_EMBED[name])).decode("utf-8")
    return (_REPO / _FILES[name]).read_text(encoding="utf-8")


def models_json():
    out = {}
    for k in ["bacterial", "fungal"]:
        spec = json.loads(get(f"models/{k}_model.json"))
        out[k] = {"spec": spec, "hmm": {f: get("models/" + f) for f in spec["hmm_files"]}}
    return out


def assemble_html(results_json=None):
    js = {k: get(k) for k in ["src/inafinder-core.js", "src/inafinder-io.js", "src/ui/worker.js", "src/ui/app.js"]}
    for k, v in js.items():
        assert "</script" not in v.lower(), k
    res = f'<script type="application/json" id="if-results">{results_json.replace("<", chr(92) + "u003c")}</script>\n' if results_json else ""
    return ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            "<title>INpro FINDER</title>\n<style>\n" + get("src/fonts.css") + "\n" + get("src/ui/style.css") + "\n</style>\n</head>\n" + get("src/ui/body.html") + "\n"
            f'<script type="application/json" id="if-models">{json.dumps(models_json()).replace("</", "<" + chr(92) + "/")}</script>\n' + res +
            '<script id="src-core">\n' + js["src/inafinder-core.js"] + "\n</script>\n"
            '<script id="src-io">\n' + js["src/inafinder-io.js"] + "\n</script>\n"
            '<script id="src-worker">\n' + js["src/ui/worker.js"] + "\n</script>\n"
            '<script id="src-app">\n' + js["src/ui/app.js"] + "\n</script>\n</body>\n</html>\n")


def strip_ext(name):
    n = re.sub(r"\.gz$", "", name, flags=re.I)
    return re.sub(r"\.(fa|fna|fasta|fsa|fas|ffn|seq|gff3?|gbff|gbk|gb|genbank)$", "", n, flags=re.I)


def stem(name):
    return re.sub(r"_genomic$", "", strip_ext(name), flags=re.I).lower()


def load_model(arg):
    if arg in ("bacterial", "fungal"):
        spec = json.loads(get(f"models/{arg}_model.json")); hmm = "".join(get("models/" + f) for f in spec["hmm_files"])
    else:
        p = Path(arg); spec = json.loads(p.read_text())
        hmm = "".join((p.parent / f).read_text() for f in spec["hmm_files"])
    return {"name": spec["name"], "hmmText": hmm, "roles": spec["roles"], "rules": spec["rules"]}, spec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("viewer", help="write the browser tool"); v.add_argument("-o", "--out", default="inafinder.html")
    x = sub.add_parser("extract", help="write the carried source and model files"); x.add_argument("dir")
    s = sub.add_parser("scan", help="scan genome FASTA files")
    s.add_argument("genomes", nargs="+"); s.add_argument("-o", "--outdir", required=True)
    s.add_argument("--model", default="bacterial"); s.add_argument("--annotation", action="append", default=[], metavar="GENOME=ANNOTATION")
    s.add_argument("--annotation-dir"); s.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    s.add_argument("--exploratory", action="store_true"); s.add_argument("--flank", type=int, default=20000)
    s.add_argument("--node", default=shutil.which("node") or "node")
    a = ap.parse_args()
    if a.cmd == "viewer":
        Path(a.out).write_text(assemble_html(), encoding="utf-8"); print(f"wrote {a.out} ({os.path.getsize(a.out) / 1e6:.2f} MB)"); return
    if a.cmd == "extract":
        for k in _FILES:
            p = Path(a.dir) / k; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(get(k), encoding="utf-8"); print("wrote", p)
        return
    model, spec = load_model(a.model)
    ann = {}
    for pair in a.annotation:
        g, _, f = pair.partition("=")
        if not f: sys.exit(f"--annotation expects GENOME=ANNOTATION, got {pair}")
        ann[os.path.abspath(g)] = os.path.abspath(f)
    if a.annotation_dir:
        cand = {stem(p.name): str(p.resolve()) for p in Path(a.annotation_dir).iterdir() if re.search(r"\.(gff3?|gbff|gbk|gb)(\.gz)?$", p.name, re.I)}
        for g in a.genomes:
            if os.path.abspath(g) not in ann and stem(os.path.basename(g)) in cand: ann[os.path.abspath(g)] = cand[stem(os.path.basename(g))]
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    params = {"model": a.model if a.model in ("bacterial", "fungal") else "custom", "model_name": model["name"], "roles": model["roles"], "rules": model["rules"],
              "search": {"F1": 0.02, "F2": 1e-3, "F3": 1e-5, "E": 1e-5, "domE": 1e-3, "omega": 1 / 256}, "exploratory": a.exploratory, "started": time.strftime("%Y-%m-%dT%H:%M:%S"), "interface": "command line"}
    with tempfile.TemporaryDirectory() as td:
        for k in ["src/inafinder-core.js", "src/inafinder-io.js", "cli/runner.js"]:
            (Path(td) / Path(k).name).write_text(get(k), encoding="utf-8")
        (Path(td) / "package.json").write_text('{"type": "commonjs"}')  # stops Node searching parent folders
        runner = str(Path(td) / "runner.js")
        jobs = []; used = set()
        for i, g in enumerate(a.genomes):
            name = strip_ext(os.path.basename(g))
            if name.lower() in ("genomic", "assembly", "contigs", "scaffolds", "genome", "sequence"):  # e.g. NCBI datasets: <accession>/genomic.fna
                name = f"{Path(g).resolve().parent.name}_{name}"
            base, k = name, 2
            while name in used: name = f"{base}_{k}"; k += 1
            used.add(name)
            cfg = {"name": name, "fasta": os.path.abspath(g), "annotation": ann.get(os.path.abspath(g)), "model": model,
                                                          "search": params["search"], "exploratory": a.exploratory, "flank": a.flank, "out": str(Path(td) / f"g{i}.json")}
            cp = Path(td) / f"g{i}.cfg.json"; cp.write_text(json.dumps(cfg)); jobs.append((name, cp, cfg["out"]))

        def run_one(j):
            name, cp, res = j
            try: r = subprocess.run([a.node, runner, "scan", str(cp)], capture_output=True, text=True)
            except FileNotFoundError: sys.exit("Node.js not found. Install it (https://nodejs.org or `conda install -c conda-forge nodejs`), or use `python inafinder.py viewer` and the browser tool.")
            line = (r.stdout.strip().splitlines() or ['{"ok": false, "error": "no output"}'])[-1]
            info = json.loads(line) if line.startswith("{") else {"ok": False, "error": r.stderr.strip()[-300:]}
            if not info.get("ok"): Path(res).write_text(json.dumps({"job": name, "error": info.get("error") or r.stderr[-300:]}))
            print(f"  {name}: " + (f"{info['loci']} locus/loci, {info['orfs']} ORFs, {info['seconds']:.1f} s" if info.get("ok") else "ERROR " + str(info.get("error"))), flush=True)
            return res
        print(f"scanning {len(jobs)} genome(s) with model '{model['name']}', {min(a.jobs, len(jobs))} in parallel")
        with ThreadPoolExecutor(max_workers=max(1, min(a.jobs, len(jobs)))) as ex: results = list(ex.map(run_one, jobs))
        fin = Path(td) / "fin.json"; rdata = Path(td) / "report_data.json"
        fin.write_text(json.dumps({"inputs": results, "outdir": str(out.resolve()), "params": params, "report_data": str(rdata)}))
        r = subprocess.run([a.node, runner, "finalize", str(fin)], capture_output=True, text=True)
        info = json.loads(r.stdout.strip().splitlines()[-1]) if r.stdout.strip() else {"ok": False, "error": r.stderr}
        if not info.get("ok"): sys.exit("error: " + str(info.get("error")))
        (out / "report.html").write_text(assemble_html(rdata.read_text(encoding="utf-8")), encoding="utf-8")
    print(f"wrote {out}/: loci.tsv, proteins.faa, genes.fna, ina.gff3, loci_20kb.fna, genbank/ ({info['genbank']} files), run.json, report.html; {info['loci']} locus/loci in {info['genomes']} genome(s)")


if __name__ == "__main__":
    main()
