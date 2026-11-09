"""Compare JS scan results (validation/js_runs/*.json) with the hmmsearch catalogue (ina_structural_loci.tsv)."""
import csv, json, sys, collections
from pathlib import Path
A = Path("/Users/tommasolamesta/Desktop/PROJECTS/PSAIN/ina_locus_pipeline/analyses/four_genera_complete_ncbi_20260930")
jobs = [l.split("\t")[0] for l in open(sys.argv[1]).read().splitlines() if l]
out = sys.argv[2]
RUNS = sys.argv[3] if len(sys.argv) > 3 else "validation/js_runs"  # exploratory loci (v1.1.0) are ignored
cat = collections.defaultdict(list)
for r in csv.DictReader(open(A / "structural_catalog/ina_structural_loci.tsv"), delimiter="\t"):
    cat[r["genome_id"]].append(r)
rows = []
missing = []
for g in jobs:
    p = Path(f"{RUNS}/{g}.json")
    if not p.exists(): missing.append(g); continue
    js = [l for l in json.load(open(p))["loci"] if not l.get("exploratory")]
    used = set()
    for r in cat.get(g, []):
        a0, a1 = int(r["locus_start"]), int(r["locus_end"])
        best, bo = None, 0
        for i, l in enumerate(js):
            if i in used or l["contig"] != r["contig"] or l["strand"] != r["strand"]: continue
            o = min(a1, l["locus_end"]) - max(a0, l["locus_start"]) + 1
            if o > bo: best, bo = i, o
        row = dict(genome_id=g, ref_contig=r["contig"], ref_strand=r["strand"], ref_locus_start=a0, ref_locus_end=a1, ref_class=r["structural_class"],
                   ref_cds=f'{r["cds_start"]}-{r["cds_end"]}', ref_len=int(r["protein_length_aa"]), ref_rep=float(r["repeat_identity"]),
                   ref_N=r["N_canonical"] == "True", ref_C=r["C_canonical"] == "True", ref_Nev=r["N_ievalue"], ref_Cev=r["C_ievalue"],
                   ref_annot=r["annotation_state"], ref_start_source=r["start_source"])
        if best is None:
            row.update(js_class="MISSING")
        else:
            used.add(best); l = js[best]
            row.update(js_class=l["structural_class"], js_cds=f'{l["cds_start"]}-{l["cds_end"]}', js_len=l["protein_length_aa"], js_rep=l["repeat_identity"],
                       js_N=l["N_canonical"], js_C=l["C_canonical"], js_Nev=l["N_ievalue"], js_Cev=l["C_ievalue"], js_annot=l["annotation_state"], js_start_source=l["start_source"])
        rows.append(row)
    for i, l in enumerate(js):
        if i not in used:
            rows.append(dict(genome_id=g, ref_contig=l["contig"], ref_strand=l["strand"], ref_locus_start=l["locus_start"], ref_locus_end=l["locus_end"], ref_class="ABSENT",
                             js_class=l["structural_class"], js_cds=f'{l["cds_start"]}-{l["cds_end"]}', js_len=l["protein_length_aa"], js_rep=l["repeat_identity"],
                             js_N=l["N_canonical"], js_C=l["C_canonical"], js_Nev=l["N_ievalue"], js_Cev=l["C_ievalue"], js_annot=l["annotation_state"]))
cols = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if rows and k in rows[0] else 99)
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
cm = collections.Counter((r["ref_class"], r["js_class"]) for r in rows)
agree = sum(v for (a, b), v in cm.items() if a == b); tot = sum(cm.values())
print("genomes", len(jobs), "missing runs", len(missing), "loci rows", tot, "class agreement", agree)
for (a, b), v in sorted(cm.items(), key=lambda x: -x[1]):
    if a != b: print("  DIFF", a, "->", b, v)
same_cds = sum(1 for r in rows if r.get("js_cds") and r["ref_class"] not in ("ABSENT",) and r["js_cds"] == r["ref_cds"])
print("identical CDS coords", same_cds, "of", sum(1 for r in rows if r.get("js_cds") and r["ref_class"] != "ABSENT"))
ann = collections.Counter((r.get("ref_annot"), r.get("js_annot")) for r in rows if r["ref_class"] not in ("ABSENT",) and r["js_class"] != "MISSING")
print("annotation state pairs", dict(ann))
