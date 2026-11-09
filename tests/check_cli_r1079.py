"""Check the command-line outputs for R10.79 (+ GFF3): python tests/check_cli_r1079.py OUTDIR"""
import csv, json, sys, zipfile
from pathlib import Path
d = Path(sys.argv[1]); ok = True
def check(cond, msg):
    global ok
    print(("ok   " if cond else "FAIL ") + msg); ok = ok and cond
rows = list(csv.DictReader(open(d / "loci.tsv"), delimiter="\t"))
check(len(rows) == 1, f"one locus ({len(rows)})")
r = rows[0]
check(r["structural_class"] == "intact_full" and r["protein_length_aa"] == "1342", f"intact_full 1342 aa ({r['structural_class']} {r['protein_length_aa']})")
check(r["annotation_state"] == "annotated_intact_CDS" and r["annot_locus_tag"] == "MGGFHP_01089", f"annotated as MGGFHP_01089 ({r['annot_locus_tag']})")
check((r["contig"], r["strand"], r["cds_start"], r["cds_end"]) == ("contig_1", "-", "1271643", "1275671"), f"CDS contig_1:1271643..1275671 (-) incl. stop ({r['cds_start']}..{r['cds_end']})")
for f in ["proteins.faa", "genes.fna", "ina.gff3", "loci_20kb.fna", "run.json", "report.html"]:
    check((d / f).exists() and (d / f).stat().st_size > 0, f"{f} written")
gb = sorted((d / "genbank").glob("*.gbk")); check(len(gb) == 1, f"one GenBank file ({len(gb)})")
try:
    from Bio import SeqIO
    rec = SeqIO.read(gb[0], "genbank"); cds = [f for f in rec.features if f.type == "CDS"]
    check(len(cds) > 5 and len(cds[0].qualifiers["translation"][0]) == 1342, f"GenBank parses with Biopython ({len(cds)} CDS)")
except ImportError:
    print("skip Biopython not installed")
run = json.load(open(d / "run.json")); check(run["version"] == "1.1.0", f"run.json version {run['version']}")
sys.exit(0 if ok else 1)
