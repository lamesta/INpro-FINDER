#!/bin/bash
# INpro FINDER test suite. Run from the package root:  bash tests/run_all_tests.sh [--full]
# Needs Node.js >= 18 and Python 3 (Biopython optional). Paths to reference data are those of the PSAIN project folder;
# tests whose inputs are missing are skipped. --full adds the 34-genome regression against stored results (about 1 min of CPU per 8 genomes).
cd "$(dirname "$0")/.."
PS=/Users/tommasolamesta/Desktop/PROJECTS/PSAIN
R=$PS/Genomes/P_syringae_R10.79
ESL=/Users/tommasolamesta/.claude-science/conda/envs/inaevo/bin/esl-translate
OUT=tests/out/suite; rm -rf $OUT; mkdir -p $OUT
pass=0; fail=0; skip=0
t() { local name=$1; shift; if "$@" > $OUT/$name.log 2>&1; then echo "PASS $name"; pass=$((pass+1)); else echo "FAIL $name (see $OUT/$name.log)"; fail=$((fail+1)); fi; }
s() { echo "SKIP $1 ($2)"; skip=$((skip+1)); }

# 1. six-frame ORFs identical to esl-translate (HMMER/Easel)
if [ -x "$ESL" ] && [ -f "$R/assembly.fasta" ]; then
  t orfs_vs_esl_translate bash -c "$ESL -l 50 $R/assembly.fasta > $OUT/esl.faa && node tests/test_orfs.js $R/assembly.fasta $OUT/esl.faa | tee /dev/stderr | grep -q '\"diffSeq\":0,\"onlyMine\":\[\],\"nOnlyMine\":0,\"onlyRef\":\[\]'"
else s orfs_vs_esl_translate "esl-translate or R10.79 not found"; fi

# 2. 16-aa repeat test identical to the Python pipeline
F=/Users/tommasolamesta/Desktop/PROJECTS/INA_EVOLUTION/01_fungal_ina/fungal_ina_proteins.faa
B=$PS/ina_locus_pipeline/analyses/four_genera_complete_ncbi_20260930/hmm/train_all.faa
if [ -f "$F" ] && [ -f "$B" ]; then t repeat_profile_parity python3 tests/test_repeat_parity.py $F $B; else s repeat_profile_parity "reference proteins not found"; fi

# 3. single-file command-line tool on R10.79 with its GFF3
if [ -f "$R/assembly.fasta" ]; then
  t cli_r1079 bash -c "python3 dist/inafinder.py scan $R/assembly.fasta --annotation $R/assembly.fasta=$R/assembly.fasta.gff3 -o $OUT/cli_r1079 --jobs 1 && python3 tests/check_cli_r1079.py $OUT/cli_r1079"
else s cli_r1079 "R10.79 not found"; fi

# 4. the single-file .py carries exactly the browser tool
t single_file_viewer bash -c "python3 dist/inafinder.py viewer -o $OUT/viewer.html && cmp $OUT/viewer.html dist/inafinder.html"

# 5. browser tool in jsdom: load, scan, downloads, saved report reopened
J=../genome_map_viewer/tests/jsdom_env/node_modules/jsdom
if [ -d "$J" ] && [ -f "$R/assembly.fasta" ]; then
  t browser_r1079 bash -c "node tests/ui_test.js dist/inafinder.html $OUT/ui_r1079 $R/assembly.fasta $R/assembly.fasta.gff3 > /dev/null && python3 -c \"
import json,sys; d=json.load(open('$OUT/ui_r1079/ui_result.json'))
l=d['runs'][0]['loci']; ok=not d['errors'] and not d['report']['errors'] and d['selftest']['ok'] and len(l)==1 and l[0][4]=='intact_full' and l[0][5]==1342 and l[0][7]=='MGGFHP_01089' and d['report']['locRows']==1 and all(d['dl_'+k]>0 for k in ['tsv','faa','fna','gff','regions','gbk','json'])
print(json.dumps({'errors':d['errors'],'loci':l,'report_rows':d['report']['locRows']})); sys.exit(0 if ok else 1)\""
else s browser_r1079 "jsdom (genome_map_viewer/tests/jsdom_env) or R10.79 not found"; fi

# 6. optional: 34-genome regression against stored 1.1.0 results (profile-based loci)
if [ "$1" == "--full" ] && [ -f validation/pilot2_jobs.txt ]; then
  t regression_pilot34 bash -c "bash validation/run_batch.sh validation/pilot2_jobs.txt $OUT/pilot bacterial - 8 $OUT/pilot_run.log && python3 -c \"
import json,os,sys
keys=['contig','strand','locus_start','locus_end','structural_class','cds_start','cds_end','protein_length_aa','repeat_identity','N_ievalue','C_ievalue','annotation_state','protein']
bad=0
for f in os.listdir('validation/v110_pilot'):
    a=[l for l in json.load(open('$OUT/pilot/'+f))['loci'] if not l.get('exploratory')]; b=[l for l in json.load(open('validation/v110_pilot/'+f))['loci'] if not l.get('exploratory')]
    bad+= len(a)!=len(b) or any(x.get(k)!=y.get(k) for x,y in zip(a,b) for k in keys)
print('genomes with differences', bad); sys.exit(1 if bad else 0)\""
fi
echo "passed $pass, failed $fail, skipped $skip"
[ $fail -eq 0 ]
