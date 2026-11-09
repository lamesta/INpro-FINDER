"""JS repeatProfile vs the Python repeat_profile of the published pipeline (ina_locus/core.py).
python tests/test_repeat_parity.py proteins.faa [more.faa ...]"""
import json, subprocess, sys
from pathlib import Path
sys.path.insert(0, "/Users/tommasolamesta/Desktop/PROJECTS/PSAIN/ina_locus_pipeline/src")
from ina_locus.core import repeat_profile
HERE = Path(__file__).resolve().parent
seqs = []
for f in sys.argv[1:]:
    name, buf = None, []
    for line in open(f):
        if line.startswith(">"):
            if name: seqs.append((name, "".join(buf).replace("*", "")))
            name, buf = line[1:].split()[0], []
        else: buf.append(line.strip())
    if name: seqs.append((name, "".join(buf).replace("*", "")))
js = subprocess.run(["node", "-e", "const C=require(process.argv[1]);const s=JSON.parse(require('fs').readFileSync(0,'utf8'));"
                     "console.log(JSON.stringify(s.map(x=>C.repeatProfile(x))))", str(HERE.parent / "src/inafinder-core.js")],
                    input=json.dumps([s for _, s in seqs]), capture_output=True, text=True, check=True)
J = json.loads(js.stdout); bad = 0
for (n, s), j in zip(seqs, J):
    p = repeat_profile(s)
    if (p["mean_consensus_identity"], p["repeat_count"], p["repeat_consensus"]) != (j["mean_consensus_identity"], j["repeat_count"], j["repeat_consensus"]):
        bad += 1; print("DIFF", n, p["mean_consensus_identity"], j["mean_consensus_identity"])
print(json.dumps({"proteins": len(seqs), "identical": len(seqs) - bad}))
sys.exit(1 if bad else 0)
