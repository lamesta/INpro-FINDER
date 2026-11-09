#!/usr/bin/env python3
"""Assemble dist/inafinder.html (single self-contained file) from src/ and models/."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
FONTS = ROOT.parent / "genome_map_viewer" / "src" / "fonts.css"   # Open Sauce Sans (patched subset, see genome_map_viewer)

def models_json():
    out = {}
    for k in ["bacterial", "fungal"]:
        spec = json.loads((ROOT / "models" / f"{k}_model.json").read_text())
        out[k] = {"spec": spec, "hmm": {f: (ROOT / "models" / f).read_text() for f in spec["hmm_files"]}}
    return json.dumps(out).replace("</", "<\\/")

def assemble(results_json=None):
    rd = lambda p: (ROOT / p).read_text(encoding="utf-8")
    for p in ["src/inafinder-core.js", "src/inafinder-io.js", "src/ui/worker.js", "src/ui/app.js"]:
        assert "</script" not in rd(p).lower(), p
    fonts = FONTS.read_text(encoding="utf-8") if FONTS.exists() else ""
    res = f'<script type="application/json" id="if-results">{results_json}</script>\n' if results_json else ""
    return ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            "<title>INpro FINDER</title>\n<style>\n" + fonts + "\n" + rd("src/ui/style.css") + "\n</style>\n</head>\n" + rd("src/ui/body.html") + "\n"
            f'<script type="application/json" id="if-models">{models_json()}</script>\n' + res +
            '<script id="src-core">\n' + rd("src/inafinder-core.js") + "\n</script>\n"
            '<script id="src-io">\n' + rd("src/inafinder-io.js") + "\n</script>\n"
            '<script id="src-worker">\n' + rd("src/ui/worker.js") + "\n</script>\n"
            '<script id="src-app">\n' + rd("src/ui/app.js") + "\n</script>\n</body>\n</html>\n")

if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "dist" / "inafinder.html"
    out.write_text(assemble(), encoding="utf-8"); print(f"wrote {out} ({out.stat().st_size/1e6:.2f} MB)")
