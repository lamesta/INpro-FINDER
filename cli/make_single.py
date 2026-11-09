"""Write dist/inafinder.py: the command-line tool with every source and model file embedded (zlib + base64)."""
import base64, re, zlib, importlib.util
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("inaf", ROOT / "cli" / "inafinder.py"); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
src = (ROOT / "cli" / "inafinder.py").read_text(encoding="utf-8")
blobs = "\n".join(f'    "{k}": """\n{base64.encodebytes(zlib.compress((ROOT / v).read_bytes(), 9)).decode()}""",' for k, v in m._FILES.items())
assert src.count("_EMBED = None") == 1
src = src.replace("_EMBED = None  # filled in the single-file edition (dist/inafinder.py); otherwise files are read from the repository", "_EMBED = {\n" + blobs + "\n}  # single-file edition: sources and models carried inside (zlib + base64)")
out = ROOT / "dist" / "inafinder.py"; out.write_text(src, encoding="utf-8"); print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
