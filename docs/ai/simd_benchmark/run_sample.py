"""Run one frozen binary against the 10,000-precursor sample.

Usage: run_sample.py LABEL BINARY
"""
import hashlib, json, os, pathlib, shutil, subprocess, sys, time

label, binary = sys.argv[1], sys.argv[2]
root = pathlib.Path("target/simd-benchmark").resolve()
m = json.loads((root / "manifest.json").read_text())
argv = m["argv"].copy()
out = root / label
if out.exists():
    shutil.rmtree(out)
out.mkdir()
argv[0] = str(pathlib.Path(binary).resolve())
argv[argv.index("--output_directory") + 1] = str(out)
argv[argv.index("--precursors") + 1] = str(root / "random-10000.parquet")
if len(sys.argv) > 3:
    cfg = sys.argv[3]
    argv = [cfg if a.endswith("sage_config.json") else a for a in argv]

digest = hashlib.sha256(pathlib.Path(argv[0]).read_bytes()).hexdigest()
(out / "command.json").write_text(json.dumps(argv, indent=2))

start = time.monotonic()
with (out / "run.log").open("w") as log:
    subprocess.run(
        argv,
        env=dict(os.environ, RAYON_NUM_THREADS="16"),
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
    )
elapsed = time.monotonic() - start

search = None
for line in (out / "run.log").read_text().splitlines():
    if "- search:" in line:
        search = line.split("- search:")[1].strip()
(out / "summary.json").write_text(
    json.dumps({"label": label, "sha256": digest, "search": search, "wall_s": elapsed}, indent=2)
)
print(f"{label}: sha={digest[:12]} search={search} wall={elapsed:.1f}s")
