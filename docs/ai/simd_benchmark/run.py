"""Run one frozen binary against the full F9477 fixture.

Usage: run.py LABEL BINARY [CONFIG]

CONFIG overrides the recorded sage config (used to sweep `bucket_size`).
"""

import hashlib, json, os, pathlib, re, subprocess, sys, time

root = pathlib.Path("target/simd-benchmark").resolve()
label, binary = sys.argv[1], sys.argv[2]
manifest = json.loads((root / "manifest.json").read_text())
out = root / label
out.mkdir()
argv = manifest["argv"].copy()
argv[0] = str(pathlib.Path(binary).resolve())
argv[argv.index("--output_directory") + 1] = str(out)
if len(sys.argv) > 3:
    argv = [sys.argv[3] if a.endswith("sage_config.json") else a for a in argv]
(out / "command.json").write_text(json.dumps(argv, indent=2))
env = dict(os.environ, RAYON_NUM_THREADS=str(manifest["threads"]))
binary_sha256 = hashlib.sha256(pathlib.Path(argv[0]).read_bytes()).hexdigest()
(out / "binary.sha256").write_text(f"{binary_sha256}  {argv[0]}\n")
start = time.monotonic()
with (out / "run.log").open("w") as log:
    result = subprocess.run(
        ["/usr/bin/time", "-v", "-o", str(out / "time.txt"), *argv],
        env=env,
        stdout=log,
        stderr=subprocess.STDOUT,
    )
log = (out / "run.log").read_text()
summary = {
    "label": label,
    "exit_code": result.returncode,
    "binary_sha256": binary_sha256,
    "wall_seconds": time.monotonic() - start,
    "file_io_ms": re.findall(r"- file IO:\s+(\d+) ms", log),
    "search_ms": re.findall(r"- search:\s+(\d+) ms", log),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
if result.returncode:
    print(log[-4000:])
sys.exit(result.returncode)
