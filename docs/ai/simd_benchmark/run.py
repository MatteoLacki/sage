import json, os, pathlib, re, subprocess, sys, time

root = pathlib.Path("target/simd-benchmark").resolve()
label, binary = sys.argv[1:]
manifest = json.loads((root / "manifest.json").read_text())
out = root / label
out.mkdir()
argv = manifest["argv"].copy()
argv[0] = str(pathlib.Path(binary).resolve())
argv[argv.index("--output_directory") + 1] = str(out)
(out / "command.json").write_text(json.dumps(argv, indent=2))
env = dict(os.environ, RAYON_NUM_THREADS=str(manifest["threads"]))
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
    "wall_seconds": time.monotonic() - start,
    "search_ms": re.findall(r"- search:\s+(\d+) ms", log),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
if result.returncode:
    print(log[-4000:])
sys.exit(result.returncode)
