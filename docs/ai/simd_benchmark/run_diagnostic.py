import json, os, pathlib, subprocess

root = pathlib.Path("target/simd-benchmark").resolve()
m = json.loads((root / "manifest.json").read_text())
argv = m["argv"].copy()
out = root / "diagnostic-10000"
out.mkdir()
argv[0] = str(root / "sage-instrumented")
argv[argv.index("--output_directory") + 1] = str(out)
argv[argv.index("--precursors") + 1] = str(root / "random-10000.parquet")
(out / "command.json").write_text(json.dumps(argv, indent=2))
with (out / "run.log").open("w") as log:
    subprocess.run(
        argv,
        env=dict(os.environ, RAYON_NUM_THREADS="16"),
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
    )
for line in (out / "run.log").read_text().splitlines():
    if "PAGE_SCAN" in line or "- search:" in line:
        print(line)
