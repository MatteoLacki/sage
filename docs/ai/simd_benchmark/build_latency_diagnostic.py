"""Build the latency-instrumented sage binary, then restore clean sources."""
from pathlib import Path
import shutil, subprocess

root = Path("target/simd-benchmark").resolve()
subprocess.run(
    ["python3", "docs/ai/simd_benchmark/instrument_latency.py"], check=True
)
try:
    with (root / "latency-build.log").open("w") as log:
        subprocess.run(
            ["cargo", "build", "--release", "--bin", "sage", "--offline"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    shutil.copy2("target/release/sage", root / "sage-latency")
    print("built", root / "sage-latency")
finally:
    shutil.copyfile(root / "database.unlatency.rs", "crates/sage/src/database.rs")
    shutil.copyfile(root / "runner.unlatency.rs", "crates/sage-cli/src/runner.rs")
    print("sources restored")
