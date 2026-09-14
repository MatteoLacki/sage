from pathlib import Path
import shutil, subprocess

root = Path("target/simd-benchmark").resolve()
subprocess.run(["python3", "docs/ai/simd_benchmark/instrument.py"], check=True)
try:
    with (root / "instrumented-build.log").open("w") as log:
        subprocess.run(
            ["cargo", "build", "--release", "--bin", "sage", "--offline"],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    shutil.copy2("target/release/sage", root / "sage-instrumented")
finally:
    shutil.copyfile(root / "database.uninstrumented.rs", "crates/sage/src/database.rs")
    shutil.copyfile(root / "runner.uninstrumented.rs", "crates/sage-cli/src/runner.rs")
    shutil.copy2(root / "sage-simd", "target/release/sage")
