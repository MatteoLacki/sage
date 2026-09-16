from pathlib import Path
import json
from pandas_ops.io import read_df

root = Path("target/simd-benchmark").resolve()
m = json.loads((root / "manifest.json").read_text())
p = m["argv"][m["argv"].index("--precursors") + 1]
df = read_df(p)
assert len(df) == 716614, len(df)
df.sample(n=10000, random_state=9477).sort_index().to_parquet(
    root / "random-10000.parquet", index=False
)
print("Uniform sample: 10000 /", len(df))
