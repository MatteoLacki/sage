import hashlib, json, shlex, tomllib
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

root = Path("target/simd-benchmark").resolve()
root.mkdir(parents=True, exist_ok=True)
node = Path(
    "../../nodes/run_sage/19b746debf1e961d99207116c6132e96c259413eca5d6d1505ce11c90673273d"
).resolve()
record = tomllib.loads((node / ".rip/dependencies.toml").read_text())
tokens = shlex.split(record["command"]["realized"])
start = tokens.index("&&") + 1
argv = tokens[start : tokens.index("&&", start)]
argv.append("--disable-telemetry-i-dont-want-to-improve-sage")
m = {"base_commit": "42510e4", "source_node": str(node), "threads": 16, "argv": argv}
dump = Path(
    "../../nodes/dump_peptides/e3c5b424dfcf90e436277be062ea0500e5b70bd493a2a501660691565f60712c/peptides.parquet"
).resolve()
d = pq.read_table(dump, columns=["peptide", "monoisotopic"]).to_pandas()
keys = pd.Index(d.peptide)
assert keys.is_unique
fingerprint = hashlib.sha256(
    np.asarray(d.monoisotopic, dtype="<f4").tobytes()
).hexdigest()
metadata = {"dumped_peptides_sha256": fingerprint}
rt_path = Path(argv[argv.index("--predicted-rt") + 1])
rt = pq.read_table(rt_path).to_pandas().set_index("sequence")
assert rt.index.is_unique
rt = rt.reindex(keys)
assert rt.rt.notna().all()
pq.write_table(
    pa.table({"rt": rt.rt.to_numpy(dtype=np.float64)}).replace_schema_metadata(
        metadata
    ),
    root / "predicted_rt.parquet",
)
frag_path = Path(argv[argv.index("--predicted-fragment-intensity-index") + 1])
starts = np.full(len(keys) * 3, -1, dtype=np.int64)
ends = starts.copy()
matched = 0
for batch in pq.ParquetFile(frag_path).iter_batches(batch_size=500000):
    b = batch.to_pandas()
    ix = keys.get_indexer(b.sequence)
    assert (ix >= 0).all()
    assert b.charge.between(2, 4).all()
    slots = ix * 3 + b.charge.to_numpy() - 2
    assert (starts[slots] == -1).all()
    starts[slots] = b.start
    ends[slots] = b.end
    matched += len(b)
assert (starts >= 0).sum() == matched
pq.write_table(
    pa.table({"start": starts, "end": ends}).replace_schema_metadata(metadata),
    root / "fragment_index.parquet",
)
for flag, name in [
    ("--predicted-rt", "predicted_rt.parquet"),
    ("--predicted-fragment-intensity-index", "fragment_index.parquet"),
]:
    argv[argv.index(flag) + 1] = str(root / name)
m["prediction_conversion"] = {
    "dump": str(dump),
    "fingerprint": fingerprint,
    "original_rt": str(rt_path),
    "original_fragment_index": str(frag_path),
    "populated_fragment_slots": matched,
}
(root / "manifest.json").write_text(json.dumps(m, indent=2))
print(m["prediction_conversion"], flush=True)
