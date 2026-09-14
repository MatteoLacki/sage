import json, sys
from pathlib import Path
import pandas as pd

root = Path("target/simd-benchmark").resolve()
left, right = sys.argv[1:]
keys = ["scannr", "peptide", "label", "charge", "rank", "isotope_error"]
a = (
    pd.read_csv(
        root / left / "results.sage.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    .set_index(keys)
    .sort_index()
)
b = (
    pd.read_csv(
        root / right / "results.sage.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    .set_index(keys)
    .sort_index()
)
assert a.index.is_unique and b.index.is_unique
assert a.index.equals(b.index), "PSM identities differ"
assert a.columns.equals(b.columns)
differences = {c: int((a[c] != b[c]).sum()) for c in a.columns if not a[c].equals(b[c])}
summary = {
    "left": left,
    "right": right,
    "psms": len(a),
    "columns": len(a.columns) + len(keys),
    "differing_columns": differences,
    "exact_except_psm_id": set(differences) <= {"psm_id"},
}
(root / f"compare-{left}-{right}.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
assert summary["exact_except_psm_id"]
