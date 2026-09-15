import json
import sys
from pathlib import Path

import pandas as pd

root = Path("target/simd-benchmark").resolve()
left, right = sys.argv[1:]
keys = ["scannr", "peptide", "label", "charge", "rank", "isotope_error"]


def read_fragments(label):
    psms = pd.read_csv(
        root / label / "results.sage.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=["psm_id", "peptide_len", *keys],
    ).set_index("psm_id")
    assert psms.index.is_unique
    fragments = pd.read_csv(
        root / label / "matched_fragments.sage.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    ).join(psms, on="psm_id", validate="many_to_one")
    assert fragments[keys].notna().all().all(), "Fragment references unknown PSM"
    ordinal = fragments["fragment_ordinals"].astype(int)
    valid = (ordinal > 0) & (ordinal < fragments["peptide_len"].astype(int))
    assert valid.all(), f"{label}: {int((~valid).sum())} invalid fragment ordinals"
    fragments = fragments.drop(columns=["psm_id", "peptide_len"])
    return fragments.sort_values(list(fragments.columns)).reset_index(drop=True)


a, b = read_fragments(left), read_fragments(right)
summary = {
    "left": left,
    "right": right,
    "left_fragments": len(a),
    "right_fragments": len(b),
    "valid_ordinals": True,
    "exact_except_psm_id_and_row_order": a.equals(b),
}
(root / f"compare-fragments-{left}-{right}.json").write_text(
    json.dumps(summary, indent=2)
)
print(json.dumps(summary, indent=2))
assert summary["exact_except_psm_id_and_row_order"]
