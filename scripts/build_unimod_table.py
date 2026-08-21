#!/usr/bin/env python3
"""Fetch unimod.obo and build crates/sage/data/unimod.csv (id,name,mono_mass).

Same source Koina's own IM2Deep_Preprocess_AC/1/modifications.py's Unimod
class downloads and parses (https://www.unimod.org/obo/unimod.obo) -- using
the identical source guarantees SAGE's resolved masses agree with what
Koina itself resolves for the same UNIMOD:i reference. Re-run this script to
refresh the embedded table against a newer Unimod release; the CSV is
checked in and embedded into the SAGE binary via include_str!, not fetched
at runtime.
"""
from __future__ import annotations

import csv
import re
import sys
import urllib.request
from pathlib import Path

UNIMOD_OBO_URL = "https://www.unimod.org/obo/unimod.obo"
OUTPUT_PATH = Path(__file__).resolve().parent.parent / "crates" / "sage" / "data" / "unimod.csv"

ID_RE = re.compile(r"^id: UNIMOD:(\d+)$", re.MULTILINE)
NAME_RE = re.compile(r'^name: (.+)$', re.MULTILINE)
MASS_RE = re.compile(r'^xref: delta_mono_mass "([-\d.]+)"$', re.MULTILINE)


def parse_obo(text: str) -> list[tuple[int, str, float]]:
    term_list = text.split("[Term]")
    term_list.pop(0)  # header/version block, not a term

    rows = []
    for term in term_list:
        id_match = ID_RE.search(term)
        mass_match = MASS_RE.search(term)
        if not id_match or not mass_match:
            continue  # e.g. UNIMOD:0 root node has no delta_mono_mass
        name_match = NAME_RE.search(term)
        name = name_match.group(1).strip() if name_match else ""
        rows.append((int(id_match.group(1)), name, float(mass_match.group(1))))
    return rows


def main() -> None:
    print(f"fetching {UNIMOD_OBO_URL}", file=sys.stderr)
    with urllib.request.urlopen(UNIMOD_OBO_URL) as response:
        text = response.read().decode("utf-8")

    rows = parse_obo(text)
    rows.sort(key=lambda r: r[0])
    print(f"parsed {len(rows)} modifications with a delta_mono_mass", file=sys.stderr)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "name", "mono_mass"])
        writer.writerows(rows)

    print(f"wrote {OUTPUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
