"""Stage 1 - Extract binding constants (Kd) from free text.

Scans a text column (e.g. abstract sentences or curated notes) for
dissociation constants, captures the value and unit, and normalises
everything to nanomolar so entries can be compared.

Usage:
    python pipeline/s01_extract_binding.py --input data/example/entries.csv \
        --output out/01_binding.csv --text-column binding_text
"""
import argparse
import re

import pandas as pd

# Molar unit -> factor to convert into nM
UNIT_TO_NM = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3, "μm": 1e3, "mm": 1e6, "m": 1e9}

KD_PATTERN = re.compile(
    r"(?:\bK\s*[dD]\b|\bKD\b|dissociation constant)"   # the term
    r"(?:\s*\([^)]*\))?"                                # optional "(Kd)"
    r"[^0-9]{0,25}?"                                    # short gap: "of", "=", "~"
    r"(?P<value>\d+(?:\.\d+)?)\s*"
    r"(?P<unit>[pnuµμm]?M)\b",
)


def extract_kd(text):
    """Return (value, unit, value_nM) for the first Kd found, else (None, None, None)."""
    if not isinstance(text, str):
        return None, None, None
    match = KD_PATTERN.search(text)
    if not match:
        return None, None, None
    value = float(match.group("value"))
    unit = match.group("unit")
    factor = UNIT_TO_NM.get(unit.lower())
    if factor is None:
        return value, unit, None
    return value, unit, round(value * factor, 4)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--text-column", default="binding_text")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    parsed = df[args.text_column].apply(extract_kd)
    df["kd_value"] = [p[0] for p in parsed]
    df["kd_unit"] = [p[1] for p in parsed]
    df["kd_nM"] = [p[2] for p in parsed]
    df["kd_status"] = df["kd_nM"].apply(lambda v: "extracted" if pd.notna(v) else "not_found")

    df.to_csv(args.output, index=False)
    found = (df["kd_status"] == "extracted").sum()
    print(f"[stage 1] Kd extracted for {found}/{len(df)} entries -> {args.output}")


if __name__ == "__main__":
    main()
