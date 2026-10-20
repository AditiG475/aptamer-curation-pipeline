#!/usr/bin/env bash
# Run all six stages on the example data (or pass your own CSV as the first argument).
set -euo pipefail

INPUT="${1:-data/example/entries.csv}"
OUT="out"
mkdir -p "$OUT"

python pipeline/s01_extract_binding.py    --input "$INPUT"                 --output "$OUT/01_binding.csv"
python pipeline/s02_doi_lookup.py         --input "$OUT/01_binding.csv"    --output "$OUT/02_doi.csv"
python pipeline/s03_pubchem_annotate.py   --input "$OUT/02_doi.csv"        --output "$OUT/03_pubchem.csv"
python pipeline/s04_classify.py           --input "$OUT/03_pubchem.csv"    --output "$OUT/04_classified.csv"
python pipeline/s05_structure_analysis.py --input "$OUT/04_classified.csv" --output "$OUT/05_structures.csv"
python pipeline/s06_export_excel.py       --input "$OUT/05_structures.csv" --output "$OUT/curated_dataset.xlsx"

echo "Done. Final workbook: $OUT/curated_dataset.xlsx"
