#!/usr/bin/env python3
"""Join every stage's output into one table, keyed on entry_id.

Usage:
    merge_results.py --typed typed.csv --secondary secondary.csv --zscore zscore.csv \
        [--tertiary tertiary_summary.csv] [--annotation annotated.csv] --output final.csv
"""
import argparse
import os

import pandas as pd


def read(path):
    return pd.read_csv(path, keep_default_na=False) if path and os.path.exists(path) and os.path.getsize(path) else None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--typed", required=True)
    parser.add_argument("--secondary", required=True)
    parser.add_argument("--zscore", required=True)
    parser.add_argument("--tertiary")
    parser.add_argument("--annotation")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = read(args.typed)
    annotation = read(args.annotation)
    if annotation is not None:
        new_cols = ["entry_id"] + [c for c in annotation.columns if c not in df.columns]
        df = df.merge(annotation[new_cols], on="entry_id", how="left")
    for path in (args.secondary, args.zscore, args.tertiary):
        part = read(path)
        if part is not None:
            part = part.drop(columns=[c for c in ("detected_type", "fold_sequence") if c in part])
            df = df.merge(part, on="entry_id", how="left")
    df.to_csv(args.output, index=False)
    print(f"[merge] {len(df)} entries x {len(df.columns)} columns -> {args.output}")


if __name__ == "__main__":
    main()
