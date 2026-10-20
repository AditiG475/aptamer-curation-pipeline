#!/usr/bin/env python3
"""Fold-stability Z-score against shuffled sequences.

Z = (MFE_real - mean(MFE_shuffled)) / sd(MFE_shuffled)

A clearly negative Z (e.g. < -2) means the aptamer folds more stably than
random sequences with the same composition. Shuffles keep the exact
nucleotide composition (mononucleotide shuffle) and use the same
molecule-specific energy parameters as the real fold. A fixed random seed
makes results reproducible.

Usage:
    shuffle_zscore.py --input typed.csv --output zscore.csv --shuffles 30 --seed 42
"""
import argparse
import math
import random

import pandas as pd
import RNA

MIN_SHUFFLES = 5
FOLDABLE = {"RNA": set("ACGUN"), "DNA": set("ACGTN")}


def mfe(seq, molecule_type, temperature):
    if molecule_type == "DNA":
        RNA.params_load_DNA_Mathews2004()
    else:
        RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature = temperature
    return float(RNA.fold_compound(seq, md).mfe()[1])


def zscore(seq, molecule_type, n_shuffles, temperature, rng):
    real = mfe(seq, molecule_type, temperature)
    shuffled_mfes = []
    for _ in range(n_shuffles):
        bases = list(seq)
        rng.shuffle(bases)
        try:
            shuffled_mfes.append(mfe("".join(bases), molecule_type, temperature))
        except Exception:
            continue
    n = len(shuffled_mfes)
    if n < MIN_SHUFFLES:
        return {"z_score": None, "shuffle_mean_mfe": None, "shuffle_sd_mfe": None,
                "shuffles_used": n, "zscore_status": "too_few_shuffles"}
    mean = sum(shuffled_mfes) / n
    sd = math.sqrt(sum((x - mean) ** 2 for x in shuffled_mfes) / (n - 1))
    if sd == 0:
        return {"z_score": None, "shuffle_mean_mfe": round(mean, 4), "shuffle_sd_mfe": 0.0,
                "shuffles_used": n, "zscore_status": "zero_variance"}
    return {"z_score": round((real - mean) / sd, 4), "shuffle_mean_mfe": round(mean, 4),
            "shuffle_sd_mfe": round(sd, 4), "shuffles_used": n, "zscore_status": "ok"}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--shuffles", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--temperature", type=float, default=25.0)
    args = parser.parse_args()

    df = pd.read_csv(args.input, keep_default_na=False)
    rng = random.Random(args.seed)
    rows = []
    for _, row in df.iterrows():
        seq, mol = row["fold_sequence"], row["detected_type"]
        if row.get("clean_status", "ok") != "ok" or not seq or set(seq) - FOLDABLE[mol]:
            rows.append({"zscore_status": "skipped"})
            continue
        rows.append(zscore(seq, mol, args.shuffles, args.temperature, rng))
    RNA.params_load_RNA_Turner2004()

    out = pd.concat([df[["entry_id"]], pd.DataFrame(rows)], axis=1)
    out.to_csv(args.output, index=False)
    print(f"[zscore] {out['zscore_status'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
