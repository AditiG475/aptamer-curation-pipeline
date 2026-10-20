#!/usr/bin/env python3
"""Decide whether each cleaned sequence is DNA or RNA.

Decision rules (in order):
  1. Only T present (no U)          -> DNA        (evidence: sequence)
  2. Only U present (no T)          -> RNA        (evidence: sequence)
  3. Both T and U present           -> majority base wins, flagged MIXED for review
  4. Neither T nor U (e.g. GGCGCC)  -> use the declared `molecule_type` column if the
                                       input has one, otherwise --default-type; flagged
If a declared type exists and disagrees with the sequence evidence, the row is
flagged `declared_mismatch` so it can be checked against the source paper.

The sequence is then written in its folding alphabet: U for RNA, T for DNA.

Usage:
    detect_molecule_type.py --input cleaned.csv --output typed.csv [--default-type RNA]
"""
import argparse

import pandas as pd


def detect(seq, declared, default_type):
    t_count, u_count = seq.count("T"), seq.count("U")
    declared = str(declared).strip().upper() if isinstance(declared, str) and declared.strip() else None
    if declared not in (None, "DNA", "RNA"):
        declared = None

    if t_count and not u_count:
        mol, evidence, flag = "DNA", f"T={t_count}, U=0", ""
    elif u_count and not t_count:
        mol, evidence, flag = "RNA", f"U={u_count}, T=0", ""
    elif t_count and u_count:
        mol = "DNA" if t_count > u_count else "RNA"
        evidence, flag = f"T={t_count}, U={u_count} (majority)", "mixed_T_and_U"
    elif declared:
        mol, evidence, flag = declared, "no T/U; used declared type", "type_from_metadata"
    else:
        mol, evidence, flag = default_type, f"no T/U and no declared type; default {default_type}", "type_assumed"

    if declared and declared != mol and flag not in ("type_from_metadata",):
        flag = ";".join(filter(None, [flag, f"declared_mismatch({declared})"]))
    return mol, evidence, flag or "ok"


def to_alphabet(seq, molecule_type):
    return seq.replace("U", "T") if molecule_type == "DNA" else seq.replace("T", "U")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--declared-column", default="molecule_type")
    parser.add_argument("--default-type", choices=["DNA", "RNA"], default="RNA")
    args = parser.parse_args()

    df = pd.read_csv(args.input, keep_default_na=False)
    declared = df[args.declared_column] if args.declared_column in df else [None] * len(df)
    results = [detect(s, d, args.default_type) for s, d in zip(df["sequence"], declared)]
    df["detected_type"] = [r[0] for r in results]
    df["type_evidence"] = [r[1] for r in results]
    df["type_flag"] = [r[2] for r in results]
    df["fold_sequence"] = [to_alphabet(s, m) for s, m in zip(df["sequence"], df["detected_type"])]

    df.to_csv(args.output, index=False)
    print(f"[detect] types {df['detected_type'].value_counts().to_dict()} | "
          f"flags {df['type_flag'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
