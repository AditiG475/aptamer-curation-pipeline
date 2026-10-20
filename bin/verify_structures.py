#!/usr/bin/env python3
"""Quality checks on the secondary-structure and Z-score results.

Hard checks (a failure means a result is wrong):
  C01 every input entry has a fold row and a Z-score row
  C02 structure length == sequence length
  C03 dot-bracket is balanced
  C04 structure uses only . ( )
  C05 free_bases == number of '.'
  C06 0 <= loop_count <= free_bases
  C07 MFE is finite and <= 0
  C08 base_pairs == number of '('
  C09 Z-score finite where status is ok
  C10 shuffles_used == requested shuffles where status is ok
Warnings (worth a manual look):
  W01 rows that were skipped or failed
  W02 MFE below --mfe-floor (suspiciously low)
  W03 Z-score > 0 (folds less stably than shuffled sequences)
  W04 duplicate sequences

Writes a CSV of every check and prints a summary. Use --strict to exit with
an error when a hard check fails.

Usage:
    verify_structures.py --typed typed.csv --secondary secondary.csv --zscore zscore.csv \
        --shuffles 30 --output qc_checks.csv
"""
import argparse
import math
import sys

import pandas as pd


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def balanced(structure):
    depth = 0
    for char in structure:
        depth += {"(": 1, ")": -1}.get(char, 0)
        if depth < 0:
            return False
    return depth == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--typed", required=True)
    parser.add_argument("--secondary", required=True)
    parser.add_argument("--zscore", required=True)
    parser.add_argument("--shuffles", type=int, default=30)
    parser.add_argument("--mfe-floor", type=float, default=-150.0)
    parser.add_argument("--output", required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    typed = pd.read_csv(args.typed, keep_default_na=False)
    sec = pd.read_csv(args.secondary, keep_default_na=False)
    zs = pd.read_csv(args.zscore, keep_default_na=False)
    folded = sec[sec["fold_status"] == "folded"].copy()
    for column in ("base_pairs", "free_bases", "loop_count"):
        folded[column] = pd.to_numeric(folded[column]).astype(float).astype(int)
    z_ok = zs[zs["zscore_status"] == "ok"]

    def ids(mask, frame):
        return list(frame.loc[mask, "entry_id"])

    checks = [
        ("C01", "every entry has fold and Z-score rows", "error",
         sorted((set(typed.entry_id) - set(sec.entry_id)) | (set(typed.entry_id) - set(zs.entry_id)))),
        ("C02", "structure length equals sequence length", "error",
         ids(folded.structure.str.len() != folded.fold_sequence.str.len(), folded)),
        ("C03", "dot-bracket is balanced", "error", ids(~folded.structure.apply(balanced), folded)),
        ("C04", "structure uses only . ( )", "error",
         ids(folded.structure.apply(lambda s: bool(set(s) - set(".()"))), folded)),
        ("C05", "free_bases equals dot count", "error",
         ids(folded.free_bases != folded.structure.str.count(r"\."), folded)),
        ("C06", "0 <= loop_count <= free_bases", "error",
         ids((folded.loop_count < 0) | (folded.loop_count > folded.free_bases), folded)),
        ("C07", "MFE finite and <= 0", "error",
         ids(~folded.mfe_kcal_mol.apply(lambda v: finite(v) and float(v) <= 0), folded)),
        ("C08", "base_pairs equals '(' count", "error",
         ids(folded.base_pairs != folded.structure.str.count(r"\("), folded)),
        ("C09", "Z-score finite where status ok", "error", ids(~z_ok.z_score.apply(finite), z_ok)),
        ("C10", f"{args.shuffles} shuffles used where status ok", "error",
         ids(pd.to_numeric(z_ok.shuffles_used).astype(float).astype(int) != args.shuffles, z_ok)),
        ("W01", "skipped or failed folds", "warning", ids(sec.fold_status != "folded", sec)),
        ("W02", f"MFE below {args.mfe_floor} kcal/mol", "warning",
         ids(folded.mfe_kcal_mol.astype(float) < args.mfe_floor, folded)),
        ("W03", "Z-score > 0 (less stable than shuffled)", "warning",
         ids(z_ok.z_score.astype(float) > 0, z_ok)),
        ("W04", "duplicate sequences", "warning",
         ids(typed.get("duplicate_of", pd.Series([""] * len(typed))) != "", typed)),
    ]

    report = pd.DataFrame([{
        "check": cid, "description": desc, "level": level,
        "result": "pass" if not bad else ("fail" if level == "error" else "warn"),
        "n_entries": len(bad), "entries": ";".join(map(str, bad[:50])),
    } for cid, desc, level, bad in checks])
    report.to_csv(args.output, index=False)

    for _, r in report.iterrows():
        mark = {"pass": "PASS", "fail": "FAIL", "warn": "WARN"}[r.result]
        extra = f" ({r.n_entries}: {r.entries})" if r.n_entries else ""
        print(f"[verify] {mark} {r.check} {r.description}{extra}")
    failed = (report.result == "fail").sum()
    print(f"[verify] {len(report) - failed - (report.result == 'warn').sum()} passed, "
          f"{failed} failed, {(report.result == 'warn').sum()} warnings")
    if failed and args.strict:
        sys.exit(1)


if __name__ == "__main__":
    main()
