#!/usr/bin/env python3
"""Secondary structure prediction with the ViennaRNA Python API.

Uses molecule-specific energy parameters (Turner 2004 for RNA, Mathews 2004
for DNA) instead of converting every sequence to RNA. All values come
straight from the library; nothing is parsed from text output.

Per sequence:
  structure (dot-bracket), MFE (kcal/mol), base pairs, GC %,
  FreeBases (unpaired nt), LoopCount (runs of unpaired nt, incl. 5'/3' tails),
  TerminalUnpaired (unpaired nt in the 5'/3' tails), ensemble diversity

Usage:
    fold_secondary.py --input typed.csv --output secondary.csv --temperature 25
"""
import argparse
import re

import pandas as pd
import RNA

FOLDABLE = {"RNA": set("ACGUN"), "DNA": set("ACGTN")}


def load_params(molecule_type):
    if molecule_type == "DNA":
        RNA.params_load_DNA_Mathews2004()
    else:
        RNA.params_load_RNA_Turner2004()


def loop_stats(structure):
    runs = re.findall(r"\.+", structure)
    lead = len(structure) - len(structure.lstrip("."))
    trail = len(structure) - len(structure.rstrip(".")) if structure.strip(".") else 0
    return {
        "free_bases": structure.count("."),
        "loop_count": len(runs),
        "terminal_unpaired": lead + trail,
    }


def fold(seq, molecule_type, temperature):
    load_params(molecule_type)
    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(seq, md)
    structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    return {
        "structure": structure,
        "mfe_kcal_mol": round(float(mfe), 4),
        "base_pairs": structure.count("("),
        "gc_percent": round(100 * sum(b in "GC" for b in seq) / len(seq), 2),
        **loop_stats(structure),
        "ensemble_diversity": round(fc.mean_bp_distance(), 3),
    }


def analyse(row, temperature):
    seq, mol = row["fold_sequence"], row["detected_type"]
    if row.get("clean_status", "ok") != "ok":
        return {"fold_status": f"skipped ({row['clean_status']})"}
    ambiguous = set(seq) - FOLDABLE[mol]
    if ambiguous:
        return {"fold_status": f"skipped (ambiguous IUPAC {''.join(sorted(ambiguous))})"}
    try:
        return {**fold(seq, mol, temperature), "fold_status": "folded"}
    except Exception as exc:  # keep the run going, record the reason
        return {"fold_status": f"failed ({exc.__class__.__name__})"}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--temperature", type=float, default=25.0)
    args = parser.parse_args()

    df = pd.read_csv(args.input, keep_default_na=False)
    rows = [analyse(r, args.temperature) for _, r in df.iterrows()]
    out = pd.concat([df[["entry_id", "detected_type", "fold_sequence"]], pd.DataFrame(rows)], axis=1)
    out["temperature_c"] = args.temperature
    out["energy_params"] = out["detected_type"].map({"DNA": "Mathews2004 (DNA)", "RNA": "Turner2004 (RNA)"})
    RNA.params_load_RNA_Turner2004()

    out.to_csv(args.output, index=False)
    print(f"[fold] {out['fold_status'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
