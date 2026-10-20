"""Stage 5 - Secondary structure analysis with ViennaRNA.

For every aptamer sequence: minimum free energy (MFE) structure in
dot-bracket notation, MFE, length, GC content, number of base pairs and
ensemble diversity. RNA uses the default Turner 2004 parameters; DNA
sequences are folded with the Mathews 2004 DNA parameter set.

Usage:
    python pipeline/s05_structure_analysis.py --input out/04_classified.csv --output out/05_structures.csv
"""
import argparse

import pandas as pd
import RNA

VALID = {"RNA": set("ACGU"), "DNA": set("ACGT")}


def clean_sequence(seq, molecule_type):
    seq = "".join(str(seq).upper().split())
    if molecule_type == "RNA":
        seq = seq.replace("T", "U")
    return seq


def fold(seq, molecule_type, temperature=25.0):
    """Fold one sequence and return a dict of structural descriptors."""
    md = RNA.md()
    md.temperature = temperature
    if molecule_type == "DNA":
        RNA.params_load_DNA_Mathews2004()
    else:
        RNA.params_load_RNA_Turner2004()

    fc = RNA.fold_compound(seq, md)
    structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    return {
        "length": len(seq),
        "gc_content": round(100 * sum(base in "GC" for base in seq) / len(seq), 1),
        "mfe_structure": structure,
        "mfe_kcal_mol": round(mfe, 2),
        "base_pairs": structure.count("("),
        "ensemble_diversity": round(fc.mean_bp_distance(), 2),
    }


def analyse_row(seq, molecule_type, temperature):
    molecule_type = str(molecule_type).upper()
    if molecule_type not in VALID:
        return {"structure_status": f"unknown molecule_type '{molecule_type}'"}
    seq = clean_sequence(seq, molecule_type)
    invalid = set(seq) - VALID[molecule_type]
    if not seq or invalid:
        return {"structure_status": f"invalid characters {sorted(invalid)}"}
    return {**fold(seq, molecule_type, temperature), "structure_status": "folded"}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--temperature", type=float, default=25.0,
                        help="Folding temperature in Celsius (match your binding assay)")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    rows = [analyse_row(s, m, args.temperature) for s, m in zip(df["sequence"], df["molecule_type"])]
    out = pd.concat([df, pd.DataFrame(rows)], axis=1)
    RNA.params_load_RNA_Turner2004()  # restore defaults

    out.to_csv(args.output, index=False)
    print(f"[stage 5] Structure status: {out['structure_status'].value_counts().to_dict()} -> {args.output}")


if __name__ == "__main__":
    main()
