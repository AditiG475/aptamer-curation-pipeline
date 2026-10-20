#!/usr/bin/env python3
"""Tertiary (3D) structure prediction for one aptamer with RhoFold+.

Wraps RhoFold+ single-sequence inference (https://github.com/ml4bio/RhoFold)
and writes a one-row summary CSV next to the predicted PDB.

RhoFold+ is trained on RNA only. DNA aptamers are therefore modelled as
their RNA analogue (T -> U) and the row is flagged `rna_analogue_of_dna`;
treat those models as approximate.

Usage:
    predict_tertiary.py --entry-id APT001 --sequence GGCGAU... --molecule-type RNA \
        --rhofold-dir /path/to/RhoFold --ckpt /path/to/RhoFold_pretrained.pt \
        --outdir APT001_3d --device cpu
"""
import argparse
import os
import shutil
import subprocess
import sys

import pandas as pd


def write_fasta(path, entry_id, seq):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f">{entry_id}\n{seq}\n")


def pick_model(outdir):
    """Prefer the relaxed model, fall back to the unrelaxed one."""
    files = sorted(os.listdir(outdir)) if os.path.isdir(outdir) else []
    relaxed = [f for f in files if f.startswith("relaxed") and f.endswith(".pdb")]
    unrelaxed = [f for f in files if f.startswith("unrelaxed") and f.endswith(".pdb")]
    chosen = (relaxed or unrelaxed or [None])[0]
    return os.path.join(outdir, chosen) if chosen else None


def count_atoms(pdb_path):
    with open(pdb_path, encoding="utf-8") as handle:
        return sum(1 for line in handle if line.startswith(("ATOM", "HETATM")))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--entry-id", required=True)
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--molecule-type", choices=["DNA", "RNA"], required=True)
    parser.add_argument("--rhofold-dir", required=True, help="Folder containing RhoFold's inference.py")
    parser.add_argument("--ckpt", required=True, help="RhoFold+ pretrained checkpoint (.pt)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--device", default="cpu", help="cpu or cuda:0")
    parser.add_argument("--relax-steps", type=int, default=1000)
    parser.add_argument("--max-length", type=int, default=200)
    parser.add_argument("--python", default=sys.executable, help="Python of the RhoFold environment")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    seq = args.sequence.upper().replace("T", "U")
    flag = "rna_analogue_of_dna" if args.molecule_type == "DNA" else "ok"
    summary = {"entry_id": args.entry_id, "model_sequence": seq, "tertiary_tool": "RhoFold+",
               "tertiary_flag": flag, "pdb_file": None, "atom_count": None}

    if len(seq) > args.max_length:
        summary["tertiary_status"] = f"skipped (length {len(seq)} > {args.max_length})"
    else:
        fasta = os.path.join(args.outdir, f"{args.entry_id}.fasta")
        write_fasta(fasta, args.entry_id, seq)
        command = [args.python, os.path.join(args.rhofold_dir, "inference.py"),
                   "--input_fas", fasta, "--single_seq_pred", "True",
                   "--output_dir", args.outdir, "--ckpt", args.ckpt,
                   "--device", args.device, "--relax_steps", str(args.relax_steps)]
        result = subprocess.run(command, capture_output=True, text=True)
        with open(os.path.join(args.outdir, "rhofold.log"), "w", encoding="utf-8") as log:
            log.write(result.stdout + "\n" + result.stderr)
        model = pick_model(args.outdir)
        if result.returncode != 0 or model is None:
            summary["tertiary_status"] = f"failed (exit {result.returncode}; see rhofold.log)"
        else:
            final = os.path.join(args.outdir, f"{args.entry_id}_3d.pdb")
            shutil.copyfile(model, final)
            summary.update(pdb_file=os.path.basename(final), atom_count=count_atoms(final),
                           tertiary_status="predicted")

    pd.DataFrame([summary]).to_csv(os.path.join(args.outdir, f"{args.entry_id}_3d_summary.csv"), index=False)
    print(f"[3d] {args.entry_id}: {summary['tertiary_status']} ({flag})")


if __name__ == "__main__":
    main()
