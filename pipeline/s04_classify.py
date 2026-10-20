"""Stage 4 - Rule-based classification of ligand role and journal tier.

Ligand role: separates real biological targets from buffer/salt components
and crystallographic additives that appear in structure or assay records.
Journal tier: groups source journals into editorial tiers.

Both rule sets live in config/classification_rules.json so they can be
edited without touching code.

Usage:
    python pipeline/s04_classify.py --input out/03_pubchem.csv --output out/04_classified.csv
"""
import argparse
import json
import re

import pandas as pd


def load_rules(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def classify_ligand(target, rules):
    """Keyword match on whole words, so 'sodium' does not match inside another word."""
    text = str(target).lower()
    for role, keywords in rules["ligand_roles"].items():
        for keyword in keywords:
            if re.search(rf"\b{re.escape(keyword)}\b", text):
                return role
    return rules["default_ligand_role"]


def classify_journal(journal, rules):
    """Exact match on the normalised journal name (avoids 'Science' matching 'Food Science')."""
    name = " ".join(str(journal).lower().replace("&", "and").split())
    for tier, journals in rules["journal_tiers"].items():
        if name in journals:
            return tier
    return rules["default_journal_tier"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--rules", default="config/classification_rules.json")
    args = parser.parse_args()

    rules = load_rules(args.rules)
    df = pd.read_csv(args.input)
    df["ligand_role"] = df["target"].apply(lambda t: classify_ligand(t, rules))
    df["journal_tier"] = df["journal"].apply(lambda j: classify_journal(j, rules))

    df.to_csv(args.output, index=False)
    print(f"[stage 4] Ligand roles: {df['ligand_role'].value_counts().to_dict()}")
    print(f"[stage 4] Journal tiers: {df['journal_tier'].value_counts().to_dict()} -> {args.output}")


if __name__ == "__main__":
    main()
