"""Stage 3 - Annotate target compounds with PubChem (PUG-REST).

Looks up each target name and adds the PubChem CID, molecular formula,
molecular weight, SMILES and InChIKey. Names PubChem cannot resolve are
kept and marked "not_found" so they can be curated by hand.

Usage:
    python pipeline/s03_pubchem_annotate.py --input out/02_doi.csv --output out/03_pubchem.csv
"""
import argparse
import time
from urllib.parse import quote

import pandas as pd
import requests

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/{props}/JSON"
PROPERTIES = ["MolecularFormula", "MolecularWeight", "SMILES", "InChIKey"]


def fetch_properties(name, session=None, timeout=20):
    """Return a dict of PubChem properties for a compound name, or None if not found."""
    session = session or requests.Session()
    url = PUG.format(name=quote(str(name)), props=",".join(PROPERTIES))
    response = session.get(url, timeout=timeout)
    if response.status_code == 404:
        return None
    response.raise_for_status()
    records = response.json().get("PropertyTable", {}).get("Properties", [])
    return records[0] if records else None


def annotate(name, session):
    empty = {"pubchem_cid": None, **{p.lower(): None for p in PROPERTIES}}
    try:
        record = fetch_properties(name, session)
    except requests.RequestException as exc:
        return {**empty, "pubchem_status": f"error: {exc.__class__.__name__}"}
    if record is None:
        return {**empty, "pubchem_status": "not_found"}
    row = {"pubchem_cid": record.get("CID")}
    row.update({p.lower(): record.get(p) for p in PROPERTIES})
    row["pubchem_status"] = "found"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--name-column", default="target")
    parser.add_argument("--delay", type=float, default=0.25,
                        help="Seconds between requests (PubChem allows max 5 requests/second)")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    session = requests.Session()
    cache, rows = {}, []
    for name in df[args.name_column]:
        if name not in cache:
            cache[name] = annotate(name, session)
            time.sleep(args.delay)
        rows.append(cache[name])

    out = pd.concat([df, pd.DataFrame(rows)], axis=1)
    out.to_csv(args.output, index=False)
    print(f"[stage 3] PubChem status counts: {out['pubchem_status'].value_counts().to_dict()} -> {args.output}")


if __name__ == "__main__":
    main()
