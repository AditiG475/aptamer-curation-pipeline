#!/usr/bin/env python3
"""Annotate target compounds with physicochemical properties from PubChem.

Identifier cascade (first one that resolves wins, and is recorded in
`pubchem_source`):  CID -> SMILES -> InChIKey -> compound name.
Requests are retried with back-off; "not found" and network errors are
kept as status values instead of being dropped.

Usage:
    pubchem_annotate.py --input entries.csv --output annotated.csv \
        [--cid-column cid] [--smiles-column smiles] [--inchikey-column inchikey] [--name-column target]
"""
import argparse
import time
from urllib.parse import quote

import pandas as pd
import requests

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPERTIES = [
    "MolecularFormula", "MolecularWeight", "SMILES", "InChIKey", "XLogP", "TPSA",
    "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount",
    "Charge", "Complexity", "DefinedAtomStereoCount", "UndefinedAtomStereoCount",
    "CovalentUnitCount",
]


def request_json(session, method, url, data=None, retries=3, delay=0.25, timeout=20):
    """GET/POST with retries. Returns parsed JSON, None for 'not found'; raises on repeated errors."""
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = session.request(method, url, data=data, timeout=timeout)
            if response.status_code == 404:
                return None
            if response.status_code == 200:
                return response.json()
            last_error = requests.HTTPError(f"HTTP {response.status_code}")
        except requests.RequestException as exc:
            last_error = exc
        time.sleep(delay * attempt)
    raise last_error


def lookup(session, namespace, identifier):
    props = ",".join(PROPERTIES)
    if namespace == "smiles":  # POST avoids URL-escaping problems with SMILES
        data = request_json(session, "POST", f"{BASE}/smiles/property/{props}/JSON", data={"smiles": identifier})
    else:
        data = request_json(session, "GET", f"{BASE}/{namespace}/{quote(str(identifier), safe='')}/property/{props}/JSON")
    records = (data or {}).get("PropertyTable", {}).get("Properties", [])
    return records[0] if records else None


def annotate(row, columns, session):
    errors = []
    for namespace, column in columns:
        value = row.get(column) if column else None
        if value is None or (isinstance(value, float) and pd.isna(value)) or str(value).strip() == "":
            continue
        value = str(int(float(value))) if namespace == "cid" else str(value).strip()
        try:
            record = lookup(session, namespace, value)
        except requests.RequestException as exc:
            errors.append(f"{namespace}:{exc.__class__.__name__}")
            continue
        if record:
            out = {"pubchem_cid": record.get("CID"), "pubchem_source": namespace, "pubchem_status": "found"}
            out.update({f"pc_{p}": record.get(p) for p in PROPERTIES})
            return out
    status = f"error ({', '.join(errors)})" if errors else "not_found"
    return {"pubchem_cid": None, "pubchem_source": None, "pubchem_status": status}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cid-column", default="cid")
    parser.add_argument("--smiles-column", default="smiles")
    parser.add_argument("--inchikey-column", default="inchikey")
    parser.add_argument("--name-column", default="target")
    parser.add_argument("--delay", type=float, default=0.25, help="PubChem allows at most 5 requests/second")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    columns = [(ns, col) for ns, col in [("cid", args.cid_column), ("smiles", args.smiles_column),
                                         ("inchikey", args.inchikey_column), ("name", args.name_column)]
               if col in df.columns]
    session, cache, rows = requests.Session(), {}, []
    for _, row in df.iterrows():
        key = tuple(str(row.get(c)) for _, c in columns)
        if key not in cache:
            cache[key] = annotate(row, columns, session)
            time.sleep(args.delay)
        rows.append(cache[key])

    out = pd.concat([df, pd.DataFrame(rows)], axis=1)
    out.to_csv(args.output, index=False)
    print(f"[pubchem] {out['pubchem_status'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
