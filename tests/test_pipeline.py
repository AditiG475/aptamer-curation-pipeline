"""Offline tests: network stages are tested with mocked API responses."""
import json
from unittest.mock import MagicMock

import pytest

from pipeline.s01_extract_binding import extract_kd
from pipeline.s02_doi_lookup import lookup_row, title_similarity
from pipeline.s03_pubchem_annotate import annotate
from pipeline.s04_classify import classify_journal, classify_ligand, load_rules
from pipeline.s05_structure_analysis import analyse_row

RULES = load_rules("config/classification_rules.json")


@pytest.mark.parametrize("text, expected_nm", [
    ("bound with a dissociation constant (Kd) of 0.32 µM", 320.0),
    ("giving Kd = 6 uM", 6000.0),
    ("a Kd of 117 nM", 117.0),
    ("apparent KD 2.5 mM", 2_500_000.0),
])
def test_extract_kd(text, expected_nm):
    assert extract_kd(text)[2] == pytest.approx(expected_nm)


def test_extract_kd_absent():
    assert extract_kd("no binding constant reported") == (None, None, None)


def fake_session(payload, status=200):
    response = MagicMock(status_code=status)
    response.json.return_value = payload
    response.raise_for_status.return_value = None
    session = MagicMock()
    session.get.return_value = response
    return session


def test_doi_verified():
    title = "A DNA aptamer that binds adenosine and ATP"
    session = fake_session({"message": {"items": [{"DOI": "10.0000/example", "title": [title]}]}})
    row = lookup_row(title, None, session, 0.85)
    assert row["doi_status"] == "verified" and row["doi"] == "10.0000/example"


def test_doi_needs_review():
    session = fake_session({"message": {"items": [{"DOI": "10.0000/other", "title": ["Unrelated paper on soil microbes"]}]}})
    assert lookup_row("A DNA aptamer that binds ATP", None, session, 0.85)["doi_status"] == "needs_review"


def test_title_similarity_ignores_case_and_hyphens():
    assert title_similarity("RNA-aptamer motifs", "rna aptamer motifs") == 1.0


def test_pubchem_found_and_missing():
    record = {"CID": 2153, "MolecularFormula": "C7H8N4O2", "MolecularWeight": "180.16", "SMILES": "CN1C2=C(C(=O)N(C1=O)C)NC=N2", "InChIKey": "ZFXYFBGIUFBOJW-UHFFFAOYSA-N"}
    assert annotate("theophylline", fake_session({"PropertyTable": {"Properties": [record]}}))["pubchem_cid"] == 2153
    assert annotate("not-a-compound", fake_session({}, status=404))["pubchem_status"] == "not_found"


@pytest.mark.parametrize("target, role", [
    ("Theophylline", "Biological target"),
    ("Magnesium chloride", "Buffer / salt component"),
    ("Glycerol", "Crystallographic additive"),
])
def test_ligand_roles(target, role):
    assert classify_ligand(target, RULES) == role


def test_journal_exact_match():
    assert classify_journal("Science", RULES) == "High-impact"
    assert classify_journal("Food Science and Technology", RULES) == "Unclassified"


def test_fold_rna_and_dna():
    rna = analyse_row("GGCGAUACCAGCCGAAAGGCCCUUGGCAGCGUC", "RNA", 25.0)
    dna = analyse_row("ACCTGGGGGAGTATTGCGGAGGAAGGT", "DNA", 25.0)
    assert rna["structure_status"] == dna["structure_status"] == "folded"
    assert rna["mfe_kcal_mol"] < 0 and len(rna["mfe_structure"]) == 33


def test_invalid_sequence_flagged():
    assert analyse_row("ACGTXX", "DNA", 25.0)["structure_status"].startswith("invalid")
