"""Offline tests. Network stages (CrossRef, PubChem) use mocked responses."""
from unittest.mock import MagicMock

import pytest

from classify import classify_journal, classify_ligand, load_rules
from clean_sequences import clean_one
from detect_molecule_type import detect, to_alphabet
from doi_lookup import lookup_row, title_similarity
from extract_binding import extract_kd
from fold_secondary import analyse, fold, loop_stats
from pubchem_annotate import annotate
from shuffle_zscore import mfe, zscore
from verify_structures import balanced

RULES = load_rules("config/classification_rules.json")


# --- cleaning ---------------------------------------------------------------
@pytest.mark.parametrize("raw, expected", [
    ("ggcg auac\tcagc\n", "GGCGAUACCAGC"),
    ("5'-ACGTACGT-3'", "ACGTACGT"),
    ("1 GGAUCC 6", "GGAUCC"),
    ("AC-GU.AC", "ACGUAC"),
    ("ACGU*XACGU", "ACGUACGU"),
])
def test_clean_sequence(raw, expected):
    assert clean_one(raw)[0] == expected


def test_clean_reports_changes():
    _, log = clean_one("acg u")
    assert any("whitespace" in c for c in log) and any("uppercase" in c for c in log)


# --- DNA / RNA detection ----------------------------------------------------
@pytest.mark.parametrize("seq, declared, mol, flag", [
    ("ACGTTGCA", None, "DNA", "ok"),
    ("ACGUUGCA", None, "RNA", "ok"),
    ("ACGTTGUA", None, "DNA", "mixed_T_and_U"),
    ("GGGCCCAG", "DNA", "DNA", "type_from_metadata"),
    ("GGGCCCAG", None, "RNA", "type_assumed"),
])
def test_detect(seq, declared, mol, flag):
    result = detect(seq, declared, "RNA")
    assert result[0] == mol and result[2] == flag


def test_declared_mismatch_flagged():
    assert "declared_mismatch(RNA)" in detect("ACGTTGCA", "RNA", "RNA")[2]


def test_alphabet_conversion():
    assert to_alphabet("ACGU", "DNA") == "ACGT" and to_alphabet("ACGT", "RNA") == "ACGU"


# --- secondary structure & Z-score -----------------------------------------
def test_dna_and_rna_use_different_parameters():
    dna = "ACCTGGGGGAGTATTGCGGAGGAAGGT"
    rna = dna.replace("T", "U")
    assert fold(dna, "DNA", 25.0)["mfe_kcal_mol"] != fold(rna, "RNA", 25.0)["mfe_kcal_mol"]
    # switching back and forth must give identical energies (no parameter leakage)
    assert mfe(dna, "DNA", 25.0) == mfe(dna, "DNA", 25.0)
    assert fold(dna, "DNA", 25.0)["mfe_kcal_mol"] == pytest.approx(mfe(dna, "DNA", 25.0), abs=1e-3)


def test_loop_stats():
    assert loop_stats("..((..))..") == {"free_bases": 6, "loop_count": 3, "terminal_unpaired": 4}


def test_skip_ambiguous_and_short():
    assert analyse({"fold_sequence": "ACGRACGU", "detected_type": "RNA", "clean_status": "ok"}, 25)["fold_status"].startswith("skipped")
    assert analyse({"fold_sequence": "ACGU", "detected_type": "RNA", "clean_status": "too_short(<8)"}, 25)["fold_status"].startswith("skipped")


def test_zscore_reproducible():
    import random
    seq = "GGAUCCCGACUGGCGAGAGCCAGGUAACGAAUGGAUCC"
    a = zscore(seq, "RNA", 20, 25.0, random.Random(1))
    b = zscore(seq, "RNA", 20, 25.0, random.Random(1))
    assert a == b and a["zscore_status"] == "ok" and a["z_score"] < 0


def test_balanced():
    assert balanced("((..))") and not balanced("(()") and not balanced(")(")


# --- annotation -------------------------------------------------------------
@pytest.mark.parametrize("text, expected_nm", [
    ("a dissociation constant (Kd) of 0.32 µM", 320.0),
    ("giving Kd = 6 uM", 6000.0),
    ("a Kd of 117 nM", 117.0),
])
def test_extract_kd(text, expected_nm):
    assert extract_kd(text)[2] == pytest.approx(expected_nm)


def fake_session(payload, status=200):
    response = MagicMock(status_code=status)
    response.json.return_value = payload
    response.raise_for_status.return_value = None
    session = MagicMock()
    session.get.return_value = response
    session.request.return_value = response
    return session


def test_doi_verified_and_review():
    title = "A DNA aptamer that binds adenosine and ATP"
    ok = fake_session({"message": {"items": [{"DOI": "10.1021/bi00002a033", "title": [title]}]}})
    bad = fake_session({"message": {"items": [{"DOI": "10.0/x", "title": ["Soil microbes"]}]}})
    assert lookup_row(title, None, ok, 0.85)["doi_status"] == "verified"
    assert lookup_row(title, None, bad, 0.85)["doi_status"] == "needs_review"
    assert title_similarity("RNA-aptamer motifs", "rna aptamer motifs") == 1.0


def test_pubchem_cascade():
    record = {"CID": 2153, "MolecularWeight": "180.16", "SMILES": "CN1C2=C(C(=O)N(C1=O)C)NC=N2"}
    found = annotate({"target": "theophylline"}, [("name", "target")], fake_session({"PropertyTable": {"Properties": [record]}}))
    assert found["pubchem_cid"] == 2153 and found["pubchem_source"] == "name"
    missing = annotate({"target": "unknown"}, [("name", "target")], fake_session({}, status=404))
    assert missing["pubchem_status"] == "not_found"


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
