"""Unit tests for the CoreMetabolome m/z annotator (pure, offline)."""

import pytest

from lipidmaps.data.annotation.mz_annotator import (
    MetaboliteDatabase,
    default_adducts,
    formula_monoisotopic_mass,
    POSITIVE_ADDUCTS,
    NEGATIVE_ADDUCTS,
)


def test_formula_monoisotopic_mass_known_values():
    assert formula_monoisotopic_mass("H2O") == pytest.approx(18.0106, abs=1e-3)
    assert formula_monoisotopic_mass("C6H12O6") == pytest.approx(180.0634, abs=1e-3)
    # Single-atom counts and multi-letter elements.
    assert formula_monoisotopic_mass("CH4") == pytest.approx(16.0313, abs=1e-3)
    assert formula_monoisotopic_mass("NaCl") == pytest.approx(57.9586, abs=1e-3)


def test_formula_parse_errors():
    with pytest.raises(ValueError):
        formula_monoisotopic_mass("")
    with pytest.raises(KeyError):
        formula_monoisotopic_mass("Xy2")  # unknown element


def test_adduct_mz_and_defaults():
    m = formula_monoisotopic_mass("C6H12O6")
    proton = POSITIVE_ADDUCTS[0]
    assert proton.name == "[M+H]+"
    assert proton.mz(m) == pytest.approx(181.0707, abs=1e-3)
    assert default_adducts("positive") == POSITIVE_ADDUCTS
    assert default_adducts("negative mode") == NEGATIVE_ADDUCTS
    assert default_adducts(None) == POSITIVE_ADDUCTS


def test_database_groups_isomers_by_formula():
    db = MetaboliteDatabase(
        {
            "C6H12O6": {
                "mass": formula_monoisotopic_mass("C6H12O6"),
                "candidates": [
                    {"name": "Glucose", "id": "HMDB1"},
                    {"name": "Fructose", "id": "HMDB2"},
                ],
            }
        }
    )
    db.build_index(POSITIVE_ADDUCTS)
    matches = db.annotate([181.0707], ppm_tol=5)
    assert len(matches) == 1
    m = matches[0]
    assert m.formula == "C6H12O6"
    assert m.adduct == "[M+H]+"
    assert m.ppm_error < 5
    assert {c["name"] for c in m.candidates} == {"Glucose", "Fructose"}
    assert m.label == "C6H12O6 [M+H]+"


def test_annotate_respects_tolerance():
    db = MetaboliteDatabase(
        {"C6H12O6": {"mass": formula_monoisotopic_mass("C6H12O6"), "candidates": [{"name": "Glc", "id": ""}]}}
    )
    db.build_index(POSITIVE_ADDUCTS)
    # 50 ppm off the [M+H]+ target -> no match at 5 ppm.
    target = POSITIVE_ADDUCTS[0].mz(formula_monoisotopic_mass("C6H12O6"))
    off = target * (1 + 50e-6)
    assert db.annotate([off], ppm_tol=5) == []
    assert db.annotate([off], ppm_tol=100) != []


def test_load_real_coremetabolome(tmp_path):
    from pathlib import Path

    db_path = Path(__file__).resolve().parent / "data" / "core_metabolome_v3.csv"
    db = MetaboliteDatabase.load(db_path)
    # Thousands of unique formulas from the 11k-row DB.
    assert len(db.by_formula) > 1000
    db.build_index(default_adducts("positive"))
    matches = db.annotate([181.0707], ppm_tol=5)
    formulas = {m.formula for m in matches}
    assert "C6H12O6" in formulas
