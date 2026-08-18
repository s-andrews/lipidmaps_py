"""Tests for shunt / alternate-route pathway detection."""

from lipidmaps.data.utils.shunt import is_shunt_pathway, reaction_is_shunt


def test_is_shunt_pathway_from_name():
    assert is_shunt_pathway({"pathway_name": "Cholesterol biosynthesis (Shunt)"})
    assert is_shunt_pathway({"name": "Sterol alternate route"})


def test_is_shunt_pathway_from_type():
    assert is_shunt_pathway({"pathway_type": ["Lipid biosynthesis", "Shunt"]})
    assert is_shunt_pathway({"type": ["shunt"]})


def test_is_shunt_pathway_from_structured_flag():
    assert is_shunt_pathway({"is_shunt": True})
    assert is_shunt_pathway({"shunt": True})


def test_is_shunt_pathway_negative():
    assert not is_shunt_pathway({"pathway_name": "Phosphatidylcholine turnover"})
    assert not is_shunt_pathway({})
    assert not is_shunt_pathway(None)


def test_structured_false_overrides_name_fallback():
    # A curator's explicit is_shunt=False is authoritative even if the name
    # happens to contain a shunt marker.
    assert not is_shunt_pathway(
        {"is_shunt": False, "pathway_name": "Cholesterol biosynthesis (Shunt)"}
    )


def test_reaction_is_shunt_any_pathway():
    pathways = [
        {"pathway_name": "Cholesterol biosynthesis"},
        {"pathway_name": "Cholesterol biosynthesis (Shunt)"},
    ]
    assert reaction_is_shunt(pathways)


def test_reaction_is_shunt_none_match():
    pathways = [
        {"pathway_name": "Cholesterol biosynthesis"},
        {"pathway_name": "Phosphatidylcholine turnover"},
    ]
    assert not reaction_is_shunt(pathways)
    assert not reaction_is_shunt([])
    assert not reaction_is_shunt(None)
