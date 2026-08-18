"""Tests for the named-compound (sterol / shunt) reaction matcher."""

from lipidmaps.data.matching import match_pathway_reactions
from lipidmaps.data.matching.base import create_matcher_context
from lipidmaps.data.matching.same_structure import SameStructureMatcher
from lipidmaps.data.models.species_reaction import ClassReaction, ReactionType


def _sterol_reaction(is_shunt: bool = True) -> ClassReaction:
    return ClassReaction(
        reactant_class="lanosterol",
        product_class="cholesterol",
        is_shunt=is_shunt,
        is_molspecies=True,
    )


def test_sterol_reaction_routes_to_named_compound():
    assert _sterol_reaction().reaction_type == ReactionType.NAMED_COMPOUND
    # Also routes without the shunt flag, because both classes are in the
    # sterol/named allow-set.
    assert _sterol_reaction(is_shunt=False).reaction_type == ReactionType.NAMED_COMPOUND


def test_named_compound_pairs_when_both_present():
    rxn = _sterol_reaction()
    result_set = match_pathway_reactions(
        lipid_names=["lanosterol", "cholesterol"], reactions=[rxn]
    )
    result = result_set.get_result("lanosterol", "cholesterol")
    assert result is not None
    assert result.pairs_matched == 1


def test_named_compound_no_pair_when_one_side_missing():
    rxn = _sterol_reaction()
    result_set = match_pathway_reactions(
        lipid_names=["lanosterol"], reactions=[rxn]
    )
    result = result_set.get_result("lanosterol", "cholesterol")
    assert result is not None
    assert result.pairs_matched == 0


def test_ordinary_reaction_unaffected():
    # A normal glycerophospholipid reaction must still be SAME_STRUCTURE and
    # match by composition, not identity.
    pc_pa = ClassReaction(reactant_class="PC", product_class="PA")
    assert pc_pa.reaction_type == ReactionType.SAME_STRUCTURE
    result_set = match_pathway_reactions(
        lipid_names=["PC 34:1", "PA 34:1", "PA 36:2"], reactions=[pc_pa]
    )
    result = result_set.get_result("PC", "PA")
    # Only the matching composition (34:1) pairs.
    assert result.pairs_matched == 1


def test_shunt_flag_survives_json_round_trip():
    rxn = _sterol_reaction()
    restored = ClassReaction.model_validate(rxn.model_dump(mode="json"))
    assert restored.is_shunt is True
    assert restored.reaction_type == ReactionType.NAMED_COMPOUND


def test_same_structure_matcher_is_not_the_intended_path_for_sterols():
    # Sterols parse MOLECULAR -> (0,0), so SAME_STRUCTURE would spuriously
    # collide unrelated sterols at (0,0); NAMED_COMPOUND is the intentional
    # matcher. This documents why routing diverts these reactions.
    rxn = _sterol_reaction()
    ctx = create_matcher_context(["lanosterol", "cholesterol"])
    ss_result = SameStructureMatcher().match(rxn, ctx)
    # The degenerate (0,0) collision pairs them by accident, not by design.
    assert ss_result.pairs_matched == 1
