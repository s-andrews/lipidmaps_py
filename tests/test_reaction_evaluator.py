import pytest

from lipidmaps.data.utils.reaction_evaluator import ReactionEvaluator
from lipidmaps.data.models.reaction import ReactionData, CompoundComponent


def make_comp(headgroup=None, abbrev=None, full_struct=None, abbrev_chains=None):
    return CompoundComponent(
        compound_headgroup=headgroup,
        compound_abbrev=abbrev,
        compound_full_struct=full_struct,
        compound_abbrev_chains=abbrev_chains,
    )


def test_linkage_mismatch_detected():
    # Reactant: ester PC, Product: plasmalogen LPC (vinyl) -> should be linkage mismatch
    r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
    p = make_comp(headgroup="LPC P-", abbrev="LPC P-(18:1)")
    reaction = ReactionData(reactants=[r], products=[p])

    ev = ReactionEvaluator()
    res = ev.evaluate_reaction(reaction)
    assert res["possible"] is False
    assert "linkage_mismatch" in res["explanation"]


def test_sphingoid_requirement_enforced():
    # Cer -> SM require sphingoid evidence; without it should fail
    r = make_comp(headgroup="Cer", abbrev="Cer(16:0)")
    p = make_comp(headgroup="SM", abbrev="SM(16:0)")
    reaction = ReactionData(reactants=[r], products=[p])

    ev = ReactionEvaluator()
    res = ev.evaluate_reaction(reaction)
    assert res["possible"] is False
    assert "missing_sphingoid" in res["explanation"]

    # provide sphingoid evidence (d18:1) in reactant
    r2 = make_comp(headgroup="Cer", full_struct="d18:1/16:0")
    reaction2 = ReactionData(reactants=[r2], products=[p])
    res2 = ev.evaluate_reaction(reaction2)
    assert res2["possible"] is True


def test_acyl_chain_requirement_ok():
    # PC -> LPC requires a 1-chain fatty acid (lyso) evidence; either side can show it
    r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
    p = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
    reaction = ReactionData(reactants=[r], products=[p])

    ev = ReactionEvaluator()
    res = ev.evaluate_reaction(reaction)
    assert res["possible"] is True


def test_total_carbon_conservation_pass():
    # PC(16:0/18:1) -> LPC(18:1) + FA(16:0) should preserve total carbons (34)
    r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
    p1 = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
    p2 = CompoundComponent(compound_headgroup="FA", compound_full_struct="16:0")
    reaction = ReactionData(reactants=[r], products=[p1, p2])

    ev = ReactionEvaluator()
    res = ev.evaluate_reaction(reaction)
    assert res["possible"] is True


def test_total_carbon_conservation_mismatch_and_unknown():
    # Mismatch: PC(16:0/18:1) -> LPC(18:1) + FA(14:0) should fail carbon check
    r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
    p1 = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
    p2 = CompoundComponent(compound_headgroup="FA", compound_full_struct="14:0")
    reaction = ReactionData(reactants=[r], products=[p1, p2])

    ev = ReactionEvaluator()
    res = ev.evaluate_reaction(reaction)
    assert res["possible"] is False
    assert "carbon_mismatch" in res["explanation"]

    # Unknown chains: reactant missing chain info should lead to missing_chain_info
    r_unknown = CompoundComponent(compound_headgroup="PC", compound_name="PC")
    p1 = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
    p2 = CompoundComponent(compound_headgroup="FA", compound_full_struct="16:0")
    reaction2 = ReactionData(reactants=[r_unknown], products=[p1, p2])
    res2 = ev.evaluate_reaction(reaction2)
    import pytest

    from lipidmaps.data.utils.reaction_evaluator import ReactionEvaluator
    from lipidmaps.data.models.reaction import ReactionData, CompoundComponent


    def make_comp(headgroup=None, abbrev=None, full_struct=None, abbrev_chains=None):
        return CompoundComponent(
            compound_headgroup=headgroup,
            compound_abbrev=abbrev,
            compound_full_struct=full_struct,
            compound_abbrev_chains=abbrev_chains,
        )


    def test_linkage_mismatch_detected():
        # Reactant: ester PC, Product: plasmalogen LPC (vinyl) -> should be linkage mismatch
        r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
        p = make_comp(headgroup="LPC P-", abbrev="LPC P-(18:1)")
        reaction = ReactionData(reactants=[r], products=[p])

        ev = ReactionEvaluator()
        res = ev.evaluate_reaction(reaction)
        assert res["possible"] is False
        assert "linkage_mismatch" in res["explanation"]


    def test_sphingoid_requirement_enforced():
        # Cer -> SM require sphingoid evidence; without it should fail
        r = make_comp(headgroup="Cer", abbrev="Cer(16:0)")
        p = make_comp(headgroup="SM", abbrev="SM(16:0)")
        reaction = ReactionData(reactants=[r], products=[p])

        ev = ReactionEvaluator()
        res = ev.evaluate_reaction(reaction)
        assert res["possible"] is False
        assert "missing_sphingoid" in res["explanation"]

        # provide sphingoid evidence (d18:1) in reactant
        r2 = make_comp(headgroup="Cer", full_struct="d18:1/16:0")
        reaction2 = ReactionData(reactants=[r2], products=[p])
        res2 = ev.evaluate_reaction(reaction2)
        assert res2["possible"] is True


    def test_acyl_chain_requirement_ok():
        # PC -> LPC requires a 1-chain fatty acid (lyso) evidence; either side can show it
        r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
        p = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
        reaction = ReactionData(reactants=[r], products=[p])

        ev = ReactionEvaluator()
        res = ev.evaluate_reaction(reaction)
        assert res["possible"] is True


    def test_total_carbon_conservation_pass():
        # PC(16:0/18:1) -> LPC(18:1) + FA(16:0) should preserve total carbons (34)
        r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
        p1 = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
        p2 = CompoundComponent(compound_headgroup="FA", compound_full_struct="16:0")
        reaction = ReactionData(reactants=[r], products=[p1, p2])

        ev = ReactionEvaluator()
        res = ev.evaluate_reaction(reaction)
        assert res["possible"] is True


    def test_total_carbon_conservation_mismatch_and_unknown():
        # Mismatch: PC(16:0/18:1) -> LPC(18:1) + FA(14:0) should fail carbon check
        r = make_comp(headgroup="PC", abbrev="PC(16:0/18:1)")
        p1 = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
        p2 = CompoundComponent(compound_headgroup="FA", compound_full_struct="14:0")
        reaction = ReactionData(reactants=[r], products=[p1, p2])

        ev = ReactionEvaluator()
        res = ev.evaluate_reaction(reaction)
        assert res["possible"] is False
        assert "carbon_mismatch" in res["explanation"]

        # Unknown chains: reactant missing chain info should lead to missing_chain_info
        r_unknown = CompoundComponent(compound_headgroup="PC", compound_name="PC")
        p1 = make_comp(headgroup="LPC", abbrev="LPC(18:1)")
        p2 = CompoundComponent(compound_headgroup="FA", compound_full_struct="16:0")
        reaction2 = ReactionData(reactants=[r_unknown], products=[p1, p2])
        res2 = ev.evaluate_reaction(reaction2)
        assert res2["possible"] is False
        assert "missing_chain_info" in res2["explanation"]
