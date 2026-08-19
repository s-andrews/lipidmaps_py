from lipidmaps.data.models.reaction import CompoundComponent, ReactionData
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
from lipidmaps.data.utils.reaction_evaluator import ReactionEvaluator


def test_reaction_evaluator_matches_generic_components_and_fa_delta():
    dataset = LipidDataset(
        samples=[
            SampleMetadata(sample_name="ctrl_1", group="control"),
            SampleMetadata(sample_name="ctrl_2", group="control"),
        ],
        lipids=[
            QuantifiedLipid(
                input_name="PC(16:0_18:1)",
                standardized_name="PC 16:0_18:1",
                generic_lm_id="LMGP01010000",
                values={"ctrl_1": 1.0, "ctrl_2": 2.0},
            ),
            QuantifiedLipid(
                input_name="LPC(18:1)",
                standardized_name="LPC 18:1",
                generic_lm_id="LMGP01050000",
                values={"ctrl_1": 1.0, "ctrl_2": 2.0},
            ),
            QuantifiedLipid(
                input_name="FA(16:0)",
                standardized_name="FA 16:0",
                generic_lm_id="LMFA01010000",
                values={"ctrl_1": 1.0, "ctrl_2": 2.0},
            ),
        ],
    )

    reaction = ReactionData(
        reactants=[
            CompoundComponent(
                compound_type="lm_main",
                compound_generic_lm_id="LMGP01010000",
            )
        ],
        products=[
            CompoundComponent(
                compound_type="lm_main",
                compound_generic_lm_id="LMGP01050000",
            )
        ],
    )

    result = ReactionEvaluator().evaluate_reaction(reaction, dataset=dataset)

    assert result["possible"] is True
    assert result["pairs_info"] == ["PC 16:0_18:1 -> LPC 18:1"]


def test_reaction_data_nonfa_helpers_include_generic_lm_ids():
    reaction = ReactionData(
        reactants=[CompoundComponent(compound_type="lm_main", compound_generic_lm_id="LMGP01010000")],
        products=[CompoundComponent(compound_type="lm_main", compound_generic_lm_id="LMGP01050000")],
    )

    assert reaction.list_nonfa_noncoa_reactant_lm_ids() == ["LMGP01010000"]
    assert reaction.list_nonfa_noncoa_product_lm_ids() == ["LMGP01050000"]