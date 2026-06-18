"""Regression tests for BioPAN FA/CoA pool construction and ratio scoring.

These lock in two fixes:

1. The FA / acyl-CoA pools are the *union* of reaction-extracted, dataset-named
   and dataset-inferred species, so a sparse reaction-derived set (e.g. only
   FA 16:0) can no longer short-circuit the richer inferred set and starve
   FA-release / FACoA-addition matching. This previously inflated subclass
   z-scores (e.g. PS->LPS) by scoring from only a "lose-a-palmitate" subset.

2. `_compute_ratio_zscore` follows the legacy R rules for missing/zero values:
   any missing per-sample aggregate makes the reaction unscorable (z = 0), and
   zero-reactant samples are dropped (with R's union-index behaviour for
   balanced unpaired designs and the ``l < n - 1`` guard).
"""

from lipidmaps.data import BioPANPathwayExporter
from lipidmaps.data.models.reaction import CompoundComponent, ReactionData
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
from lipidmaps.data.models.species_reaction import ClassReaction, CompoundRequirement


def _fa_reaction(abbrev: str, lm_id: str = "LMFA01010001") -> ReactionData:
    """A reaction whose product is an LM-ID'd fatty acid (feeds extract_fa)."""
    return ReactionData(
        reactants=[CompoundComponent(compound_name="PS", compound_headgroup="PS")],
        products=[CompoundComponent(compound_abbrev=abbrev, compound_name=abbrev, compound_lm_id=lm_id)],
    )


def _ps_lps_dataset() -> LipidDataset:
    """Sum-level PS reactants + LPS products across several acyl chains.

    The only reaction-derived FA is the sparse FA 16:0; the inferred pool from
    the LPS lyso species is much richer (18:0, 18:1, 20:4, 22:6).
    """
    samples = [SampleMetadata(sample_name=f"ctrl_{i}", group="control") for i in range(1, 5)]
    samples += [SampleMetadata(sample_name=f"case_{i}", group="case") for i in range(1, 5)]

    def lipid(name, ctrl, case):
        values = {f"ctrl_{i}": ctrl[i - 1] for i in range(1, 5)}
        values.update({f"case_{i}": case[i - 1] for i in range(1, 5)})
        return QuantifiedLipid(input_name=name, values=values)

    lipids = [
        lipid("PS(34:0)", [20, 22, 19, 21], [18, 17, 19, 16]),
        lipid("PS(34:1)", [30, 31, 29, 30], [25, 26, 24, 25]),
        lipid("PS(36:4)", [40, 42, 39, 41], [35, 34, 36, 33]),
        lipid("PS(38:6)", [15, 16, 14, 15], [12, 13, 11, 12]),
        lipid("LPS(18:0)", [2.0, 2.2, 1.8, 2.1], [3.5, 3.8, 3.2, 3.6]),
        lipid("LPS(18:1)", [1.5, 1.6, 1.4, 1.5], [2.4, 2.5, 2.3, 2.4]),
        lipid("LPS(20:4)", [1.0, 1.1, 0.9, 1.0], [1.6, 1.7, 1.5, 1.6]),
        lipid("LPS(22:6)", [0.8, 0.9, 0.7, 0.8], [1.3, 1.4, 1.2, 1.3]),
    ]
    dataset = LipidDataset(samples=samples, lipids=lipids)
    dataset.reactions = [_fa_reaction("FA 16:0")]
    return dataset


def test_fa_pool_is_union_not_short_circuited():
    dataset = _ps_lps_dataset()
    exporter = BioPANPathwayExporter(dataset=dataset)

    _, fa_names, _ = exporter._get_matching_inputs(dataset)

    # Sparse reaction FA is present, but did NOT prevent the richer inferred set.
    assert "FA 16:0" in fa_names
    for inferred in ("FA 18:0", "FA 18:1", "FA 20:4", "FA 22:6"):
        assert inferred in fa_names, f"{inferred} missing from unioned FA pool"
    assert len(fa_names) > 1


def test_release_reaction_scores_full_matched_species_set():
    dataset = _ps_lps_dataset()
    exporter = BioPANPathwayExporter(
        dataset=dataset,
        class_reactions=[ClassReaction(reactant_class="PS", product_class="LPS", compound_require=CompoundRequirement.FA)],
    )
    result_set, _ = exporter.build_reaction_match_set(dataset)
    diag = exporter.diagnose_reaction(
        "PS,LPS", disease_group="case", control_group="control", level="class", mode="active", result_set=result_set,
    )

    # All four LPS products participate, not just the FA-16:0 subset.
    assert diag["product_species"] == ["LPS(18:0)", "LPS(18:1)", "LPS(20:4)", "LPS(22:6)"]
    assert diag["pairs_matched"] >= 4
    assert diag["n_disease_used"] == 4 and diag["n_control_used"] == 4
    # LPS rises relative to PS in the case group -> active (positive) z.
    assert diag["z_score"] > 0


def test_missing_aggregate_makes_reaction_unscorable():
    exporter = BioPANPathwayExporter()
    # A None product sum on one sample -> R's na_values guard -> z 0.
    z = exporter._compute_ratio_zscore(
        disease_products=[2.0, None, 2.2, 2.1],
        disease_reactants=[10.0, 10.0, 11.0, 10.5],
        control_products=[1.0, 1.1, 0.9, 1.0],
        control_reactants=[10.0, 10.0, 10.0, 10.0],
        alt="greater",
        paired=False,
    )
    assert z == 0.0


def test_zero_reactant_samples_dropped_union_index():
    exporter = BioPANPathwayExporter()
    # ctrl position 0 has a zero reactant -> that position is dropped from BOTH
    # groups (R union-index behaviour), leaving 3 paired-by-position samples.
    selected = exporter._select_ratio_vectors(
        disease_products=[2.0, 2.1, 2.2, 2.3],
        disease_reactants=[10.0, 10.0, 10.0, 10.0],
        control_products=[1.0, 1.1, 0.9, 1.0],
        control_reactants=[0.0, 10.0, 10.0, 10.0],
        paired=False,
    )
    assert selected is not None
    disease_ratios, control_ratios = selected
    assert len(disease_ratios) == 3 and len(control_ratios) == 3


def test_sphingomyelin_synthase_decomposed_to_cer_sm_edge():
    """Multi-substrate SM synthase yields a Cer->SM backbone class reaction.

    'PC; Ceramide -> 1,2-DG; Sphingomyelin' is a 2-in/2-out reaction that the
    1:1 extraction rule skips; the sphingoid backbone passes through unchanged,
    so a same-structure Cer->SM edge must be derived from it. dhCer->dhSM is
    derived likewise.
    """
    dataset = LipidDataset(samples=[SampleMetadata(sample_name="s1", group="g")], lipids=[])
    dataset.reactions = [
        ReactionData(
            reactants=[
                CompoundComponent(compound_name="PC", compound_headgroup="PC"),
                CompoundComponent(compound_name="Ceramide", compound_headgroup="Cer"),
            ],
            products=[
                CompoundComponent(compound_name="1,2-DG", compound_headgroup="DG"),
                CompoundComponent(compound_name="Sphingomyelin", compound_headgroup="SM"),
            ],
        ),
        ReactionData(
            reactants=[
                CompoundComponent(compound_name="PC", compound_headgroup="PC"),
                CompoundComponent(compound_name="Dihydroceramide", compound_headgroup="dhCer"),
            ],
            products=[
                CompoundComponent(compound_name="1,2-DG", compound_headgroup="DG"),
                CompoundComponent(compound_name="Dihydrosphingomyelin", compound_headgroup="dhSM"),
            ],
        ),
    ]
    exporter = BioPANPathwayExporter(dataset=dataset)
    class_reactions, _ = exporter._extract_class_reactions(dataset)
    keys = {(cr.reactant_class, cr.product_class) for cr in class_reactions}
    assert ("Cer", "SM") in keys
    assert ("dhCer", "dhSM") in keys
    # Backbone passes through unchanged -> same-structure (no compound required).
    cer_sm = next(cr for cr in class_reactions if (cr.reactant_class, cr.product_class) == ("Cer", "SM"))
    assert cer_sm.compound_require == CompoundRequirement.NONE


def test_multisubstrate_without_single_sphingo_pair_is_skipped():
    """A 2-in/2-out reaction with no clean single sphingo pairing adds no edge."""
    dataset = LipidDataset(samples=[SampleMetadata(sample_name="s1", group="g")], lipids=[])
    dataset.reactions = [
        ReactionData(
            reactants=[
                CompoundComponent(compound_name="PC", compound_headgroup="PC"),
                CompoundComponent(compound_name="MG", compound_headgroup="MG"),
            ],
            products=[
                CompoundComponent(compound_name="1,2-DG", compound_headgroup="DG"),
                CompoundComponent(compound_name="LPC", compound_headgroup="LPC"),
            ],
        ),
    ]
    exporter = BioPANPathwayExporter(dataset=dataset)
    class_reactions, _ = exporter._extract_class_reactions(dataset)
    assert class_reactions == []


def test_too_many_zero_reactants_unscorable():
    exporter = BioPANPathwayExporter()
    # 3 of 4 positions have a zero reactant -> l >= n - 1 -> unscorable.
    selected = exporter._select_ratio_vectors(
        disease_products=[2.0, 2.1, 2.2, 2.3],
        disease_reactants=[0.0, 0.0, 0.0, 10.0],
        control_products=[1.0, 1.1, 0.9, 1.0],
        control_reactants=[10.0, 10.0, 10.0, 10.0],
        paired=False,
    )
    assert selected is None
