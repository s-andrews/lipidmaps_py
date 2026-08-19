"""Regression tests for BioPAN FA/CoA pool construction and ratio scoring.

These lock in two behaviours:

1. The FA / acyl-CoA pools are the *union* of reaction-extracted, dataset-named
   and dataset-inferred species, so a sparse reaction-derived set (e.g. only
   FA 16:0) can no longer short-circuit the richer set inferred from the
   dataset's full-structure acyl chains. Mono-acyl (lyso) self-inference is
   gated on the dataset having no measured FA/FaCoA, mirroring legacy BioPAN
   (parse_data.r): lyso/CE/sphingo chains seed the pool only when no fatty
   acids are measured. When FAs *are* measured, an LPS species must not seed
   the FA that then validates its own PS->LPS release (which would over-match
   release pairs and inflate the subclass z-score).

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


def test_lyso_seeds_fa_pool_when_no_measured_fa():
    """Legacy BioPAN synthesises an FA pool from lyso/CE/sphingo chains when the
    dataset has no measured FA. _ps_lps_dataset carries no FA species, so the
    LPS chains DO seed the pool."""
    dataset = _ps_lps_dataset()
    exporter = BioPANPathwayExporter(dataset=dataset)

    _, fa_names, _ = exporter._get_matching_inputs(dataset)

    for seeded in ("FA 18:0", "FA 18:1", "FA 20:4", "FA 22:6"):
        assert seeded in fa_names, f"{seeded} should be seeded from LPS when no FA is measured"


def test_lyso_excluded_from_fa_pool_when_fa_measured():
    """Once any FA is measured, legacy BioPAN stops synthesising from lyso
    chains, so an LPS product no longer seeds the FA that validates its own
    release. The measured FA stays; the lyso-only chains drop out."""
    dataset = _ps_lps_dataset()
    # Add a measured FA species (present across all samples).
    values = {s.sample_name: 1.0 for s in dataset.samples}
    dataset.lipids.append(QuantifiedLipid(input_name="FA(16:0)", values=values))
    exporter = BioPANPathwayExporter(dataset=dataset)

    _, fa_names, _ = exporter._get_matching_inputs(dataset)

    assert "FA 16:0" in fa_names  # measured FA retained
    for excluded in ("FA 18:1", "FA 20:4", "FA 22:6"):
        assert excluded not in fa_names, f"{excluded} should not be lyso-seeded once FA is measured"


def test_union_still_includes_full_structure_inferred_fa():
    """The union must still surface FAs inferred from full-structure species,
    independent of the mono-acyl gate. With a measured FA present (so mono-acyl
    inference is off), FA 18:1 from a full-structure PC is still pooled."""
    samples = [SampleMetadata(sample_name="s1", group="g")]
    lipids = [
        QuantifiedLipid(input_name="PC 16:0_18:1", values={"s1": 1.0}),
        QuantifiedLipid(input_name="FA(16:0)", values={"s1": 1.0}),  # measured FA -> gate off
    ]
    dataset = LipidDataset(samples=samples, lipids=lipids)
    exporter = BioPANPathwayExporter(dataset=dataset)

    _, fa_names, _ = exporter._get_matching_inputs(dataset)

    assert "FA 16:0" in fa_names  # measured
    assert "FA 18:1" in fa_names  # inferred from the full-structure PC despite gate off


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


def test_legacy_substrate_consumption_drops_orphan_reactant():
    """Legacy greedy consumption (get_reaction_fa_coa) assigns each product to a
    single reactant in ascending (carbons, double_bonds) order and drops
    reactants left with none. Default many-to-many keeps both."""
    from lipidmaps.data.models.species_reaction import (
        PathwayReactionSet,
        ReactionMatchResult,
        SpeciesReactionPair,
    )
    from lipidmaps.data.utils.chain_parser import ChainParser

    parser = ChainParser()
    cr = ClassReaction(reactant_class="LPC", product_class="PC")

    def pair(reactant: str, product: str) -> SpeciesReactionPair:
        return SpeciesReactionPair(
            reactant=parser.parse(reactant), product=parser.parse(product), class_reaction=cr
        )

    # Both LPC 16:0 and LPC 18:0 can form the single product PC 34:0.
    def fresh_pairs():
        return [pair("LPC 16:0", "PC 34:0"), pair("LPC 18:0", "PC 34:0")]

    exporter = BioPANPathwayExporter(
        dataset=LipidDataset(samples=[SampleMetadata(sample_name="s", group="g")], lipids=[]),
        legacy_substrate_consumption=True,
    )

    result = ReactionMatchResult(class_reaction=cr, pairs=fresh_pairs())
    rs = PathwayReactionSet(results={"lpc,pc": result})
    exporter._apply_legacy_substrate_consumption(rs)
    surviving = {(p.reactant.total_carbons, p.reactant.total_double_bonds) for p in result.pairs}
    # LPC 16:0 (sorted first) claims PC 34:0; LPC 18:0 is left with no product.
    assert surviving == {(16, 0)}

    # With two products, both reactants survive (one product each).
    result2 = ReactionMatchResult(
        class_reaction=cr,
        pairs=[pair("LPC 16:0", "PC 34:0"), pair("LPC 18:0", "PC 34:0"), pair("LPC 18:0", "PC 36:0")],
    )
    rs2 = PathwayReactionSet(results={"lpc,pc": result2})
    exporter._apply_legacy_substrate_consumption(rs2)
    surviving2 = {(p.reactant.total_carbons, p.reactant.total_double_bonds) for p in result2.pairs}
    assert surviving2 == {(16, 0), (18, 0)}


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
