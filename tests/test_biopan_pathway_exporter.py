import json

from lipidmaps.data import BioPANPathwayExporter
from lipidmaps.data.models.reaction import CompoundComponent, ReactionData
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
from lipidmaps.data.models.species_reaction import ClassReaction


def make_reaction_dataset() -> LipidDataset:
    samples = [
        SampleMetadata(sample_name="ctrl_1", group="control"),
        SampleMetadata(sample_name="ctrl_2", group="control"),
        SampleMetadata(sample_name="ctrl_3", group="control"),
        SampleMetadata(sample_name="case_1", group="case"),
        SampleMetadata(sample_name="case_2", group="case"),
        SampleMetadata(sample_name="case_3", group="case"),
    ]
    lipids = [
        QuantifiedLipid(input_name="PC(34:1)", values={"ctrl_1": 12.0, "ctrl_2": 11.0, "ctrl_3": 10.0, "case_1": 6.0, "case_2": 6.5, "case_3": 7.0}),
        QuantifiedLipid(input_name="PA(34:1)", values={"ctrl_1": 3.0, "ctrl_2": 2.5, "ctrl_3": 2.0, "case_1": 12.0, "case_2": 11.0, "case_3": 10.0}),
        QuantifiedLipid(input_name="PC(36:2)", values={"ctrl_1": 11.0, "ctrl_2": 10.5, "ctrl_3": 10.0, "case_1": 5.5, "case_2": 6.0, "case_3": 6.5}),
        QuantifiedLipid(input_name="PA(36:2)", values={"ctrl_1": 2.5, "ctrl_2": 2.0, "ctrl_3": 1.5, "case_1": 11.0, "case_2": 10.0, "case_3": 9.0}),
    ]
    return LipidDataset(samples=samples, lipids=lipids)


def make_pathway_dataset() -> LipidDataset:
    dataset = make_reaction_dataset()
    dataset.reactions = [
        ReactionData(
            reactants=[CompoundComponent(compound_name="PC", compound_headgroup="PC")],
            products=[CompoundComponent(compound_name="PA", compound_headgroup="PA")],
            genes=[{"gene_symbol": "PLA2G15"}],
            pathways=[
                {
                    "pathway_name": "Phosphatidylcholine turnover",
                    "pathway_type": ["Glycerolipids and Glycerophospholipids"],
                }
            ],
        )
    ]
    return dataset


def make_suppressed_pathway_dataset() -> LipidDataset:
    dataset = LipidDataset(
        samples=[
            SampleMetadata(sample_name="ctrl_1", group="control"),
            SampleMetadata(sample_name="ctrl_2", group="control"),
            SampleMetadata(sample_name="ctrl_3", group="control"),
            SampleMetadata(sample_name="case_1", group="case"),
            SampleMetadata(sample_name="case_2", group="case"),
            SampleMetadata(sample_name="case_3", group="case"),
        ],
        lipids=[
            QuantifiedLipid(input_name="PC(34:1)", values={"ctrl_1": 6.0, "ctrl_2": 6.5, "ctrl_3": 7.0, "case_1": 12.0, "case_2": 11.0, "case_3": 10.0}),
            QuantifiedLipid(input_name="PA(34:1)", values={"ctrl_1": 12.0, "ctrl_2": 11.0, "ctrl_3": 10.0, "case_1": 3.0, "case_2": 2.5, "case_3": 2.0}),
            QuantifiedLipid(input_name="PC(36:2)", values={"ctrl_1": 5.5, "ctrl_2": 6.0, "ctrl_3": 6.5, "case_1": 11.0, "case_2": 10.0, "case_3": 9.0}),
            QuantifiedLipid(input_name="PA(36:2)", values={"ctrl_1": 11.0, "ctrl_2": 10.0, "ctrl_3": 9.0, "case_1": 2.5, "case_2": 2.0, "case_3": 1.5}),
        ],
    )
    dataset.reactions = [
        ReactionData(
            reactants=[CompoundComponent(compound_name="PC", compound_headgroup="PC")],
            products=[CompoundComponent(compound_name="PA", compound_headgroup="PA")],
            genes=[{"gene_symbol": "PLA2G15"}],
            pathways=[
                {
                    "pathway_name": "Phosphatidylcholine turnover",
                    "pathway_type": ["Glycerolipids and Glycerophospholipids"],
                }
            ],
        )
    ]
    return dataset


def make_multistep_pathway_dataset() -> LipidDataset:
    samples = [
        SampleMetadata(sample_name="ctrl_1", group="control"),
        SampleMetadata(sample_name="ctrl_2", group="control"),
        SampleMetadata(sample_name="ctrl_3", group="control"),
        SampleMetadata(sample_name="case_1", group="case"),
        SampleMetadata(sample_name="case_2", group="case"),
        SampleMetadata(sample_name="case_3", group="case"),
    ]
    lipids = [
        QuantifiedLipid(input_name="DG(36:0)", values={"ctrl_1": 10.0, "ctrl_2": 10.0, "ctrl_3": 10.0, "case_1": 5.0, "case_2": 5.0, "case_3": 5.0}),
        QuantifiedLipid(input_name="PE(36:0)", values={"ctrl_1": 3.0, "ctrl_2": 3.0, "ctrl_3": 3.0, "case_1": 10.0, "case_2": 10.0, "case_3": 10.0}),
        QuantifiedLipid(input_name="PC(36:0)", values={"ctrl_1": 1.0, "ctrl_2": 1.0, "ctrl_3": 1.0, "case_1": 9.0, "case_2": 9.0, "case_3": 9.0}),
    ]
    dataset = LipidDataset(samples=samples, lipids=lipids)
    dataset.reactions = [
        ReactionData(
            reactants=[CompoundComponent(compound_name="DG", compound_headgroup="DG")],
            products=[CompoundComponent(compound_name="PE", compound_headgroup="PE")],
            genes=[{"gene_symbol": "CEPT1"}],
            pathways=[
                {
                    "pathway_name": "Biosynthesis of PC",
                    "pathway_type": ["Glycerolipids and Glycerophospholipids"],
                }
            ],
        ),
        ReactionData(
            reactants=[CompoundComponent(compound_name="PE", compound_headgroup="PE")],
            products=[CompoundComponent(compound_name="PC", compound_headgroup="PC")],
            genes=[{"gene_symbol": "PEMT"}],
            proteins=[{"predicted_gene_symbol": "PEMT2"}],
            pathways=[
                {
                    "pathway_name": "Biosynthesis of PC",
                    "pathway_type": ["Glycerolipids and Glycerophospholipids"],
                }
            ],
        ),
    ]
    return dataset


def test_reaction_exporter_builds_initial_reaction_assets(tmp_path):
    exporter = BioPANPathwayExporter(
        dataset=make_reaction_dataset(),
        class_reactions=[ClassReaction(reactant_class="PC", product_class="PA", genes=["PLA2G15"])],
    )

    written = exporter.export_reaction_files(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        disease_group="case",
        control_group="control",
        threshold=0.05,
    )

    assert "lp_class_reaction.json" in written
    assert "lp_species_reaction_case_control_active_notpaired.json" in written
    assert "pc,pa.json" in written

    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_reaction.json").open(encoding="utf-8") as handle:
        tree = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_reaction_case_control_active_notpaired.json").open(encoding="utf-8") as handle:
        graph = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_reaction_case_control_active_0.05_notpaired_tbl.json").open(encoding="utf-8") as handle:
        table = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "pc,pa.json").open(encoding="utf-8") as handle:
        edge = json.load(handle)

    assert tree == [{"text": "Matched reactions", "children": [{"text": "PA"}, {"text": "PC"}]}]
    assert {node["data"]["label"] for node in graph["nodes"]} == {"PC", "PA"}
    assert len(graph["edges"]) == 1
    assert graph["edges"][0]["data"]["weight"] > 1.645
    assert table["pathways"][0]["data"]["pathway"] == "PC&#8594;PA"
    assert table["pathways"][0]["data"]["gene"] == "PLA2G15"
    assert edge["selected_re"] == "PC(34:1),PC(36:2)"
    assert edge["selected_pro"] == "PA(34:1),PA(36:2)"


def test_reaction_exporter_builds_species_tree_and_highlight_payload():
    exporter = BioPANPathwayExporter(
        dataset=make_reaction_dataset(),
        class_reactions=[ClassReaction(reactant_class="PC", product_class="PA", genes=["PLA2G15"])],
    )

    result_set, _ = exporter.build_reaction_match_set()
    species_tree = exporter.build_reaction_tree(level="species", result_set=result_set)
    highlight = exporter.build_reaction_highlight(
        disease_group="case",
        control_group="control",
        threshold=0.05,
        level="class",
        mode="active",
        result_set=result_set,
        reaction_lookup={},
    )

    assert species_tree[0]["children"][0]["text"] == "PA"
    assert species_tree[0]["children"][1]["text"] == "PC"
    assert "#pc" in highlight["nodes"]
    assert "#pa" in highlight["nodes"]
    assert "#PCPA" in highlight["edges"]


def test_reaction_exporter_builds_pathway_assets(tmp_path):
    exporter = BioPANPathwayExporter(dataset=make_pathway_dataset())

    written = exporter.export_reaction_files(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        disease_group="case",
        control_group="control",
        threshold=0.05,
    )

    assert "lp_class_pathway.json" in written
    assert "lp_species_pathway_case_control_active_notpaired.json" in written
    assert "lp_class_pathway_case_control_active_0.05_notpaired_tbl.json" in written

    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_pathway.json").open(encoding="utf-8") as handle:
        tree = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_pathway_case_control_active_notpaired.json").open(encoding="utf-8") as handle:
        graph = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_pathway_case_control_active_0.05_notpaired.json").open(encoding="utf-8") as handle:
        highlight = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_pathway_case_control_active_0.05_notpaired_tbl.json").open(encoding="utf-8") as handle:
        table = json.load(handle)

    assert tree == [{"text": "Glycerolipids and Glycerophospholipids", "children": [{"text": "PA"}, {"text": "PC"}]}]
    assert {node["data"]["label"] for node in graph["nodes"]} == {"PC", "PA"}
    assert graph["edges"][0]["data"]["id"] == "PCPA"
    assert highlight["nodes"] == "#pa,#pc"
    assert highlight["edges"] == "#PCPA"
    assert table["pathways"][0]["data"]["pathway"] == "PC&#8594;PA"
    assert table["pathways"][0]["data"]["class"] == "Phosphatidylcholine turnover (Glycerolipids and Glycerophospholipids)"


def test_pathway_exporter_uses_legacy_one_sided_scores_for_active_and_suppressed_modes():
    active_exporter = BioPANPathwayExporter(dataset=make_pathway_dataset())
    active_reaction_rows = active_exporter.build_reaction_table(
        disease_group="case",
        control_group="control",
        threshold=0.05,
        level="class",
        mode="active",
    )["pathways"]
    active_pathway_rows = active_exporter.build_pathway_table(
        disease_group="case",
        control_group="control",
        threshold=0.05,
        level="class",
        mode="active",
    )["pathways"]

    assert active_reaction_rows[0]["data"]["score"] > 1.645
    assert active_pathway_rows[0]["data"]["score"] > 1.645

    suppressed_exporter = BioPANPathwayExporter(dataset=make_suppressed_pathway_dataset())
    suppressed_reaction_rows = suppressed_exporter.build_reaction_table(
        disease_group="case",
        control_group="control",
        threshold=0.05,
        level="class",
        mode="suppressed",
    )["pathways"]
    suppressed_pathway_rows = suppressed_exporter.build_pathway_table(
        disease_group="case",
        control_group="control",
        threshold=0.05,
        level="class",
        mode="suppressed",
    )["pathways"]

    assert suppressed_reaction_rows[0]["data"]["score"] > 1.645
    assert suppressed_pathway_rows[0]["data"]["score"] > 1.645


def test_pathway_exporter_builds_multistep_chain_and_gene_aggregation(tmp_path):
    exporter = BioPANPathwayExporter(dataset=make_multistep_pathway_dataset())

    written = exporter.export_reaction_files(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        disease_group="case",
        control_group="control",
        threshold=0.05,
    )

    assert "lp_species_pathway_case_control_active_0.05_notpaired_tbl.json" in written

    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_species_pathway_case_control_active_0.05_notpaired_tbl.json").open(encoding="utf-8") as handle:
        table = json.load(handle)
    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_species_pathway_case_control_active_0.05_notpaired.json").open(encoding="utf-8") as handle:
        highlight = json.load(handle)

    assert table["pathways"][0]["data"]["pathway"] == "DG(36:0)&#8594;PE(36:0)&#8594;PC(36:0)"
    assert table["pathways"][0]["data"]["gene"] == "CEPT1,PEMT,PEMT2"
    assert table["pathways"][0]["data"]["class"] == "Biosynthesis of PC (Glycerolipids and Glycerophospholipids)"
    assert table["pathways"][0]["data"]["score"] > 1.645
    assert highlight["edges"] == "#DG360PE360,#PE360PC360"