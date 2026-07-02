import json
from types import SimpleNamespace

from lipidmaps.data import BioPANPathwayExporter
from lipidmaps.data.models.reaction import CompoundComponent, ReactionData
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
from lipidmaps.data.models.species_reaction import ClassReaction


def load_comparison_payload(output_root, key):
    with (output_root / "biopan" / "comparison_bundle.json").open(encoding="utf-8") as handle:
        bundle = json.load(handle)
    return bundle["payloads"][key]


def load_edge_detail_payload(output_root, key):
    with (output_root / "biopan" / "edge_details_bundle.json").open(encoding="utf-8") as handle:
        bundle = json.load(handle)
    return bundle["payloads"][key]


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
    def test_pathway_exporter_reads_object_gene_entries_and_fallback():
        exporter = BioPANPathwayExporter()
        exporter._reaction_gene_lookup = {"411": ["HSD3B2", "HSD3B1", "CYP11A1"]}

        reaction = ReactionData(
            reaction_id=411,
            genes=[SimpleNamespace(gene_name="HSD3B2")],
            proteins=[SimpleNamespace(protein_gene_name="HSD3B1")],
        )

        assert exporter._get_reaction_gene_symbols(reaction) == ["HSD3B2", "HSD3B1", "CYP11A1"]
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
    assert "comparison_bundle.json" in written
    assert "edge_details_bundle.json" in written

    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_reaction.json").open(encoding="utf-8") as handle:
        tree = json.load(handle)
    graph = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_class_reaction_case_control_active_notpaired.json",
    )
    table = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_class_reaction_case_control_active_0.05_notpaired_tbl.json",
    )
    edge = load_edge_detail_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "pc,pa.json",
    )

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
    assert "comparison_bundle.json" in written

    with (tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "lp_class_pathway.json").open(encoding="utf-8") as handle:
        tree = json.load(handle)
    graph = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_class_pathway_case_control_active_notpaired.json",
    )
    highlight = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_class_pathway_case_control_active_0.05_notpaired.json",
    )
    table = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_class_pathway_case_control_active_0.05_notpaired_tbl.json",
    )

    assert tree == [{"text": "Glycerolipids and Glycerophospholipids", "children": [{"text": "PA"}, {"text": "PC"}]}]
    assert {node["data"]["label"] for node in graph["nodes"]} == {"PC", "PA"}
    assert graph["edges"][0]["data"]["id"] == "PCPA"
    assert highlight["nodes"] == "#pa,#pc"
    assert highlight["edges"] == "#PCPA"
    assert table["pathways"][0]["data"]["pathway"] == "PC&#8594;PA"
    assert table["pathways"][0]["data"]["class"] == "Phosphatidylcholine turnover (Glycerolipids and Glycerophospholipids)"


def make_fa_dataset() -> LipidDataset:
    samples = [
        SampleMetadata(sample_name="ctrl_1", group="control"),
        SampleMetadata(sample_name="ctrl_2", group="control"),
        SampleMetadata(sample_name="ctrl_3", group="control"),
        SampleMetadata(sample_name="case_1", group="case"),
        SampleMetadata(sample_name="case_2", group="case"),
        SampleMetadata(sample_name="case_3", group="case"),
    ]
    lipids = [
        # FA(16:0) consumed; FA(18:0) and FA(18:1) produced in the case group, so
        # the elongation/desaturation steps read as "active".
        QuantifiedLipid(input_name="FA(16:0)", values={"ctrl_1": 10.0, "ctrl_2": 11.0, "ctrl_3": 12.0, "case_1": 5.0, "case_2": 6.0, "case_3": 7.0}),
        QuantifiedLipid(input_name="FA(18:0)", values={"ctrl_1": 2.0, "ctrl_2": 2.5, "ctrl_3": 3.0, "case_1": 10.0, "case_2": 11.0, "case_3": 12.0}),
        # Space-style name (no parentheses) to exercise FA name canonicalisation.
        QuantifiedLipid(input_name="FA 18:1", values={"ctrl_1": 1.0, "ctrl_2": 1.5, "ctrl_3": 2.0, "case_1": 9.0, "case_2": 10.0, "case_3": 11.0}),
    ]
    return LipidDataset(samples=samples, lipids=lipids)


def test_fa_graph_and_table_build_for_measured_fatty_acids():
    exporter = BioPANPathwayExporter(dataset=make_fa_dataset())

    graph = exporter.build_fa_graph(disease_group="case", control_group="control", mode="active")
    node_labels = {node["data"]["label"] for node in graph["nodes"]}
    # FA species use shorthand nomenclature (FA 16:0, not FA(16:0)).
    assert {"FA 16:0", "FA 18:0", "FA 18:1"} <= node_labels
    assert not any("(" in label for label in node_labels)
    # FA nodes render as triangles (matches the legend).
    assert all(node["data"]["shape"] == "triangle" for node in graph["nodes"])
    # The elongation/desaturation edges are present and scored.
    edge_pairs = {(e["data"]["source"], e["data"]["target"]) for e in graph["edges"]}
    assert (exporter._safe_node_id("FA(16:0)"), exporter._safe_node_id("FA(18:0)")) in edge_pairs

    table = exporter.build_fa_table(disease_group="case", control_group="control", threshold=0.05, mode="active")
    rows = table["pathways"]
    assert rows, "expected significant FA reactions"
    top = rows[0]
    assert "&#8594;" in top["data"]["pathway"]
    assert "(" not in top["data"]["pathway"]  # shorthand FA names
    assert top["data"]["score"] > 1.645
    # Enzyme genes from the curated FA network flow through to the table.
    assert top["data"]["gene"] != "NA"
    # FA enzyme genes are tagged human for the frontend's taxonomy flag.
    assert "9606" in table["gene_taxonomy"].get(top["data"]["gene"].split(",")[0], [])


def test_export_reaction_files_includes_fa_payloads(tmp_path):
    exporter = BioPANPathwayExporter(dataset=make_fa_dataset())
    exporter.export_reaction_files(
        tmp_path / "FAsession00000001",
        disease_group="case",
        control_group="control",
        threshold=0.05,
        paired=False,
        lazy=False,
    )
    graph_payload = load_comparison_payload(tmp_path / "FAsession00000001", "fa_case_control_active_notpaired.json")
    assert graph_payload["edges"], "FA graph payload should contain edges"
    table_payload = load_comparison_payload(tmp_path / "FAsession00000001", "fa_case_control_active_0.05_notpaired_tbl.json")
    assert "pathways" in table_payload


def test_fa_view_does_not_alter_lipid_scoring():
    # Building the FA payloads must not change the lipid reaction scores.
    baseline = BioPANPathwayExporter(dataset=make_pathway_dataset()).build_reaction_table(
        disease_group="case", control_group="control", threshold=0.05, level="class", mode="active",
    )["pathways"][0]["data"]["score"]

    exporter = BioPANPathwayExporter(dataset=make_pathway_dataset())
    exporter._fa_payload_items("case", "control", 0.05, False, make_pathway_dataset())
    after = exporter.build_reaction_table(
        disease_group="case", control_group="control", threshold=0.05, level="class", mode="active",
    )["pathways"][0]["data"]["score"]

    assert baseline == after


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

    assert "comparison_bundle.json" in written

    table = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_species_pathway_case_control_active_0.05_notpaired_tbl.json",
    )
    highlight = load_comparison_payload(
        tmp_path / "WSMhyRzAGt5IoYEJ",
        "lp_species_pathway_case_control_active_0.05_notpaired.json",
    )

    assert table["pathways"][0]["data"]["pathway"] == "DG(36:0)&#8594;PE(36:0)&#8594;PC(36:0)"
    assert table["pathways"][0]["data"]["gene"] == "CEPT1,PEMT,PEMT2"
    assert table["pathways"][0]["data"]["class"] == "Biosynthesis of PC (Glycerolipids and Glycerophospholipids)"
    assert table["pathways"][0]["data"]["score"] > 1.645
    assert highlight["edges"] == "#DG360PE360,#PE360PC360"


def test_pathway_exporter_reads_gene_name_fields_and_reaction_gene_fallback():
    exporter = BioPANPathwayExporter()
    exporter._reaction_gene_lookup = {"411": ["HSD3B2", "HSD3B1"]}

    reaction = ReactionData(
        reaction_id=411,
        genes=[{"gene_name": "HSD3B2"}],
        proteins=[{"protein_gene_name": "HSD3B1"}],
    )

    assert exporter._get_reaction_gene_symbols(reaction) == ["HSD3B2", "HSD3B1"]


def test_pathway_exporter_uses_reaction_gene_lookup_when_payload_omits_gene_names():
    exporter = BioPANPathwayExporter()
    exporter._reaction_gene_lookup = {"411": ["HSD3B2", "HSD3B1"]}

    reaction = ReactionData(reaction_id=411)

    assert exporter._get_reaction_gene_symbols(reaction) == ["HSD3B2", "HSD3B1"]


def test_pathway_exporter_prefers_api_genes_over_uniprot_lookup():
    # The reactions API now serves genes from a ready table, so API-provided
    # genes are used directly and the live UniProt lookup is not consulted.
    exporter = BioPANPathwayExporter()
    exporter._uniprot_gene_lookup = {"10273": ["SHOULD_NOT_BE_USED"]}

    reaction = ReactionData(
        genes=[
            {"gene_name": "LCAT", "uniprot_id": "P04180", "taxonomy_id": 9606},
            {"gene_name": "PLA2G15", "uniprot_id": "Q8NCC3", "taxonomy_id": 10090},
        ],
        curations=[{"database_name": "rhea", "database_id": "10273"}],
    )

    assert exporter._get_reaction_gene_symbols(reaction) == ["LCAT", "PLA2G15"]
    # The gene -> UniProt accession map is populated for the frontend links.
    assert exporter._gene_uniprot_lookup["LCAT"] == "P04180"
    assert exporter._gene_uniprot_lookup["PLA2G15"] == "Q8NCC3"
    # The gene -> taxonomy map is populated for the human flag / hover.
    assert exporter._gene_taxonomy_lookup["LCAT"] == ["9606"]
    assert exporter._gene_taxonomy_lookup["PLA2G15"] == ["10090"]


def test_pathway_exporter_records_gene_organism_names_from_api_genes():
    # The reactions API serves the organism scientific/common name alongside each
    # gene; the exporter records a taxonomy name (scientific_name preferred, then
    # common_name) so the frontend can display the organism on hover.
    exporter = BioPANPathwayExporter()
    reaction = ReactionData(
        genes=[
            {
                "taxonomy_id": 9606,
                "uniprot_id": "P48449",
                "gene_name": "LSS",
                "scientific_name": "Homo sapiens",
                "common_name": "human",
            },
            {
                "taxonomy_id": 559292,
                "uniprot_id": "P38604",
                "gene_name": "ERG7",
                "scientific_name": "Saccharomyces cerevisiae S288C",
                "common_name": "",
            },
            {
                # No scientific_name -> falls back to common_name.
                "taxonomy_id": 10090,
                "uniprot_id": "Q99N16",
                "gene_name": "Lss",
                "common_name": "mouse",
            },
        ],
    )

    assert exporter._get_reaction_gene_symbols(reaction) == ["LSS", "ERG7", "Lss"]
    assert exporter._gene_taxonomy_name_lookup["LSS"] == ["Homo sapiens"]
    assert exporter._gene_taxonomy_name_lookup["ERG7"] == ["Saccharomyces cerevisiae S288C"]
    assert exporter._gene_taxonomy_name_lookup["Lss"] == ["mouse"]


def test_pathway_exporter_falls_back_to_uniprot_lookup_when_no_api_genes():
    # When the API has no genes for a reaction yet, fall back to the live
    # UniProt lookup keyed by rhea_id.
    exporter = BioPANPathwayExporter()
    exporter._uniprot_gene_lookup = {"10273": ["LCAT", "PLA2G15"]}

    reaction = ReactionData(
        curations=[{"database_name": "rhea", "database_id": "10273"}],
    )

    assert exporter._get_reaction_gene_symbols(reaction) == ["LCAT", "PLA2G15"]


def test_reaction_data_lists_normalized_rhea_ids():
    reaction = ReactionData(
        curations=[
            {"database_name": "rhea", "database_id": "RHEA:10273"},
            {"database_name": "kegg", "database_id": "C00001"},
            {"database_name": "rhea", "database_id": 15421},
            {"database_name": "rhea", "database_id": " RHEA:20311 "},
        ]
    )

    assert reaction.list_rhea_ids() == ["10273", "15421", "20311"]


def test_reaction_data_preserves_explicit_rhea_id():
    reaction = ReactionData(
        rhea_id="RHEA:18766",
        curations=[{"database_name": "rhea", "database_id": "10273"}],
    )

    assert reaction.list_rhea_ids() == ["18766"]
