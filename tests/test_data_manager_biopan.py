import json

from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models.reaction import CompoundComponent, ReactionData
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleConditions, SampleMetadata


def make_biopan_dataset() -> LipidDataset:
    samples = [
        SampleMetadata(sample_name="ctrl_1", group="control"),
        SampleMetadata(sample_name="ctrl_2", group="control"),
        SampleMetadata(sample_name="case_1", group="case"),
    ]
    lipids = [
        QuantifiedLipid(
            input_name="PC(34:1)",
            standardized_name="PC 34:1",
            generic_lm_id="LMGP01010000",
            reactions=[ReactionData(reaction_id=1, reaction_name="PC to LPC")],
            values={"ctrl_1": 10.0, "ctrl_2": 12.0, "case_1": 14.0},
        ),
        QuantifiedLipid(
            input_name="CE(18:1)",
            standardized_name="CE 18:1",
            generic_lm_id="LMST01020000",
            values={"ctrl_1": 4.0, "ctrl_2": 5.0, "case_1": 6.0},
        ),
        QuantifiedLipid(
            input_name="FA(18:1)",
            standardized_name="FA 18:1",
            generic_lm_id="LMFA01010000",
            reactions=[ReactionData(reaction_id=2, reaction_name="FA activation")],
            values={"ctrl_1": 2.0, "ctrl_2": 3.0, "case_1": 4.0},
        ),
        QuantifiedLipid(
            input_name="O-TG(50:0)",
            values={"ctrl_1": 1.0, "ctrl_2": 1.5, "case_1": 2.0},
        ),
    ]
    return LipidDataset(
        samples=samples,
        lipids=lipids,
        column_info={"empty_columns": []},
    )


def test_build_biopan_summary_classifies_species():
    manager = DataManager(dataset=make_biopan_dataset())

    summary = manager.build_biopan_summary()

    assert summary["total"] == 4
    assert summary["groups"] == {
        "ctrl_1": "control",
        "ctrl_2": "control",
        "case_1": "case",
    }
    assert summary["pathway"]["processed"] == ["PC(34:1)", "FA(18:1)"]
    assert summary["pathway"]["unprocessed"] == ["CE(18:1)"]
    assert summary["undef"] == ["O-TG(50:0)"]
    assert summary["processed_dataset"] == {
        "PC(34:1)": " PC 34:1",
        "FA(18:1)": " FA 18:1",
    }
    assert summary["unprocessed_dataset"] == {"CE(18:1)": " CE 18:1"}


def test_get_sample_conditions_returns_pydantic_mapping_for_auto_detected_groups():
    manager = DataManager()
    dataset = LipidDataset(
        samples=manager.extract_sample_metadata(["ctrl_1", "ctrl_2", "case_1"]),
        lipids=[],
    )

    conditions = manager.get_sample_conditions(dataset)

    assert isinstance(conditions, SampleConditions)
    assert conditions.model_dump() == {
        "ctrl_1": "ctrl",
        "ctrl_2": "ctrl",
        "case_1": "case",
    }


def test_set_sample_conditions_updates_user_selected_conditions():
    manager = DataManager(dataset=make_biopan_dataset())

    updated = manager.set_sample_conditions({"ctrl_1": "baseline", "case_1": "treated"})

    assert updated.model_dump() == {
        "ctrl_1": "baseline",
        "ctrl_2": "control",
        "case_1": "treated",
    }
    assert [sample.group for sample in manager.dataset.samples] == ["baseline", "control", "treated"]


def test_set_sample_conditions_rejects_unknown_samples():
    dataset = make_biopan_dataset()

    try:
        dataset.set_sample_conditions({"missing_sample": "treated"})
    except ValueError as exc:
        assert str(exc) == "Unknown sample names in condition mapping: missing_sample"
    else:
        raise AssertionError("Expected ValueError for unknown sample names")


def test_build_biopan_msg1_and_msg2_capture_display_state():
    manager = DataManager(dataset=make_biopan_dataset())

    msg1 = manager.build_biopan_msg1()
    msg2 = manager.build_biopan_msg2()

    assert msg1["na_values"] == [False]
    assert msg1["empty_rows"] == [False]
    assert msg1["empty_cols"] == [False]
    assert msg1["has_fa"] == ["facoa"]
    assert msg2["valid"] == {"groups": ["control"], "freqs": [2]}
    assert msg2["notvalid"] == {"groups": ["case"], "freqs": [1]}
    assert msg2["reaction"]["lp"] == [True]
    assert msg2["reaction"]["fa"] == [True]
    assert msg2["subset"]["reaction"] == [True]
    assert msg2["subset"]["pathway"] == [True]


def test_export_biopan_display_files_writes_session_biopan_folder(tmp_path):
    manager = DataManager(dataset=make_biopan_dataset())

    written = manager.export_biopan_display_files(tmp_path / "WSMhyRzAGt5IoYEJ")

    assert set(written) == {"msg1", "summary", "msg2"}
    summary_path = tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "summary.json"
    msg2_path = tmp_path / "WSMhyRzAGt5IoYEJ" / "biopan" / "msg2.json"
    assert summary_path.exists()
    assert msg2_path.exists()

    with summary_path.open(encoding="utf-8") as handle:
        summary = json.load(handle)
    with msg2_path.open(encoding="utf-8") as handle:
        msg2 = json.load(handle)

    assert summary["pathway"]["processed"] == ["PC(34:1)", "FA(18:1)"]
    assert msg2["valid"]["groups"] == ["control"]


def test_build_biopan_summary_uses_empty_object_for_empty_undef():
    dataset = make_biopan_dataset()
    dataset.lipids[-1].generic_lm_id = "LMGL03010000"
    manager = DataManager(dataset=dataset)

    summary = manager.build_biopan_summary(lipidlynxx="yes")

    assert summary["lipidlynxx"] == "yes"
    assert summary["undef"] == {}


def test_get_lipids_for_component_matches_compound_generic_lm_id():
    dataset = make_biopan_dataset()

    component = CompoundComponent(
        compound_type="lm_main",
        compound_generic_lm_id="LMGP01010000",
    )

    matched = dataset.get_lipids_for_component(component)

    assert [lipid.input_name for lipid in matched] == ["PC(34:1)"]