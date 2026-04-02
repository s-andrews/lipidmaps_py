from lipidmaps.biopan_cli import build_parser, _apply_group_overrides, _parse_sample_group_overrides, _resolve_groups, _resolve_input_path
from lipidmaps.data import BioPANExporter
from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata


def make_dataset() -> LipidDataset:
    return LipidDataset(
        samples=[
            SampleMetadata(sample_name="S1", group="alpha"),
            SampleMetadata(sample_name="S2", group="alpha"),
            SampleMetadata(sample_name="S3", group="beta"),
        ],
        lipids=[
            QuantifiedLipid(input_name="PC(34:1)", values={"S1": 1.0, "S2": 2.0, "S3": 3.0}),
        ],
    )


def test_parse_sample_group_overrides_accepts_multiple_entries():
    overrides = _parse_sample_group_overrides(["S1=control", "S3=treated"])

    assert overrides == {"S1": "control", "S3": "treated"}


def test_parse_sample_group_overrides_rejects_invalid_format():
    try:
        _parse_sample_group_overrides(["S1-control"])
    except ValueError as exc:
        assert "Expected SAMPLE=GROUP" in str(exc)
    else:
        raise AssertionError("Expected ValueError for malformed sample-group override")


def test_apply_group_overrides_updates_exported_groups():
    dataset = make_dataset()

    _apply_group_overrides(dataset, {"S2": "treated", "S3": "treated"})

    msg2 = BioPANExporter(dataset=dataset).build_msg2()

    assert [sample.group for sample in dataset.samples] == ["alpha", "treated", "treated"]
    assert msg2["valid"] == {"groups": ["treated"], "freqs": [2]}
    assert msg2["notvalid"] == {"groups": ["alpha"], "freqs": [1]}


def test_resolve_groups_uses_first_seen_group_order_by_default():
    manager = DataManager()
    manager.dataset = LipidDataset(
        samples=[
            SampleMetadata(sample_name="S1", group="dr_young"),
            SampleMetadata(sample_name="S2", group="dr_young"),
            SampleMetadata(sample_name="S3", group="al_young"),
            SampleMetadata(sample_name="S4", group="al_young"),
        ],
        lipids=[],
    )

    disease_group, control_group = _resolve_groups(manager, None, None)

    assert disease_group == "dr_young"
    assert control_group == "al_young"


def test_build_parser_supports_has_labels_flag():
    parser = build_parser()

    args = parser.parse_args(["/tmp/session", "--has-labels"])

    assert args.has_labels is True


def test_resolve_input_path_prefers_existing_tsv_file(tmp_path):
    session_dir = tmp_path / "session"
    input_dir = session_dir / "input"
    input_dir.mkdir(parents=True)
    tsv_path = input_dir / "input.tsv"
    tsv_path.write_text("lipid\tsample\nPC 34:1\t1\n")

    resolved = _resolve_input_path(session_dir, None)

    assert resolved == tsv_path