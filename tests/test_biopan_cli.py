from lipidmaps import biopan_cli
from lipidmaps.biopan_cli import build_parser, _apply_group_overrides, _load_cached_dataset, _parse_sample_group_overrides, _resolve_groups, _resolve_input_path, _write_cached_dataset
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


def test_resolve_groups_defaults_match_legacy_order():
    # Legacy BioPAN / frontend default: control = first group, condition of
    # interest (disease) = second group.
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

    assert disease_group == "al_young"
    assert control_group == "dr_young"


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


def test_cached_dataset_round_trip(tmp_path):
    session_dir = tmp_path / "session"
    dataset = make_dataset()

    cache_path = _write_cached_dataset(session_dir, dataset)
    loaded = _load_cached_dataset(session_dir)

    assert cache_path == session_dir / "config" / "processed_dataset.json"
    assert loaded is not None
    assert [sample.sample_name for sample in loaded.samples] == ["S1", "S2", "S3"]
    assert loaded.lipids[0].input_name == "PC(34:1)"


def test_main_uses_cached_dataset_without_reprocessing_csv(tmp_path, monkeypatch):
    session_dir = tmp_path / "session"
    input_dir = session_dir / "input"
    input_dir.mkdir(parents=True)
    (input_dir / "input.csv").write_text("lipid,S1,S2,S3,S4\nPC(34:1),1,2,3,4\n")

    cached_dataset = LipidDataset(
        samples=[
            SampleMetadata(sample_name="S1", group="alpha"),
            SampleMetadata(sample_name="S2", group="alpha"),
            SampleMetadata(sample_name="S3", group="beta"),
            SampleMetadata(sample_name="S4", group="beta"),
        ],
        lipids=[
            QuantifiedLipid(input_name="PC(34:1)", values={"S1": 1.0, "S2": 2.0, "S3": 3.0, "S4": 4.0}),
        ],
    )
    _write_cached_dataset(session_dir, cached_dataset)

    exported = {}

    class DummyManager:
        def __init__(self, *args, **kwargs):
            self.dataset = None

        def process_csv(self, csv_path):
            raise AssertionError("process_csv should not be called when a cache is available")

        def export_biopan_display_files(self, output_path, dataset=None, **kwargs):
            exported["display"] = dataset or self.dataset
            return {}

        def export_biopan_reaction_files(self, output_path, disease_group, control_group, dataset=None, **kwargs):
            exported["reaction"] = {
                "dataset": dataset or self.dataset,
                "disease_group": disease_group,
                "control_group": control_group,
            }
            return {}

    monkeypatch.setattr(biopan_cli, "DataManager", DummyManager)
    monkeypatch.setattr(
        biopan_cli,
        "configure_logging",
        lambda: tmp_path / "logs",
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "biopan_cli",
            str(session_dir),
            "--sample-group",
            "S1=control",
            "--sample-group",
            "S2=control",
            "--sample-group",
            "S3=treated",
            "--sample-group",
            "S4=treated",
        ],
    )

    biopan_cli.main()

    assert [sample.group for sample in exported["display"].samples] == ["control", "control", "treated", "treated"]
    # Legacy default order: control = first group, disease = second group.
    assert exported["reaction"]["disease_group"] == "treated"
    assert exported["reaction"]["control_group"] == "control"