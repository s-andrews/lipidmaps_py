from lipidmaps.biopan_cli import _apply_group_overrides, _parse_sample_group_overrides
from lipidmaps.data import BioPANExporter
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