import logging
from typing import Any

import pytest

from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models.sample import QuantifiedLipid, SampleMetadata, LipidDataset
from lipidmaps.data.models.refmet import RefMetResult


def make_dataset(names):
    samples = [SampleMetadata(sample_id="S1", group="g1")]
    lipids = [QuantifiedLipid(input_name=n, values={"S1": 1.0}) for n in names]
    return LipidDataset(samples=samples, lipids=lipids)


def test_run_lmsd_fill_and_report_updates(monkeypatch, caplog):
    caplog.set_level(logging.INFO)

    mgr = DataManager()
    ds = make_dataset(["A", "B"])

    # Ensure starting state has no lm_ids
    assert ds.lipids[0].lm_id is None
    assert ds.lipids[1].lm_id is None

    # Mock LMSD.get_lm_ids_by_name to return two lm_id entries
    def fake_get(names):
        assert names == ["A", "B"]
        return [
            {"input_name": "A", "lm_id": "LM_A", "matched_field": "name"},
            {"input_name": "B", "lm_id": "LM_B", "matched_field": "abbrev"},
        ]

    from lipidmaps.data.models import lmsd

    monkeypatch.setattr(lmsd.LMSD, "get_lm_ids_by_name", staticmethod(fake_get))

    updated = mgr.run_lmsd_fill_and_report(ds)

    assert updated == 2
    assert ds.lipids[0].lm_id == "LM_A"
    assert ds.lipids[0].lm_id_found_by == "LMSD"
    assert ds.lipids[0].matched_field == "name"
    assert ds.lipids[1].lm_id == "LM_B"
    assert ds.lipids[1].lm_id_found_by == "LMSD"


def test_annotate_lipids_with_refmet_sets_provenance(monkeypatch):
    mgr = DataManager()
    ds = make_dataset(["X", "Y"])

    # Create RefMetResult objects to return
    r1 = RefMetResult(input_name="X", standardized_name="StdX", lm_id="LMX")
    r2 = RefMetResult(input_name="Y", standardized_name=None, lm_id=None)

    from lipidmaps.data.models import refmet

    monkeypatch.setattr(refmet.RefMet, "validate_metabolite_names", staticmethod(lambda names: [r1, r2]))

    mgr.annotate_lipids_with_refmet(ds.lipids)

    assert ds.lipids[0].standardized_name == "StdX"
    assert ds.lipids[0].standardized_by == "RefMet"
    assert ds.lipids[0].lm_id == "LMX"
    assert ds.lipids[0].lm_id_found_by == "RefMet"

    # second lipid should not have provenance set
    assert ds.lipids[1].standardized_name is None
    assert ds.lipids[1].standardized_by is None
    assert ds.lipids[1].lm_id is None
    assert ds.lipids[1].lm_id_found_by is None
