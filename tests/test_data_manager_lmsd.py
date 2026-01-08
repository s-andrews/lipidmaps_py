import pytest

from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.models.sample import QuantifiedLipid, LipidDataset
from lipidmaps.data.models.lmsd import LMSD


def test_fill_missing_lm_ids_uses_standardized_name(monkeypatch):
    q = QuantifiedLipid(input_name='Butyrylcarnitine', standardized_name='Butyrylcarnitine', values={})
    mgr = DataManager()
    mgr.dataset = LipidDataset(samples=[], lipids=[q])

    def fake_get(names):
        assert names == ['Butyrylcarnitine']
        return [{'input_name': 'Butyrylcarnitine', 'lm_id': 'LMFA07070054', 'matched_field': 'name'}]

    monkeypatch.setattr(LMSD, 'get_lm_ids_by_name', fake_get)

    updated = mgr.fill_missing_lm_ids_from_lmsd()
    assert updated == 1
    assert mgr.dataset.lipids[0].lm_id == 'LMFA07070054'


def test_fill_missing_lm_ids_falls_back_to_input_name(monkeypatch):
    q = QuantifiedLipid(input_name='Cholesterol', standardized_name=None, values={})
    mgr = DataManager()
    mgr.dataset = LipidDataset(samples=[], lipids=[q])

    def fake_get(names):
        assert names == ['Cholesterol']
        return [{'input_name': 'Cholesterol', 'lm_id': 'LMST01010001'}]

    monkeypatch.setattr(LMSD, 'get_lm_ids_by_name', fake_get)

    updated = mgr.fill_missing_lm_ids_from_lmsd()
    assert updated == 1
    assert mgr.dataset.lipids[0].lm_id == 'LMST01010001'


def test_fill_missing_lm_ids_does_not_overwrite_existing(monkeypatch):
    q1 = QuantifiedLipid(input_name='HasID', standardized_name='HasID', lm_id='EXIST', values={})
    q2 = QuantifiedLipid(input_name='NoID', standardized_name='NoID', values={})
    mgr = DataManager()
    mgr.dataset = LipidDataset(samples=[], lipids=[q1, q2])

    def fake_get(names):
        # only queried for the missing one
        assert names == ['NoID']
        return [{'input_name': 'NoID', 'lm_id': 'LMNO0001'}]

    monkeypatch.setattr(LMSD, 'get_lm_ids_by_name', fake_get)

    updated = mgr.fill_missing_lm_ids_from_lmsd()
    assert updated == 1
    assert mgr.dataset.lipids[0].lm_id == 'EXIST'
    assert mgr.dataset.lipids[1].lm_id == 'LMNO0001'


def test_fill_missing_lm_ids_handles_lmsd_error(monkeypatch):
    q = QuantifiedLipid(input_name='X', values={})
    mgr = DataManager()
    mgr.dataset = LipidDataset(samples=[], lipids=[q])

    def fake_get(names):
        return {'error': 'service unavailable'}

    monkeypatch.setattr(LMSD, 'get_lm_ids_by_name', fake_get)

    updated = mgr.fill_missing_lm_ids_from_lmsd()
    assert updated == 0
    assert mgr.dataset.lipids[0].lm_id is None
