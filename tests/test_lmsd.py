import json
import requests
import pytest

from lipidmaps.data.models.lmsd import LMSD


class FakeResponse:
    def __init__(self, status_code=200, text='', headers=None, json_data=None, raise_exc=None):
        self.status_code = status_code
        self._text = text
        self.headers = headers or {'Content-Type': 'application/json'}
        self._json_data = json_data
        self._raise_exc = raise_exc

    def json(self):
        if isinstance(self._json_data, Exception):
            raise self._json_data
        return self._json_data

    def raise_for_status(self):
        if self._raise_exc:
            raise self._raise_exc

    @property
    def text(self):
        return self._text


def test_get_lm_ids_by_name_json(monkeypatch):
    sample = [
        {
            "input_name": "Butyrylcarnitine",
            "matched_field": "name",
            "name": "Butyrylcarnitine",
            "sys_name": "3-(butanoyloxy)-4-(trimethylazaniumyl)butanoate",
            "abbrev": "CAR 4:0",
            "abbrev_chains": None,
            "lm_id": "LMFA07070054",
        },
        {
            "input_name": "Cholesterol",
            "matched_field": "name",
            "name": "Cholesterol",
            "sys_name": "cholest-5-en-3beta-ol",
            "abbrev": "ST 27:1;O",
            "abbrev_chains": None,
            "lm_id": "LMST01010001",
        },
    ]

    def fake_post(*args, **kwargs):
        return FakeResponse(status_code=200, text=json.dumps(sample), json_data=sample)

    monkeypatch.setattr(requests, 'post', fake_post)

    res = LMSD.get_lm_ids_by_name(["Butyrylcarnitine", "Cholesterol"])
    assert isinstance(res, list)
    assert len(res) == 2
    assert res[0]['lm_id'] == 'LMFA07070054'
    assert res[1]['lm_id'] == 'LMST01010001'


def test_get_lm_ids_by_name_tsv_fallback(monkeypatch):
    header = '\t'.join(['input_name', 'matched_field', 'name', 'sys_name', 'abbrev', 'abbrev_chains', 'lm_id'])
    row1 = '\t'.join(['Butyrylcarnitine', 'name', 'Butyrylcarnitine', '3-(butanoyloxy)-4-(trimethylazaniumyl)butanoate', 'CAR 4:0', '', 'LMFA07070054'])
    row2 = '\t'.join(['Cholesterol', 'name', 'Cholesterol', 'cholest-5-en-3beta-ol', 'ST 27:1;O', '', 'LMST01010001'])
    tsv = header + '\n' + row1 + '\n' + row2 + '\n'

    def fake_post(*args, **kwargs):
        return FakeResponse(status_code=200, text=tsv, json_data=ValueError('no json'))

    monkeypatch.setattr(requests, 'post', fake_post)

    res = LMSD.get_lm_ids_by_name(["Butyrylcarnitine", "Cholesterol"])
    assert isinstance(res, list)
    assert len(res) == 2
    assert res[0]['lm_id'] == 'LMFA07070054'
    assert res[1]['lm_id'] == 'LMST01010001'


def test_get_lm_ids_by_name_empty_response(monkeypatch):
    def fake_post(*args, **kwargs):
        return FakeResponse(status_code=200, text='', json_data=ValueError('no json'))

    monkeypatch.setattr(requests, 'post', fake_post)

    res = LMSD.get_lm_ids_by_name(['NoMatch'])
    assert res == []


def test_get_lm_ids_by_name_request_exception(monkeypatch):
    def fake_post(*args, **kwargs):
        raise requests.RequestException('boom')

    monkeypatch.setattr(requests, 'post', fake_post)

    res = LMSD.get_lm_ids_by_name(['X'])
    assert isinstance(res, dict)
    assert 'error' in res
    assert 'boom' in res['error']
