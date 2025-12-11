import json
import requests

from lipidmaps.data.data_manager import DataManager


def test_get_reactions_for_lm_ids(monkeypatch):
    # Minimal synthetic API response mimicking the Laravel endpoint structure
    sample_raw = [
        {
            "id": 1,
            "reactants": [
                {
                    "compound_type": "lm_main",
                    "compound_name": "FA(22:5)",
                    "compound_lm_id": "LMFA04000049",
                }
            ],
            "products": [
                {
                    "compound_type": "lm_main",
                    "compound_name": "PC(16:0/18:1)",
                    "compound_lm_id": "LMGP06010000",
                }
            ],
        }
    ]

    class DummyResp:
        def __init__(self, json_data, status_code=200):
            self._json = json_data
            self.status_code = status_code

        def json(self):
            return self._json

        def raise_for_status(self):
            if not (200 <= self.status_code < 300):
                raise requests.HTTPError()

    def fake_post(url, json=None, timeout=None):
        return DummyResp(sample_raw)

    # Patch requests.post used by reaction_checker
    monkeypatch.setattr("lipidmaps.data.reaction_checker.requests.post", fake_post)

    mgr = DataManager()
    result = mgr.get_reactions_for_lm_ids(["LMFA04000049", "LMGP06010000"], base_url="http://localhost")

    assert isinstance(result, dict)
    assert "reactions" in result
    assert isinstance(result["reactions"], list)
    assert len(result["reactions"]) == 1
    reaction = result["reactions"][0]
    assert "reactants" in reaction and "products" in reaction
    assert len(reaction["reactants"]) == 1
    assert len(reaction["products"]) == 1
