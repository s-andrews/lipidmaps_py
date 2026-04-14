import logging
import sys
from pathlib import Path

import requests


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidmaps.data.models.refmet import RefMet


class FakeResponse:
    def __init__(self, text: str, status_code: int = 200, raise_exc=None):
        self._text = text
        self.status_code = status_code
        self._raise_exc = raise_exc

    @property
    def text(self):
        return self._text

    @property
    def content(self):
        return self._text.encode("utf-8")

    def raise_for_status(self):
        if self._raise_exc:
            raise self._raise_exc


def test_refmet_fetch_logging_and_results(monkeypatch, caplog):
    caplog.set_level(logging.INFO)

    response_text = "\t".join(["Input name", "Standardized name", "LM_ID"]) + "\n"
    response_text += "\t".join(["Alpha", "Alpha standardized", "LMFA00000001"]) + "\n"

    calls = []

    def fake_post(*args, **kwargs):
        calls.append((args, kwargs))
        return FakeResponse(response_text)

    monkeypatch.setattr(requests, "post", fake_post)

    results = RefMet.validate_metabolite_names(["Alpha"])

    assert len(calls) == 1
    assert len(results) == 1
    assert results[0].standardized_name == "Alpha standardized"
    assert results[0].lm_id == "LMFA00000001"
    assert "Fetching RefMet annotations for 1 metabolites from API" in caplog.text
    assert "RefMet API fetched" in caplog.text