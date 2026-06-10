import sys
from pathlib import Path

import requests
import pytest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidmaps.data.models.reaction import ReactionChecker


class FakeResponse:
    status_code = 200
    text = "[]"

    def raise_for_status(self):
        return None

    def json(self):
        return []


class PayloadResponse:
    status_code = 200

    def __init__(self, payload):
        self._payload = payload
        self.text = str(payload)

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def test_check_reactions_preserves_genes_for_generic_reactions(monkeypatch):
    payload = [
        {
            "reaction_id": 55,
            "genes": [
                {
                    "taxonomy_id": 9606,
                    "uniprot_id": "P21964",
                    "gene_name": "COMT",
                }
            ],
            "reactants": [
                {
                    "compound_type": "lm_main",
                    "compound_lm_id": "LMSP03010000",
                    "compound_name": "Sphingomyelin",
                },
                {
                    "compound_type": "lm_notlipids",
                    "compound_name": "water",
                },
            ],
            "products": [
                {
                    "compound_type": "lm_main",
                    "compound_lm_id": "LMSP02010000",
                    "compound_name": "Ceramide",
                },
                {
                    "compound_type": "lm_notlipids",
                    "compound_name": "phosphocholine",
                },
            ],
        }
    ]

    def fake_post(url, json, timeout):
        assert json["generic_reactions"] is True
        assert json["search_type"] == "all"
        assert json["lmsd_ids"] == ["LMSP03010000"]
        return PayloadResponse(payload)

    monkeypatch.setattr(requests, "post", fake_post)

    checker = ReactionChecker(base_url="http://localhost")
    response = checker.check_reactions(
        lm_ids=["LMSP03010000"],
        search_type="all",
        generic_reactions=True,
        only_lipid_components=True,
    )

    assert len(response.reactions) == 1
    reaction = response.reactions[0]
    assert reaction.reaction_id == 55
    assert reaction.genes == [
        {
            "taxonomy_id": 9606,
            "uniprot_id": "P21964",
            "gene_name": "COMT",
        }
    ]
    assert [component.compound_name for component in reaction.reactants] == ["Sphingomyelin"]
    assert [component.compound_name for component in reaction.products] == ["Ceramide"]