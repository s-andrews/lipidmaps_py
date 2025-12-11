import os
import requests
import pytest

REACTIONS_URL = os.environ.get("REACTIONS_URL", "http://localhost/api/reactions")


def test_reactions_laravel_endpoint():
    payload = {
        "search_mode": "default",
        "search_type": "lipids",
        "generic_reactions": True,
        "lmsd_ids": ["LMGP06010000", "LMFA0103000", "LMFA010100150"],
    }
    try:
        r = requests.post(REACTIONS_URL, json=payload, timeout=10)
    except requests.RequestException as exc:
        pytest.skip(f"Cannot reach reactions endpoint {REACTIONS_URL}: {exc}")

    assert r.status_code == 200, f"Unexpected status {r.status_code}: {r.text}"
    data = r.json()
    assert isinstance(data, list), f"Expected list, got {type(data)}"
    assert data, "Empty response list"
    first = data[0]
    assert "reaction_id" in first, "Missing reaction_id in first item"
    assert "reactants" in first and "products" in first, "Missing reactants/products"