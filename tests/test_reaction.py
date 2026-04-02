import sys
from pathlib import Path

import requests


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidmaps.data.models.reaction import ReactionChecker


class FakeResponse:
    status_code = 200
    text = "[]"

    def raise_for_status(self):
        return None

    def json(self):
        return []