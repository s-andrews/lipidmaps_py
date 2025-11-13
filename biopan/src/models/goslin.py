import logging
from typing import List, Dict, Optional, Any

import requests
from pydantic import BaseModel
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)


class GoslinResult(BaseModel):
    input_name: Optional[str] = None
    normalized_name: Optional[str] = None
    lipid_string: Optional[str] = None
    lm_id: Optional[str] = None
    raw: Dict[str, Any] = {}

    def to_dict(self) -> Dict[str, Any]:
        return self.dict()


class Goslin:
    """Client wrapper for the GOSLIN REST validate endpoint.

    Example request body:
    {
      "lipidNames": ["Cholesterol"],
      "grammars": ["GOSLIN"],
      "skipInvalid": true
    }

    The service returns a JSON object with a `results` list. We extract a
    normalized name and any LM reference (databaseElementId) when present.
    """

    BASE_URL = "https://apps.lifs-tools.org/goslin/rest/validate"

    @staticmethod
    def _session_with_retries(retries: int = 3, backoff: float = 0.3) -> requests.Session:
        s = requests.Session()
        retry = Retry(
            total=retries,
            read=retries,
            connect=retries,
            backoff_factor=backoff,
            status_forcelist=(500, 502, 503, 504),
            allowed_methods=frozenset(['POST', 'GET']),
        )
        adapter = HTTPAdapter(max_retries=retry)
        s.mount('https://', adapter)
        s.mount('http://', adapter)
        return s

    @staticmethod
    def annotate_samples(samples: List[Any], grammars: Optional[List[str]] = None, skip_invalid: bool = True) -> Dict[str, Dict[str, Any]]:
        """Annotate a list of sample-like objects.

        Each sample is expected to have a `sample_name` attribute (string).
        This function will attach a `goslin_result` attribute (a dict) to
        each sample and return a lookup mapping input_name -> result dict.
        """
        grammars = grammars or ["GOSLIN"]
        names = [getattr(s, 'sample_name', '') for s in samples]
        payload = {
            'lipidNames': names,
            'grammars': grammars,
            'skipInvalid': bool(skip_invalid),
        }

        session = Goslin._session_with_retries()
        try:
            logger.info(f"Sending request to GOSLIN API {Goslin.BASE_URL}")
            resp = session.post(Goslin.BASE_URL, json=payload, timeout=15)
            resp.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"GOSLIN API call failed: {e}")
            # Attach fallback results: no normalization
            for sample in samples:
                key = getattr(sample, 'sample_name', None)
                sample.goslin_result = {
                    'input_name': key,
                    'normalized_name': key or '',
                    'lipid_string': None,
                    'lm_id': None,
                    'raw': {},
                }
            return {}

        try:
            data = resp.json()
        except ValueError:
            logger.error('GOSLIN returned non-JSON response')
            return {}

        results = data.get('results', []) if isinstance(data, dict) else []

        lookup: Dict[str, Dict[str, Any]] = {}
        for r in results:
            input_name = r.get('lipidName') or r.get('input') or None
            normalized = r.get('normalizedName') or r.get('lipidString') or None
            lipid_string = r.get('lipidString') or None
            lm_id = None
            lm_refs = r.get('lipidMapsReferences') or []
            if lm_refs and isinstance(lm_refs, list) and len(lm_refs) > 0:
                # Prefer the first LM reference that has a databaseElementId
                for ref in lm_refs:
                    dbid = ref.get('databaseElementId') if isinstance(ref, dict) else None
                    if dbid and dbid.upper().startswith('LM'):
                        lm_id = dbid
                        break

            result = GoslinResult(
                input_name=input_name,
                normalized_name=normalized,
                lipid_string=lipid_string,
                lm_id=lm_id,
                raw=r,
            )

            if input_name:
                lookup[input_name] = result.dict()

        # Attach to samples; fallback to identity if no match
        for sample in samples:
            key = getattr(sample, 'sample_name', None)
            if key in lookup:
                sample.goslin_result = lookup[key]
            else:
                sample.goslin_result = {
                    'input_name': key,
                    'normalized_name': key or '',
                    'lipid_string': None,
                    'lm_id': None,
                    'raw': {},
                }
            logger.info(f"Annotated sample {key} with Goslin normalized '{sample.goslin_result.get('normalized_name')}' and lm_id {sample.goslin_result.get('lm_id')}")

        return lookup

    @staticmethod
    def get_lm_ids(samples: List[Any]) -> List[str]:
        lm_ids = []
        for s in samples:
            r = getattr(s, 'goslin_result', None)
            lm_id = None
            if isinstance(r, GoslinResult):
                lm_id = r.lm_id
            elif isinstance(r, dict):
                lm_id = r.get('lm_id')
            if lm_id:
                lm_ids.append(lm_id)
        return list(set(lm_ids))

    @staticmethod
    def get_unmatched_results(samples: List[Any]) -> List[str]:
        unmatched = []
        for s in samples:
            r = getattr(s, 'goslin_result', None)
            normalized = None
            lm_id = None
            if isinstance(r, GoslinResult):
                normalized = r.normalized_name
                lm_id = r.lm_id
            elif isinstance(r, dict):
                normalized = r.get('normalized_name')
                lm_id = r.get('lm_id')
            if normalized and lm_id is None:
                unmatched.append(normalized)
        return list(set(unmatched))
