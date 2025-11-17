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
        return self.model_dump()


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
    def _session_with_retries(
        retries: int = 3, backoff: float = 0.3
    ) -> requests.Session:
        s = requests.Session()
        retry = Retry(
            total=retries,
            read=retries,
            connect=retries,
            backoff_factor=backoff,
            status_forcelist=(500, 502, 503, 504),
            allowed_methods=frozenset(["POST", "GET"]),
        )
        adapter = HTTPAdapter(max_retries=retry)
        s.mount("https://", adapter)
        s.mount("http://", adapter)
        return s

    @staticmethod
    def validate_lipid_names(
        lipid_names: List[str],
        grammars: Optional[List[str]] = None,
        skip_invalid: bool = True,
    ) -> List[GoslinResult]:
        """Validate lipid names using GOSLIN API and return GoslinResult objects.

        Args:
            lipid_names: List of lipid name strings to validate
            grammars: List of grammar types to use (default: ["GOSLIN"])
            skip_invalid: Whether to skip invalid names (default: True)

        Returns:
            List of GoslinResult objects, one per input name
        """
        grammars = grammars or ["GOSLIN"]
        payload = {
            "lipidNames": lipid_names,
            "grammars": grammars,
            "skipInvalid": bool(skip_invalid),
        }

        session = Goslin._session_with_retries()
        try:
            logger.info(f"Sending request to GOSLIN API {Goslin.BASE_URL}")
            resp = session.post(Goslin.BASE_URL, json=payload, timeout=15)
            resp.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"GOSLIN API call failed: {e}")
            # Return fallback results
            return [
                GoslinResult(
                    input_name=name,
                    normalized_name=name,
                    lipid_string=None,
                    lm_id=None,
                    raw={},
                )
                for name in lipid_names
            ]

        try:
            data = resp.json()
        except ValueError:
            logger.error("GOSLIN returned non-JSON response")
            return []

        results = data.get("results", []) if isinstance(data, dict) else []

        goslin_results: List[GoslinResult] = []
        lookup: Dict[str, GoslinResult] = {}
        
        for r in results:
            input_name = r.get("lipidName") or r.get("input") or None
            normalized = r.get("normalizedName") or r.get("lipidString") or None
            lipid_string = r.get("lipidString") or None
            lm_id = None
            lm_refs = r.get("lipidMapsReferences") or []
            if lm_refs and isinstance(lm_refs, list) and len(lm_refs) > 0:
                # Prefer the first LM reference that has a databaseElementId
                for ref in lm_refs:
                    dbid = (
                        ref.get("databaseElementId") if isinstance(ref, dict) else None
                    )
                    if dbid and dbid.upper().startswith("LM"):
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
                lookup[input_name] = result

        # Ensure we have a result for each input name (fallback to identity if no match)
        for name in lipid_names:
            if name in lookup:
                goslin_results.append(lookup[name])
            else:
                goslin_results.append(
                    GoslinResult(
                        input_name=name,
                        normalized_name=name,
                        lipid_string=None,
                        lm_id=None,
                        raw={},
                    )
                )
            logger.info(
                f"Validated '{name}' -> normalized: '{goslin_results[-1].normalized_name}', lm_id: {goslin_results[-1].lm_id}"
            )

        return goslin_results

    @staticmethod
    def attach_results_to_samples(
        samples: List[Any], results: List[GoslinResult]
    ) -> None:
        """Attach GoslinResult objects to sample objects as goslin_result attribute.

        Args:
            samples: List of sample objects with sample_name attribute
            results: List of GoslinResult objects (must match length of samples)
        """
        if len(samples) != len(results):
            raise ValueError(
                f"Length mismatch: {len(samples)} samples vs {len(results)} results"
            )

        for sample, result in zip(samples, results):
            sample.goslin_result = result.model_dump()
            logger.debug(
                f"Attached result to sample {getattr(sample, 'sample_name', 'unknown')}"
            )

    @staticmethod
    def annotate_samples(
        samples: List[Any],
        grammars: Optional[List[str]] = None,
        skip_invalid: bool = True,
    ) -> List[GoslinResult]:
        """Convenience method: validate lipid names and attach results to samples.

        This is a convenience wrapper that calls validate_lipid_names and
        attach_results_to_samples. For more control, use those methods separately.

        Args:
            samples: List of sample objects with sample_name attribute
            grammars: List of grammar types to use (default: ["GOSLIN"])
            skip_invalid: Whether to skip invalid names (default: True)

        Returns:
            List of GoslinResult objects
        """
        names = [getattr(s, "sample_name", "") for s in samples]
        results = Goslin.validate_lipid_names(names, grammars, skip_invalid)
        Goslin.attach_results_to_samples(samples, results)
        return results

    @staticmethod
    def get_lm_ids(samples: List[Any]) -> List[str]:
        lm_ids = []
        for s in samples:
            r = getattr(s, "goslin_result", None)
            lm_id = None
            if isinstance(r, GoslinResult):
                lm_id = r.lm_id
            elif isinstance(r, dict):
                lm_id = r.get("lm_id")
            if lm_id:
                lm_ids.append(lm_id)
        return list(set(lm_ids))

    @staticmethod
    def get_unmatched_results(samples: List[Any]) -> List[str]:
        unmatched = []
        for s in samples:
            r = getattr(s, "goslin_result", None)
            normalized = None
            lm_id = None
            if isinstance(r, GoslinResult):
                normalized = r.normalized_name
                lm_id = r.lm_id
            elif isinstance(r, dict):
                normalized = r.get("normalized_name")
                lm_id = r.get("lm_id")
            if normalized and lm_id is None:
                unmatched.append(normalized)
        return list(set(unmatched))
