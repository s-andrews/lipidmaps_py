import logging
from typing import List, Dict, Optional, Any
import requests
from pydantic import BaseModel

logger = logging.getLogger(__name__)

"""
This template module provides functionality to check and standardize lipid names
using the LipidMaps database API. Work in progress.
"""


class LMNameResult(BaseModel):
    input_name: Optional[str] = None
    standardized_name: str
    lm_id: Optional[str] = None
    sub_class: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()  # Pydantic v2 method instead of .dict()


""" Refmet LMSD mapping is not providing enough information. We should also check LipidMaps by name"""


class LipidMapsNameChecker:
    LMBaseURL = "https://www.lipidmaps.org/api/lm_lookup"

    @staticmethod
    def annotate_samples(samples: List[Any]) -> Dict[str, Dict[str, Any]]:
        """
        Query LipidMaps for a list of samples and attach a LMNameResult (as dict) to each sample
        in sample.lm_result. Returns a lookup mapping input_name -> result dict.
        Samples are expected to have a 'sample_name' attribute.
        """
        names = [s.sample_name for s in samples]
        data = {"metabolite_name": "\n".join(names)}
        try:
            logger.info("Sending request to LipidMaps API")
            response = requests.post(
                LipidMapsNameChecker.LMBaseURL, data=data, verify=False, timeout=10
            )
            response.raise_for_status()
        except requests.RequestException as e:
            logger.error("LipidMaps API call failed: %s", e)
            return {}

        lines = [ln for ln in response.text.splitlines() if ln.strip()]
        if not lines:
            logger.warning("LipidMaps returned empty response")
            return {}

        header = lines[0].split("\t")

        def idx(name: str) -> Optional[int]:
            return header.index(name) if name in header else None

        input_idx = idx("Input name")
        standardized_idx = idx("Standardized name")
        lm_id_idx = idx("LM_ID")
        sub_class_idx = idx("Sub class")

        # store plain dicts so callers that expect dicts (.get) keep working
        lm_lookup: Dict[str, Dict[str, Any]] = {}
        for line in lines[1:]:
            fields = line.split("\t")
            input_name = (
                fields[input_idx]
                if input_idx is not None and len(fields) > input_idx
                else None
            )
            standardized = (
                fields[standardized_idx]
                if standardized_idx is not None and len(fields) > standardized_idx
                else None
            )
            lm_id = (
                fields[lm_id_idx]
                if lm_id_idx is not None and len(fields) > lm_id_idx
                else None
            )
            sub_class = (
                fields[sub_class_idx]
                if sub_class_idx is not None and len(fields) > sub_class_idx
                else None
            )

            # Refmet returns "-" if no value is found
            # Convert dash values to None
            input_name = None if input_name == "-" else input_name
            standardized = None if standardized == "-" else standardized
            lm_id = None if lm_id == "-" else lm_id
            sub_class = None if sub_class == "-" else sub_class

            result = LMNameResult(
                input_name=input_name,
                standardized_name=standardized or "",
                lm_id=lm_id,
                sub_class=sub_class,
            )
            lm_lookup[input_name] = result.to_dict()

        return lm_lookup
