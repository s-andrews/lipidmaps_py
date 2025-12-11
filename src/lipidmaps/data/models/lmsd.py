import logging
from typing import List, Dict, Optional, Union, Any
import requests
from pydantic import BaseModel

logger = logging.getLogger(__name__)


class LMSDResult(BaseModel):
    input_name: Optional[str] = None
    name: Optional[str] = None
    lm_id: Optional[str] = None
    sys_name: Optional[str] = None
    abbrev: Optional[str] = None
    abbrev_chains: Optional[Union[str, float]] = None
    matched_field: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

class LMSD:
    LMSDNameUrl = "http://localhost/api/reactions/names"

    @staticmethod
    def get_lm_ids_by_name(lipid_names: List[str]) -> Union[List[Dict[str, Any]], Dict[str, Any], None]:
        """Return lm_id's and associated names using LMSD API.

        Args:
            lipid_names: List of lipid name strings to validate

        Returns:
            List of dictionaries (serialized LMSDResult) one per input name,
            or an error dictionary with an `error` key on failure.
        """
        data = {"names": lipid_names}
        try:
            logger.info("Sending request to LMSD API")
            response = requests.post(
                LMSD.LMSDNameUrl, json=data, verify=False, timeout=20
            )
            response.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"LMSD API call failed: {e}")
            # Return error dict on failure
            return {"error": str(e)}
        # Try to parse JSON first (the LMSD endpoint often returns JSON)
        try:
            json_data = response.json()
        except ValueError:
            json_data = None

        if json_data is not None:
            # If the response is a dict that looks like an error, return it
            if isinstance(json_data, dict):
                if 'error' in json_data and len(json_data) == 1:
                    return json_data
                json_list = [json_data]
            elif isinstance(json_data, list):
                json_list = json_data
            else:
                json_list = None

            if json_list is not None:
                results: List[Dict[str, Any]] = []

                for item in json_list:  

                    res = LMSDResult(
                        input_name=item.get('input_name'),
                        matched_field=item.get('matched_field'),
                        name=item.get('name'),
                        sys_name=item.get('sys_name'),
                        abbrev=item.get('abbrev'),
                        abbrev_chains=item.get('abbrev_chains'),
                        lm_id=item.get('lm_id'),
                    )
                    results.append(res.to_dict())

                return results

        # Fallback: treat the response as TSV/text (legacy behaviour)
        lines = [ln for ln in response.text.splitlines() if ln.strip()]

        if not lines:
            logger.info("LMSD returned empty response")
            return []

        header = lines[0].split("\t")

        def idx(name: str) -> Optional[int]:
            return header.index(name) if name in header else None

        # Build a case-insensitive header map for flexible matching
        hdr = [h.strip() for h in header]
        hdr_map = {h.lower(): i for i, h in enumerate(hdr)}

        def find(*candidates: str) -> Optional[int]:
            for cand in candidates:
                cand_l = cand.lower()
                if cand_l in hdr_map:
                    return hdr_map[cand_l]

            return None

        input_idx = find('input_name')
        matched_idx = find('matched_field')
        name_idx = find('name')
        sys_name_idx = find('sys_name')
        abbrev_idx = find('abbrev')
        abbrev_chains_idx = find('abbrev_chains')
        lm_id_idx = find('lm_id')

        results: List[Dict[str, Any]] = []

        for ln in lines[1:]:
            cols = ln.split('\t')

            def get(i: Optional[int]) -> Optional[str]:
                if i is None:
                    return None
                if i < 0 or i >= len(cols):
                    return None
                val = cols[i].strip()
                return val if val != '' else None

            abbrev_chains_val = get(abbrev_chains_idx)
            try:
                abbrev_chains = float(abbrev_chains_val) if abbrev_chains_val is not None else None
            except (ValueError, TypeError):
                abbrev_chains = None

            res = LMSDResult(
                input_name=get(input_idx),
                matched_field=get(matched_idx),
                name=get(name_idx),
                sys_name=get(sys_name_idx),
                abbrev=get(abbrev_idx),
                abbrev_chains=abbrev_chains,
                lm_id=get(lm_id_idx),
            )

            results.append(res.to_dict())

        return results