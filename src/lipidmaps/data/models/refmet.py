import logging
from typing import List, Dict, Optional, Any
import requests
from pydantic import BaseModel

logger = logging.getLogger(__name__)


class RefMetResult(BaseModel):
    input_name: Optional[str] = None
    standardized_name: Optional[str] = None
    lm_id: Optional[str] = None
    sub_class: Optional[str] = None
    formula: Optional[str] = None
    exact_mass: Optional[float] = None
    super_class: Optional[str] = None
    main_class: Optional[str] = None
    chebi_id: Optional[str] = None
    kegg_id: Optional[str] = None
    refmet_id: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()


class RefMet:
    MWBaseURL = "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php"

    @staticmethod
    def validate_metabolite_names(metabolite_names: List[str]) -> List[RefMetResult]:
        """Validate metabolite names using RefMet API and return RefMetResult objects.

        Args:
            metabolite_names: List of metabolite name strings to validate

        Returns:
            List of RefMetResult objects, one per input name
        """
        data = {"metabolite_name": "\n".join(metabolite_names)}
        try:
            logger.info("Sending request to RefMet API")
            response = requests.post(
                RefMet.MWBaseURL, data=data, verify=True, timeout=20
            )
            response.raise_for_status()
        except requests.RequestException as e:
            logger.error("RefMet API call failed: %s", e)
            # Return empty list on failure
            return {"error": str(e)}

        lines = [ln for ln in response.text.splitlines() if ln.strip()]
        if not lines:
            logger.warning("RefMet returned empty response")
            return []

        header = lines[0].split("\t")

        def idx(name: str) -> Optional[int]:
            return header.index(name) if name in header else None

        input_idx = idx("Input name")
        standardized_idx = idx("Standardized name")
        lm_id_idx = idx("LM_ID")
        formula_idx = idx("Formula")
        exact_mass_idx = idx("Exact mass")
        super_class_idx = idx("Super class")
        main_class_idx = idx("Main class")
        sub_class_idx = idx("Sub class")
        chebi_id_idx = idx("ChEBI_ID")
        kegg_id_idx = idx("KEGG_ID")
        refmet_id_idx = idx("RefMet_ID")

        refmet_results: List[RefMetResult] = []
        lookup: Dict[str, RefMetResult] = {}

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
            formula = (
                fields[formula_idx]
                if formula_idx is not None and len(fields) > formula_idx
                else None
            )
            exact_mass = (
                fields[exact_mass_idx]
                if exact_mass_idx is not None and len(fields) > exact_mass_idx
                else None
            )
            super_class = (
                fields[super_class_idx]
                if super_class_idx is not None and len(fields) > super_class_idx
                else None
            )
            main_class = (
                fields[main_class_idx]
                if main_class_idx is not None and len(fields) > main_class_idx
                else None
            )
            sub_class = (
                fields[sub_class_idx]
                if sub_class_idx is not None and len(fields) > sub_class_idx
                else None
            )
            chebi_id = (
                fields[chebi_id_idx]
                if chebi_id_idx is not None and len(fields) > chebi_id_idx
                else None
            )
            kegg_id = (
                fields[kegg_id_idx]
                if kegg_id_idx is not None and len(fields) > kegg_id_idx
                else None
            )
            refmet_id = (
                fields[refmet_id_idx]
                if refmet_id_idx is not None and len(fields) > refmet_id_idx
                else None
            )

            # Convert dash values to None
            input_name = None if input_name == "-" else input_name
            standardized = None if standardized == "-" else standardized
            lm_id = None if lm_id == "-" else lm_id
            sub_class = None if sub_class == "-" else sub_class
            formula = None if formula == "-" else formula
            exact_mass = None if exact_mass == "-" else exact_mass
            super_class = None if super_class == "-" else super_class
            main_class = None if main_class == "-" else main_class
            chebi_id = None if chebi_id == "-" else chebi_id
            kegg_id = None if kegg_id == "-" else kegg_id
            refmet_id = None if refmet_id == "-" else refmet_id

            result = RefMetResult(
                input_name=input_name,
                standardized_name=standardized,
                lm_id=(lm_id if lm_id else None),
                sub_class=(sub_class if sub_class else None),
                formula=(formula if formula else None),
                exact_mass=(exact_mass if exact_mass else None),
                super_class=(super_class if super_class else None),
                main_class=(main_class if main_class else None),
                chebi_id=(chebi_id if chebi_id else None),
                kegg_id=(kegg_id if kegg_id else None),
                refmet_id=(refmet_id if refmet_id else None),
            )

            logger.info(f"RefMet result: {result.standardized_name}")
            if input_name:
                lookup[input_name] = result

        # Ensure we have a result for each input name (fallback to identity if no match)
        for name in metabolite_names:
            if name in lookup:
                refmet_results.append(lookup[name])
            else:
                refmet_results.append(
                    RefMetResult(input_name=name, standardized_name=name)
                )
            logger.info(
                f"Validated '{name}' -> standardized: '{refmet_results[-1].standardized_name}', lm_id: {refmet_results[-1].lm_id}"
            )

        logger.info(f"Annotated {len(refmet_results)} metabolites via RefMet")
        return refmet_results

    @staticmethod
    def attach_results_to_samples(
        samples: List[Any], results: List[RefMetResult]
    ) -> None:
        """Attach RefMetResult objects to sample objects as refmet_result attribute.

        Args:
            samples: List of sample objects with sample_name attribute
            results: List of RefMetResult objects (must match length of samples)
        """
        if len(samples) != len(results):
            raise ValueError(
                f"Length mismatch: {len(samples)} samples vs {len(results)} results"
            )

        for sample, result in zip(samples, results):
            sample.refmet_result = result.model_dump()
            logger.debug(
                f"Attached RefMet result to sample {getattr(sample, 'sample_name', 'unknown')}"
            )

    @staticmethod
    def annotate_samples(samples: List[Any]) -> List[RefMetResult]:
        """Convenience method: validate metabolite names and attach results to samples.

        This is a convenience wrapper that calls validate_metabolite_names and
        attach_results_to_samples. For more control, use those methods separately.

        Args:
            samples: List of sample objects with sample_name attribute

        Returns:
            List of RefMetResult objects
        """
        names = [s.sample_name for s in samples]
        results = RefMet.validate_metabolite_names(names)
        RefMet.attach_results_to_samples(samples, results)
        return results

    @staticmethod
    def get_lm_ids(samples: List[Any]) -> List[str]:
        """Return unique LM_IDs from samples that have a refmet_result with lm_id starting with 'LM'."""
        lm_ids = []
        for s in samples:
            r = getattr(s, "refmet_result", None)
            lm_id = None
            if isinstance(r, RefMetResult):
                lm_id = r.lm_id
            elif isinstance(r, dict):
                lm_id = r.get("lm_id")
            if lm_id and lm_id.startswith("LM"):
                lm_ids.append(lm_id)
        return list(set(lm_ids))

    @staticmethod
    def get_unmatched_results(samples: List[Any]) -> List[str]:
        """Return unique standardized names for samples that were not matched (lm_id is None)."""
        unmatched = []
        for s in samples:
            r = getattr(s, "refmet_result", None)
            standardized = None
            lm_id = None
            if isinstance(r, RefMetResult):
                standardized = r.standardized_name
                lm_id = r.lm_id
            elif isinstance(r, dict):
                standardized = r.get("standardized_name")
                lm_id = r.get("lm_id")
            if standardized and lm_id is None:
                unmatched.append(standardized)
        # logger.info("Unmatched samples: %s", unmatched)
        return list(set(unmatched))
