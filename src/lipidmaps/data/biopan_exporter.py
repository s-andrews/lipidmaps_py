import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from pydantic import ConfigDict, Field

from .models.base import LipidmapsBaseModel
from .models.sample import LipidDataset, QuantifiedLipid


def _ordered_unique_strings(values: List[Any]) -> List[str]:
    """Return non-empty string values once, preserving input order."""
    seen = set()
    ordered = []
    for value in values:
        text = str(value).strip() if value is not None else ""
        if not text or text in seen:
            continue
        seen.add(text)
        ordered.append(text)
    return ordered


class BioPANExporter(LipidmapsBaseModel):
    """Build and write the BioPAN JSON files consumed by the PHP frontend."""

    dataset: Optional[Any] = Field(default=None)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def _get_dataset(self, dataset: Optional[LipidDataset] = None) -> LipidDataset:
        resolved_dataset = dataset or self.dataset
        if resolved_dataset is None:
            raise ValueError("BioPAN export requires a populated dataset")
        return resolved_dataset

    @staticmethod
    def _get_output_dir(output_path: Union[str, Path]) -> Path:
        output_dir = Path(output_path)
        if output_dir.name != "biopan":
            output_dir = output_dir / "biopan"
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    @staticmethod
    def _is_fa_like(name: str) -> bool:
        fa_name = name.replace(" ", "").upper()
        return fa_name.startswith("FA(") or fa_name.startswith("FACOA(")

    def _get_has_fa_status(self, dataset: LipidDataset) -> str:
        names = [lipid.input_name for lipid in dataset.lipids if lipid.input_name]
        has_fa = any(name.replace(" ", "").upper().startswith("FA(") for name in names)
        has_facoa = any(name.replace(" ", "").upper().startswith("FACOA(") for name in names)

        if has_fa and not has_facoa:
            return "facoa"
        if not has_fa and not has_facoa:
            return "none"
        return ""

    def _classify_species(self, dataset: LipidDataset) -> Dict[str, List[str]]:
        dataset.fill_generic_lm_ids_from_headgroups()

        processed = []
        unprocessed = []
        undefined = []

        for lipid in dataset.lipids:
            lipid_name = (lipid.input_name or "").strip()
            if not lipid_name:
                continue

            has_reactions = bool(getattr(lipid, "reactions", None))
            is_known_species = bool(getattr(lipid, "generic_lm_id", None) or getattr(lipid, "lm_id", None))

            if has_reactions:
                processed.append(lipid_name)
            elif is_known_species:
                unprocessed.append(lipid_name)
            else:
                undefined.append(lipid_name)

        return {
            "processed": _ordered_unique_strings(processed),
            "unprocessed": _ordered_unique_strings(unprocessed),
            "undef": _ordered_unique_strings(undefined),
        }

    @staticmethod
    def _get_dataset_display_name(lipid: QuantifiedLipid) -> str:
        standardized_name = (getattr(lipid, "standardized_name", None) or "").strip()
        if standardized_name:
            return f" {standardized_name}"
        return ""

    def _build_dataset_name_map(
        self,
        dataset: LipidDataset,
        selected_names: List[str],
    ) -> Dict[str, str]:
        selected_lookup = set(selected_names)
        mapped_names: Dict[str, str] = {}
        for lipid in dataset.lipids:
            lipid_name = (lipid.input_name or "").strip()
            if not lipid_name or lipid_name not in selected_lookup or lipid_name in mapped_names:
                continue
            mapped_names[lipid_name] = self._get_dataset_display_name(lipid)
        return mapped_names

    def build_summary(
        self,
        dataset: Optional[LipidDataset] = None,
        lipidlynxx: str = "no",
    ) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        classified = self._classify_species(resolved_dataset)
        groups = resolved_dataset.get_sample_conditions().model_dump()
        undef = classified["undef"]
        processed = classified["processed"]
        unprocessed = classified["unprocessed"]

        return {
            "total": len(resolved_dataset.lipids),
            "groups": groups,
            "lipidlynxx": lipidlynxx,
            "undef": undef if undef else {},
            "undef_dataset": [],
            "processed_dataset": self._build_dataset_name_map(resolved_dataset, processed),
            "unprocessed_dataset": self._build_dataset_name_map(resolved_dataset, unprocessed),
            "pathway": {
                "name": "Pathway anlysis",
                "processed": processed,
                "unprocessed": unprocessed,
            },
        }

    def build_msg1(
        self,
        dataset: Optional[LipidDataset] = None,
        lipidlynxx_error: bool = False,
    ) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        column_info = getattr(resolved_dataset, "column_info", None) or {}

        duplicates = [
            lipid_name
            for lipid_name, count in Counter(
                lipid.input_name for lipid in resolved_dataset.lipids if lipid.input_name
            ).items()
            if count > 1
        ]
        has_missing_values = any(
            value is None
            for lipid in resolved_dataset.lipids
            for value in lipid.values.values()
        )

        return {
            "na_values": [has_missing_values],
            "empty_rows": [False],
            "empty_cols": [bool(column_info.get("empty_columns"))],
            "dup_lp": _ordered_unique_strings(duplicates) or {},
            "has_fa": [self._get_has_fa_status(resolved_dataset)],
            "lipidlynxx_error": [lipidlynxx_error],
        }

    def build_msg2(self, dataset: Optional[LipidDataset] = None) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        classified = self._classify_species(resolved_dataset)

        group_counts: Dict[str, int] = {}
        for sample in resolved_dataset.samples:
            if not sample.group:
                continue
            group_counts.setdefault(sample.group, 0)
            group_counts[sample.group] += 1

        valid_groups = [group for group, count in group_counts.items() if count > 1]
        valid_freqs = [group_counts[group] for group in valid_groups]
        invalid_groups = [group for group, count in group_counts.items() if count <= 1]
        invalid_freqs = [group_counts[group] for group in invalid_groups]

        processed_fa = [name for name in classified["processed"] if self._is_fa_like(name)]
        processed_lipids = [name for name in classified["processed"] if not self._is_fa_like(name)]
        has_processed_content = bool(classified["processed"])

        notvalid: Union[List[Any], Dict[str, Any]]
        if invalid_groups:
            notvalid = {"groups": invalid_groups, "freqs": invalid_freqs}
        else:
            notvalid = []

        return {
            "valid": {"groups": valid_groups, "freqs": valid_freqs},
            "notvalid": notvalid,
            "reaction": {
                "lp": [bool(processed_lipids)],
                "fa": [bool(processed_fa)],
            },
            "subset": {
                "reaction": [has_processed_content],
                "pathway": [has_processed_content],
            },
        }

    def export_display_files(
        self,
        output_path: Union[str, Path],
        dataset: Optional[LipidDataset] = None,
        lipidlynxx: str = "no",
        lipidlynxx_error: bool = False,
        include_msg1: bool = True,
        include_summary: bool = True,
        include_msg2: bool = True,
    ) -> Dict[str, str]:
        resolved_dataset = self._get_dataset(dataset)
        output_dir = self._get_output_dir(output_path)
        written_files: Dict[str, str] = {}

        file_payloads = []
        if include_msg1:
            file_payloads.append(("msg1", output_dir / "msg1.json", self.build_msg1(resolved_dataset, lipidlynxx_error=lipidlynxx_error)))
        if include_summary:
            file_payloads.append(("summary", output_dir / "summary.json", self.build_summary(resolved_dataset, lipidlynxx=lipidlynxx)))
        if include_msg2:
            file_payloads.append(("msg2", output_dir / "msg2.json", self.build_msg2(resolved_dataset)))

        for key, file_path, payload in file_payloads:
            with file_path.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            written_files[key] = str(file_path)

        return written_files