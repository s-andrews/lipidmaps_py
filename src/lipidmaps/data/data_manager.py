import csv
import logging
import re
from typing import List, Tuple, Dict, Any, Union, Optional
from pathlib import Path
import pandas as pd

from pydantic import BaseModel, Field, field_validator

# import the data models we will produce
from .models.sample import SampleMetadata, QuantifiedLipid, LipidDataset
from .models.refmet import RefMet
from .models.lmsd import LMSD

# import new ingestion and validation modules
from .ingestion.csv_reader import CSVIngestion, CSVFormat
from .validation.data_validator import DataValidator, ValidationReport


logger = logging.getLogger(__name__)


class DataManager(BaseModel):
    """Pydantic v2 DataManager: reads CSVs into LipidDataset objects.

    Now uses CSVIngestion for file reading and DataValidator for quality checks.

    Usage:
        # Default behavior - first column is lipid names, rest are samples
        mgr = DataManager()
        dataset = mgr.process_csv("tests/inputs/quantified_test_file.csv")

        # Specify custom columns:
        mgr = DataManager(
            lipid_name_column=0,  # or column name as string
            sample_columns=[1, 2, 3],  # or column names as list of strings
        )
        dataset = mgr.process_csv("data.csv")

        # Specify group-to-sample mapping:
        mgr = DataManager(
            group_mapping={
                "Control": ["Sample1", "Sample2"],
                "Treatment": ["Sample3", "Sample4"]
            }
        )
        dataset = mgr.process_csv("data.csv")

        # With validation:
        mgr = DataManager(validate_data=True)
        dataset = mgr.process_csv("tests/inputs/quantified_test_file.csv")
        if mgr.validation_report and not mgr.validation_report.passed:
            mgr.validation_report.print_report()

    The CSV is expected to have a first column with lipid name (e.g. NAME)
    and subsequent columns one column per sample (sample ids as headers).
    """

    dataset: Optional[LipidDataset] = Field(default=None)
    lipid_species: List[Any] = Field(default_factory=list)
    validation_report: Optional[ValidationReport] = Field(default=None)

    # Configuration for ingestion and validation
    validate_data: bool = Field(default=False)
    csv_format: CSVFormat = Field(default=CSVFormat.AUTO)
    has_labels: bool = Field(default=False)

    # User-specified column configuration
    lipid_name_column: Optional[Union[int, str]] = Field(
        default=0,
        description="Column index (0-based) or name for lipid names. Default: first column (0)",
    )
    sample_columns: Optional[Union[List[int], List[str]]] = Field(
        default=None,
        description="List of column indices or names for sample data. Default: all columns after lipid column",
    )

    # Group-to-sample mapping
    group_mapping: Optional[Dict[str, List[str]]] = Field(
        default=None,
        description="Map group names to sample IDs. Example: {'Control': ['S1', 'S2'], 'Treatment': ['S3', 'S4']}",
    )

    group_label: Optional[str] = Field(
        default=None,
        description="Label to use for groups - Usually derived from group_mapping or second row in CSV.",
    )

    model_config = {"arbitrary_types_allowed": True}

    @field_validator("sample_columns", mode="before")
    @classmethod
    def validate_sample_columns(cls, v):
        """Ensure sample_columns is a list if provided."""
        if v is None:
            return None
        if isinstance(v, (int, str)):
            return [v]
        return v

    def model_post_init(self, __context: dict) -> None:
        logger.info(
            f"Initialized DataManager (validation={self.validate_data}, lipid_col={self.lipid_name_column}, "
            f"sample_cols={len(self.sample_columns) if self.sample_columns else 'auto'}, "
            f"groups={len(self.group_mapping) if self.group_mapping else 'auto'})"
        )

    def process_csv(self, csv_path: Union[str, Path]) -> LipidDataset:
        """Read CSV and populate SampleMetadata, QuantifiedLipid and LipidDataset.

        Now uses CSVIngestion for reading and DataValidator for quality checks.

        Args:
            csv_path: Path to CSV file

        Returns:
            LipidDataset with processed data
        """
        csv_path = Path(csv_path)
        logger.info(f"Loading CSV file: {csv_path}")

        # Use CSVIngestion to read file
        ingestion = CSVIngestion(has_labels=self.has_labels)
        raw_df = ingestion.read_csv(csv_path, format_type=self.csv_format)
        column_info = ingestion.get_column_info(raw_df)
        logger.info(
            f"CSV ingested: {raw_df.row_count} rows x {raw_df.column_count} columns."
        )

        # {
        #     "column_count": len(raw_df.fieldnames),
        #     "columns": raw_df.fieldnames,
        #     "empty_columns": [],
        #     "column_types": {},
        # }

        # Validate data if requested
        if self.validate_data:
            validator = DataValidator()
            self.validation_report = validator.validate(raw_df)

            if not self.validation_report.passed:
                logger.warning(
                    f"Validation found {len(self.validation_report.issues)} issues"
                )
                # Optionally print report
                # self.validation_report.print_report()

        # Process the raw data
        if raw_df.is_empty():
            ds = LipidDataset(samples=[], lipids=[])
            self.dataset = ds
            return ds

        # Determine lipid name column
        name_col = self._resolve_lipid_column(raw_df.fieldnames)

        # Determine sample columns
        sample_ids = self._resolve_sample_columns(raw_df.fieldnames, name_col)
        labels = raw_df.labels if hasattr(raw_df, "labels") else []
        # Create sample metadata with group mapping if provided
        samples_meta = self.extract_sample_metadata(sample_ids, labels=labels)

        # logger.info(f"column_info: {column_info.get('column_types', {})}")
        # logger.info(f"Empty: {column_info.get('empty_columns', [])}")
        # Extract quantified lipids
        quantified = self.extract_quantified_lipids(raw_df.rows, name_col, sample_ids, column_info)
        self.annotate_lipids_with_refmet(quantified)

        dataset = LipidDataset(samples=samples_meta, lipids=quantified, column_info=column_info)
        self.dataset = dataset
        logger.info(
            f"Created LipidDataset: {len(samples_meta)} samples, {len(quantified)} lipids"
        )
        return dataset

    def _resolve_lipid_column(self, fieldnames: List[str]) -> str:
        """Resolve the lipid name column from user specification or default.

        Args:
            fieldnames: List of column names from CSV

        Returns:
            Column name to use for lipid names
        """
        if self.lipid_name_column is None:
            # Default: first column
            return fieldnames[0]

        if isinstance(self.lipid_name_column, int):
            # Column index specified
            if 0 <= self.lipid_name_column < len(fieldnames):
                return fieldnames[self.lipid_name_column]
            else:
                raise ValueError(
                    f"lipid_name_column index {self.lipid_name_column} out of range. "
                    f"CSV has {len(fieldnames)} columns."
                )

        # Column name specified
        if self.lipid_name_column in fieldnames:
            return self.lipid_name_column
        else:
            raise ValueError(
                f"lipid_name_column '{self.lipid_name_column}' not found in CSV. "
                f"Available columns: {fieldnames}"
            )

    def _resolve_sample_columns(
        self, fieldnames: List[str], lipid_col: str
    ) -> List[str]:
        """Resolve sample columns from user specification or default.

        Args:
            fieldnames: List of column names from CSV
            lipid_col: The lipid name column (to exclude)

        Returns:
            List of column names to use for sample data
        """
        if self.sample_columns is None:
            # Default: all columns except lipid column
            return [
                col for col in fieldnames if col != lipid_col and col and col.strip()
            ]

        resolved = []
        for spec in self.sample_columns:
            if isinstance(spec, int):
                # Column index specified
                if 0 <= spec < len(fieldnames):
                    resolved.append(fieldnames[spec])
                else:
                    raise ValueError(
                        f"sample_columns index {spec} out of range. "
                        f"CSV has {len(fieldnames)} columns."
                    )
            else:
                # Column name specified
                if spec in fieldnames:
                    resolved.append(spec)
                else:
                    raise ValueError(
                        f"sample_columns name '{spec}' not found in CSV. "
                        f"Available columns: {fieldnames}"
                    )

        return [col for col in resolved if col and col.strip()]


    def extract_sample_metadata(self, sample_ids: List[str], labels: Optional[List[str]]=None) -> List[SampleMetadata]:
        """Create SampleMetadata for each sample id.

        If group_mapping is provided, uses it to assign groups.
        Otherwise, extracts group from sample ID using pattern matching.
        """
        # Build reverse mapping: sample_id -> group_name
        sample_to_group = {}
        if self.group_mapping:
            for group_name, samples in self.group_mapping.items():
                for sample_id in samples:
                    sample_to_group[sample_id] = group_name

        def extract_group(sample_id: str) -> str:
            # First check explicit mapping
            if sample_id in sample_to_group:
                return sample_to_group[sample_id]

            # Fall back to pattern extraction
            if not sample_id or not sample_id.strip():
                return "unknown"
            match = re.match(r"^(\D+)", sample_id)
            if match:
                group = match.group(1).strip("_")
                return group if group else "unknown"
            return "unknown"

        if labels:
            # Use labels to assign groups if available
            label_map = {sid: lbl for sid, lbl in zip(sample_ids, labels[1:]) if lbl}
            samples = [
                SampleMetadata(
                    sample_id=sid,
                    group=extract_group(sid),
                    label=label_map.get(sid)
                )
                for sid in sample_ids
            ]
        else:
            samples = [
                SampleMetadata(sample_id=sid, group=extract_group(sid))
                for sid in sample_ids
            ]
        if self.group_mapping:
            logger.info(
                f"Applied group_mapping: {len(self.group_mapping)} groups, "
                f"mapped {len(sample_to_group)} samples"
            )

        return samples

    def extract_quantified_lipids(
        self, rows: List[Dict], name_col: str, sample_ids: List[str], column_info: Optional[Dict[str, Any]] = None
    ) -> List[QuantifiedLipid]:
        """Extract QuantifiedLipid objects from CSV rows."""
        logger.info(f"Extracting quantified lipids using name_col='{name_col}' and {sample_ids} samples")
        quantified = []
        skipped_rows = 0
        empty_columns = []
        non_numeric_columns = []

        if column_info is not None:
            empty_columns = column_info.get("empty_columns", [])
            non_numeric_columns = [col for col, ctype in column_info.get("column_types", {}).items() if ctype != "numeric"]

        
        for row_idx, row in enumerate(rows, start=1):
            lipid_name = (row.get(name_col) or "").strip()
            if not lipid_name:
                skipped_rows += 1
                logger.info(f"Skipping row {row_idx}: empty lipid name")
                continue


            values = {}
            skipped_values = 0
            for sid in sample_ids:
                if sid in empty_columns or sid in non_numeric_columns:
                    skipped_values += 1
                    continue
                raw = (row.get(sid) or "").strip()
                if raw == "":
                    skipped_values += 1
                    continue
                try:
                    values[sid] = float(raw)
                except ValueError:
                    skipped_values += 1
                    # logger.warning(
                    #     f"Non-numeric value for sample {sid} at row {row_idx}: {raw!r}"
                    # )
                    continue
            if values:
                quantified.append(QuantifiedLipid(input_name=lipid_name, values=values))
            else:
                skipped_rows += 1
                logger.info(f"Skipping row {row_idx}: no valid values found")
        if skipped_rows > 0:
            logger.info(f"Total skipped rows: {skipped_rows}")
        return quantified

    def annotate_lipids_with_refmet(self, quantified: List[Any]) -> None:
        """Annotate QuantifiedLipid objects with RefMet data."""
        try:
            # Extract lipid names
            lipid_names = [q.input_name for q in quantified]

            # Call RefMet API to get results
            refmet_results = RefMet.validate_metabolite_names(lipid_names)
            logger.info(f"RefMet returned {len(refmet_results)} results")

            # Apply results to quantified lipids
            for q, result in zip(quantified, refmet_results):
                q.standardized_name = result.standardized_name
                # record that the standardized name came from RefMet when present
                try:
                    if getattr(result, "standardized_name", None):
                        q.standardized_by = "RefMet"
                except Exception:
                    pass

                q.lm_id = result.lm_id
                # if RefMet supplied the lm_id, record that source
                try:
                    if getattr(result, "lm_id", None):
                        q.lm_id_found_by = "RefMet"
                except Exception:
                    pass
                q.sub_class = result.sub_class
                q.formula = result.formula
                q.mass = result.exact_mass
                q.super_class = result.super_class
                q.main_class = result.main_class
                q.chebi_id = result.chebi_id
                q.kegg_id = result.kegg_id
                q.refmet_id = result.refmet_id
        except Exception:
            logger.exception(
                "RefMet annotation failed; continuing without standardized names"
            )

    def fill_missing_lm_ids_from_lmsd(
        self, quantified: Optional[List[Any]] = None, use_standardized_name: bool = True
    ) -> int:
        """Fill missing `lm_id` fields on QuantifiedLipid objects using LMSD.

        Args:
            quantified: Optional list of QuantifiedLipid objects to operate on. If
                not provided, uses `self.dataset.lipids`.
            use_standardized_name: If True, prefer `standardized_name` when calling
                the LMSD API; fall back to `input_name` when absent.

        Returns:
            Number of lipids updated with an `lm_id`.
        """
        if quantified is None:
            if self.dataset is None:
                logger.info("No dataset available to update lm_ids")
                return 0
            quantified = self.dataset.lipids

        # Collect indices and names for lipids missing lm_id
        missing_indices: List[int] = []
        query_names: List[str] = []
        for i, q in enumerate(quantified):
            current = getattr(q, "lm_id", None)
            if not current:
                name = None
                if use_standardized_name:
                    name = getattr(q, "standardized_name", None)
                if not name:
                    name = getattr(q, "input_name", None)
                # ensure we have a name to query
                if name:
                    missing_indices.append(i)
                    query_names.append(name)

        if not query_names:
            logger.info("No missing lm_id entries to update via LMSD")
            return 0

        try:
            logger.info(f"Querying LMSD for {len(query_names)} names to fill missing lm_id fields")
            resp = LMSD.get_lm_ids_by_name(query_names)
        except Exception:
            logger.info("LMSD lookup failed")
            return 0

        # If LMSD returned an error dict, log and exit
        if isinstance(resp, dict) and resp.get("error"):
            logger.info(f"LMSD returned error: {resp.get('error')}")
            return 0

        if not isinstance(resp, list):
            logger.info(f"Unexpected LMSD response type: {type(resp)}")
            return 0

        updated = 0
        # Map results back to quantified list
        for idx, item in zip(missing_indices, resp):
            if not isinstance(item, dict):
                continue
            lm_id = item.get("lm_id")
            if lm_id:
                try:
                    quantified[idx].lm_id = lm_id
                    # Optionally record which field matched (if model supports it)
                    if hasattr(quantified[idx], "matched_field"):
                        setattr(quantified[idx], "matched_field", item.get("matched_field"))
                    # record that this lm_id was found via LMSD
                    try:
                        quantified[idx].lm_id_found_by = "LMSD"
                    except Exception:
                        pass
                except Exception:
                    logger.exception(f"Failed to set lm_id on lipid at index {idx}")
                    continue
                updated += 1

        logger.info(f"Updated {updated} lm_id fields using LMSD")
        return updated

    def run_lmsd_fill_and_report(self, dataset: Optional[Any] = None) -> int:
        """Run LMSD fill for `dataset` and return number updated.

        This wraps `fill_missing_lm_ids_from_lmsd` and records/report changes
        so it can be called from external code (API/CLI).
        """
        # Accept either an explicit dataset or use the manager's current dataset
        if dataset is None:
            dataset = getattr(self, "dataset", None)

        # Record pre-existing lm_id state if we have a dataset
        pre_lm: Dict[str, Any] = {}
        if dataset is not None and getattr(dataset, "lipids", None) is not None:
            pre_lm = {q.input_name: q.lm_id for q in dataset.lipids}

        # If we have an explicit dataset, pass its lipids list to the fill helper
        if dataset is not None and getattr(dataset, "lipids", None) is not None:
            updated_count = self.fill_missing_lm_ids_from_lmsd(quantified=dataset.lipids)
        else:
            updated_count = self.fill_missing_lm_ids_from_lmsd()

        # Logging/reporting: mirror previous CLI behavior
        logger = logging.getLogger(__name__)
        logger.info(f"LMSD fill completed: {updated_count} updated")
        if updated_count and dataset is not None and getattr(dataset, "lipids", None) is not None:
            for q in dataset.lipids:
                prev = pre_lm.get(q.input_name)
                if not prev and q.lm_id:
                    logger.info(
                        f"  {q.input_name} -> {q.lm_id} (matched_field={getattr(q, 'matched_field', None)})"
                    )

        return updated_count

    def print_report(self) -> None:
        """Print the most recent validation report if available."""
        if not self.validation_report:
            logger.info(
                "No validation report available; run with validate_data=True to generate one"
            )
            return
        self.validation_report.print_report()

    def dataset_dict(self) -> Dict[str, Any]:
        """Serialize the dataset to plain dict for JSON output or downstream analysis."""
        if self.dataset is None:
            return {}
        return (
            self.dataset.model_dump()
            if hasattr(self.dataset, "model_dump")
            else self.dataset.dict()
        )

    # small helper to compute dataframe for quick analysis
    def dataset_as_dataframe(self) -> pd.DataFrame:
        """Return a pandas DataFrame with lipids as rows and samples as columns."""
        if self.dataset is None:
            return pd.DataFrame()
        records = []
        for q in self.dataset.lipids:
            rec = {"lipid": q.input_name}
            rec.update(q.values)
            records.append(rec)
        df = pd.DataFrame.from_records(records).set_index("lipid")
        return df

    def add_lipid_species(self, lipid: Any) -> None:
        """Add a lipid species to the manager (legacy helper used by tests).

        Keeps a simple list `lipid_species` and also appends to dataset.lipids when present.
        """
        try:
            self.lipid_species.append(lipid)
        except Exception:
            self.lipid_species = getattr(self, "lipid_species", []) + [lipid]

        if self.dataset is None:
            # create minimal dataset if necessary
            try:
                self.dataset = LipidDataset(samples=[], lipids=[lipid])
            except Exception:
                pass
        else:
            try:
                self.dataset.lipids.append(lipid)
            except Exception:
                pass

    # STATISTICS AND COMPARISONS

    def get_group_statistics(self) -> Dict[str, Dict[str, Any]]:
        """Calculate statistics for each group across all lipids.

        Returns:
            Dict mapping group names to their statistics including:
            - sample_count: number of samples in group
            - lipid_coverage: how many lipids have data for this group
            - mean_values: dict of lipid -> mean value in this group
            - std_values: dict of lipid -> std dev in this group
        """
        if self.dataset is None:
            return {}

        import numpy as np

        # Group samples by their group attribute
        groups: Dict[str, List[SampleMetadata]] = {}
        for sample in self.dataset.samples:
            groups.setdefault(sample.group, []).append(sample)

        group_stats = {}
        for group_name, samples in groups.items():
            sample_ids = [s.sample_id for s in samples]

            lipid_means = {}
            lipid_stds = {}
            lipid_coverage = 0

            for lipid in self.dataset.lipids:
                # Extract values for this group's samples
                group_values = [
                    lipid.values.get(sid) for sid in sample_ids if sid in lipid.values
                ]

                if group_values:
                    lipid_coverage += 1
                    lipid_means[lipid.input_name] = float(np.mean(group_values))
                    lipid_stds[lipid.input_name] = float(np.std(group_values))

            group_stats[group_name] = {
                "sample_count": len(samples),
                "lipid_coverage": lipid_coverage,
                "mean_values": lipid_means,
                "std_values": lipid_stds,
            }

        return group_stats


if __name__ == "__main__":
    """Quick demo showcasing CSVIngestion + DataValidator + DataManager."""
    import sys

    logging.basicConfig(level=logging.INFO)
    csv_path = (
        Path(__file__).parents[3]
        / "tests"
        / "data"
        / "inputs"
        / "file_structure_negative.csv"
    )
    print(f"Using CSV: {csv_path}")

    # --- CSV ingestion -----------------------------------------------------
    ingestion = CSVIngestion()
    try:
        raw_df = ingestion.read_csv(csv_path)
    except Exception as exc:
        logger.exception(f"CSV ingestion failed: {exc}")
        sys.exit(2)

    print(
        f"Ingested {raw_df.row_count} rows x {raw_df.column_count} columns "
        f"(format={raw_df.format_type.value})"
    )
    if raw_df.fieldnames:
        preview = ", ".join(raw_df.fieldnames[:5])
        if len(raw_df.fieldnames) > 5:
            preview += ", ..."
        print(f"Columns: {preview}")

    # --- Data validation ---------------------------------------------------
    validator = DataValidator()
    validation_report = validator.validate(raw_df)
    status = "PASSED" if validation_report.passed else "FAILED"
    print(
        f"Validation status: {status} "
        f"({len(validation_report.issues)} issues, warnings={validation_report.has_warnings})"
    )
    if validation_report.issues:
        sample_issues = validation_report.issues[:3]
        print("Sample issues:")
        for issue in sample_issues:
            print(f"  - {issue}")
        remaining = len(validation_report.issues) - len(sample_issues)
        if remaining:
            print(f"  ...and {remaining} more")

    # --- Dataset processing ------------------------------------------------
    mgr = DataManager(validate_data=True)
    try:
        ds = mgr.process_csv(csv_path)
    except Exception as exc:
        logger.exception(f"Demo failed to process CSV: {exc}")
        sys.exit(2)

    print("\nDataset summary:")
    print(f"  samples: {len(ds.samples)}")
    print(f"  lipids: {len(ds.lipids)}")

    # Group analysis
    print("\n" + "=" * 60)
    print("GROUP-LEVEL ANALYSIS")
    print("=" * 60)

    # Show group statistics
    group_stats = mgr.get_group_statistics()
    print(f"\nFound {len(group_stats)} groups:")
    for group_name, stats in group_stats.items():
        print(
            f"  {group_name}: {stats['sample_count']} samples, {stats['lipid_coverage']} lipids with data"
        )
