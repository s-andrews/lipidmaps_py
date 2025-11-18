import csv
import logging
import re
from typing import List, Tuple, Dict, Any, Union, Optional
from pathlib import Path
import pandas as pd

from pydantic import BaseModel, Field

# import the data models we will produce
from .models.sample import SampleMetadata, QuantifiedLipid, LipidDataset
from .models.refmet import RefMet

# import new ingestion and validation modules
from .ingestion.csv_reader import CSVIngestion, RawDataFrame, CSVFormat
from .validation.data_validator import DataValidator, ValidationReport


logger = logging.getLogger(__name__)


class DataManager(BaseModel):
    """Pydantic v2 DataManager: reads CSVs into LipidDataset objects.

    Now uses CSVIngestion for file reading and DataValidator for quality checks.

    Usage:
        mgr = DataManager()
        dataset = mgr.process_csv("tests/inputs/quantified_test_file.csv")
        
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

    model_config = {"arbitrary_types_allowed": True}

    def model_post_init(self, __context: dict) -> None:
        logger.info("Initialized DataManager (validation=%s)", self.validate_data)

    def process_csv(self, csv_path: Union[str, Path]) -> LipidDataset:
        """Read CSV and populate SampleMetadata, QuantifiedLipid and LipidDataset.
        
        Now uses CSVIngestion for reading and DataValidator for quality checks.
        
        Args:
            csv_path: Path to CSV file
            
        Returns:
            LipidDataset with processed data
        """
        csv_path = Path(csv_path)
        logger.info("Loading CSV file: %s", csv_path)
        
        # Use CSVIngestion to read file
        ingestion = CSVIngestion()
        raw_df = ingestion.read_csv(csv_path, format_type=self.csv_format)
        
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

        name_col = raw_df.fieldnames[0]
        sample_ids = [sid for sid in raw_df.fieldnames[1:] if sid and sid.strip()]
        samples_meta = self.extract_sample_metadata(sample_ids)
        quantified = self.extract_quantified_lipids(raw_df.rows, name_col, sample_ids)
        self.annotate_lipids_with_refmet(quantified)

        dataset = LipidDataset(samples=samples_meta, lipids=quantified)
        self.dataset = dataset
        logger.info(
            f"Created LipidDataset: {len(samples_meta)} samples, {len(quantified)} lipids"
        )
        return dataset

    def read_csv_rows(self, csv_path: Path) -> Tuple[List[Dict], List[str]]:
        """Read CSV and return rows and fieldnames.
        
        DEPRECATED: Use CSVIngestion directly instead.
        Kept for backward compatibility.
        """
        logger.debug("Using legacy read_csv_rows (consider using CSVIngestion)")
        with csv_path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
            fieldnames = reader.fieldnames or []
        return rows, fieldnames

    def extract_sample_metadata(self, sample_ids: List[str]) -> List[SampleMetadata]:
        """Create SampleMetadata for each sample id."""

        def extract_group(sample_id: str) -> str:
            if not sample_id or not sample_id.strip():
                return "unknown"
            match = re.match(r"^(\D+)", sample_id)
            if match:
                group = match.group(1).strip("_")
                return group if group else "unknown"
            return "unknown"

        return [
            SampleMetadata(sample_id=sid, group=extract_group(sid))
            for sid in sample_ids
        ]

    def extract_quantified_lipids(
        self, rows: List[Dict], name_col: str, sample_ids: List[str]
    ) -> List[QuantifiedLipid]:
        """Extract QuantifiedLipid objects from CSV rows."""
        quantified = []
        skipped_rows = 0
        for row_idx, row in enumerate(rows, start=1):
            lipid_name = (row.get(name_col) or "").strip()
            if not lipid_name:
                skipped_rows += 1
                logger.debug(f"Skipping row {row_idx}: empty lipid name")
                continue
            values = {}
            skipped_values = 0
            for sid in sample_ids:
                raw = (row.get(sid) or "").strip()
                if raw == "":
                    skipped_values += 1
                    continue
                try:
                    values[sid] = float(raw)
                except ValueError:
                    skipped_values += 1
                    logger.warning(
                        f"Non-numeric value for sample {sid} at row {row_idx}: {raw!r}"
                    )
                    continue
            if values:
                quantified.append(QuantifiedLipid(input_name=lipid_name, values=values))
                logger.debug(
                    f"Quantified lipid: {lipid_name}, number of values: {len(values)}"
                )
                logger.debug(
                    f"Skipped rows: {skipped_rows}, skipped values: {skipped_values}"
                )
            else:
                skipped_rows += 1
                logger.debug(f"Skipping row {row_idx}: no valid values found")
        if skipped_rows > 0:
            logger.debug(f"Total skipped rows: {skipped_rows}")
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
                q.lm_id = result.lm_id
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
    """Quick demo: run this module directly to parse the example CSV and print a short summary."""
    import sys

    logging.basicConfig(level=logging.INFO)
    mgr = DataManager()
    csv_path = Path(__file__).parents[1] / "tests" / "inputs" / "biopan_test_input.csv"
    print(f"Using CSV: {csv_path}")
    try:
        ds = mgr.process_csv(csv_path)
    except Exception as exc:
        logger.exception("Demo failed to process CSV: %s", exc)
        sys.exit(2)

    print("\nDataset summary:")
    print(f"  samples: {len(ds.samples)}")
    print(f"  lipids: {len(ds.lipids)}")

    # Group analysis
    print("\n" "=" * 60)
    print("GROUP-LEVEL ANALYSIS")
    print("=" * 60)

    # Show group statistics
    group_stats = mgr.get_group_statistics()
    print(f"\nFound {len(group_stats)} groups:")
    for group_name, stats in group_stats.items():
        print(
            f"  {group_name}: {stats['sample_count']} samples, {stats['lipid_coverage']} lipids with data"
        )
