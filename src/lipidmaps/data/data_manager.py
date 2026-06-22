import csv
import logging
import re
import tempfile
from typing import List, Dict, Any, Union, Optional
from pathlib import Path
import pandas as pd
# import networkx as nx
from pydantic import ConfigDict, Field, field_validator
from .models.base import LipidmapsBaseModel
from .biopan_exporter import BioPANExporter
from .biopan_pathway_exporter import BioPANPathwayExporter
# import the data models we will produce
from .models.sample import SampleMetadata, QuantifiedLipid, LipidDataset, LipidAnnotation, SampleConditions
from .models.refmet import RefMet
from .models.lmsd import LMSD, LMSDResult
# from .reaction_checker import ReactionChecker, ReactionData
from .models.reaction import ReactionData

# import new ingestion and validation modules
from .ingestion.csv_reader import CSVIngestion, CSVFormat
from .validation.data_validator import DataValidator, ValidationReport


logger = logging.getLogger(__name__)

class DataManager(LipidmapsBaseModel):

    """DataManager: reads CSVs into LipidDataset objects.
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

    # use_refmet: bool = Field(default=True, description="Whether to use RefMet for annotation during CSV processing")
    use_refmet: bool = Field(default=True, description="Whether to use RefMet for annotation during CSV processing")
    use_headgroups: bool = Field(default=True, description="Whether to use headgroup mapping for filling missing LM IDs")
    fetch_reactions: bool = Field(default=True, description="Whether to fetch reactions by LM ID after processing CSV")
    taxonomy_group: Optional[str] = Field(default="all", description="Taxonomy group filter for reaction fetching (e.g. 'bacteria', 'mammalia', 'all')")
    transpose_file: bool = Field(
        default=False,
        description="If True, transpose the input CSV before ingestion (useful when lipids are columns and samples are rows)."
    )
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

    model_config = ConfigDict(arbitrary_types_allowed=True)

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

        Uses CSVIngestion for reading and DataValidator for quality checks.

        Args:
            csv_path: Path to CSV file

        Returns:
            LipidDataset with processed data
        """
        csv_path = Path(csv_path)
        logger.info(f"Loading CSV file: {csv_path}")

        # Optional manual transpose: write a transposed temp CSV and use that file
        if self.transpose_file:
            try:
                logger.info("transpose_file=True: transposing CSV so lipids become rows")
                # read with first column as index if present, otherwise fall back to default read
                try:
                    df = pd.read_csv(csv_path, header=0, index_col=0)
                except Exception:
                    df = pd.read_csv(csv_path, header=0)
                transposed = df.T.reset_index()
                tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
                transposed.to_csv(tmp.name, index=False)
                csv_path = Path(tmp.name)
                logger.info(f"Transposed CSV written to temporary file: {csv_path}")
            except Exception:
                logger.exception("Failed to transpose CSV; proceeding with original file")


        # Use CSVIngestion to read file
        ingestion = CSVIngestion(has_labels=self.has_labels)
        raw_df = ingestion.read_csv(csv_path, format_type=self.csv_format)
        column_info = ingestion.get_column_info(raw_df)
        logger.info(
            f"CSV ingested: {raw_df.row_count} rows x {raw_df.column_count} columns."
        )

        # Process the raw data
        if raw_df.is_empty():
            ds = LipidDataset(samples=[], lipids=[])
            self.dataset = ds
            return ds

        # Determine lipid name column
        name_col = self._resolve_lipid_column(raw_df.fieldnames)

        # Determine sample columns
        sample_names = self._resolve_sample_columns(raw_df.fieldnames, name_col)
        labels = raw_df.labels if hasattr(raw_df, "labels") else []
        # Create sample metadata with group mapping if provided
        samples_meta = self.extract_sample_metadata(sample_names, labels=labels)

        # Extract quantified lipids
        quantified = self.extract_quantified_lipids(raw_df.rows, name_col, sample_names, column_info)

        # If refmet is enabled, annotate lipids with refmet results before creating the dataset, so that standardized names and lm_ids are available on the QuantifiedLipid objects within the dataset from the start
        # Default is true since it provides standardized names and can improve LMSD matching downstream, but can be disabled if users want to skip that step or handle annotation separately
        refmet_failed = False
        if self.use_refmet:
            refmet_success = self.annotate_lipids_with_refmet(quantified)
            refmet_failed = not refmet_success

        dataset = LipidDataset(samples=samples_meta, lipids=quantified, column_info=column_info, refmet_failed=refmet_failed)
        self.dataset = dataset
        logger.info(
            f"Created LipidDataset: {len(samples_meta)} samples, {len(quantified)} lipids"
        )

        if self.use_headgroups:
            headgroup_updates = self.dataset.fill_generic_lm_ids_from_headgroups()
            logger.info(f"Filled missing LM IDs using headgroup mapping: {headgroup_updates} updated")

        if self.fetch_reactions:
            reaction_updates = self.dataset.fetch_reactions_by_lm_id(taxonomy_group=self.taxonomy_group)
            try:
                count = len(reaction_updates) if reaction_updates is not None else 0
            except Exception:
                count = 0
            logger.info(f"Fetched reactions for lipids: {count} reactions")

        if self.validate_data:
            validator = DataValidator()
            self.validation_report = validator.validate(raw_df)
            
            if not self.validation_report.passed:
                logger.warning(
                    f"Validation found {len(self.validation_report.issues)} issues"
                )
            # After processing, attach the validation report to the dataset for downstream use
            self.dataset.validation_report = self.validation_report

        # Below we'll fetch lmid details for all lipids with an lm_id and annotate them, which will be used by downstream analysis and reporting
        lm_ids = dataset.list_lm_ids()
        molecules = LMSD.get_molecules_by_lm_id(lm_ids)  # Fetch details for first 5 LM IDs as a check
        logger.info(f"LMSD molecules fetched: {molecules[:5] if isinstance(molecules, list) else molecules}")
        self.annotate_lipids_with_lmsd_details(dataset=dataset, molecules=molecules)
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

    def _resolve_sample_columns(self, fieldnames: List[str], lipid_col: str) -> List[str]:
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


    def extract_sample_metadata(self, sample_names: List[str], labels: Optional[List[str]]=None) -> List[SampleMetadata]:
        """Create SampleMetadata for each sample id.

        If group_mapping is provided, uses it to assign groups.
        Otherwise, extracts group from sample ID using pattern matching.
        """
        # Build reverse mapping: sample_name -> group_name
        sample_to_group = {}
        if self.group_mapping:
            for group_name, samples in self.group_mapping.items():
                for sample_name in samples:
                    sample_to_group[sample_name] = group_name

        def extract_group(sample_name: str) -> str:
            # First check explicit mapping
            if sample_name in sample_to_group:
                return sample_to_group[sample_name]

            # Fall back to pattern extraction
            if not sample_name or not sample_name.strip():
                return "unknown"
            match = re.match(r"^(\D+)", sample_name)
            if match:
                group = match.group(1).strip("_")
                return group if group else "unknown"
            return "unknown"

        if labels:
            # Use labels to assign groups if available
            label_map = {sid: lbl for sid, lbl in zip(sample_names, labels[1:]) if lbl}
            samples = [
                SampleMetadata(
                    sample_name=sid,
                    group=extract_group(sid),
                    label=label_map.get(sid)
                )
                for sid in sample_names
            ]
        else:
            samples = [
                SampleMetadata(sample_name=sid, group=extract_group(sid))
                for sid in sample_names
            ]
        if self.group_mapping:
            logger.info(
                f"Applied group_mapping: {len(self.group_mapping)} groups, "
                f"mapped {len(sample_to_group)} samples"
            )

        return samples

    def extract_quantified_lipids(self, rows: List[Dict], name_col: str, sample_names: List[str], column_info: Optional[Dict[str, Any]] = None) -> List[QuantifiedLipid]:
        """Extract QuantifiedLipid objects from CSV rows."""
        logger.info(f"Extracting quantified lipids using name_col='{name_col}' and {len(sample_names)} samples")
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
            for name in sample_names:
                if name in empty_columns or name in non_numeric_columns:
                    skipped_values += 1
                    continue
                raw = (row.get(name) or "").strip()
                if raw == "":
                    skipped_values += 1
                    continue
                try:
                    values[name] = float(raw)
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

    def annotate_lipids_with_refmet(self, quantified: List[Any]) -> bool:
        """Annotate QuantifiedLipid objects with RefMet data.
        
        Returns:
            True if annotation succeeded, False if it failed.
        """
        try:
            # Extract lipid names
            lipid_names = [q.input_name for q in quantified]

            # Call RefMet API to get results
            refmet_results = RefMet.validate_metabolite_names(lipid_names)
            if isinstance(refmet_results, dict):
                logger.error(
                    "RefMet returned an error payload; continuing without standardized names: %s",
                    refmet_results.get("error", "unknown error"),
                )
                return False
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
                
                # Create annotation object with classification and IDs
                q.annotation = LipidAnnotation(
                    sub_class=result.sub_class,
                    main_class=result.main_class,
                    super_class=result.super_class,
                    chebi_id=result.chebi_id,
                    kegg_id=result.kegg_id,
                    refmet_id=result.refmet_id,
                    formula=result.formula,
                    mass=result.exact_mass,
                )
            return True
        except Exception:
            logger.exception(
                "RefMet annotation failed; continuing without standardized names"
            )
            return False

    #NOTE: This method is separate from the main CSV processing flow since we get lm_id results 
    # from RefMet annotation, and later we fill generic lm_id's from headgroups.
    #NOTE: This method updates input list of quantified object and if no dataset is provided, 
    # it will update the manager's current dataset.lipids list if available. 
    # This allows it to be used as a standalone method for filling lm_id details on any list of 
    # QuantifiedLipid objects, or as part of the main CSV processing flow.

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


    def fill_generic_lm_ids_from_headgroups(self, dataset: Optional[LipidDataset] = None) -> int:
        """
        Fill missing lm_id fields on QuantifiedLipid objects using headgroup mapping from headgroups.py.
        Args:
            dataset: Optional LipidDataset to operate on. If not provided, uses self.dataset.
        Returns:
            Number of lipids updated with an lm_id.
        """
        if dataset is None:
            dataset = self.dataset
        if dataset is None or not hasattr(dataset, "lipids"):
            logger.warning("No dataset or lipids to fill with headgroup mapping.")
            return 0
        # Call the new method on LipidDataset
        return dataset.fill_generic_lm_ids_from_headgroups()

    def annotate_lipids_with_lmsd_details(self, dataset: Optional[LipidDataset] = None, molecules: Optional[List[LMSDResult]] = None) -> int:
        """Annotate quantified lipids with additional LMSD details (abbrev, generic_lm_id).

        This collects LM IDs present on the provided `quantified` list (or
        manager's dataset), queries `LMSD.get_molecules_by_lm_ids`, and maps
        returned fields (`abbrev`, `generic_lm_id`, `abbrev_chains`, `smiles`)
        back onto the QuantifiedLipid objects.

        Returns number of lipids updated.
        """

        # If called with no explicit list, try to use the manager's dataset
        if dataset is None:
            if self.dataset is None:
                logger.info("No dataset available to annotate with LMSD details")
                return 0
            # Delegate to LipidDataset method if available
            dataset = self.dataset


        if not molecules:
            logger.info("LMSD output doesn't exist yet; no molecules provided")
            return 0

        mapping: Dict[Optional[str], LMSDResult] = {}
        for m in (molecules or []):
            if m is None:
                continue
            if isinstance(m, LMSDResult):
                lm = m.lm_id
                mapping[lm] = m
            else:
                logger.warning(f"Unexpected molecule type in LMSD results: {type(m)}")

        updated = 0
        for q in dataset.lipids:
            lm = getattr(q, "lm_id", None)
            if not lm or lm not in mapping:
                continue
            item = mapping[lm]
            try:
                # if hasattr(q, "abbrev") and item.get("abbrev"):
                #     q.abbrev = item.get("abbrev")
                if hasattr(q, "generic_lm_id") and item.generic_lm_id:
                    q.generic_lm_id = item.generic_lm_id
                # if hasattr(q, "abbrev_chains") and item.get("abbrev_chains"):
                #     q.abbrev_chains = item.get("abbrev_chains")
                # if hasattr(q, "smiles") and item.get("smiles"):
                #     q.smiles = item.get("smiles")
                updated += 1
            except Exception:
                logger.exception(f"Failed to annotate lipid {getattr(q, 'input_name', None)} with LMSD details")

        logger.info(f"Annotated {updated} lipids with LMSD details")
        return updated

    # Reaction fetching has been moved to `LipidDataset.fetch_reactions_by_lm_id()`
    # to keep dataset-level operations within the dataset object.
    def annotate_lipids_with_reactions(self, dataset: LipidDataset, reactions: Optional[list[ReactionData]] = None) -> None:
        """
        For each QuantifiedLipid in the dataset, find all reactions where the lipid's lm_id or generic_lm_id is a reactant or product (compound_lm_id),
        and update the QuantifiedLipid's 'reactions' field with a list of unique ReactionData objects. Also, assign all unique reactions to dataset.reactions.
        """
        if reactions is None:
            logger.warning("No reactions provided.")
            return

        if dataset is None or not hasattr(dataset, "lipids"):
            logger.warning("No dataset or lipids to annotate with reactions.")
            return

        # Deduplicate reactions by reaction_id (or id)
        reaction_dict = {}
        for reaction in reactions:
            reaction_id = getattr(reaction, "reaction_id", None)
            if reaction_id is None:
                reaction_id = getattr(reaction, "id", None)
            if reaction_id is not None:
                reaction_id = str(reaction_id)
                if reaction_id not in reaction_dict:
                    reaction_dict[reaction_id] = reaction
        # Assign all unique reactions to dataset.reactions
        dataset.reactions = list(reaction_dict.values())

        # For each lipid, collect matching unique reactions (by id)
        for lipid in dataset.lipids:
            lipid_ids = set()
            if getattr(lipid, "lm_id", None):
                lipid_ids.add(lipid.lm_id)
            if getattr(lipid, "generic_lm_id", None):
                lipid_ids.add(lipid.generic_lm_id)
            matched_ids = set()
            matched_reactions = []
            for reaction in dataset.reactions:
                for role in ["reactants", "products"]:
                    items = getattr(reaction, role, [])
                    for item in items:
                        compound_lm_id = None
                        if isinstance(item, dict):
                            compound_lm_id = item.get("compound_lm_id")
                        elif hasattr(item, "compound_lm_id"):
                            compound_lm_id = getattr(item, "compound_lm_id", None)
                        if compound_lm_id and compound_lm_id in lipid_ids:
                            rid = getattr(reaction, "reaction_id", None) or getattr(reaction, "id", None)
                            if rid is not None:
                                rid = str(rid)
                                if rid not in matched_ids:
                                    matched_reactions.append(reaction)
                                    matched_ids.add(rid)
                            break  # Only need to add once per reaction
                    else:
                        continue
                    break
            lipid.reactions = matched_reactions if matched_reactions else None

    def print_report(self) -> None:
        """Print the most recent validation report if available."""
        if not self.validation_report:
            logger.info(
                "No validation report available; run with validate_data=True to generate one"
            )
            return
        self.validation_report.print_report()


    def summarize_csv(self, dataset: LipidDataset) -> Dict[str, Any]:
        """
        Summarize a quantitative CSV file and return a JSON-style dict with key stats and structure.
        Includes: number of samples, number of lipids, sample IDs, lipid names, and column info.
        """
        summary = {
            "num_samples": len(dataset.list_sample_names()),
            "num_lipids": len(dataset.list_lipid_names()),
            "num_lm_ids": len(dataset.list_lipids_with_lmid()),
            # "sample_names": dataset.list_sample_names() if hasattr(dataset, 'list_sample_names') else [s.sample_name for s in dataset.samples],
            # "lipid_names": dataset.list_lipid_names() if hasattr(dataset, 'list_lipid_names') else [l.input_name for l in dataset.lipids],
            "column_info": getattr(dataset, 'column_info', None),
        }

        return summary

    def get_biopan_exporter(self, dataset: Optional[LipidDataset] = None) -> BioPANExporter:
        return BioPANExporter(dataset=dataset or self.dataset)

    def get_sample_conditions(
        self,
        dataset: Optional[LipidDataset] = None,
    ) -> SampleConditions:
        resolved_dataset = dataset or self.dataset
        if resolved_dataset is None:
            raise ValueError("Sample conditions require a populated dataset")
        return resolved_dataset.get_sample_conditions()

    def set_sample_conditions(
        self,
        conditions: Union[SampleConditions, Dict[str, str]],
        dataset: Optional[LipidDataset] = None,
        *,
        strict: bool = True,
    ) -> SampleConditions:
        resolved_dataset = dataset or self.dataset
        if resolved_dataset is None:
            raise ValueError("Sample conditions require a populated dataset")
        return resolved_dataset.set_sample_conditions(conditions, strict=strict)

    def get_biopan_pathway_exporter(self, dataset: Optional[LipidDataset] = None) -> BioPANPathwayExporter:
        return BioPANPathwayExporter(dataset=dataset or self.dataset)

    def build_biopan_summary(
        self,
        dataset: Optional[LipidDataset] = None,
        lipidlynxx: str = "no",
    ) -> Dict[str, Any]:
        """Build the subset of BioPAN summary data needed by the PHP frontend."""
        return self.get_biopan_exporter(dataset).build_summary(lipidlynxx=lipidlynxx)

    def build_biopan_msg1(
        self,
        dataset: Optional[LipidDataset] = None,
        lipidlynxx_error: bool = False,
    ) -> Dict[str, Any]:
        """Build BioPAN parse-stage status data used by the summary page."""
        return self.get_biopan_exporter(dataset).build_msg1(lipidlynxx_error=lipidlynxx_error)

    def build_biopan_msg2(self, dataset: Optional[LipidDataset] = None) -> Dict[str, Any]:
        """Build BioPAN processing-stage status data used by the pathway page."""
        return self.get_biopan_exporter(dataset).build_msg2()

    def export_biopan_display_files(
        self,
        output_path: Union[str, Path],
        dataset: Optional[LipidDataset] = None,
        lipidlynxx: str = "no",
        lipidlynxx_error: bool = False,
        include_msg1: bool = True,
        include_summary: bool = True,
        include_msg2: bool = True,
    ) -> Dict[str, str]:
        """Write BioPAN display JSON files into a session directory or its biopan subdir."""
        written_files = self.get_biopan_exporter(dataset).export_display_files(
            output_path=output_path,
            lipidlynxx=lipidlynxx,
            lipidlynxx_error=lipidlynxx_error,
            include_msg1=include_msg1,
            include_summary=include_summary,
            include_msg2=include_msg2,
        )
        logger.info("Exported BioPAN display files to %s", output_path)
        return written_files

    def export_biopan_reaction_files(
        self,
        output_path: Union[str, Path],
        disease_group: str,
        control_group: str,
        threshold: float = 0.05,
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        lazy: bool = False,
    ) -> Dict[str, str]:
        written_files = self.get_biopan_pathway_exporter(dataset).export_reaction_files(
            output_path=output_path,
            disease_group=disease_group,
            control_group=control_group,
            threshold=threshold,
            paired=paired,
            lazy=lazy,
        )
        logger.info("Exported BioPAN reaction files to %s", output_path)
        return written_files
    
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
            sample_names = [s.sample_name for s in samples]

            lipid_means = {}
            lipid_stds = {}
            lipid_coverage = 0

            for lipid in self.dataset.lipids:
                # Extract values for this group's samples
                group_values = [
                    lipid.values.get(name) for name in sample_names if name in lipid.values
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

    def selected(self, n: int = 10) -> List[QuantifiedLipid]:
        """Return the first `n` QuantifiedLipid objects from the current dataset.

        This is a small convenience wrapper used by CLI/debug helpers to inspect
        a subset of the dataset. If no dataset is present, returns an empty list.
        """
        if self.dataset is None or not hasattr(self.dataset, "lipids"):
            return []
        try:
            return list(self.dataset.lipids)[: max(0, int(n))]
        except Exception:
            # Defensive fallback
            return self.dataset.lipids[:n] if hasattr(self.dataset, "lipids") else []
