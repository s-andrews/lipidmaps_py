"""
Data import functions for LIPID MAPS

Provides high-level API for importing and processing lipid quantification data.
Built on top of DataManager for robust CSV processing and RefMet integration.
"""

import logging
from typing import Union, Optional, List, Dict, Any
from pathlib import Path

from pydantic import BaseModel, Field, computed_field, ConfigDict

from .data.data_manager import DataManager
from .data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
from .data.models.reaction import ReactionData, CompoundComponent
from .data.quantitation import (
    QuantitationAnalyzer,
    QuantitationConfig,
    NormalizationMethod,
    QuantitationUnit,
    create_analyzer,
)

logger = logging.getLogger(__name__)


class LipidData(BaseModel):
    """High-level interface for lipid data imported from CSV files.

    This class wraps DataManager and LipidDataset to provide a simple API
    for accessing imported lipid data, with backward compatibility for
    legacy code.

    Attributes:
        dataset: The underlying LipidDataset object
        manager: The DataManager instance used for processing
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    dataset: LipidDataset = Field(
        default_factory=lambda: LipidDataset(samples=[], lipids=[])
    )
    manager: Optional[DataManager] = Field(default=None)

    def model_post_init(self, __context: Any) -> None:
        """Initialize manager if not provided."""
        if self.manager is None:
            object.__setattr__(self, "manager", DataManager(dataset=self.dataset))

    @computed_field  # type: ignore[misc]
    @property
    def lipids_list(self) -> List[QuantifiedLipid]:
        """Return list of QuantifiedLipid objects."""
        return self.dataset.lipids

    @computed_field  # type: ignore[misc]
    @property
    def failed_lipids(self) -> List[str]:
        """Return list of lipid names that failed to import or annotate."""
        # Lipids without standardized names could be considered "failed"
        return [
            lipid.input_name
            for lipid in self.dataset.lipids
            if lipid.standardized_name is None
        ]

    @computed_field  # type: ignore[misc]
    @property
    def sample_names(self) -> List[str]:
        """Return list of sample IDs."""
        return [s.sample_name for s in self.dataset.samples]

    def successful_import_count(self) -> int:
        """Return count of successfully imported lipids."""
        return len(self.dataset.lipids)

    def failed_import_count(self) -> int:
        """Return count of failed lipid imports."""
        return len(self.failed_lipids)

    def failed_import_names(self) -> List[str]:
        """Return list of names that failed to import."""
        return self.failed_lipids

    def samples(self) -> List[str]:
        """Return list of sample names."""
        return self.sample_names

    def lipids(self) -> List[QuantifiedLipid]:
        """Return list of QuantifiedLipid objects."""
        return self.dataset.lipids

    def get_value_for_lipid(
        self, lipid: Union[str, QuantifiedLipid], sample: str
    ) -> Optional[float]:
        """Get quantitation value for a specific lipid in a sample.

        Args:
            lipid: Lipid name (str) or QuantifiedLipid object
            sample: Sample ID

        Returns:
            Quantitation value or None if not found
        """
        if isinstance(lipid, str):
            # Find lipid by name
            lipid_obj = next(
                (
                    sample_lipid
                    for sample_lipid in self.dataset.lipids
                    if sample_lipid.input_name == lipid or sample_lipid.standardized_name == lipid
                ),
                None,
            )
            if lipid_obj is None:
                return None
            return lipid_obj.values.get(sample)
        else:
            return lipid.values.get(sample)

    def get_lipid_by_name(self, name: str) -> Optional[QuantifiedLipid]:
        """Get a QuantifiedLipid by input name or standardized name.

        Args:
            name: Lipid name to search for

        Returns:
            QuantifiedLipid object or None if not found
        """
        return next(
            (
                sample_lipid
                for sample_lipid in self.dataset.lipids
                if sample_lipid.input_name == name or sample_lipid.standardized_name == name
            ),
            None,
        )

    def get_lipids_by_class(self, lipid_class: str) -> List[QuantifiedLipid]:
        """Get all lipids belonging to a specific class.

        Args:
            lipid_class: Lipid class name (e.g., 'PC', 'TAG')

        Returns:
            List of QuantifiedLipid objects
        """
        return [
            lipid
            for lipid in self.dataset.lipids
            if lipid.annotation and (
                lipid.annotation.main_class == lipid_class or 
                lipid.annotation.sub_class == lipid_class
            )
        ]

    def get_lm_ids(self) -> List[str]:
        """Get all unique LIPID MAPS IDs from the dataset.

        Returns:
            List of LM IDs
        """
        return list(
            set(
                lipid.lm_id
                for lipid in self.dataset.lipids
                if lipid.lm_id and lipid.lm_id.startswith("LM")
            )
        )

    def as_dataframe(self):
        """Return dataset as pandas DataFrame with lipids as rows and samples as columns.

        Returns:
            pandas.DataFrame
        """
        return self.manager.dataset_as_dataframe()

    def get_group_statistics(self) -> Dict[str, Dict[str, Any]]:
        """Calculate statistics for each sample group across all lipids.

        Returns:
            Dict mapping group names to their statistics including:
            - sample_count: number of samples in group
            - lipid_coverage: how many lipids have data for this group
            - mean_values: dict of lipid -> mean value in this group
            - std_values: dict of lipid -> std dev in this group
        """
        return self.manager.get_group_statistics()

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the dataset to a dictionary.

        Returns:
            Dictionary representation of the dataset
        """
        if self.manager:
            return self.manager.dataset_dict()
        return self.model_dump()

    # TODO: Future methods for reactions integration
    def get_reactions(self, species: str = "human", complete: bool = True):
        """
        Retrieve reactions for imported lipids from LIPID MAPS API.

        Args:
            species: Species to filter reactions (default: "human")
            complete: Only include complete reactions (default: True)

        Returns:
            Reactions object

        Note:
            This method requires LIPID MAPS API integration (future implementation)
        """
        raise NotImplementedError(
            "Reactions integration not yet implemented. "
            "This will call LIPID MAPS API to retrieve reaction data."
        )

    def get_lipids_for_reaction_component(
        self, component: Union[str, CompoundComponent]
    ) -> List[QuantifiedLipid]:
        """Get all lipids matching a reaction component.

        Args:
            component: Component LM ID (str) or CompoundComponent object

        Returns:
            List of QuantifiedLipid objects matching the component
        """
        # If reactions integration is not enabled for this LipidData instance,
        # signal via NotImplementedError to match the higher-level API contract.
        if not getattr(self, "_reactions_available", True):
            raise NotImplementedError(
                "Reaction helper methods are not available for this LipidData instance."
            )

        if hasattr(component, "compound_lm_id"):
            comp_id = component.compound_lm_id
        else:
            comp_id = str(component)

        if not comp_id:
            return []

        comp_id_lower = comp_id.lower()
        return [
            l
            for l in self.dataset.lipids
            if (l.lm_id and l.lm_id.lower() == comp_id_lower)
            or (l.generic_lm_id and l.generic_lm_id.lower() == comp_id_lower)
        ]

    def get_value_for_reaction_component(
        self,
        component: Union[str, CompoundComponent],
        sample: Optional[Union[str, SampleMetadata]] = None,
        method: str = "sum",
    ) -> Optional[float]:
        """Get quantitation for a reaction component.

        A reaction component may map to multiple lipids in the dataset.
        This method aggregates values across all matching lipids.

        Args:
            component: Component LM ID (str) or CompoundComponent object
            sample: Sample name (str) or SampleMetadata object
            method: Aggregation method ('sum', 'mean', 'max', 'min')

        Returns:
            Aggregated quantitation value, or None if no matches
        """
        if not getattr(self, "_reactions_available", True):
            raise NotImplementedError(
                "Reaction helper methods are not available for this LipidData instance."
            )

        if sample is None:
            # Mirror previous behavior where callers may omit sample; raise TypeError
            # to indicate misuse rather than silently proceeding.
            raise TypeError("get_value_for_reaction_component() missing required 'sample' argument")

        return self.analyzer.get_value_for_reaction_component(component, sample, method)

    # =========================================================================
    # QUANTITATION ANALYSIS
    # =========================================================================

    @property
    def analyzer(self) -> QuantitationAnalyzer:
        """Get the quantitation analyzer for this dataset."""
        if not hasattr(self, "_analyzer") or self._analyzer is None:
            object.__setattr__(self, "_analyzer", create_analyzer(self.dataset))
        return self._analyzer

    @property
    def quantitation_config(self) -> QuantitationConfig:
        """Get the current quantitation configuration."""
        return self.analyzer.config

    def set_quantitation_config(
        self,
        unit: Optional[str] = None,
        method: Optional[str] = None,
        notes: Optional[str] = None,
    ) -> None:
        """Set quantitation configuration (unit, measurement method, notes).

        Args:
            unit: Unit of measurement (e.g., 'pmol', 'ng', 'area')
            method: Measurement method (e.g., 'LC-MS', 'GC-MS')
            notes: Additional notes about the quantitation
        """
        self.analyzer.set_config(unit=unit, method=method, notes=notes)

    def get_values(self, lipid_name: str) -> Optional[Dict[str, float]]:
        """Get quantified values for a lipid by name.

        Args:
            lipid_name: Lipid name (input_name or standardized_name)

        Returns:
            Dict mapping sample_name -> value, or None if not found
        """
        return self.dataset.get_values(lipid_name)

    # =========================================================================
    # NORMALIZATION METHODS
    # =========================================================================

    def normalize(
        self,
        method: Union[str, NormalizationMethod],
        internal_standard: Optional[str] = None,
        log_base: int = 2,
        log_offset: float = 1.0,
        total_lipid_scale: float = 1e6,
        in_place: bool = False,
    ) -> Dict[str, Dict[str, float]]:
        """Apply normalization to the dataset.

        Args:
            method: Normalization method ('total_lipid', 'internal_standard',
                   'log2', 'log10', 'median_center', 'zscore')
            internal_standard: Required if method is 'internal_standard'
            log_base: Base for log transformation (2 or 10)
            log_offset: Offset for log transformation to handle zeros
            total_lipid_scale: Scale factor for total lipid normalization
            in_place: If True, update lipid values in the dataset

        Returns:
            Dict mapping lipid_name -> {sample_name: normalized_value}
        """
        if isinstance(method, str):
            method = NormalizationMethod(method)
        return self.analyzer.apply_normalization(
            method=method,
            internal_standard=internal_standard,
            log_base=log_base,
            log_offset=log_offset,
            total_lipid_scale=total_lipid_scale,
            in_place=in_place,
        )

    def normalize_total_lipid(
        self, scale_factor: float = 1e6
    ) -> Dict[str, Dict[str, float]]:
        """Normalize by total lipid content per sample.

        Each sample's values are divided by the sum of all lipid values.

        Args:
            scale_factor: Factor to multiply normalized values (default: 1e6 for ppm)

        Returns:
            Dict mapping lipid_name -> {sample_name: normalized_value}
        """
        return self.analyzer.normalize_total_lipid(scale_factor)

    def normalize_log(
        self, base: int = 2, offset: float = 1.0
    ) -> Dict[str, Dict[str, float]]:
        """Apply log transformation.

        Args:
            base: Logarithm base (2 or 10)
            offset: Value added before log to handle zeros

        Returns:
            Dict mapping lipid_name -> {sample_name: log_value}
        """
        return self.analyzer.normalize_log(base, offset)

    def normalize_median_center(self) -> Dict[str, Dict[str, float]]:
        """Center values by subtracting sample median.

        Returns:
            Dict mapping lipid_name -> {sample_name: centered_value}
        """
        return self.analyzer.normalize_median_center()

    # =========================================================================
    # STATISTICAL ANALYSIS
    # =========================================================================

    def calculate_fold_change(
        self, group1: str, group2: str, log2: bool = True
    ) -> Dict[str, float]:
        """Calculate fold change between two groups.

        Args:
            group1: First group name (numerator)
            group2: Second group name (denominator/reference)
            log2: If True, return log2 fold change

        Returns:
            Dict mapping lipid_name -> fold_change
        """
        return self.analyzer.calculate_fold_change(group1, group2, log2)

    def calculate_pvalue(
        self,
        group1: str,
        group2: str,
        test: str = "ttest",
        paired: bool = False,
    ) -> Dict[str, float]:
        """Calculate p-values comparing two groups.

        Args:
            group1: First group name
            group2: Second group name
            test: Statistical test ('ttest' or 'mannwhitney')
            paired: If True, use paired test

        Returns:
            Dict mapping lipid_name -> p_value
        """
        return self.analyzer.calculate_pvalue(group1, group2, test, paired)

    def calculate_cv(
        self, by_group: bool = False
    ) -> Union[Dict[str, float], Dict[str, Dict[str, float]]]:
        """Calculate coefficient of variation for each lipid.

        CV = (std / mean) * 100

        Args:
            by_group: If True, calculate CV per group

        Returns:
            Dict mapping lipid_name -> CV (or group -> {lipid_name -> CV})
        """
        return self.analyzer.calculate_cv(by_group)

    def differential_analysis(
        self,
        group1: str,
        group2: str,
        fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05,
        test: str = "ttest",
    ) -> List[Dict[str, Any]]:
        """Perform differential analysis between two groups.

        Args:
            group1: First group name
            group2: Second group name (reference)
            fc_threshold: Absolute log2 fold change threshold
            pvalue_threshold: P-value significance threshold
            test: Statistical test to use

        Returns:
            List of dicts with lipid analysis results, sorted by p-value
        """
        return self.analyzer.differential_analysis(
            group1, group2, fc_threshold, pvalue_threshold, test
        )

    # =========================================================================
    # GROUP-LEVEL QUANTITATION
    # =========================================================================

    def get_group_means(self) -> Dict[str, Dict[str, float]]:
        """Get mean values per group for each lipid.

        Returns:
            Dict mapping group_name -> {lipid_name: mean_value}
        """
        return self.analyzer.get_group_means()

    def get_group_stds(self) -> Dict[str, Dict[str, float]]:
        """Get standard deviation per group for each lipid.

        Returns:
            Dict mapping group_name -> {lipid_name: std_value}
        """
        return self.analyzer.get_group_stds()

    def get_group_summary(self) -> Dict[str, Dict[str, Dict[str, float]]]:
        """Get comprehensive summary statistics per group.

        Returns:
            Dict mapping group_name -> {lipid_name: {mean, std, min, max, median, n}}
        """
        return self.analyzer.get_group_summary()

    def compare_groups(self, group1: str, group2: str) -> Dict[str, Dict[str, Any]]:
        """Compare two groups with full statistics.

        Args:
            group1: First group name
            group2: Second group name

        Returns:
            Dict mapping lipid_name -> {group1_mean, group1_std, group2_mean,
                                        group2_std, log2_fc, p_value, significant}
        """
        return self.analyzer.compare_groups(group1, group2)


def import_data(
    filename: Union[str, Path],
    lipid_col: Optional[Union[int, str]] = None,
    sample_cols: Optional[Union[List[int], List[str]]] = None,
    group_mapping: Optional[Dict[str, List[str]]] = None,
    validate: bool = False,
) -> LipidData:
    """
    Import lipid data from a CSV file with flexible column specification.

    The CSV file should have lipid names in one column (default: first column)
    and quantitation values in other columns (one per sample).

    The function automatically:
    - Validates and parses the CSV structure
    - Calls RefMet API to standardize lipid names and retrieve metadata
    - Creates sample metadata with automatic or user-specified group assignments
    - Returns a LipidData object with full annotation

    Args:
        filename: Path to CSV file
        lipid_col: Column index (0-based) or column name for lipid names (default: first column)
        sample_cols: List of column indices or names for sample data (default: all columns after lipid_col)
        group_mapping: Dict mapping group names to lists of sample IDs.
            Example: {"Control": ["S1", "S2"], "Treatment": ["S3", "S4"]}
            If not provided, groups are auto-extracted from sample IDs.
        validate: Whether to run data quality validation (default: False)

    Returns:
        LipidData object containing the imported and annotated data

    Examples:
        >>> # Basic import (auto-detect columns)
        >>> data = import_data("lipids.csv")

        >>> # Specify columns by index
        >>> data = import_data("lipids.csv", lipid_col=0, sample_cols=[1, 2, 3])

        >>> # Specify columns by name
        >>> data = import_data("lipids.csv", lipid_col="Name", sample_cols=["Control1", "Control2"])

        >>> # Specify group mapping
        >>> data = import_data(
        ...     "lipids.csv",
        ...     group_mapping={
        ...         "Control": ["Sample1", "Sample2"],
        ...         "Treatment": ["Sample3", "Sample4"]
        ...     }
        ... )

        >>> # With validation
        >>> data = import_data("lipids.csv", validate=True)
        >>> if data.manager.validation_report and not data.manager.validation_report.passed:
        ...     data.manager.validation_report.print_report()
    """
    logger.info(f"Importing lipid data from {filename}")

    # Create DataManager with configuration
    manager = DataManager(
        lipid_name_column=lipid_col,
        sample_columns=sample_cols,
        group_mapping=group_mapping,
        validate_data=validate,
    )
    dataset = manager.process_csv(filename)

    # Wrap in LipidData for high-level API
    lipid_data = LipidData(dataset=dataset, manager=manager)

    # Mark imported LipidData instances as not providing built-in reaction helpers
    # until explicit reaction-fetching/integration is implemented.
    try:
        object.__setattr__(lipid_data, "_reactions_available", False)
    except Exception:
        lipid_data._reactions_available = False

    logger.info(
        f"Import complete: {lipid_data.successful_import_count()} lipids, "
        f"{len(lipid_data.samples())} samples"
    )

    return lipid_data


def import_msdial(filename: Union[str, Path]) -> LipidData:
    """
    Import MS-DIAL formatted lipid data.

    MS-DIAL is a popular lipidomics software that exports data in a specific format.
    This function handles the MS-DIAL output format and converts it to LipidData.

    Args:
        filename: Path to MS-DIAL output file

    Returns:
        LipidData object

    Note:
        Current implementation treats MS-DIAL files as standard CSV.
        Future versions may add MS-DIAL-specific parsing logic.
    """
    logger.info(f"Importing MS-DIAL data from {filename}")

    # For now, treat as standard CSV
    # TODO: Add MS-DIAL-specific parsing logic if needed
    return import_data(filename)
