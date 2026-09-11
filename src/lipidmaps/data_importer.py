"""
Data import functions for LIPID MAPS

Provides high-level API for importing and processing lipid quantification data.
Built on top of DataManager for robust CSV processing and RefMet integration.
"""

import logging
from typing import Union, Optional, List, Dict, Any
from pathlib import Path

from pydantic import BaseModel, Field, PrivateAttr, computed_field, ConfigDict

from .data.data_manager import DataManager
from .data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata, PixelCoordinate
from .data.ingestion.imzml_reader import ImzMLIngestion, IonAnnotation, pixel_name
from .data.quantitation import (
    QuantitationAnalyzer,
    QuantitationConfig,
    NormalizationMethod,
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
    _analyzer: Optional[QuantitationAnalyzer] = PrivateAttr(default=None)

    def model_post_init(self, __context: Any) -> None:
        """Initialize manager if not provided."""
        if self.manager is None:
            object.__setattr__(self, "manager", DataManager(dataset=self.dataset))

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
    use_refmet: bool = True,
    use_headgroups: bool = True,
    fetch_reactions: bool = True,
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
        use_refmet=use_refmet,
        use_headgroups=use_headgroups,
        fetch_reactions=fetch_reactions,
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
    )

    return lipid_data


def import_imzml(
    imzml_path: Union[str, Path],
    annotation_path: Optional[Union[str, Path, List[IonAnnotation]]] = None,
    database: Optional[Union[str, Path, Any]] = None,
    mz_tolerance_ppm: float = 5.0,
    bbox: Optional[tuple] = None,
    adducts: Optional[list] = None,
    top_n_peaks: Optional[int] = None,
    stride: int = 1,
    progress=None,
    group: str = "tissue",
    use_refmet: Optional[bool] = None,
    use_headgroups: bool = True,
    fetch_reactions: bool = True,
) -> LipidData:
    """Import a mass-spectrometry-imaging (MSI) ``.imzML`` dataset.

    Each pixel becomes a sample (``SampleMetadata`` with ``coordinates``); each
    annotated ion becomes a ``QuantifiedLipid`` whose ``values`` hold that ion's
    per-pixel intensity. Names come from ``annotation_path`` (a companion table of
    lipid name + target m/z, e.g. a METASPACE export); the resulting names are then
    standardized and reaction-annotated through the SAME pipeline as
    :func:`import_data` (RefMet -> headgroups -> LMSD -> reactions).

    imzML holds only m/z + intensity, so lipid names must be supplied via
    ``annotation_path``; the m/z itself is never looked up against a database here.

    Args:
        imzml_path: Path to the ``.imzML`` file (its ``.ibd`` must sit alongside).
        annotation_path: Path to an annotation CSV (``name``, ``mz`` columns, plus
            optional ``adduct``/``formula``/``lm_id``) or a list of ``IonAnnotation``.
        mz_tolerance_ppm: Peak-match window when locating an ion in each pixel.
        bbox: Optional ``(x_min, y_min, x_max, y_max)`` inclusive crop to subsample.
        group: Group label assigned to every pixel sample (region grouping is a
            future extension).
        use_refmet, use_headgroups, fetch_reactions: Forwarded to the shared
            standardization/reaction pipeline.

    Returns:
        LipidData wrapping a spatial ``LipidDataset`` (``dataset.is_spatial`` True).

    Examples:
        >>> data = import_imzml("brain.imzML", "brain_annotations.csv")
        >>> data.dataset.is_spatial
        True
    """
    logger.info(f"Importing imzML data from {imzml_path}")

    ingestion = ImzMLIngestion()
    result = ingestion.read(
        imzml_path,
        annotations=annotation_path,
        mz_tolerance_ppm=mz_tolerance_ppm,
        bbox=bbox,
        database=database,
        adducts=adducts,
        top_n_peaks=top_n_peaks,
        stride=stride,
        progress=progress,
    )

    # Pixel -> SampleMetadata (carrying spatial coordinates).
    samples = [
        SampleMetadata(
            sample_name=name,
            group=group,
            coordinates=PixelCoordinate(x=x, y=y, z=z),
        )
        for (name, x, y, z) in result.pixels
    ]

    # Annotated ion -> QuantifiedLipid (per-pixel intensities keyed by pixel name).
    # Auto-annotated ions carry candidate molecule names for lm_id resolution.
    lipids = [
        QuantifiedLipid(
            input_name=ann.name,
            values=result.ion_values.get(ann.name, {}),
            annotation_candidates=(list(ann.candidates) if ann.candidates else None),
        )
        for ann in result.annotations
    ]

    dataset = LipidDataset(samples=samples, lipids=lipids)
    if result.pixel_size is not None:
        dataset.pixel_size_um = result.pixel_size

    # RefMet on formula-labeled auto-annotated ions is pointless (no names to match);
    # default it off in auto mode and rely on candidate-name resolution instead.
    if use_refmet is None:
        use_refmet = annotation_path is not None

    # Reuse the standardization + reaction annotation pipeline used by process_csv.
    manager = DataManager(
        use_refmet=use_refmet,
        use_headgroups=use_headgroups,
        fetch_reactions=fetch_reactions,
    )
    manager._annotate_and_react(dataset)

    lipid_data = LipidData(dataset=dataset, manager=manager)
    logger.info(
        f"imzML import complete: {len(lipids)} ions across {len(samples)} pixels"
    )
    return lipid_data


def import_imzml_stack(
    imzml_paths: List[Union[str, Path]],
    annotation_path: Optional[Union[str, Path, List[IonAnnotation]]] = None,
    database: Optional[Union[str, Path, Any]] = None,
    mz_tolerance_ppm: float = 5.0,
    bbox: Optional[tuple] = None,
    adducts: Optional[list] = None,
    top_n_peaks: Optional[int] = None,
    stride: int = 1,
    z_values: Optional[List[int]] = None,
    z_spacing_um: Optional[float] = None,
    progress=None,
    group: str = "tissue",
    use_refmet: Optional[bool] = None,
    use_headgroups: bool = True,
    fetch_reactions: bool = True,
) -> LipidData:
    """Stack several imzML files (one per serial tissue section) into a 3D dataset.

    Each file becomes one **z-layer** (``z_values[i]`` or the file order), overriding
    the per-file z (which is usually a meaningless constant of 1). Pixels merge into a
    single ``LipidDataset`` whose ``sample_name`` encodes (x, y, layer), so the whole
    volume shares one ion table, standardization run, and reaction set. Ions are unioned
    by their ``formula [adduct]`` (or annotation) label across files.

    Assumptions/limits: sections are treated as **already x/y-aligned** (no automatic
    image registration), and the same annotation source is applied to every file so ion
    labels line up. Provide ``z_spacing_um`` (section thickness) for a true physical z.

    Returns a ``LipidData`` whose ``dataset.is_3d`` is True (``z_slice_count`` == number
    of files).
    """
    paths = [Path(p).expanduser() for p in imzml_paths]
    if not paths:
        raise ValueError("import_imzml_stack requires at least one imzML path.")
    layers = z_values if z_values is not None else list(range(len(paths)))
    if len(layers) != len(paths):
        raise ValueError("z_values must have one entry per imzML path.")

    samples: List[SampleMetadata] = []
    merged_values: Dict[str, Dict[str, Optional[float]]] = {}
    candidates: Dict[str, list] = {}
    ordered_labels: List[str] = []
    pixel_size = None
    ingestion = ImzMLIngestion()

    for i, path in enumerate(paths):
        layer = int(layers[i])

        def _p(phase, done, total, _layer=layer, _i=i):
            if progress:
                progress(f"file {_i + 1}/{len(paths)}:{phase}", done, total)

        result = ingestion.read(
            path, annotations=annotation_path, mz_tolerance_ppm=mz_tolerance_ppm,
            bbox=bbox, database=database, adducts=adducts, top_n_peaks=top_n_peaks,
            stride=stride, progress=_p,
        )
        if pixel_size is None and result.pixel_size is not None:
            pixel_size = result.pixel_size
        for ann in result.annotations:
            if ann.name not in merged_values:
                merged_values[ann.name] = {}
                ordered_labels.append(ann.name)
                candidates[ann.name] = list(ann.candidates) if ann.candidates else None
        for (name, x, y, _z) in result.pixels:
            sn = pixel_name(x, y, layer)
            samples.append(SampleMetadata(
                sample_name=sn, group=group, coordinates=PixelCoordinate(x=x, y=y, z=layer),
            ))
            for label in result.ion_values:
                merged_values[label][sn] = result.ion_values[label].get(name)

    lipids = [
        QuantifiedLipid(
            input_name=label, values=merged_values[label],
            annotation_candidates=candidates.get(label),
        )
        for label in ordered_labels
    ]
    dataset = LipidDataset(samples=samples, lipids=lipids)
    if pixel_size is not None:
        dataset.pixel_size_um = pixel_size
    dataset.z_spacing_um = z_spacing_um

    if use_refmet is None:
        use_refmet = annotation_path is not None
    manager = DataManager(
        use_refmet=use_refmet, use_headgroups=use_headgroups, fetch_reactions=fetch_reactions,
    )
    manager._annotate_and_react(dataset)

    lipid_data = LipidData(dataset=dataset, manager=manager)
    logger.info(
        "imzML stack import complete: %d ions across %d pixels in %d z-slices",
        len(lipids), len(samples), dataset.z_slice_count,
    )
    return lipid_data


#TODO to implement in future
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
