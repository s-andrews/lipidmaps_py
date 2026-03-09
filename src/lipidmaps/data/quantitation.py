"""
Quantitation analysis module for lipid datasets.

Provides normalization methods, statistical analysis, and group-level quantitation
functions for lipidomics data.
"""

import logging
import math
from typing import Dict, List, Optional, Tuple, Union, Any, TYPE_CHECKING
from enum import Enum

import numpy as np
from scipy import stats

from .models.base import LipidmapsBaseModel
from pydantic import PrivateAttr

if TYPE_CHECKING:
    from .models.sample import LipidDataset, QuantifiedLipid, SampleMetadata
    from .models.reaction import ReactionData, CompoundComponent

logger = logging.getLogger(__name__)


class NormalizationMethod(str, Enum):
    """Available normalization methods."""
    NONE = "none"
    TOTAL_LIPID = "total_lipid"
    INTERNAL_STANDARD = "internal_standard"
    LOG2 = "log2"
    LOG10 = "log10"
    MEDIAN_CENTER = "median_center"
    ZSCORE = "zscore"
    QUANTILE = "quantile"


class QuantitationUnit(str, Enum):
    """Common quantitation units for lipidomics."""
    PMOL = "pmol"
    NMOL = "nmol"
    NG = "ng"
    UG = "ug"
    MG = "mg"
    AREA = "area"  # Raw peak area
    RELATIVE = "relative"  # Normalized/relative abundance
    PERCENT = "percent"  # Percentage
    ARBITRARY = "AU"  # Arbitrary units


class QuantitationConfig(LipidmapsBaseModel):
    """Configuration for quantitation tracking and normalization."""
    
    unit: Optional[str] = None
    method: Optional[str] = None  # e.g., 'LC-MS', 'GC-MS'
    normalization: NormalizationMethod = NormalizationMethod.NONE
    internal_standard_lipid: Optional[str] = None  # Lipid name used as internal standard
    notes: Optional[str] = None


class QuantitationAnalyzer(LipidmapsBaseModel):
    """
    Analyzer class for performing quantitation operations on lipid datasets.

    Provides methods for:
    - Normalization (total lipid, internal standard, log, median centering)
    - Statistical analysis (fold change, p-values, CV)
    - Group-level aggregations
    - Reaction component quantitation
    """

    dataset: Any
    model_config = {"arbitrary_types_allowed": True}
    _config: QuantitationConfig = PrivateAttr(default_factory=QuantitationConfig)

    def __init__(self, dataset: Any = None, **data):
        """Allow positional-style construction `QuantitationAnalyzer(dataset)` for tests.

        If `dataset` is provided positionally, inject it into kwargs forwarded
        to the Pydantic BaseModel initializer.
        """
        if dataset is not None:
            data.setdefault("dataset", dataset)
        super().__init__(**data)

    def __repr__(self) -> str:
        return f"QuantitationAnalyzer(dataset={getattr(self.dataset, 'name', None)})"

    @property
    def config(self) -> QuantitationConfig:
        """Get the current quantitation configuration."""
        return self._config

    def set_config(
        self,
        unit: Optional[str] = None,
        method: Optional[str] = None,
        normalization: Optional[NormalizationMethod] = None,
        internal_standard_lipid: Optional[str] = None,
        notes: Optional[str] = None,
    ) -> None:
        """Update quantitation configuration."""
        if unit is not None:
            self._config.unit = unit
        if method is not None:
            self._config.method = method
        if normalization is not None:
            self._config.normalization = normalization
        if internal_standard_lipid is not None:
            self._config.internal_standard_lipid = internal_standard_lipid
        if notes is not None:
            self._config.notes = notes

    # =========================================================================
    # NORMALIZATION METHODS
    # =========================================================================

    def normalize_total_lipid(self, scale_factor: float = 1e6) -> Dict[str, Dict[str, float]]:
        """
        Normalize values by total lipid content per sample.

        Each sample's values are divided by the sum of all lipid values in that sample,
        then optionally scaled (e.g., to parts per million).

        Args:
            scale_factor: Factor to multiply normalized values by (default: 1e6 for ppm)

        Returns:
            Dict mapping lipid input_name -> {sample_name: normalized_value}
        """
        normalized: Dict[str, Dict[str, float]] = {}

        # Calculate total per sample
        sample_totals: Dict[str, float] = {}
        for sample in self.dataset.samples:
            total = 0.0
            for lipid in self.dataset.lipids:
                val = lipid.values.get(sample.sample_name)
                if val is not None and not math.isnan(val):
                    total += val
            sample_totals[sample.sample_name] = total if total > 0 else 1.0

        # Normalize each lipid
        for lipid in self.dataset.lipids:
            normalized[lipid.input_name] = {}
            for sample_name, val in lipid.values.items():
                if val is not None and not math.isnan(val):
                    normalized[lipid.input_name][sample_name] = (
                        val / sample_totals.get(sample_name, 1.0)
                    ) * scale_factor
                else:
                    normalized[lipid.input_name][sample_name] = float('nan')

        logger.info(f"Total lipid normalization applied (scale_factor={scale_factor})")
        return normalized

    def normalize_internal_standard(
        self,
        internal_standard: Union[str, "QuantifiedLipid"],
    ) -> Dict[str, Dict[str, float]]:
        """
        Normalize values by an internal standard lipid.

        Each sample's values are divided by the internal standard's value in that sample.

        Args:
            internal_standard: Lipid name (str) or QuantifiedLipid object to use as standard

        Returns:
            Dict mapping lipid input_name -> {sample_name: normalized_value}

        Raises:
            ValueError: If internal standard is not found in dataset
        """
        # Find the internal standard lipid
        if isinstance(internal_standard, str):
            std_lipid = next(
                (l for l in self.dataset.lipids
                 if l.input_name.lower() == internal_standard.lower()
                 or (l.standardized_name and l.standardized_name.lower() == internal_standard.lower())),
                None
            )
            if std_lipid is None:
                raise ValueError(f"Internal standard '{internal_standard}' not found in dataset")
        else:
            std_lipid = internal_standard

        normalized: Dict[str, Dict[str, float]] = {}

        for lipid in self.dataset.lipids:
            normalized[lipid.input_name] = {}
            for sample_name, val in lipid.values.items():
                std_val = std_lipid.values.get(sample_name)
                if std_val is not None and std_val != 0 and not math.isnan(std_val):
                    if val is not None and not math.isnan(val):
                        normalized[lipid.input_name][sample_name] = val / std_val
                    else:
                        normalized[lipid.input_name][sample_name] = float('nan')
                else:
                    normalized[lipid.input_name][sample_name] = float('nan')

        self._config.internal_standard_lipid = std_lipid.input_name
        logger.info(f"Internal standard normalization applied (standard={std_lipid.input_name})")
        return normalized

    def normalize_log(
        self,
        base: int = 2,
        offset: float = 1.0,
    ) -> Dict[str, Dict[str, float]]:
        """
        Apply log transformation to values.

        Args:
            base: Logarithm base (2 or 10)
            offset: Value to add before log transform to handle zeros (default: 1.0)

        Returns:
            Dict mapping lipid input_name -> {sample_name: log_transformed_value}
        """
        if base not in (2, 10):
            raise ValueError("Log base must be 2 or 10")

        log_func = math.log2 if base == 2 else math.log10
        normalized: Dict[str, Dict[str, float]] = {}

        for lipid in self.dataset.lipids:
            normalized[lipid.input_name] = {}
            for sample_name, val in lipid.values.items():
                if val is not None and not math.isnan(val):
                    # Add offset to handle zeros
                    adjusted_val = val + offset
                    if adjusted_val <= 0:
                        # Can't take log of zero or negative
                        normalized[lipid.input_name][sample_name] = float('-inf')
                    else:
                        normalized[lipid.input_name][sample_name] = log_func(adjusted_val)
                else:
                    normalized[lipid.input_name][sample_name] = float('nan')

        logger.info(f"Log{base} transformation applied (offset={offset})")
        return normalized

    def normalize_median_center(self) -> Dict[str, Dict[str, float]]:
        """
        Center values by subtracting the median per sample.

        Each sample's values are centered by subtracting the median value across
        all lipids in that sample.

        Returns:
            Dict mapping lipid input_name -> {sample_name: centered_value}
        """
        normalized: Dict[str, Dict[str, float]] = {}

        # Calculate median per sample
        sample_medians: Dict[str, float] = {}
        for sample in self.dataset.samples:
            sample_vals = []
            for lipid in self.dataset.lipids:
                val = lipid.values.get(sample.sample_name)
                if val is not None and not math.isnan(val):
                    sample_vals.append(val)
            if sample_vals:
                sample_medians[sample.sample_name] = float(np.median(sample_vals))
            else:
                sample_medians[sample.sample_name] = 0.0

        # Center each lipid
        for lipid in self.dataset.lipids:
            normalized[lipid.input_name] = {}
            for sample_name, val in lipid.values.items():
                if val is not None and not math.isnan(val):
                    normalized[lipid.input_name][sample_name] = (
                        val - sample_medians.get(sample_name, 0.0)
                    )
                else:
                    normalized[lipid.input_name][sample_name] = float('nan')

        logger.info("Median centering applied")
        return normalized

    def apply_normalization(
        self,
        method: NormalizationMethod,
        internal_standard: Optional[str] = None,
        log_base: int = 2,
        log_offset: float = 1.0,
        total_lipid_scale: float = 1e6,
        in_place: bool = False,
    ) -> Dict[str, Dict[str, float]]:
        """
        Apply a normalization method to the dataset.

        Args:
            method: Normalization method to apply
            internal_standard: Required if method is INTERNAL_STANDARD
            log_base: Base for log transformation (2 or 10)
            log_offset: Offset for log transformation
            total_lipid_scale: Scale factor for total lipid normalization
            in_place: If True, update lipid values in place

        Returns:
            Dict mapping lipid input_name -> {sample_name: normalized_value}
        """
        if method == NormalizationMethod.NONE:
            return {l.input_name: l.values.copy() for l in self.dataset.lipids}
        elif method == NormalizationMethod.TOTAL_LIPID:
            result = self.normalize_total_lipid(scale_factor=total_lipid_scale)
        elif method == NormalizationMethod.INTERNAL_STANDARD:
            if internal_standard is None:
                raise ValueError("internal_standard required for INTERNAL_STANDARD normalization")
            result = self.normalize_internal_standard(internal_standard)
        elif method == NormalizationMethod.LOG2:
            result = self.normalize_log(base=2, offset=log_offset)
        elif method == NormalizationMethod.LOG10:
            result = self.normalize_log(base=10, offset=log_offset)
        elif method == NormalizationMethod.MEDIAN_CENTER:
            result = self.normalize_median_center()
        elif method == NormalizationMethod.ZSCORE:
            result = {}
            for lipid in self.dataset.lipids:
                result[lipid.input_name] = lipid.zscore()
        else:
            raise ValueError(f"Unknown normalization method: {method}")

        self._config.normalization = method

        if in_place:
            for lipid in self.dataset.lipids:
                if lipid.input_name in result:
                    lipid.values = result[lipid.input_name]
            logger.info(f"Applied {method.value} normalization in place")

        return result

    # =========================================================================
    # STATISTICAL ANALYSIS
    # =========================================================================

    def calculate_fold_change(
        self,
        group1: str,
        group2: str,
        log2: bool = True,
    ) -> Dict[str, float]:
        """
        Calculate fold change between two groups for each lipid.

        Args:
            group1: Name of first group (numerator)
            group2: Name of second group (denominator/reference)
            log2: If True, return log2 fold change

        Returns:
            Dict mapping lipid input_name -> fold_change
        """
        # Get sample names for each group
        group1_samples = [s.sample_name for s in self.dataset.samples if s.group == group1]
        group2_samples = [s.sample_name for s in self.dataset.samples if s.group == group2]

        if not group1_samples or not group2_samples:
            logger.warning(f"Empty groups: {group1}={len(group1_samples)}, {group2}={len(group2_samples)}")
            return {}

        fold_changes: Dict[str, float] = {}

        for lipid in self.dataset.lipids:
            # Get means for each group
            g1_vals = [lipid.values.get(s) for s in group1_samples if lipid.values.get(s) is not None]
            g2_vals = [lipid.values.get(s) for s in group2_samples if lipid.values.get(s) is not None]

            if g1_vals and g2_vals:
                mean1 = np.nanmean(g1_vals)
                mean2 = np.nanmean(g2_vals)

                if mean2 != 0 and not math.isnan(mean2):
                    fc = mean1 / mean2
                    if log2:
                        fc = math.log2(fc) if fc > 0 else float('nan')
                    fold_changes[lipid.input_name] = fc
                else:
                    fold_changes[lipid.input_name] = float('nan')
            else:
                fold_changes[lipid.input_name] = float('nan')

        logger.info(f"Calculated fold change: {group1} vs {group2} (log2={log2})")
        return fold_changes

    def calculate_pvalue(
        self,
        group1: str,
        group2: str,
        test: str = "ttest",
        paired: bool = False,
    ) -> Dict[str, float]:
        """
        Calculate p-values comparing two groups for each lipid.

        Args:
            group1: Name of first group
            group2: Name of second group
            test: Statistical test to use ('ttest' or 'mannwhitney')
            paired: If True, use paired test (only for ttest)

        Returns:
            Dict mapping lipid input_name -> p_value
        """
        group1_samples = [s.sample_name for s in self.dataset.samples if s.group == group1]
        group2_samples = [s.sample_name for s in self.dataset.samples if s.group == group2]

        if not group1_samples or not group2_samples:
            return {}

        pvalues: Dict[str, float] = {}

        for lipid in self.dataset.lipids:
            g1_vals = [lipid.values.get(s) for s in group1_samples
                       if lipid.values.get(s) is not None and not math.isnan(lipid.values.get(s))]
            g2_vals = [lipid.values.get(s) for s in group2_samples
                       if lipid.values.get(s) is not None and not math.isnan(lipid.values.get(s))]

            if len(g1_vals) >= 2 and len(g2_vals) >= 2:
                try:
                    if test == "ttest":
                        if paired and len(g1_vals) == len(g2_vals):
                            _, pval = stats.ttest_rel(g1_vals, g2_vals)
                        else:
                            _, pval = stats.ttest_ind(g1_vals, g2_vals)
                    elif test == "mannwhitney":
                        _, pval = stats.mannwhitneyu(g1_vals, g2_vals, alternative='two-sided')
                    else:
                        pval = float('nan')
                    pvalues[lipid.input_name] = pval
                except Exception:
                    pvalues[lipid.input_name] = float('nan')
            else:
                pvalues[lipid.input_name] = float('nan')

        logger.info(f"Calculated p-values: {group1} vs {group2} (test={test}, paired={paired})")
        return pvalues

    def calculate_cv(self, by_group: bool = False) -> Union[Dict[str, float], Dict[str, Dict[str, float]]]:
        """
        Calculate coefficient of variation (CV) for each lipid.

        CV = (standard deviation / mean) * 100

        Args:
            by_group: If True, calculate CV separately for each group

        Returns:
            If by_group=False: Dict mapping lipid input_name -> CV (%)
            If by_group=True: Dict mapping group -> {lipid input_name -> CV (%)}
        """
        if not by_group:
            cvs: Dict[str, float] = {}
            for lipid in self.dataset.lipids:
                vals = [v for v in lipid.values.values() if v is not None and not math.isnan(v)]
                if len(vals) >= 2:
                    mean = np.mean(vals)
                    std = np.std(vals, ddof=1)  # Sample std
                    if mean != 0:
                        cvs[lipid.input_name] = (std / abs(mean)) * 100
                    else:
                        cvs[lipid.input_name] = float('nan')
                else:
                    cvs[lipid.input_name] = float('nan')
            return cvs

        # By group
        group_cvs: Dict[str, Dict[str, float]] = {}
        groups = set(s.group for s in self.dataset.samples)

        for group in groups:
            group_samples = [s.sample_name for s in self.dataset.samples if s.group == group]
            group_cvs[group] = {}

            for lipid in self.dataset.lipids:
                vals = [lipid.values.get(s) for s in group_samples
                        if lipid.values.get(s) is not None and not math.isnan(lipid.values.get(s))]
                if len(vals) >= 2:
                    mean = np.mean(vals)
                    std = np.std(vals, ddof=1)
                    if mean != 0:
                        group_cvs[group][lipid.input_name] = (std / abs(mean)) * 100
                    else:
                        group_cvs[group][lipid.input_name] = float('nan')
                else:
                    group_cvs[group][lipid.input_name] = float('nan')

        return group_cvs

    def differential_analysis(
        self,
        group1: str,
        group2: str,
        fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05,
        test: str = "ttest",
    ) -> List[Dict[str, Any]]:
        """
        Perform differential analysis between two groups.

        Args:
            group1: Name of first group
            group2: Name of second group (reference)
            fc_threshold: Absolute log2 fold change threshold
            pvalue_threshold: P-value significance threshold
            test: Statistical test to use

        Returns:
            List of dicts with lipid analysis results, sorted by significance
        """
        fold_changes = self.calculate_fold_change(group1, group2, log2=True)
        pvalues = self.calculate_pvalue(group1, group2, test=test)

        results = []
        for lipid in self.dataset.lipids:
            name = lipid.input_name
            fc = fold_changes.get(name, float('nan'))
            pval = pvalues.get(name, float('nan'))

            is_significant = (
                not math.isnan(fc) and
                not math.isnan(pval) and
                abs(fc) >= fc_threshold and
                pval <= pvalue_threshold
            )

            direction = "up" if fc > 0 else "down" if fc < 0 else "unchanged"

            results.append({
                "lipid_name": name,
                "standardized_name": lipid.standardized_name,
                "lm_id": lipid.lm_id,
                "log2_fold_change": fc,
                "p_value": pval,
                "significant": is_significant,
                "direction": direction,
            })

        # Sort by p-value
        results.sort(key=lambda x: x["p_value"] if not math.isnan(x["p_value"]) else 1.0)

        significant_count = sum(1 for r in results if r["significant"])
        logger.info(f"Differential analysis: {significant_count} significant lipids")

        return results

    # =========================================================================
    # GROUP-LEVEL QUANTITATION
    # =========================================================================

    def get_group_means(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate mean values per group for each lipid.

        Returns:
            Dict mapping group_name -> {lipid_name: mean_value}
        """
        groups: Dict[str, List[str]] = {}
        for sample in self.dataset.samples:
            groups.setdefault(sample.group, []).append(sample.sample_name)

        result: Dict[str, Dict[str, float]] = {}
        for group, sample_names in groups.items():
            result[group] = {}
            for lipid in self.dataset.lipids:
                vals = [lipid.values.get(s) for s in sample_names
                        if lipid.values.get(s) is not None and not math.isnan(lipid.values.get(s))]
                if vals:
                    result[group][lipid.input_name] = float(np.mean(vals))
                else:
                    result[group][lipid.input_name] = float('nan')

        return result

    def get_group_stds(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate standard deviation per group for each lipid.

        Returns:
            Dict mapping group_name -> {lipid_name: std_value}
        """
        groups: Dict[str, List[str]] = {}
        for sample in self.dataset.samples:
            groups.setdefault(sample.group, []).append(sample.sample_name)

        result: Dict[str, Dict[str, float]] = {}
        for group, sample_names in groups.items():
            result[group] = {}
            for lipid in self.dataset.lipids:
                vals = [lipid.values.get(s) for s in sample_names
                        if lipid.values.get(s) is not None and not math.isnan(lipid.values.get(s))]
                if len(vals) >= 2:
                    result[group][lipid.input_name] = float(np.std(vals, ddof=1))
                else:
                    result[group][lipid.input_name] = float('nan')

        return result

    def get_group_summary(self) -> Dict[str, Dict[str, Dict[str, float]]]:
        """
        Get comprehensive summary statistics per group.

        Returns:
            Dict mapping group_name -> {
                lipid_name: {mean, std, min, max, median, n}
            }
        """
        groups: Dict[str, List[str]] = {}
        for sample in self.dataset.samples:
            groups.setdefault(sample.group, []).append(sample.sample_name)

        result: Dict[str, Dict[str, Dict[str, float]]] = {}
        for group, sample_names in groups.items():
            result[group] = {}
            for lipid in self.dataset.lipids:
                vals = [lipid.values.get(s) for s in sample_names
                        if lipid.values.get(s) is not None and not math.isnan(lipid.values.get(s))]
                if vals:
                    result[group][lipid.input_name] = {
                        "mean": float(np.mean(vals)),
                        "std": float(np.std(vals, ddof=1)) if len(vals) >= 2 else 0.0,
                        "min": float(np.min(vals)),
                        "max": float(np.max(vals)),
                        "median": float(np.median(vals)),
                        "n": len(vals),
                    }
                else:
                    result[group][lipid.input_name] = {
                        "mean": float('nan'),
                        "std": float('nan'),
                        "min": float('nan'),
                        "max": float('nan'),
                        "median": float('nan'),
                        "n": 0,
                    }

        return result

    def compare_groups(
        self,
        group1: str,
        group2: str,
    ) -> Dict[str, Dict[str, Any]]:
        """
        Compare two groups with full statistics.

        Args:
            group1: First group name
            group2: Second group name

        Returns:
            Dict mapping lipid_name -> {
                group1_mean, group1_std, group2_mean, group2_std,
                log2_fc, p_value, significant
            }
        """
        means = self.get_group_means()
        stds = self.get_group_stds()
        fold_changes = self.calculate_fold_change(group1, group2, log2=True)
        pvalues = self.calculate_pvalue(group1, group2)

        result: Dict[str, Dict[str, Any]] = {}
        for lipid in self.dataset.lipids:
            name = lipid.input_name
            fc = fold_changes.get(name, float('nan'))
            pval = pvalues.get(name, float('nan'))

            result[name] = {
                f"{group1}_mean": means.get(group1, {}).get(name, float('nan')),
                f"{group1}_std": stds.get(group1, {}).get(name, float('nan')),
                f"{group2}_mean": means.get(group2, {}).get(name, float('nan')),
                f"{group2}_std": stds.get(group2, {}).get(name, float('nan')),
                "log2_fold_change": fc,
                "p_value": pval,
                "significant": (not math.isnan(pval) and pval < 0.05 and
                               not math.isnan(fc) and abs(fc) >= 1.0),
            }

        return result

    # =========================================================================
    # REACTION COMPONENT QUANTITATION
    # =========================================================================

    def get_value_for_reaction_component(
        self,
        component: Union[str, "CompoundComponent"],
        sample: Union[str, "SampleMetadata"],
        method: str = "sum",
    ) -> Optional[float]:
        """
        Get quantitation value for a reaction component.

        A reaction component may map to multiple lipids in the dataset.
        This method aggregates values across all matching lipids.

        Args:
            component: Component LM ID (str) or CompoundComponent object
            sample: Sample name (str) or SampleMetadata object
            method: Aggregation method ('sum', 'mean', 'max', 'min')

        Returns:
            Aggregated quantitation value, or None if no matches
        """
        # Resolve component ID
        if hasattr(component, 'compound_lm_id'):
            comp_id = component.compound_lm_id
        else:
            comp_id = str(component)

        if not comp_id:
            return None

        # Resolve sample name
        if hasattr(sample, 'sample_name'):
            sample_name = sample.sample_name
        else:
            sample_name = str(sample)

        # Find matching lipids
        comp_id_lower = comp_id.lower()
        matching_lipids = [
            l for l in self.dataset.lipids
            if (l.lm_id and l.lm_id.lower() == comp_id_lower) or
               (l.generic_lm_id and l.generic_lm_id.lower() == comp_id_lower)
        ]

        if not matching_lipids:
            return None

        # Collect values
        values = []
        for lipid in matching_lipids:
            val = lipid.values.get(sample_name)
            if val is not None and not math.isnan(val):
                values.append(val)

        if not values:
            return None

        # Aggregate
        if method == "sum":
            return sum(values)
        elif method == "mean":
            return float(np.mean(values))
        elif method == "max":
            return max(values)
        elif method == "min":
            return min(values)
        else:
            raise ValueError(f"Unknown aggregation method: {method}")

    def get_values_for_reaction_component(
        self,
        component: Union[str, "CompoundComponent"],
        method: str = "sum",
    ) -> Dict[str, float]:
        """
        Get quantitation values for a reaction component across all samples.

        Args:
            component: Component LM ID (str) or CompoundComponent object
            method: Aggregation method ('sum', 'mean', 'max', 'min')

        Returns:
            Dict mapping sample_name -> aggregated_value
        """
        result: Dict[str, float] = {}
        for sample in self.dataset.samples:
            val = self.get_value_for_reaction_component(component, sample, method)
            if val is not None:
                result[sample.sample_name] = val
        return result

    def get_reaction_flux_estimate(
        self,
        reaction: "ReactionData",
        sample: Union[str, "SampleMetadata"],
        method: str = "ratio",
    ) -> Optional[float]:
        """
        Estimate relative flux for a reaction based on reactant/product ratios.

        This is a simplified estimate based on abundance ratios.

        Args:
            reaction: ReactionData object
            sample: Sample name or SampleMetadata object
            method: 'ratio' (products/reactants) or 'difference' (products-reactants)

        Returns:
            Estimated flux value, or None if insufficient data
        """
        reactant_sum = 0.0
        product_sum = 0.0

        for reactant in reaction.reactants:
            val = self.get_value_for_reaction_component(reactant, sample, method="sum")
            if val is not None:
                reactant_sum += val

        for product in reaction.products:
            val = self.get_value_for_reaction_component(product, sample, method="sum")
            if val is not None:
                product_sum += val

        if method == "ratio":
            if reactant_sum > 0:
                return product_sum / reactant_sum
            return None
        elif method == "difference":
            if reactant_sum > 0 or product_sum > 0:
                return product_sum - reactant_sum
            return None
        else:
            raise ValueError(f"Unknown method: {method}")


# Ensure forward-referenced models are resolved when possible (pydantic v2)
try:
    # LipidDataset and related models live in models.sample; import to ensure
    # QuantitationAnalyzer's forward references are resolved.
    from .models.sample import LipidDataset  # noqa: F401

    # Rebuild pydantic models so forward refs are resolved (safe no-op if already built)
    try:
        QuantitationAnalyzer.model_rebuild()
    except Exception:
        # If model_rebuild is not available or fails, ignore to avoid breaking imports
        pass
except Exception:
    # If import fails due to import-order issues, skip — the caller can rebuild later.
    pass


def create_analyzer(dataset: "LipidDataset") -> QuantitationAnalyzer:
    """
    Factory function to create a QuantitationAnalyzer for a dataset.

    Args:
        dataset: LipidDataset to analyze

    Returns:
        QuantitationAnalyzer instance
    """
    return QuantitationAnalyzer(dataset)
