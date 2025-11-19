"""
Data import functions for LIPID MAPS

Provides high-level API for importing and processing lipid quantification data.
Built on top of DataManager for robust CSV processing and RefMet integration.
"""
import logging
from typing import Union, Optional, List, Dict, Any
from pathlib import Path

from pydantic import BaseModel, Field, computed_field, ConfigDict

from .biopan.data_manager import DataManager
from .biopan.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata

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
    
    dataset: LipidDataset = Field(default_factory=lambda: LipidDataset(samples=[], lipids=[]))
    manager: Optional[DataManager] = Field(default=None)
    
    def model_post_init(self, __context: Any) -> None:
        """Initialize manager if not provided."""
        if self.manager is None:
            object.__setattr__(self, 'manager', DataManager(dataset=self.dataset))

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
        return [s.sample_id for s in self.dataset.samples]

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

    def get_value_for_lipid(self, lipid: Union[str, QuantifiedLipid], sample: str) -> Optional[float]:
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
                (l for l in self.dataset.lipids if l.input_name == lipid or l.standardized_name == lipid),
                None
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
            (l for l in self.dataset.lipids if l.input_name == name or l.standardized_name == name),
            None
        )

    def get_lipids_by_class(self, lipid_class: str) -> List[QuantifiedLipid]:
        """Get all lipids belonging to a specific class.
        
        Args:
            lipid_class: Lipid class name (e.g., 'PC', 'TAG')
            
        Returns:
            List of QuantifiedLipid objects
        """
        return [
            lipid for lipid in self.dataset.lipids
            if lipid.main_class == lipid_class or lipid.sub_class == lipid_class
        ]

    def get_lm_ids(self) -> List[str]:
        """Get all unique LIPID MAPS IDs from the dataset.
        
        Returns:
            List of LM IDs
        """
        return list(set(
            lipid.lm_id for lipid in self.dataset.lipids
            if lipid.lm_id and lipid.lm_id.startswith('LM')
        ))

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

    def get_lipids_for_reaction_component(self, component: str):
        """Get valid lipids for a reaction component.
        
        Note:
            Future implementation - requires reactions integration
        """
        raise NotImplementedError("Reactions integration not yet implemented.")

    def get_value_for_reaction_component(self, component: str, method: str = "sum"):
        """Get quantitation for a reaction component.
        
        Note:
            Future implementation - requires reactions integration
        """
        raise NotImplementedError("Reactions integration not yet implemented.")


def import_data(
    filename: Union[str, Path],
    lipid_col: Optional[Union[int, str]] = None,
    sample_cols: Optional[Union[List[int], List[str]]] = None,
    group_mapping: Optional[Dict[str, List[str]]] = None,
    validate: bool = False
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
        validate_data=validate
    )
    dataset = manager.process_csv(filename)
    
    # Wrap in LipidData for high-level API
    lipid_data = LipidData(dataset=dataset, manager=manager)
    
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
