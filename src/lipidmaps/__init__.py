"""
LIPID MAPS Python API
=====================

A Python package for importing, processing, and analyzing lipid data
using LIPID MAPS resources.

Main functions:
    - import_data: Import lipid data from CSV files
    - process_csv: Convenience function to process CSV directly

Main classes:
    - QuantifiedLipid: Typed model representing a quantified lipid, its per-sample values, and annotation metadata.
    - LipidDataset: Container for samples and `QuantifiedLipid` objects with helpers for querying, aggregation, normalization, and reaction annotation.
    - RefMet: Utilities for RefMet name standardization and LM ID lookup for input lipid names.
    - Reactions: Tools and models for fetching, filtering, and inspecting LIPID MAPS reactions.
    - QuantitationAnalyzer: Analysis tools for normalization and quantitation workflows.


"""

# Convenience wrapper for direct CSV processing
from .data.data_manager import DataManager
from .data.biopan_exporter import BioPANExporter
from .data.biopan_pathway_exporter import BioPANPathwayExporter


def process_csv(csv_path, **kwargs):
    """
    Directly process a CSV file into a LipidDataset using DataManager.
    Accepts all DataManager keyword arguments.
    Example:
        from lipidmaps import process_csv
        dataset = process_csv("your_file.csv", validate_data=True)
    """
    manager = DataManager(**kwargs)
    return manager.process_csv(csv_path)

from .data_importer import import_data, import_msdial, LipidData

# Quantitation analysis exports
from .data.quantitation import (
    QuantitationAnalyzer,
    QuantitationConfig,
    NormalizationMethod,
    QuantitationUnit,
    create_analyzer,
)

# Import subpackages to make them accessible
from . import data

__version__ = "0.1.0"
__all__ = [
    "import_data",
    "process_csv",
    "DataManager",
    "BioPANExporter",
    "BioPANPathwayExporter",
    "QuantitationAnalyzer",
    "QuantitationConfig",
    "NormalizationMethod",
    "QuantitationUnit",
    "create_analyzer",
    "data",
]
