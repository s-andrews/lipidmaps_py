"""
Data module for LIPID MAPS Python API.

This module provides data models, analysis tools, and utilities for
working with lipidomics datasets.

Main components:
    - DataManager: Main class for processing CSV files
    - LipidDataset: Container for samples and lipids
    - QuantifiedLipid: Model for a lipid with quantitation values
    - SampleMetadata: Model for sample information
    - QuantitationAnalyzer: Analysis tools for quantitation
"""

from .data_manager import DataManager
from .models.sample import (
    LipidDataset,
    QuantifiedLipid,
    SampleMetadata,
    Quantitation,
)
from .models.reaction import (
    ReactionData,
    CompoundComponent,
    ReactionChecker,
    ReactionResponse,
)
from .quantitation import (
    QuantitationAnalyzer,
    QuantitationConfig,
    NormalizationMethod,
    QuantitationUnit,
    create_analyzer,
)

__all__ = [
    "DataManager",
    "LipidDataset",
    "QuantifiedLipid",
    "SampleMetadata",
    "Quantitation",
    "ReactionData",
    "CompoundComponent",
    "ReactionChecker",
    "ReactionResponse",
    "QuantitationAnalyzer",
    "QuantitationConfig",
    "NormalizationMethod",
    "QuantitationUnit",
    "create_analyzer",
]
