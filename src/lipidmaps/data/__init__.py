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
from .biopan_exporter import BioPANExporter
from .biopan_pathway_exporter import BioPANPathwayExporter
from .models.sample import (
    LipidAnnotation,
    LipidDataset,
    QuantifiedLipid,
    SampleConditions,
    SampleMetadata,
    Quantitation,
)
from .utils.chain_parser import (
    LipidStructure,
    parse_lipid,
)
from .models.refmet import (
    RefMet
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

# Pathway analysis
from .matching import (
    MatcherContext,
    MatcherRegistry,
    match_pathway_reactions,
    create_matcher_context,
)
from .models.species_reaction import (
    ClassReaction,
    CompoundRequirement,
    PathwayReactionSet,
    ReactionMatchResult,
    ReactionType,
    SpeciesReactionPair,
)
from .utils.chain_parser import (
    AcylChain,
    ChainParser,
    LipidStructure,
    SphingoBackbone,
    # FA/FACoA utilities
    get_common_fa_names,
    get_common_facoa_names,
    get_extended_fa_names,
    get_extended_facoa_names,
    infer_fa_from_lipids,
    infer_facoa_from_lipids,
)

__all__ = [
    "DataManager",
    "BioPANExporter",
    "BioPANPathwayExporter",
    "LipidDataset",
    "QuantifiedLipid",
    "SampleConditions",
    "SampleMetadata",
    "Quantitation",
    "RefMet",
    "ReactionData",
    "CompoundComponent",
    "ReactionChecker",
    "ReactionResponse",
    "QuantitationAnalyzer",
    "QuantitationConfig",
    "NormalizationMethod",
    "QuantitationUnit",
    "create_analyzer",
    # Matching
    "MatcherContext",
    "MatcherRegistry",
    "match_pathway_reactions",
    "create_matcher_context",
    "ClassReaction",
    "CompoundRequirement",
    "PathwayReactionSet",
    "ReactionMatchResult",
    "ReactionType",
    "SpeciesReactionPair",
    "AcylChain",
    "ChainParser",
    "LipidStructure",
    "SphingoBackbone",
    # FA/FACoA utilities
    "get_common_fa_names",
    "get_common_facoa_names",
    "get_extended_fa_names",
    "get_extended_facoa_names",
    "infer_fa_from_lipids",
    "infer_facoa_from_lipids",
]

