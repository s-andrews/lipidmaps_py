"""
Base class for reaction matchers (Strategy pattern).

Each matcher implements the logic for matching reactant-product pairs
for a specific reaction type (same structure, FA release, etc.).
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Set, Tuple

from pydantic import BaseModel, Field

from ..models.species_reaction import (
    ClassReaction,
    CompoundRequirement,
    ReactionMatchResult,
    ReactionType,
    SpeciesReactionPair,
)
from ..utils.chain_parser import AcylChain, ChainParser, LipidStructure


class MatcherContext(BaseModel):
    """Context data passed to matchers for species matching.
    
    Contains the available lipid species and compounds from the dataset.
    """
    # Species grouped by headgroup
    species_by_class: Dict[str, List[LipidStructure]] = Field(default_factory=dict)
    
    # Available FA and FACoA compounds
    fa_species: Dict[str, AcylChain] = Field(default_factory=dict,
        description="Available FA keyed by chain string (e.g., '16:0')")
    facoa_species: Dict[str, AcylChain] = Field(default_factory=dict,
        description="Available FACoA keyed by chain string (e.g., '18:1')")
    
    # All lipid names for quick lookup
    all_lipid_names: Set[str] = Field(default_factory=set)
    
    class Config:
        arbitrary_types_allowed = True
    
    def get_species(self, lipid_class: str) -> List[LipidStructure]:
        """Get all species for a given lipid class."""
        # Handle class name variations
        normalized = lipid_class.upper()
        for key in self.species_by_class:
            if key.upper() == normalized:
                return self.species_by_class[key]
        return []
    
    def has_fa(self, chain: AcylChain) -> bool:
        """Check if FA with this chain is available."""
        key = chain.notation
        return key in self.fa_species
    
    def has_facoa(self, chain: AcylChain) -> bool:
        """Check if FACoA with this chain is available."""
        key = chain.notation
        return key in self.facoa_species
    
    def fa_for_difference(self, carbon_diff: int, bond_diff: int) -> Optional[AcylChain]:
        """Find FA matching a carbon:bond difference."""
        key = f"{carbon_diff}:{bond_diff}"
        return self.fa_species.get(key)
    
    def facoa_for_difference(self, carbon_diff: int, bond_diff: int) -> Optional[AcylChain]:
        """Find FACoA matching a carbon:bond difference."""
        key = f"{carbon_diff}:{bond_diff}"
        return self.facoa_species.get(key)


class ReactionMatcher(ABC):
    """Abstract base class for reaction matching strategies.
    
    Each concrete implementation handles one reaction type.
    Implements the Strategy pattern for extensible reaction matching.
    """
    
    @property
    @abstractmethod
    def reaction_type(self) -> ReactionType:
        """The reaction type this matcher handles."""
        ...
    
    @abstractmethod
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match species-level pairs for a class reaction.
        
        Args:
            class_reaction: The class-level reaction to match
            context: Dataset context with available species and compounds
            
        Returns:
            ReactionMatchResult with all valid species pairs
        """
        ...
    
    def _create_result(
        self,
        class_reaction: ClassReaction,
        reactants: List[LipidStructure],
        products: List[LipidStructure],
    ) -> ReactionMatchResult:
        """Helper to create a base result with metadata."""
        return ReactionMatchResult(
            class_reaction=class_reaction,
            reactant_species_found=len(reactants),
            product_species_found=len(products),
        )
    
    def _create_pair(
        self,
        reactant: LipidStructure,
        product: LipidStructure,
        class_reaction: ClassReaction,
        required_compound: Optional[AcylChain] = None,
    ) -> SpeciesReactionPair:
        """Helper to create a species pair."""
        return SpeciesReactionPair(
            reactant=reactant,
            product=product,
            class_reaction=class_reaction,
            required_compound=required_compound,
        )
    
    def _get_chain_difference(
        self,
        larger: LipidStructure,
        smaller: LipidStructure,
    ) -> Tuple[int, int]:
        """Calculate the carbon/bond difference between two species.
        
        Returns (carbon_diff, bond_diff) where larger - smaller.
        """
        carbon_diff = larger.total_carbons - smaller.total_carbons
        bond_diff = larger.total_double_bonds - smaller.total_double_bonds
        return carbon_diff, bond_diff
    
    def _chains_match(
        self,
        species1: LipidStructure,
        species2: LipidStructure,
    ) -> bool:
        """Check if two species have the same total chain composition."""
        return (species1.total_carbons == species2.total_carbons and
                species1.total_double_bonds == species2.total_double_bonds)


def create_matcher_context(
    lipid_names: List[str],
    fa_names: Optional[List[str]] = None,
    facoa_names: Optional[List[str]] = None,
) -> MatcherContext:
    """Create a MatcherContext from lipid names.
    
    Args:
        lipid_names: List of lipid species names (e.g., ['PC 34:1', 'PE 36:2'])
        fa_names: Optional list of FA names (e.g., ['FA 16:0', 'FA 18:1'])
        facoa_names: Optional list of FACoA/CAR names (e.g., ['CAR 16:0'])
        
    Returns:
        MatcherContext populated with parsed species
    """
    parser = ChainParser()
    
    # Group species by class
    species_by_class: Dict[str, List[LipidStructure]] = {}
    all_names = set()
    for name in lipid_names:
        all_names.add(name)
        species = parser.parse(name)
        if species:
            headgroup = species.headgroup
            if headgroup not in species_by_class:
                species_by_class[headgroup] = []
            species_by_class[headgroup].append(species)
    
    # Parse FA species (handles both "FA 16:0" and "FA(16:0)" formats)
    fa_dict: Dict[str, AcylChain] = {}
    if fa_names:
        for name in fa_names:
            species = parser.parse(name)
            if species and species.headgroup == "FA" and species.chains:
                chain = species.chains[0]
                fa_dict[chain.notation] = chain
    
    # Parse FACoA/CoA species (handles multiple formats)
    facoa_dict: Dict[str, AcylChain] = {}
    if facoa_names:
        for name in facoa_names:
            # Try parsing directly (handles "CAR 16:0" format)
            species = parser.parse(name)
            if species and species.headgroup in ("CoA", "FACoA", "FaCoA", "FACOA") and species.chains:
                chain = species.chains[0]
                facoa_dict[chain.notation] = chain
            elif "FaCoA" in name or "FACoA" in name:
                # Legacy format: FACoA(16:0) -> parse as FA
                species = parser.parse(name.replace("FACoA", "FA").replace("FaCoA", "FA"))
                if species and species.chains:
                    chain = species.chains[0]
                    facoa_dict[chain.notation] = chain
    
    return MatcherContext(
        species_by_class=species_by_class,
        fa_species=fa_dict,
        facoa_species=facoa_dict,
        all_lipid_names=all_names,
    )
