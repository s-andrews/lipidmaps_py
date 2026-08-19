"""
Matcher registry and factory.

Provides automatic selection of the appropriate reaction matcher
based on reaction type (strategy pattern).
"""

from typing import Dict, List, Optional

from ..models.species_reaction import (
    ClassReaction,
    PathwayReactionSet,
    ReactionMatchResult,
    ReactionType,
)
from .base import MatcherContext, ReactionMatcher, create_matcher_context
from .same_structure import DesaturationMatcher, SameStructureMatcher
from .fa_compound import FACompoundMatcher, FACompoundMatcherFullStructure
from .facoa_compound import FACoACompoundMatcher, FACoACompoundMatcherFullStructure
from .sphingolipid import (
    SphingoFACoAAdditionMatcher,
    SphingoFAReleaseMatcher,
)
from .cardiolipin import CardiolipinMatcher


class MatcherRegistry:
    """Registry of reaction matchers indexed by reaction type.
    
    The registry allows registering custom matchers and provides
    automatic matcher selection based on reaction classification.
    """
    
    def __init__(self, use_full_structure: bool = False):
        """Initialize registry with default matchers.
        
        Args:
            use_full_structure: If True, use structure-aware matchers
                when available for FA/FACoA reactions.
        """
        self._matchers: Dict[ReactionType, ReactionMatcher] = {}
        self._use_full_structure = use_full_structure
        self._register_defaults()
    
    def _register_defaults(self) -> None:
        """Register the default set of matchers."""
        # Same structure
        self.register(ReactionType.SAME_STRUCTURE, SameStructureMatcher())
        self.register(ReactionType.DESATURATION, DesaturationMatcher())
        
        # FA/FACoA compound reactions
        if self._use_full_structure:
            self.register(ReactionType.FA_RELEASE, FACompoundMatcherFullStructure())
            self.register(ReactionType.FACOA_ADDITION, FACoACompoundMatcherFullStructure())
        else:
            self.register(ReactionType.FA_RELEASE, FACompoundMatcher())
            self.register(ReactionType.FACOA_ADDITION, FACoACompoundMatcher())
        
        # Sphingolipid reactions
        self.register(ReactionType.SPHINGO_FA_RELEASE, SphingoFAReleaseMatcher())
        self.register(ReactionType.SPHINGO_FACOA_ADD, SphingoFACoAAdditionMatcher())
        
        # Cardiolipin
        self.register(ReactionType.CARDIOLIPIN, CardiolipinMatcher())
    
    def register(self, reaction_type: ReactionType, matcher: ReactionMatcher) -> None:
        """Register a matcher for a reaction type.
        
        Args:
            reaction_type: The reaction type to handle
            matcher: The matcher implementation
        """
        self._matchers[reaction_type] = matcher
    
    def get_matcher(self, reaction_type: ReactionType) -> Optional[ReactionMatcher]:
        """Get the registered matcher for a reaction type."""
        return self._matchers.get(reaction_type)
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match species pairs using the appropriate matcher.
        
        Automatically selects the matcher based on reaction type.
        
        Args:
            class_reaction: The class-level reaction to match
            context: Dataset context with available species
            
        Returns:
            ReactionMatchResult with all valid pairs
        """
        reaction_type = class_reaction.reaction_type
        matcher = self.get_matcher(reaction_type)
        
        if matcher is None:
            # Return empty result if no matcher registered
            return ReactionMatchResult(
                class_reaction=class_reaction,
                pairs=[],
            )
        
        return matcher.match(class_reaction, context)
    
    def match_all(
        self,
        reactions: List[ClassReaction],
        context: MatcherContext,
    ) -> PathwayReactionSet:
        """Match all class reactions against the dataset.
        
        Args:
            reactions: List of class-level reactions to process
            context: Dataset context with available species
            
        Returns:
            PathwayReactionSet containing all results
        """
        result_set = PathwayReactionSet(
            reactions_checked=len(reactions),
            fa_species=list(context.fa_species.keys()),
            facoa_species=list(context.facoa_species.keys()),
        )
        
        for class_reaction in reactions:
            result = self.match(class_reaction, context)
            result_set.add_result(result)
        
        return result_set


def create_default_registry(use_full_structure: bool = False) -> MatcherRegistry:
    """Create a matcher registry with default configuration.
    
    Args:
        use_full_structure: If True, use structure-aware matchers
            for more precise FA/FACoA matching.
            
    Returns:
        Configured MatcherRegistry
    """
    return MatcherRegistry(use_full_structure=use_full_structure)


def match_pathway_reactions(
    lipid_names: List[str],
    reactions: List[ClassReaction],
    fa_names: Optional[List[str]] = None,
    facoa_names: Optional[List[str]] = None,
    use_full_structure: bool = False,
) -> PathwayReactionSet:
    """Convenience function to match all reactions for a dataset.
    
    This is the main entry point for pathway analysis.
    
    Args:
        lipid_names: List of lipid species in the dataset
        reactions: List of class-level reactions to check
        fa_names: Optional list of FA species
        facoa_names: Optional list of FACoA species
        use_full_structure: Whether to use structure-aware matching
        
    Returns:
        PathwayReactionSet with all match results
        
    Example:
        >>> lipids = ["PC(34:1)", "PA(34:1)", "PC(36:2)", "PA(36:2)"]
        >>> reactions = [ClassReaction(reactant_class="PC", product_class="PA")]
        >>> results = match_pathway_reactions(lipids, reactions)
        >>> print(results.summary())
    """
    context = create_matcher_context(
        lipid_names=lipid_names,
        fa_names=fa_names,
        facoa_names=facoa_names,
    )
    
    registry = create_default_registry(use_full_structure=use_full_structure)
    
    result_set = registry.match_all(reactions, context)
    result_set.dataset_lipid_count = len(lipid_names)
    
    return result_set
