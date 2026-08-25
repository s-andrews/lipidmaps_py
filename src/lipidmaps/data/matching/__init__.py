"""
Reaction matching package.

Provides strategy-pattern matchers for molecular species
reaction pairing.
"""

from .base import MatcherContext, ReactionMatcher, create_matcher_context
from .same_structure import SameStructureMatcher, DesaturationMatcher
from .fa_compound import FACompoundMatcher, FACompoundMatcherFullStructure
from .facoa_compound import FACoACompoundMatcher, FACoACompoundMatcherFullStructure
from .sphingolipid import (
    SphingoSameStructureMatcher,
    SphingoFAReleaseMatcher,
    SphingoFACoAAdditionMatcher,
)
from .cardiolipin import CardiolipinMatcher, CardiolipinCompositionMatcher
from .named_compound import NamedCompoundMatcher
from .registry import (
    MatcherRegistry,
    create_default_registry,
    match_pathway_reactions,
)

__all__ = [
    # Base
    "MatcherContext",
    "ReactionMatcher",
    "create_matcher_context",
    # Same structure
    "SameStructureMatcher",
    "DesaturationMatcher",
    # FA compound
    "FACompoundMatcher",
    "FACompoundMatcherFullStructure",
    # FACoA compound
    "FACoACompoundMatcher",
    "FACoACompoundMatcherFullStructure",
    # Sphingolipid
    "SphingoSameStructureMatcher",
    "SphingoFAReleaseMatcher",
    "SphingoFACoAAdditionMatcher",
    # Cardiolipin
    "CardiolipinMatcher",
    "CardiolipinCompositionMatcher",
    # Named compound (sterol / shunt)
    "NamedCompoundMatcher",
    # Registry
    "MatcherRegistry",
    "create_default_registry",
    "match_pathway_reactions",
]
