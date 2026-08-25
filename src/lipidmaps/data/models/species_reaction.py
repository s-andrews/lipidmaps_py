"""
Species-level reaction models for pathway analysis.

These models represent molecular species reaction pairs generated from
class-level reaction rules.
"""

from enum import Enum
from typing import Any, Dict, List, Optional, Set

from pydantic import Field, computed_field

from .base import LipidmapsBaseModel
from ..utils.chain_parser import LipidStructure, AcylChain


class CompoundRequirement(str, Enum):
    """Type of compound required for a reaction."""
    NONE = "none"       # No compound required (same structure)
    FA = "fa"           # Requires free fatty acid (release)
    FACOA = "facoa"     # Requires fatty acyl-CoA (addition)


class ReactionType(str, Enum):
    """Category of reaction for matching strategy selection."""
    SAME_STRUCTURE = "same_structure"       # No FA/FACoA involved
    FA_RELEASE = "fa_release"               # Releases FA (e.g., PC -> LPC)
    FACOA_ADDITION = "facoa_addition"       # Adds FACoA (e.g., LPC -> PC)
    SPHINGO_FA_RELEASE = "sphingo_fa_release"     # Sphingolipid FA release (Cer -> SPB)
    SPHINGO_FACOA_ADD = "sphingo_facoa_addition"  # Sphingolipid FACoA add (SPB -> Cer)
    CARDIOLIPIN = "cardiolipin"             # Special CL handling (all-to-all)
    DESATURATION = "desaturation"           # dhCer -> Cer type reactions
    NAMED_COMPOUND = "named_compound"        # Pair by compound identity, not composition


class ClassReaction(LipidmapsBaseModel):
    """Class-level reaction definition for pathway analysis.
    
    Represents a reaction between lipid classes (e.g., PC -> PA)
    with metadata about the reaction mechanism.
    """
    reactant_class: str = Field(..., description="Reactant lipid class (e.g., 'PC')")
    product_class: str = Field(..., description="Product lipid class (e.g., 'PA')")
    
    # Reaction metadata
    reaction_class: str = Field(default="", description="Classification (e.g., 'Glycerophospholipids')")
    sub_class: str = Field(default="", description="Sub-classification")
    
    # Mechanism
    compound_require: CompoundRequirement = Field(default=CompoundRequirement.NONE)
    acyl_add: bool = Field(default=False, description="True if reaction adds acyl chain")
    is_molspecies: bool = Field(default=False, description="True if described at molecular species level")

    # Explicit override: when set, force NAMED_COMPOUND routing regardless of the
    # class names. Used when extraction has already established (by lm_id) that
    # both endpoints are specific measured compounds with no generic form -- e.g.
    # sterol intermediates or oxylipins whose measured names ("27-Hydroxy-
    # cholesterol") are not in NAMED_COMPOUND_CLASSES and would otherwise fall
    # through to composition-based matching and mis-pair at the degenerate level.
    named_compound: bool = Field(default=False, description="Force identity-based (named-compound) matching")

    # Specific LIPID MAPS IDs of the endpoints for named-compound reactions. These
    # drive node shape/coloring and structure links, which the class-headgroup
    # lookup cannot supply for specific compounds.
    reactant_lm_id: Optional[str] = Field(default=None, description="Reactant LM ID (named-compound reactions)")
    product_lm_id: Optional[str] = Field(default=None, description="Product LM ID (named-compound reactions)")

    # Associated genes
    genes: List[str] = Field(default_factory=list)
    
    @computed_field
    @property
    def reaction_key(self) -> str:
        """Reaction key (e.g., 'PC,PA')."""
        return f"{self.reactant_class},{self.product_class}"
    
    @computed_field
    @property
    def reaction_type(self) -> ReactionType:
        """Determine the reaction type for matcher selection."""
        # Extraction may have already proven (by lm_id) that both endpoints are
        # specific measured named compounds; honour that first so routing does
        # not depend on the endpoint names appearing in NAMED_COMPOUND_CLASSES.
        if self.named_compound:
            return ReactionType.NAMED_COMPOUND

        # Check for cardiolipin
        if "CL" in self.reactant_class or "CL" in self.product_class:
            return ReactionType.CARDIOLIPIN
        
        # Check for sphingoid base reactions
        sphingo_bases = {"SPB", "SPBP", "dhSPB", "dhSPBP"}
        re_is_sphingo = self.reactant_class in sphingo_bases
        pro_is_sphingo = self.product_class in sphingo_bases
        
        if re_is_sphingo or pro_is_sphingo:
            if self.compound_require == CompoundRequirement.FA:
                return ReactionType.SPHINGO_FA_RELEASE
            elif self.compound_require == CompoundRequirement.FACOA:
                return ReactionType.SPHINGO_FACOA_ADD
        
        # Check for desaturation (dh -> non-dh)
        if (self.reactant_class.startswith("dh") and 
            not self.product_class.startswith("dh") and
            self.reactant_class[2:] == self.product_class):
            return ReactionType.DESATURATION
        
        # Check compound requirements
        if self.compound_require == CompoundRequirement.NONE:
            return ReactionType.SAME_STRUCTURE
        elif self.compound_require == CompoundRequirement.FA:
            return ReactionType.FA_RELEASE
        else:  # FACOA
            return ReactionType.FACOA_ADDITION


class SpeciesReactionPair(LipidmapsBaseModel):
    """A single reactant-product pair at the molecular species level.
    
    Represents one specific reaction, e.g., PC(34:1) -> PA(34:1).
    """
    reactant: LipidStructure = Field(..., description="Reactant molecular species")
    product: LipidStructure = Field(..., description="Product molecular species")
    
    # Parent class reaction
    class_reaction: ClassReaction = Field(..., description="Parent class-level reaction")
    
    # Validation status
    is_valid: bool = Field(default=True, description="Whether this pair is biochemically valid")
    validation_notes: List[str] = Field(default_factory=list)
    
    # Required compound (if any)
    required_compound: Optional[AcylChain] = Field(default=None, 
        description="FA/FACoA required for this specific reaction")
    
    @computed_field
    @property
    def pair_id(self) -> str:
        """Unique identifier for this pair."""
        def clean(s: str) -> str:
            return s.replace("(", "").replace(")", "").replace(":", "")
        return f"{clean(self.reactant.full_name)}{clean(self.product.full_name)}"


class ReactionMatchResult(LipidmapsBaseModel):
    """Result of matching a class reaction to molecular species pairs.
    
    Contains all valid species-level reaction pairs for a given class reaction.
    """
    class_reaction: ClassReaction = Field(..., description="The class-level reaction")
    pairs: List[SpeciesReactionPair] = Field(default_factory=list)
    
    # Matching metadata
    reactant_species_found: int = Field(default=0, description="Number of reactant species in dataset")
    product_species_found: int = Field(default=0, description="Number of product species in dataset")
    pairs_matched: int = Field(default=0, description="Number of valid pairs generated")
    
    # Available compounds (for FA/FACoA reactions)
    available_fa: Set[str] = Field(default_factory=set, description="Available FA acyl chains")
    available_facoa: Set[str] = Field(default_factory=set, description="Available FACoA acyl chains")
    
    @computed_field
    @property
    def reaction_key(self) -> str:
        """Key for file naming."""
        return self.class_reaction.reaction_key.lower()
    
    @computed_field
    @property
    def has_pairs(self) -> bool:
        """Whether any valid pairs were found."""
        return len(self.pairs) > 0
    
    def get_grouped_by_reactant(self) -> List[tuple[LipidStructure, List[LipidStructure]]]:
        """Group product species by their reactant species.
        
        Returns list of (reactant, [products]) tuples.
        """
        grouped: Dict[str, tuple[LipidStructure, List[LipidStructure]]] = {}
        for pair in self.pairs:
            key = pair.reactant.full_name
            if key not in grouped:
                grouped[key] = (pair.reactant, [])
            grouped[key][1].append(pair.product)
        return list(grouped.values())


class PathwayReactionSet(LipidmapsBaseModel):
    """Collection of all reaction match results for a dataset.
    
    This is the top-level container for pathway analysis.
    """
    results: Dict[str, ReactionMatchResult] = Field(default_factory=dict,
        description="Match results keyed by reaction key (e.g., 'pc,pa')")
    
    # Dataset metadata
    dataset_lipid_count: int = Field(default=0)
    reactions_checked: int = Field(default=0)
    reactions_with_pairs: int = Field(default=0)
    
    # Available compounds
    fa_species: List[str] = Field(default_factory=list)
    facoa_species: List[str] = Field(default_factory=list)
    
    def add_result(self, result: ReactionMatchResult) -> None:
        """Add a match result to the set."""
        self.results[result.reaction_key] = result
        if result.has_pairs:
            self.reactions_with_pairs += 1
    
    def get_result(self, reactant: str, product: str) -> Optional[ReactionMatchResult]:
        """Get result for a specific class reaction."""
        key = f"{reactant},{product}".lower()
        return self.results.get(key)
    
    def get_all_pairs(self) -> List[SpeciesReactionPair]:
        """Get all species pairs across all reactions."""
        pairs = []
        for result in self.results.values():
            pairs.extend(result.pairs)
        return pairs
    
    def summary(self) -> Dict[str, Any]:
        """Generate summary statistics."""
        return {
            "dataset_lipids": self.dataset_lipid_count,
            "reactions_checked": self.reactions_checked,
            "reactions_with_matches": self.reactions_with_pairs,
            "total_species_pairs": sum(len(r.pairs) for r in self.results.values()),
            "fa_species_count": len(self.fa_species),
            "facoa_species_count": len(self.facoa_species),
        }
