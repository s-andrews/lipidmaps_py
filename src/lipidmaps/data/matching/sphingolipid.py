"""
Sphingolipid reaction matchers.

Handles sphingolipid-specific reactions involving sphingoid backbones:
- SPB <-> SPBP (phosphorylation, same structure)
- Cer -> SPB (FA release from ceramide)
- SPB -> Cer (FACoA addition to sphingoid base)
"""

from typing import List

from ..models.species_reaction import (
    ClassReaction,
    ReactionMatchResult,
    ReactionType,
)
from ..utils.chain_parser import AcylChain, LipidStructure
from .base import MatcherContext, ReactionMatcher


# Sphingoid base classes
# TODO: Check if these are the correct class names in the dataset (SPB vs SPH, etc.)
# TODO: Check if they have LM_ID's and they have their own group extractable from lm_id
SPHINGOID_BASES = {"SPB", "SPBP", "dhSPB", "dhSPBP", "S1P", "dhS1P"}

# Ceramide and related classes
CERAMIDE_CLASSES = {
    "Cer", "dhCer", "SM", "dhSM", "HexCer", "dhHexCer",
    "Hex2Cer", "GlcCer", "GalCer", "LacCer"
}


class SphingoSameStructureMatcher(ReactionMatcher):
    """Matches sphingolipid reactions with same structure.
    
    SPB <-> SPBP: Phosphorylation adds/removes phosphate
    dhSPB <-> SPB: Desaturation (handled by DesaturationMatcher)
    
    For SPB(d18:1) -> SPBP(d18:1), the backbone stays the same.
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.SAME_STRUCTURE
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match sphingolipids by backbone composition."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        
        if not reactants or not products:
            return result
        
        # Index products by backbone composition
        products_by_backbone: dict[tuple[int, int], List[LipidStructure]] = {}
        for prod in products:
            # Use backbone if available, else total composition
            if prod.backbone:
                key = (prod.backbone.carbons, prod.backbone.double_bonds)
            else:
                key = (prod.total_carbons, prod.total_double_bonds)
            
            if key not in products_by_backbone:
                products_by_backbone[key] = []
            products_by_backbone[key].append(prod)
        
        # Match by backbone
        for reactant in reactants:
            if reactant.backbone:
                key = (reactant.backbone.carbons, reactant.backbone.double_bonds)
            else:
                key = (reactant.total_carbons, reactant.total_double_bonds)
            
            matching_products = products_by_backbone.get(key, [])
            for product in matching_products:
                pair = self._create_pair(reactant, product, class_reaction)
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result


class SphingoFAReleaseMatcher(ReactionMatcher):
    """Matches Cer -> SPB type reactions (FA release).
    
    For Cer(d18:1/16:0) -> SPB(d18:1):
    - Ceramide loses its N-acyl FA chain
    - FA(16:0) must exist in dataset
    - Product backbone matches ceramide backbone
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.SPHINGO_FA_RELEASE
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match ceramide -> sphingoid base reactions."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        result.available_fa = set(context.fa_species.keys())
        
        if not reactants or not products:
            return result
        
        for reactant in reactants:
            matches = self._find_sphingo_fa_release_matches(
                reactant, products, context
            )
            for product, fa in matches:
                pair = self._create_pair(
                    reactant, product, class_reaction,
                    required_compound=fa
                )
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result
    
    def _find_sphingo_fa_release_matches(
        self,
        ceramide: LipidStructure,
        spb_products: List[LipidStructure],
        context: MatcherContext,
    ) -> List[tuple[LipidStructure, AcylChain]]:
        """Find SPB products from ceramide FA release."""
        matches = []
        
        # If ceramide has full structure Cer(d18:1/16:0)
        if ceramide.backbone and ceramide.chains:
            # The N-acyl chain is released as FA
            for chain in ceramide.chains:
                if not context.has_fa(chain):
                    continue
                
                # Product should match the backbone
                for product in spb_products:
                    backbone_match = False
                    if product.backbone:
                        backbone_match = (
                            product.backbone.carbons == ceramide.backbone.carbons and
                            product.backbone.double_bonds == ceramide.backbone.double_bonds
                        )
                    else:
                        # Match by total if no backbone info
                        backbone_match = (
                            product.total_carbons == ceramide.backbone.carbons and
                            product.total_double_bonds == ceramide.backbone.double_bonds
                        )
                    
                    if backbone_match:
                        matches.append((product, chain))
        
        else:
            # Sum composition matching for sphingolipids
            # Cer(34:1) -> SPB + FA: need to find valid backbone/FA splits
            for product in spb_products:
                prod_carbons = product.total_carbons
                prod_bonds = product.total_double_bonds
                
                fa_carbons = ceramide.total_carbons - prod_carbons
                fa_bonds = ceramide.total_double_bonds - prod_bonds
                
                if fa_carbons > 0 and fa_bonds >= 0:
                    fa = context.fa_for_difference(fa_carbons, fa_bonds)
                    if fa is not None:
                        matches.append((product, fa))
        
        return matches


class SphingoFACoAAdditionMatcher(ReactionMatcher):
    """Matches SPB -> Cer type reactions (FACoA addition).
    
    For SPB(d18:1) + FACoA(16:0) -> Cer(d18:1/16:0):
    - Sphingoid base gains N-acyl chain
    - FACoA(16:0) must exist in dataset
    - Product backbone matches SPB
    """
    
    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.SPHINGO_FACOA_ADD
    
    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        """Match sphingoid base -> ceramide reactions."""
        
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)
        
        result = self._create_result(class_reaction, reactants, products)
        result.available_facoa = set(context.facoa_species.keys())
        
        if not reactants or not products:
            return result
        
        for reactant in reactants:
            matches = self._find_sphingo_facoa_addition_matches(
                reactant, products, context
            )
            for product, facoa in matches:
                pair = self._create_pair(
                    reactant, product, class_reaction,
                    required_compound=facoa
                )
                result.pairs.append(pair)
        
        result.pairs_matched = len(result.pairs)
        return result
    
    def _find_sphingo_facoa_addition_matches(
        self,
        spb: LipidStructure,
        ceramide_products: List[LipidStructure],
        context: MatcherContext,
    ) -> List[tuple[LipidStructure, AcylChain]]:
        """Find ceramide products from SPB + FACoA."""
        matches = []
        
        spb_backbone = spb.backbone
        spb_carbons = spb_backbone.carbons if spb_backbone else spb.total_carbons
        spb_bonds = spb_backbone.double_bonds if spb_backbone else spb.total_double_bonds
        
        for ceramide in ceramide_products:
            # If ceramide has full structure
            if ceramide.backbone and ceramide.chains:
                # Check backbone matches SPB
                if (ceramide.backbone.carbons == spb_carbons and
                    ceramide.backbone.double_bonds == spb_bonds):
                    # N-acyl chain should match available FACoA
                    for chain in ceramide.chains:
                        if context.has_facoa(chain):
                            matches.append((ceramide, chain))
                            break  # One match per ceramide
            else:
                # Sum composition matching
                facoa_carbons = ceramide.total_carbons - spb_carbons
                facoa_bonds = ceramide.total_double_bonds - spb_bonds
                
                if facoa_carbons > 0 and facoa_bonds >= 0:
                    facoa = context.facoa_for_difference(facoa_carbons, facoa_bonds)
                    if facoa is not None:
                        matches.append((ceramide, facoa))
        
        return matches
