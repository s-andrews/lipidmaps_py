"""
Chain parser for lipid species nomenclature.

Parses lipid species names into structured components for reaction matching.
Handles sum composition (e.g., PC 34:1) and full structure (e.g., PC 16:0_18:1).
"""

import re
from enum import Enum
from typing import List, Optional, Tuple, TYPE_CHECKING

from pydantic import BaseModel, Field, field_validator, computed_field
from .headgroups import lipidmaps_headgroups

if TYPE_CHECKING:
    from lipidmaps.data.models.reaction import ReactionData


class StructureLevel(str, Enum):
    """Level of structural detail in lipid nomenclature."""
    SUM = "sum"           # Total carbons:bonds (e.g., PC 34:1)
    FULL = "full"         # Individual chains (e.g., PC 16:0_18:1)
    MOLECULAR = "molecular"  # Named species without structure (e.g., ST, cholesterol)

class AcylChain(BaseModel):
    """Represents a single fatty acyl chain."""
    carbons: int = Field(..., ge=0, description="Number of carbon atoms")
    double_bonds: int = Field(..., ge=0, description="Number of double bonds")
    modifier: Optional[str] = Field(default=None, description="Chain modifier (O-, P-, d, t, etc.)")
    
    @computed_field
    @property
    def notation(self) -> str:
        """Return standard notation (e.g., '16:0', '18:1')."""
        prefix = f"{self.modifier}" if self.modifier else ""
        return f"{prefix}{self.carbons}:{self.double_bonds}"
    
    def __str__(self) -> str:
        return self.notation
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AcylChain):
            return False
        return (self.carbons == other.carbons and 
                self.double_bonds == other.double_bonds and
                self.modifier == other.modifier)
    
    def __hash__(self) -> int:
        return hash((self.carbons, self.double_bonds, self.modifier))


class SphingoBackbone(BaseModel):
    """Represents a sphingoid base backbone (e.g., d18:1, t18:0)."""
    prefix: str = Field(..., description="Backbone type (d=dihydroxy, t=trihydroxy)")
    carbons: int = Field(..., ge=0)
    double_bonds: int = Field(..., ge=0)
    hydroxyl_info: Optional[str] = Field(default=None, description="Hydroxylation info (e.g., ;O2)")
    
    @computed_field
    @property
    def notation(self) -> str:
        suffix = self.hydroxyl_info or ""
        return f"{self.prefix}{self.carbons}:{self.double_bonds}{suffix}"
    
    def __str__(self) -> str:
        return self.notation


class LipidStructure(BaseModel):
    """Parsed structural components of a lipid name."""
    input_name: str = Field(..., description="Original input name")
    headgroup: str = Field(..., description="Lipid class/headgroup (e.g., PC, PE, Cer)")
    
    # Structure components
    chains: List[AcylChain] = Field(default_factory=list, description="Individual acyl chains")
    backbone: Optional[SphingoBackbone] = Field(default=None, description="Sphingoid backbone if present")
    
    # Sum composition (always populated)
    total_carbons: int = Field(default=0, ge=0)
    total_double_bonds: int = Field(default=0, ge=0)
    
    # Metadata
    level: StructureLevel = Field(default=StructureLevel.SUM)
    linkage_type: Optional[str] = Field(default=None, description="ester, ether_alkyl, ether_vinyl, amide")
    
    @computed_field
    @property
    def sum_composition(self) -> str:
        """Return sum composition notation (e.g., '34:1')."""
        return f"{self.total_carbons}:{self.total_double_bonds}"
    
    @computed_field
    @property 
    def full_name(self) -> str:
        """Reconstruct full species name."""
        if self.level == StructureLevel.MOLECULAR:
            return self.headgroup
        
        if self.backbone and self.chains:
            chain_str = "/".join(str(c) for c in self.chains)
            return f"{self.headgroup} {self.backbone}/{chain_str}"
        elif self.chains and self.level == StructureLevel.FULL:
            chain_str = "_".join(str(c) for c in self.chains)
            return f"{self.headgroup} {chain_str}"
        else:
            return f"{self.headgroup} {self.sum_composition}"
    
    def matches_sum_composition(self, other: "LipidStructure") -> bool:
        """Check if two structures have matching sum composition."""
        return (self.total_carbons == other.total_carbons and
                self.total_double_bonds == other.total_double_bonds)
    
    def get_acyl_chain_notation(self) -> str:
        """Get the acyl chain portion without headgroup."""
        if self.level == StructureLevel.FULL and self.chains:
            return "_".join(str(c) for c in self.chains)
        return self.sum_composition


class ChainParser:
    """Parser for lipid species nomenclature."""
    
    # Headgroup patterns - order matters (longer matches first)
    HEADGROUPS = list(lipidmaps_headgroups.keys())

    # Sphingoid backbone prefixes
    SPHINGO_PREFIXES = {"d", "t", "m"}
    
    # Linkage type mappings
    LINKAGE_MAP = {
        "O-": "ether_alkyl",
        "P-": "ether_vinyl",
    }
    
    # Sphingolipid headgroups (use amide linkage)
    SPHINGO_HEADGROUPS = {
        "Cer", "dhCer", "SM", "dhSM", "Cer1P", "dhCer1P",
        "Hex-Cer", "Hex-dhCer", "Glc-Cer", "Gal-Cer",
        "Lac-Cer", "Lac-dhCer",
    }
    
    @classmethod
    def parse(cls, name: str) -> LipidStructure:
        """Parse a lipid species name into structured components.
        
        Examples:
            parse("PC 34:1") -> sum composition
            parse("PC 16:0_18:1") -> full structure
            parse("Cer d18:1/18:0") -> sphingolipid with backbone
            parse("ST") -> molecular species
        """
        name = name.strip()
        
        # Try to extract headgroup
        headgroup = cls._extract_headgroup(name)
        if not headgroup:
            # Might be a molecular species without parentheses
            return LipidStructure(
                input_name=name,
                headgroup=name,
                level=StructureLevel.MOLECULAR,
            )
        
        # Determine linkage type
        linkage = cls._determine_linkage(headgroup)
        
        # Extract structure portion (inside parentheses)
        structure = cls._extract_structure(name, headgroup)
        
        if not structure:
            # No structure info (e.g., just "PC" or "SPB")
            return LipidStructure(
                input_name=name,
                headgroup=headgroup,
                level=StructureLevel.MOLECULAR,
                linkage_type=linkage,
            )
        
        # Parse the structure
        backbone, chains, level = cls._parse_structure(structure, headgroup)
        
        # Calculate totals
        total_c, total_db = cls._calculate_totals(backbone, chains, structure)
        
        return LipidStructure(
            input_name=name,
            headgroup=headgroup,
            chains=chains,
            backbone=backbone,
            total_carbons=total_c,
            total_double_bonds=total_db,
            level=level,
            linkage_type=linkage,
        )
    
    @classmethod
    def _extract_headgroup(cls, name: str) -> Optional[str]:
        """Extract the headgroup from a lipid name."""
        # Sort by length descending to match longer patterns first
        if "coa" in name.lower():
            # Handle CoA species (e.g., "CoA 16:0", "FACoA(16:0)")
            return "acyl CoA"
        for hg in sorted(cls.HEADGROUPS, key=len, reverse=True):
            if name.startswith(hg):
                # Accept direct structural notation after prefixed headgroups such as
                # 'PC P-16:0_18:1' or 'PC O-34:1' in addition to spaced/parenthetical forms.
                remainder = name[len(hg):]
                if (
                    not remainder
                    or remainder.startswith("(")
                    or remainder.startswith(" ")
                    or remainder[0].isdigit()
                ):
                    return hg
        return None
    
    @classmethod
    def _extract_structure(cls, name: str, headgroup: Optional[str] = None) -> Optional[str]:
        """Extract the structure portion from parentheses or space-separated.
        
        Handles both formats:
        - PC 34:1 -> 34:1
        - FA 16:0 -> 16:0
        - PC(34:1) -> 34:1 (legacy format)
        """
        # First try parentheses format
        match = re.search(r"\(([^)]+)\)", name)
        if match:
            return match.group(1)

        # If we already know the headgroup, strip it directly so prefixed forms like
        # 'PC P-16:0_16:1' keep the correct headgroup and only expose the chain text.
        if headgroup and name.startswith(headgroup):
            remainder = name[len(headgroup):].strip()
            if remainder:
                return remainder
        
        # Try space-separated format (e.g., "FA 16:0")
        match = re.search(r"^[A-Za-z-]+\s+(\S+)$", name)
        if match:
            return match.group(1)
        
        return None
    
    @classmethod
    def _determine_linkage(cls, headgroup: str) -> str:
        """Determine linkage type from headgroup."""
        stripped_headgroup = headgroup.strip()

        if stripped_headgroup.startswith("O-") or stripped_headgroup.endswith(" O-"):
            return "ether_alkyl"
        if stripped_headgroup.startswith("P-") or stripped_headgroup.endswith(" P-"):
            return "ether_vinyl"
        
        # Strip prefix for sphingolipid check
        base_hg = stripped_headgroup.removeprefix("O-").removeprefix("P-")
        if base_hg.endswith(" O-") or base_hg.endswith(" P-"):
            base_hg = base_hg[:-3]
        if base_hg in cls.SPHINGO_HEADGROUPS:
            return "amide"
        
        return "ester"
    
    @classmethod
    def _parse_structure(cls, structure: str, headgroup: str) -> Tuple[Optional[SphingoBackbone], List[AcylChain], StructureLevel]:
        """Parse structure string into backbone and chains."""
        backbone = None
        chains: List[AcylChain] = []
        level = StructureLevel.SUM
        
        # Check for sphingoid backbone (d18:1, t18:0, etc.)
        sphingo_match = re.match(r"^([dDtTmM])(\d+):(\d+)(;O\d+)?", structure)
        if sphingo_match:
            backbone = SphingoBackbone(
                prefix=sphingo_match.group(1).lower(),
                carbons=int(sphingo_match.group(2)),
                double_bonds=int(sphingo_match.group(3)),
                hydroxyl_info=sphingo_match.group(4),
            )
            # Get the remainder after backbone
            remainder = structure[sphingo_match.end():]
            if remainder.startswith("/"):
                remainder = remainder[1:]
                level = StructureLevel.FULL
                chains = cls._parse_chains(remainder)
            return backbone, chains, level
        
        # Check for full structure (contains _ or / separators)
        if "_" in structure or "/" in structure:
            level = StructureLevel.FULL
            chains = cls._parse_chains(structure)
        else:
            # Sum composition
            level = StructureLevel.SUM
            chain = cls._parse_single_chain(structure)
            if chain:
                chains = [chain]
        
        return backbone, chains, level
    
    @classmethod
    def _parse_chains(cls, structure: str) -> List[AcylChain]:
        """Parse multiple chains from structure string."""
        chains = []
        # Split by _ or /
        parts = re.split(r"[_/]", structure)
        for part in parts:
            chain = cls._parse_single_chain(part.strip())
            if chain:
                chains.append(chain)
        return chains
    
    @classmethod
    def _parse_single_chain(cls, chain_str: str) -> Optional[AcylChain]:
        """Parse a single chain notation (e.g., '16:0', 'O-16:0')."""
        # Handle modifiers (O-, P-, etc.)
        modifier = None
        match = re.match(r"^([OPop]-)?(\d+):(\d+)", chain_str)
        if match:
            if match.group(1):
                modifier = match.group(1).upper().rstrip("-") + "-"
            return AcylChain(
                carbons=int(match.group(2)),
                double_bonds=int(match.group(3)),
                modifier=modifier,
            )
        return None
    
    @classmethod
    def _calculate_totals(cls, backbone: Optional[SphingoBackbone], 
                          chains: List[AcylChain], 
                          structure: str) -> Tuple[int, int]:
        """Calculate total carbons and double bonds."""
        total_c = 0
        total_db = 0
        
        if backbone:
            total_c += backbone.carbons
            total_db += backbone.double_bonds
        
        for chain in chains:
            total_c += chain.carbons
            total_db += chain.double_bonds
        
        # If we couldn't parse chains, try to extract from sum composition
        if not chains and not backbone:
            match = re.search(r"(\d+):(\d+)", structure)
            if match:
                total_c = int(match.group(1))
                total_db = int(match.group(2))
        
        return total_c, total_db


# Convenience function
def parse_lipid(name: str) -> LipidStructure:
    """Parse a lipid species name into structured components."""
    if name is None:
        return None
    return ChainParser.parse(name)





# =============================================================================
# Common Fatty Acids
# =============================================================================

# Standard fatty acids found in biological membranes
# Format: (carbons, double_bonds, common_name)
COMMON_FATTY_ACIDS: List[Tuple[int, int, str]] = [
    # Saturated fatty acids
    (12, 0, "lauric acid"),
    (14, 0, "myristic acid"),
    (16, 0, "palmitic acid"),
    (18, 0, "stearic acid"),
    (20, 0, "arachidic acid"),
    (22, 0, "behenic acid"),
    (24, 0, "lignoceric acid"),
    # Monounsaturated fatty acids
    (16, 1, "palmitoleic acid"),
    (18, 1, "oleic acid"),
    (20, 1, "gondoic acid"),
    (22, 1, "erucic acid"),
    (24, 1, "nervonic acid"),
    # Polyunsaturated fatty acids (omega-6)
    (18, 2, "linoleic acid"),
    (18, 3, "gamma-linolenic acid"),
    (20, 3, "dihomo-gamma-linolenic acid"),
    (20, 4, "arachidonic acid"),
    (22, 4, "adrenic acid"),
    # Polyunsaturated fatty acids (omega-3)
    (18, 3, "alpha-linolenic acid"),
    (20, 5, "eicosapentaenoic acid (EPA)"),
    (22, 5, "docosapentaenoic acid (DPA)"),
    (22, 6, "docosahexaenoic acid (DHA)"),
]

# Extended set including less common but biologically relevant FAs
EXTENDED_FATTY_ACIDS: List[Tuple[int, int, str]] = COMMON_FATTY_ACIDS + [
    # Short chain
    (8, 0, "caprylic acid"),
    (10, 0, "capric acid"),
    # Very long chain saturated
    (26, 0, "cerotic acid"),
    (28, 0, "montanic acid"),
    # Odd chain (bacterial, dietary)
    (15, 0, "pentadecanoic acid"),
    (17, 0, "margaric acid"),
    (17, 1, "heptadecenoic acid"),
    # Additional unsaturated
    (14, 1, "myristoleic acid"),
    (20, 2, "eicosadienoic acid"),
    (22, 2, "docosadienoic acid"),
]


def get_common_fa_chains() -> List[AcylChain]:
    """Get common fatty acids as AcylChain objects.
    
    Returns the most common fatty acids found in biological systems.
    """
    seen = set()
    chains = []
    for carbons, double_bonds, _ in COMMON_FATTY_ACIDS:
        key = (carbons, double_bonds)
        if key not in seen:
            seen.add(key)
            chains.append(AcylChain(carbons=carbons, double_bonds=double_bonds))
    return chains


def get_extended_fa_chains() -> List[AcylChain]:
    """Get extended set of fatty acids as AcylChain objects.
    
    Includes less common but biologically relevant fatty acids.
    """
    seen = set()
    chains = []
    for carbons, double_bonds, _ in EXTENDED_FATTY_ACIDS:
        key = (carbons, double_bonds)
        if key not in seen:
            seen.add(key)
            chains.append(AcylChain(carbons=carbons, double_bonds=double_bonds))
    return chains


def get_common_fa_names() -> List[str]:
    """Get common fatty acid names in LIPID MAPS shorthand format.
    
    Example: ['FA 12:0', 'FA 14:0', 'FA 16:0', ...]
    """
    seen = set()
    names = []
    for carbons, double_bonds, _ in COMMON_FATTY_ACIDS:
        name = f"FA {carbons}:{double_bonds}"
        if name not in seen:
            seen.add(name)
            names.append(name)
    return names


def get_common_facoa_names() -> List[str]:
    """Get common fatty acyl-CoA names in LIPID MAPS shorthand format.
    
    Example: ['CoA 12:0', 'CoA 14:0', 'CoA 16:0', ...]
    
    Note: Uses CoA as the LIPID MAPS shorthand for acyl-CoA.
    """
    seen = set()
    names = []
    for carbons, double_bonds, _ in COMMON_FATTY_ACIDS:
        name = f"CoA {carbons}:{double_bonds}"
        if name not in seen:
            seen.add(name)
            names.append(name)
    return names


def get_extended_fa_names() -> List[str]:
    """Get extended fatty acid names in LIPID MAPS shorthand format."""
    seen = set()
    names = []
    for carbons, double_bonds, _ in EXTENDED_FATTY_ACIDS:
        name = f"FA {carbons}:{double_bonds}"
        if name not in seen:
            seen.add(name)
            names.append(name)
    return names


def get_extended_facoa_names() -> List[str]:
    """Get extended fatty acyl-CoA names in LIPID MAPS shorthand format."""
    seen = set()
    names = []
    for carbons, double_bonds, _ in EXTENDED_FATTY_ACIDS:
        name = f"CAR {carbons}:{double_bonds}"
        if name not in seen:
            seen.add(name)
            names.append(name)
    return names


# Lipid classes whose sum-composition residue is itself a single fatty acyl
# chain. Datasets reported at sum level (e.g. 'LPI 18:0', 'CE 18:1') still carry
# a usable acyl chain for these classes, so they can seed the FA / acyl-CoA
# pools. This mirrors the legacy BioPAN R behaviour (parse_data.r), which
# derived FA/FaCoA acyl chains from lyso, CE and sphingolipid species.
SINGLE_ACYL_HEADGROUPS = {
    "LPA", "LPC", "LPE", "LPG", "LPI", "LPS",  # lyso-glycerophospholipids
    "MG",                                        # monoacylglycerol
    "CE",                                        # cholesteryl ester
}


def _strip_linkage_suffix(headgroup: str) -> str:
    """Strip a trailing ether linkage marker (' O-' / ' P-') from a headgroup."""
    base = headgroup.strip()
    for suffix in (" O-", " P-"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def infer_fa_from_lipids(lipid_names: List[str], include_mono_acyl: bool = True) -> List[str]:
    """Infer possible FA species from the acyl chains in a lipid dataset.

    Two sources are used:
    1. Every individual chain of full-structure species (e.g. PC 16:0_18:1).
       For sphingolipids the backbone is excluded by the parser, so only the
       N-acyl chain contributes (e.g. Cer d18:1/16:0 -> FA 16:0).
    2. The single residual chain of sum-level species whose class is itself
       mono-acyl (lyso-glycerophospholipids, MG, CE; see SINGLE_ACYL_HEADGROUPS).
       Disable this source with ``include_mono_acyl=False``.

    Sum-composition species of multi-acyl classes (e.g. PC 34:1) are NOT used,
    since their sum is not a fatty acid.

    Args:
        lipid_names: List of lipid species names
        include_mono_acyl: When True (default) the single chain of sum-level
            mono-acyl species (source 2) is treated as a candidate FA. Set to
            False to match legacy BioPAN, which does not let a lyso product seed
            the FA pool that in turn validates its own release reaction (e.g. the
            LPS species in PS->LPS); this prevented over-matching of release
            pairs. See ``BioPANPathwayExporter._get_matching_inputs``.

    Returns:
        List of FA names inferred from the dataset

    Example:
        >>> infer_fa_from_lipids(['PC 16:0_18:1', 'PE 18:0_20:4', 'LPI 18:0'])
        ['FA 16:0', 'FA 18:1', 'FA 18:0', 'FA 20:4', 'FA 18:0']
        >>> infer_fa_from_lipids(['PC 16:0_18:1', 'LPI 18:0'], include_mono_acyl=False)
        ['FA 16:0', 'FA 18:1']
    """
    parser = ChainParser()
    seen = set()
    fa_names = []

    def _add(chain: AcylChain) -> None:
        fa_name = f"FA {chain.notation}"
        if fa_name not in seen:
            seen.add(fa_name)
            fa_names.append(fa_name)

    for name in lipid_names:
        species = parser.parse(name)
        if not species or not species.chains:
            continue
        # Full structure: every individual chain is a candidate FA.
        if species.level == StructureLevel.FULL:
            for chain in species.chains:
                _add(chain)
            continue
        # Sum level: only mono-acyl classes expose a usable fatty acyl chain.
        if not include_mono_acyl:
            continue
        base_hg = _strip_linkage_suffix(species.headgroup)
        if base_hg in SINGLE_ACYL_HEADGROUPS:
            _add(species.chains[0])

    return fa_names


def infer_facoa_from_lipids(lipid_names: List[str], include_mono_acyl: bool = True) -> List[str]:
    """Infer possible FACoA species from the acyl chains in a lipid dataset.

    Same as infer_fa_from_lipids but returns CoA (acyl-CoA) names
    in LIPID MAPS shorthand format. See ``infer_fa_from_lipids`` for the
    meaning of ``include_mono_acyl``.
    """
    fa_names = infer_fa_from_lipids(lipid_names, include_mono_acyl=include_mono_acyl)
    return [name.replace("FA ", "CoA ") for name in fa_names]


# =============================================================================
# LM ID-based extraction from LIPID MAPS reactions
# =============================================================================

# LIPID MAPS FA subclass prefixes
# See: https://www.lipidmaps.org/databases/lmsd/overview#FA
FA_LMID_PREFIX = "LMFA0101"       # Fatty acids and conjugates
FACOA_LMID_PREFIX = "LMFA0705"   # Fatty acyl coenzyme A
# CAR_LMID_PREFIX = "LMFA0707"     # Fatty acyl carnitines


def extract_fa_from_reactions(reactions: List["ReactionData"]) -> List[str]:
    """Extract unique FA species names from LIPID MAPS reactions.
    
    Uses compound LM IDs to identify fatty acids (LMFA0101xxxx).
    Returns shorthand names (e.g., 'FA 16:0').
    
    Args:
        reactions: List of ReactionData objects from LIPID MAPS
        
    Returns:
        List of unique FA names in LIPID MAPS shorthand format
        
    Example:
        >>> reactions = [...]  # from LIPID MAPS API
        >>> extract_fa_from_reactions(reactions)
        ['FA 16:0', 'FA 18:1', 'FA 20:4']
    """
    seen = set()
    fa_names = []
    parser = ChainParser()
    
    for reaction in reactions:
        all_compounds = reaction.reactants + reaction.products
        for compound in all_compounds:
            lm_id = compound.compound_lm_id or ""
            if not lm_id.startswith(FA_LMID_PREFIX):
                continue
                
            # Try compound_abbrev first (e.g., "FA 16:0"), then name
            name = compound.compound_abbrev or compound.compound_name
            if not name:
                continue
                
            # Parse to standardize the name
            species = parser.parse(name)
            if species and species.headgroup == "FA" and species.chains:
                standardized = f"FA {species.chains[0].notation}"
                if standardized not in seen:
                    seen.add(standardized)
                    fa_names.append(standardized)
    
    return fa_names


def extract_facoa_from_reactions(reactions: List["ReactionData"]) -> List[str]:
    """Extract unique FACoA/CoA species names from LIPID MAPS reactions.
    
    Uses compound LM IDs to identify:
    - Fatty acyl coenzyme A (LMFA0705xxxx)  
    
    TODO: Fatty acyl carnitines (LMFA0707xxxx) ?
    
    Returns standardized shorthand names (e.g., 'CoA 16:0').
    
    Args:
        reactions: List of ReactionData objects from LIPID MAPS
        
    Returns:
        List of unique CoA names in LIPID MAPS shorthand format
        
    Example:
        >>> reactions = [...]  # from LIPID MAPS API
        >>> extract_facoa_from_reactions(reactions)
        ['CoA 4:0', 'CoA 16:0', 'CoA 18:1']
    """
    seen = set()
    facoa_names = []
    parser = ChainParser()
    
    for reaction in reactions:
        all_compounds = reaction.reactants + reaction.products
        for compound in all_compounds:
            lm_id = compound.compound_lm_id or ""
            # Match both CoA and carnitine subclasses
            if not (lm_id.startswith(FACOA_LMID_PREFIX)):
                continue
                
            # Try compound_abbrev first (e.g., "CoA 4:0"), then name
            name = compound.compound_abbrev or compound.compound_name
            if not name:
                continue
                
            # Parse to standardize the name
            species = parser.parse(name)
            if species and species.headgroup in ("CoA", "FACoA", "FaCoA", "FACOA") and species.chains:
                standardized = f"CoA {species.chains[0].notation}"
                if standardized not in seen:
                    seen.add(standardized)
                    facoa_names.append(standardized)
    
    return facoa_names


def extract_fa_facoa_from_reactions(
    reactions: List["ReactionData"]
) -> Tuple[List[str], List[str]]:
    """Extract both FA and FACoA species from LIPID MAPS reactions.
    
    Convenience function that returns both FA and FACoA names.
    
    Args:
        reactions: List of ReactionData objects from LIPID MAPS
        
    Returns:
        Tuple of (fa_names, facoa_names) in LIPID MAPS shorthand format
    """
    return (
        extract_fa_from_reactions(reactions),
        extract_facoa_from_reactions(reactions)
    )
