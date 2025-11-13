from typing import Optional, List, Dict, Any, Tuple
from enum import Enum
import re
import logging

from pydantic import BaseModel, Field, model_validator

logger = logging.getLogger(__name__)
""" Early version of lipid species model using Pydantic for validation and parsing."""
""" Will be used to update sample.py """
class LipidNotationType(Enum):
    """Enum for different lipid notation types."""
    SUM_COMPOSITION = "sum"           # PC 36:4
    SPECIES_LEVEL = "species"         # PC 18:2/18:2
    MOLECULAR_SPECIES = "molecular"   # PC 18:2(9Z,12Z)/18:2(9Z,12Z)
    UNKNOWN = "unknown"


class FattyAcidChain(BaseModel):
    """Represents a single fatty acid chain."""
    carbons: int
    double_bonds: int
    stereochemistry: Optional[str] = None  # e.g., "9Z,12Z"

    def __str__(self) -> str:
        base = f"{self.carbons}:{self.double_bonds}"
        return f"{base}({self.stereochemistry})" if self.stereochemistry else base

    def to_dict(self) -> Dict[str, Any]:
        return {
            "carbons": self.carbons,
            "double_bonds": self.double_bonds,
            "stereochemistry": self.stereochemistry
        }


class LipidSpecies(BaseModel):
    """
    Pydantic-based lipid species model.
    - name: original input name
    - lm_id, subclass: optional metadata
    - lipid_class, notation_type, fatty_acid_chains: parsed fields
    """
    name: str
    lm_id: Optional[str] = None
    subclass: Optional[str] = None

    lipid_class: Optional[str] = None
    notation_type: LipidNotationType = LipidNotationType.UNKNOWN
    fatty_acid_chains: List[FattyAcidChain] = Field(default_factory=list)

    # preserve original_name for compatibility
    original_name: Optional[str] = None

    @model_validator(mode="before")
    def parse_name_and_chains(cls, values: Dict[str, Any]) -> Dict[str, Any]:
        raw_name = values.get("name", "") or ""
        name = raw_name.strip()
        values["original_name"] = raw_name
        values["name"] = name

        # extract lipid class prefix e.g., "PC", "PE", "TG"
        m = re.match(r'^([A-Z]{1,4})\s*(.*)$', name)
        if m:
            values["lipid_class"] = m.group(1)
            body = m.group(2).strip()
        else:
            values["lipid_class"] = None
            body = name

        # Determine notation type
        notation = LipidNotationType.UNKNOWN
        if '(' in body and ')' in body:
            notation = LipidNotationType.MOLECULAR_SPECIES
        elif '/' in body or '_' in body:
            notation = LipidNotationType.SPECIES_LEVEL
        elif ':' in body and '/' not in body and '_' not in body:
            notation = LipidNotationType.SUM_COMPOSITION
        values["notation_type"] = notation

        chains: List[FattyAcidChain] = []
        if notation == LipidNotationType.SUM_COMPOSITION:
            mm = re.search(r'(\d+):(\d+)', body)
            if mm:
                chains.append(FattyAcidChain(carbons=int(mm.group(1)), double_bonds=int(mm.group(2))))
        elif notation in (LipidNotationType.SPECIES_LEVEL, LipidNotationType.MOLECULAR_SPECIES):
            parts = [p.strip() for p in re.split(r'\/', body) if p.strip()]
            for part in parts:
                mm = re.match(r'(\d+):(\d+)(?:\(([^)]+)\))?', part)
                if mm:
                    chains.append(FattyAcidChain(
                        carbons=int(mm.group(1)),
                        double_bonds=int(mm.group(2)),
                        stereochemistry=mm.group(3)
                    ))
        values["fatty_acid_chains"] = chains
        return values

    # convenience properties / methods (preserve API)
    @property
    def total_carbons(self) -> Optional[int]:
        if not self.fatty_acid_chains:
            return None
        if self.notation_type == LipidNotationType.SUM_COMPOSITION:
            return self.fatty_acid_chains[0].carbons
        return sum(c.carbons for c in self.fatty_acid_chains)

    @property
    def total_double_bonds(self) -> Optional[int]:
        if not self.fatty_acid_chains:
            return None
        if self.notation_type == LipidNotationType.SUM_COMPOSITION:
            return self.fatty_acid_chains[0].double_bonds
        return sum(c.double_bonds for c in self.fatty_acid_chains)

    def is_sum_composition(self) -> bool:
        return self.notation_type == LipidNotationType.SUM_COMPOSITION

    def is_species_level(self) -> bool:
        return self.notation_type == LipidNotationType.SPECIES_LEVEL

    def is_molecular_species(self) -> bool:
        return self.notation_type == LipidNotationType.MOLECULAR_SPECIES

    def get_sum_composition(self) -> str:
        if self.lipid_class and self.total_carbons is not None and self.total_double_bonds is not None:
            return f"{self.lipid_class} {self.total_carbons}:{self.total_double_bonds}"
        return self.name

    def get_chain_count(self) -> int:
        if self.notation_type == LipidNotationType.SUM_COMPOSITION:
            return 1
        return len(self.fatty_acid_chains)

    def __str__(self) -> str:
        return self.name

    def to_dict(self) -> Dict[str, Any]:
        """Return dict compatible with previous implementation."""
        return {
            "name": self.name,
            "original_name": self.original_name,
            "lm_id": self.lm_id,
            "subclass": self.subclass,
            "lipid_class": self.lipid_class,
            "notation_type": self.notation_type.value,
            "total_carbons": self.total_carbons,
            "total_double_bonds": self.total_double_bonds,
            "fatty_acid_chains": [c.to_dict() for c in self.fatty_acid_chains],
            "sum_composition": self.get_sum_composition(),
            "chain_count": self.get_chain_count()
        }


# Factory function for creating lipid species (keeps previous signature)
def create_lipid_species(name: str, lm_id: Optional[str] = None, subclass: Optional[str] = None) -> LipidSpecies:
    return LipidSpecies(name=name, lm_id=lm_id, subclass=subclass)


# Example usage and quick smoke tests
if __name__ == "__main__":
    test_lipids = [
        "PC 36:4",
        "PC 18:2/18:2",
        "PC 18:2(9Z,12Z)/18:2(9Z,12Z)",
        "TG 54:6",
        "PE 16:0/20:4",
    ]

    for lipid_name in test_lipids:
        lipid = create_lipid_species(lipid_name)
        print(f"\n{lipid}")
        print(f"  Type: {lipid.notation_type.value}")
        print(f"  Class: {lipid.lipid_class}")
        print(f"  Total C:DB: {lipid.total_carbons}:{lipid.total_double_bonds}")
        print(f"  Chains: {len(lipid.fatty_acid_chains)}")