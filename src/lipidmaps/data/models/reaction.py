from typing import List, Optional, Union, Any, Dict
import logging

from pydantic import BaseModel, Field, field_validator


logger = logging.getLogger(__name__)

""" IN TEMPLATE PHASE"""


class Reaction(BaseModel):
    """
    - reaction_id: identifier for the reaction (int or str)
    - reaction_name: optional human-readable name
    - reaction_level: numeric or string level/score
    - proteins: list of {ec_number, description} dicts
    - genes: list of gene identifiers or dicts
    - curations: list of curation dicts (database_name, database_id, citation_id, info)
    - reactants / products: lists of dicts describing lipid species (keeps LMSD/Reaction API shape)
    """

    reaction_id: Union[int, str]
    reaction_name: Optional[str] = None
    reaction_level: Optional[Union[int, float, str]] = None

    # biological context / provenance
    proteins: List[Dict[str, Any]] = Field(default_factory=list)
    genes: List[Union[str, Dict[str, Any]]] = Field(default_factory=list)
    curations: List[Dict[str, Any]] = Field(default_factory=list)

    # reaction participants: keep as flexible dicts to match API payloads
    reactants: List[Dict[str, Any]] = Field(default_factory=list)
    products: List[Dict[str, Any]] = Field(default_factory=list)

    # metadata
    type: Optional[str] = None  # "species-level" or "class-level"
    pathway_id: Optional[str] = None
    enzyme_id: Optional[str] = None

    # Pydantic v2 config
    model_config = {"arbitrary_types_allowed": True}

    @field_validator("proteins", "genes", "curations", "reactants", "products", mode="before")
    def _ensure_list(cls, v):
        # Accept None -> empty list; wrap single item into list
        if v is None:
            return []
        if isinstance(v, (dict, str)):
            return [v]
        if not isinstance(v, list):
            return list(v)
        return v

    def __init__(self, **data):
        super().__init__(**data)
        logger.info(f"Created Reaction: {self.reaction_id}")

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize reaction to a plain dict. Preserve proteins, genes, curations,
        reactants and products as provided by upstream APIs.
        """

        def _serialize_item(item):
            if hasattr(item, "dict"):
                return item.dict()
            if hasattr(item, "to_dict"):
                return item.to_dict()
            return item

        return {
            "reaction_id": self.reaction_id,
            "reaction_name": self.reaction_name,
            "reaction_level": self.reaction_level,
            "proteins": [_serialize_item(p) for p in self.proteins],
            "genes": [_serialize_item(g) for g in self.genes],
            "curations": [_serialize_item(c) for c in self.curations],
            "reactants": [_serialize_item(r) for r in self.reactants],
            "products": [_serialize_item(p) for p in self.products],
            "type": self.type,
            "pathway_id": self.pathway_id,
            "enzyme_id": self.enzyme_id,
        }