from typing import List, Optional, Union, Any, Dict
import logging

from pydantic import BaseModel, Field, field_validator


logger = logging.getLogger(__name__)

""" IN TEMPLATE PHASE"""


class Reaction(BaseModel):
    """
    - reaction_id: identifier for the reaction
    - reaction_level: numeric or string level/score
    - reactants / products: lists of dicts describing lipid species or raw items convertible to dicts
    """

    reaction_id: str
    reaction_name: str
    reaction_level: Union[int, float, str]
    reactants: List[Union[Dict[str, Any], Any]] = Field(default_factory=list)
    products: List[Union[Dict[str, Any], Any]] = Field(default_factory=list)
    type: str  # "species-level" or "class-level"
    pathway_id: Optional[str]
    enzyme_id: Optional[str]

    # Pydantic v2 config
    model_config = {"arbitrary_types_allowed": True}

    @field_validator("reactants", "products", mode="before")
    def _ensure_list(cls, v):
        # Accept None -> empty list
        if v is None:
            return []
        # If single item passed, wrap in list
        if not isinstance(v, list):
            return [v]
        return v

    def __init__(self, **data):
        super().__init__(**data)
        logger.info("Created Reaction: %s", self.reaction_id)

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize reaction to a plain dict. Entries that expose .dict()/.to_dict()
        are converted automatically.
        """

        def _serialize_item(item):
            if hasattr(item, "dict"):
                return item.dict()
            if hasattr(item, "to_dict"):
                return item.to_dict()
            return item

        return {
            "reaction_id": self.reaction_id,
            "reaction_level": self.reaction_level,
            "reactants": [_serialize_item(s) for s in self.reactants],
            "products": [_serialize_item(p) for p in self.products],
        }
