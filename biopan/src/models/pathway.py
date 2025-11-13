from typing import List, Optional, Union, Any, Dict
from .lipid_species import LipidSpecies
import logging

from pydantic import BaseModel, Field, field_validator
""" IN TEMPLATE PHASE"""

logger = logging.getLogger(__name__)

class Pathway(BaseModel):
    name: str
    reactions: List[Any]  # Replace Any with the actual Reaction type when defined

    @field_validator('name')
    def name_must_not_be_empty(cls, v):
        if not v:
            raise ValueError('Pathway name must not be empty')
        return v

    def __init__(self, **data):
        super().__init__(**data)
        logger.info(f"Created Pathway: {self.name}")
