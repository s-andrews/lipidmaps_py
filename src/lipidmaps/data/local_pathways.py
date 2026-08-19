"""Local curated pathway definitions for BioPAN exports."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from typing import Any, Dict, List, Sequence

from pydantic import Field, computed_field

from .models.base import LipidmapsBaseModel


class CuratedPathwayStep(LipidmapsBaseModel):
    """A single ordered class-to-class step within a curated pathway."""

    reactant_class: str = Field(...)
    product_class: str = Field(...)

    @computed_field
    @property
    def reaction_key(self) -> str:
        return f"{self.reactant_class},{self.product_class}".lower()


class CuratedPathwayDefinition(LipidmapsBaseModel):
    """Ordered curated pathway definition used for BioPAN pathway export."""

    pathway_id: str = Field(...)
    pathway_name: str = Field(...)
    pathway_type: List[str] = Field(default_factory=list)
    description: str = Field(default="")
    steps: List[CuratedPathwayStep] = Field(default_factory=list)

    @computed_field
    @property
    def step_keys(self) -> List[str]:
        return [step.reaction_key for step in self.steps]


class CuratedPathwayCatalog(LipidmapsBaseModel):
    """Container for packaged curated pathway definitions."""

    pathways: List[CuratedPathwayDefinition] = Field(default_factory=list)


@lru_cache(maxsize=1)
def load_curated_biopan_pathways() -> List[CuratedPathwayDefinition]:
    """Load packaged curated pathway definitions for BioPAN exports."""

    config_dir = resources.files("lipidmaps.data.config")
    with config_dir.joinpath("biopan_pathways.json").open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    catalog = CuratedPathwayCatalog.model_validate(payload)
    return catalog.pathways


def build_pathway_step_lookup(
    pathways: Sequence[CuratedPathwayDefinition],
) -> Dict[str, List[Dict[str, Any]]]:
    """Map reaction keys to pathway metadata for edges participating in curated pathways."""

    lookup: Dict[str, List[Dict[str, Any]]] = {}
    for pathway in pathways:
        pathway_entry = {
            "pathway_name": pathway.pathway_name,
            "pathway_type": list(pathway.pathway_type),
        }
        for step in pathway.steps:
            entries = lookup.setdefault(step.reaction_key, [])
            if pathway_entry not in entries:
                entries.append(pathway_entry)
    return lookup