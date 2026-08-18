
import logging
import requests
from typing import List, Optional, Union, Dict, Any

from pydantic import ConfigDict, Field, field_validator
from .base import LipidmapsBaseModel
from ..utils.shunt import reaction_is_shunt
logger = logging.getLogger(__name__)


class CompoundComponent(LipidmapsBaseModel):
    """Represents a compound in a reaction (reactant or product)."""

    compound_type: Optional[str] = None
    compound_name: Optional[str] = None
    compound_lm_id: Optional[str] = None
    compound_sys_name: Optional[str] = None
    compound_synonyms: Optional[str] = None
    compound_generic_lm_id: Optional[str] = None
    compound_abbrev: Optional[str] = None
    compound_abbrev_chains: Optional[str] = None
    compound_headgroup: Optional[str] = None
    compound_full_struct: Optional[str] = None
    compound_smiles: Optional[str] = None
    def display_name(self) -> str:
        """Get the best available name for this compound."""
        return (
            self.compound_name
            or self.compound_lm_id
            or self.compound_generic_lm_id
            or "Unknown"
        )


class ReactionData(LipidmapsBaseModel):
    """Represents a single reaction with reactants and products."""

    reactants: List[CompoundComponent] = Field(default_factory=list)
    products: List[CompoundComponent] = Field(default_factory=list)
    genes: List[Dict[str, Any]] = Field(default_factory=list)
    proteins: List[Dict[str, Any]] = Field(default_factory=list)
    curations: List[Dict[str, Any]] = Field(default_factory=list)
    pathways: List[Dict[str, Any]] = Field(default_factory=list)
    rhea_id: Optional[Union[str, List[str]]] = None
    reaction_name: Optional[str] = None
    reaction_id: Optional[int] = None
    reaction_type: Optional[str] = None
    # True when this reaction belongs to a shunt / alternate-route pathway.
    # Derived from ``pathways`` in ``model_post_init`` (or honoured if the API
    # sends it directly). See ``lipidmaps.data.utils.shunt``.
    is_shunt: bool = False
    # Detailed evaluation produced by `ReactionEvaluator.evaluate_reaction()`.
    # Contains keys like `possible` (bool), `explanation` (str), and any
    # additional diagnostic `details` the evaluator may include.
    evaluation: Optional[Dict[str, Any]] = None

    @staticmethod
    def _coerce_rhea_id(value: Any) -> Optional[str]:
        if value is None:
            return None
        text = str(value).strip()
        if not text:
            return None
        if text.lower().startswith("rhea:"):
            text = text.split(":", 1)[1].strip()
        return text or None

    @classmethod
    def _extract_rhea_ids_from_curations(cls, curations: List[Dict[str, Any]]) -> List[str]:
        ordered: List[str] = []
        seen = set()

        for curation in curations or []:
            database_name = str(curation.get("database_name", "") or "").strip().lower()
            if database_name != "rhea":
                continue

            rhea_id = cls._coerce_rhea_id(curation.get("database_id"))
            if rhea_id and rhea_id not in seen:
                ordered.append(rhea_id)
                seen.add(rhea_id)

        return ordered

    @field_validator("rhea_id", mode="before")
    @classmethod
    def _normalize_rhea_id_field(cls, value: Any) -> Optional[Union[str, List[str]]]:
        if value is None:
            return None

        values = value if isinstance(value, list) else [value]
        ordered: List[str] = []
        seen = set()
        for item in values:
            rhea_id = cls._coerce_rhea_id(item)
            if rhea_id and rhea_id not in seen:
                ordered.append(rhea_id)
                seen.add(rhea_id)

        if not ordered:
            return None
        if len(ordered) == 1:
            return ordered[0]
        return ordered

    def list_rhea_ids(self) -> List[str]:
        """Return normalized Rhea identifiers found on this reaction."""

        ordered: List[str] = []
        seen = set()

        values = self.rhea_id if isinstance(self.rhea_id, list) else [self.rhea_id]
        for value in values:
            rhea_id = self._coerce_rhea_id(value)
            if rhea_id and rhea_id not in seen:
                ordered.append(rhea_id)
                seen.add(rhea_id)
        return ordered

    def model_post_init(self, __context: Any) -> None:
        """Populate normalized Rhea identifiers and the shunt flag."""

        # Derive the shunt flag from pathway metadata (honouring any structured
        # value the API already provided).
        self.is_shunt = self.is_shunt or reaction_is_shunt(self.pathways)

        if self.list_rhea_ids():
            return

        curated_rhea_ids = self._extract_rhea_ids_from_curations(self.curations)
        if not curated_rhea_ids:
            return

        if len(curated_rhea_ids) == 1:
            self.rhea_id = curated_rhea_ids[0]
        else:
            self.rhea_id = curated_rhea_ids

    @property
    def organisms(self) -> List[str]:
        """Return list of organisms where proteins in this reaction are present."""
        organism_fields = ['bacteria', 'archaea', 'fungi', 'viridiplantae', 
                          'mammalia', 'arthropoda', 'eukaryota']
        organisms_set = set()
        
        # Iterate through all proteins and collect organisms
        for protein in self.proteins:
            for org in organism_fields:
                if protein.get(org, False):
                    organisms_set.add(org)
        
        return sorted(list(organisms_set))
    
    # Allow additional fields from API response
    model_config = ConfigDict(extra="allow")

    @staticmethod
    def _component_identifier(comp: CompoundComponent) -> Optional[str]:
        return comp.compound_lm_id or comp.compound_generic_lm_id

    def list_generic_lm_ids(self) -> List[str]:
        """List all generic_lm_ids from reactants and products."""
        ids = set()
        for comp in self.reactants + self.products:
            if comp.compound_generic_lm_id or comp.compound_lm_id:
                ids.add(comp.compound_generic_lm_id or comp.compound_lm_id)
        return sorted(ids)
    
    def list_lm_ids(self) -> List[str]:
        """List all lm_ids from reactants and products."""
        ids = set()
        for comp in self.reactants + self.products:
            if comp.compound_lm_id:
                ids.add(comp.compound_lm_id)
        return sorted(ids)
    
    def list_nonfa_noncoa_reactant_lm_ids(self) -> List[str]:
        """List all lm_ids from reactants."""
        ids = set()
        for component in self.reactants:
            component_id = self._component_identifier(component)
            if component_id and not (component_id == "LMFA01010000" or component_id.startswith("LMFA0705")):
                ids.add(component_id)
        return sorted(ids)
    
    def list_nonfa_noncoa_product_lm_ids(self) -> List[str]:
        """List all lm_ids from products."""
        ids = set()
        for component in self.products:
            component_id = self._component_identifier(component)
            if component_id and not (component_id == "LMFA01010000" or component_id.startswith("LMFA0705")):
                ids.add(component_id)
        return sorted(ids)

    def has_lm_main_components(self) -> bool:
        """Check if reaction has any lm_main compounds."""
        all_components = self.reactants + self.products
        return any(c.compound_type == "lm_main" for c in all_components)

    def filter_reaction(self, only_lipid_components: bool = True) -> "ReactionData":
        """Return a new ReactionData with only lm_main components."""
        if only_lipid_components:
            filtered_reactants = [c for c in self.reactants if c.compound_type == "lm_main"]
            filtered_products = [c for c in self.products if c.compound_type == "lm_main"]
        else:
            filtered_reactants = self.reactants
            filtered_products = self.products
        # Generate reaction name from filtered components by using lm_main names
        reactant_names = "; ".join(c.display_name() for c in filtered_reactants if c.compound_type == "lm_main")
        product_names = "; ".join(c.display_name() for c in filtered_products if c.compound_type == "lm_main")
        reaction_name = (
            f"{reactant_names} -> {product_names}"
            if (reactant_names or product_names)
            else None
        )

        return ReactionData(
            reaction_id=self.reaction_id,
            reaction_name=reaction_name,
            reactants=filtered_reactants,
            products=filtered_products,
            genes=self.genes,
            proteins=self.proteins,
            curations=self.curations,
            rhea_id=self.rhea_id,
            pathways=self.pathways,
            reaction_type=self.reaction_type,
            is_shunt=self.is_shunt,
        )


class ReactionResponse(LipidmapsBaseModel):
    """Response from reaction checking API."""

    reactions: List[ReactionData] = Field(default_factory=list)
    error: Optional[str] = None


class ReactionChecker(LipidmapsBaseModel):
    """Handle reaction checking via LIPID MAPS API with Pydantic validation."""

    base_url: str = Field(..., description="Base URL for the reaction API")
    endpoint: str = Field(default="/api/reactions", description="API endpoint path")
    timeout: int = Field(
        default=120, description="Request timeout in seconds", ge=1, le=600
    )

    # Computed field for full API URL
    @property
    def api_url(self) -> str:
        """Get the full API URL."""
        return self.base_url.rstrip("/") + self.endpoint

    def model_post_init(self, __context: Any) -> None:
        """Initialize after model creation."""
        logger.info(f"Initialized ReactionChecker with URL: {self.api_url}")

    def check_reactions(
        self,
        lm_ids: Union[List[str], str],
        search_type: str = "lipids",
        search_mode: str = "default",
        generic_reactions: bool = True,
        only_lipid_components: bool = True,
        taxonomy_group: Optional[str] = "all",
        reaction_type: Optional[str] = None,
    ) -> ReactionResponse:
        """Check reactions for given LIPID MAPS IDs.

        Builds a request payload compatible with the LIPID MAPS reaction API.

        Args:
            lm_ids: List of LIPID MAPS IDs to check
            search_type: Type of search (e.g. "lipids")
            search_mode: Mode string for the API (default: "default")
            generic_reactions: Whether to request generic reactions
            only_lipid_components: Whether to include only lipid components in the response
            taxonomy_group: Optional taxonomy group to filter reactions
        Returns:
            ReactionResponse with filtered reactions containing only lm_main components if specified
        """
        # Build payload matching the API used in the HTTP example (keys like "lmsd_ids").
        # Include `search_type` only when requesting lipid-only components (backwards-compatible).
        payload: Dict[str, Any] = {
            "search_mode": search_mode,
            "generic_reactions": generic_reactions,
            "lmsd_ids": lm_ids,
            "search_source": "lipidmaps_py",
        }
        if only_lipid_components:
            payload["search_type"] = search_type
        if taxonomy_group != "all":
            payload["taxonomy_group"] = taxonomy_group

        logger.debug("Reaction check payload: %s", payload)
        try:
            logger.info(f"Sending reaction check request for {len(lm_ids)} LM IDs to {self.api_url}")
            response = requests.post(
                self.api_url,
                json=payload,
                timeout=self.timeout,
            )
            try:
                response.raise_for_status()
            except requests.HTTPError:
                # Log response body where available to help debugging 400/500 errors
                try:
                    body = response.text
                except Exception:
                    body = "<unavailable>"
                logger.error(
                    "Reaction API returned HTTP %s: %s",
                    response.status_code,
                    body,
                )
                raise

            raw_data = response.json()

            # Parse raw data into ReactionData objects
            reactions = []
            for raw_reaction in raw_data:
                try:
                    # Parse reactants and products
                    reactants = [
                        CompoundComponent(**comp)
                        for comp in raw_reaction.get("reactants", [])
                    ]
                    products = [
                        CompoundComponent(**comp)
                        for comp in raw_reaction.get("products", [])
                    ]

                    reaction = ReactionData(
                        reactants=reactants,
                        products=products,
                        **{
                            k: v
                            for k, v in raw_reaction.items()
                            if k not in ["reactants", "products"]
                        },
                        reaction_type=reaction_type,
                    )

                    # Filter to keep only lm_main components
                    filtered_reaction = reaction.filter_reaction(only_lipid_components=only_lipid_components)

                    # Only include if it has lm_main components
                    if filtered_reaction.reactants or filtered_reaction.products:
                        reactions.append(filtered_reaction)

                except Exception as e:
                    logger.warning(f"Failed to parse reaction: {e}")
                    continue

            logger.info(f"Successfully retrieved {len(reactions)} reactions")
            return ReactionResponse(reactions=reactions)

        except requests.RequestException as e:
            logger.error("Reaction API call failed: %s", e)
            return ReactionResponse(reactions=[], error=str(e))
