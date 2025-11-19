import logging
import requests
from typing import List, Dict, Any, Optional

from pydantic import BaseModel, Field, HttpUrl, field_validator

logger = logging.getLogger(__name__)


class CompoundComponent(BaseModel):
    """Represents a compound in a reaction (reactant or product)."""
    compound_type: Optional[str] = None
    compound_name: Optional[str] = None
    compound_lm_id: Optional[str] = None
    compound_generic_id: Optional[str] = None
    
    def display_name(self) -> str:
        """Get the best available name for this compound."""
        return (
            self.compound_name 
            or self.compound_lm_id 
            or self.compound_generic_id 
            or "Unknown"
        )


class ReactionData(BaseModel):
    """Represents a single reaction with reactants and products."""
    reactants: List[CompoundComponent] = Field(default_factory=list)
    products: List[CompoundComponent] = Field(default_factory=list)
    reaction_name: Optional[str] = None
    # Allow additional fields from API response
    model_config = {"extra": "allow"}
    
    def has_lm_main_components(self) -> bool:
        """Check if reaction has any lm_main compounds."""
        all_components = self.reactants + self.products
        return any(c.compound_type == "lm_main" for c in all_components)
    
    def filter_lm_main(self) -> "ReactionData":
        """Return a new ReactionData with only lm_main components."""
        filtered_reactants = [
            c for c in self.reactants if c.compound_type == "lm_main"
        ]
        filtered_products = [
            c for c in self.products if c.compound_type == "lm_main"
        ]
        
        # Generate reaction name from filtered components
        reactant_names = "; ".join(c.display_name() for c in filtered_reactants)
        product_names = "; ".join(c.display_name() for c in filtered_products)
        reaction_name = f"{reactant_names} -> {product_names}" if (reactant_names or product_names) else None
        
        return ReactionData(
            reactants=filtered_reactants,
            products=filtered_products,
            reaction_name=reaction_name
        )


class ReactionResponse(BaseModel):
    """Response from reaction checking API."""
    reactions: List[ReactionData] = Field(default_factory=list)
    error: Optional[str] = None


class ReactionChecker(BaseModel):
    """Handle reaction checking via LIPID MAPS API with Pydantic validation."""
    
    base_url: str = Field(..., description="Base URL for the reaction API")
    endpoint: str = Field(default="/api/reactions", description="API endpoint path")
    timeout: int = Field(default=10, description="Request timeout in seconds", ge=1, le=300)
    
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
        lm_ids: List[str], 
        search_type: str = "lm_id"
    ) -> ReactionResponse:
        """Check reactions for given LIPID MAPS IDs.
        
        Args:
            lm_ids: List of LIPID MAPS IDs to check
            search_type: Type of search (default: "lm_id")
            
        Returns:
            ReactionResponse with filtered reactions containing only lm_main components
        """
        payload = {
            "search_source": "biopan",
            "search_type": search_type,
            "lm_ids": lm_ids,
        }
        
        try:
            logger.info(f"Sending reaction check request for {len(lm_ids)} LM IDs")
            response = requests.post(
                self.api_url, 
                json=payload, 
                timeout=self.timeout
            )
            response.raise_for_status()
            
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
                        **{k: v for k, v in raw_reaction.items() 
                           if k not in ["reactants", "products"]}
                    )
                    
                    # Filter to keep only lm_main components
                    filtered_reaction = reaction.filter_lm_main()
                    
                    # Only include if it has lm_main components
                    if filtered_reaction.reactants or filtered_reaction.products:
                        reactions.append(filtered_reaction)
                        
                except Exception as e:
                    logger.warning(f"Failed to parse reaction: {e}")
                    continue
            
            logger.info(f"Successfully retrieved {len(reactions)} reactions")
            return ReactionResponse(reactions=reactions)
            
        except requests.RequestException as e:
            logger.error(f"Reaction API call failed: {e}")
            return ReactionResponse(reactions=[], error=str(e))
