"""
Pydantic models for Cytoscape.js graph data structures.

These models provide type-safe representations of nodes, edges, and complete graphs
that can be serialized to JSON for Cytoscape.js visualization.
"""

from typing import Dict, List, Any, Optional, Union, Literal
from pydantic import BaseModel, Field, field_validator, ConfigDict
from datetime import datetime


class NodePosition(BaseModel):
    """2D position for a node in the graph."""

    x: float = Field(description="X coordinate")
    y: float = Field(description="Y coordinate")


class NodeData(BaseModel):
    """Data attributes for a Cytoscape node.

    This represents the 'data' object within a Cytoscape node element.
    """

    id: str = Field(description="Unique identifier for the node")
    name: Optional[str] = Field(None, description="Display name for the node")
    type: Optional[str] = Field(
        None, description="Node type (e.g., 'lipid', 'reaction', 'pathway')"
    )

    # Lipid-specific fields
    lm_id: Optional[str] = Field(None, description="LIPID MAPS ID")
    standardized_name: Optional[str] = Field(
        None, description="RefMet standardized name"
    )
    formula: Optional[str] = Field(None, description="Chemical formula")
    mass: Optional[Union[float, str]] = Field(None, description="Molecular mass")

    # Reaction-specific fields
    reaction_id: Optional[str] = Field(None, description="Reaction identifier")
    pathway: Optional[str] = Field(None, description="Associated pathway name")
    enzyme: Optional[str] = Field(None, description="Enzyme catalyzing the reaction")

    # Quantification data (for lipids)
    fold_change: Optional[float] = Field(None, description="Fold change value")
    pvalue: Optional[float] = Field(None, description="Statistical p-value")
    significant: Optional[bool] = Field(
        None, description="Whether change is significant"
    )

    # Additional custom attributes
    parent: Optional[str] = Field(
        None, description="Parent node ID for hierarchical graphs"
    )
    label: Optional[str] = Field(None, description="Alternative label for display")

    model_config = ConfigDict(extra="allow")  # Allow additional custom fields


class CytoscapeNode(BaseModel):
    """A Cytoscape node element.

    Represents a single node in the graph with its data and optional position.
    """

    data: NodeData = Field(description="Node data attributes")
    position: Optional[NodePosition] = Field(
        None, description="Node position in 2D space"
    )
    selected: bool = Field(False, description="Whether node is selected")
    selectable: bool = Field(True, description="Whether node can be selected")
    locked: bool = Field(False, description="Whether node position is locked")
    grabbed: bool = Field(False, description="Whether node is currently grabbed")
    grabbable: bool = Field(True, description="Whether node can be grabbed")

    # CSS classes for styling
    classes: Optional[str] = Field(None, description="Space-separated CSS class names")

    model_config = ConfigDict(extra="allow")


class EdgeData(BaseModel):
    """Data attributes for a Cytoscape edge.

    This represents the 'data' object within a Cytoscape edge element.
    """

    id: str = Field(description="Unique identifier for the edge")
    source: str = Field(description="ID of the source node")
    target: str = Field(description="ID of the target node")

    label: Optional[str] = Field(
        None, description="Edge label (e.g., 'substrate', 'product')"
    )
    interaction: Optional[str] = Field(None, description="Type of interaction")
    pathway: Optional[str] = Field(None, description="Pathway this edge belongs to")

    # Reaction-specific fields
    stoichiometry: Optional[float] = Field(
        None, description="Stoichiometric coefficient"
    )
    reversible: Optional[bool] = Field(
        None, description="Whether reaction is reversible"
    )

    # Edge weight/confidence
    weight: Optional[float] = Field(None, description="Edge weight or confidence score")

    model_config = ConfigDict(extra="allow")


class CytoscapeEdge(BaseModel):
    """A Cytoscape edge element.

    Represents a connection between two nodes in the graph.
    """

    data: EdgeData = Field(description="Edge data attributes")
    selected: bool = Field(False, description="Whether edge is selected")
    selectable: bool = Field(True, description="Whether edge can be selected")

    # CSS classes for styling
    classes: Optional[str] = Field(None, description="Space-separated CSS class names")

    model_config = ConfigDict(extra="allow")

    @field_validator("data")
    @classmethod
    def validate_source_target(cls, v: EdgeData) -> EdgeData:
        """Ensure source and target are different."""
        if v.source == v.target:
            raise ValueError("Edge source and target cannot be the same")
        return v


class GraphElements(BaseModel):
    """Collection of nodes and edges in a Cytoscape graph."""

    nodes: List[CytoscapeNode] = Field(
        default_factory=list, description="List of nodes"
    )
    edges: List[CytoscapeEdge] = Field(
        default_factory=list, description="List of edges"
    )

    def add_node(self, node: CytoscapeNode) -> None:
        """Add a node to the graph."""
        # Check for duplicate IDs
        if any(n.data.id == node.data.id for n in self.nodes):
            raise ValueError(f"Node with ID '{node.data.id}' already exists")
        self.nodes.append(node)

    def add_edge(self, edge: CytoscapeEdge) -> None:
        """Add an edge to the graph."""
        # Verify source and target nodes exist
        node_ids = {n.data.id for n in self.nodes}
        if edge.data.source not in node_ids:
            raise ValueError(f"Source node '{edge.data.source}' does not exist")
        if edge.data.target not in node_ids:
            raise ValueError(f"Target node '{edge.data.target}' does not exist")

        self.edges.append(edge)

    def get_node_by_id(self, node_id: str) -> Optional[CytoscapeNode]:
        """Retrieve a node by its ID."""
        return next((n for n in self.nodes if n.data.id == node_id), None)

    def get_edges_for_node(self, node_id: str) -> List[CytoscapeEdge]:
        """Get all edges connected to a node."""
        return [
            e
            for e in self.edges
            if e.data.source == node_id or e.data.target == node_id
        ]


class GraphMetadata(BaseModel):
    """Metadata about the graph."""

    generated_by: str = Field(
        default="BioPAN Python", description="Tool that generated the graph"
    )
    timestamp: datetime = Field(
        default_factory=datetime.now, description="Generation timestamp"
    )
    lipid_count: int = Field(default=0, description="Number of lipid nodes")
    reaction_count: int = Field(default=0, description="Number of reaction nodes")
    sample_count: int = Field(default=0, description="Number of samples in dataset")
    pathway_count: int = Field(default=0, description="Number of unique pathways")
    description: Optional[str] = Field(None, description="Graph description")

    model_config = ConfigDict(extra="allow")


class LayoutOptions(BaseModel):
    """Layout algorithm options for Cytoscape."""

    name: Literal["cose", "circle", "grid", "breadthfirst", "concentric", "random"] = (
        Field("cose", description="Layout algorithm name")
    )
    animate: bool = Field(True, description="Whether to animate the layout")
    animationDuration: int = Field(
        500, description="Animation duration in milliseconds"
    )
    fit: bool = Field(True, description="Fit the graph to viewport")
    padding: int = Field(30, description="Padding around the graph")

    # COSE-specific options
    nodeRepulsion: Optional[int] = Field(
        None, description="Node repulsion force (COSE)"
    )
    idealEdgeLength: Optional[int] = Field(None, description="Ideal edge length (COSE)")
    edgeElasticity: Optional[int] = Field(None, description="Edge elasticity (COSE)")

    # Hierarchical options
    directed: Optional[bool] = Field(
        None, description="Use directed layout (breadthfirst)"
    )
    spacingFactor: Optional[float] = Field(None, description="Spacing between nodes")

    model_config = ConfigDict(extra="allow")


class StyleSelector(BaseModel):
    """CSS-like selector for Cytoscape styling."""

    selector: str = Field(
        description="CSS selector (e.g., 'node', 'edge', '.highlighted')"
    )
    style: Dict[str, Any] = Field(description="Style properties")


class CytoscapeGraph(BaseModel):
    """Complete Cytoscape.js graph structure.

    This is the top-level model that can be directly serialized to JSON
    and loaded by Cytoscape.js.
    """

    elements: GraphElements = Field(description="Graph nodes and edges")
    metadata: Optional[GraphMetadata] = Field(None, description="Optional metadata")
    layout: Optional[LayoutOptions] = Field(None, description="Layout configuration")
    style: Optional[List[StyleSelector]] = Field(
        None, description="Visual style definitions"
    )

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "elements": {
                    "nodes": [
                        {
                            "data": {
                                "id": "LMFA01010001",
                                "name": "Palmitic acid",
                                "type": "lipid",
                                "lm_id": "LMFA01010001",
                                "formula": "C16H32O2",
                            }
                        }
                    ],
                    "edges": [
                        {
                            "data": {
                                "id": "e1",
                                "source": "LMFA01010001",
                                "target": "rxn_123",
                                "label": "substrate",
                            }
                        }
                    ],
                },
                "metadata": {
                    "generated_by": "BioPAN Python",
                    "lipid_count": 1,
                    "reaction_count": 0,
                },
            }
        }
    )

    def to_json(self, **kwargs) -> str:
        """Export graph as JSON string for Cytoscape.js."""
        return self.model_dump_json(exclude_none=True, **kwargs)

    def add_lipid_node(self, lipid_id: str, name: str, **kwargs) -> CytoscapeNode:
        """Helper method to add a lipid node."""
        node_data = NodeData(id=lipid_id, name=name, type="lipid", **kwargs)
        node = CytoscapeNode(data=node_data)
        self.elements.add_node(node)
        return node

    def add_reaction_node(self, reaction_id: str, name: str, **kwargs) -> CytoscapeNode:
        """Helper method to add a reaction node."""
        node_data = NodeData(
            id=f"rxn_{reaction_id}",
            name=name,
            type="reaction",
            reaction_id=reaction_id,
            **kwargs,
        )
        node = CytoscapeNode(data=node_data)
        self.elements.add_node(node)
        return node

    def add_interaction_edge(
        self, source: str, target: str, interaction_type: str, **kwargs
    ) -> CytoscapeEdge:
        """Helper method to add an interaction edge."""
        edge_data = EdgeData(
            id=f"edge_{source}_to_{target}",
            source=source,
            target=target,
            label=interaction_type,
            interaction=interaction_type,
            **kwargs,
        )
        edge = CytoscapeEdge(data=edge_data)
        self.elements.add_edge(edge)
        return edge

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CytoscapeGraph":
        """Create graph from dictionary (e.g., loaded from JSON)."""
        return cls.model_validate(data)


class GraphBuilder(BaseModel):
    """Helper class for building Cytoscape graphs programmatically."""

    graph: CytoscapeGraph = Field(
        default_factory=lambda: CytoscapeGraph(elements=GraphElements())
    )

    def build_lipid_reaction_network(
        self, lipids: List[Dict[str, Any]], reactions: List[Dict[str, Any]]
    ) -> CytoscapeGraph:
        """Build a network from lipids and reactions.

        Args:
            lipids: List of lipid dictionaries with 'lm_id', 'name', etc.
            reactions: List of reaction dictionaries with 'reaction_id', 'reactants', 'products'

        Returns:
            CytoscapeGraph instance
        """
        # Add lipid nodes
        for lipid in lipids:
            self.graph.add_lipid_node(
                lipid_id=lipid.get("lm_id", lipid.get("id")),
                name=lipid.get("name", "Unknown"),
                lm_id=lipid.get("lm_id"),
                standardized_name=lipid.get("standardized_name"),
                formula=lipid.get("formula"),
                mass=lipid.get("mass"),
                fold_change=lipid.get("fold_change"),
                pvalue=lipid.get("pvalue"),
            )

        # Add reaction nodes and edges
        for reaction in reactions:
            reaction_id = reaction.get("reaction_id") or reaction.get("id")
            if not reaction_id:
                continue

            # Add reaction node
            self.graph.add_reaction_node(
                reaction_id=reaction_id,
                name=reaction.get("name", f"Reaction {reaction_id}"),
                pathway=reaction.get("pathway"),
                enzyme=reaction.get("enzyme"),
            )

            # Add edges from reactants to reaction
            for reactant in reaction.get("reactants", []):
                lm_id = reactant.get("compound_lm_id")
                if lm_id:
                    try:
                        self.graph.add_interaction_edge(
                            source=lm_id,
                            target=f"rxn_{reaction_id}",
                            interaction_type="substrate",
                            stoichiometry=reactant.get("stoichiometry"),
                        )
                    except ValueError:
                        # Node doesn't exist, skip
                        pass

            # Add edges from reaction to products
            for product in reaction.get("products", []):
                lm_id = product.get("compound_lm_id")
                if lm_id:
                    try:
                        self.graph.add_interaction_edge(
                            source=f"rxn_{reaction_id}",
                            target=lm_id,
                            interaction_type="product",
                            stoichiometry=product.get("stoichiometry"),
                        )
                    except ValueError:
                        # Node doesn't exist, skip
                        pass

        # Update metadata
        self.graph.metadata = GraphMetadata(
            lipid_count=len(
                [n for n in self.graph.elements.nodes if n.data.type == "lipid"]
            ),
            reaction_count=len(
                [n for n in self.graph.elements.nodes if n.data.type == "reaction"]
            ),
            pathway_count=len(
                set(n.data.pathway for n in self.graph.elements.nodes if n.data.pathway)
            ),
        )

        return self.graph
