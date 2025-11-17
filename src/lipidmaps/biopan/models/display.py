import logging
from typing import Optional

from pydantic import BaseModel


logger = logging.getLogger(__name__)


class CytoscapeOutput(BaseModel):
    nodes: list[dict]
    edges: list[dict]


class CytoscapeNode(BaseModel):
    data: dict
    position: dict
    selected: bool = False
    locked: bool = False


class CytoscapeEdge(BaseModel):
    data: dict
    source: str
    target: str
    selected: bool = False
    locked: bool = False


class CytoscapeNodeData(BaseModel):
    id: str
    name: str
    type: Optional[str] = None  # compound / reaction / pathway etc
    lm_id: Optional[str] = None
    formula: Optional[str] = None
    mass: Optional[float] = None


class CytoscapeEdgeData(BaseModel):
    id: str
    source: str  # source node id
    target: str  # target node id
    direction: Optional[str] = None  # direction of the edge
    pathway: Optional[str] = None  # pathway id
