import os
import json
import pytest
from biopan.src.models.cytoscape import CytoscapeGraph, GraphBuilder

@pytest.fixture
def graph_builder():
    return GraphBuilder()

def test_cytoscape_graph_creation_and_export(graph_builder, tmp_path):
    # Add nodes
    graph_builder.graph.add_lipid_node(
        lipid_id="LMFA01010001",
        name="Palmitic acid",
        formula="C16H32O2",
        mass=256.24
    )
    graph_builder.graph.add_reaction_node(
        reaction_id="123",
        name="Fatty acid synthesis",
        pathway="Lipid metabolism"
    )
    # Add edge
    graph_builder.graph.add_interaction_edge(
        source="LMFA01010001",
        target="rxn_123",
        interaction_type="substrate"
    )

    # Export to temporary file
    graph_json_path = tmp_path / "graph.json"
    with open(graph_json_path, 'w') as f:
        f.write(graph_builder.graph.to_json(indent=2))

    # Load and validate graph
    with open(graph_json_path, 'r') as f:
        data = json.load(f)
    graph = CytoscapeGraph.from_dict(data)
    assert len(graph.elements.nodes) > 0
    assert any(node.data.name == "Palmitic acid" for node in graph.elements.nodes)
    assert any(node.data.name == "Fatty acid synthesis" for node in graph.elements.nodes)

def test_cytoscape_graph_loading(tmp_path):
    # Prepare a sample graph file
    graph_data = {
        "elements": {
            "nodes": [
                {"data": {"id": "LMFA01010001", "lm_id": "LMFA01010001", "name": "Palmitic acid"}},
                {"data": {"id": "rxn_123", "reaction_id": "123", "name": "Fatty acid synthesis"}}
            ],
            "edges": [
                {"data": {"id": "edge_1", "source": "LMFA01010001", "target": "rxn_123", "interaction_type": "substrate"}}
            ]
        },
                        "metadata": {
                    "generated_by": "BioPAN Python",
                    "lipid_count": 1,
                    "reaction_count": 0
                }
    }
    graph_json_path = tmp_path / "lipid_graph.json"
    with open(graph_json_path, 'w') as f:
        json.dump(graph_data, f)

    # Load and validate existing graph
    with open(graph_json_path, 'r') as f:
        data = json.load(f)
    graph = CytoscapeGraph.from_dict(data)
    assert len(graph.elements.nodes) == 2
    assert graph.elements.nodes[0].data.name == "Palmitic acid"
    assert graph.elements.nodes[1].data.name == "Fatty acid synthesis"
    assert graph.to_json(indent=2) is not None
    # Prepare a sample graph file
    graph_data = {
        "elements": {
            "nodes": [
                {"data": {"id": "LMFA01010001","lm_id": "LMFA01010001", "name": "Palmitic acid"}},
                {"data": {"id": "rxn_123", "reaction_id": "123", "name": "Fatty acid synthesis"}}
            ],
            "edges": [
                {"data": {"id": "edge_1", "source": "LMFA01010001", "target": "rxn_123", "interaction_type": "substrate"}}
            ]
        }
    }
    graph_json_path = tmp_path / "lipid_graph.json"
    with open(graph_json_path, 'w') as f:
        json.dump(graph_data, f)

    # Load and validate existing graph
    with open(graph_json_path, 'r') as f:
        data = json.load(f)
    graph = CytoscapeGraph.from_dict(data)
    assert len(graph.elements.nodes) == 2
    assert graph.elements.nodes[0].data.name == "Palmitic acid"
    assert graph.elements.nodes[1].data.name == "Fatty acid synthesis"
    assert graph.to_json(indent=2) is not None