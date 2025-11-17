import logging
from typing import Any, List, Dict, Optional
from pydantic import BaseModel
import numpy as np

logger = logging.getLogger(__name__)


class SampleMetadata(BaseModel):
    sample_id: str
    group: str  # e.g., "Control", "WT"
    condition: Optional[str] = None  # e.g., "Fasted", "Fed"


class QuantifiedLipid(BaseModel):
    input_name: str
    values: Dict[str, float]  # sample_id -> value
    pathway_ids: Optional[List[str]] = None  # e.g., KEGG or Reactome IDs
    pathway_names: Optional[List[str]] = None  # Human-readable names
    enzyme_ids: Optional[List[str]] = None  # e.g., EC numbers or UniProt IDs
    # RefMet annotations
    standardized_name: Optional[str] = None
    lm_id: Optional[str] = None
    generic_lm_id: Optional[str] = None
    sub_class: Optional[str] = None
    super_class: Optional[str] = None
    main_class: Optional[str] = None
    chebi_id: Optional[str] = None
    kegg_id: Optional[str] = None
    refmet_id: Optional[str] = None
    formula: Optional[str] = None
    mass: Optional[float] = None
    reactions: Optional[List[Dict[str, Any]]] = (
        None  # e.g., {"reactant": "PC", "product": "LPC", "type": "class-level"}
    )
    weight: Optional[float] = None  # For species or class-level reaction

    def zscore(self) -> Dict[str, float]:
        vals = np.array(list(self.values.values()))
        mean = np.mean(vals)
        std = np.std(vals)
        return {
            k: (v - mean) / std if std != 0 else 0.0 for k, v in self.values.items()
        }


class LipidDataset(BaseModel):
    samples: List[SampleMetadata]
    lipids: List[QuantifiedLipid]

    def get_grouped_data(self) -> Dict[str, List[QuantifiedLipid]]:
        grouped = {}
        for sample in self.samples:
            grouped.setdefault(sample.group, []).append(sample.sample_id)
        result = {}
        for group, sample_ids in grouped.items():
            result[group] = [
                QuantifiedLipid(
                    input_name=lipid.input_name,
                    values={
                        sid: lipid.values[sid]
                        for sid in sample_ids
                        if sid in lipid.values
                    },
                )
                for lipid in self.lipids
            ]
        return result

    # def normalize(self) -> None:
    #     # Normalize concentrations across samples
    #     totals = {sid: sum(lipid.values.get(sid, 0) for lipid in self.lipids) for sid in self.sample_ids()}
    #     for lipid in self.lipids:
    #         lipid.values = {sid: val / totals[sid] if totals[sid] != 0 else 0 for sid, val in lipid.values.items()}

    # def sample_ids(self) -> List[str]:
    #     return [s.sample_id for s in self.samples]

    # def compute_weights(self, reactions: List[Dict[str, str]]) -> None:
    #     # reactions: [{"reactant": "PC", "product": "LPC", "type": "class-level"}]
    #     for r in reactions:
    #         reactant_sum = sum(l.values.get(sid, 0) for l in self.lipids if r["reactant"] in l.input_name for sid in self.sample_ids())
    #         product_sum = sum(l.values.get(sid, 0) for l in self.lipids if r["product"] in l.input_name for sid in self.sample_ids())
    #         weight = product_sum / reactant_sum if reactant_sum != 0 else None
    #         # Store weight at dataset or reaction level


class Sample(BaseModel):
    sample_name: str
    conditions: Dict[str, Any]
    level: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
    refmet_result: Optional[Dict[str, Any]] = None
    processed: bool = False

    def model_post_init(self, __context: dict) -> None:
        logger.info("Initialized Sample: %s", self.sample_name)

    @property
    def recognized(self):
        return bool(self.refmet_result and self.refmet_result.get("lm_id"))

    def to_dict(self):
        return {
            "sample_name": self.sample_name,
            "conditions": self.conditions,
            "level": self.level,
            "recognized": self.recognized,
            "metadata": self.metadata,
            "refmet_result": self.refmet_result,
        }


if __name__ == "__main__":
    # Example usage
    sample = Sample(
        sample_name="Sample1",
        conditions={"group": "Control", "condition": "Fasted"},
        level="species",
        metadata={"experiment_date": "2024-06-01"},
    )
    print(sample.to_dict())

    lipid = QuantifiedLipid(
        input_name="PC(16:0/18:1)",
        values={"Sample1": 10.2, "Sample2": 11.3, "Sample3": 9.8},
        pathway_ids=["R-HSA-1483257"],
        pathway_names=["Glycerophospholipid metabolism"],
        enzyme_ids=["EC 2.3.1.51"],
    )
    zscores = lipid.zscore()
    print("Z-scores:", zscores, lipid)
