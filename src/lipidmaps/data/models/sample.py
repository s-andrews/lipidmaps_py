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


if __name__ == "__main__":

    lipid = QuantifiedLipid(
        input_name="PC(16:0/18:1)",
        values={"Sample1": 10.2, "Sample2": 11.3, "Sample3": 9.8},
        pathway_ids=["R-HSA-1483257"],
        pathway_names=["Glycerophospholipid metabolism"],
        enzyme_ids=["EC 2.3.1.51"],
    )
    zscores = lipid.zscore()
    print("Z-scores:", zscores, lipid)
