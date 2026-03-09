
import logging
from typing import Any, List, Dict, Optional, Union, Callable
import numpy as np
import re
from ..utils.headgroups import lipidmaps_headgroups
from .query import Query, from_callable, attr_eq, attr_in, attr_contains, attr_gt, has_attr
from .base import LipidmapsBaseModel
from pydantic import Field
from .reaction import ReactionData, CompoundComponent, ReactionChecker
from ...config import LMSD_REACTIONS_BASE_URL
from ..validation.data_validator import DataValidator, ValidationReport


logger = logging.getLogger(__name__)


class SampleMetadata(LipidmapsBaseModel):
    sample_name: str
    group: str  # e.g., "Control", "WT"
    label: Optional[str] = None  # e.g., "Fasted", "Fed"
    values: Optional[Dict[str, float]] = None  # lipid input_name -> value (optional cache)

    def get_value_for_lipid(self, lipid: "QuantifiedLipid") -> Optional[float]:
        """
        Retrieve the quantitation value for a given lipid.
        Reads directly from the lipid's values dict (no data duplication).
        Returns None if not found.
        """
        return lipid.values.get(self.sample_name)

    def mean_value_for_lipids(
        self,
        lipid: Union["QuantifiedLipid", str, List[Union["QuantifiedLipid", str]]],
        dataset: Optional["LipidDataset"] = None,
        skip_missing: bool = True,
    ) -> Optional[float]:
        """
        Compute the mean quantitation value for this sample across a list of lipids.

        Parameters:
        - lipids: list of `QuantifiedLipid` objects or lipid `input_name` strings.
        - dataset: optional `LipidDataset` used to resolve string names to lipid objects.
        - skip_missing: if True, missing values are ignored; if False, missing values are treated as NaN.

        Returns the mean as a float, or `None` if no valid numeric values were found.
        """
        if dataset is None:
            raise ValueError("dataset is required to resolve lipid names; call dataset.mean_value_for_lipids(...) instead")
        
        lipid_list = lipid if isinstance(lipid, list) else [lipid]
        return dataset.mean_value_for_lipids(
            sample=self,
            lipids=lipid_list,
            skip_missing=skip_missing,
        )


class QuantifiedLipid(LipidmapsBaseModel):
    input_name: str
    values: Dict[str, float]  # sample_name -> value
    pathway_ids: Optional[List[str]] = None  # e.g., KEGG or Reactome IDs
    pathway_names: Optional[List[str]] = None  # Human-readable names
    enzyme_ids: Optional[List[str]] = None  
    # RefMet annotations
    standardized_name: Optional[str] = None
    standardized_by: Optional[str] = None # e.g., "RefMet"
    lm_id: Optional[str] = None
    lm_id_found_by: Optional[str] = None  # e.g., "LMSD", "RefMet"
    matched_field: Optional[str] = None
    generic_lm_id: Optional[str] = None
    sub_class: Optional[str] = None
    super_class: Optional[str] = None
    main_class: Optional[str] = None
    chebi_id: Optional[str] = None
    kegg_id: Optional[str] = None
    refmet_id: Optional[str] = None
    formula: Optional[str] = None
    mass: Optional[float] = None
    reactions: Optional[List[ReactionData]] = None # List of associated reactions
    weight: Optional[float] = None  # For species or class-level reaction

    def get_value_for_sample(self, sample: "SampleMetadata") -> Optional[float]:
        """
        Retrieve the quantitation value for a given SampleMetadata object.
        Returns None if not found.
        """
        return self.values.get(sample.sample_name)

    
    @property
    def recognized(self) -> bool:
        return self.standardized_name is not None
    
    def zscore(self) -> Dict[str, float]:
        vals = np.array(list(self.values.values()))
        mean = np.mean(vals)
        std = np.std(vals)
        return {
            k: (v - mean) / std if std != 0 else 0.0 for k, v in self.values.items()
        }

class LipidDataset(LipidmapsBaseModel):
    samples: List[SampleMetadata]
    lipids: List[QuantifiedLipid]
    column_info: Optional[Dict[str, Any]] = None  # Metadata about CSV columns
    reactions: List[ReactionData] = Field(default_factory=list)  # All reactions in dataset
    validation_report: Optional[ValidationReport] = Field(default=None)

    def list_sample_names(self) -> List[str]:
        return [s.sample_name for s in self.samples]

    def list_lipid_names(self) -> List[str]:
        return [l.input_name for l in self.lipids]
    
    def list_reactions(self) -> Optional[List[str]]:
        return [f"{r.reaction_name} {r.reaction_id}" for r in (self.reactions)]
    
    def list_lipids_with_lmid(self) -> List[str]:
        return [l.input_name for l in self.lipids if l.lm_id is not None]
    
    def list_lm_ids(self) -> List[str]:
        # Collect lm_ids, skip missing values, and preserve original order
        lm_ids = [l.lm_id for l in self.lipids if l.lm_id is not None]
        seen = set()
        unique = []
        for lid in lm_ids:
            if lid not in seen:
                seen.add(lid)
                unique.append(lid)
        return unique

    def list_lipids_with_reactions(self) -> List[str]:
        return [l.input_name for l in self.lipids if l.reactions is not None and len(l.reactions) > 0]
    
    def get_lipids_with_reactions(self) -> List[QuantifiedLipid]:
        return [l for l in self.lipids if l.reactions is not None and len(l.reactions) > 0]

    def query_lipids(self, *preds: Union[Query, Callable[["QuantifiedLipid"], bool]], combine: str = "and") -> List[QuantifiedLipid]:
        """
        Query lipids using one or more predicates. Predicates may be `Query` objects
        or plain callables that accept a `QuantifiedLipid` and return truthy/falsy.

        combine: "and" (default) or "or" - how to combine multiple predicates.
        """
        if not preds:
            return list(self.lipids)

        qobjs: List[Query] = []
        for p in preds:
            if isinstance(p, Query):
                qobjs.append(p)
            elif callable(p):
                qobjs.append(from_callable(p))
            else:
                raise TypeError("Predicates must be Query or callable")

        combined = qobjs[0]
        for q in qobjs[1:]:
            combined = combined & q if combine == "and" else combined | q

        return [l for l in self.lipids if combined.matches(l)]

    def fetch_reactions_by_lm_id(self, reaction_type: Optional[str] = None, only_lipid_components: bool = True, taxonomy_group: Optional[str] = "all") -> List[ReactionData]:
        """
        Fetch reactions for LM IDs present in this dataset using the ReactionChecker API.
        Attaches the fetched reactions to `self.reactions` and annotates lipids in-place.

        Returns the list of fetched ReactionData objects (may be empty).
        """
        lm_ids = [lipid.lm_id for lipid in self.lipids if lipid.lm_id]
        generic_lm_ids = [lipid.generic_lm_id for lipid in self.lipids if lipid.generic_lm_id]
        lm_ids = list(set(lm_ids).union(set(generic_lm_ids)))

        if not lm_ids:
            logger.info("No LM IDs provided for reaction fetching.")
            self.reactions = []
            return []

        try:
            checker = ReactionChecker(base_url=LMSD_REACTIONS_BASE_URL)
            response = checker.check_reactions(
                lm_ids=lm_ids,
                generic_reactions=False,
                reaction_type=("class-level" if reaction_type is None else reaction_type),
                only_lipid_components=only_lipid_components,
                taxonomy_group=taxonomy_group
            )
            reactions = getattr(response, "reactions", []) or []
            logger.info(f"Fetched {len(reactions)} reactions from ReactionChecker API.")
            logger.debug(f"Reaction fetching parameters: lm_ids={lm_ids}, reaction_type={reaction_type}, only_lipid_components={only_lipid_components}, taxonomy_group={taxonomy_group}")
            # Deduplicate reactions by id and assign to dataset
            reaction_dict = {}
            for reaction in reactions:
                rid = getattr(reaction, "reaction_id", None) or getattr(reaction, "id", None)
                if rid is None:
                    continue
                rid = str(rid)
                if rid not in reaction_dict:
                    reaction_dict[rid] = reaction

            self.reactions = list(reaction_dict.values())

            # Annotate lipids with matching reactions
            for lipid in self.lipids:
                lipid_ids = set()
                if getattr(lipid, "lm_id", None):
                    lipid_ids.add(lipid.lm_id)
                if getattr(lipid, "generic_lm_id", None):
                    lipid_ids.add(lipid.generic_lm_id)

                matched_ids = set()
                matched_reactions = []
                for reaction in self.reactions:
                    for role in ["reactants", "products"]:
                        items = getattr(reaction, role, [])
                        for item in items:
                            compound_lm_id = None
                            if isinstance(item, dict):
                                compound_lm_id = item.get("compound_lm_id")
                            elif hasattr(item, "compound_lm_id"):
                                compound_lm_id = getattr(item, "compound_lm_id", None)
                            if compound_lm_id and compound_lm_id in lipid_ids:
                                rid = getattr(reaction, "reaction_id", None) or getattr(reaction, "id", None)
                                if rid is not None:
                                    rid = str(rid)
                                    if rid not in matched_ids:
                                        matched_reactions.append(reaction)
                                        matched_ids.add(rid)
                                break
                        else:
                            continue
                        break

                lipid.reactions = matched_reactions if matched_reactions else None

            return self.reactions
        except Exception:
            logger.exception("Failed to fetch reactions from ReactionChecker API.")
            self.reactions = []
            return []

    def get_lipids_for_component(self, component: CompoundComponent) -> List[QuantifiedLipid]:
        """
        Return all QuantifiedLipid objects where the component matches lm_id (case-insensitive).
        """
        comp_names = set()

        if hasattr(component, "compound_lm_id") and component.compound_lm_id:
            comp_names.add(component.compound_lm_id.lower())

        return [
            l for l in self.lipids
            if (l.lm_id and l.lm_id.lower() in comp_names)
            or (l.generic_lm_id and l.generic_lm_id.lower() in comp_names)
        ]

    def find_lipids(self, query: str) -> List[QuantifiedLipid]:
        q = query.lower()
        return [
            l for l in self.lipids
            if q in (l.input_name or "").lower() or (l.standardized_name and q in l.standardized_name.lower())
        ]

    def _find_headgroup_from_name(self, name: str) -> Optional[str]:
        """
        Find headgroup string using headgroup mapping.
        Returns:
           Headgroup string if a match is found, otherwise None.
        """
        if name.lower().startswith("fa ") or name.lower().startswith("fa("):
            return "FA"
        elif " O-" in name or " P-" in name:
            dash_index = name.index("-")
            return name[:dash_index+1]
        else:
            match = re.match(r"^([A-Za-z0-9\-]+)", name)
            if match:
                return match.group(1)
        return None
    
    def _find_generic_lm_id_from_name(self, name: str) -> Optional[str]:
        """
        Find Generic LM_ID string using headgroup mapping.
        Returns:
           Generic LM_ID string if a headgroup match is found, otherwise None.
        """
        headgroup = self._find_headgroup_from_name(name)
        
        if headgroup:
            lm_ids = lipidmaps_headgroups.get(headgroup)
            if lm_ids and lm_ids[0]:
                return lm_ids[0]
        return None
    
    def fill_generic_lm_ids_from_headgroups(self) -> int:
        """
        Fill missing lm_id fields on QuantifiedLipid objects using headgroup mapping from headgroups.py.
        Returns:
            Number of lipids updated with an lm_id.
        """
        updated = 0
        for lipid in self.lipids:
            if not getattr(lipid, "generic_lm_id", None):
                generic_lm_id = None
                if lipid.standardized_name:
                    generic_lm_id = self._find_generic_lm_id_from_name(lipid.standardized_name)
                if not generic_lm_id and (lipid.input_name.lower().startswith("fa ") or lipid.input_name.lower().startswith("fa(")):
                    generic_lm_id = self._find_generic_lm_id_from_name(lipid.input_name)
                if generic_lm_id:
                    lipid.generic_lm_id = generic_lm_id
                    lipid.lm_id_found_by = "headgroup"
                    updated += 1

        logger = logging.getLogger(__name__)
        logger.info(f"Updated {updated} generic_lm_id fields using headgroup mapping (via LipidDataset)")
        return updated

    def get_value(self, sample: "SampleMetadata", lipid: "QuantifiedLipid") -> Optional[float]:
        """
        Retrieve the quantitation value for a given sample and lipid object.
        Returns None if not found.
        """
        return lipid.values.get(sample.sample_name)

    def mean_value_for_lipids(
        self,
        sample: Union["SampleMetadata", str],
        lipids: List[Union["QuantifiedLipid", str]],
        skip_missing: bool = True,
    ) -> Optional[float]:
        """
        Compute the mean quantitation value for `sample` across a list of lipids.

        Parameters:
        - sample: `SampleMetadata` or sample_name string.
        - lipids: list of `QuantifiedLipid` objects or lipid `input_name` strings.
        - skip_missing: if True, missing values are ignored; if False, missing values are treated as NaN.

        Returns the mean as a float, or `None` if no valid numeric values were found.
        """
        # resolve sample name
        sample_name = sample.sample_name if hasattr(sample, "sample_name") else sample

        collected_values: List[float] = []
        for item in lipids:
            lipid_obj = None
            if isinstance(item, str):
                name_lc = item.lower()
                lipid_obj = next(
                    (
                        l
                        for l in self.lipids
                        if (l.input_name or "").lower() == name_lc
                        or (l.standardized_name and l.standardized_name.lower() == name_lc)
                    ),
                    None,
                )
                if lipid_obj is None:
                    continue
            else:
                lipid_obj = item

            value = lipid_obj.values.get(sample_name)
            if value is None:
                if skip_missing:
                    continue
                collected_values.append(np.nan)
            else:
                collected_values.append(value)

        if not collected_values:
            return None
        return float(np.nanmean(collected_values))

    def get_values(self, lipid_name: str) -> Optional[Dict[str, float]]:
        """
        Return the values dict for a lipid, matching input_name case-insensitively.
        Returns None if not found.
        """
        name_lc = lipid_name.lower()
        for l in self.lipids:
            if (l.input_name or '').lower() == name_lc:
                return l.values
        return None
        
    def get_grouped_data(self) -> Dict[str, List[QuantifiedLipid]]:
        grouped = {}
        for sample in self.samples:
            grouped.setdefault(sample.group, []).append(sample.sample_name)
        result = {}
        for group, sample_names in grouped.items():
            result[group] = [
                QuantifiedLipid(
                    input_name=lipid.input_name,
                    values={
                        sid: lipid.values[sid]
                        for sid in sample_names
                        if sid in lipid.values
                    },
                )
                for lipid in self.lipids
            ]
        return result
    
    def get_lipids_by_generic_lm_id(self, generic_lm_id: str) -> List[QuantifiedLipid]:
        """Return all QuantifiedLipid objects with the given generic lm_id."""
        return [l for l in self.lipids if l.generic_lm_id == generic_lm_id]

    def get_lipid_values_for_samples(self, sample_name: str) -> List[Dict[str, Optional[float]]]:
        """
        Return a list of objects for a given `sample_name` where each object contains:
            - `input_name`: the lipid's input name
            - `value`: the quantitation value for the provided sample_name (or None)

        This is useful for building per-sample plots or tables.
        """
        # maintain backward compatibility: accept either a SampleMetadata object or a sample name string
        if hasattr(sample_name, "sample_name"):
            sid = sample_name.sample_name
        else:
            sid = sample_name

        return [
            {"input_name": l.input_name, "value": l.values.get(sid)}
            for l in self.lipids
        ]

    def get_lipids_for_reaction(self, reaction_or_id, role: str = "reactant") -> List[QuantifiedLipid]:
        """
        Return QuantifiedLipid objects that are reactants or products in a given reaction.
        Accepts either a ReactionData object or a reaction_id (str).
        role: "reactant" or "product"
        """
        # Determine if input is a ReactionData object or an id
        if isinstance(reaction_or_id, ReactionData):
            reaction = reaction_or_id
        else:
            reaction = next((r for r in self.reactions if getattr(r, 'reaction_id', None) == reaction_or_id), None)
        if not reaction:
            return []
        lm_ids = getattr(reaction, 'reactant_lm_ids', []) if role == "reactant" else getattr(reaction, 'product_lm_ids', [])
        # If reaction uses CompoundComponent, adjust accordingly
        if not lm_ids and hasattr(reaction, 'reactants') and role == "reactant":
            lm_ids = [c.compound_lm_id for c in getattr(reaction, 'reactants', []) if hasattr(c, 'compound_lm_id')]
        if not lm_ids and hasattr(reaction, 'products') and role == "product":
            lm_ids = [c.compound_lm_id for c in getattr(reaction, 'products', []) if hasattr(c, 'compound_lm_id')]
        return [l for l in self.lipids if l.lm_id in lm_ids]



class Quantitation(LipidmapsBaseModel):
    lipid: "QuantifiedLipid"  # Reference to a QuantifiedLipid object
    sample_values: Dict[str, float]  # sample_name -> value
    method: Optional[str] = None  # e.g., 'LC-MS', 'GC-MS'
    unit: Optional[str] = None  # e.g., 'pmol', 'ng'
    notes: Optional[str] = None
    # Add more fields as needed

    def get_value_for_sample(self, sample_name: str) -> Optional[float]:
        return self.sample_values.get(sample_name)
    
if __name__ == "__main__":

    lipid = QuantifiedLipid(
        input_name="PC(16:0/18:1)",
        values={"Sample1": 10.2, "Sample2": 11.3, "Sample3": 9.8},
    )
    zscores = lipid.zscore()
    print(f"Z-scores: {zscores} {lipid}")
