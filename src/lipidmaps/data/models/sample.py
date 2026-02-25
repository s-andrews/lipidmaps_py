
import logging
from typing import Any, List, Dict, Optional, Union, Callable
import numpy as np
import re
from ..utils.headgroups import lipidmaps_headgroups
from .query import Query, from_callable, attr_eq, attr_in, attr_contains, attr_gt, has_attr
from .base import LipidmapsBaseModel
from pydantic import Field
from enum import Enum
from datetime import datetime
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

class MissingType(str, Enum):
    MCAR = "Missing Completely At Random"
    MAR = "Missing At Random"
    MNAR = "Missing Not At Random"
    OTHER = "OTHER"

class MissingInfo(LipidmapsBaseModel):
    type: MissingType
    note: Optional[str] = None
    imputed: bool = False

    model_config = {"extra": "forbid"}


class NormalizedResult(LipidmapsBaseModel):
    """Typed container for normalized values and provenance metadata."""

    method_key: str
    values: Dict[str, Optional[float]] = Field(default_factory=dict)
    meta: Dict[str, Any] = Field(default_factory=dict)
    created_iso: Optional[str] = None

    def model_post_init(self, __context: Any) -> None:  # type: ignore[misc]
        if self.created_iso is None:
            self.created_iso = datetime.utcnow().isoformat() + "Z"

class QuantifiedLipid(LipidmapsBaseModel):
    input_name: str
    values: Dict[str, Optional[float]] = Field(default_factory=dict)
    missing: Dict[str, MissingInfo] = Field(default_factory=dict)
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
    # Quantitation metadata
    unit: Optional[str] = None  # e.g., 'pmol', 'ng', 'area'
    measurement_method: Optional[str] = None  # e.g., 'LC-MS', 'GC-MS'
    quantitation_notes: Optional[str] = None
    # Store typed normalized results keyed by method_key -> NormalizedResult
    normalized: Dict[str, "NormalizedResult"] = Field(default_factory=dict)

    def get_value_for_sample(self, sample: SampleMetadata) -> Optional[float]:
        """
        Retrieve the quantitation value for a given SampleMetadata object.
        Returns None if not found.
        """
        return self.values.get(sample.sample_name)

    def set_missing_for_sample(self, sample: SampleMetadata, mtype: MissingType, note: Optional[str] = None, imputed: bool = False) -> None:
        sample_name = sample.sample_name
        self.missing[sample_name] = MissingInfo(type=mtype, note=note, imputed=imputed)
        # ensure explicit None value for missing entries
        if sample_name not in self.values:
            self.values[sample_name] = None

    def is_missing(self, sample: SampleMetadata) -> bool:
        return sample.sample_name in self.missing

    def impute_missing(self, strategy: str = "mean", fill_value: Optional[float] = None) -> Dict[str, float]:
        """
        Impute missing values in self.values.

        strategy: "mean", "median", "zero", or "value" (requires fill_value).
        Returns a dict of sample_name -> imputed_value for entries that were imputed.
        """
        # collect existing numeric values (ignore None and NaN)
        present = [v for v in self.values.values() if v is not None and not (isinstance(v, float) and np.isnan(v))]
        if strategy == "mean":
            fill = float(np.mean(present)) if present else (0.0 if fill_value is None else float(fill_value))
        elif strategy == "median":
            fill = float(np.median(present)) if present else (0.0 if fill_value is None else float(fill_value))
        elif strategy == "zero":
            fill = 0.0
        elif strategy == "value":
            if fill_value is None:
                raise ValueError("fill_value required for strategy='value'")
            fill = float(fill_value)
        else:
            raise ValueError(f"Unsupported imputation strategy: {strategy}")

        imputed: Dict[str, float] = {}
        for sid, val in list(self.values.items()):
            if val is None or (isinstance(val, float) and np.isnan(val)):
                self.values[sid] = fill
                mi = self.missing.get(sid)
                if mi:
                    mi.imputed = True
                    if mi.note is None:
                        mi.note = f"Imputed ({strategy})"
                else:
                    # Mark as imputed missing if no MissingInfo existed
                    self.missing[sid] = MissingInfo(type=MissingType.OTHER, note=f"Imputed ({strategy})", imputed=True)
                imputed[sid] = fill
        return imputed
        
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

    def set_normalized(self, entry: Union["NormalizedResult", str], vals: Optional[Dict[str, Optional[float]]] = None, meta: Optional[Dict[str, Any]] = None) -> None:
        """Store normalized values. Accepts either a NormalizedResult or (method_key, vals[, meta])."""
        if hasattr(entry, "method_key") and hasattr(entry, "values"):
            nr: "NormalizedResult" = entry  # type: ignore[var-annotated]
            self.normalized[nr.method_key] = nr
            return

        method_key = entry
        if vals is None:
            raise ValueError("vals required when passing method_key string")
        nr = NormalizedResult(method_key=method_key, values=vals, meta=meta or {})
        self.normalized[method_key] = nr

    def get_normalized(self, method_key: str) -> Optional["NormalizedResult"]:
        """Return a NormalizedResult for `method_key` or None if not present."""
        return self.normalized.get(method_key)

class LipidDataset(LipidmapsBaseModel):
    samples: List[SampleMetadata]
    lipids: List[QuantifiedLipid]
    column_info: Optional[Dict[str, Any]] = None  # Metadata about CSV columns
    reactions: List[ReactionData] = Field(default_factory=list)  # All reactions in dataset
    validation_report: Optional[ValidationReport] = Field(default=None)
    # Quantitation metadata for the entire dataset
    quantitation_unit: Optional[str] = None  # e.g., 'pmol', 'ng', 'area'
    quantitation_method: Optional[str] = None  # e.g., 'LC-MS', 'GC-MS'
    quantitation_notes: Optional[str] = None

    def set_quantitation_info(
        self,
        unit: Optional[str] = None,
        method: Optional[str] = None,
        notes: Optional[str] = None,
        apply_to_lipids: bool = False,
    ) -> None:
        """Set quantitation metadata for the dataset.

        Args:
            unit: Unit of measurement (e.g., 'pmol', 'ng', 'area')
            method: Measurement method (e.g., 'LC-MS', 'GC-MS')
            notes: Additional notes about the quantitation
            apply_to_lipids: If True, also set these values on all lipids
        """
        if unit is not None:
            self.quantitation_unit = unit
        if method is not None:
            self.quantitation_method = method
        if notes is not None:
            self.quantitation_notes = notes
        if apply_to_lipids:
            for lipid in self.lipids:
                if unit is not None:
                    lipid.unit = unit
                if method is not None:
                    lipid.measurement_method = method
                if notes is not None:
                    lipid.quantitation_notes = notes

    def list_sample_names(self) -> List[str]:
        return [s.sample_name for s in self.samples]

    def list_lipid_names(self) -> List[str]:
        return [l.input_name for l in self.lipids]

    def normalized_dataframe(self, method_key: str):
        """Return a pandas DataFrame of normalized values for `method_key`.

        Rows are lipids (input_name) and columns are sample names.
        """
        try:
            import pandas as pd
        except Exception:
            raise

        data = {}
        for l in self.lipids:
            nr = None
            try:
                nr = l.get_normalized(method_key)
            except Exception:
                nr = None

            if nr is None:
                # fallback to old-style dict stored directly on object
                try:
                    raw = getattr(l, 'normalized', {})
                    if isinstance(raw, dict) and method_key in raw:
                        data[l.input_name] = raw[method_key]
                except Exception:
                    pass
            else:
                # nr may be a NormalizedResult instance
                try:
                    data[l.input_name] = nr.values if hasattr(nr, 'values') else nr
                except Exception:
                    data[l.input_name] = nr

        # If no per-lipid normalized entries found, check dataset-level store
        if not data:
            store = getattr(self, '_normalized_store', None)
            if store and method_key in store:
                df = pd.DataFrame.from_dict(store[method_key], orient='index')
            else:
                df = pd.DataFrame.from_dict(data, orient='index')
        else:
            df = pd.DataFrame.from_dict(data, orient='index')
        # order columns if possible
        cols = self.list_sample_names() if hasattr(self, 'list_sample_names') else None
        if cols is not None and len(cols) > 0:
            try:
                df = df[cols]
            except Exception:
                pass
        return df
    
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

    def _find_generic_lm_id_from_headgroup(self, name: str) -> Optional[str]:
        """
        Find Generic LM_ID string using headgroup mapping.
        Returns:
           Generic LM_ID string if a headgroup match is found, otherwise None.
        """
        headgroup = None
        if name.lower().startswith("fa ") or name.lower().startswith("fa("):
            headgroup = "FA"
        elif " O-" in name or " P-" in name:
            dash_index = name.index("-")
            headgroup = name[:dash_index+1]
        else:
            match = re.match(r"^([A-Za-z0-9\-]+)", name)
            if match:
                headgroup = match.group(1)
        if headgroup:
            lm_ids = lipidmaps_headgroups.get(headgroup)
            if lm_ids and lm_ids[0]:
                return lm_ids[0]
        return None
    
    def fill_missing_lm_ids_from_headgroups(self) -> int:
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
                    generic_lm_id = self._find_generic_lm_id_from_headgroup(lipid.standardized_name)
                if not generic_lm_id and (lipid.input_name.lower().startswith("fa ") or lipid.input_name.lower().startswith("fa(")):
                    generic_lm_id = self._find_generic_lm_id_from_headgroup(lipid.input_name)
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

    def possible_reactions(self, reactions: Optional[List[Any]] = None) -> List[Any]:
        """Return reactions from `reactions` that are plausible for this dataset.

        Wraps the rule-based `reactions_possible_in_dataset` helper which inspects
        headgroups present in this dataset and filters based on known headgroup
        conversion rules.
        """
        try:
            if reactions is None:
                reactions = getattr(self, 'reactions', []) or []
            # local import to avoid top-level dependency cycles
            from ..utils.lipid_reaction_rules import reactions_possible_in_dataset

            return reactions_possible_in_dataset(self, reactions)
        except Exception:
            logger.exception("Failed to compute plausible reactions for dataset.")
            return []

    def impute_missing_values(self, strategy: str = "mean", fill_value: Optional[float] = None) -> Dict[str, Dict[str, float]]:
        """
        Impute missing values across all lipids in the dataset.
        Returns a mapping: lipid_input_name -> {sample_name: imputed_value, ...}
        """
        summary: Dict[str, Dict[str, float]] = {}
        for lipid in self.lipids:
            imputed = lipid.impute_missing(strategy=strategy, fill_value=fill_value)
            if imputed:
                summary[lipid.input_name] = imputed
        return summary



class Quantitation(LipidmapsBaseModel):
    lipid: Union["QuantifiedLipid", str]  # Reference to a QuantifiedLipid object or its input_name
    sample_values: Dict[str, Optional[float]]  # sample_name -> value (snapshot)
    method: Optional[str] = None  # e.g., 'LC-MS', 'GC-MS'
    unit: Optional[str] = None  # e.g., 'pmol', 'ng'
    notes: Optional[str] = None
    # Add more fields as needed

    @classmethod
    def from_lipid(cls, lipid: "QuantifiedLipid", samples: Optional[List[str]] = None, method: Optional[str] = None, unit: Optional[str] = None, notes: Optional[str] = None) -> "Quantitation":
        """
        Create a Quantitation snapshot from a QuantifiedLipid.
        If samples is provided, only include those sample names.
        """
        sample_values = {sid: val for sid, val in lipid.values.items() if (samples is None or sid in samples)}
        return cls(lipid=lipid.input_name, sample_values=sample_values, method=method, unit=unit, notes=notes)

    def get_value_for_sample(self, sample_name: str) -> Optional[float]:
        # If lipid was stored as object elsewhere, prefer snapshot; allow object resolution if needed
        if hasattr(self.lipid, "values"):
            return getattr(self.lipid, "values").get(sample_name)  # type: ignore
        return self.sample_values.get(sample_name)
    
if __name__ == "__main__":

    lipid = QuantifiedLipid(
        input_name="PC(16:0/18:1)",
        values={"Sample1": 10.2, "Sample2": 11.3, "Sample3": 9.8},
    )
    zscores = lipid.zscore()
    print(f"Z-scores: {zscores} {lipid}")
