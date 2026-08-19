
import logging
from typing import Any, List, Dict, Optional, Union, Callable
import numpy as np
import re
from ..utils.headgroups import lipidmaps_headgroups
from ..utils.lipid_reaction_rules import lipid_reaction_rules
from ..utils.chain_parser import LipidStructure, parse_lipid
from .query import Query, from_callable
from .base import LipidmapsBaseModel
from pydantic import ConfigDict, Field, PrivateAttr, RootModel, field_validator
from enum import Enum
from datetime import datetime
from .reaction import ReactionData, CompoundComponent, ReactionChecker
from ...config import LMSD_REACTIONS_BASE_URL
from ..validation.data_validator import ValidationReport



logger = logging.getLogger(__name__)


class SampleConditions(RootModel[Dict[str, str]]):
    """Mapping of sample_name to condition/group name."""

    root: Dict[str, str] = Field(default_factory=dict)

    @field_validator("root")
    @classmethod
    def validate_root(cls, value: Dict[str, str]) -> Dict[str, str]:
        normalized: Dict[str, str] = {}
        for sample_name, condition_name in value.items():
            sample_text = str(sample_name).strip() if sample_name is not None else ""
            condition_text = str(condition_name).strip() if condition_name is not None else ""
            if not sample_text:
                raise ValueError("Condition mapping keys must be non-empty sample names")
            if not condition_text:
                raise ValueError(
                    f"Condition name for sample '{sample_text}' must be a non-empty string"
                )
            normalized[sample_text] = condition_text
        return normalized

    def as_dict(self) -> Dict[str, str]:
        return dict(self.root)


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

    model_config = ConfigDict(extra="forbid")


class NormalizedResult(LipidmapsBaseModel):
    """Typed container for normalized values and provenance metadata."""

    method_key: str
    values: Dict[str, Optional[float]] = Field(default_factory=dict)
    meta: Dict[str, Any] = Field(default_factory=dict)
    created_iso: Optional[str] = None

    def model_post_init(self, __context: Any) -> None:  # type: ignore[misc]
        if self.created_iso is None:
            self.created_iso = datetime.utcnow().isoformat() + "Z"


class LipidAnnotation(LipidmapsBaseModel):
    """External database annotations and classification for a lipid.
    
    Contains classification hierarchy and database identifiers from sources
    like RefMet, CHEBI, KEGG, and LIPID MAPS.
    """
    # Classification hierarchy
    sub_class: Optional[str] = None  # e.g., "PC", "TAG"
    main_class: Optional[str] = None  # e.g., "Glycerophosphocholines"
    super_class: Optional[str] = None  # e.g., "Glycerophospholipids"
    
    # Database identifiers
    chebi_id: Optional[str] = None
    kegg_id: Optional[str] = None
    refmet_id: Optional[str] = None
    
    # Chemical properties
    formula: Optional[str] = None
    mass: Optional[float] = None


class QuantifiedLipid(LipidmapsBaseModel):
    input_name: str
    #TODO: We should have a headgroup class to refer here rather than duplication of headgroup string
    #Delegating headgroups to LipidStructure class
    #headgroup: Optional[str] = None
    standardized_name: Optional[str] = None # RefMet name if standardized
    standardized_by: Optional[str] = None # e.g., "RefMet"
    lm_id: Optional[str] = None
    lm_id_found_by: Optional[str] = None  # e.g., "LMSD", "RefMet"
    abbrev: Optional[str] = None  # e.g., "PC 34:1"
    abbrev_chains: Optional[str] = None  # e.g., "PC 16:0_18:1"
    values: Dict[str, Optional[float]] = Field(default_factory=dict)
    missing: Dict[str, MissingInfo] = Field(default_factory=dict) # If lipid value is missing for a sample, store info about the missingness here keyed by sample_name
    pathway_ids: Optional[List[str]] = None  # e.g., KEGG or Reactome IDs
    pathway_names: Optional[List[str]] = None  # Human-readable names
    enzyme_ids: Optional[List[str]] = None  


    matched_field: Optional[str] = None
    generic_lm_id: Optional[str] = None

    reactions: Optional[List[ReactionData]] = None # List of associated reactions
    weight: Optional[float] = None  # For species or class-level reaction
    # Quantitation metadata
    unit: Optional[str] = None  # e.g., 'pmol', 'ng', 'area'
    measurement_method: Optional[str] = None  # e.g., 'LC-MS', 'GC-MS'
    quantitation_notes: Optional[str] = None
    # Store typed normalized results keyed by method_key -> NormalizedResult
    normalized: Dict[str, "NormalizedResult"] = Field(default_factory=dict)


    # Composed structural data (parsed lazily or eagerly)
    _structure: Optional[LipidStructure] = PrivateAttr(default=None)
    # External annotations (populated from RefMet, LMSD, etc.)
    _annotation: Optional[LipidAnnotation] = PrivateAttr(default=None)
    
    @property
    def structure(self) -> LipidStructure:
        if self._structure is not None:
            return self._structure
        std = self.standardized_name
        std_lower = std.lower() if std is not None else ""
        if std is not None and "FA" not in self.input_name and "-coa" not in std_lower:
            self._structure = parse_lipid(std)
        elif "fa" in self.input_name.lower() or "-coa" in std_lower:
            self._structure = parse_lipid(self.input_name)
        return self._structure
    
    @property
    def annotation(self) -> Optional[LipidAnnotation]:
        """External database annotations (classification, IDs, properties)."""
        return self._annotation
    
    @annotation.setter
    def annotation(self, value: Optional[LipidAnnotation]) -> None:
        self._annotation = value
    
    #### METHODS FOR ACCESSING AND MANIPULATING QUANTITATION DATA ####

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
        def _present(v: Optional[float]) -> bool:
            return v is not None and not (isinstance(v, float) and np.isnan(v))

        present = [v for v in self.values.values() if _present(v)]
        if len(present) < 2:
            return {k: 0.0 for k in self.values}
        mean = float(np.mean(present))
        std = float(np.std(present, ddof=1))
        if std == 0:
            return {k: 0.0 for k in self.values}
        return {
            k: ((v - mean) / std) if _present(v) else 0.0 for k, v in self.values.items()
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
    # Status flags for API calls
    refmet_failed: bool = Field(default=False, description="True if RefMet API call failed during processing")
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

    def get_sample_conditions(self) -> SampleConditions:
        """Return the current sample_name -> condition mapping."""
        return SampleConditions.model_validate(
            {
                sample.sample_name: sample.group
                for sample in self.samples
                if sample.sample_name and sample.group
            }
        )

    def set_sample_conditions(
        self,
        conditions: Union[SampleConditions, Dict[str, str]],
        *,
        strict: bool = True,
    ) -> SampleConditions:
        """Update sample conditions from a partial or full sample_name -> condition mapping."""
        condition_map = (
            conditions
            if isinstance(conditions, SampleConditions)
            else SampleConditions.model_validate(conditions)
        )
        sample_lookup = {
            sample.sample_name: sample
            for sample in self.samples
            if sample.sample_name
        }
        unknown_samples = sorted(set(condition_map.root) - set(sample_lookup))
        if strict and unknown_samples:
            missing_text = ", ".join(unknown_samples)
            raise ValueError(f"Unknown sample names in condition mapping: {missing_text}")

        for sample_name, condition_name in condition_map.root.items():
            sample = sample_lookup.get(sample_name)
            if sample is not None:
                sample.group = condition_name

        return self.get_sample_conditions()

    def list_lipid_names(self) -> List[str]:
        return [l.input_name for l in self.lipids]

    def list_headgroups(self) -> List[str]:
        headgroups = set()
        for lipid in self.lipids:
            hg = None
            try:
                if hasattr(lipid, "structure") and lipid.structure and getattr(lipid.structure, "headgroup", None):
                    hg = lipid.structure.headgroup
            except Exception:
                hg = None
            if hg is not None:
                headgroups.add(hg)
        return sorted(headgroups)
    
    def list_generic_lm_ids(self) -> List[str]:
        generic_lm_ids = set()        
        for lipid in self.lipids:
            glid = getattr(lipid, "generic_lm_id", None)
            if glid is not None:
                generic_lm_ids.add(glid)
        return list(generic_lm_ids)
    
    def generic_lm_ids_exist(self, generic_lm_ids: List[str]) -> bool:
        existing = self.list_generic_lm_ids()
        return all(glid in existing for glid in generic_lm_ids)
    
    def headgroups_exist(self, headgroups: List[str]) -> bool:
        existing = self.list_headgroups()
        return all(hg in existing for hg in headgroups)
    
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

    def get_structures(self) -> List[LipidStructure]:
        """Return all parsed LipidStructure objects (triggers lazy parsing).
        
        Returns:
            List of LipidStructure objects for all lipids in the dataset.
        """
        structs: List[LipidStructure] = []
        for lipid in self.lipids:
            try:
                s = lipid.structure
            except Exception:
                s = None
            if s is not None:
                structs.append(s)
        return structs

    def get_structures_by_headgroup(self, headgroup: Optional[str] = None) -> Dict[str, List[LipidStructure]]:
        """Group LipidStructure objects by headgroup, or filter by headgroup if supplied."""
        result: Dict[str, List[LipidStructure]] = {}
        for lipid in self.lipids:
            try:
                s = lipid.structure
            except Exception:
                s = None
            if s is None:
                continue
            hg = getattr(s, "headgroup", None) or "Unknown"
            if headgroup is not None:
                if hg == headgroup:
                    result.setdefault(hg, []).append(s)
            else:
                result.setdefault(hg, []).append(s)
        if headgroup is not None:
            # Only return the dict for the requested headgroup
            return {headgroup: result.get(headgroup, [])}
        return result

    def check_structure_exist(self, headgroup: Optional[str] = None, total_carbons: Optional[int] = None, total_double_bonds: Optional[int] = None) -> bool:
        """LipidStructure object by headgroup, total_carbons, total_double_bonds"""
        if headgroup is None or total_carbons is None or total_double_bonds is None:
            return False
        result = False
        for lipid in self.lipids:
            try:
                s = lipid.structure
            except Exception:
                s = None
            if s is None:
                continue
            headgroup_origin = getattr(s, "headgroup", None)
            total_carbons_origin = getattr(s, "total_carbons", None)
            total_double_bonds_origin = getattr(s, "total_double_bonds", None)
            if headgroup_origin == headgroup and total_carbons_origin == total_carbons and total_double_bonds_origin == total_double_bonds:
                return True
            else:
                continue

        return result
        
    def fetch_reactions(self, lm_ids: List[str], reaction_type: Optional[str] = None, only_lipid_components: bool = True, taxonomy_group: Optional[str] = "all") -> List[ReactionData]:
        """
        Fetch reactions for a list of LM IDs using the ReactionChecker API.

        This is a focused helper which only performs the network call and
        returns a deduplicated list of ReactionData objects. It does NOT
        mutate the dataset or annotate lipids.
        """
        if not lm_ids:
            return []

        try:
            checker = ReactionChecker(base_url=LMSD_REACTIONS_BASE_URL)
            response = checker.check_reactions(
                lm_ids=lm_ids,
                generic_reactions=False,
                reaction_type=("class-level" if reaction_type is None else reaction_type),
                only_lipid_components=only_lipid_components,
                taxonomy_group=taxonomy_group,
            )
            reactions = getattr(response, "reactions", []) or []

            # Deduplicate by id
            reaction_dict = {}
            for reaction in reactions:
                rid = getattr(reaction, "reaction_id", None) or getattr(reaction, "id", None)
                if rid is None:
                    continue
                rid = str(rid)
                if rid not in reaction_dict:
                    reaction_dict[rid] = reaction

            return list(reaction_dict.values())
        except Exception:
            logger.exception("Failed to fetch reactions from ReactionChecker API.")
            return []



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

    def fetch_reactions_by_lm_id(self, reaction_type: Optional[str] = None, only_lipid_components: bool = True, taxonomy_group: Optional[str] = "all", annotate_lipids: bool = True, annotate_using_evaluator: bool = True) -> List[ReactionData]:
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

        # Use the focused helper to perform the network fetch and deduplication
        reactions = self.fetch_reactions(lm_ids=lm_ids, reaction_type=reaction_type, only_lipid_components=only_lipid_components, taxonomy_group=taxonomy_group)
        logger.info(f"Fetched {len(reactions)} reactions from ReactionChecker API.")

        # Assign to dataset
        self.reactions = reactions

        if annotate_lipids and reactions:
            # Annotate lipids with matching reactions (simple LM ID matching)
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

        if annotate_using_evaluator and self.reactions:
            try:
                from ..utils.reaction_evaluator import ReactionEvaluator

                evaluator = ReactionEvaluator(lipid_reaction_rules)
                evaluator.annotate_reactions(self.reactions, dataset=self)
            except Exception:
                logger.exception("ReactionEvaluator failed; leaving reactions unannotated")

        return self.reactions

    def get_lipids_for_component(self, component: CompoundComponent) -> List[QuantifiedLipid]:
        """
        Return all QuantifiedLipid objects where the component matches lm_id (case-insensitive).
        """
        comp_names = set()

        if hasattr(component, "compound_lm_id") and component.compound_lm_id:
            comp_names.add(component.compound_lm_id.lower())
        if hasattr(component, "compound_generic_lm_id") and component.compound_generic_lm_id:
            comp_names.add(component.compound_generic_lm_id.lower())

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

    def _find_linkage_from_name(self, name: str) -> Optional[str]:
        """
        Detect linkage type from lipid name patterns.
        Returns one of: 'ester', 'ether_alkyl', 'ether_vinyl', 'amide', or None.
        """
        if not name:
            return None
        n = name.strip()
        # plasmalogen marker (P-)
        if " P-" in n or n.startswith("P-"):
            return "ether_vinyl"
        # alkyl ether marker (O-)
        if " O-" in n or n.startswith("O-"):
            return "ether_alkyl"
        # sphingolipids often have 'Cer' or 'SM' in name
        if n.startswith("Cer") or n.startswith("SM") or n.lower().startswith("dhcer"):
            return "amide"
        # default: assume ester for glycerophospholipids
        return "ester"

    def _get_headgroup_from_compound(self, comp: CompoundComponent) -> Optional[str]:
        if not comp:
            return None
        # Prefer explicit field
        if getattr(comp, "compound_headgroup", None):
            return comp.compound_headgroup
        # Fall back to generic id or name parsing
        if getattr(comp, "compound_generic_lm_id", None):
            return comp.compound_generic_lm_id
        if getattr(comp, "compound_name", None):
            return self._find_headgroup_from_name(comp.compound_name)
        return None

    def _get_linkage_from_compound(self, comp: CompoundComponent) -> Optional[str]:
        # Prefer explicit abbrev markers
        if not comp:
            return None
        if getattr(comp, "compound_abbrev", None):
            return self._find_linkage_from_name(comp.compound_abbrev)
        if getattr(comp, "compound_headgroup", None):
            return self._find_linkage_from_name(comp.compound_headgroup)
        if getattr(comp, "compound_name", None):
            return self._find_linkage_from_name(comp.compound_name)
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

    def fill_headgroups_from_names(self) -> int:
        """
        Count lipids that have a headgroup parsed from their structure.
        
        Note: headgroup is now a computed property from LipidStructure,
        parsed automatically from input_name. This method counts lipids
        with successfully parsed headgroups rather than setting them.
        
        Returns:
            Number of lipids with a parsed headgroup.
        """
        count = 0
        for lipid in self.lipids:
            # Access structure.headgroup to check parsing
            if lipid.structure.headgroup:
                count += 1

        logger = logging.getLogger(__name__)
        logger.info(f"Found {count} lipids with parsed headgroup (via LipidDataset)")
        return count

    def fill_compound_headgroups_from_lipids(self) -> int:
        """
        Fill missing `compound_headgroup` fields on `CompoundComponent` objects
        inside `self.reactions` by matching `compound_lm_id` (or
        `compound_generic_lm_id`) against `QuantifiedLipid.lm_id` /
        `QuantifiedLipid.generic_lm_id` and copying `QuantifiedLipid.headgroup`.

        Only updates components where `compound_type == 'lm_main'` and
        `compound_headgroup` is not already set.

        Returns:
            Number of components updated.
        """
        updated = 0
        if not getattr(self, 'reactions', None):
            return 0

        # build quick lookup by lm_id and generic_lm_id
        lookup = {}
        for l in self.lipids:
            if getattr(l, 'lm_id', None):
                lookup.setdefault(str(l.lm_id), []).append(l)
            if getattr(l, 'generic_lm_id', None):
                lookup.setdefault(str(l.generic_lm_id), []).append(l)

        for rx in self.reactions:
            for comp in (getattr(rx, 'reactants', []) or []) + (getattr(rx, 'products', []) or []):
                try:
                    ctype = getattr(comp, 'compound_type', None)
                    if ctype != 'lm_main':
                        continue
                    if getattr(comp, 'compound_headgroup', None):
                        continue
                    lm = getattr(comp, 'compound_lm_id', None) or getattr(comp, 'compound_generic_lm_id', None)
                    if not lm:
                        continue
                    # exact-string lookup (use str keys)
                    matches = lookup.get(str(lm)) or []
                    if not matches:
                        continue
                    # prefer first match with headgroup set
                    for ml in matches:
                        try:
                            structure = getattr(ml, 'structure', None)
                            headgroup = getattr(structure, 'headgroup', None)
                        except Exception:
                            headgroup = None
                        if headgroup:
                            comp.compound_headgroup = headgroup
                            updated += 1
                            break
                except Exception:
                    continue

        logger = logging.getLogger(__name__)
        logger.info(f"Updated {updated} compound headgroups from quantified lipids")
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
        """Return reactions from `reactions` that are plausible for this dataset."""
        try:
            if reactions is None:
                reactions = getattr(self, 'reactions', []) or []

            return [
                reaction
                for reaction in reactions
                if getattr(reaction, 'evaluation', None)
                and reaction.evaluation.get('possible') is True
            ]
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
