import json
import math
import re
from collections import OrderedDict
from pathlib import Path
from typing import Any, ClassVar, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from pydantic import ConfigDict, Field, PrivateAttr
from scipy import stats

from .matching import match_pathway_reactions
from .models.base import LipidmapsBaseModel
from .models.reaction import ReactionData
from .models.sample import LipidDataset, QuantifiedLipid
from .models.species_reaction import ClassReaction, CompoundRequirement, PathwayReactionSet, ReactionMatchResult, SpeciesReactionPair
from .models.uniprot import UniProtRheaClient
from .utils.chain_parser import (
    extract_facoa_from_reactions,
    extract_fa_from_reactions,
    get_common_facoa_names,
    get_common_fa_names,
    infer_facoa_from_lipids,
    infer_fa_from_lipids,
)
from .utils.headgroups import lipidmaps_headgroups, lm_id_to_headgroup


class BioPANPathwayExporter(LipidmapsBaseModel):
    """Build BioPAN reaction graph, table, and edge-detail assets from a dataset."""

    COMPARISON_BUNDLE_NAME: ClassVar[str] = "comparison_bundle.json"
    EDGE_DETAILS_BUNDLE_NAME: ClassVar[str] = "edge_details_bundle.json"

    # Sphingolipid backbone classes whose base + N-acyl chain pass through
    # multi-substrate reactions (e.g. SM synthase) unchanged, so a same-structure
    # backbone edge (Cer->SM, dhCer->dhSM) can be derived from them.
    SPHINGO_BACKBONE_CLASSES: ClassVar[frozenset] = frozenset({"Cer", "dhCer", "SM", "dhSM"})

    dataset: Optional[Any] = Field(default=None)
    class_reactions: Optional[List[Any]] = Field(default=None)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    _edge_cache: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = PrivateAttr(default_factory=dict)
    _reaction_table_cache: Dict[Tuple[Any, ...], Dict[str, Any]] = PrivateAttr(default_factory=dict)
    _pathway_table_cache: Dict[Tuple[Any, ...], Dict[str, Any]] = PrivateAttr(default_factory=dict)
    _reaction_gene_lookup: Optional[Dict[str, List[str]]] = PrivateAttr(default=None)
    _uniprot_gene_lookup: Dict[str, List[str]] = PrivateAttr(default_factory=dict)
    _uniprot_client: UniProtRheaClient = PrivateAttr(default_factory=UniProtRheaClient)

    def _get_dataset(self, dataset: Optional[LipidDataset] = None) -> LipidDataset:
        resolved_dataset = dataset or self.dataset
        if resolved_dataset is None:
            raise ValueError("BioPAN pathway export requires a populated dataset")
        return resolved_dataset

    @staticmethod
    def _get_output_dir(output_path: Union[str, Path]) -> Path:
        output_dir = Path(output_path)
        if output_dir.name != "biopan":
            output_dir = output_dir / "biopan"
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    @staticmethod
    def _paired_suffix(paired: bool) -> str:
        return "paired" if paired else "notpaired"

    @staticmethod
    def _safe_node_id(value: str) -> str:
        return re.sub(r"[^A-Za-z0-9]+", "", value).lower()

    @staticmethod
    def _safe_edge_id(source: str, target: str) -> str:
        return re.sub(r"[^A-Za-z0-9]+", "", f"{source}{target}").upper()

    @staticmethod
    def _structure_label(value: Any) -> str:
        if value is None:
            return ""
        headgroup = getattr(value, "headgroup", None)
        if not headgroup:
            return str(value)
        if getattr(value, "backbone", None) and getattr(value, "chains", None):
            chain_text = "/".join(str(chain) for chain in value.chains)
            return f"{headgroup}({value.backbone}/{chain_text})"
        if getattr(value, "chains", None) and getattr(value, "level", None) and str(value.level).lower().endswith("full"):
            chain_text = "_".join(str(chain) for chain in value.chains)
            return f"{headgroup}({chain_text})"
        if getattr(value, "total_carbons", None) is not None and getattr(value, "total_double_bonds", None) is not None:
            return f"{headgroup}({value.total_carbons}:{value.total_double_bonds})"
        return str(value)

    @staticmethod
    def _get_lipid_display_name(lipid: QuantifiedLipid) -> str:
        if getattr(lipid, "standardized_name", None):
            return lipid.standardized_name
        if lipid.input_name:
            return lipid.input_name
        structure = getattr(lipid, "structure", None)
        return BioPANPathwayExporter._structure_label(structure)

    def _display_name_for_structure(
        self,
        structure: Any,
        lipid_lookup: Dict[str, QuantifiedLipid],
    ) -> str:
        structure_label = self._structure_label(structure)
        lipid = lipid_lookup.get(structure_label)
        if lipid is not None:
            return self._get_lipid_display_name(lipid)
        return structure_label

    def _build_lipid_lookup(self, dataset: LipidDataset) -> Dict[str, QuantifiedLipid]:
        lookup: Dict[str, QuantifiedLipid] = {}
        for lipid in dataset.lipids:
            display_name = self._get_lipid_display_name(lipid)
            if display_name:
                lookup.setdefault(display_name, lipid)
            structure = getattr(lipid, "structure", None)
            structure_label = self._structure_label(structure)
            if structure_label:
                lookup.setdefault(structure_label, lipid)
        return lookup

    def _load_reaction_gene_lookup(self) -> Dict[str, List[str]]:
        if self._reaction_gene_lookup is not None:
            return self._reaction_gene_lookup

        lookup: Dict[str, List[str]] = {}
        sql_path = Path(__file__).resolve().parents[3] / "reactions_genes.sql"
        if not sql_path.exists():
            self._reaction_gene_lookup = lookup
            return lookup

        in_copy_block = False
        with sql_path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not in_copy_block:
                    if line.startswith("COPY public.reactions_genes "):
                        in_copy_block = True
                    continue

                if line == "\\.":
                    break

                parts = line.split("\t")
                if len(parts) < 5:
                    continue

                reaction_id = parts[1].strip()
                gene_name = parts[4].strip()
                if not reaction_id or not gene_name:
                    continue

                lookup.setdefault(reaction_id, [])
                if gene_name not in lookup[reaction_id]:
                    lookup[reaction_id].append(gene_name)

        self._reaction_gene_lookup = lookup
        return lookup

    @staticmethod
    def _reaction_identifier(reaction: ReactionData) -> Optional[str]:
        for key in ("reaction_id", "id"):
            value = getattr(reaction, key, None)
            if value is not None:
                return str(value).strip()
        return None

    @staticmethod
    def _entry_value(entry: Any, key: str) -> Any:
        if isinstance(entry, dict):
            return entry.get(key)
        return getattr(entry, key, None)

    def _get_reaction_gene_symbols(self, reaction: ReactionData) -> List[str]:
        rhea_ids = reaction.list_rhea_ids()
        if rhea_ids:
            names: List[str] = []
            seen = set()
            for rhea_id in rhea_ids:
                genes = self._uniprot_gene_lookup.get(rhea_id)
                if genes is None:
                    genes = self._uniprot_client.fetch_gene_symbols(rhea_id)
                    self._uniprot_gene_lookup[rhea_id] = genes
                for gene in genes:
                    if gene and gene not in seen:
                        names.append(gene)
                        seen.add(gene)
            return names

        names: List[str] = []
        for entry in getattr(reaction, "genes", []) or []:
            for key in (
                "gene_symbol",
                "gene_name",
                "gene",
                "symbol",
                "name",
                "predicted_gene_symbol",
                "predicted_gene",
                "predicted_symbol",
            ):
                value = self._entry_value(entry, key)
                if value:
                    names.append(str(value).strip())
                    break
        for entry in getattr(reaction, "proteins", []) or []:
            for key in (
                "gene_symbol",
                "gene_name",
                "gene",
                "symbol",
                "protein_gene_name",
                "protein_gene_symbol",
                "predicted_gene_symbol",
                "predicted_gene",
                "predicted_symbol",
            ):
                value = self._entry_value(entry, key)
                if value:
                    names.append(str(value).strip())
                    break

        reaction_id = self._reaction_identifier(reaction)
        if reaction_id:
            names.extend(self._load_reaction_gene_lookup().get(reaction_id, []))

        ordered: List[str] = []
        seen = set()
        for name in names:
            if name and name not in seen:
                ordered.append(name)
                seen.add(name)
        return ordered

    def _component_headgroup(self, component: Any) -> Optional[str]:
        if component is None:
            return None
        explicit_headgroup = getattr(component, "compound_headgroup", None)
        if explicit_headgroup:
            return explicit_headgroup
        for lm_key in ("compound_generic_lm_id", "compound_lm_id"):
            lm_id = getattr(component, lm_key, None)
            if lm_id and lm_id in lm_id_to_headgroup:
                return lm_id_to_headgroup[lm_id]
        for name_key in ("compound_abbrev", "compound_name"):
            value = getattr(component, name_key, None)
            if value:
                text = str(value).strip()
                if text.lower().startswith("coa") or "facoa" in text.lower():
                    return "acyl CoA"
                if text.startswith("FA(") or text.startswith("FA "):
                    return "FA"
                match = re.match(r"^([A-Za-z0-9\- ]+?)(?:\(|\s|$)", text)
                if match:
                    return match.group(1).strip()
        return None

    def _infer_compound_requirement(self, reaction: ReactionData) -> Tuple[CompoundRequirement, bool]:
        reactant_headgroups = {self._component_headgroup(component) for component in reaction.reactants}
        product_headgroups = {self._component_headgroup(component) for component in reaction.products}
        if "acyl CoA" in reactant_headgroups:
            return CompoundRequirement.FACOA, True
        if "FA" in product_headgroups:
            return CompoundRequirement.FA, False
        return CompoundRequirement.NONE, False

    def _extract_class_reactions(self, dataset: LipidDataset) -> Tuple[List[ClassReaction], Dict[str, List[ReactionData]]]:
        if self.class_reactions:
            empty_lookup = {reaction.reaction_key.lower(): [] for reaction in self.class_reactions}
            return self.class_reactions, empty_lookup

        reaction_lookup: Dict[str, List[ReactionData]] = {}
        reaction_map: Dict[str, ClassReaction] = {}

        def register(reaction, reactant_class, product_class, compound_require, acyl_add):
            key = f"{reactant_class},{product_class}".lower()
            if key not in reaction_map:
                reaction_map[key] = ClassReaction(
                    reactant_class=reactant_class,
                    product_class=product_class,
                    reaction_class="Matched reactions",
                    compound_require=compound_require,
                    acyl_add=acyl_add,
                    genes=self._get_reaction_gene_symbols(reaction),
                )
            else:
                merged = reaction_map[key].genes + self._get_reaction_gene_symbols(reaction)
                seen = set()
                reaction_map[key].genes = [gene for gene in merged if gene and not (gene in seen or seen.add(gene))]
            reaction_lookup.setdefault(key, []).append(reaction)

        for reaction in getattr(dataset, "reactions", []) or []:
            reactants = [component for component in reaction.reactants if self._component_headgroup(component) not in {None, "FA", "acyl CoA"}]
            products = [component for component in reaction.products if self._component_headgroup(component) not in {None, "FA", "acyl CoA"}]

            if len(reactants) == 1 and len(products) == 1:
                reactant_class = self._component_headgroup(reactants[0])
                product_class = self._component_headgroup(products[0])
                if not reactant_class or not product_class:
                    continue
                compound_require, acyl_add = self._infer_compound_requirement(reaction)
                register(reaction, reactant_class, product_class, compound_require, acyl_add)
                continue

            # Multi-substrate reactions such as sphingomyelin synthase
            # ('PC; Ceramide -> 1,2-DG; Sphingomyelin') are skipped by the 1:1
            # rule above, so the sphingolipid backbone edge BioPAN tracks
            # (Cer->SM, dhCer->dhSM) is otherwise lost. The sphingoid base and
            # N-acyl chain pass through unchanged, so derive that same-structure
            # backbone edge. Fire only with exactly one sphingolipid on each side
            # to avoid mis-pairing.
            sphingo_reactants = [hg for hg in (self._component_headgroup(c) for c in reactants) if hg in self.SPHINGO_BACKBONE_CLASSES]
            sphingo_products = [hg for hg in (self._component_headgroup(c) for c in products) if hg in self.SPHINGO_BACKBONE_CLASSES]
            if len(sphingo_reactants) == 1 and len(sphingo_products) == 1 and sphingo_reactants[0] != sphingo_products[0]:
                register(reaction, sphingo_reactants[0], sphingo_products[0], CompoundRequirement.NONE, False)

        return list(reaction_map.values()), reaction_lookup

    def _get_matching_inputs(
        self,
        dataset: LipidDataset,
    ) -> Tuple[List[str], List[str], List[str]]:
        lipid_names = [self._get_lipid_display_name(lipid) for lipid in dataset.lipids if self._get_lipid_display_name(lipid)]

        # Union all available sources rather than taking the first non-empty
        # one. Reaction extraction often returns only a sparse set of LM-ID'd
        # fatty acids (e.g. FA 16:0 + 2:0), which would otherwise short-circuit
        # the richer pool inferable from the dataset's own acyl chains and limit
        # FA-release / FACoA-addition matching to that sparse subset. Common
        # defaults are used only when nothing can be sourced from the dataset.
        fa_names = self._ordered_union(
            extract_fa_from_reactions(dataset.reactions) if dataset.reactions else [],
            [name for name in lipid_names if name.startswith("FA(") or name.startswith("FA ")],
            infer_fa_from_lipids(lipid_names),
        ) or get_common_fa_names()

        facoa_names = self._ordered_union(
            extract_facoa_from_reactions(dataset.reactions) if dataset.reactions else [],
            [name for name in lipid_names if name.startswith("CoA ") or name.startswith("FACoA") or name.startswith("FaCoA")],
            infer_facoa_from_lipids(lipid_names),
        ) or get_common_facoa_names()

        return lipid_names, fa_names, facoa_names

    @staticmethod
    def _ordered_union(*groups: Sequence[str]) -> List[str]:
        """Concatenate name groups, preserving first-seen order and dropping dups."""
        seen = set()
        ordered: List[str] = []
        for group in groups:
            for name in group:
                if name and name not in seen:
                    seen.add(name)
                    ordered.append(name)
        return ordered

    def build_reaction_match_set(
        self,
        dataset: Optional[LipidDataset] = None,
    ) -> Tuple[PathwayReactionSet, Dict[str, List[ReactionData]]]:
        resolved_dataset = self._get_dataset(dataset)
        class_reactions, reaction_lookup = self._extract_class_reactions(resolved_dataset)
        lipid_names, fa_names, facoa_names = self._get_matching_inputs(resolved_dataset)
        result_set = match_pathway_reactions(
            lipid_names=lipid_names,
            reactions=class_reactions,
            fa_names=fa_names,
            facoa_names=facoa_names,
            use_full_structure=True,
        )
        result_set.dataset_lipid_count = len(resolved_dataset.lipids)
        return result_set, reaction_lookup

    def _group_samples(self, dataset: LipidDataset, group: str) -> List[str]:
        return [sample.sample_name for sample in dataset.samples if sample.group == group and sample.sample_name]

    @staticmethod
    def _round_zscore(z_score: float) -> float:
        if z_score == 0:
            return 0.0
        return float(f"{round(z_score, 3):.3f}")

    @staticmethod
    def _critical_zscore(p_value: float) -> float:
        return float(stats.norm.ppf(1 - p_value))

    @staticmethod
    def _mode_rank_value(score: float, mode: str) -> float:
        return score

    @staticmethod
    def _log_ratio(product_value: Optional[float], reactant_value: Optional[float]) -> Optional[float]:
        if product_value is None or reactant_value is None:
            return None
        numerator = max(product_value, 0.0) + 1e-9
        denominator = max(reactant_value, 0.0) + 1e-9
        return math.log2(numerator / denominator)

    @staticmethod
    def _clean_series(values: Sequence[Optional[float]]) -> List[float]:
        cleaned: List[float] = []
        for value in values:
            if value is None:
                continue
            if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
                continue
            cleaned.append(float(value))
        return cleaned

    def _sum_lipid_values(
        self,
        lipid_names: Sequence[str],
        sample_names: Sequence[str],
        lipid_lookup: Dict[str, QuantifiedLipid],
    ) -> List[Optional[float]]:
        series: List[Optional[float]] = []
        for sample_name in sample_names:
            total = 0.0
            found = False
            for lipid_name in lipid_names:
                lipid = lipid_lookup.get(lipid_name)
                if lipid is None:
                    continue
                value = lipid.values.get(sample_name)
                if value is None:
                    continue
                if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
                    continue
                total += float(value)
                found = True
            series.append(total if found else None)
        return series

    @staticmethod
    def _is_missing(value: Optional[float]) -> bool:
        """True for None / NaN / Inf aggregates."""
        if value is None:
            return True
        return isinstance(value, float) and (math.isnan(value) or math.isinf(value))

    def _select_ratio_vectors(
        self,
        disease_products: Sequence[Optional[float]],
        disease_reactants: Sequence[Optional[float]],
        control_products: Sequence[Optional[float]],
        control_reactants: Sequence[Optional[float]],
        paired: bool,
    ) -> Optional[Tuple[List[float], List[float]]]:
        """Build the disease/control ratio vectors fed to the t-test.

        Mirrors the legacy R BioPAN rules (lib_pathway_analysis.r,
        ``get_lipid_pathway_zscore`` / ``extract_zscore``):

        1. If *any* per-sample product or reactant sum is missing (None / NaN /
           Inf), the reaction is unscorable and ``None`` is returned (R's
           ``na_values == 0`` guard -> z-score 0).
        2. Samples whose reactant sum is 0 are dropped. For a balanced, unpaired
           design R removes the *union* of zero-reactant positions from both
           groups; this is reproduced when the two groups have equal length.
           For unequal groups (where R's cross-indexing is undefined) each group
           is filtered independently.
        3. If the number of dropped positions is >= n - 1 the reaction is
           unscorable (R's ``l < length(...) - 1`` guard).

        Returns ``(disease_ratios, control_ratios)`` or ``None`` to signal a
        zero z-score.
        """
        vectors = (disease_products, disease_reactants, control_products, control_reactants)
        if any(self._is_missing(value) for vector in vectors for value in vector):
            return None

        if not paired and len(disease_reactants) == len(control_reactants):
            size = len(disease_reactants)
            dropped = {
                index
                for index in range(size)
                if disease_reactants[index] == 0 or control_reactants[index] == 0
            }
            if dropped and len(dropped) >= size - 1:
                return None
            disease_ratios = [
                float(disease_products[i]) / float(disease_reactants[i])
                for i in range(size) if i not in dropped
            ]
            control_ratios = [
                float(control_products[i]) / float(control_reactants[i])
                for i in range(size) if i not in dropped
            ]
            return disease_ratios, control_ratios

        if paired:
            size = min(len(disease_reactants), len(control_reactants))
            indices = [
                index for index in range(size)
                if disease_reactants[index] != 0 and control_reactants[index] != 0
            ]
            disease_ratios = [float(disease_products[i]) / float(disease_reactants[i]) for i in indices]
            control_ratios = [float(control_products[i]) / float(control_reactants[i]) for i in indices]
            return disease_ratios, control_ratios

        disease_ratios = [
            float(product) / float(reactant)
            for product, reactant in zip(disease_products, disease_reactants) if reactant != 0
        ]
        control_ratios = [
            float(product) / float(reactant)
            for product, reactant in zip(control_products, control_reactants) if reactant != 0
        ]
        return disease_ratios, control_ratios

    def _compute_ratio_zscore(
        self,
        disease_products: Sequence[Optional[float]],
        disease_reactants: Sequence[Optional[float]],
        control_products: Sequence[Optional[float]],
        control_reactants: Sequence[Optional[float]],
        alt: str,
        paired: bool,
    ) -> float:
        selected = self._select_ratio_vectors(
            disease_products, disease_reactants, control_products, control_reactants, paired
        )
        if selected is None:
            return 0.0
        disease_values, control_values = selected

        if len(disease_values) <= 1 or len(control_values) <= 1:
            return 0.0

        try:
            if paired:
                test_result = stats.ttest_rel(disease_values, control_values, alternative=alt)
            else:
                test_result = stats.ttest_ind(disease_values, control_values, equal_var=False, alternative=alt)
        except Exception:
            return 0.0

        p_value = getattr(test_result, "pvalue", None)
        if p_value is None or math.isnan(p_value) or math.isinf(p_value):
            return 0.0

        z_score = float(stats.norm.ppf(1 - p_value))
        if math.isnan(z_score):
            return 0.0
        return self._round_zscore(z_score)

    def _score_pair(
        self,
        pair: SpeciesReactionPair,
        lipid_lookup: Dict[str, QuantifiedLipid],
        disease_samples: Sequence[str],
        control_samples: Sequence[str],
        alt: str,
        paired: bool,
    ) -> float:
        reactant_name = self._display_name_for_structure(pair.reactant, lipid_lookup)
        product_name = self._display_name_for_structure(pair.product, lipid_lookup)
        disease_products = self._sum_lipid_values([product_name], disease_samples, lipid_lookup)
        disease_reactants = self._sum_lipid_values([reactant_name], disease_samples, lipid_lookup)
        control_products = self._sum_lipid_values([product_name], control_samples, lipid_lookup)
        control_reactants = self._sum_lipid_values([reactant_name], control_samples, lipid_lookup)
        return self._compute_ratio_zscore(
            disease_products,
            disease_reactants,
            control_products,
            control_reactants,
            alt=alt,
            paired=paired,
        )

    def _score_result(
        self,
        result: ReactionMatchResult,
        lipid_lookup: Dict[str, QuantifiedLipid],
        disease_samples: Sequence[str],
        control_samples: Sequence[str],
        alt: str,
        paired: bool,
    ) -> float:
        if not result.pairs:
            return 0.0
        reactant_names = sorted({self._display_name_for_structure(pair.reactant, lipid_lookup) for pair in result.pairs})
        product_names = sorted({self._display_name_for_structure(pair.product, lipid_lookup) for pair in result.pairs})
        disease_products = self._sum_lipid_values(product_names, disease_samples, lipid_lookup)
        disease_reactants = self._sum_lipid_values(reactant_names, disease_samples, lipid_lookup)
        control_products = self._sum_lipid_values(product_names, control_samples, lipid_lookup)
        control_reactants = self._sum_lipid_values(reactant_names, control_samples, lipid_lookup)
        return self._compute_ratio_zscore(
            disease_products,
            disease_reactants,
            control_products,
            control_reactants,
            alt=alt,
            paired=paired,
        )

    def _diagnose_ratio(
        self,
        reactant_names: Sequence[str],
        product_names: Sequence[str],
        lipid_lookup: Dict[str, QuantifiedLipid],
        disease_samples: Sequence[str],
        control_samples: Sequence[str],
        alt: str,
        paired: bool,
    ) -> Dict[str, Any]:
        """Expose every intermediate behind a single ratio z-score.

        Mirrors `_compute_ratio_zscore` step for step so the per-sample product
        and reactant sums, the per-sample ratios, the cleaned vectors actually
        fed to the t-test, and the resulting t-statistic / p-value / z-score can
        all be diffed against the legacy R BioPAN output for the same reaction.
        """
        reactant_names = sorted(set(reactant_names))
        product_names = sorted(set(product_names))

        disease_products = self._sum_lipid_values(product_names, disease_samples, lipid_lookup)
        disease_reactants = self._sum_lipid_values(reactant_names, disease_samples, lipid_lookup)
        control_products = self._sum_lipid_values(product_names, control_samples, lipid_lookup)
        control_reactants = self._sum_lipid_values(reactant_names, control_samples, lipid_lookup)

        def _ratios(products: Sequence[Optional[float]], reactants: Sequence[Optional[float]]) -> List[Optional[float]]:
            ratios: List[Optional[float]] = []
            for product, reactant in zip(products, reactants):
                if product is None or reactant is None or reactant == 0:
                    ratios.append(None)
                else:
                    ratios.append(float(product) / float(reactant))
            return ratios

        disease_ratios = _ratios(disease_products, disease_reactants)
        control_ratios = _ratios(control_products, control_reactants)

        # Use the exact selection routine that scoring uses, so the reported
        # t-test inputs match the z-score below. None signals an unscorable
        # reaction (missing aggregate or too many zero-reactant samples).
        selected = self._select_ratio_vectors(
            disease_products, disease_reactants, control_products, control_reactants, paired
        )
        disease_used, control_used = selected if selected is not None else ([], [])

        diagnostic: Dict[str, Any] = {
            "reactant_species": reactant_names,
            "product_species": product_names,
            "per_sample": {
                "disease": [
                    {"sample": sample, "product_sum": product, "reactant_sum": reactant, "ratio": ratio}
                    for sample, product, reactant, ratio in zip(disease_samples, disease_products, disease_reactants, disease_ratios)
                ],
                "control": [
                    {"sample": sample, "product_sum": product, "reactant_sum": reactant, "ratio": ratio}
                    for sample, product, reactant, ratio in zip(control_samples, control_products, control_reactants, control_ratios)
                ],
            },
            "disease_ratios_used": disease_used,
            "control_ratios_used": control_used,
            "n_disease_used": len(disease_used),
            "n_control_used": len(control_used),
            "t_statistic": None,
            "p_value": None,
        }

        if len(disease_used) > 1 and len(control_used) > 1:
            try:
                if paired:
                    test_result = stats.ttest_rel(disease_used, control_used, alternative=alt)
                else:
                    test_result = stats.ttest_ind(disease_used, control_used, equal_var=False, alternative=alt)
                diagnostic["t_statistic"] = float(getattr(test_result, "statistic", float("nan")))
                diagnostic["p_value"] = float(getattr(test_result, "pvalue", float("nan")))
            except Exception as exc:  # pragma: no cover - defensive
                diagnostic["t_error"] = str(exc)

        diagnostic["z_score"] = self._compute_ratio_zscore(
            disease_products,
            disease_reactants,
            control_products,
            control_reactants,
            alt=alt,
            paired=paired,
        )
        return diagnostic

    def diagnose_reaction(
        self,
        reaction_key: str,
        disease_group: str,
        control_group: str,
        *,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
    ) -> Dict[str, Any]:
        """Dump the inputs and intermediate values behind one reaction z-score.

        Use this to diff dev against the legacy R BioPAN tool for a single
        reaction (e.g. ``diagnose_reaction("PS,LPS", disease, control)``).

        At ``level='class'`` the diagnostic mirrors `_score_result` (sums every
        matched reactant/product species). At ``level='species'`` it returns one
        block per matched pair, mirroring `_score_pair`.
        """
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None:
            result_set, _ = self.build_reaction_match_set(resolved_dataset)

        key = reaction_key.lower()
        result = result_set.results.get(key)
        if result is None:
            return {
                "reaction_key": reaction_key,
                "error": "reaction not found in match set",
                "available_keys": sorted(result_set.results),
            }

        lipid_lookup = self._build_lipid_lookup(resolved_dataset)
        disease_samples = self._group_samples(resolved_dataset, disease_group)
        control_samples = self._group_samples(resolved_dataset, control_group)
        alt = "greater" if mode in {"active", "most_active"} else "less"

        summary: Dict[str, Any] = {
            "reaction_key": reaction_key,
            "level": level,
            "mode": mode,
            "alt": alt,
            "paired": paired,
            "disease_group": disease_group,
            "control_group": control_group,
            "disease_samples": list(disease_samples),
            "control_samples": list(control_samples),
            "pairs_matched": len(result.pairs),
            "fa_pool": sorted(result_set.fa_species),
            "facoa_pool": sorted(result_set.facoa_species),
        }

        if level == "species":
            pair_diags: List[Dict[str, Any]] = []
            for pair in result.pairs:
                reactant_name = self._display_name_for_structure(pair.reactant, lipid_lookup)
                product_name = self._display_name_for_structure(pair.product, lipid_lookup)
                diag = self._diagnose_ratio(
                    [reactant_name], [product_name], lipid_lookup,
                    disease_samples, control_samples, alt, paired,
                )
                pair_diags.append(diag)
            summary["pairs"] = pair_diags
            return summary

        reactant_names = sorted({self._display_name_for_structure(pair.reactant, lipid_lookup) for pair in result.pairs})
        product_names = sorted({self._display_name_for_structure(pair.product, lipid_lookup) for pair in result.pairs})
        summary.update(
            self._diagnose_ratio(
                reactant_names, product_names, lipid_lookup,
                disease_samples, control_samples, alt, paired,
            )
        )
        return summary

    def _gene_text(self, result: ReactionMatchResult, reaction_lookup: Dict[str, List[ReactionData]]) -> str:
        genes = list(result.class_reaction.genes)
        for reaction in reaction_lookup.get(result.reaction_key, []):
            genes.extend(self._get_reaction_gene_symbols(reaction))
        ordered: List[str] = []
        seen = set()
        for gene in genes:
            if gene and gene not in seen:
                ordered.append(gene)
                seen.add(gene)
        return ",".join(ordered) if ordered else "NA"

    @staticmethod
    def _unique_pathway_type(value: Any) -> List[str]:
        if value is None:
            return []
        if isinstance(value, (list, tuple, set)):
            values = value
        else:
            values = [value]
        ordered: List[str] = []
        seen = set()
        for item in values:
            text = str(item).strip()
            if text and text not in seen:
                ordered.append(text)
                seen.add(text)
        return ordered

    def _get_result_pathways(
        self,
        result: ReactionMatchResult,
        reaction_lookup: Dict[str, List[ReactionData]],
    ) -> List[Dict[str, Any]]:
        ordered: List[Dict[str, Any]] = []
        seen = set()
        for reaction in reaction_lookup.get(result.reaction_key, []):
            for entry in getattr(reaction, "pathways", []) or []:
                if not isinstance(entry, dict):
                    continue
                pathway_name = str(entry.get("pathway_name") or entry.get("name") or "").strip()
                pathway_types = self._unique_pathway_type(entry.get("pathway_type") or entry.get("type") or entry.get("classification"))
                if not pathway_name and not pathway_types:
                    continue
                key = (pathway_name, tuple(pathway_types))
                if key in seen:
                    continue
                seen.add(key)
                ordered.append({
                    "pathway_name": pathway_name,
                    "pathway_type": pathway_types,
                })
        return ordered

    @staticmethod
    def _format_pathway_labels(pathways: Sequence[Dict[str, Any]]) -> str:
        labels: List[str] = []
        seen = set()
        for entry in pathways:
            pathway_name = str(entry.get("pathway_name") or "").strip()
            categories = [str(item).strip() for item in entry.get("pathway_type", []) if str(item).strip()]
            if pathway_name and categories:
                label = f"{pathway_name} ({'/'.join(categories)})"
            else:
                label = pathway_name or "/".join(categories)
            if label and label not in seen:
                labels.append(label)
                seen.add(label)
        return ", ".join(labels) if labels else "NA"

    def _get_pathway_keys(self, pathways: Sequence[Dict[str, Any]]) -> List[str]:
        keys: List[str] = []
        seen = set()
        for entry in pathways:
            pathway_name = str(entry.get("pathway_name") or "").strip()
            categories = [str(item).strip() for item in entry.get("pathway_type", []) if str(item).strip()]
            label = pathway_name or " / ".join(categories)
            if label and label not in seen:
                keys.append(label)
                seen.add(label)
        return keys

    @staticmethod
    def _merge_genes(gene_lists: Sequence[Sequence[str]]) -> str:
        ordered: List[str] = []
        seen = set()
        for gene_list in gene_lists:
            for gene in gene_list:
                gene_text = str(gene).strip()
                if gene_text and gene_text != "NA" and gene_text not in seen:
                    ordered.append(gene_text)
                    seen.add(gene_text)
        return ",".join(ordered) if ordered else "NA"

    def _edge_gene_list(self, result: ReactionMatchResult, reaction_lookup: Dict[str, List[ReactionData]]) -> List[str]:
        gene_text = self._gene_text(result, reaction_lookup)
        if not gene_text or gene_text == "NA":
            return []
        return [entry.strip() for entry in gene_text.split(",") if entry.strip()]

    def _build_reaction_edges(
        self,
        disease_group: str,
        control_group: str,
        level: str,
        dataset: LipidDataset,
        result_set: PathwayReactionSet,
        reaction_lookup: Dict[str, List[ReactionData]],
        paired: bool,
        mode: str,
    ) -> List[Dict[str, Any]]:
        lipid_lookup = self._build_lipid_lookup(dataset)
        disease_samples = self._group_samples(dataset, disease_group)
        control_samples = self._group_samples(dataset, control_group)
        alt = "greater" if mode in {"active", "most_active"} else "less"
        edges: List[Dict[str, Any]] = []

        for result in result_set.results.values():
            if not result.has_pairs:
                continue
            edge_genes = self._edge_gene_list(result, reaction_lookup)
            edge_pathways = self._get_result_pathways(result, reaction_lookup)
            pathway_keys = self._get_pathway_keys(edge_pathways)
            pathway_label = self._format_pathway_labels(edge_pathways)

            if level == "class":
                reactant = result.class_reaction.reactant_class
                product = result.class_reaction.product_class
                score = self._score_result(result, lipid_lookup, disease_samples, control_samples, alt=alt, paired=paired)
                edges.append({
                    "source_label": reactant,
                    "target_label": product,
                    "source_id": self._safe_node_id(reactant),
                    "target_id": self._safe_node_id(product),
                    "source_lm_id": lipidmaps_headgroups.get(reactant, [""])[0] if reactant in lipidmaps_headgroups else "",
                    "target_lm_id": lipidmaps_headgroups.get(product, [""])[0] if product in lipidmaps_headgroups else "",
                    "edge_id": self._safe_edge_id(reactant, product),
                    "score": score,
                    "genes": edge_genes,
                    "pathways": edge_pathways,
                    "pathway_keys": pathway_keys,
                    "pathway_label": pathway_label,
                })
                continue

            for pair in result.pairs:
                reactant_label = self._display_name_for_structure(pair.reactant, lipid_lookup)
                product_label = self._display_name_for_structure(pair.product, lipid_lookup)
                score = self._score_pair(pair, lipid_lookup, disease_samples, control_samples, alt=alt, paired=paired)
                edges.append({
                    "source_label": reactant_label,
                    "target_label": product_label,
                    "source_id": self._safe_node_id(reactant_label),
                    "target_id": self._safe_node_id(product_label),
                    "source_lm_id": "",
                    "target_lm_id": "",
                    "edge_id": self._safe_edge_id(reactant_label, product_label),
                    "score": score,
                    "genes": edge_genes,
                    "pathways": edge_pathways,
                    "pathway_keys": pathway_keys,
                    "pathway_label": pathway_label,
                })
        return edges

    def _get_reaction_edges_cached(
        self,
        disease_group: str,
        control_group: str,
        level: str,
        dataset: LipidDataset,
        result_set: PathwayReactionSet,
        reaction_lookup: Dict[str, List[ReactionData]],
        paired: bool,
        mode: str,
    ) -> List[Dict[str, Any]]:
        cache_key = (
            id(dataset),
            disease_group,
            control_group,
            level,
            paired,
            mode,
        )
        cached = self._edge_cache.get(cache_key)
        if cached is None:
            cached = self._build_reaction_edges(
                disease_group=disease_group,
                control_group=control_group,
                level=level,
                dataset=dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
                paired=paired,
                mode=mode,
            )
            self._edge_cache[cache_key] = cached
        return cached

    def _select_significant_edges(
        self,
        edges: Sequence[Dict[str, Any]],
        p_value: float,
        mode: str,
    ) -> List[Dict[str, Any]]:
        cutoff = self._critical_zscore(p_value)
        return [
            edge
            for edge in edges
            if self._mode_rank_value(edge["score"], mode) > cutoff
        ]

    def _chain_score(self, scores: Sequence[float]) -> float:
        if not scores:
            return 0.0
        return self._round_zscore(sum(scores) / math.sqrt(len(scores)))

    def _pathway_chain_rows(
        self,
        edges: Sequence[Dict[str, Any]],
        subset: str,
        p_value: float,
        mode: str,
        max_depth: int = 6,
        max_rows: int = 2000,
        beam_width: int = 24,
    ) -> List[Dict[str, Any]]:
        cutoff = self._critical_zscore(p_value)
        adjacency: Dict[str, List[Dict[str, Any]]] = {}
        label_lookup: Dict[str, str] = {}
        state_scores: Dict[Tuple[str, Tuple[str, ...]], List[float]] = {}
        for edge in edges:
            adjacency.setdefault(edge["source_id"], []).append(edge)
            label_lookup[edge["source_id"]] = edge["source_label"]
            label_lookup[edge["target_id"]] = edge["target_label"]

        for edge_list in adjacency.values():
            edge_list.sort(key=lambda edge: self._mode_rank_value(edge["score"], mode), reverse=True)

        chains: Dict[str, Dict[str, Any]] = {}

        def walk(
            node_id: str,
            node_labels: List[str],
            chain_edges: List[Dict[str, Any]],
            shared_keys: Optional[set[str]],
            visited: set[str],
        ) -> None:
            if len(chains) >= max_rows:
                return
            for edge in adjacency.get(node_id, []):
                target_id = edge["target_id"]
                if target_id in visited:
                    continue

                edge_keys = set(edge["pathway_keys"])
                next_shared = set(edge_keys) if shared_keys is None else set(shared_keys) & edge_keys
                if subset == "pathway" and not next_shared:
                    continue

                next_edges = chain_edges + [edge]
                next_labels = node_labels + [edge["target_label"]]
                score = self._chain_score([item["score"] for item in next_edges])
                if self._mode_rank_value(score, mode) <= cutoff:
                    continue

                state_key = (target_id, tuple(sorted(next_shared)))
                prior_scores = state_scores.setdefault(state_key, [])
                ranked_score = self._mode_rank_value(score, mode)
                if len(prior_scores) >= beam_width and ranked_score <= min(prior_scores):
                    continue
                prior_scores.append(ranked_score)
                prior_scores.sort(reverse=True)
                if len(prior_scores) > beam_width:
                    del prior_scores[beam_width:]

                chain_key = "&#8594;".join(next_labels)
                prefix_key = "&#8594;".join(node_labels)
                if prefix_key in chains:
                    del chains[prefix_key]

                shared_pathways = [pathway for pathway in edge["pathways"] if (pathway.get("pathway_name") or " / ".join(pathway.get("pathway_type", []))) in next_shared]
                pathway_label = self._format_pathway_labels(shared_pathways)
                genes = self._merge_genes([item["genes"] for item in next_edges])
                chains[chain_key] = {
                    "data": {
                        "pathway": chain_key,
                        "score": score,
                        "gene": genes,
                        "class": pathway_label,
                    },
                    "edge_ids": [item["edge_id"] for item in next_edges],
                    "node_ids": [item["source_id"] for item in next_edges] + [next_edges[-1]["target_id"]],
                    "edge_scores": [item["score"] for item in next_edges],
                }

                if len(next_edges) >= max_depth or len(chains) >= max_rows:
                    continue

                walk(
                    target_id,
                    next_labels,
                    next_edges,
                    next_shared,
                    visited | {target_id},
                )

        for start_id in sorted(adjacency):
            walk(start_id, [label_lookup[start_id]], [], None, {start_id})

        rows = list(chains.values())
        rows.sort(key=lambda row: self._mode_rank_value(row["data"]["score"], mode), reverse=True)
        return rows

    def _most_significant_rows(self, rows: Sequence[Dict[str, Any]], mode: str) -> List[Dict[str, Any]]:
        selected: Dict[str, Dict[str, Any]] = {}
        for row in rows:
            node_ids = row.get("node_ids", [])
            if len(node_ids) < 2:
                continue
            root = node_ids[0]
            current = selected.get(root)
            if current is None:
                selected[root] = row
                continue

            current_nodes = current.get("node_ids", [])
            current_scores = current.get("edge_scores", [])
            new_scores = row.get("edge_scores", [])
            replace = False
            same_length = min(len(current_nodes), len(node_ids))
            for index in range(1, same_length):
                if node_ids[index] != current_nodes[index]:
                    new_score = new_scores[index - 1] if index - 1 < len(new_scores) else float("-inf")
                    current_score = current_scores[index - 1] if index - 1 < len(current_scores) else float("-inf")
                    replace = self._mode_rank_value(new_score, mode) > self._mode_rank_value(current_score, mode)
                    break
            if not replace and self._mode_rank_value(row["data"]["score"], mode) > self._mode_rank_value(current["data"]["score"], mode):
                replace = True
            if replace:
                selected[root] = row

        result = list(selected.values())
        result.sort(key=lambda row: self._mode_rank_value(row["data"]["score"], mode), reverse=True)
        return result

    def _rows_to_highlight(self, rows: Sequence[Dict[str, Any]]) -> Dict[str, str]:
        node_ids: List[str] = []
        edge_ids: List[str] = []
        for row in rows:
            node_ids.extend(row.get("node_ids", []))
            edge_ids.extend(row.get("edge_ids", []))
        unique_nodes = sorted(set(node_ids))
        unique_edges = sorted(set(edge_ids))
        return {
            "nodes": ",".join(f"#{node_id}" for node_id in unique_nodes),
            "edges": ",".join(f"#{edge_id}" for edge_id in unique_edges),
        }

    def _iter_pathway_results(
        self,
        result_set: PathwayReactionSet,
        reaction_lookup: Dict[str, List[ReactionData]],
    ) -> Iterable[Tuple[ReactionMatchResult, List[Dict[str, Any]]]]:
        for result in result_set.results.values():
            if not result.has_pairs:
                continue
            pathways = self._get_result_pathways(result, reaction_lookup)
            if pathways:
                yield result, pathways

    def build_reaction_tree(
        self,
        level: str = "class",
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
    ) -> List[Dict[str, Any]]:
        resolved_dataset = self._get_dataset(dataset)
        resolved_result_set = result_set or self.build_reaction_match_set(resolved_dataset)[0]
        category = {"text": "Matched reactions", "children": []}

        if level == "class":
            seen_classes = set()
            for result in resolved_result_set.results.values():
                if not result.has_pairs:
                    continue
                for class_name in (result.class_reaction.reactant_class, result.class_reaction.product_class):
                    if class_name not in seen_classes:
                        category["children"].append({"text": class_name})
                        seen_classes.add(class_name)
            category["children"].sort(key=lambda item: item["text"])
            return [category]

        lipid_lookup = self._build_lipid_lookup(resolved_dataset)
        grouped_species: Dict[str, List[str]] = {}
        for result in resolved_result_set.results.values():
            if not result.has_pairs:
                continue
            for pair in result.pairs:
                grouped_species.setdefault(pair.reactant.headgroup, []).append(self._display_name_for_structure(pair.reactant, lipid_lookup))
                grouped_species.setdefault(pair.product.headgroup, []).append(self._display_name_for_structure(pair.product, lipid_lookup))

        for class_name in sorted(grouped_species):
            unique_species = sorted(set(grouped_species[class_name]))
            category["children"].append({
                "text": class_name,
                "children": [{"text": species} for species in unique_species],
            })
        return [category]

    def build_pathway_tree(
        self,
        level: str = "class",
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
    ) -> List[Dict[str, Any]]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)

        lipid_lookup = self._build_lipid_lookup(resolved_dataset)
        grouped: "OrderedDict[str, OrderedDict[str, set[str]]]" = OrderedDict()
        for result, pathways in self._iter_pathway_results(result_set, reaction_lookup):
            categories: List[str] = []
            for pathway in pathways:
                categories.extend(pathway.get("pathway_type", []))
            ordered_categories: List[str] = []
            seen = set()
            for category in categories:
                if category and category not in seen:
                    ordered_categories.append(category)
                    seen.add(category)
            for category in ordered_categories:
                category_bucket = grouped.setdefault(category, OrderedDict())
                if level == "class":
                    for class_name in (result.class_reaction.reactant_class, result.class_reaction.product_class):
                        category_bucket.setdefault(class_name, set())
                else:
                    for pair in result.pairs:
                        category_bucket.setdefault(pair.reactant.headgroup, set()).add(self._display_name_for_structure(pair.reactant, lipid_lookup))
                        category_bucket.setdefault(pair.product.headgroup, set()).add(self._display_name_for_structure(pair.product, lipid_lookup))

        tree: List[Dict[str, Any]] = []
        for category, classes in grouped.items():
            children: List[Dict[str, Any]] = []
            for class_name in sorted(classes):
                species = classes[class_name]
                if level == "class":
                    children.append({"text": class_name})
                else:
                    children.append({
                        "text": class_name,
                        "children": [{"text": name} for name in sorted(species)],
                    })
            tree.append({"text": category, "children": children})
        return tree

    def build_reaction_graph(
        self,
        disease_group: str,
        control_group: str,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
    ) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        edge_rows = self._get_reaction_edges_cached(
            disease_group=disease_group,
            control_group=control_group,
            level=level,
            dataset=resolved_dataset,
            result_set=result_set,
            reaction_lookup=reaction_lookup,
            paired=paired,
            mode=mode,
        )

        node_map: Dict[str, Dict[str, Any]] = {}
        edges: List[Dict[str, Any]] = []
        for row in edge_rows:
            source_name = row["source_label"] if level == "class" else ""
            target_name = row["target_label"] if level == "class" else ""
            node_map.setdefault(row["source_id"], {"data": {"id": row["source_id"], "lm_id": row["source_lm_id"], "label": row["source_label"], "name": source_name, "shape": "ellipse"}})
            node_map.setdefault(row["target_id"], {"data": {"id": row["target_id"], "lm_id": row["target_lm_id"], "label": row["target_label"], "name": target_name, "shape": "ellipse"}})
            edges.append({
                "data": {
                    "id": row["edge_id"],
                    "source": row["source_id"],
                    "target": row["target_id"],
                    "weight": row["score"],
                    "color": "#24a19c" if row["score"] > 0 else "#A634C7",
                }
            })
        return {"nodes": list(node_map.values()), "edges": edges}

    def build_pathway_graph(
        self,
        disease_group: str,
        control_group: str,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
    ) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        edge_rows = [
            edge for edge in self._get_reaction_edges_cached(
                disease_group=disease_group,
                control_group=control_group,
                level=level,
                dataset=resolved_dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
                paired=paired,
                mode=mode,
            )
            if edge["pathway_keys"]
        ]

        node_map: Dict[str, Dict[str, Any]] = {}
        edges: List[Dict[str, Any]] = []
        for row in edge_rows:
            source_name = row["source_label"] if level == "class" else {}
            target_name = row["target_label"] if level == "class" else {}
            node_map.setdefault(row["source_id"], {"data": {"id": row["source_id"], "lm_id": row["source_lm_id"], "label": row["source_label"], "name": source_name, "shape": "ellipse"}})
            node_map.setdefault(row["target_id"], {"data": {"id": row["target_id"], "lm_id": row["target_lm_id"], "label": row["target_label"], "name": target_name, "shape": "ellipse"}})
            edges.append({
                "data": {
                    "id": row["edge_id"],
                    "source": row["source_id"],
                    "target": row["target_id"],
                    "weight": row["score"],
                    "color": "#24a19c" if row["score"] > 0 else "#A634C7",
                }
            })
        return {"nodes": list(node_map.values()), "edges": edges}

    def build_reaction_highlight(
        self,
        disease_group: str,
        control_group: str,
        threshold: float,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
    ) -> Dict[str, str]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        rows = self.build_reaction_table(
            disease_group=disease_group,
            control_group=control_group,
            threshold=threshold,
            level=level,
            mode=mode,
            paired=paired,
            dataset=resolved_dataset,
            result_set=result_set,
            reaction_lookup=reaction_lookup,
        )
        return self._rows_to_highlight(rows["pathways"])

    def build_reaction_table(
        self,
        disease_group: str,
        control_group: str,
        threshold: float,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
        limit: Optional[int] = None,
    ) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        cache_key = (
            id(resolved_dataset),
            disease_group,
            control_group,
            threshold,
            level,
            mode,
            paired,
        )
        cached = self._reaction_table_cache.get(cache_key)
        if cached is None:
            edges = self._get_reaction_edges_cached(
                disease_group=disease_group,
                control_group=control_group,
                level=level,
                dataset=resolved_dataset,
                result_set=result_set,
                reaction_lookup=reaction_lookup,
                paired=paired,
                mode=mode,
            )
            selected_edges = self._select_significant_edges(edges, threshold, mode)
            rows: List[Dict[str, Any]] = []
            for edge in selected_edges:
                rows.append({
                    "data": {
                        "pathway": f"{edge['source_label']}&#8594;{edge['target_label']}",
                        "score": edge["score"],
                        "gene": self._merge_genes([edge["genes"]]),
                    },
                    "node_ids": [edge["source_id"], edge["target_id"]],
                    "edge_ids": [edge["edge_id"]],
                    "edge_scores": [edge["score"]],
                })

            rows.sort(key=lambda row: self._mode_rank_value(row["data"]["score"], mode), reverse=True)
            if mode in {"most_active", "most_suppressed"}:
                rows = self._most_significant_rows(rows, mode)
            if limit is not None:
                rows = rows[:limit]
            cached = {"pathways": rows}
            self._reaction_table_cache[cache_key] = cached
        rows = list(cached["pathways"])
        if limit is not None:
            rows = rows[:limit]
        return {"pathways": rows}

    def build_pathway_highlight(
        self,
        disease_group: str,
        control_group: str,
        threshold: float,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
    ) -> Dict[str, str]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        rows = self.build_pathway_table(
            disease_group=disease_group,
            control_group=control_group,
            threshold=threshold,
            level=level,
            mode=mode,
            paired=paired,
            dataset=resolved_dataset,
            result_set=result_set,
            reaction_lookup=reaction_lookup,
        )
        return self._rows_to_highlight(rows["pathways"])

    def build_pathway_table(
        self,
        disease_group: str,
        control_group: str,
        threshold: float,
        level: str = "class",
        mode: str = "active",
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
        reaction_lookup: Optional[Dict[str, List[ReactionData]]] = None,
        limit: Optional[int] = None,
    ) -> Dict[str, Any]:
        resolved_dataset = self._get_dataset(dataset)
        if result_set is None or reaction_lookup is None:
            result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        cache_key = (
            id(resolved_dataset),
            disease_group,
            control_group,
            threshold,
            level,
            mode,
            paired,
        )
        cached = self._pathway_table_cache.get(cache_key)
        if cached is None:
            edges = [
                edge
                for edge in self._get_reaction_edges_cached(
                    disease_group=disease_group,
                    control_group=control_group,
                    level=level,
                    dataset=resolved_dataset,
                    result_set=result_set,
                    reaction_lookup=reaction_lookup,
                    paired=paired,
                    mode=mode,
                )
                if edge["pathway_keys"]
            ]
            rows = self._pathway_chain_rows(edges, subset="pathway", p_value=threshold, mode=mode)
            if mode in {"most_active", "most_suppressed"}:
                rows = self._most_significant_rows(rows, mode)
            cached = {"pathways": rows}
            self._pathway_table_cache[cache_key] = cached
        rows = list(cached["pathways"])
        if limit is not None:
            rows = rows[:limit]
        return {"pathways": rows}

    def build_edge_details(
        self,
        level: str = "class",
        dataset: Optional[LipidDataset] = None,
        result_set: Optional[PathwayReactionSet] = None,
    ) -> Dict[str, Dict[str, str]]:
        resolved_dataset = self._get_dataset(dataset)
        resolved_result_set = result_set
        if resolved_result_set is None:
            resolved_result_set, _ = self.build_reaction_match_set(resolved_dataset)
        lipid_lookup = self._build_lipid_lookup(resolved_dataset)
        by_headgroup = resolved_dataset.get_structures_by_headgroup()
        details: Dict[str, Dict[str, str]] = {}

        for result in resolved_result_set.results.values():
            if not result.has_pairs:
                continue

            selected_reactants = sorted({self._display_name_for_structure(pair.reactant, lipid_lookup) for pair in result.pairs})
            selected_products = sorted({self._display_name_for_structure(pair.product, lipid_lookup) for pair in result.pairs})
            reactant_pool = sorted({self._display_name_for_structure(structure, lipid_lookup) for structure in by_headgroup.get(result.class_reaction.reactant_class, [])})
            product_pool = sorted({self._display_name_for_structure(structure, lipid_lookup) for structure in by_headgroup.get(result.class_reaction.product_class, [])})
            edge_key = ",".join(
                [
                    self._safe_node_id(result.class_reaction.reactant_class),
                    self._safe_node_id(result.class_reaction.product_class),
                ]
            )
            details[edge_key] = {
                "non_selected_re": ",".join([name for name in reactant_pool if name not in selected_reactants]),
                "selected_re": ",".join(selected_reactants),
                "non_selected_pro": ",".join([name for name in product_pool if name not in selected_products]),
                "selected_pro": ",".join(selected_products),
            }

            if level == "species":
                for pair in result.pairs:
                    reactant_label = self._display_name_for_structure(pair.reactant, lipid_lookup)
                    product_label = self._display_name_for_structure(pair.product, lipid_lookup)
                    pair_key = ",".join(
                        [
                            self._safe_node_id(reactant_label),
                            self._safe_node_id(product_label),
                        ]
                    )
                    details[pair_key] = {
                        "non_selected_re": "",
                        "selected_re": reactant_label,
                        "non_selected_pro": "",
                        "selected_pro": product_label,
                    }

        return details

    def _build_comparison_payloads(
        self,
        disease_group: str,
        control_group: str,
        threshold: float,
        paired: bool,
        dataset: LipidDataset,
        result_set: PathwayReactionSet,
        reaction_lookup: Dict[str, List[ReactionData]],
    ) -> Dict[str, Any]:
        paired_suffix = self._paired_suffix(paired)
        return {
            f"lp_class_reaction_{disease_group}_{control_group}_active_{paired_suffix}.json": self.build_reaction_graph(disease_group, control_group, level="class", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_suppressed_{paired_suffix}.json": self.build_reaction_graph(disease_group, control_group, level="class", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_active_{paired_suffix}.json": self.build_reaction_graph(disease_group, control_group, level="species", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_suppressed_{paired_suffix}.json": self.build_reaction_graph(disease_group, control_group, level="species", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_active_{paired_suffix}.json": self.build_pathway_graph(disease_group, control_group, level="class", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_suppressed_{paired_suffix}.json": self.build_pathway_graph(disease_group, control_group, level="class", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_active_{paired_suffix}.json": self.build_pathway_graph(disease_group, control_group, level="species", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_suppressed_{paired_suffix}.json": self.build_pathway_graph(disease_group, control_group, level="species", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="class", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="class", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="class", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="class", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="species", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="species", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="species", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}.json": self.build_reaction_highlight(disease_group, control_group, threshold, level="species", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="class", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="class", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="class", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="class", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="species", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="species", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="species", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}.json": self.build_pathway_highlight(disease_group, control_group, threshold, level="species", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="class", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="class", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_reaction_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="class", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_class_reaction_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="class", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_species_reaction_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="species", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="species", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_reaction_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="species", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_species_reaction_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_reaction_table(disease_group, control_group, threshold, level="species", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_class_pathway_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="class", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="class", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_class_pathway_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="class", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_class_pathway_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="class", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_species_pathway_{disease_group}_{control_group}_active_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="species", mode="active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="species", mode="suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            f"lp_species_pathway_{disease_group}_{control_group}_most_active_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="species", mode="most_active", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
            f"lp_species_pathway_{disease_group}_{control_group}_most_suppressed_{threshold}_{paired_suffix}_tbl.json": self.build_pathway_table(disease_group, control_group, threshold, level="species", mode="most_suppressed", paired=paired, dataset=dataset, result_set=result_set, reaction_lookup=reaction_lookup, limit=10),
        }

    def _build_edge_detail_payloads(
        self,
        dataset: LipidDataset,
        result_set: PathwayReactionSet,
    ) -> Dict[str, Any]:
        payloads: Dict[str, Any] = {}
        payloads.update({
            f"{name}.json": payload
            for name, payload in self.build_edge_details(level="class", dataset=dataset, result_set=result_set).items()
            if not name.endswith(".json")
        })
        payloads.update({
            f"{name}.json": payload
            for name, payload in self.build_edge_details(level="species", dataset=dataset, result_set=result_set).items()
            if not name.endswith(".json")
        })
        return payloads

    def _cleanup_generated_json_files(self, output_dir: Path) -> None:
        preserved = {"summary.json", "msg1.json", "msg2.json"}
        for file_path in output_dir.glob("*.json"):
            if file_path.name in preserved:
                continue
            file_path.unlink(missing_ok=True)

    def export_reaction_files(
        self,
        output_path: Union[str, Path],
        disease_group: str,
        control_group: str,
        threshold: float = 0.05,
        paired: bool = False,
        dataset: Optional[LipidDataset] = None,
    ) -> Dict[str, str]:
        resolved_dataset = self._get_dataset(dataset)
        output_dir = self._get_output_dir(output_path)
        result_set, reaction_lookup = self.build_reaction_match_set(resolved_dataset)
        written_files: Dict[str, str] = {}

        self._cleanup_generated_json_files(output_dir)

        payloads: Dict[str, Any] = {
            "lp_class_reaction.json": self.build_reaction_tree(level="class", dataset=resolved_dataset, result_set=result_set),
            "lp_species_reaction.json": self.build_reaction_tree(level="species", dataset=resolved_dataset, result_set=result_set),
            "lp_class_pathway.json": self.build_pathway_tree(level="class", dataset=resolved_dataset, result_set=result_set, reaction_lookup=reaction_lookup),
            "lp_species_pathway.json": self.build_pathway_tree(level="species", dataset=resolved_dataset, result_set=result_set, reaction_lookup=reaction_lookup),
        }

        edge_detail_payloads = self._build_edge_detail_payloads(
            dataset=resolved_dataset,
            result_set=result_set,
        )
        comparison_payloads = self._build_comparison_payloads(
            disease_group=disease_group,
            control_group=control_group,
            threshold=threshold,
            paired=paired,
            dataset=resolved_dataset,
            result_set=result_set,
            reaction_lookup=reaction_lookup,
        )

        for file_name, payload in payloads.items():
            file_path = output_dir / file_name
            with file_path.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            written_files[file_name] = str(file_path)

        bundle_path = output_dir / self.COMPARISON_BUNDLE_NAME
        with bundle_path.open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "metadata": {
                        "disease_group": disease_group,
                        "control_group": control_group,
                        "threshold": threshold,
                        "paired": paired,
                    },
                    "payloads": comparison_payloads,
                },
                handle,
                indent=2,
            )
        written_files[self.COMPARISON_BUNDLE_NAME] = str(bundle_path)

        edge_bundle_path = output_dir / self.EDGE_DETAILS_BUNDLE_NAME
        with edge_bundle_path.open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "payloads": edge_detail_payloads,
                },
                handle,
                indent=2,
            )
        written_files[self.EDGE_DETAILS_BUNDLE_NAME] = str(edge_bundle_path)

        return written_files
