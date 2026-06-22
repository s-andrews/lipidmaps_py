from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .data_manager import DataManager
from .models.lmsd import LMSD
from .ingestion.csv_reader import CSVFormat
from .matching import match_pathway_reactions, create_matcher_context
from .models.species_reaction import ClassReaction, CompoundRequirement
from .utils.chain_parser import (
    get_common_fa_names,
    get_common_facoa_names,
    infer_fa_from_lipids,
    infer_facoa_from_lipids,
    extract_fa_from_reactions,
    extract_facoa_from_reactions,
)
from ..logging_utils import configure_logging

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:

    """Construct CLI argument parser for quick package smoke testing."""
    parser = argparse.ArgumentParser(
        description="Process a lipidomics CSV using the lipidmaps_py DataManager"
    )
    parser.add_argument(
        "csv",
        metavar="CSV_PATH",
        help="Path to the quantified lipid CSV file",
    )
    parser.add_argument(
        "--format",
        choices=[fmt.value for fmt in CSVFormat],
        default=CSVFormat.AUTO.value,
        help="Optional CSV format override (default: auto-detect)",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Run DataValidator checks during processing",
    )
    parser.add_argument(
        "--has-labels",
        action="store_true",
        help="Indicate if the second row in the CSV contains labels",
    )
    parser.add_argument(
        "--fill-lmsd",
        dest="fill_lmsd",
        action="store_true",
        help="Query LMSD to fill missing LM IDs after RefMet annotation",
    )
    parser.add_argument(
        "--fill-headgroups",
        dest="fill_headgroups",
        action="store_true",
        help="Fill missing LM IDs using headgroup mapping after LMSD fill",
    )
    parser.add_argument(
        "--groups",
        nargs="*",
        metavar="GROUP=S1,S2",
        help="Optional group mapping entries (e.g. Control=S1,S2 Treatment=S3,S4)",
    )
    parser.add_argument(
        "--test-matching",
        dest="test_matching",
        action="store_true",
        help="Test species-level reaction matching on the dataset",
    )
    return parser


def parse_group_mapping(raw_groups: list[str] | None) -> dict[str, list[str]] | None:
    """Convert GROUP=S1,S2 CLI inputs into mapping dict."""
    if not raw_groups:
        return None
    mapping: dict[str, list[str]] = {}
    for entry in raw_groups:
        if "=" not in entry:
            logger.warning(f"Skipping malformed group mapping '{entry}'")
            continue
        name, samples_str = entry.split("=", 1)
        samples = [s.strip() for s in samples_str.split(",") if s.strip()]
        if not samples:
            logger.warning(f"Group '{name}' has no samples; skipping")
            continue
        mapping[name.strip()] = samples
    return mapping or None


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    log_dir = configure_logging()

    csv_path = Path(args.csv).expanduser()
    if not csv_path.exists():
        parser.error(f"CSV file not found: {csv_path}")

    logger.info(f"CLI logging to: {log_dir}")
    logger.info(f"Processing CSV: {csv_path}")

    group_mapping = parse_group_mapping(args.groups)
    manager = DataManager(
        validate_data=args.validate,
        has_labels=args.has_labels,
        csv_format=CSVFormat(args.format),
        group_mapping=group_mapping,
    )

    # (reaction annotation will be performed after CSV processing)

    try:
        dataset = manager.process_csv(csv_path)
    except Exception as exc:
        logger.exception(f"Failed to process CSV: {exc}")
        raise SystemExit(2) from exc

    # Dataset is now ready; perform reaction fetching/annotation and display selected lipids
    logger.info(f"Dataset ready: {len(dataset.samples)} samples, {len(dataset.lipids)} lipids")
    logger.info(f"Sample column info: {dataset.samples[:3]}")
    logger.info(f"Samples: {dataset.list_sample_names()[:3]}")

    # Get structures grouped by headgroup from the dataset
    try:
        structures_by_headgroup = dataset.get_structures_by_headgroup()
        logger.info(f"Found {len(structures_by_headgroup)} unique headgroups from dataset")
        logger.info(f"Headgroup names: {list(structures_by_headgroup.keys())[:5]}")
        for hg_name, structures in list(structures_by_headgroup.items())[:10]:
            structure_names = [s.full_name for s in structures]
            print(f"Headgroup: {hg_name} | Structures ({len(structures)}): {structure_names[:5]}")
    except Exception:
        logger.exception("Failed to get structures by headgroup")

    lmids = dataset.list_lm_ids()
    print(f"LM IDs in dataset: {LMSD.get_molecules_by_lm_id(lmids)}")  # Print first 5 LM IDs
    # Optionally fill missing LM IDs using LMSD and report what changed
    if getattr(args, "fill_lmsd", False):
        # Use DataManager helper to run LMSD fill and report updates
        updated_count = manager.run_lmsd_fill_and_report(dataset)
        logger.info(f"Filled {updated_count} missing LM IDs using LMSD")

    # Optionally fill missing LM IDs using headgroup mapping
    if getattr(args, "fill_headgroups", False):
        updated_count = manager.fill_generic_lm_ids_from_headgroups(dataset)
        logger.info(f"Filled {updated_count} missing LM IDs using headgroup mapping")
        
    group_stats = manager.get_group_statistics()
    logger.info(f"Computed statistics for {len(group_stats)} groups")
    for group_name, stats in group_stats.items():
        logger.info(
            f"{group_name} -> {stats['sample_count']} samples, {stats['lipid_coverage']} lipids"
        )

    # Fetch reactions for all LM IDs in the dataset
    reactions = dataset.fetch_reactions_by_lm_id(reaction_type="class-level", only_lipid_components=True, taxonomy_group="bacteria")
    # print(reactions[:1])  # Print first 1 reaction for brevity
    logger.info(f"Lipids with reactions: {dataset.list_lipids_with_reactions()[:1]}")
    
    lipid_objects_with_reactions = dataset.get_lipids_with_reactions()
    logger.info(f"Total lipids with reactions: {len(lipid_objects_with_reactions)}")
    logger.info(f"Total reaction object in dataset: {len(dataset.reactions)}")
    # for reaction_name in dataset.list_reactions():  # Print first 1 reaction name
    #     logger.info(f"{reaction_name}")
    # for reaction in dataset.reactions[:3]:  # Print first 3 lipid with reactions
    #     print(f"Reaction name: {reaction.reaction_name} {reaction.reaction_id}\n"
    #           f"Reactant names: {[(r.compound_name, r.compound_lm_id) for r in (reaction.reactants or [])]}\n"

    #           f"Product names: {[(p.compound_name, p.compound_lm_id) for p in (reaction.products or [])]}")

    # Test species-level reaction matching
    if getattr(args, "test_matching", False):
        test_species_matching(dataset)


def test_species_matching(dataset) -> None:
    """Test species-level reaction matching on the dataset."""
    logger.info("=" * 60)
    logger.info("Testing species-level reaction matching")
    logger.info("=" * 60)
    
    # Get all lipid names from dataset
    lipid_names = [lipid.standardized_name or lipid.input_name or None for lipid in dataset.lipids]
    logger.info(f"Dataset has {len(lipid_names)} lipids")
    
    # --- FA / FACoA pool construction ---
    # Union the available sources rather than taking the first non-empty one:
    # reaction extraction often yields only a sparse, non-representative set
    # (e.g. just FA 16:0 + 2:0 from LM IDs), which would otherwise mask the
    # richer set of acyl chains inferable from the dataset's own lipids and
    # starve FA-release / FACoA-addition matching. Only fall back to common
    # defaults when nothing at all can be sourced.

    def _ordered_union(*groups: List[str]) -> List[str]:
        seen = set()
        ordered: List[str] = []
        for group in groups:
            for name in group:
                if name and name not in seen:
                    seen.add(name)
                    ordered.append(name)
        return ordered

    fa_from_reactions = extract_fa_from_reactions(dataset.reactions) if dataset.reactions else []
    fa_from_names = [name for name in lipid_names if name.startswith("FA ") or name.startswith("FA(")]
    facoa_from_names = [name for name in lipid_names
                        if name.startswith("CoA ") or name.startswith("FACoA") or name.startswith("FaCoA")]
    # Gate mono-acyl (lyso) self-inference on the dataset having no measured
    # FA/FaCoA, mirroring legacy BioPAN (parse_data.r): lyso/CE/sphingo chains
    # only seed the pool when the input carries no fatty acids of its own, so a
    # lyso product (e.g. LPS in PS->LPS) does not validate its own release.
    infer_mono_acyl = not (fa_from_names or facoa_from_names)
    fa_inferred = infer_fa_from_lipids(lipid_names, include_mono_acyl=infer_mono_acyl)
    fa_names = _ordered_union(fa_from_reactions, fa_from_names, fa_inferred)
    if fa_names:
        logger.info(
            f"FA pool: {len(fa_names)} species "
            f"(reactions={len(fa_from_reactions)}, dataset={len(fa_from_names)}, inferred={len(fa_inferred)})"
        )
    else:
        fa_names = get_common_fa_names()
        logger.info(f"Using {len(fa_names)} common FA species (no FA in dataset)")

    facoa_from_reactions = extract_facoa_from_reactions(dataset.reactions) if dataset.reactions else []
    facoa_inferred = infer_facoa_from_lipids(lipid_names, include_mono_acyl=infer_mono_acyl)
    facoa_names = _ordered_union(facoa_from_reactions, facoa_from_names, facoa_inferred)
    if facoa_names:
        logger.info(
            f"CoA pool: {len(facoa_names)} species "
            f"(reactions={len(facoa_from_reactions)}, dataset={len(facoa_from_names)}, inferred={len(facoa_inferred)})"
        )
    else:
        facoa_names = get_common_facoa_names()
        logger.info(f"Using {len(facoa_names)} common CoA species (no CoA in dataset)")
    
    # Define some common class reactions to test
    test_reactions = [
        ClassReaction(reactant_class="PC", product_class="PA"),
        ClassReaction(reactant_class="PE", product_class="PC"),
        ClassReaction(reactant_class="PA", product_class="DG"),
        ClassReaction(reactant_class="DG", product_class="TG", compound_require=CompoundRequirement.FACOA, acyl_add=True),
        ClassReaction(reactant_class="PC", product_class="LPC", compound_require=CompoundRequirement.FA),
        ClassReaction(reactant_class="Cer", product_class="SM"),
        ClassReaction(reactant_class="PG", product_class="CL"),
    ]
    
    # Run matching
    results = match_pathway_reactions(
        lipid_names=lipid_names,
        reactions=test_reactions,
        fa_names=fa_names,
        facoa_names=facoa_names,
    )
    
    # Report results
    logger.info(f"\nMatching Summary:")
    logger.info(f"  Reactions checked: {results.reactions_checked}")
    logger.info(f"  Reactions with matches: {results.reactions_with_pairs}")
    total_pairs = sum(len(r.pairs) for r in results.results.values())
    logger.info(f"  Total species pairs: {total_pairs}")
    
    # Show details for each reaction
    for reaction_key, result in results.results.items():
        if result.has_pairs:
            logger.info(f"\n{result.class_reaction.reactant_class} -> {result.class_reaction.product_class}:")
            logger.info(f"  Reactant species: {result.reactant_species_found}")
            logger.info(f"  Product species: {result.product_species_found}")
            logger.info(f"  Valid pairs: {len(result.pairs)}")
            
            # Show first few pairs
            for pair in result.pairs[:5]:
                compound_info = ""
                if pair.required_compound:
                    compound_info = f" (via {pair.required_compound.notation})"
                logger.info(f"    {pair.reactant.full_name} -> {pair.product.full_name}{compound_info}")
            
            if len(result.pairs) > 5:
                logger.info(f"    ... and {len(result.pairs) - 5} more pairs")


if __name__ == "__main__":
    main()
