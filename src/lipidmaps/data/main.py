from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .data_manager import DataManager
from .models.lmsd import LMSD
from .ingestion.csv_reader import CSVFormat


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("lipidmaps_py.log"), logging.StreamHandler()],
)

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

    csv_path = Path(args.csv).expanduser()
    if not csv_path.exists():
        parser.error(f"CSV file not found: {csv_path}")

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
    print(reactions[:1])  # Print first 1 reaction for brevity
    logger.info(f"Lipids with reactions: {dataset.list_lipids_with_reactions()[:1]}")
    
    lipid_objects_with_reactions = dataset.get_lipids_with_reactions()
    logger.info(f"Total lipids with reactions: {len(lipid_objects_with_reactions)}")
    logger.info(f"Total reaction object in dataset: {len(dataset.reactions)}")
    for reaction_name in dataset.list_reactions():  # Print first 1 reaction name
        logger.info(f"{reaction_name}")
    for reaction in dataset.reactions[:3]:  # Print first 3 lipid with reactions
        print(f"Reaction name: {reaction.reaction_name} {reaction.reaction_id}\n"
              f"Reactant names: {[(r.compound_name, r.compound_lm_id) for r in (reaction.reactants or [])]}\n"

              f"Product names: {[(p.compound_name, p.compound_lm_id) for p in (reaction.products or [])]}")

if __name__ == "__main__":
    main()
