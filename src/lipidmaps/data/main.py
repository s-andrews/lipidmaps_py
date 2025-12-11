from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from .data_manager import DataManager
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
        "--use-headgroups",
        dest="use_headgroups",
        action="store_true",
        help="Query LMSD to fill missing LM IDs after RefMet annotation",
    )
    parser.add_argument(
        "--groups",
        nargs="*",
        metavar="GROUP=S1,S2",
        help="Optional group mapping entries (e.g. Control=S1,S2 Treatment=S3,S4)",
    )
    parser.add_argument(
        "--reactions",
        action="store_true",
        help="Query reaction API for LM IDs found in the dataset",
    )
    parser.add_argument(
        "--reaction-api",
        dest="reaction_api",
        default="http://localhost",
        help="Base URL for the reaction API (default: http://localhost)",
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

    # if manager.validation_report and manager.validation_report.passed:
    #     logger.info("Data validation passed successfully.")
    # else:
    #     logger.warning("Data validation failed.")
    #     # exit()

    try:
        dataset = manager.process_csv(csv_path)
    except Exception as exc:
        logger.exception(f"Failed to process CSV: {exc}")
        raise SystemExit(2) from exc

    logger.info(f"Dataset ready: {len(dataset.samples)} samples, {len(dataset.lipids)} lipids")
    logger.info(f"Sample column info: {dataset.samples[:4]}")

    # Optionally fill missing LM IDs using LMSD and report what changed
    if getattr(args, "fill_lmsd", False):
        # Use DataManager helper to run LMSD fill and report updates
        updated_count = manager.run_lmsd_fill_and_report(dataset)
        logger.info(f"Filled {updated_count} missing LM IDs using LMSD")

    # if getattr(args, "use_headgroups", False):
    #     # Use DataManager helper to run LM_ID filling by headgroups and report updates

    group_stats = manager.get_group_statistics()
    logger.info(f"Computed statistics for {len(group_stats)} groups")
    for group_name, stats in group_stats.items():
        logger.info(
            f"{group_name} -> {stats['sample_count']} samples, {stats['lipid_coverage']} lipids"
        )

    # Optionally query reactions for LM IDs present in the dataset
    if getattr(args, "reactions", False):
        # collect all non-empty lm_id values from the dataset
        lm_ids = [q.lm_id for q in dataset.lipids if getattr(q, "lm_id", None)]
        # remove duplicates while preserving order
        seen = set()
        unique_lm = []
        for lid in lm_ids:
            if lid not in seen:
                seen.add(lid)
                unique_lm.append(lid)

        if not unique_lm:
            logger.info("No LM IDs found in dataset to query reactions")
        else:
            logger.info(f"Querying reaction API for {len(unique_lm)} LM IDs")
            reactions = manager.get_reactions_for_lm_ids(unique_lm, base_url=args.reaction_api)
            # Log a compact summary and the full JSON at debug level
            try:
                count = len(reactions.get("reactions", [])) if isinstance(reactions, dict) else 0
            except Exception:
                count = 0
            logger.info(f"Reaction API returned {count} reactions")
            logger.debug("Reaction API full response:\n%s", json.dumps(reactions, indent=2))


if __name__ == "__main__":
    main()
