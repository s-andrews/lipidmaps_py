from __future__ import annotations

import argparse
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
            logger.warning("Skipping malformed group mapping '%s'", entry)
            continue
        name, samples_str = entry.split("=", 1)
        samples = [s.strip() for s in samples_str.split(",") if s.strip()]
        if not samples:
            logger.warning("Group '%s' has no samples; skipping", name)
            continue
        mapping[name.strip()] = samples
    return mapping or None


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    csv_path = Path(args.csv).expanduser()
    if not csv_path.exists():
        parser.error(f"CSV file not found: {csv_path}")

    logger.info("Processing CSV: %s", csv_path)

    group_mapping = parse_group_mapping(args.groups)
    manager = DataManager(
        validate_data=args.validate,
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
        logger.exception("Failed to process CSV: %s", exc)
        raise SystemExit(2) from exc

    logger.info(
        "Dataset ready: %s samples, %s lipids",
        len(dataset.samples),
        len(dataset.lipids),
    )

    group_stats = manager.get_group_statistics()
    logger.info("Computed statistics for %s groups", len(group_stats))
    for group_name, stats in group_stats.items():
        logger.info(
            "%s -> %s samples, %s lipids",
            group_name,
            stats["sample_count"],
            stats["lipid_coverage"],
        )


if __name__ == "__main__":
    main()
