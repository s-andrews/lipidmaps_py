from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .data import DataManager


logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Regenerate BioPAN session assets using lipidmaps_py"
    )
    parser.add_argument("session_dir", help="BioPAN session directory, e.g. /lipidmaps/temp/biopan/<session_id>")
    parser.add_argument("--csv", dest="csv_path", help="Optional explicit input CSV path. Defaults to <session_dir>/input/input.csv")
    parser.add_argument("--taxonomy-group", default="all", help="Taxonomy group for reaction fetching")
    parser.add_argument("--validate-data", action="store_true", help="Enable CSV validation during import")
    parser.add_argument("--use-refmet", action="store_true", default=True, help="Use RefMet annotation")
    parser.add_argument("--no-use-refmet", dest="use_refmet", action="store_false", help="Disable RefMet annotation")
    parser.add_argument("--use-headgroups", action="store_true", default=True, help="Fill generic LM IDs from headgroups")
    parser.add_argument("--no-use-headgroups", dest="use_headgroups", action="store_false", help="Disable headgroup fill")
    parser.add_argument("--fetch-reactions", action="store_true", default=True, help="Fetch reactions from LIPID MAPS")
    parser.add_argument("--no-fetch-reactions", dest="fetch_reactions", action="store_false", help="Disable reaction fetching")
    parser.add_argument("--paired", action="store_true", help="Mark pathway calculation as paired")
    parser.add_argument("--threshold", type=float, default=0.05, help="Threshold used for highlighted and ranked BioPAN outputs")
    parser.add_argument("--disease-group", help="Condition of interest for reaction asset export")
    parser.add_argument("--control-group", help="Control condition for reaction asset export")
    parser.add_argument(
        "--sample-group",
        action="append",
        default=[],
        metavar="SAMPLE=GROUP",
        help="Override a sample's group assignment before exporting BioPAN assets",
    )
    parser.add_argument("--summary-only", action="store_true", help="Only write summary/msg assets")
    parser.add_argument("--reaction-only", action="store_true", help="Only write reaction assets")
    return parser


def _parse_sample_group_overrides(entries: list[str]) -> dict[str, str]:
    overrides: dict[str, str] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"Invalid --sample-group value '{entry}'. Expected SAMPLE=GROUP")
        sample_name, group_name = entry.split("=", 1)
        sample_name = sample_name.strip()
        group_name = group_name.strip()
        if not sample_name or not group_name:
            raise ValueError(f"Invalid --sample-group value '{entry}'. Expected SAMPLE=GROUP")
        overrides[sample_name] = group_name
    return overrides


def _apply_group_overrides(dataset, overrides: dict[str, str]) -> None:
    if not overrides:
        return

    sample_lookup = {sample.sample_name: sample for sample in dataset.samples}
    missing_samples = sorted(set(overrides) - set(sample_lookup))
    if missing_samples:
        missing_text = ", ".join(missing_samples)
        raise ValueError(f"Unknown sample names in --sample-group: {missing_text}")

    for sample_name, group_name in overrides.items():
        sample_lookup[sample_name].group = group_name


def _resolve_groups(manager: DataManager, disease_group: str | None, control_group: str | None) -> tuple[str, str]:
    if disease_group and control_group:
        return disease_group, control_group

    groups = []
    for sample in manager.dataset.samples if manager.dataset else []:
        if sample.group and sample.group not in groups:
            groups.append(sample.group)

    if len(groups) < 2:
        raise ValueError("At least two sample groups are required to export BioPAN reaction assets")

    return disease_group or groups[1], control_group or groups[0]


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    session_dir = Path(args.session_dir).expanduser().resolve()
    csv_path = Path(args.csv_path).expanduser().resolve() if args.csv_path else session_dir / "input" / "input.csv"
    if not csv_path.exists():
        parser.error(f"Input CSV not found: {csv_path}")

    manager = DataManager(
        validate_data=args.validate_data,
        use_refmet=args.use_refmet,
        use_headgroups=args.use_headgroups,
        fetch_reactions=args.fetch_reactions,
        taxonomy_group=args.taxonomy_group,
    )
    dataset = manager.process_csv(csv_path)
    try:
        group_overrides = _parse_sample_group_overrides(args.sample_group)
        _apply_group_overrides(dataset, group_overrides)
    except ValueError as exc:
        parser.error(str(exc))

    written: dict[str, str] = {}

    if not args.reaction_only:
        written.update(manager.export_biopan_display_files(session_dir, dataset=dataset))

    if not args.summary_only:
        disease_group, control_group = _resolve_groups(manager, args.disease_group, args.control_group)
        written.update(
            manager.export_biopan_reaction_files(
                session_dir,
                disease_group=disease_group,
                control_group=control_group,
                threshold=args.threshold,
                paired=args.paired,
                dataset=dataset,
            )
        )

    print({"written": written})


if __name__ == "__main__":
    main()