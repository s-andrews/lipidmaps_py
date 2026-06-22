from __future__ import annotations

import argparse
import logging
from pathlib import Path

from .data import BioPANPathwayExporter, DataManager
from .data.models.sample import LipidDataset
from .logging_utils import configure_logging


logger = logging.getLogger(__name__)

SESSION_DATASET_CACHE_NAME = "processed_dataset.json"


def _resolve_input_path(session_dir: Path, explicit_path: str | None) -> Path:
    if explicit_path:
        return Path(explicit_path).expanduser().resolve()

    supported_paths = [
        session_dir / "input" / "input.csv",
        session_dir / "input" / "input.tsv",
    ]
    for candidate in supported_paths:
        if candidate.exists():
            return candidate

    return supported_paths[0]


def _get_dataset_cache_path(session_dir: Path) -> Path:
    return session_dir / "config" / SESSION_DATASET_CACHE_NAME


def _load_cached_dataset(session_dir: Path) -> LipidDataset | None:
    cache_path = _get_dataset_cache_path(session_dir)
    if not cache_path.exists():
        return None

    try:
        return LipidDataset.model_validate_json(cache_path.read_text(encoding="utf-8"))
    except Exception:
        logger.warning("Failed to load cached BioPAN dataset from %s", cache_path, exc_info=True)
        return None


def _write_cached_dataset(session_dir: Path, dataset: LipidDataset) -> Path:
    cache_path = _get_dataset_cache_path(session_dir)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(dataset.model_dump_json(indent=2), encoding="utf-8")
    return cache_path


def _invalidate_match_set_cache(session_dir: Path) -> None:
    cache_path = session_dir / "config" / BioPANPathwayExporter.MATCH_SET_CACHE_NAME
    try:
        cache_path.unlink(missing_ok=True)
    except Exception:
        logger.warning("Failed to invalidate reaction match set cache at %s", cache_path, exc_info=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Regenerate BioPAN session assets using lipidmaps_py"
    )
    parser.add_argument("session_dir", help="BioPAN session directory, e.g. /lipidmaps/temp/biopan/<session_id>")
    parser.add_argument("--csv", dest="csv_path", help="Optional explicit input table path. Defaults to <session_dir>/input/input.csv or input.tsv")
    parser.add_argument("--taxonomy-group", default="all", help="Taxonomy group for reaction fetching")
    parser.add_argument("--validate-data", action="store_true", help="Enable CSV validation during import")
    parser.add_argument("--has-labels", action="store_true", help="Treat the second row of the CSV as sample labels")
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
    parser.add_argument("--lazy-bundle", action="store_true", help="Write the comparison bundle as a skeleton (metadata + empty payloads); build per-view payloads on demand with --build-view")
    parser.add_argument("--build-view", action="store_true", help="Lazily build only one view's comparison payloads and merge them into the existing bundle")
    parser.add_argument("--scope", choices=["graph", "tables"], default="graph", help="With --build-view: which payloads to build")
    parser.add_argument("--family", choices=["reaction", "pathway"], default="reaction", help="With --build-view: reaction or pathway family")
    parser.add_argument("--level", choices=["class", "species"], default="class", help="With --build-view: class or species level")
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
    try:
        dataset.set_sample_conditions(overrides, strict=True)
    except ValueError as exc:
        message = str(exc)
        if message.startswith("Unknown sample names in condition mapping:"):
            missing_text = message.split(":", 1)[1].strip()
            raise ValueError(f"Unknown sample names in --sample-group: {missing_text}") from exc
        raise


def _resolve_groups(manager: DataManager, disease_group: str | None, control_group: str | None) -> tuple[str, str]:
    if disease_group and control_group:
        return disease_group, control_group

    groups = []
    for sample in manager.dataset.samples if manager.dataset else []:
        if sample.group and sample.group not in groups:
            groups.append(sample.group)

    if len(groups) < 2:
        raise ValueError("At least two sample groups are required to export BioPAN reaction assets")

    # Default order matches the legacy BioPAN tool and the frontend dropdown
    # defaults (graph.js): control = first group, condition of interest
    # (disease) = second group.
    return disease_group or groups[1], control_group or groups[0]


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    log_dir = configure_logging()

    session_dir = Path(args.session_dir).expanduser().resolve()
    csv_path = _resolve_input_path(session_dir, args.csv_path)
    if not csv_path.exists():
        parser.error(f"Input file not found: {csv_path}")

    logger.info("BioPAN CLI logging to %s", log_dir)
    logger.info("Starting BioPAN export for session %s", session_dir)

    try:
        manager = DataManager(
            validate_data=args.validate_data,
            has_labels=args.has_labels,
            use_refmet=args.use_refmet,
            use_headgroups=args.use_headgroups,
            fetch_reactions=args.fetch_reactions,
            taxonomy_group=args.taxonomy_group,
        )
        dataset = None
        if args.csv_path is None:
            dataset = _load_cached_dataset(session_dir)
            if dataset is not None:
                logger.info("Loaded cached BioPAN dataset from %s", _get_dataset_cache_path(session_dir))

        if dataset is None:
            dataset = manager.process_csv(csv_path)
            cache_path = _write_cached_dataset(session_dir, dataset)
            logger.info("Cached processed BioPAN dataset at %s", cache_path)
            # The reaction match set is derived from the dataset, so a freshly
            # processed dataset invalidates any cached match set.
            _invalidate_match_set_cache(session_dir)

        manager.dataset = dataset
        group_overrides = _parse_sample_group_overrides(args.sample_group)
        _apply_group_overrides(dataset, group_overrides)
    except ValueError as exc:
        parser.error(str(exc))
    except Exception:
        logger.exception("BioPAN export failed for session %s", session_dir)
        raise

    written: dict[str, str] = {}

    # Lazy per-view build: only build the requested view's comparison payloads
    # and merge them into the existing bundle. Reuses the cached match set.
    if args.build_view:
        try:
            disease_group, control_group = _resolve_groups(manager, args.disease_group, args.control_group)
            exporter = manager.get_biopan_pathway_exporter(dataset)
            built = exporter.build_and_merge_view(
                session_dir,
                disease_group,
                control_group,
                args.threshold,
                args.paired,
                scope=args.scope,
                family=args.family,
                level=args.level,
                dataset=dataset,
            )
        except ValueError as exc:
            parser.error(str(exc))
        logger.info("BioPAN built %s lazy view payloads for session %s", len(built), session_dir)
        print({"built_view": sorted(built.keys())})
        return

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
                lazy=args.lazy_bundle,
            )
        )

    logger.info("BioPAN export completed for session %s", session_dir)
    logger.info("Wrote %s assets", len(written))

    print({"written": written})


if __name__ == "__main__":
    main()