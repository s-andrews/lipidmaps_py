from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

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
    base = BioPANPathwayExporter.MATCH_SET_CACHE_NAME
    # Both the default and the legacy-pairing cache variants must be cleared.
    for name in (base, base.replace(".json", "_legacy.json")):
        cache_path = session_dir / "config" / name
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
    parser.add_argument("--legacy-substrate-consumption", action="store_true", help="Reproduce legacy BioPAN greedy substrate-consumption pairing (matches old-tool z-scores)")
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


@dataclass
class RunParams:
    """A single BioPAN export request. Mirrors the CLI arguments so that the CLI
    and the warm service (`biopan_service`) drive the exact same code path."""

    session_dir: Path
    csv_path: Optional[str] = None
    taxonomy_group: str = "all"
    validate_data: bool = False
    has_labels: bool = False
    use_refmet: bool = True
    use_headgroups: bool = True
    fetch_reactions: bool = True
    paired: bool = False
    legacy_substrate_consumption: bool = False
    threshold: float = 0.05
    disease_group: Optional[str] = None
    control_group: Optional[str] = None
    sample_group: list[str] = field(default_factory=list)
    summary_only: bool = False
    reaction_only: bool = False
    lazy_bundle: bool = False
    build_view: bool = False
    scope: str = "graph"
    family: str = "reaction"
    level: str = "class"

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "RunParams":
        return cls(
            session_dir=Path(args.session_dir).expanduser().resolve(),
            csv_path=args.csv_path,
            taxonomy_group=args.taxonomy_group,
            validate_data=args.validate_data,
            has_labels=args.has_labels,
            use_refmet=args.use_refmet,
            use_headgroups=args.use_headgroups,
            fetch_reactions=args.fetch_reactions,
            paired=args.paired,
            legacy_substrate_consumption=args.legacy_substrate_consumption,
            threshold=args.threshold,
            disease_group=args.disease_group,
            control_group=args.control_group,
            sample_group=list(args.sample_group or []),
            summary_only=args.summary_only,
            reaction_only=args.reaction_only,
            lazy_bundle=args.lazy_bundle,
            build_view=args.build_view,
            scope=args.scope,
            family=args.family,
            level=args.level,
        )


@dataclass
class RunResult:
    """Outcome of `run_session`. `exporter` and `dataset` are returned so a warm
    caller (the service registry) can keep them between requests; `reprocessed`
    signals that the dataset was rebuilt from CSV and any warm cache is stale."""

    dataset: LipidDataset
    exporter: Optional[BioPANPathwayExporter]
    reprocessed: bool
    output: dict[str, Any]


def run_session(
    params: RunParams,
    *,
    dataset: Optional[LipidDataset] = None,
    exporter: Optional[BioPANPathwayExporter] = None,
) -> RunResult:
    """Run one BioPAN export request.

    This is the shared core for both the CLI (`main`) and the warm service. A
    warm caller may inject an already-loaded ``dataset`` and ``exporter`` (whose
    in-memory match set / table caches survive between requests); when omitted
    the dataset is loaded from the on-disk cache or reprocessed from CSV exactly
    as the standalone CLI does. ``ValueError`` is raised for user-facing input
    errors (missing file, bad groups) so callers can map it to an exit / 400.
    """
    session_dir = params.session_dir
    csv_path = _resolve_input_path(session_dir, params.csv_path)
    if not csv_path.exists():
        raise ValueError(f"Input file not found: {csv_path}")

    logger.info("Starting BioPAN export for session %s", session_dir)

    manager = DataManager(
        validate_data=params.validate_data,
        has_labels=params.has_labels,
        use_refmet=params.use_refmet,
        use_headgroups=params.use_headgroups,
        fetch_reactions=params.fetch_reactions,
        taxonomy_group=params.taxonomy_group,
        legacy_substrate_consumption=params.legacy_substrate_consumption,
    )

    reprocessed = False
    if dataset is None and params.csv_path is None:
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
        reprocessed = True
        # A reprocess produces a different dataset object, so any warm exporter
        # bound to the previous dataset must not be reused.
        exporter = None

    manager.dataset = dataset
    # A warm exporter built with a different legacy-pairing setting must not be
    # reused, or a toggled comparison would serve stale (other-mode) z-scores.
    if exporter is not None and getattr(
        exporter, "legacy_substrate_consumption", False
    ) != params.legacy_substrate_consumption:
        exporter = None
    _apply_group_overrides(dataset, _parse_sample_group_overrides(params.sample_group))

    # Lazy per-view build: only build the requested view's comparison payloads
    # and merge them into the existing bundle. Reuses the cached match set.
    if params.build_view:
        disease_group, control_group = _resolve_groups(manager, params.disease_group, params.control_group)
        if exporter is None:
            exporter = manager.get_biopan_pathway_exporter(dataset)
        built = exporter.build_and_merge_view(
            session_dir,
            disease_group,
            control_group,
            params.threshold,
            params.paired,
            scope=params.scope,
            family=params.family,
            level=params.level,
            dataset=dataset,
        )
        logger.info("BioPAN built %s lazy view payloads for session %s", len(built), session_dir)
        return RunResult(dataset=dataset, exporter=exporter, reprocessed=reprocessed, output={"built_view": built})

    written: dict[str, str] = {}

    if not params.reaction_only:
        written.update(manager.export_biopan_display_files(session_dir, dataset=dataset))

    if not params.summary_only:
        disease_group, control_group = _resolve_groups(manager, params.disease_group, params.control_group)
        written.update(
            manager.export_biopan_reaction_files(
                session_dir,
                disease_group=disease_group,
                control_group=control_group,
                threshold=params.threshold,
                paired=params.paired,
                dataset=dataset,
                lazy=params.lazy_bundle,
            )
        )

    logger.info("BioPAN export completed for session %s", session_dir)
    logger.info("Wrote %s assets", len(written))

    return RunResult(dataset=dataset, exporter=exporter, reprocessed=reprocessed, output={"written": written})


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    log_dir = configure_logging()
    logger.info("BioPAN CLI logging to %s", log_dir)

    params = RunParams.from_args(args)

    try:
        result = run_session(params)
    except ValueError as exc:
        parser.error(str(exc))
    except Exception:
        logger.exception("BioPAN export failed for session %s", params.session_dir)
        raise

    if params.build_view:
        print({"built_view": sorted(result.output.get("built_view", {}).keys())})
    else:
        print({"written": result.output.get("written", {})})


if __name__ == "__main__":
    main()