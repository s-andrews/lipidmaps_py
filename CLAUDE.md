# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`lipidmaps_py` ingests, normalizes, validates, and analyzes mass-spectrometry lipidomics datasets, and interfaces with LIPID MAPS web resources (RefMet standardization, LM ID lookup, reactions). It is a `src/`-layout package (`src/lipidmaps`) installed via `setup.py`. A key consumer is the downstream BioPAN tool at `/lipidmaps/lipidmaps/www/html/resources/tools/biopan`, which imports this package (top-level helpers *and* some internal modules) while its legacy R processing is retired — changes here can break it even though it lives outside this repo.

## Commands

```bash
pip install -e .[dev]          # dev install (pydantic, pandas, numpy, requests, scipy + pytest/black/flake8/mypy/streamlit/plotly)
pytest                         # full suite; pytest.ini adds importlib mode, coverage (htmlcov/), and report.html
pytest tests/data/test_csv_ingestion.py::TestName::test_case   # single test
PYTEST_ADDOPTS="" pytest -q     # fast loop, no coverage/html reporting add-ons
flake8                         # config in .flake8 (max-line-length 100; E203/E501/W503 ignored)
mkdocs serve                   # preview API docs (mkdocstrings-generated from docstrings)
```

- Always run tests through `pytest`, never by ad-hoc `python -m` execution — the suite depends on the configured `--import-mode=importlib`.
- SQLite is required for tests. `ModuleNotFoundError: No module named '_sqlite3'` means Python was built without SQLite; use a system Python or rebuild with `sqlite-devel`.

## Architecture

**User-facing entry points** live at the top level:
- `process_csv(path, **kwargs)` (`src/lipidmaps/__init__.py`) — one-call CSV → `LipidDataset`; forwards all `DataManager` kwargs (`validate_data`, `use_refmet`, `use_headgroups`, `taxonomy_group`, `fetch_reactions`, ...).
- `import_data()` / `import_msdial()` (`data_importer.py`) — lower-level importers returning `LipidData`.

**Processing pipeline** — `DataManager.process_csv()` (`data/data_manager.py`) orchestrates: ingest CSV (`data/ingestion/`) → optional validation (`data/validation/data_validator.py`) → RefMet name standardization + LM ID lookup (`data/models/refmet.py`) → optional headgroup-based generic LM ID fill (`data/utils/headgroups.py`) → optional reaction fetch/annotation. Output is a `LipidDataset`.

**Data models** (`data/models/`, Pydantic v2) are the typed spine used everywhere:
- `sample.py` — `LipidDataset` (the central container: `samples`, `lipids`, plus query/aggregation/normalization/reaction-annotation helpers), `QuantifiedLipid`, `SampleMetadata`, `SampleConditions`, `Quantitation`, `LipidAnnotation`. Reaction fetching lives here as `LipidDataset.fetch_reactions_by_lm_id()`, not on `DataManager`.
- `reaction.py` / `species_reaction.py` — `ReactionData`, `ReactionChecker`, and species-level `ClassReaction`/`ReactionType`/`CompoundRequirement` used by matching.
- `refmet.py`, `lmsd.py`, `uniprot.py` — external-API clients.

**Reaction matching** (active area — `reaction_matching` branch):
- `data/matching/` — Strategy pattern. `registry.py` (`MatcherRegistry`) selects a matcher by `ReactionType`; matchers (`same_structure`, `fa_compound`, `facoa_compound`, `sphingolipid`, `cardiolipin`) inherit `base.ReactionMatcher`. `use_full_structure=True` swaps in structure-aware FA/FACoA matchers.
- `data/utils/reaction_evaluator.py` (`ReactionEvaluator`) — evaluates whether a class reaction is possible for given species, driven by `lipid_reaction_rules.py`; used by `fetch_reactions_by_lm_id(annotate_using_evaluator=True)`.
- `data/utils/chain_parser.py` — parses lipid names into `LipidStructure`/`AcylChain`; foundational to matching.

**BioPAN integration**:
- `biopan_cli.py` → `lipidmaps-biopan` console script (entry point in `setup.py`); orchestration in `data/main.py`. Regenerates BioPAN session assets (`msg1.json`, `summary.json`, `msg2.json`, reaction/pathway graph + table payloads) into a session dir. Supports group comparison (`--disease-group`/`--control-group`/`--threshold`), per-sample group overrides (`--sample-group A=control`), lazy per-view building (`--lazy-bundle`/`--build-view`), and `--legacy-substrate-consumption` for parity with the old R tool's z-scores.
- `data/biopan_exporter.py` (`BioPANExporter`) and `data/biopan_pathway_exporter.py` (`BioPANPathwayExporter`) write the JSON the BioPAN PHP frontend expects. `DataManager.export_biopan_display_files()` is a compatibility wrapper.

**External endpoints** default to `https://dev.lipidmaps.org` (RefMet, reactions, LMSD, molecules). The reactions base URL is overridable via the `LMSD_REACTIONS_BASE_URL` env var (`src/lipidmaps/config.py`).

**Demo/scripts** (`scripts/`) — Streamlit demo (`app.py` → `scripts/streamlit_demo.py`, deployed via `ecosystem.config.js`/pm2 on port 8501), plus headless eval and live-comparison scripts. These are not part of the installed package API.

## Conventions

- **Pydantic v2 only.** Use `model_config = ConfigDict(...)` / dict `model_config`, `field_validator`, `computed_field`, `model_post_init`, `model_dump()`. Never introduce v1 idioms (`class Config`, `.dict()`). Most models subclass `LipidmapsBaseModel` (`data/models/base.py`), which auto-generates UUIDs — tests must not depend on fixed object IDs.
- Prefer top-level APIs (`process_csv`, `import_data`) for new call sites; keep them easy for downstream callers.
- Preserve object-based access patterns (sample↔lipid helpers on `LipidDataset`) when editing models.
- Be conservative renaming/moving anything under `src/lipidmaps/data/` — downstream BioPAN scripts import internal modules directly. Preserve the `lipidmaps-biopan` CLI surface and export flows.
- Docs to update rather than duplicate: `README.md` (install/usage), `INSTALL.md` (setup), `docs/api.md`, `docs/models.md`, `docs/custom_columns_guide.md`.
