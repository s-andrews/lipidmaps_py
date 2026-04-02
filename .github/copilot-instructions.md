# Project Guidelines

## Code Style

- This repository is a Python package under `src/lipidmaps` built with `setup.py`; keep changes consistent with the existing package layout and import style.
- Use the shared virtual environment at `/lipidmaps/lipidmaps/www/venv` when running Python commands for this repo.
- Treat Pydantic as v2-only. Follow existing patterns such as `model_config = ConfigDict(...)` or dict-based `model_config`, `field_validator`, `computed_field`, `model_post_init`, and `model_dump()`.
- Keep existing model patterns in place rather than rewriting them to a different style. Representative files: `src/lipidmaps/data/models/base.py`, `src/lipidmaps/data/models/sample.py`, and `src/lipidmaps/data_importer.py`.

## Architecture

- `src/lipidmaps/data_importer.py` and `src/lipidmaps/__init__.py` expose the main user-facing entry points such as `process_csv()` and `import_data()`.
- `src/lipidmaps/data/` contains the core ingestion, validation, quantitation, export, and reaction logic.
- `src/lipidmaps/data/models/` contains the typed data layer used across the package.
- `src/lipidmaps/biopan_cli.py` is the compatibility CLI used to generate BioPAN session assets.
- An adjacent downstream project at `/lipidmaps/lipidmaps/www/html/resources/tools/biopan` consumes this package while R-based processing is being removed. Changes here can affect that tool even though it lives outside this workspace.

## Build And Test

- Install for development with `pip install -e .[dev]`.
- Run tests with `pytest`. The repo config already enables `--import-mode=importlib`, coverage output, `htmlcov/`, and `report.html`.
- Use targeted pytest runs when changing a narrow area, then run broader tests if the change affects shared models, ingestion, exporters, or BioPAN integration.
- Build docs with `mkdocs build` or preview with `mkdocs serve` when changing documentation-oriented code.

## Conventions

- Prefer the high-level package APIs unless a task specifically requires internal plumbing. Keep `process_csv()` and related top-level flows easy for downstream callers to use.
- Preserve compatibility for BioPAN-facing outputs and entry points, especially `lipidmaps-biopan`, dataset export flows, and modules imported by external scripts.
- When editing models or dataset access patterns, preserve object-based access conventions such as sample-to-lipid and lipid-to-sample helpers.
- Do not introduce Pydantic v1 idioms such as `class Config` or `.dict()` in new code.
- Link to existing docs instead of duplicating them: `README.md` for installation and usage, `INSTALL.md` for setup, `docs/api.md` for API coverage, `docs/models.md` for model references, and `WORKING_NOTES.txt` for design intent and active direction.

## Pitfalls

- Some downstream Biopan scripts import internal modules directly, not just top-level helpers. Be careful when changing names or locations under `src/lipidmaps/data/`.
- Tests should be run through `pytest`, not by ad hoc direct module execution, because the repo relies on the configured import mode and reporting defaults.
- Base models generate UUIDs automatically; tests should not depend on fixed object IDs.