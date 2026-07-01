# LIPID MAPS Python Package

Python API for importing, processing, and analyzing lipid data using LIPID MAPS resources.

## Installation

### Development Installation (Editable)

For development, install in editable mode so changes are immediately reflected:

```bash
cd /lipidmaps/lipidmaps_py
pip install -e .
```

### Standard Installation

For production use:

```bash
cd /lipidmaps/lipidmaps_py
pip install .
```

### Installation from Requirements File

If you have a requirements.txt that references this package:

```bash
# In requirements.txt
-e /lipidmaps/lipidmaps_py

# Or for production
/lipidmaps/lipidmaps_py
```

## Usage

After installation, you can import the package as shown in `lipidmaps_reactions_api.py`:

```python
import lipidmaps
from lipidmaps import data

# Import data
lipid_data = lipidmaps.import_data("mydata.csv", lipid_col=1, sample_cols=[4,5,6,7])

# Or import MS-DIAL format
lipid_data = lipidmaps.import_msdial("mydata_msdial.csv")

# Access imported data
print(f"Imported Lipids: {lipid_data.successful_import_count()}")
print(f"Unrecognised Lipids: {lipid_data.failed_import_count()}")

# Fetch reactions for the imported lipids from LIPID MAPS
reactions = lipid_data.dataset.fetch_reactions_by_lm_id(taxonomy_group="human")
print(f"Reactions: {len(reactions)}")
```

## Project Structure

```
lipidmaps_py/
├── setup.py                    # Package configuration
├── src/
│   └── lipidmaps/              # Main package
│       ├── __init__.py         # Top-level API: process_csv, import_data, ...
│       ├── data_importer.py    # import_data / import_msdial -> LipidData
│       ├── biopan_cli.py       # `lipidmaps-biopan` console entry point
│       └── data/               # Core data framework (subpackage)
│           ├── data_manager.py       # DataManager pipeline orchestration
│           ├── quantitation.py       # QuantitationAnalyzer (normalization, stats)
│           ├── biopan_exporter.py            # BioPAN display JSON
│           ├── biopan_pathway_exporter.py    # BioPAN pathway/reaction z-scores
│           ├── ingestion/            # CSV/TSV + metadata readers
│           ├── validation/           # DataValidator
│           ├── matching/             # Species-level reaction matchers
│           ├── models/               # Pydantic v2 data models
│           ├── utils/                # chain_parser, reaction_evaluator, headgroups
│           └── config/               # Bundled JSON (biopan_pathways.json)
└── lipidmaps_reactions_api.py  # API usage examples
```

## Development

The package is structured to allow:
- Easy imports: `import lipidmaps`
- Modular design with separate modules for different functionality
- Integration with LIPID MAPS web APIs

## Dependencies

- pandas
- numpy
- requests

These will be installed automatically when you install the package.
