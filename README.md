# lipidmaps_py

A Python package providing tools to ingest, normalize, validate, and manage lipidomics datasets and to interface with LIPID MAPS resources.

This project is intended for researchers and developers working with mass-spectrometry lipidomics data who need reproducible preprocessing (ingestion, normalization, and QC), integration with RefMet identifiers, and programmatic access to dataset management utilities.

## Purpose

- Provide robust CSV/TSV ingestion with flexible column handling and format detection.
- Normalize lipid names to RefMet where possible so downstream analyses work with standardized identifiers.
- Validate datasets and generate concise QC reports highlighting missing values, format inconsistencies, and common data issues.
- Offer a `DataManager` abstraction for working with quantified lipids, samples, and simple cohort metadata.
- Lay the groundwork for LIPID MAPS API integration (LM ID lookup) and reaction-analysis features.

## Development Status

### ✅ Complete
- **Data Import & Validation**: CSV/TSV data ingestion with format detection
- **Name Standardization**: RefMet standardization of user provided metabolite names
- **Quality Control**: Data validation and issue reporting
- **Data Management**: DataManager for handling quantified lipid datasets

### 🚧 In Progress
- **LIPID MAPS API Integration**: LM ID lookup and validation
- **Reaction Analysis**: Integration with LIPID MAPS reactions database
- **Sample Metadata**: Support for experimental metadata and conditions

## Installation

### Prerequisites
- Python 3.9 or higher

> **Note**: Test functions require sqlite3 support in python. If you encounter `ModuleNotFoundError: No module named '_sqlite3'`, your Python installation was built without SQLite support. Either:
> - Use your system's Python (e.g., `python3` instead of a custom-built Python)
> - Rebuild Python with SQLite development libraries installed (`sudo dnf install sqlite-devel` on AlmaLinux/Fedora/RHEL, then rebuild Python)

### Install from Source

1. Clone the repository:
```bash
git clone https://github.com/s-andrews/lipidmaps_py.git
cd lipidmaps_py
```

2. Create a virtual environment (recommended):
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package:
```bash
pip install .
```

or to include dev requirements
```bash
pip install -e .[dev]
```

4. Install development dependencies (optional and if you haven't already used .[dev], for running tests):
```bash
pip install pytest pytest-cov pytest-html black flake8 mypy
```

### Verify Installation

Test that the package is correctly installed:
```bash
python -c "import lipidmaps; print('Installation successful')"
```

### API Documentation

By using mkdocs package, you can view API documentation locally
```bash
pip install mkdocs mkdocs-material mkdocstrings[python]
mkdocs build
mkdocs serve
```

## Running Tests & Reports

The package ships with a comprehensive pytest configuration (`pytest.ini`) that automatically:

- Runs tests in `importlib` mode (fixes relative-import issues)
- Collects coverage for the `lipidmaps` package
- Writes a terminal coverage summary and an HTML coverage site (`htmlcov/index.html`)
- Generates a standalone HTML test report (`report.html`) you can archive or share

Running the full suite is therefore as simple as:

```bash
# Inside the repo (venv recommended)
pytest
```

### Targeted Test Runs

```bash
# Verbose output
pytest -v

# Specific directory or file
pytest tests/data/
pytest tests/data/test_csv_ingestion.py

# Disable reporting add-ons if you need a quicker loop
PYTEST_ADDOPTS="" pytest -q
```

### Viewing Reports

- **Test results**: open `report.html` in any browser for per-test details, logs, and attachments.
- **Coverage**: open `htmlcov/index.html` for annotated source along with percentage metrics.
- Both artifacts are produced on every `pytest` run locally and in CI (uploaded as workflow artifacts).

## Quick Start

### Basic Usage

```python
from lipidmaps import process_csv

# Load a CSV file. The package includes sample datasets in the `tests/data/inputs/` directory:
dataset = process_csv("path/to/your/data.csv")
# process_csv accepts many parameters that DataManager class accepts.
dataset = process_csv("path/to/your/data.csv", validate_data=True, use_refmet=True, use_headgroups=True, taxonomy_group="mammalia")

# Csv file is processed into an object with iterable samples and lipids data
# samples - the list of SampleMetadata type objects with sample_name, group and label attributes
print(dataset.samples[:1]) 

# lipids - the list of QuantifiedLipid type objects with input_name, standardized_name, lm_id, recognized and values object "sample_name": "value"
print(dataset.lipids[:1]) 

# List first 5 sample names 
print(f"Samples: {dataset.list_sample_names()[:5]}")

# List first 5 lipid names
print(f"Lipids: {dataset.list_lipid_names()[:5]}")

```

## Refmet Name Standardization
- [Refmet](https://www.metabolomicsworkbench.org/databases/refmet/index.php) standardization is used by default and it can be dectivated with ```use_refmet=False``` option. 
- Lipid names can also be validated separately by importing RefMet class and providing an array of names to validate_metabolite_names() function
```python

from lipidmaps.data import RefMet
lipid_names = ['PC(12:1)','Chol']
refmet_results = RefMet.validate_metabolite_names(lipid_names)

# Output will be RefMetResult object that contains standardized_name and lm_id if recognized 
print(refmet_results)

```

## Adding LIPID MAPS Ids
- RefMet will return lipid maps id's if they are mapped to refmet entries
- Defining molecular species in lipidomics assays are difficult and we rely on generic lipid structures where R groups are not specified. We use headgroup matching for generic structure ID's if RefMet doesn't return lm_id's.  
```python
# Update LIPID MAPS ids by headgroups
# fill_generic_lm_ids_from_headgroups(dataset) will assign headgroup LIPID MAPS ids to lipids as generic lm_id and return the updated count.
updated_count = dataset.fill_generic_lm_ids_from_headgroups()

# List lipid names where an lm id is assigned
print(f"Lipid names with assigned lm ids: {dataset.list_lipids_with_lmid()[:5]}")
```
## Quantitation

- Quantitation values are stored on `QuantifiedLipid` objects (accessible via `dataset.lipids`) as a `values` mapping from sample name to numeric value. The package provides several small helpers to read, query and aggregate those values.

Examples

```python
# Iterate lipids and samples (sample may be a SampleMetadata object or a sample id string)
for lipid in dataset.lipids[:5]:
      for sample in dataset.samples[:5]:
            val = lipid.get_value_for_sample(sample)
            print(f"{sample.sample_name}\t{lipid.input_name}\t{val}")

# Equivalent access via dataset helper (sample can be id or SampleMetadata)
print(dataset.get_value(dataset.samples[0], dataset.lipids[0]))

# Find lipids by name (returns a list of lipid objects where the query matches input or standardized name)
matches = dataset.find_lipids("cardiolipin")

# Get quantified values for a lipid name (returns dict: {"sample_name": value, ...})
values = dataset.get_values("PC(16:0/18:1)")
print(values)

# Per-sample list of lipid values (useful for plotting or tables)
vals = dataset.get_lipid_values_for_samples('sample_01')
# returns something like: [{"input_name": "orig name", "value": 123.4}, ...]

# Aggregations — compute mean across a set of lipid objects for a given sample
class_lipids = dataset.get_lipids_by_generic_lm_id('LMGP01010000')  # example generic LM ID
mean = dataset.mean_value_for_lipids('sample_01', class_lipids, skip_missing=True)

# Fetch reactions by LM ID (attaches ReactionData objects to the dataset)
reactions = dataset.fetch_reactions_by_lm_id(reaction_type="class-level", only_lipid_components=False)
```

### Notes
- `QuantifiedLipid.values` is a plain mapping — you can access values directly (`lipid.values.get('sample_01')`) but prefer the helper methods above to handle missing samples or alternative sample identifiers.
- The Streamlit demo (`scripts/streamlit_demo.py`) uses these helpers to build per-sample charts, class-level summaries and the interactive search UI.


## Reactions

- LIPID MAPS reactions can be fetched by LM ID using `dataset.fetch_reactions_by_lm_id(...)`. The method returns a list of
`ReactionData` objects and also attaches the fetched reactions to the `LipidDataset` instance so downstream code (or the Streamlit demo) can inspect or render them.

Example

```python
# Fetch reactions (optionally restrict by taxonomy_group)
reactions = dataset.fetch_reactions_by_lm_id(
   reaction_type='class-level',
   only_lipid_components=False,
   taxonomy_group='mammalia',
)

# Inspect a few reactions
for rx in reactions[:5]:
   print(rx.reaction_id)
   print('Reactants:', [c.lm_id for c in getattr(rx, 'reactants', [])])
   print('Products: ', [c.lm_id for c in getattr(rx, 'products', [])])

# Get lipid objects participating in a reaction (helper accepts ReactionData or reaction id)
lipids_in_rx = dataset.get_lipids_for_reaction(reactions[0])
```
## Normalization

- Normalized quantitation outputs are stored separately from raw measurements so original
`QuantifiedLipid.values` are never overwritten. Normalized results are kept as list of
`NormalizedResult` objects in `QuantifiedLipid.normalized`.

Example — total-lipid normalization and storing results per-lipid:

```python
from lipidmaps.data.quantitation import QuantitationAnalyzer

# build an analyzer for the dataset
qa = QuantitationAnalyzer(dataset=dataset)

# perform total-lipid normalization (scale factor e.g. 1e6 for ppm)
method_key = "total_lipid:scale=1e6"
norm_map = qa.normalize_total_lipid(scale_factor=1e6)  # returns {input_name: {sample: value, ...}, ...}

# store normalized values on each QuantifiedLipid (non-destructive)
for l in dataset.lipids:
   if l.input_name in norm_map:
      l.set_normalized(method_key, norm_map[l.input_name])

# produce a DataFrame of normalized values (rows=lipids, cols=samples)
df_norm = dataset.normalized_dataframe(method_key)
print(df_norm.head())

# save CSV if desired
df_norm.to_csv("normalized.csv")
```

Accessing stored normalized results for a single lipid:

```python
# returns a NormalizedResult or None
nr = dataset.lipids[0].get_normalized(method_key)
if nr:
   print(nr.values)         # per-sample mapping
   print(nr.meta)           # any provenance metadata
   print(nr.created_iso)    # ISO timestamp when created
```

Note: The Streamlit demo (`scripts/streamlit_demo.py`) includes UI to apply
normalization methods and download the resulting CSV without modifying raw values.

Notes
- Reaction responses may include non-lipid entities; reactants/products can be generic or species-level LM IDs.
- `ReactionData` typically exposes fields like `reaction_id`, `proteins`, `genes`, `curations`, `reactants` and `products`.
- Fetching attaches reactions to the dataset and will annotate matching lipids in-place.

## Queries

The package provides a composable query API for filtering `QuantifiedLipid` objects via `dataset.query_lipids(*preds, combine='and')`.
Predicates can be `Query` objects from `lipidmaps.data.models.query` (e.g. `attr_eq`, `attr_contains`) or plain callables that
accept a `QuantifiedLipid` and return a boolean.

Examples

```python
from lipidmaps.data.models.query import attr_eq, attr_contains

# Find cardiolipins by main class or name
q = attr_eq('main_class', 'Cardiolipins') | attr_contains('input_name', 'cardiolipin')
results = dataset.query_lipids(q, combine='or')

# Combine Query with a callable to filter by a sample value
q2 = attr_eq('main_class', 'Cardiolipins') & (lambda l: l.values.get('Sample1', 0) > 100)
results2 = dataset.query_lipids(q2)
```

Tips
- Use `combine='or'` to broaden matches across predicates.
- Callables are useful for numeric comparisons (sample values, mass ranges) that are awkward to express as attribute queries.
- The README and `scripts/streamlit_demo.py` contain examples showing how to use queries to power interactive selection and plotting.


## Streamlit demo

We have a basic streamlit demo script that you can try.

```python
pip install streamlit
streamlit run scripts/streamlit_demo.py
```

In this demo, you can either use existing csv and tsv files in demo folder or upload your own.
You can preview your data, see validation report and associated LIPID MAPS reactions.

## New / Updated API Methods

The codebase recently added a number of convenience helpers on the dataset and model objects for working with quantitation values and for grouping lipids. These are useful when writing downstream analyses or powering the Streamlit demo.

- `LipidDataset.mean_value_for_lipids(sample, lipids, skip_missing=True)`: compute the mean quantitation for a given sample across a list of lipid objects. `sample` may be a `SampleMetadata` object or a sample id string. Example:

```python
# sample can be a SampleMetadata or sample id
mean = dataset.mean_value_for_lipids(sample, class_lipids, skip_missing=True)
```

- `LipidDataset.get_lipid_values_for_samples(sample_name)`: returns a list of objects describing each lipid's reported value for the given sample, typically in the form `[{"input_name": "orig name", "value": 123.4}, ...]`. Example:

```python
vals = dataset.get_lipid_values_for_samples('sample_01')
```

- `LipidDataset.get_lipids_by_generic_lm_id(generic_lm_id)`: return lipid objects that share a generic LM ID (useful for aggregations by headgroup).

- `LipidDataset.get_lipids_for_reaction(reaction_or_id)`: return lipid objects participating in a reaction (accepts a `ReactionData` object or a reaction id).

- `QuantifiedLipid.get_value_for_sample(sample_or_id)`: convenience on the lipid object to fetch a single sample's value from the lipid's internal `values` mapping.

These helpers are used by `scripts/streamlit_demo.py` to build per-sample and per-lipid views, compute mean values by class, and annotate reactions. 

- `LipidDataset.query_lipids(*preds, combine='and')`: composable query API for filtering lipids. Predicates may be
   `Query` objects from `lipidmaps.data.models.query` or plain callables that accept a `QuantifiedLipid`.

Example:
```python
from lipidmaps.data.models.query import attr_eq, attr_contains

# find cardiolipins or any lipid with "cardiolipin" in the input name
q = attr_eq('main_class', 'Cardiolipins') | attr_contains('input_name', 'cardiolipin')
results = dataset.query_lipids(q, combine='or')
print(len(results))

# combine a Query with a callable to filter by a sample value
q2 = attr_eq('main_class', 'Cardiolipins') & (lambda l: l.values.get('Sample1', 0) > 100)
results2 = dataset.query_lipids(q2)
```
## Example Datasets

Sample datasets are available in `tests/data/inputs/`:
- `small_demo.csv`: Small example dataset for quick testing
- `large_demo.csv`: Larger dataset for comprehensive testing

## Documentation

For more detailed documentation, see:
- `docs/custom_columns_guide.md`: Guide for working with custom data columns
- `INSTALL.md`: Detailed installation instructions

## Project Structure

```
lipidmaps_py/
├── src/lipidmaps/           # Main package code
│   ├── data/                # Data analysis module
│   │   ├── models/         # Data models
│   │   ├── ingestion/      # Data import
│   │   ├── validation/     # Data validation
│   │   └── config/         # Configuration
│   └── tools/              # Utility tools
├── tests/                   # Test suite
│   └── data/               # Data module tests
│       └── inputs/         # Sample datasets
├── docs/                    # Documentation
└── scripts/                 # Demo scripts
```

## Troubleshooting

### SQLite3 Module Not Found

If you get `ModuleNotFoundError: No module named '_sqlite3'`:

1. **Use system Python** instead of custom-built Python:
   ```bash
   /usr/bin/python3 -m venv venv # or /bin/python3 
   source venv/bin/activate
   pip install -e .
   ```

2. **Or rebuild Python with SQLite support**:
   ```bash
   sudo dnf install sqlite-devel  # CentOS/RHEL
   # Then rebuild and reinstall Python from source
   ```

### Import Errors

If you get import errors, make sure the package is installed:
```bash
pip install -e .
```

### Test Failures

If tests fail, ensure you have all dependencies:
```bash
pip install pytest pandas numpy requests
```

## Contributing

Contributions are welcome! Please ensure:
1. All tests pass: `pytest`
2. Code follows the project style
3. New features include tests
4. Documentation is updated

## License

See LICENSE file for details.

## Support

For issues, questions, or contributions, please visit the [GitHub repository](https://github.com/s-andrews/lipidmaps_py).
