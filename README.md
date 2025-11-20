# lipidmaps_py

A Python package to interface with the LIPID MAPS database.

## Development Status

### ✅ Complete
- **Data Import & Validation**: CSV/TSV data ingestion with format detection
- **Data Normalization**: RefMet standardization
- **Quality Control**: Data validation and issue reporting
- **Data Management**: DataManager for handling quantified lipid datasets
- **Sample Metadata**: Support for experimental metadata and conditions

### 🚧 In Progress
- **LIPID MAPS API Integration**: LM ID lookup and validation
- **Reaction Analysis**: Integration with LIPID MAPS reactions database

## Installation

### Prerequisites
- Python 3.8 or higher (with SQLite3 support if test coverage report is used)
- pip (Python package installer)

> **Note**: If you encounter `ModuleNotFoundError: No module named '_sqlite3'`, your Python installation was built without SQLite support. Either:
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
# Use system Python if your custom Python lacks SQLite support
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the package in development mode:
```bash
pip install -e .
```

4. Install development dependencies (optional, for running tests):
```bash
pip install pytest pytest-cov black flake8 mypy
```

### Verify Installation

Test that the package is correctly installed:
```bash
python -c "import lipidmaps; print('Installation successful')"
```

## Running Tests

The package includes a comprehensive test suite. To run all tests:

```bash
# Run all tests
pytest

# Run tests with verbose output
pytest -v

# Run tests for specific module
pytest tests/data/

# Run specific test file
pytest tests/data/test_data_manager.py

# Run with coverage report (requires SQLite support)
pytest --cov=lipidmaps --cov-report=html

# If coverage fails due to SQLite issues, run tests without coverage
pytest -v
```

### Test Coverage

View test coverage by opening `htmlcov/index.html` in your browser after running with coverage.

## Quick Start

### Basic Usage

```python
from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.ingestion import CSVIngestion

# Load a CSV file
ingestion = CSVIngestion("path/to/your/data.csv")
raw_data = ingestion.load()

# Create a DataManager instance
dm = DataManager(raw_data)

# Access quantified lipids
for lipid in dm.get_quantified_lipids():
    print(f"{lipid.name}: {lipid.abundance}")
```

### Working with Sample Data

The package includes sample datasets in the `tests/data/inputs/` directory:

```python
from lipidmaps.data.data_manager import DataManager
from lipidmaps.data.ingestion import CSVIngestion

# Load the demo dataset
ingestion = CSVIngestion("tests/data/inputs/small_demo.csv")
raw_data = ingestion.load()

# Create DataManager
dm = DataManager(raw_data)

# Get dataset information
print(f"Number of lipids: {len(dm.get_quantified_lipids())}")
print(f"Number of samples: {len(dm.get_samples())}")
```

### Data Validation

```python
from lipidmaps.data.ingestion import CSVIngestion
from lipidmaps.data.validation import DataValidator

# Load and validate data
ingestion = CSVIngestion("path/to/data.csv")
raw_data = ingestion.load()

validator = DataValidator(raw_data)
report = validator.validate()

# Print validation issues
for issue in report.get_all_issues():
    print(f"{issue.severity}: {issue.message}")
```

### RefMet Standardization

```python
from lipidmaps.data.models.refmet import RefMet

# Standardize lipid names
refmet = RefMet()
result = refmet.standardize("PC(16:0/18:1)")

if result.success:
    print(f"Standardized name: {result.standardized_name}")
    print(f"RefMet ID: {result.refmet_id}")
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
│   ├── data/                # Data analysis module (formerly biopan)
│   │   ├── models/         # Data models
│   │   ├── ingestion/      # Data import
│   │   ├── validation/     # Data validation
│   │   └── config/         # Configuration
│   └── tools/              # Utility tools
├── tests/                   # Test suite
│   └── data/               # Data module tests
│       └── inputs/         # Sample datasets
└── docs/                    # Documentation
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
