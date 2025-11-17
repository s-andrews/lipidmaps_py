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
from lipidmaps import biopan, tools

# Import data
lipid_data = lipidmaps.import_data("mydata.csv", lipid_col=1, sample_cols=[4,5,6,7])

# Or import MS-DIAL format
lipid_data = lipidmaps.import_msdial("mydata_msdial.csv")

# Access imported data
print(f"Imported Lipids: {lipid_data.successful_import_count()}")
print(f"Unrecognised Lipids: {lipid_data.failed_import_count()}")

# Get reactions
reactions = lipid_data.get_reactions(species="human", complete=True)

# Use biopan subpackage
import lipidmaps.biopan as biopan
# Or: from lipidmaps.biopan import data_manager
```

## Project Structure

```
lipidmaps_py/
├── setup.py              # Package configuration
├── src/
│   └── lipidmaps/       # Main package
│       ├── __init__.py
│       ├── data_importer.py
│       ├── biopan/      # BioPAN framework (subpackage)
│       │   ├── __init__.py
│       │   ├── data_manager.py
│       │   ├── reaction_checker.py
│       │   ├── config/
│       │   ├── models/
│       │   └── utils/
│       └── tools/       # Utility tools (subpackage)
│           └── __init__.py
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
