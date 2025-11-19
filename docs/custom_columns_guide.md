# Custom Column Specification and Group Mapping

## Overview

The `import_data()` function and `DataManager` now support flexible CSV column specification and custom group-to-sample mappings. This allows users to:

1. **Specify which column contains lipid names** (by index or name)
2. **Specify which columns contain sample data** (by index or name)
3. **Define custom groups** and associate them with specific samples

## Basic Usage

### Default Behavior

By default, the first column is treated as lipid names, and all other columns are treated as samples:

```python
from lipidmaps.data_importer import import_data

# Auto-detect: first column = lipids, rest = samples
data = import_data("lipids.csv")
```

### Specify Columns by Index (0-based)

```python
# Lipid names in column 2, sample data in columns 3-6
data = import_data(
    "lipids.csv",
    lipid_col=2,
    sample_cols=[3, 4, 5, 6]
)
```

### Specify Columns by Name

```python
# Use column names from CSV header
data = import_data(
    "lipids.csv",
    lipid_col="LipidName",
    sample_cols=["Control1", "Control2", "Treatment1", "Treatment2"]
)
```

### Mix Index and Name Specification

```python
# Can mix indices and names
data = import_data(
    "lipids.csv",
    lipid_col=2,  # Column index
    sample_cols=["Control1", "Control2"]  # Column names
)
```

## Group Mapping

### Custom Group Assignment

Instead of auto-extracting groups from sample names, you can explicitly define groups:

```python
data = import_data(
    "lipids.csv",
    group_mapping={
        "Control": ["Sample1", "Sample2", "Sample3"],
        "Treatment": ["Sample4", "Sample5", "Sample6"],
        "HighDose": ["Sample7", "Sample8"]
    }
)
```

### Accessing Group Information

```python
# View sample-to-group assignments
for sample in data.dataset.samples:
    print(f"{sample.sample_id} -> {sample.group}")

# Get statistics by group
stats = data.get_group_statistics()
for group_name, group_stats in stats.items():
    print(f"{group_name}: {group_stats['sample_count']} samples")
```

### Combined: Custom Columns + Group Mapping

```python
data = import_data(
    "lipids.csv",
    lipid_col="LipidName",
    sample_cols=["S1", "S2", "S3", "S4"],
    group_mapping={
        "Control": ["S1", "S2"],
        "Treatment": ["S3", "S4"]
    }
)
```

## Data Validation

Enable validation to check data quality:

```python
data = import_data("lipids.csv", validate=True)

# Check validation results
if data.manager.validation_report:
    if not data.manager.validation_report.passed:
        data.manager.validation_report.print_report()
```

## Using DataManager Directly

For more control, use `DataManager` directly:

```python
from lipidmaps.biopan.data_manager import DataManager

manager = DataManager(
    lipid_name_column="LipidName",
    sample_columns=["C1", "C2", "T1", "T2"],
    group_mapping={
        "Control": ["C1", "C2"],
        "Treatment": ["T1", "T2"]
    },
    validate_data=True
)

dataset = manager.process_csv("lipids.csv")
```

## Example CSV Formats

### Standard Format (Default)
```csv
Lipid,Sample1,Sample2,Sample3
PC(16:0/18:1),100,110,105
TAG(54:3),200,210,205
```

### Custom Format with Extra Columns
```csv
ID,Metadata,LipidName,Control1,Control2,Treatment1,Extra
1,foo,PC(16:0/18:1),100,110,150,ignore
2,bar,TAG(54:3),200,210,180,ignore
```

Use with:
```python
data = import_data(
    "custom.csv",
    lipid_col="LipidName",
    sample_cols=["Control1", "Control2", "Treatment1"]
)
```

## Error Handling

The system validates column specifications and provides clear error messages:

```python
# Invalid column index
data = import_data("lipids.csv", lipid_col=99)
# ValueError: lipid_name_column index 99 out of range. CSV has 5 columns.

# Invalid column name
data = import_data("lipids.csv", lipid_col="NonExistent")
# ValueError: lipid_name_column 'NonExistent' not found in CSV. Available columns: [...]
```

## Tips

1. **Use column names** when possible - more readable and resilient to column reordering
2. **Validate your data** with `validate=True` to catch quality issues early
3. **Define explicit groups** when sample names don't follow a pattern
4. **Check group assignments** after import to ensure they match your expectations
5. **Use group statistics** to verify data distribution across experimental groups
