import os
import csv
import unittest

from biopan.src.data_manager import DataManager
from biopan.src.models.sample import SampleMetadata, QuantifiedLipid

class TestPopulateManager(unittest.TestCase):
    def test_populate_manager_from_csv(self):
        """
        Read tests/inputs/val_subset_test.csv (expected format: first column = lipid id/name,
        remaining columns = sample ids with numeric values) and populate a DataManager
        with SampleMetadata and QuantifiedLipid instances.
        """
    csv_path = os.path.join(os.path.dirname(__file__), "inputs", "val_subset_test.csv")
    assert os.path.exists(csv_path), f"CSV not found: {csv_path}"

    with open(csv_path, newline='') as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    assert rows, "Input CSV is empty"
    assert len(fieldnames) >= 2, "CSV must have at least one lipid column and one sample column"

    first_col = fieldnames[0]
    sample_ids = fieldnames[1:]

    # build sample metadata
    samples_meta = [SampleMetadata(sample_id=sid, group="unknown") for sid in sample_ids]

    # build quantified lipids
    lipid_data = []
    for row in rows:
        lipid_species = (row.get(first_col) or "").strip()
        if not lipid_species:
            continue
        values = {}
        for sid in sample_ids:
            raw = (row.get(sid) or "").strip()
            if raw == "":
                continue
            try:
                values[sid] = float(raw)
            except ValueError:
                # skip non-numeric cells
                continue
        lipid_data.append(QuantifiedLipid(input_name=lipid_species, values=values))

    # Use DataManager to process the CSV and get the dataset
    manager = DataManager()
    dataset = manager.process_csv(csv_path)

    # basic assertions
    assert len(dataset.samples) == len(sample_ids)
    assert len(dataset.lipids) == len(lipid_data)
    # spot-check a QuantifiedLipid shape
    assert isinstance(dataset.lipids[0].values, dict)
    # ensure DataManager initialized and dataset populated
    assert hasattr(manager, "process_csv")
