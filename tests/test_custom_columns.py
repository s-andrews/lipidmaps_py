"""
Tests for custom column specification and group mapping in DataManager and import_data.
"""

import unittest
import tempfile
from pathlib import Path

from lipidmaps.data_importer import import_data
from lipidmaps.data.data_manager import DataManager


class TestCustomColumns(unittest.TestCase):
    """Test custom column specification features."""

    def setUp(self):
        """Create temporary test CSV files."""
        # Standard CSV with default layout
        self.standard_csv = tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".csv"
        )
        self.standard_csv.write("Lipid,Sample1,Sample2,Sample3\n")
        self.standard_csv.write("PC(16:0/18:1),100,110,105\n")
        self.standard_csv.write("TAG(54:3),200,210,205\n")
        self.standard_csv.close()

        # CSV with custom column order
        self.custom_csv = tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".csv"
        )
        self.custom_csv.write("ID,Extra,LipidName,Value1,Value2,Ignore\n")
        self.custom_csv.write("1,foo,PC(16:0/18:1),100,110,bar\n")
        self.custom_csv.write("2,baz,TAG(54:3),200,210,qux\n")
        self.custom_csv.close()

    def tearDown(self):
        """Clean up temporary files."""
        Path(self.standard_csv.name).unlink(missing_ok=True)
        Path(self.custom_csv.name).unlink(missing_ok=True)

    def test_default_column_behavior(self):
        """Test default behavior: first column is lipid, rest are samples."""
        data = import_data(self.standard_csv.name)

        self.assertEqual(len(data.lipids()), 2)
        self.assertEqual(len(data.samples()), 3)
        self.assertIn("Sample1", data.sample_names)
        self.assertIn("Sample2", data.sample_names)
        self.assertIn("Sample3", data.sample_names)

    def test_custom_lipid_column_by_index(self):
        """Test specifying lipid column by index."""
        data = import_data(
            self.custom_csv.name,
            lipid_col=2,  # LipidName column
            sample_cols=[3, 4],  # Value1, Value2
        )

        self.assertEqual(len(data.lipids()), 2)
        self.assertEqual(len(data.samples()), 2)
        self.assertIn("Value1", data.sample_names)
        self.assertIn("Value2", data.sample_names)

    def test_custom_lipid_column_by_name(self):
        """Test specifying lipid column by name."""
        data = import_data(
            self.custom_csv.name,
            lipid_col="LipidName",
            sample_cols=["Value1", "Value2"],
        )

        self.assertEqual(len(data.lipids()), 2)
        self.assertEqual(len(data.samples()), 2)
        self.assertIn("Value1", data.sample_names)
        self.assertIn("Value2", data.sample_names)

    def test_mixed_column_specification(self):
        """Test mixing index and name specifications."""
        data = import_data(
            self.custom_csv.name,
            lipid_col=2,  # Index for LipidName
            sample_cols=["Value1", "Value2"],  # Names for samples
        )

        self.assertEqual(len(data.lipids()), 2)
        self.assertEqual(len(data.samples()), 2)

    def test_group_mapping(self):
        """Test explicit group-to-sample mapping."""
        data = import_data(
            self.standard_csv.name,
            group_mapping={"Control": ["Sample1", "Sample2"], "Treatment": ["Sample3"]},
        )

        self.assertEqual(len(data.lipids()), 2)
        self.assertEqual(len(data.samples()), 3)

        # Check group assignments
        samples_dict = {s.sample_id: s.group for s in data.dataset.samples}
        self.assertEqual(samples_dict["Sample1"], "Control")
        self.assertEqual(samples_dict["Sample2"], "Control")
        self.assertEqual(samples_dict["Sample3"], "Treatment")

    def test_group_mapping_with_custom_columns(self):
        """Test group mapping combined with custom column specification."""
        data = import_data(
            self.custom_csv.name,
            lipid_col="LipidName",
            sample_cols=["Value1", "Value2"],
            group_mapping={"GroupA": ["Value1"], "GroupB": ["Value2"]},
        )

        self.assertEqual(len(data.lipids()), 2)
        self.assertEqual(len(data.samples()), 2)

        samples_dict = {s.sample_id: s.group for s in data.dataset.samples}
        self.assertEqual(samples_dict["Value1"], "GroupA")
        self.assertEqual(samples_dict["Value2"], "GroupB")

    def test_invalid_lipid_column_index(self):
        """Test error handling for invalid lipid column index."""
        with self.assertRaises(ValueError) as ctx:
            import_data(self.standard_csv.name, lipid_col=99)
        self.assertIn("out of range", str(ctx.exception))

    def test_invalid_lipid_column_name(self):
        """Test error handling for invalid lipid column name."""
        with self.assertRaises(ValueError) as ctx:
            import_data(self.standard_csv.name, lipid_col="NonExistent")
        self.assertIn("not found", str(ctx.exception))

    def test_invalid_sample_column_index(self):
        """Test error handling for invalid sample column index."""
        with self.assertRaises(ValueError) as ctx:
            import_data(self.standard_csv.name, sample_cols=[99])
        self.assertIn("out of range", str(ctx.exception))

    def test_invalid_sample_column_name(self):
        """Test error handling for invalid sample column name."""
        with self.assertRaises(ValueError) as ctx:
            import_data(self.standard_csv.name, sample_cols=["NonExistent"])
        self.assertIn("not found", str(ctx.exception))

    def test_data_manager_direct_usage(self):
        """Test DataManager with custom configuration directly."""
        manager = DataManager(
            lipid_name_column="LipidName",
            sample_columns=["Value1", "Value2"],
            group_mapping={"Control": ["Value1"], "Treatment": ["Value2"]},
        )

        dataset = manager.process_csv(self.custom_csv.name)

        self.assertEqual(len(dataset.lipids), 2)
        self.assertEqual(len(dataset.samples), 2)

        # Verify group assignments
        samples_dict = {s.sample_id: s.group for s in dataset.samples}
        self.assertEqual(samples_dict["Value1"], "Control")
        self.assertEqual(samples_dict["Value2"], "Treatment")

    def test_partial_group_mapping(self):
        """Test group mapping with some samples unmapped (should fall back to auto-detection)."""
        data = import_data(
            self.standard_csv.name,
            group_mapping={
                "Control": ["Sample1"],
                # Sample2 and Sample3 not mapped
            },
        )

        samples_dict = {s.sample_id: s.group for s in data.dataset.samples}
        self.assertEqual(samples_dict["Sample1"], "Control")
        # Sample2 and Sample3 should use auto-detection
        self.assertIn(samples_dict["Sample2"], ["Sample", "unknown"])
        self.assertIn(samples_dict["Sample3"], ["Sample", "unknown"])

    def test_group_statistics_with_mapping(self):
        """Test that group statistics work correctly with custom group mapping."""
        data = import_data(
            self.standard_csv.name,
            group_mapping={"Control": ["Sample1", "Sample2"], "Treatment": ["Sample3"]},
        )

        stats = data.get_group_statistics()

        self.assertIn("Control", stats)
        self.assertIn("Treatment", stats)
        self.assertEqual(stats["Control"]["sample_count"], 2)
        self.assertEqual(stats["Treatment"]["sample_count"], 1)


if __name__ == "__main__":
    unittest.main()
