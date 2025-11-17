"""
Tests for lipidmaps.import_data() and LipidData high-level API
"""
import unittest
import logging
from pathlib import Path

import lipidmaps
from lipidmaps import LipidData


logging.basicConfig(
    level=logging.WARNING, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class TestDataImporter(unittest.TestCase):
    """Test the high-level import_data API and LipidData class."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.test_file = Path(__file__).parent / "biopan" / "inputs" / "biopan_small_demo.csv"
        if not cls.test_file.exists():
            raise FileNotFoundError(f"Test file not found: {cls.test_file}")

    def test_import_data_basic(self):
        """Test basic import_data functionality."""
        # Import CSV data
        data = lipidmaps.import_data(str(self.test_file))

        # Verify data was imported
        self.assertIsInstance(data, LipidData)
        self.assertGreater(data.successful_import_count(), 0)
        self.assertIsInstance(data.samples(), list)
        self.assertGreater(len(data.samples()), 0)

    def test_lipid_data_counts(self):
        """Test LipidData count methods."""
        data = lipidmaps.import_data(str(self.test_file))

        # Test count methods
        successful = data.successful_import_count()
        failed = data.failed_import_count()
        total_lipids = len(data.lipids())

        self.assertEqual(successful, total_lipids)
        self.assertGreaterEqual(failed, 0)
        self.assertIsInstance(data.failed_import_names(), list)

    def test_lipid_data_samples(self):
        """Test sample access methods."""
        data = lipidmaps.import_data(str(self.test_file))

        # Get samples
        samples = data.samples()
        self.assertIsInstance(samples, list)
        self.assertTrue(all(isinstance(s, str) for s in samples))

        # Verify sample_names property
        self.assertEqual(samples, data.sample_names)

    def test_get_lm_ids(self):
        """Test extraction of LIPID MAPS IDs."""
        data = lipidmaps.import_data(str(self.test_file))

        # Get LM IDs
        lm_ids = data.get_lm_ids()
        self.assertIsInstance(lm_ids, list)

        # Verify all IDs start with 'LM'
        for lm_id in lm_ids:
            self.assertTrue(lm_id.startswith("LM"), f"Invalid LM ID: {lm_id}")

        # Verify uniqueness
        self.assertEqual(len(lm_ids), len(set(lm_ids)))

    def test_get_lipid_by_name(self):
        """Test retrieving specific lipids by name."""
        data = lipidmaps.import_data(str(self.test_file))

        if data.lipids():
            # Get first lipid
            first_lipid = data.lipids()[0]
            input_name = first_lipid.input_name

            # Retrieve by input name
            retrieved = data.get_lipid_by_name(input_name)
            self.assertIsNotNone(retrieved)
            self.assertEqual(retrieved.input_name, input_name)

            # Try retrieving by standardized name if available
            if first_lipid.standardized_name:
                retrieved_by_std = data.get_lipid_by_name(first_lipid.standardized_name)
                if retrieved_by_std:
                    self.assertEqual(retrieved_by_std.input_name, input_name)

    def test_get_value_for_lipid(self):
        """Test getting quantitation values for specific lipid/sample pairs."""
        data = lipidmaps.import_data(str(self.test_file))

        if data.lipids() and data.samples():
            first_lipid = data.lipids()[0]
            lipid_name = first_lipid.input_name

            # Get value for each sample
            for sample in data.samples():
                value = data.get_value_for_lipid(lipid_name, sample)
                if value is not None:
                    self.assertIsInstance(value, float)

            # Test with QuantifiedLipid object directly
            if data.samples():
                sample = data.samples()[0]
                value = data.get_value_for_lipid(first_lipid, sample)
                if value is not None:
                    self.assertIsInstance(value, float)

    def test_get_lipids_by_class(self):
        """Test filtering lipids by class."""
        data = lipidmaps.import_data(str(self.test_file))

        # Find a lipid with a class
        lipid_with_class = None
        for lipid in data.lipids():
            if lipid.sub_class or lipid.main_class:
                lipid_with_class = lipid
                break

        if lipid_with_class:
            # Get lipids by sub_class
            if lipid_with_class.sub_class:
                lipids_in_class = data.get_lipids_by_class(lipid_with_class.sub_class)
                self.assertIsInstance(lipids_in_class, list)
                self.assertGreater(len(lipids_in_class), 0)
                # Verify the original lipid is in the results
                self.assertIn(lipid_with_class, lipids_in_class)

    def test_as_dataframe(self):
        """Test pandas DataFrame export."""
        data = lipidmaps.import_data(str(self.test_file))

        # Export to DataFrame
        df = data.as_dataframe()

        # Verify DataFrame structure
        import pandas as pd
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), len(data.lipids()))
        # Columns should be samples
        self.assertTrue(any(s in df.columns for s in data.samples()))

    def test_get_group_statistics(self):
        """Test group-level statistical analysis."""
        data = lipidmaps.import_data(str(self.test_file))

        # Get group statistics
        stats = data.get_group_statistics()

        self.assertIsInstance(stats, dict)
        # Should have at least one group
        if stats:
            for group_name, group_stats in stats.items():
                self.assertIsInstance(group_name, str)
                self.assertIn("sample_count", group_stats)
                self.assertIn("lipid_coverage", group_stats)
                self.assertIn("mean_values", group_stats)
                self.assertIn("std_values", group_stats)
                self.assertIsInstance(group_stats["mean_values"], dict)
                self.assertIsInstance(group_stats["std_values"], dict)

    def test_to_dict(self):
        """Test dictionary serialization."""
        data = lipidmaps.import_data(str(self.test_file))

        # Serialize to dict
        data_dict = data.to_dict()

        self.assertIsInstance(data_dict, dict)
        self.assertIn("samples", data_dict)
        self.assertIn("lipids", data_dict)

    def test_import_msdial(self):
        """Test import_msdial function."""
        # For now, import_msdial should work the same as import_data
        data = lipidmaps.import_msdial(str(self.test_file))

        self.assertIsInstance(data, LipidData)
        self.assertGreater(data.successful_import_count(), 0)

    def test_usage_example_from_docs(self):
        """Test the usage example from documentation."""
        # Import CSV data
        data = lipidmaps.import_data(str(self.test_file))

        # Access data
        imported_count = data.successful_import_count()
        self.assertGreater(imported_count, 0)

        samples = data.samples()
        self.assertIsInstance(samples, list)

        lm_ids = data.get_lm_ids()
        self.assertIsInstance(lm_ids, list)

        # Get specific lipid (if we have lipids with standardized names)
        for lipid in data.lipids()[:5]:  # Check first 5 lipids
            retrieved = data.get_lipid_by_name(lipid.input_name)
            if retrieved:
                self.assertEqual(retrieved.input_name, lipid.input_name)

                # Get value for first sample if available
                if samples:
                    value = data.get_value_for_lipid(lipid.input_name, samples[0])
                    # Value might be None if this sample doesn't have data for this lipid
                    if value is not None:
                        self.assertIsInstance(value, float)

        # Get all lipids of a specific class (if any exist)
        for lipid in data.lipids():
            if lipid.sub_class:
                lipids_in_class = data.get_lipids_by_class(lipid.sub_class)
                self.assertIsInstance(lipids_in_class, list)
                break

        # Export to pandas
        df = data.as_dataframe()
        import pandas as pd
        self.assertIsInstance(df, pd.DataFrame)

        # Group statistics
        stats = data.get_group_statistics()
        self.assertIsInstance(stats, dict)

    def test_reactions_not_implemented(self):
        """Test that reaction methods raise NotImplementedError."""
        data = lipidmaps.import_data(str(self.test_file))

        with self.assertRaises(NotImplementedError):
            data.get_reactions()

        with self.assertRaises(NotImplementedError):
            data.get_lipids_for_reaction_component("component")

        with self.assertRaises(NotImplementedError):
            data.get_value_for_reaction_component("component")


if __name__ == "__main__":
    unittest.main()
