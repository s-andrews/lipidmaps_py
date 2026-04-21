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
        cls.test_file = Path(__file__).parent / "data" / "inputs" / "small_demo.csv"
        if not cls.test_file.exists():
            raise FileNotFoundError(f"Test file not found: {cls.test_file}")

    def test_import_data_basic(self):
        """Test basic import_data functionality."""
        # Import CSV data
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

        # Verify data was imported
        self.assertIsInstance(data, LipidData)
        self.assertGreater(data.successful_import_count(), 0)
        self.assertIsInstance(data.dataset.samples, list)
        self.assertGreater(len(data.dataset.samples), 0)

    def test_lipid_data_counts(self):
        """Test LipidData count methods."""
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

        # Test count methods
        successful = data.successful_import_count()
        failed = data.failed_import_count()
        total_lipids = len(data.dataset.lipids)

        self.assertEqual(successful, total_lipids)
        self.assertGreaterEqual(failed, 0)

    def test_lipid_data_samples(self):
        """Test sample access methods."""
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

        # Get samples
        samples = data.dataset.samples
        self.assertIsInstance(samples, list)

    def test_as_dataframe(self):
        """Test pandas DataFrame export."""
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

        # Export to DataFrame
        df = data.as_dataframe()

        # Verify DataFrame structure
        import pandas as pd

        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), len(data.dataset.lipids))

    def test_get_group_statistics(self):
        """Test group-level statistical analysis."""
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

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
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

        # Serialize to dict
        data_dict = data.to_dict()

        self.assertIsInstance(data_dict, dict)
        self.assertIn("samples", data_dict)
        self.assertIn("lipids", data_dict)


    def test_usage_example_from_docs(self):
        """Test the usage example from documentation."""
        # Import CSV data
        data = lipidmaps.import_data(str(self.test_file), use_refmet=False, use_headgroups=False, fetch_reactions=False)

        # Access data
        imported_count = data.successful_import_count()
        self.assertGreater(imported_count, 0)

        samples = data.dataset.samples
        self.assertIsInstance(samples, list)

        # Export to pandas
        df = data.as_dataframe()
        import pandas as pd

        self.assertIsInstance(df, pd.DataFrame)

        # Group statistics
        stats = data.get_group_statistics()
        self.assertIsInstance(stats, dict)


if __name__ == "__main__":
    unittest.main()
