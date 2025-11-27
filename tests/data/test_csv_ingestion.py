"""
Tests for CSVIngestion class
"""

import unittest
import tempfile
import logging
from pathlib import Path

from lipidmaps.data.ingestion import CSVIngestion, RawDataFrame, CSVFormat
from lipidmaps.data.validation.data_validator import DataValidator, IssueSeverity


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class TestCSVIngestion(unittest.TestCase):
    """Test CSVIngestion class for CSV file reading."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path(__file__).parent / "inputs"
        self.ingestion = CSVIngestion()

    def test_read_standard_csv(self):
        """Test reading a standard CSV file."""
        test_file = self.test_data_dir / "small_demo.csv"

        raw_df = self.ingestion.read_csv(test_file, CSVFormat.STANDARD)

        self.assertIsInstance(raw_df, RawDataFrame)
        self.assertGreater(raw_df.row_count, 0)
        self.assertGreater(raw_df.column_count, 0)
        self.assertEqual(raw_df.format_type, CSVFormat.STANDARD)
        self.assertIn("source_file", raw_df.metadata)

    def test_auto_detect_format(self):
        """Test automatic format detection."""
        test_file = self.test_data_dir / "small_demo.csv"

        raw_df = self.ingestion.read_csv(test_file, CSVFormat.AUTO)

        self.assertIsInstance(raw_df, RawDataFrame)
        self.assertIn(raw_df.format_type, [CSVFormat.STANDARD, CSVFormat.MSDIAL])

    def test_read_nonexistent_file(self):
        """Test reading a file that doesn't exist."""
        nonexistent = Path("/nonexistent/file.csv")

        with self.assertRaises(FileNotFoundError):
            self.ingestion.read_csv(nonexistent)

    def test_empty_csv(self):
        """Test reading an empty CSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("lipid,sample1\n")  # Header only
            temp_path = Path(f.name)

        try:
            raw_df = self.ingestion.read_csv(temp_path)
            self.assertTrue(raw_df.is_empty())
            self.assertEqual(raw_df.row_count, 0)
        finally:
            temp_path.unlink()

    def test_delimiter_detection(self):
        """Test automatic delimiter detection."""
        # Create CSV with tab delimiter
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("lipid\tsample1\tsample2\n")
            f.write("PC(16:0/18:1)\t123.45\t678.90\n")
            temp_path = Path(f.name)

        try:
            ingestion = CSVIngestion()  # No delimiter specified
            raw_df = ingestion.read_csv(temp_path)

            self.assertEqual(raw_df.column_count, 3)
            self.assertEqual(raw_df.row_count, 1)
        finally:
            temp_path.unlink()

    def test_raw_dataframe_properties(self):
        """Test RawDataFrame properties."""
        test_file = self.test_data_dir / "small_demo.csv"
        raw_df = self.ingestion.read_csv(test_file)

        # Test properties
        self.assertEqual(raw_df.row_count, len(raw_df.rows))
        self.assertEqual(raw_df.column_count, len(raw_df.fieldnames))
        self.assertFalse(raw_df.is_empty())
        self.assertIsInstance(raw_df.metadata, dict)

    def test_get_column_info(self):
        """Test getting column information."""
        test_file = self.test_data_dir / "small_demo.csv"
        raw_df = self.ingestion.read_csv(test_file)

        col_info = self.ingestion.get_column_info(raw_df)

        self.assertIn("column_count", col_info)
        self.assertIn("columns", col_info)
        self.assertIn("column_types", col_info)
        self.assertEqual(col_info["column_count"], len(raw_df.fieldnames))

    def test_read_batch(self):
        """Test reading multiple files at once."""
        test_files = [
            self.test_data_dir / "small_demo.csv",
        ]

        results = self.ingestion.read_batch(test_files)

        self.assertIsInstance(results, list)
        self.assertEqual(len(results), len(test_files))
        for raw_df in results:
            self.assertIsInstance(raw_df, RawDataFrame)

    def test_custom_delimiter(self):
        """Test reading CSV with custom delimiter."""
        # Create CSV with semicolon delimiter
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("lipid;sample1;sample2\n")
            f.write("PC(16:0/18:1);123.45;678.90\n")
            temp_path = Path(f.name)

        try:
            ingestion = CSVIngestion(delimiter=";")
            raw_df = ingestion.read_csv(temp_path)

            self.assertEqual(raw_df.column_count, 3)
            self.assertEqual(raw_df.row_count, 1)
        finally:
            temp_path.unlink()

    def test_csv_with_missing_values(self):
        """Test reading CSV with missing values."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("lipid,sample1,sample2\n")
            f.write("PC(16:0/18:1),123.45,\n")
            f.write("TAG(16:0/18:1/18:2),,678.90\n")
            temp_path = Path(f.name)

        try:
            raw_df = self.ingestion.read_csv(temp_path)

            self.assertEqual(raw_df.row_count, 2)
            # Check that missing values are empty strings
            self.assertEqual(raw_df.rows[0]["sample2"], "")
            self.assertEqual(raw_df.rows[1]["sample1"], "")
        finally:
            temp_path.unlink()


class TestCSVStructureValidation(unittest.TestCase):
    """Validate CSV structural rules."""

    def setUp(self):
        self.test_data_dir = Path(__file__).parent / "inputs"
        self.ingestion = CSVIngestion()
        self.validator = DataValidator(min_samples=2, min_lipids=1)

    def _validate_file(self, filename: str):
        path = self.test_data_dir / filename
        raw_df = self.ingestion.read_csv(path)
        return self.validator.validate(raw_df)

    def test_input_positive_structure_passes(self):
        report = self._validate_file("input_positive.csv")
        self.assertTrue(
            report.passed, "Expected positive input file to pass validation"
        )
        self.assertEqual(len(report.get_issues_by_severity(IssueSeverity.ERROR)), 0)
        self.assertEqual(len(report.get_issues_by_severity(IssueSeverity.CRITICAL)), 0)

    def test_input_empty_file_reports_error(self):
        report = self._validate_file("input_empty.csv")
        self.assertFalse(report.passed)
        critical_issues = report.get_issues_by_severity(IssueSeverity.CRITICAL)
        self.assertTrue(
            any("Dataset is empty" in issue.message for issue in critical_issues)
        )

    def test_missing_lipid_names_are_flagged(self):
        report = self._validate_file("input_missing_names.csv")
        name_issues = [
            issue for issue in report.issues if issue.category == "missing_data"
        ]
        self.assertGreater(
            len(name_issues), 0, "Expected missing lipid names to be reported"
        )
        self.assertFalse(report.passed)

    def test_inconsistent_field_counts_detected(self):
        report = self._validate_file("input_inconsistent_fields.csv")
        structure_errors = [
            issue for issue in report.issues if issue.category == "structure"
        ]
        self.assertGreater(len(structure_errors), 0)
        self.assertTrue(
            any(
                "missing values" in issue.message or "extra fields" in issue.message
                for issue in structure_errors
            )
        )

    def test_non_numeric_quantitation_values_flagged(self):
        report = self._validate_file("input_non_numeric.csv")
        invalid_value_errors = [
            issue for issue in report.issues if issue.category == "invalid_value"
        ]
        self.assertGreater(len(invalid_value_errors), 0)
        self.assertFalse(report.passed)


if __name__ == "__main__":
    unittest.main()
