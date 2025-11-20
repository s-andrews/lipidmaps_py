"""
Tests for DataValidator class
"""
import unittest
import tempfile
import logging
from pathlib import Path

from lipidmaps.data.ingestion import CSVIngestion, RawDataFrame
from lipidmaps.data.validation import DataValidator, ValidationReport, IssueSeverity


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


class TestDataValidator(unittest.TestCase):
    """Test DataValidator class for data quality checks."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_data_dir = Path(__file__).parent / "inputs"
        self.validator = DataValidator()
        self.ingestion = CSVIngestion()

    def test_validate_good_data(self):
        """Test validation on a good dataset."""
        test_file = self.test_data_dir / "small_demo.csv"
        raw_df = self.ingestion.read_csv(test_file)
        
        report = self.validator.validate(raw_df)
        
        self.assertIsInstance(report, ValidationReport)
        self.assertIsInstance(report.issues, list)
        self.assertIn('total_rows', report.summary)

    def test_empty_dataset_validation(self):
        """Test validation of empty dataset."""
        # Create empty CSV
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid,sample1\n")  # Header only
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = self.validator.validate(raw_df)
            
            # Should have critical issue for empty dataset
            self.assertFalse(report.passed)
            critical_issues = report.get_issues_by_severity(IssueSeverity.CRITICAL)
            self.assertGreater(len(critical_issues), 0)
        finally:
            temp_path.unlink()

    def test_missing_values_detection(self):
        """Test detection of missing values."""
        # Create CSV with missing values
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid,sample1,sample2\n")
            f.write("PC(16:0/18:1),123.45,\n")
            f.write("TAG(16:0/18:1/18:2),,678.90\n")
            f.write("CE(18:1),,\n")
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = self.validator.validate(raw_df)
            
            # Check for missing data issues
            missing_issues = report.get_issues_by_category('missing_data')
            self.assertGreater(len(missing_issues), 0)
        finally:
            temp_path.unlink()

    def test_non_numeric_values_detection(self):
        """Test detection of non-numeric values in sample columns."""
        # Create CSV with non-numeric values
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid,sample1,sample2\n")
            f.write("PC(16:0/18:1),123.45,abc\n")
            f.write("TAG(16:0/18:1/18:2),xyz,678.90\n")
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = self.validator.validate(raw_df)
            
            # Check for invalid value issues
            invalid_issues = report.get_issues_by_category('invalid_value')
            self.assertGreater(len(invalid_issues), 0)
            self.assertFalse(report.passed)
        finally:
            temp_path.unlink()

    def test_duplicate_lipid_names_detection(self):
        """Test detection of duplicate lipid names."""
        # Create CSV with duplicates
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid,sample1,sample2\n")
            f.write("PC(16:0/18:1),123.45,234.56\n")
            f.write("PC(16:0/18:1),678.90,789.01\n")
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = self.validator.validate(raw_df)
            
            # Check for consistency issues
            consistency_issues = report.get_issues_by_category('data_consistency')
            self.assertGreater(len(consistency_issues), 0)
        finally:
            temp_path.unlink()

    def test_negative_values_detection(self):
        """Test detection of negative values."""
        # Create CSV with negative values
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid,sample1,sample2\n")
            f.write("PC(16:0/18:1),-123.45,234.56\n")
            f.write("TAG(16:0/18:1/18:2),678.90,-789.01\n")
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = self.validator.validate(raw_df)
            
            # Check for negative value warnings
            invalid_issues = report.get_issues_by_category('invalid_value')
            self.assertGreater(len(invalid_issues), 0)
        finally:
            temp_path.unlink()

    def test_validation_report_properties(self):
        """Test ValidationReport properties and methods."""
        test_file = self.test_data_dir / "small_demo.csv"
        raw_df = self.ingestion.read_csv(test_file)
        
        report = self.validator.validate(raw_df)
        
        # Test properties
        self.assertIsInstance(report.passed, bool)
        self.assertIsInstance(report.has_warnings, bool)
        
        # Test methods
        warnings = report.get_issues_by_severity(IssueSeverity.WARNING)
        self.assertIsInstance(warnings, list)
        
        report_dict = report.to_dict()
        self.assertIsInstance(report_dict, dict)
        self.assertIn('passed', report_dict)
        self.assertIn('issues', report_dict)

    def test_insufficient_columns(self):
        """Test validation with insufficient columns."""
        # Create CSV with only one column
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid\n")
            f.write("PC(16:0/18:1)\n")
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = self.validator.validate(raw_df)
            
            # Should have critical or error issue
            self.assertFalse(report.passed)
            structure_issues = report.get_issues_by_category('structure')
            self.assertGreater(len(structure_issues), 0)
        finally:
            temp_path.unlink()

    def test_custom_validator_config(self):
        """Test validator with custom configuration."""
        validator = DataValidator(
            min_samples=3,
            min_lipids=5,
            max_missing_percent=10.0
        )
        
        # Create CSV that violates these limits
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("lipid,sample1\n")  # Only 1 sample (< 3)
            f.write("PC(16:0/18:1),123.45\n")
            temp_path = Path(f.name)
        
        try:
            raw_df = self.ingestion.read_csv(temp_path)
            report = validator.validate(raw_df)
            
            # Should have error for insufficient samples
            structure_issues = report.get_issues_by_category('structure')
            sample_issue = [i for i in structure_issues if 'sample' in i.message.lower()]
            self.assertGreater(len(sample_issue), 0)
        finally:
            temp_path.unlink()

    def test_validation_summary(self):
        """Test that validation summary contains expected keys."""
        test_file = self.test_data_dir / "small_demo.csv"
        raw_df = self.ingestion.read_csv(test_file)
        
        report = self.validator.validate(raw_df)
        
        expected_keys = [
            'total_rows',
            'total_columns',
            'sample_columns',
            'data_completeness_percent',
            'validation_passed'
        ]
        
        for key in expected_keys:
            self.assertIn(key, report.summary)


if __name__ == "__main__":
    unittest.main()
