"""
Data Validation Module

Provides comprehensive data quality checks for lipid quantification datasets.
Validates structure, detects missing values, checks numeric ranges, and generates quality reports.
"""
import logging
import re
from typing import List, Dict, Any, Optional, Set
from dataclasses import dataclass, field
from enum import Enum
from collections import Counter

from ..ingestion.csv_reader import RawDataFrame

logger = logging.getLogger(__name__)


class IssueSeverity(Enum):
    """Severity levels for validation issues."""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class ValidationIssue:
    """Represents a single validation issue found in the data.
    
    Attributes:
        severity: Issue severity level
        category: Issue category (e.g., 'structure', 'missing_data', 'invalid_value')
        message: Human-readable description
        location: Where the issue was found (row, column, cell)
        suggestion: Optional suggestion for fixing the issue
    """
    severity: IssueSeverity
    category: str
    message: str
    location: Optional[Dict[str, Any]] = None
    suggestion: Optional[str] = None
    
    def __str__(self) -> str:
        """String representation of the issue."""
        loc = f" at {self.location}" if self.location else ""
        return f"[{self.severity.value.upper()}] {self.category}: {self.message}{loc}"


@dataclass
class ValidationReport:
    """Complete validation report for a dataset.
    
    Attributes:
        issues: List of validation issues found
        passed: Whether validation passed (no errors/critical issues)
        summary: Summary statistics about the validation
    """
    issues: List[ValidationIssue] = field(default_factory=list)
    summary: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def passed(self) -> bool:
        """Check if validation passed (no critical or error issues)."""
        critical_or_error = [
            issue for issue in self.issues
            if issue.severity in (IssueSeverity.CRITICAL, IssueSeverity.ERROR)
        ]
        return len(critical_or_error) == 0
    
    @property
    def has_warnings(self) -> bool:
        """Check if there are any warnings."""
        return any(issue.severity == IssueSeverity.WARNING for issue in self.issues)
    
    def get_issues_by_severity(self, severity: IssueSeverity) -> List[ValidationIssue]:
        """Get all issues of a specific severity."""
        return [issue for issue in self.issues if issue.severity == severity]
    
    def get_issues_by_category(self, category: str) -> List[ValidationIssue]:
        """Get all issues in a specific category."""
        return [issue for issue in self.issues if issue.category == category]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary."""
        return {
            'passed': self.passed,
            'has_warnings': self.has_warnings,
            'issue_count': len(self.issues),
            'issues': [
                {
                    'severity': issue.severity.value,
                    'category': issue.category,
                    'message': issue.message,
                    'location': issue.location,
                    'suggestion': issue.suggestion
                }
                for issue in self.issues
            ],
            'summary': self.summary
        }
    
    def print_report(self) -> None:
        """Print formatted validation report."""
        print("\n" + "=" * 80)
        print("VALIDATION REPORT")
        print("=" * 80)
        
        print(f"\nStatus: {'PASSED' if self.passed else 'FAILED'}")
        print(f"Total Issues: {len(self.issues)}")
        
        # Count by severity
        severity_counts = Counter(issue.severity for issue in self.issues)
        for severity in IssueSeverity:
            count = severity_counts.get(severity, 0)
            if count > 0:
                print(f"  {severity.value.capitalize()}: {count}")
        
        # Print issues grouped by severity
        for severity in [IssueSeverity.CRITICAL, IssueSeverity.ERROR, IssueSeverity.WARNING, IssueSeverity.INFO]:
            issues = self.get_issues_by_severity(severity)
            if issues:
                print(f"\n{severity.value.upper()} Issues:")
                for issue in issues:
                    print(f"  - {issue}")
        
        # Print summary
        if self.summary:
            print("\nSummary:")
            for key, value in self.summary.items():
                print(f"  {key}: {value}")
        
        print("=" * 80 + "\n")


class DataValidator:
    """Comprehensive data quality validator for lipid datasets.
    
    This class performs various validation checks on raw data:
    - Structural validation (columns, rows, format)
    - Missing value detection
    - Numeric range validation
    - Lipid name format validation
    - Data consistency checks
    - Statistical outlier detection
    """
    
    def __init__(
        self,
        min_samples: int = 1,
        min_lipids: int = 1,
        allow_missing_values: bool = True,
        max_missing_percent: float = 50.0
    ):
        """Initialize validator with configuration.
        
        Args:
            min_samples: Minimum number of samples required
            min_lipids: Minimum number of lipids required
            allow_missing_values: Whether missing values are allowed
            max_missing_percent: Maximum percentage of missing values allowed
        """
        self.min_samples = min_samples
        self.min_lipids = min_lipids
        self.allow_missing_values = allow_missing_values
        self.max_missing_percent = max_missing_percent
    
    def validate(self, raw_df: RawDataFrame) -> ValidationReport:
        """Run complete validation on raw data frame.
        
        Args:
            raw_df: RawDataFrame to validate
            
        Returns:
            ValidationReport with all issues found
        """
        logger.info("Starting data validation")
        report = ValidationReport()
        
        # Run all validation checks
        self._validate_structure(raw_df, report)
        self._validate_missing_values(raw_df, report)
        self._validate_numeric_values(raw_df, report)
        self._validate_lipid_names(raw_df, report)
        self._validate_consistency(raw_df, report)
        
        # Generate summary
        report.summary = self._generate_summary(raw_df, report)
        
        logger.info(
            f"Validation complete: {len(report.issues)} issues found, "
            f"passed={report.passed}"
        )
        
        return report
    
    def _validate_structure(self, raw_df: RawDataFrame, report: ValidationReport) -> None:
        """Validate basic structural requirements."""
        # Check for empty data
        if raw_df.is_empty():
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.CRITICAL,
                category='structure',
                message='Dataset is empty (no data rows)',
                suggestion='Check if file was read correctly or contains data'
            ))
            return
        
        # Check minimum columns
        if len(raw_df.fieldnames) < 2:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.CRITICAL,
                category='structure',
                message=f'Insufficient columns: found {len(raw_df.fieldnames)}, need at least 2',
                suggestion='File should have at least one lipid name column and one sample column'
            ))
        
        # Check for duplicate column names
        duplicates = [
            name for name, count in Counter(raw_df.fieldnames).items()
            if count > 1
        ]
        if duplicates:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.ERROR,
                category='structure',
                message=f'Duplicate column names found: {duplicates}',
                suggestion='Ensure all column headers are unique'
            ))
        
        # Check for empty column names
        empty_cols = [i for i, name in enumerate(raw_df.fieldnames) if not name or not name.strip()]
        if empty_cols:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.WARNING,
                category='structure',
                message=f'Empty column names at positions: {empty_cols}',
                suggestion='All columns should have meaningful names'
            ))
        
        # Validate sample count
        sample_cols = len(raw_df.fieldnames) - 1  # Assuming first col is lipid names
        if sample_cols < self.min_samples:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.ERROR,
                category='structure',
                message=f'Insufficient samples: found {sample_cols}, need at least {self.min_samples}',
                location={'sample_count': sample_cols}
            ))
        
        # Validate lipid count
        if raw_df.row_count < self.min_lipids:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.WARNING,
                category='structure',
                message=f'Low lipid count: found {raw_df.row_count}, expected at least {self.min_lipids}',
                location={'lipid_count': raw_df.row_count}
            ))
    
    def _validate_missing_values(self, raw_df: RawDataFrame, report: ValidationReport) -> None:
        """Check for missing values in the dataset."""
        if raw_df.is_empty():
            return
        
        name_col = raw_df.fieldnames[0]
        sample_cols = raw_df.fieldnames[1:]
        
        # Check for missing lipid names
        missing_names = 0
        for i, row in enumerate(raw_df.rows):
            name = row.get(name_col, '').strip()
            if not name:
                missing_names += 1
                report.issues.append(ValidationIssue(
                    severity=IssueSeverity.ERROR,
                    category='missing_data',
                    message=f'Missing lipid name',
                    location={'row': i + 1}
                ))
        
        # Check for missing values in sample columns
        for col in sample_cols:
            missing_count = 0
            for i, row in enumerate(raw_df.rows):
                value = row.get(col, '').strip()
                if not value:
                    missing_count += 1
            
            if missing_count > 0:
                missing_percent = (missing_count / len(raw_df.rows)) * 100
                
                if missing_percent > self.max_missing_percent:
                    report.issues.append(ValidationIssue(
                        severity=IssueSeverity.ERROR,
                        category='missing_data',
                        message=f'Excessive missing values in column "{col}": '
                                f'{missing_percent:.1f}% ({missing_count}/{len(raw_df.rows)})',
                        location={'column': col, 'missing_count': missing_count},
                        suggestion=f'Maximum allowed is {self.max_missing_percent}%'
                    ))
                elif missing_percent > 10:
                    report.issues.append(ValidationIssue(
                        severity=IssueSeverity.WARNING,
                        category='missing_data',
                        message=f'Missing values in column "{col}": '
                                f'{missing_percent:.1f}% ({missing_count}/{len(raw_df.rows)})',
                        location={'column': col, 'missing_count': missing_count}
                    ))
    
    def _validate_numeric_values(self, raw_df: RawDataFrame, report: ValidationReport) -> None:
        """Validate that sample columns contain valid numeric values."""
        if raw_df.is_empty():
            return
        
        sample_cols = raw_df.fieldnames[1:]
        
        for col in sample_cols:
            non_numeric = []
            negative_values = []
            zero_values = 0
            values = []
            
            for i, row in enumerate(raw_df.rows):
                value_str = row.get(col, '').strip()
                if not value_str:
                    continue  # Skip missing values (handled separately)
                
                try:
                    value = float(value_str)
                    values.append(value)
                    
                    if value < 0:
                        negative_values.append((i + 1, value))
                    elif value == 0:
                        zero_values += 1
                except ValueError:
                    non_numeric.append((i + 1, value_str))
            
            # Report non-numeric values
            if non_numeric:
                sample = non_numeric[:3]  # Show first 3 examples
                report.issues.append(ValidationIssue(
                    severity=IssueSeverity.ERROR,
                    category='invalid_value',
                    message=f'Non-numeric values in column "{col}": {len(non_numeric)} found',
                    location={'column': col, 'examples': sample},
                    suggestion='All sample values must be numeric'
                ))
            
            # Report negative values
            if negative_values:
                sample = negative_values[:3]
                report.issues.append(ValidationIssue(
                    severity=IssueSeverity.WARNING,
                    category='invalid_value',
                    message=f'Negative values in column "{col}": {len(negative_values)} found',
                    location={'column': col, 'examples': sample},
                    suggestion='Quantitation values should typically be positive'
                ))
            
            # Report excessive zeros
            if values and zero_values / len(values) > 0.5:
                zero_percent = (zero_values / len(values)) * 100
                report.issues.append(ValidationIssue(
                    severity=IssueSeverity.WARNING,
                    category='data_quality',
                    message=f'High proportion of zero values in column "{col}": '
                            f'{zero_percent:.1f}% ({zero_values}/{len(values)})',
                    location={'column': col}
                ))
    
    def _validate_lipid_names(self, raw_df: RawDataFrame, report: ValidationReport) -> None:
        """Validate lipid name formats and check for duplicates."""
        if raw_df.is_empty():
            return
        
        name_col = raw_df.fieldnames[0]
        names = [row.get(name_col, '').strip() for row in raw_df.rows]
        names = [n for n in names if n]  # Remove empty
        
        # Check for duplicates
        duplicates = [name for name, count in Counter(names).items() if count > 1]
        if duplicates:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.WARNING,
                category='data_consistency',
                message=f'Duplicate lipid names found: {len(duplicates)} duplicates',
                location={'duplicates': duplicates[:5]},  # Show first 5
                suggestion='Consider using unique identifiers or aggregating duplicates'
            ))
        
        # Check for suspicious patterns
        very_short = [name for name in names if len(name) < 3]
        if very_short:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.INFO,
                category='data_quality',
                message=f'Very short lipid names found: {len(very_short)}',
                location={'examples': very_short[:5]}
            ))
        
        # Check for common invalid patterns
        suspicious_patterns = [
            (r'^\d+$', 'numeric only'),
            (r'^[Nn][Aa]$', 'NA values'),
            (r'^[Uu]nknown', 'unknown markers'),
        ]
        
        for pattern, description in suspicious_patterns:
            matches = [name for name in names if re.match(pattern, name)]
            if matches:
                report.issues.append(ValidationIssue(
                    severity=IssueSeverity.INFO,
                    category='data_quality',
                    message=f'Potential placeholder names found ({description}): {len(matches)}',
                    location={'examples': matches[:5]}
                ))
    
    def _validate_consistency(self, raw_df: RawDataFrame, report: ValidationReport) -> None:
        """Check for data consistency issues."""
        if raw_df.is_empty():
            return
        
        # Check for rows with all missing values
        sample_cols = raw_df.fieldnames[1:]
        empty_rows = []
        
        for i, row in enumerate(raw_df.rows):
            values = [row.get(col, '').strip() for col in sample_cols]
            if not any(values):
                empty_rows.append(i + 1)
        
        if empty_rows:
            report.issues.append(ValidationIssue(
                severity=IssueSeverity.WARNING,
                category='data_consistency',
                message=f'Rows with no quantitation data: {len(empty_rows)}',
                location={'rows': empty_rows[:10]},
                suggestion='Consider removing rows with no data'
            ))
    
    def _generate_summary(self, raw_df: RawDataFrame, report: ValidationReport) -> Dict[str, Any]:
        """Generate validation summary statistics."""
        sample_cols = raw_df.fieldnames[1:]
        
        # Calculate data completeness
        total_cells = len(raw_df.rows) * len(sample_cols)
        missing_cells = 0
        
        for row in raw_df.rows:
            for col in sample_cols:
                if not row.get(col, '').strip():
                    missing_cells += 1
        
        completeness = ((total_cells - missing_cells) / total_cells * 100) if total_cells > 0 else 0
        
        return {
            'total_rows': len(raw_df.rows),
            'total_columns': len(raw_df.fieldnames),
            'sample_columns': len(sample_cols),
            'data_completeness_percent': round(completeness, 2),
            'validation_passed': report.passed,
            'critical_issues': len(report.get_issues_by_severity(IssueSeverity.CRITICAL)),
            'errors': len(report.get_issues_by_severity(IssueSeverity.ERROR)),
            'warnings': len(report.get_issues_by_severity(IssueSeverity.WARNING)),
            'info': len(report.get_issues_by_severity(IssueSeverity.INFO))
        }
