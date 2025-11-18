"""
CSV Ingestion Module

Handles reading CSV files in various formats and provides basic structural validation.
Supports standard CSV, MS-DIAL, and auto-detection of common lipidomics formats.
"""
import csv
import logging
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional, Union
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class CSVFormat(Enum):
    """Supported CSV formats."""
    STANDARD = "standard"
    MSDIAL = "msdial"
    AUTO = "auto"


@dataclass
class RawDataFrame:
    """Container for raw CSV data before processing.
    
    Attributes:
        rows: List of dictionaries, one per data row
        fieldnames: List of column headers
        format_type: Detected or specified format type
        metadata: Additional metadata about the file
    """
    rows: List[Dict[str, str]]
    fieldnames: List[str]
    format_type: CSVFormat
    metadata: Dict[str, Any]
    
    @property
    def row_count(self) -> int:
        """Return number of data rows."""
        return len(self.rows)
    
    @property
    def column_count(self) -> int:
        """Return number of columns."""
        return len(self.fieldnames)
    
    def is_empty(self) -> bool:
        """Check if data frame has no rows."""
        return len(self.rows) == 0


class CSVIngestion:
    """Handle CSV file ingestion with format detection and basic validation.
    
    This class is responsible for:
    - Reading CSV files with various delimiters and encodings
    - Detecting common lipidomics data formats
    - Basic structural validation
    - Converting raw CSV to RawDataFrame objects
    
    Does NOT handle:
    - Data quality validation (use DataValidator)
    - Business logic transformations
    - API calls for annotation
    """
    
    SUPPORTED_DELIMITERS = [',', '\t', ';', '|']
    SUPPORTED_ENCODINGS = ['utf-8', 'latin-1', 'iso-8859-1']
    
    def __init__(self, delimiter: Optional[str] = None, encoding: str = 'utf-8'):
        """Initialize CSV reader.
        
        Args:
            delimiter: CSV delimiter (default: auto-detect)
            encoding: File encoding (default: utf-8)
        """
        self.delimiter = delimiter
        self.encoding = encoding
    
    def read_csv(
        self,
        path: Union[str, Path],
        format_type: CSVFormat = CSVFormat.AUTO
    ) -> RawDataFrame:
        """Read CSV file and return raw data frame.
        
        Args:
            path: Path to CSV file
            format_type: Format type (AUTO, STANDARD, MSDIAL)
            
        Returns:
            RawDataFrame containing raw data
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file cannot be parsed
        """
        path = Path(path)
        logger.info(f"Reading CSV file: {path}")
        
        if not path.exists():
            raise FileNotFoundError(f"CSV file not found: {path}")
        
        # Detect format if AUTO
        if format_type == CSVFormat.AUTO:
            format_type = self.detect_format(path)
            logger.info(f"Detected format: {format_type.value}")
        
        # Read based on format
        if format_type == CSVFormat.MSDIAL:
            return self.read_msdial(path)
        else:
            return self.read_standard_csv(path)
    
    def read_standard_csv(self, path: Path) -> RawDataFrame:
        """Read standard CSV file.
        
        Args:
            path: Path to CSV file
            
        Returns:
            RawDataFrame with parsed data
        """
        # Auto-detect delimiter if not specified
        delimiter = self.delimiter or self._detect_delimiter(path)
        
        rows = []
        fieldnames = []
        
        try:
            with path.open('r', encoding=self.encoding, newline='') as fh:
                reader = csv.DictReader(fh, delimiter=delimiter)
                fieldnames = reader.fieldnames or []
                rows = list(reader)
        except UnicodeDecodeError:
            # Try alternative encoding
            logger.warning(f"Failed to decode with {self.encoding}, trying latin-1")
            with path.open('r', encoding='latin-1', newline='') as fh:
                reader = csv.DictReader(fh, delimiter=delimiter)
                fieldnames = reader.fieldnames or []
                rows = list(reader)
        
        metadata = {
            'source_file': str(path),
            'delimiter': delimiter,
            'encoding': self.encoding,
            'file_size_bytes': path.stat().st_size
        }
        
        logger.info(
            f"Read {len(rows)} rows, {len(fieldnames)} columns from {path.name}"
        )
        
        return RawDataFrame(
            rows=rows,
            fieldnames=fieldnames,
            format_type=CSVFormat.STANDARD,
            metadata=metadata
        )
    
    def read_msdial(self, path: Path) -> RawDataFrame:
        """Read MS-DIAL formatted CSV file.
        
        MS-DIAL files may have special headers or metadata rows.
        For now, this treats them as standard CSV but can be extended.
        
        Args:
            path: Path to MS-DIAL file
            
        Returns:
            RawDataFrame with parsed data
        """
        logger.info(f"Reading MS-DIAL format from {path}")
        
        # TODO: Add MS-DIAL specific parsing logic
        # - Skip metadata rows
        # - Handle special column names
        # - Parse retention time, mass, etc.
        
        # For now, read as standard CSV with tab delimiter
        original_delimiter = self.delimiter
        self.delimiter = '\t'
        
        result = self.read_standard_csv(path)
        result.format_type = CSVFormat.MSDIAL
        result.metadata['format_notes'] = 'MS-DIAL format (basic parsing)'
        
        self.delimiter = original_delimiter
        
        return result
    
    def detect_format(self, path: Path) -> CSVFormat:
        """Detect CSV format by inspecting file contents.
        
        Args:
            path: Path to CSV file
            
        Returns:
            Detected CSVFormat
        """
        # Read first few lines to detect format
        with path.open('r', encoding=self.encoding) as fh:
            lines = [fh.readline() for _ in range(5)]
        
        # Check for MS-DIAL indicators
        header = lines[0].lower() if lines else ''
        
        msdial_indicators = [
            'alignment id',
            'average rt',
            'metabolite name',
            'ms-dial'
        ]
        
        if any(indicator in header for indicator in msdial_indicators):
            logger.debug("Detected MS-DIAL format")
            return CSVFormat.MSDIAL
        
        logger.debug("Defaulting to STANDARD format")
        return CSVFormat.STANDARD
    
    def _detect_delimiter(self, path: Path) -> str:
        """Detect CSV delimiter by analyzing first few lines.
        
        Args:
            path: Path to CSV file
            
        Returns:
            Most likely delimiter character
        """
        with path.open('r', encoding=self.encoding) as fh:
            sample = fh.read(8192)  # Read first 8KB
        
        # Count occurrences of each delimiter
        delimiter_counts = {
            delim: sample.count(delim)
            for delim in self.SUPPORTED_DELIMITERS
        }
        
        # Return delimiter with highest count
        detected = max(delimiter_counts, key=delimiter_counts.get)
        logger.debug(f"Detected delimiter: {repr(detected)}")
        
        return detected
    
    def read_batch(
        self,
        paths: List[Union[str, Path]],
        format_type: CSVFormat = CSVFormat.AUTO
    ) -> List[RawDataFrame]:
        """Read multiple CSV files.
        
        Args:
            paths: List of file paths
            format_type: Format type for all files
            
        Returns:
            List of RawDataFrame objects
        """
        results = []
        for path in paths:
            try:
                df = self.read_csv(path, format_type)
                results.append(df)
            except Exception as e:
                logger.error(f"Failed to read {path}: {e}")
                # Continue with other files
        
        logger.info(f"Successfully read {len(results)} of {len(paths)} files")
        return results
    
    def get_column_info(self, raw_df: RawDataFrame) -> Dict[str, Any]:
        """Get information about columns in the data frame.
        
        Args:
            raw_df: RawDataFrame to analyze
            
        Returns:
            Dictionary with column information
        """
        info = {
            'column_count': len(raw_df.fieldnames),
            'columns': raw_df.fieldnames,
            'empty_columns': [],
            'column_types': {}
        }
        
        # Check for empty columns
        for col in raw_df.fieldnames:
            values = [row.get(col, '').strip() for row in raw_df.rows]
            non_empty = [v for v in values if v]
            
            if not non_empty:
                info['empty_columns'].append(col)
            else:
                # Try to determine column type
                info['column_types'][col] = self._guess_column_type(non_empty)
        
        return info
    
    def _guess_column_type(self, values: List[str]) -> str:
        """Guess the data type of a column based on values.
        
        Args:
            values: List of non-empty string values
            
        Returns:
            Guessed type: 'numeric', 'text', 'identifier', 'mixed'
        """
        if not values:
            return 'empty'
        
        # Sample up to 100 values
        sample = values[:100]
        
        numeric_count = 0
        for val in sample:
            try:
                float(val)
                numeric_count += 1
            except ValueError:
                pass
        
        # If >80% numeric, consider it numeric
        if numeric_count / len(sample) > 0.8:
            return 'numeric'
        
        # Check if looks like identifiers (short, alphanumeric)
        avg_length = sum(len(v) for v in sample) / len(sample)
        if avg_length < 20 and all(v.replace('_', '').replace('-', '').isalnum() for v in sample[:10]):
            return 'identifier'
        
        return 'text'
