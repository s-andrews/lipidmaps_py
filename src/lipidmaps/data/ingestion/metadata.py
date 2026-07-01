"""
Metadata extraction for CSV and tabular files.

This class provides methods to extract, validate, and standardize metadata from CSV files, following industry standards (e.g., ISA-Tab, MIAPE, FAIR principles).
"""

from typing import Dict, Any, Optional, List
from pathlib import Path
import pandas as pd
from pydantic import BaseModel

class MetadataExtractor(BaseModel):
    """
    Extracts and standardizes metadata from CSV/tabular files.
    - Supports summary statistics, column types and missing values.
    - Can be extended for ISA-Tab, MIAPE, or FAIR metadata standards.
    """
    file_path: str
    file_name: Optional[str] = None
    file_size_bytes: Optional[int] = None
    file_resolved_path: Optional[str] = None
    source: Optional[str] = "user-uploaded"
    num_rows: Optional[int] = None
    num_columns: Optional[int] = None
    columns: Optional[List[str]] = None
    column_types: Optional[Dict[str, str]] = None
    missing_values: Optional[Dict[str, int]] = None
    preview: Optional[List[Dict[str, Any]]] = None

    def extract_basic_metadata(self):
        path = Path(self.file_path)
        self.file_name = path.name
        self.file_size_bytes = path.stat().st_size if path.exists() else None
        self.file_resolved_path = str(path.resolve())
        self.source = self.source or "user-uploaded"

    def extract_from_dataframe(self, df: pd.DataFrame):
        self.num_rows = df.shape[0]
        self.num_columns = df.shape[1]
        self.columns = list(df.columns)
        self.column_types = {col: str(df[col].dtype) for col in df.columns}
        self.missing_values = {col: int(df[col].isnull().sum()) for col in df.columns}
        self.preview = df.head(5).to_dict(orient="records")

    def to_json(self) -> str:
        import json
        return json.dumps(self.model_dump(), indent=2)

    def get_metadata(self) -> Dict[str, Any]:
        return self.model_dump()

    # Optionally: add ISA-Tab, MIAPE, or FAIR metadata methods here
    # def extract_isa_tab(self, ...):
    #     pass
    # def extract_fair_metadata(self, ...):
    #     pass
