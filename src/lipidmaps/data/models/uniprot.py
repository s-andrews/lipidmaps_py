"""Helpers for retrieving reaction-linked gene symbols from UniProt."""

from __future__ import annotations

import csv
import io
import logging
from typing import Dict, List

import requests
from pydantic import Field

from .base import LipidmapsBaseModel

logger = logging.getLogger(__name__)


class UniProtGeneRecord(LipidmapsBaseModel):
    """Minimal UniProt record used for reaction gene lookup."""

    accession: str = Field(default="")
    gene_primary: str = Field(default="")
    organism_id: str = Field(default="")
    ec: str = Field(default="")


class UniProtRheaClient(LipidmapsBaseModel):
    """Fetch reviewed UniProt entries associated with a Rhea reaction."""

    base_url: str = Field(default="https://rest.uniprot.org/uniprotkb/search")
    timeout: int = Field(default=30, ge=1, le=300)

    def fetch_gene_records(self, rhea_id: str) -> List[UniProtGeneRecord]:
        """Return reviewed UniProt records for a Rhea identifier."""

        response = requests.get(
            self.base_url,
            params={
                "query": f"rhea:{rhea_id} AND reviewed:true",
                "fields": "gene_primary,organism_id,ec,accession",
                "format": "tsv",
            },
            timeout=self.timeout,
        )
        response.raise_for_status()
        return self._parse_tsv(response.text)

    def fetch_gene_symbols(self, rhea_id: str) -> List[str]:
        """Return deduplicated primary gene symbols for a Rhea identifier."""

        ordered: List[str] = []
        seen = set()
        for record in self.fetch_gene_records(rhea_id):
            gene = record.gene_primary.strip()
            if gene and gene not in seen:
                ordered.append(gene)
                seen.add(gene)
        return ordered

    @staticmethod
    def _parse_tsv(payload: str) -> List[UniProtGeneRecord]:
        """Parse UniProt TSV search results."""

        if not payload.strip():
            return []

        reader = csv.DictReader(io.StringIO(payload), delimiter="\t")
        if not reader.fieldnames:
            return []

        header_map: Dict[str, str] = {}
        for header in reader.fieldnames:
            normalized = header.strip().lower()
            if normalized == "entry":
                header_map["accession"] = header
            elif "gene" in normalized and "primary" in normalized:
                header_map["gene_primary"] = header
            elif "organism" in normalized and "id" in normalized:
                header_map["organism_id"] = header
            elif normalized.startswith("ec"):
                header_map["ec"] = header

        records: List[UniProtGeneRecord] = []
        for row in reader:
            records.append(
                UniProtGeneRecord(
                    accession=(row.get(header_map.get("accession", ""), "") or "").strip(),
                    gene_primary=(row.get(header_map.get("gene_primary", ""), "") or "").strip(),
                    organism_id=(row.get(header_map.get("organism_id", ""), "") or "").strip(),
                    ec=(row.get(header_map.get("ec", ""), "") or "").strip(),
                )
            )

        logger.debug("Parsed %s UniProt records", len(records))
        return records
