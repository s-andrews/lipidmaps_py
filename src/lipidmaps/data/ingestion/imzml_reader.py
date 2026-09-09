"""Read mass-spectrometry-imaging (MSI) ``.imzML`` files into per-pixel ion series.

imzML stores only m/z + intensity per pixel -- no molecule names. Names come from
a companion *annotation* table (e.g. a METASPACE export) giving, per ion, a lipid
name and its target m/z. This reader parses the imzML with ``pyimzml``, then for
each annotated ion extracts the matching peak intensity in every pixel, yielding a
neutral intermediate the importer turns into ``SampleMetadata`` + ``QuantifiedLipid``.

``pyimzml`` is an optional dependency; import errors are raised only when a reader
is actually constructed, so the rest of the package works without it installed.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

logger = logging.getLogger(__name__)


class IonAnnotation(BaseModel):
    """One annotated ion: a lipid name plus the m/z to extract for it."""

    name: str
    mz: float
    adduct: Optional[str] = None
    formula: Optional[str] = None
    lm_id: Optional[str] = None


class ImzMLIngestionResult(BaseModel):
    """Neutral intermediate produced by :class:`ImzMLIngestion`.

    - ``pixels``: ordered ``(pixel_name, x, y, z)`` tuples, one per spectrum kept.
    - ``ion_values``: ``{ion_name: {pixel_name: intensity}}``.
    - ``annotations``: the ions that were extracted (for provenance/reporting).
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    pixels: List[Tuple[str, int, int, int]] = Field(default_factory=list)
    ion_values: Dict[str, Dict[str, float]] = Field(default_factory=dict)
    annotations: List[IonAnnotation] = Field(default_factory=list)

    @property
    def pixel_count(self) -> int:
        return len(self.pixels)


def pixel_name(x: int, y: int, z: int = 0) -> str:
    """Deterministic, sortable sample name for a pixel coordinate."""
    return f"px_x{int(x):04d}_y{int(y):04d}_z{int(z):02d}"


def parse_annotation_csv(path: str | Path) -> List[IonAnnotation]:
    """Parse a companion annotation CSV into :class:`IonAnnotation` rows.

    Expects at least ``name`` and ``mz`` columns (case-insensitive); ``adduct``,
    ``formula`` and ``lm_id`` are optional. Rows without a usable name+mz are
    skipped. Column matching is lenient so METASPACE-style headers can be mapped
    by the demo fetcher into this simple schema.
    """
    path = Path(path)
    annotations: List[IonAnnotation] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        # Case-insensitive header lookup.
        field_map = {(name or "").strip().lower(): name for name in (reader.fieldnames or [])}

        def col(*candidates: str) -> Optional[str]:
            for cand in candidates:
                if cand in field_map:
                    return field_map[cand]
            return None

        name_col = col("name", "molecule", "molecule_name", "moleculenames")
        mz_col = col("mz", "m/z", "mz_value")
        adduct_col = col("adduct", "ion")
        formula_col = col("formula", "sumformula", "ion_formula")
        lmid_col = col("lm_id", "lmid", "moleculeids", "molecule_ids")

        if name_col is None or mz_col is None:
            raise ValueError(
                f"Annotation CSV {path} must have 'name' and 'mz' columns; "
                f"found {reader.fieldnames}"
            )

        for row in reader:
            name = (row.get(name_col) or "").strip()
            raw_mz = (row.get(mz_col) or "").strip()
            if not name or not raw_mz:
                continue
            try:
                mz = float(raw_mz)
            except ValueError:
                continue
            annotations.append(
                IonAnnotation(
                    name=name,
                    mz=mz,
                    adduct=(row.get(adduct_col) or "").strip() or None if adduct_col else None,
                    formula=(row.get(formula_col) or "").strip() or None if formula_col else None,
                    lm_id=(row.get(lmid_col) or "").strip() or None if lmid_col else None,
                )
            )
    return annotations


def _match_intensity(
    mzs: np.ndarray, intensities: np.ndarray, target_mz: float, tol_ppm: float
) -> float:
    """Intensity of the peak nearest ``target_mz`` within ``tol_ppm``, else 0.0.

    Assumes ``mzs`` is sorted ascending (true for imzML m/z arrays). Absent signal
    is returned as 0.0, which is the natural "dark pixel" value for an ion image.
    """
    n = len(mzs)
    if n == 0:
        return 0.0
    idx = int(np.searchsorted(mzs, target_mz))
    best_j: Optional[int] = None
    best_diff: Optional[float] = None
    for j in (idx - 1, idx):
        if 0 <= j < n:
            diff = abs(float(mzs[j]) - target_mz)
            if best_diff is None or diff < best_diff:
                best_diff = diff
                best_j = j
    if best_j is None or best_diff is None:
        return 0.0
    if (best_diff / target_mz) * 1.0e6 <= tol_ppm:
        return float(intensities[best_j])
    return 0.0


class ImzMLIngestion:
    """Reads an imzML file and extracts per-pixel intensities for annotated ions."""

    def __init__(self) -> None:
        try:
            from pyimzml.ImzMLParser import ImzMLParser  # noqa: F401
        except ImportError as exc:  # pragma: no cover - exercised only without extra
            raise ImportError(
                "Reading imzML requires the optional 'pyimzml' dependency. "
                "Install it with: pip install lipidmaps_py[msi]"
            ) from exc

    def read(
        self,
        imzml_path: str | Path,
        annotations: List[IonAnnotation] | str | Path,
        mz_tolerance_ppm: float = 10.0,
        bbox: Optional[Tuple[int, int, int, int]] = None,
    ) -> ImzMLIngestionResult:
        """Parse ``imzml_path`` and extract intensities for each annotated ion.

        Args:
            imzml_path: Path to the ``.imzML`` (its ``.ibd`` must sit alongside).
            annotations: ``IonAnnotation`` list, or a path to an annotation CSV.
            mz_tolerance_ppm: Match window for locating an ion's peak in a pixel.
            bbox: Optional ``(x_min, y_min, x_max, y_max)`` crop (inclusive) to
                subsample a spatial region; ``None`` keeps all pixels.
        """
        from pyimzml.ImzMLParser import ImzMLParser

        if isinstance(annotations, (str, Path)):
            annotations = parse_annotation_csv(annotations)
        if not annotations:
            raise ValueError("No ion annotations provided; cannot extract ion images.")

        parser = ImzMLParser(str(imzml_path))

        result = ImzMLIngestionResult(annotations=list(annotations))
        for ann in annotations:
            result.ion_values[ann.name] = {}

        kept = 0
        for idx, coord in enumerate(parser.coordinates):
            x, y, z = (coord + (0, 0, 0))[:3]
            if bbox is not None:
                x_min, y_min, x_max, y_max = bbox
                if not (x_min <= x <= x_max and y_min <= y <= y_max):
                    continue
            name = pixel_name(x, y, z)
            mzs, intensities = parser.getspectrum(idx)
            mzs = np.asarray(mzs, dtype=float)
            intensities = np.asarray(intensities, dtype=float)
            result.pixels.append((name, int(x), int(y), int(z)))
            for ann in annotations:
                result.ion_values[ann.name][name] = _match_intensity(
                    mzs, intensities, ann.mz, mz_tolerance_ppm
                )
            kept += 1

        logger.info(
            "imzML %s: kept %d pixels, extracted %d ions (tol=%.1f ppm)",
            imzml_path,
            kept,
            len(annotations),
            mz_tolerance_ppm,
        )
        return result
