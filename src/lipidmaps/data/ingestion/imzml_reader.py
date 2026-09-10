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
    """One annotated ion: a display name plus the m/z to extract for it.

    For a supplied annotation table ``name`` is the lipid name. For auto-annotation
    from a molecule database, ions are labeled at ``formula [adduct]`` level and the
    isomeric molecule names are kept in ``candidates`` (each ``{name, id}``) so an
    lm_id can be resolved from the names downstream.
    """

    name: str
    mz: float
    adduct: Optional[str] = None
    formula: Optional[str] = None
    lm_id: Optional[str] = None
    theoretical_mz: Optional[float] = None
    ppm_error: Optional[float] = None
    candidates: List[Dict[str, str]] = Field(default_factory=list)


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
    # (x, y) pixel size in micrometres from imzML metadata, when available.
    pixel_size: Optional[Tuple[float, float]] = None

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
        annotations: Optional[List[IonAnnotation] | str | Path] = None,
        mz_tolerance_ppm: float = 5.0,
        bbox: Optional[Tuple[int, int, int, int]] = None,
        database=None,
        adducts=None,
        top_n_peaks: Optional[int] = None,
        stride: int = 1,
        progress=None,
    ) -> ImzMLIngestionResult:
        """Parse ``imzml_path`` and extract per-pixel intensities for each ion.

        Provide either an explicit annotation set or a molecule ``database`` for
        auto-annotation (matching observed peaks by formula+adduct). An explicit
        ``annotations`` argument takes precedence over the database.

        Reads spectra lazily and reports progress via the optional ``progress``
        callback ``(phase, done, total)``. For continuous files (shared m/z axis)
        the database is matched against the axis directly -- no extra pass -- and
        extraction reads fixed channel indices; processed files discover peaks from a
        pixel sample. ``stride`` keeps every Nth selected pixel to bound work/memory.

        Args:
            imzml_path: Path to the ``.imzML`` (its ``.ibd`` must sit alongside).
            annotations: ``IonAnnotation`` list or a path to an annotation CSV.
                If ``None``, auto-annotate from ``database``.
            mz_tolerance_ppm: Match window (ppm) for peak matching.
            bbox: Optional ``(x_min, y_min, x_max, y_max)`` crop (inclusive).
            database: A ``MetaboliteDatabase`` or a path to a CoreMetabolome-style
                table, used only when ``annotations`` is ``None``.
            adducts: Adduct list for auto-annotation; defaults to a polarity-based set.
            top_n_peaks: Cap the number of observed peaks/ions (ranked by intensity).
            stride: Keep every Nth selected pixel (>=1) to subsample large datasets.
            progress: Optional callback ``(phase: str, done: int, total: int)``.
        """
        from pyimzml.ImzMLParser import ImzMLParser

        self._report(progress, "parse", 0, 1)
        parser = ImzMLParser(str(imzml_path))
        self._report(progress, "parse", 1, 1)

        selected = self._selected_indices(parser, bbox, stride)
        if not selected:
            raise ValueError("No pixels selected (check --bbox / --stride).")

        if isinstance(annotations, (str, Path)):
            annotations = parse_annotation_csv(annotations)
        if not annotations:
            annotations = self._auto_annotate(
                parser, selected, database, mz_tolerance_ppm, adducts, top_n_peaks, progress
            )
        if not annotations:
            raise ValueError(
                "No ions to extract: provide annotations, or a database that yields "
                "matches for the observed peaks."
            )

        # Fast path for continuous data: extract from fixed channel indices.
        continuous, axis = self._axis_if_continuous(parser, selected)
        channel_idx = (
            self._channel_indices(axis, annotations, mz_tolerance_ppm)
            if continuous and axis is not None else None
        )

        result = ImzMLIngestionResult(
            annotations=list(annotations), pixel_size=self._pixel_size(parser)
        )
        for ann in annotations:
            result.ion_values[ann.name] = {}

        total = len(selected)
        step = max(1, total // 50)
        for n, (idx, x, y, z) in enumerate(selected):
            name = pixel_name(x, y, z)
            mzs, intensities = parser.getspectrum(idx)
            intensities = np.asarray(intensities, dtype=float)
            result.pixels.append((name, x, y, z))
            if channel_idx is not None:
                for ann, ci in zip(annotations, channel_idx):
                    result.ion_values[ann.name][name] = (
                        float(intensities[ci]) if 0 <= ci < len(intensities) else 0.0
                    )
            else:
                mzs = np.asarray(mzs, dtype=float)
                for ann in annotations:
                    result.ion_values[ann.name][name] = _match_intensity(
                        mzs, intensities, ann.mz, mz_tolerance_ppm
                    )
            if n % step == 0 or n + 1 == total:
                self._report(progress, "extract", n + 1, total)

        logger.info(
            "imzML %s: %d pixels (stride=%d), %d ions (%s, tol=%.1f ppm)",
            imzml_path, total, stride, len(annotations),
            "continuous" if channel_idx is not None else "processed", mz_tolerance_ppm,
        )
        return result

    @staticmethod
    def _pixel_size(parser) -> Optional[Tuple[float, float]]:
        """(x, y) pixel size in µm from imzML metadata, or None if not recorded."""
        meta = getattr(parser, "imzmldict", {}) or {}
        x = meta.get("pixel size (x)") or meta.get("pixel size x") or meta.get("pixel size")
        y = meta.get("pixel size (y)") or meta.get("pixel size y") or x
        try:
            return (float(x), float(y)) if x is not None else None
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _report(progress, phase: str, done: int, total: int) -> None:
        if progress is not None:
            try:
                progress(phase, done, total)
            except Exception:  # progress must never break ingestion
                pass

    @staticmethod
    def _selected_indices(parser, bbox, stride):
        """(idx, x, y, z) for pixels kept after bbox crop and every-Nth stride."""
        stride = max(1, int(stride or 1))
        out = []
        kept = 0
        for i, coord in enumerate(parser.coordinates):
            x, y, z = (tuple(coord) + (0, 0, 0))[:3]
            if bbox is not None:
                x0, y0, x1, y1 = bbox
                if not (x0 <= x <= x1 and y0 <= y <= y1):
                    continue
            if kept % stride == 0:
                out.append((i, int(x), int(y), int(z)))
            kept += 1
        return out

    @staticmethod
    def _axis_if_continuous(parser, selected):
        """Return (is_continuous, shared_axis) by comparing the first two spectra."""
        if not selected:
            return False, None
        mz0, _ = parser.getspectrum(selected[0][0])
        mz0 = np.asarray(mz0, dtype=float)
        if len(selected) == 1:
            return True, mz0
        mz1, _ = parser.getspectrum(selected[1][0])
        mz1 = np.asarray(mz1, dtype=float)
        if len(mz1) == len(mz0) and np.allclose(mz1, mz0):
            return True, mz0
        return False, None

    @staticmethod
    def _channel_indices(axis, annotations, ppm):
        """Axis index per annotation m/z (nearest within ppm), else -1."""
        idxs = []
        for ann in annotations:
            j = int(np.searchsorted(axis, ann.mz))
            best, best_diff = -1, None
            for k in (j - 1, j):
                if 0 <= k < len(axis):
                    d = abs(float(axis[k]) - ann.mz)
                    if best_diff is None or d < best_diff:
                        best_diff, best = d, k
            if best >= 0 and best_diff is not None and (best_diff / ann.mz) * 1e6 <= ppm:
                idxs.append(best)
            else:
                idxs.append(-1)
        return idxs

    def _auto_annotate(
        self, parser, selected, database, mz_tolerance_ppm, adducts, top_n_peaks, progress
    ) -> List[IonAnnotation]:
        """Build IonAnnotations by matching observed peaks to a molecule database."""
        if database is None:
            raise ValueError(
                "No annotations supplied and no molecule database given for "
                "auto-annotation. Pass annotations=... or database=<path or DB>."
            )
        from ..annotation.mz_annotator import MetaboliteDatabase, default_adducts

        if isinstance(database, (str, Path)):
            self._report(progress, "db-load", 0, 1)
            database = MetaboliteDatabase.load(database)
            self._report(progress, "db-load", 1, 1)
        if adducts is None:
            adducts = default_adducts(getattr(parser, "polarity", None))
        database.build_index(adducts)

        continuous, axis = self._axis_if_continuous(parser, selected)
        if continuous and axis is not None:
            # Match the shared axis directly -- no full discovery pass needed.
            observed = [float(v) for v in axis]
            self._report(progress, "discover", 1, 1)
        else:
            observed = self._discover_processed(parser, selected, top_n_peaks, progress)

        matches = database.annotate(observed, ppm_tol=mz_tolerance_ppm, adducts=adducts)
        if top_n_peaks and len(matches) > top_n_peaks:
            matches = self._rank_matches(
                parser, selected, axis if continuous else None, matches,
                mz_tolerance_ppm, top_n_peaks,
            )
        annotations = [
            IonAnnotation(
                name=m.label, mz=m.observed_mz, adduct=m.adduct, formula=m.formula,
                theoretical_mz=m.theoretical_mz, ppm_error=m.ppm_error, candidates=m.candidates,
            )
            for m in matches
        ]
        logger.info(
            "Auto-annotated %d ions from %d observed peaks (%d formulas in DB)",
            len(annotations), len(observed), len(database.by_formula),
        )
        return annotations

    # Cap how many pixels are read to build a representative peak list / ranking.
    _MAX_DISCOVERY_PIXELS = 2000

    @classmethod
    def _sample(cls, selected):
        """Evenly sample up to _MAX_DISCOVERY_PIXELS selected pixels."""
        n = len(selected)
        if n <= cls._MAX_DISCOVERY_PIXELS:
            return selected
        stepf = n / cls._MAX_DISCOVERY_PIXELS
        return [selected[int(i * stepf)] for i in range(cls._MAX_DISCOVERY_PIXELS)]

    @classmethod
    def _discover_processed(cls, parser, selected, top_n, progress):
        """Pool peaks from a pixel sample into ~1 mDa bins, ranked by intensity."""
        sample = cls._sample(selected)
        acc: Dict[float, float] = {}
        total = len(sample)
        step = max(1, total // 25)
        for n, (idx, _x, _y, _z) in enumerate(sample):
            mzs, inten = parser.getspectrum(idx)
            for mz, it in zip(np.asarray(mzs, dtype=float), np.asarray(inten, dtype=float)):
                key = round(float(mz), 3)
                acc[key] = acc.get(key, 0.0) + float(it)
            if n % step == 0 or n + 1 == total:
                cls._report(progress, "discover", n + 1, total)
        items = sorted(acc.items(), key=lambda kv: kv[1], reverse=True)
        if top_n:
            items = items[: max(top_n * 4, top_n)]  # over-select; ranked again after matching
        return sorted(k for k, _ in items)

    @classmethod
    def _rank_matches(cls, parser, selected, axis, matches, ppm, top_n):
        """Keep the top_n matches by summed intensity over a pixel sample."""
        sample = cls._sample(selected)
        # Precompute a channel index per match for the continuous fast path.
        channel = None
        if axis is not None:
            channel = []
            for m in matches:
                j = int(np.searchsorted(axis, m.observed_mz))
                best, best_diff = -1, None
                for k in (j - 1, j):
                    if 0 <= k < len(axis):
                        d = abs(float(axis[k]) - m.observed_mz)
                        if best_diff is None or d < best_diff:
                            best_diff, best = d, k
                channel.append(best)
        totals = np.zeros(len(matches), dtype=float)
        for (idx, _x, _y, _z) in sample:
            mzs, inten = parser.getspectrum(idx)
            inten = np.asarray(inten, dtype=float)
            if channel is not None:
                for j, ci in enumerate(channel):
                    if 0 <= ci < len(inten):
                        totals[j] += inten[ci]
            else:
                mzs = np.asarray(mzs, dtype=float)
                for j, m in enumerate(matches):
                    totals[j] += _match_intensity(mzs, inten, m.observed_mz, ppm)
        keep = np.argsort(totals)[::-1][:top_n]
        return [matches[i] for i in sorted(keep)]
