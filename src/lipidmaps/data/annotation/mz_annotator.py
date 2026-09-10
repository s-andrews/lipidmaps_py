"""Annotate observed m/z peaks against a molecule database (e.g. CoreMetabolome v3).

imzML gives only m/z; this module turns a molecule list (name + formula) into
theoretical adduct m/z and matches observed peaks within a ppm tolerance. Matches
are made at the **(formula, adduct)** level -- m/z cannot distinguish isomers -- and
every molecule name sharing that formula is retained as a *candidate* so downstream
standardization can resolve an LM ID from the names.

This is a simple monoisotopic mass match (no isotope-pattern / FDR scoring like
METASPACE); treat results as putative. Pure/offline: no network, numpy only.
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Monoisotopic masses (u) of the elements present in CoreMetabolome v3, plus the
# electron mass for charge-carrier corrections.
_ELECTRON_MASS = 0.00054858
ELEMENT_MASSES: Dict[str, float] = {
    "H": 1.0078250319,
    "C": 12.0,
    "N": 14.0030740052,
    "O": 15.9949146221,
    "P": 30.97376151,
    "S": 31.97207069,
    "Na": 22.98976928,
    "K": 38.96370649,
    "Cl": 34.96885268,
    "Br": 78.9183376,
    "F": 18.99840322,
    "I": 126.904473,
    "Co": 58.9331943,
    "Fe": 55.9349375,
    "Se": 79.9165218,
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def formula_monoisotopic_mass(formula: str) -> float:
    """Monoisotopic neutral mass of a bracket-free molecular formula (e.g. ``C6H12O6``).

    Raises ``KeyError`` for an unknown element and ``ValueError`` for a formula that
    does not parse cleanly (leftover characters).
    """
    if not formula:
        raise ValueError("empty formula")
    mass = 0.0
    consumed = 0
    for match in _FORMULA_TOKEN.finditer(formula):
        consumed += len(match.group(0))
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        mass += ELEMENT_MASSES[element] * count
    if consumed != len(formula):
        raise ValueError(f"could not fully parse formula {formula!r}")
    return mass


@dataclass(frozen=True)
class Adduct:
    """A singly/multiply charged adduct: m/z = (M + delta_mass) / charge."""

    name: str
    delta_mass: float
    charge: int = 1
    polarity: str = "positive"

    def mz(self, neutral_mass: float) -> float:
        return (neutral_mass + self.delta_mass) / abs(self.charge)


# Default adducts by polarity. Deltas include the charge carrier and electron mass.
_PROTON = 1.0072764666  # H - electron
POSITIVE_ADDUCTS: List[Adduct] = [
    Adduct("[M+H]+", _PROTON, 1, "positive"),
    Adduct("[M+Na]+", ELEMENT_MASSES["Na"] - _ELECTRON_MASS, 1, "positive"),
    Adduct("[M+K]+", ELEMENT_MASSES["K"] - _ELECTRON_MASS, 1, "positive"),
]
NEGATIVE_ADDUCTS: List[Adduct] = [
    Adduct("[M-H]-", -_PROTON, 1, "negative"),
    Adduct("[M+Cl]-", ELEMENT_MASSES["Cl"] + _ELECTRON_MASS, 1, "negative"),
]


def default_adducts(polarity: Optional[str]) -> List[Adduct]:
    """Default adduct set for a polarity string; positive when unknown/None."""
    if polarity and str(polarity).lower().startswith("neg"):
        return list(NEGATIVE_ADDUCTS)
    return list(POSITIVE_ADDUCTS)


@dataclass
class FormulaMatch:
    """A matched (formula, adduct) with its observed peak and candidate molecules."""

    formula: str
    adduct: str
    observed_mz: float
    theoretical_mz: float
    ppm_error: float
    candidates: List[Dict[str, str]] = field(default_factory=list)

    @property
    def label(self) -> str:
        return f"{self.formula} {self.adduct}"


class MetaboliteDatabase:
    """A molecule database grouped by formula, with theoretical adduct m/z indexing.

    Load from a CoreMetabolome-style table (columns ``id``, ``name``, ``formula``;
    tab- or comma-delimited). Isomers collapse under one formula; their names/ids are
    kept as candidates.
    """

    def __init__(self, by_formula: Dict[str, Dict]):
        # by_formula: {formula: {"mass": float, "candidates": [{"name","id"}, ...]}}
        self.by_formula = by_formula
        self._index_mz: Optional[np.ndarray] = None
        self._index_keys: List[Tuple[str, str]] = []  # (formula, adduct) per index row
        self._index_adducts: List[Adduct] = []

    @classmethod
    def load(cls, path: str | Path) -> "MetaboliteDatabase":
        path = Path(path)
        # Sniff delimiter: CoreMetabolome ships tab-delimited despite a .csv name.
        with path.open("r", encoding="utf-8", newline="") as handle:
            first = handle.readline()
        delimiter = "\t" if first.count("\t") >= first.count(",") else ","

        by_formula: Dict[str, Dict] = {}
        skipped = 0
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            fmap = {(k or "").strip().lower(): k for k in (reader.fieldnames or [])}
            name_col = fmap.get("name")
            formula_col = fmap.get("formula")
            id_col = fmap.get("id")
            if not name_col or not formula_col:
                raise ValueError(
                    f"Database {path} needs 'name' and 'formula' columns; got {reader.fieldnames}"
                )
            for row in reader:
                formula = (row.get(formula_col) or "").strip()
                name = (row.get(name_col) or "").strip()
                if not formula or not name:
                    continue
                entry = by_formula.get(formula)
                if entry is None:
                    try:
                        mass = formula_monoisotopic_mass(formula)
                    except (KeyError, ValueError):
                        skipped += 1
                        continue
                    entry = {"mass": mass, "candidates": []}
                    by_formula[formula] = entry
                entry["candidates"].append(
                    {"name": name, "id": (row.get(id_col) or "").strip() if id_col else ""}
                )
        logger.info(
            "Loaded metabolite DB %s: %d formulas (%d rows skipped for unparseable formula)",
            path, len(by_formula), skipped,
        )
        return cls(by_formula)

    def build_index(self, adducts: List[Adduct]) -> None:
        """Precompute a sorted array of theoretical adduct m/z for fast matching."""
        mzs: List[float] = []
        keys: List[Tuple[str, str]] = []
        for formula, entry in self.by_formula.items():
            for adduct in adducts:
                mzs.append(adduct.mz(entry["mass"]))
                keys.append((formula, adduct.name))
        order = np.argsort(mzs)
        self._index_mz = np.asarray(mzs, dtype=float)[order]
        self._index_keys = [keys[i] for i in order]
        self._index_adducts = list(adducts)

    def annotate(
        self, observed_mzs, ppm_tol: float = 5.0, adducts: Optional[List[Adduct]] = None
    ) -> List[FormulaMatch]:
        """Match observed peaks to theoretical adduct m/z within ``ppm_tol``.

        Returns one :class:`FormulaMatch` per matched ``(formula, adduct)``, keeping the
        closest observed peak, with all molecules of that formula as candidates.
        """
        if self._index_mz is None:
            self.build_index(adducts or POSITIVE_ADDUCTS)
        index_mz = self._index_mz
        assert index_mz is not None

        best: Dict[Tuple[str, str], FormulaMatch] = {}
        for obs in observed_mzs:
            obs = float(obs)
            tol = obs * ppm_tol * 1e-6
            lo = int(np.searchsorted(index_mz, obs - tol, side="left"))
            hi = int(np.searchsorted(index_mz, obs + tol, side="right"))
            for j in range(lo, hi):
                theo = float(index_mz[j])
                ppm = abs(theo - obs) / obs * 1e6
                if ppm > ppm_tol:
                    continue
                formula, adduct_name = self._index_keys[j]
                key = (formula, adduct_name)
                prev = best.get(key)
                if prev is None or ppm < prev.ppm_error:
                    best[key] = FormulaMatch(
                        formula=formula,
                        adduct=adduct_name,
                        observed_mz=obs,
                        theoretical_mz=theo,
                        ppm_error=ppm,
                        candidates=list(self.by_formula[formula]["candidates"]),
                    )
        return sorted(best.values(), key=lambda m: m.observed_mz)
