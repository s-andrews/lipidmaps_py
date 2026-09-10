"""Generate the small synthetic MSI demo fixture used by tests and the Streamlit app.

This produces an offline-safe ``.imzML``/``.ibd`` plus an ``annotations.csv`` under
``tests/data/inputs/demo/msi/``. It is **synthetic** (a tissue-like ellipse with
per-lipid spatial patterns) but uses **real lipid names** and plausible ``[M+H]+``
m/z, so the standardization + reaction pipeline behaves realistically. For a real
open-data slice, use ``scripts/fetch_msi_demo.py`` instead (it overwrites these
files); see the README written beside the fixture.

Run:  python scripts/make_synth_msi_demo.py
"""

from __future__ import annotations

from pathlib import Path

import sys

import numpy as np
from pyimzml.ImzMLWriter import ImzMLWriter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidmaps.data.annotation.mz_annotator import formula_monoisotopic_mass  # noqa: E402

OUT_DIR = Path(__file__).resolve().parents[1] / "tests" / "data" / "inputs" / "demo" / "msi"
WIDTH, HEIGHT = 40, 30
_PROTON = 1.0072764666

# Real lipid names + formulas confirmed present in CoreMetabolome v3, so the
# annotation-free demo (matching against that DB) produces near-0-ppm hits. Several
# are PC/LPC/PA (reaction-bearing) for the reaction-overlay demo. Peaks are [M+H]+.
LIPIDS = [
    # name, formula, spatial-pattern kind
    ("LysoPC(16:0)", "C24H50NO7P", "left_gradient"),
    ("PA(20:0/14:1)", "C37H71O8P", "center_blob"),
    ("SM(d18:1/16:0)", "C39H79N2O6P", "top_blob"),
    ("PC(31:1)", "C39H76NO8P", "ring"),
    ("PC(34:1)", "C42H82NO8P", "right_gradient"),
    ("PC(36:1)", "C44H84NO8P", "bottom_blob"),
]
# Attach [M+H]+ m/z computed from the formula.
LIPIDS = [(name, formula_monoisotopic_mass(formula) + _PROTON, kind, formula)
          for (name, formula, kind) in LIPIDS]


def _tissue_mask() -> np.ndarray:
    """An elliptical 'tissue' region so background pixels read as absent."""
    yy, xx = np.mgrid[0:HEIGHT, 0:WIDTH]
    cx, cy = WIDTH / 2.0, HEIGHT / 2.0
    return (((xx - cx) / (WIDTH * 0.46)) ** 2 + ((yy - cy) / (HEIGHT * 0.46)) ** 2) <= 1.0


def _pattern(kind: str) -> np.ndarray:
    yy, xx = np.mgrid[0:HEIGHT, 0:WIDTH]
    cx, cy = WIDTH / 2.0, HEIGHT / 2.0
    if kind == "left_gradient":
        field = 1.0 - xx / WIDTH
    elif kind == "right_gradient":
        field = xx / WIDTH
    elif kind == "top_blob":
        field = np.exp(-(((xx - cx) ** 2 + (yy - cy * 0.5) ** 2) / (2 * 6.0 ** 2)))
    elif kind == "bottom_blob":
        field = np.exp(-(((xx - cx) ** 2 + (yy - cy * 1.5) ** 2) / (2 * 6.0 ** 2)))
    elif kind == "center_blob":
        field = np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 5.0 ** 2)))
    elif kind == "ring":
        r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
        field = np.exp(-((r - 9.0) ** 2) / (2 * 2.5 ** 2))
    else:  # uniform
        field = np.full((HEIGHT, WIDTH), 0.7)
    return field


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    imzml_path = OUT_DIR / "demo_brain.imzML"

    mask = _tissue_mask()
    rng = np.random.default_rng(42)
    # imzML m/z arrays are always sorted ascending; write peaks in that order so
    # the reader's nearest-peak search behaves like it does on real data.
    lipids_sorted = sorted(LIPIDS, key=lambda item: item[1])
    mz_axis = np.array([mz for _, mz, _, _ in lipids_sorted], dtype=float)
    fields = [_pattern(kind) for _, _, kind, _ in lipids_sorted]

    with ImzMLWriter(str(imzml_path)) as writer:
        for y in range(HEIGHT):
            for x in range(WIDTH):
                if not mask[y, x]:
                    continue
                base = np.array([f[y, x] for f in fields], dtype=float)
                noise = rng.normal(1.0, 0.05, size=base.shape)
                intensities = np.clip(base * 1000.0 * noise, 0.0, None)
                # imzML coordinates are 1-based.
                writer.addSpectrum(mz_axis, intensities, (x + 1, y + 1, 1))

    ann_path = OUT_DIR / "annotations.csv"
    with ann_path.open("w", encoding="utf-8") as handle:
        handle.write("name,mz,adduct,formula\n")
        for name, mz, _, formula in LIPIDS:
            handle.write(f"{name},{mz:.4f},[M+H]+,{formula}\n")

    n_pixels = int(mask.sum())
    print(f"Wrote {imzml_path} ({n_pixels} pixels), {ann_path} ({len(LIPIDS)} ions)")


if __name__ == "__main__":
    main()
