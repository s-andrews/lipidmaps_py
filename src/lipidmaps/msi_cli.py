"""``lipidmaps-msi`` — process a mass-spectrometry-imaging ``.imzML`` from disk.

Reads the imzML/ibd directly (no upload, unlike the Streamlit app, so multi-GB files
are fine), auto-annotates observed peaks against a molecule database (or a supplied
annotation CSV), resolves LM IDs from candidate names, optionally fetches reactions,
and writes data files plus rendered PNG diagrams.

Usage:
    lipidmaps-msi path/to/run.imzML --database core_metabolome_v3.csv --out run_out
    lipidmaps-msi run.imzML --annotations ions.csv --no-fetch-reactions
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lipidmaps-msi",
        description="Process an imzML MSI dataset into a LipidDataset, summaries, and images.",
    )
    parser.add_argument("imzml", help="Path to the .imzML file (its .ibd must sit alongside)")
    parser.add_argument("--annotations", help="Optional annotation CSV (name, mz[, adduct]); takes precedence over --database")
    parser.add_argument("--database", help="Molecule DB (CoreMetabolome-style id,name,formula) for auto-annotation")
    parser.add_argument("--database-url", help="Download the molecule DB from this URL and cache it, if no DB is found")
    parser.add_argument("--ppm", type=float, default=5.0, help="m/z match tolerance in ppm (default: 5)")
    parser.add_argument("--bbox", nargs=4, type=int, metavar=("X0", "Y0", "X1", "Y1"), help="Inclusive pixel crop")
    parser.add_argument("--stride", type=int, default=1, help="Keep every Nth pixel (subsample large datasets; default: 1)")
    parser.add_argument("--top-n-peaks", type=int, help="Cap ions (ranked by intensity) for auto-annotation")
    parser.add_argument("--out", help="Output directory (default: <imzml_stem>_out)")
    parser.add_argument("--no-json", dest="write_json", action="store_false", default=True, help="Skip writing processed_dataset.json (large for big datasets)")
    parser.add_argument("--use-headgroups", action="store_true", default=True, help="Fill generic LM IDs from headgroups")
    parser.add_argument("--no-headgroups", dest="use_headgroups", action="store_false", help="Disable headgroup fill")
    parser.add_argument("--fetch-reactions", action="store_true", default=True, help="Fetch reactions from LIPID MAPS")
    parser.add_argument("--no-fetch-reactions", dest="fetch_reactions", action="store_false", help="Disable reaction fetching")
    parser.add_argument("--no-images", dest="render_images", action="store_false", default=True, help="Skip rendering (data files only)")
    parser.add_argument("--html", dest="render_html", action="store_true", default=False, help="Also write interactive HTML maps (hover for pixel x,y,z + value)")
    parser.add_argument("--max-images", type=int, default=20, help="Max ion maps to render, ranked by total intensity (default: 20)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def _ion_total(lipid) -> float:
    return float(sum(v for v in lipid.values.values() if v))


def _write_ions_csv(dataset, path: Path) -> int:
    rows = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        w.writerow([
            "input_name", "standardized_name", "lm_id", "generic_lm_id",
            "lm_id_found_by", "n_candidates", "n_reactions", "total_intensity",
        ])
        for lp in dataset.lipids:
            w.writerow([
                lp.input_name, lp.standardized_name or "", lp.lm_id or "",
                lp.generic_lm_id or "", lp.lm_id_found_by or "",
                len(lp.annotation_candidates or []),
                len(lp.reactions or []), f"{_ion_total(lp):.4g}",
            ])
            rows += 1
    return rows


def _write_reactions_csv(dataset, path: Path) -> int:
    # Only reactions with measured reactant+product ions (meaningful over tissue).
    from .data.utils.spatial import spatial_reactions

    spatial_rx = spatial_reactions(dataset)
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle)
        w.writerow(["reaction_name", "reactant_ions", "product_ions"])
        for rx, react_ions, prod_ions in spatial_rx:
            w.writerow([
                rx.reaction_name or f"reaction {rx.reaction_id}",
                ";".join(react_ions), ";".join(prod_ions),
            ])
    return len(spatial_rx)


def _safe(name: str) -> str:
    for a, b in (("/", "_"), (" ", "_"), ("[", ""), ("]", ""), (">", "to")):
        name = name.replace(a, b)
    return name


def _write_html(fig, path: Path) -> None:
    # include_plotlyjs="directory" writes plotly.min.js once into the folder; each
    # HTML references it, so files stay small and work offline.
    if fig is not None:
        fig.write_html(str(path), include_plotlyjs="directory", full_html=True)


def _render_images(dataset, out_dir: Path, max_images: int, render_html: bool = False) -> int:
    from .data.utils.spatial import (
        lipid_spatial_series, ratio_series, spatial_reactions, z_slices,
    )
    from .data.utils import spatial_render as sr
    html = None
    if render_html:
        from .data.utils import spatial_html as html

    ps = getattr(dataset, "pixel_size_um", None)
    images = out_dir / "images"
    ranked = sorted(dataset.lipids, key=_ion_total, reverse=True)
    written = 0
    for lp in ranked[:max_images]:
        coords, values = lipid_spatial_series(dataset, lp.input_name)
        if not coords:
            continue
        z = (z_slices(coords) or [0])[0]  # render the first acquired z-slice
        safe = _safe(lp.input_name)
        sr.save_ion_grid_png(coords, values, images / f"ion_{safe}_grid.png",
                             lp.input_name, z=z, pixel_size_um=ps)
        sr.save_voronoi_png(coords, values, images / f"ion_{safe}_voronoi.png",
                            f"{lp.input_name} (Voronoi)", z=z, pixel_size_um=ps)
        if html is not None:
            _write_html(html.grid_figure(coords, values, z, lp.input_name, pixel_size_um=ps),
                        images / f"ion_{safe}_grid.html")
            _write_html(html.voronoi_figure(coords, values, z, f"{lp.input_name} (Voronoi)",
                                            pixel_size_um=ps), images / f"ion_{safe}_voronoi.html")
        written += 1

    # Reaction ratio maps for reactions with measured reactant+product.
    for rx, react_ions, prod_ions in spatial_reactions(dataset)[:max_images]:
        coords, ratio = ratio_series(dataset, react_ions[0], prod_ions[0])
        z = (z_slices(coords) or [0])[0]
        label = _safe(rx.reaction_name or "reaction")
        title = f"{rx.reaction_name}: log2(product/reactant)"
        sr.save_ratio_png(coords, ratio, images / f"rxn_{label}_ratio.png", title, z=z)
        if html is not None:
            _write_html(
                html.grid_figure(coords, ratio, z, title, pixel_size_um=ps,
                                 cmap="RdBu_r", midpoint=0.0),
                images / f"rxn_{label}_ratio.html",
            )
        written += 1
    return written


_PHASE_LABEL = {
    "parse": "Parsing imzML metadata",
    "db-load": "Loading molecule database",
    "discover": "Scanning peaks",
    "extract": "Extracting ions",
}


def _make_progress():
    """A stderr progress printer for the reader's (phase, done, total) callback."""
    state = {"announced_extract": False}

    def progress(phase: str, done: int, total: int) -> None:
        label = _PHASE_LABEL.get(phase, phase)
        if phase == "extract" and not state["announced_extract"]:
            state["announced_extract"] = True
            print(f"  {total} pixels selected for extraction", file=sys.stderr)
        pct = (done / total * 100.0) if total else 100.0
        end = "\n" if done >= total else "\r"
        sys.stderr.write(f"  {label:<26} {done}/{total} ({pct:5.1f}%)      {end}")
        sys.stderr.flush()

    return progress


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    # Root stays quiet (hide pyimzml noise); surface lipidmaps stage logs at INFO so
    # the post-extraction pipeline (RefMet/headgroups/LMSD/reactions) is visible.
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    if not args.verbose:
        logging.getLogger("lipidmaps").setLevel(logging.INFO)

    imzml = Path(args.imzml).expanduser()
    if not imzml.exists():
        print(f"error: imzML not found: {imzml}", file=sys.stderr)
        return 2
    if not imzml.with_suffix(".ibd").exists():
        print(f"warning: no .ibd beside {imzml.name}; pyimzml may fail to read spectra.",
              file=sys.stderr)
    out_dir = Path(args.out).expanduser() if args.out else imzml.with_name(imzml.stem + "_out")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve the annotation source: explicit CSV wins; else provision a molecule DB.
    database = None
    if not args.annotations:
        from .data.annotation.db_provision import resolve_metabolome_db
        try:
            database = str(resolve_metabolome_db(args.database, download_url=args.database_url))
            print(f"  annotation DB: {database}", file=sys.stderr)
        except FileNotFoundError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

    from . import import_imzml

    print(f"Processing {imzml.name} -> {out_dir}/", file=sys.stderr)
    if args.stride > 1:
        print(f"  subsampling: every {args.stride}th pixel", file=sys.stderr)
    data = import_imzml(
        imzml,
        annotation_path=args.annotations,
        database=database,
        mz_tolerance_ppm=args.ppm,
        bbox=tuple(args.bbox) if args.bbox else None,
        top_n_peaks=args.top_n_peaks,
        stride=args.stride,
        progress=_make_progress(),
        use_headgroups=args.use_headgroups,
        fetch_reactions=args.fetch_reactions,
    )
    dataset = data.dataset

    if args.write_json:
        print("  writing processed_dataset.json ...", file=sys.stderr)
        (out_dir / "processed_dataset.json").write_text(
            dataset.model_dump_json(indent=2), encoding="utf-8"
        )
    print("  writing ions.csv / reactions.csv ...", file=sys.stderr)
    n_ions = _write_ions_csv(dataset, out_dir / "ions.csv")
    n_rx = _write_reactions_csv(dataset, out_dir / "reactions.csv")
    n_img = 0
    if args.render_images:
        kind = "PNG + interactive HTML" if args.render_html else "PNG"
        print(f"  rendering up to {args.max_images} ion/reaction diagrams ({kind}) ...", file=sys.stderr)
        n_img = _render_images(dataset, out_dir, args.max_images, render_html=args.render_html)

    resolved = sum(1 for lp in dataset.lipids if lp.lm_id or lp.generic_lm_id)
    print(
        f"Done: {n_ions} ions ({resolved} with LM IDs) across "
        f"{len(dataset.spatial_samples())} pixels; {n_rx} spatial reactions; "
        f"{n_img} images. Output in {out_dir}/",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
