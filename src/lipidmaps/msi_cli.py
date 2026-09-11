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
import json
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
    parser.add_argument("imzml", nargs="*", help="Path(s) to .imzML file(s), each with its .ibd alongside. Multiple files are stacked as z-slices (one section per file) into a 3D dataset.")
    parser.add_argument("--stack", action="append", metavar="LABEL=PATH[,PATH...]",
                        help="A labelled stack for BETWEEN-STACK reaction comparison, e.g. "
                             "--stack control=a.imzML --stack disease=b.imzML,c.imzML. Repeat for each stack; "
                             "give two or more to compare reactions across stacks.")
    parser.add_argument("--z-spacing", type=float, help="Section thickness in µm between stacked slices (3D metadata)")
    parser.add_argument("--annotations", help="Optional annotation CSV (name, mz[, adduct]); takes precedence over --database")
    parser.add_argument("--database", help="Molecule DB (CoreMetabolome-style id,name,formula) for auto-annotation")
    parser.add_argument("--database-url", help="Download the molecule DB from this URL and cache it, if no DB is found")
    parser.add_argument("--ppm", type=float, default=5.0, help="m/z match tolerance in ppm (default: 5)")
    parser.add_argument("--threshold", type=float, default=0.05, help="Significance threshold for stack comparison (default: 0.05)")
    parser.add_argument("--replicate-by", choices=["pixel", "section", "tile"], default="pixel",
                        help="Replicate unit for the between-stack t-test: pixel (default; inflates "
                             "significance), section (per z-slice), or tile (spatial grid).")
    parser.add_argument("--tiles", type=int, default=16, help="Number of tiles per stack for --replicate-by tile (default: 16)")
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
    parser.add_argument("--top-reactions", type=int, default=6, help="How many most-prominent/most-changed reactions to render as 3D HTML (default: 6)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser


def _unique_out_dir(base: Path) -> Path:
    """Return ``base`` if free, else ``base1``, ``base2``, … so runs stay distinct."""
    if not base.exists():
        return base
    i = 1
    while True:
        cand = base.with_name(f"{base.name}{i}")
        if not cand.exists():
            return cand
        i += 1


def _parse_stacks(entries):
    """Parse repeated ``--stack LABEL=path[,path...]`` into an ordered {label: [paths]}."""
    stacks = {}
    for entry in entries or []:
        if "=" not in entry:
            raise ValueError(f"--stack must be LABEL=path[,path...]; got {entry!r}")
        label, rest = entry.split("=", 1)
        paths = [p.strip() for p in rest.split(",") if p.strip()]
        if not label.strip() or not paths:
            raise ValueError(f"--stack must be LABEL=path[,path...]; got {entry!r}")
        stacks[label.strip()] = paths
    return stacks


def _run_comparison(args) -> int:
    """Load ≥2 labelled stacks and compare reactions between them (control vs condition)."""
    try:
        stacks = _parse_stacks(args.stack)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if len(stacks) < 2:
        print("error: --stack given fewer than 2 stacks; need at least two to compare.", file=sys.stderr)
        return 2

    all_paths = [Path(p).expanduser() for paths in stacks.values() for p in paths]
    for p in all_paths:
        if not p.exists():
            print(f"error: imzML not found: {p}", file=sys.stderr)
            return 2
        if not p.with_suffix(".ibd").exists():
            print(f"warning: no .ibd beside {p.name}; pyimzml may fail to read spectra.", file=sys.stderr)

    base_out = Path(args.out).expanduser() if args.out else all_paths[0].with_name(all_paths[0].stem + "_compare_out")
    out_dir = _unique_out_dir(base_out)
    out_dir.mkdir(parents=True, exist_ok=True)

    annotation_path = args.annotations
    if not annotation_path:
        maf = _find_maf(all_paths[0].parent)
        if maf is not None:
            annotation_path = str(maf)
            print(f"  using MAF annotation: {maf.name}", file=sys.stderr)
    database = None
    if not annotation_path:
        from .data.annotation.db_provision import resolve_metabolome_db
        try:
            database = str(resolve_metabolome_db(args.database, download_url=args.database_url))
        except FileNotFoundError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

    labels = list(stacks.keys())
    control, condition = labels[0], labels[1]
    print(f"Comparing reactions: control={control} vs condition={condition} "
          f"({len(stacks)} stacks) -> {out_dir}/", file=sys.stderr)

    from . import import_imzml_groups
    from .data.utils.spatial import (
        aggregate_replicates, compare_stack_reactions, reaction_effect_sizes,
    )

    data = import_imzml_groups(
        stacks, annotation_path=annotation_path, database=database,
        mz_tolerance_ppm=args.ppm, top_n_peaks=args.top_n_peaks, stride=args.stride,
        progress=_make_progress(), fetch_reactions=True,
    )
    dataset = data.dataset
    print(f"  merged: {len(dataset.lipids)} ions, {len(dataset.spatial_samples())} pixels, "
          f"groups={labels}", file=sys.stderr)

    # Replicate unit for the t-test (pixel by default; section/tile avoid pixel inflation).
    compare_ds = dataset
    if args.replicate_by != "pixel":
        compare_ds = aggregate_replicates(dataset, args.replicate_by, args.tiles)
        counts = {}
        for s in compare_ds.samples:
            counts[s.group] = counts.get(s.group, 0) + 1
        print(f"  replicate-by {args.replicate_by}: {counts}", file=sys.stderr)
        for g, n in counts.items():
            if n < 2:
                print(f"warning: group '{g}' has {n} replicate(s); z-score needs >=2 "
                      f"(try --replicate-by tile or more sections).", file=sys.stderr)

    comparison = compare_stack_reactions(compare_ds, control, condition, threshold=args.threshold)
    rows = comparison["rows"]
    with (out_dir / "reaction_comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction", "source", "target", "z_score", "direction", "significant"])
        for r in rows:
            w.writerow([r["reaction"], r["source"], r["target"], f"{r['z_score']:.4g}",
                        r["direction"], r["significant"]])
    (out_dir / "reaction_comparison_graph.json").write_text(
        json.dumps(comparison["graph"], indent=2), encoding="utf-8")

    effects = reaction_effect_sizes(compare_ds, control, condition)
    with (out_dir / "reaction_effect_sizes.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["reaction", "reactant", "product",
                    f"mean_log2ratio_{control}", f"mean_log2ratio_{condition}", "effect_size_delta"])
        for e in effects:
            w.writerow([e["reaction"], e["reactant"], e["product"],
                        e.get(f"mean_log2ratio_{control}"), e.get(f"mean_log2ratio_{condition}"),
                        e["effect_size_delta"]])

    # 3D HTML for the most-changed reactions, with genes. Per reaction we write:
    #  - <reaction>.html         : per-stack cubes, coloured relative to each section
    #  - <reaction>_tilediff.html: one cube of aligned-tile Δ between the stacks
    n_r3d = 0
    if args.top_reactions and effects:
        import math
        from .data.utils.spatial import reaction_ratio_cloud, tile_aligned_delta
        from .data.utils import spatial_html as html
        r3d_dir = out_dir / "reactions3d"
        r3d_dir.mkdir(parents=True, exist_ok=True)
        zsp = getattr(dataset, "z_spacing_um", None)
        n_tiles = max(2, int(round(math.sqrt(max(4, args.tiles)))))
        for e in effects[:args.top_reactions]:
            if e["effect_size_delta"] is None:
                continue
            base = _safe(e["reaction"])
            ca, va = reaction_ratio_cloud(dataset, e["reactant"], e["product"], group=control)
            cb, vb = reaction_ratio_cloud(dataset, e["reactant"], e["product"], group=condition)
            title = f"{e['reaction']}  (Δlog2 ratio {e['effect_size_delta']:+.2f})"
            fig = html.reaction_volume_figure(
                [(control, ca, va), (condition, cb, vb)], title,
                gene_text=e.get("genes"), z_spacing_um=zsp, color_mode="section",
            )
            if fig is not None:
                fig.write_html(str(r3d_dir / f"{base}.html"),
                               include_plotlyjs="directory", full_html=True)
                n_r3d += 1
            delta = tile_aligned_delta(dataset, e["reactant"], e["product"], control, condition, n=n_tiles)
            tfig = html.tile_delta_figure(delta, f"{title} — aligned-tile Δ ({condition}−{control})",
                                          gene_text=e.get("genes"), z_spacing_um=zsp)
            if tfig is not None:
                tfig.write_html(str(r3d_dir / f"{base}_tilediff.html"),
                                include_plotlyjs="directory", full_html=True)
                n_r3d += 1

    sig = sum(1 for r in rows if r["significant"])
    print(f"Done: {len(rows)} reaction edges ({sig} significant at p<{args.threshold}); "
          f"{len(effects)} effect sizes; {n_r3d} reaction 3D files. Output in {out_dir}/",
          file=sys.stderr)
    return 0


def _find_maf(directory: Path) -> Optional[Path]:
    """Find a MetaboLights MAF (``m_MTBLS*..._maf.tsv``) in ``directory``, if any."""
    for pattern in ("m_*_maf.tsv", "*_maf.tsv", "*maf.tsv"):
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


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


def _ion_filebase(lipid) -> str:
    """Filename stem for an ion: molecule name + formula label (keeps uniqueness)."""
    from .data.utils.spatial import ion_display_name

    name = ion_display_name(lipid)
    if name and name != lipid.input_name:
        return f"{_safe(name)}__{_safe(lipid.input_name)}"
    return _safe(lipid.input_name)


def _render_images(dataset, out_dir: Path, max_images: int, render_html: bool = False) -> int:
    from .data.utils.spatial import (
        ion_display_full, lipid_spatial_series, ratio_series, spatial_reactions, z_slices,
    )
    from .data.utils import spatial_render as sr
    html = None
    if render_html:
        from .data.utils import spatial_html as html

    ps = getattr(dataset, "pixel_size_um", None)
    zsp = getattr(dataset, "z_spacing_um", None)
    is_3d = dataset.is_3d
    images = out_dir / "images"
    ranked = sorted(dataset.lipids, key=_ion_total, reverse=True)
    written = 0
    for lp in ranked[:max_images]:
        coords, values = lipid_spatial_series(dataset, lp.input_name)
        if not coords:
            continue
        z = (z_slices(coords) or [0])[0]  # first acquired z-slice for the 2D maps
        base = _ion_filebase(lp)
        title = ion_display_full(lp)  # molecule name (formula [adduct]) when known
        sr.save_ion_grid_png(coords, values, images / f"ion_{base}_grid.png",
                             title, z=z, pixel_size_um=ps)
        sr.save_voronoi_png(coords, values, images / f"ion_{base}_voronoi.png",
                            f"{title} (Voronoi)", z=z, pixel_size_um=ps)
        if html is not None:
            _write_html(html.grid_figure(coords, values, z, title, pixel_size_um=ps),
                        images / f"ion_{base}_grid.html")
            _write_html(html.voronoi_figure(coords, values, z, f"{title} (Voronoi)",
                                            pixel_size_um=ps), images / f"ion_{base}_voronoi.html")
        if is_3d:  # volumetric point cloud across all z-slices
            sr.save_scatter3d_png(coords, values, images / f"ion_{base}_3d.png",
                                  f"{title} (3D)", z_spacing_um=zsp)
            if html is not None:
                _write_html(html.scatter3d_figure(coords, values, f"{title} (3D)", z_spacing_um=zsp),
                            images / f"ion_{base}_3d.html")
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

    # For a true 3D volume, also render the most-prominent reactions as 3D HTML (+genes).
    if is_3d and html is not None:
        from .data.utils.spatial import (
            prominent_reactions, reaction_gene_labels, reaction_ratio_cloud,
        )
        zsp = getattr(dataset, "z_spacing_um", None)
        for rx, react_ions, prod_ions, _score in prominent_reactions(dataset, top=max_images):
            coords3d, ratio3d = reaction_ratio_cloud(dataset, react_ions[0], prod_ions[0])
            rtitle = rx.reaction_name or "reaction"
            fig = html.reaction_volume_figure(
                [("log2(product/reactant)", coords3d, ratio3d)], rtitle,
                gene_text=", ".join(reaction_gene_labels(rx)), z_spacing_um=zsp,
                color_mode="section",
            )
            if fig is not None:
                _write_html(fig, images / f"rxn3d_{_safe(rtitle)}.html")
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

    # Between-stack reaction comparison mode.
    if args.stack:
        return _run_comparison(args)
    if not args.imzml:
        print("error: give an imzML file (or use --stack LABEL=... twice to compare stacks).",
              file=sys.stderr)
        return 2

    paths = [Path(p).expanduser() for p in args.imzml]
    for p in paths:
        if not p.exists():
            print(f"error: imzML not found: {p}", file=sys.stderr)
            return 2
        if not p.with_suffix(".ibd").exists():
            print(f"warning: no .ibd beside {p.name}; pyimzml may fail to read spectra.",
                  file=sys.stderr)
    stacked = len(paths) > 1
    if args.out:
        base_out = Path(args.out).expanduser()
    else:
        suffix = "_stack_out" if stacked else "_out"
        base_out = paths[0].with_name(paths[0].stem + suffix)
    out_dir = _unique_out_dir(base_out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Annotation source precedence: explicit --annotations, else a sibling MetaboLights
    # MAF next to the imzML, else auto-annotate from the provisioned molecule DB.
    annotation_path = args.annotations
    if not annotation_path:
        maf = _find_maf(paths[0].parent)
        if maf is not None:
            annotation_path = str(maf)
            print(f"  using MAF annotation: {maf.name}", file=sys.stderr)

    database = None
    if not annotation_path:
        from .data.annotation.db_provision import resolve_metabolome_db
        try:
            database = str(resolve_metabolome_db(args.database, download_url=args.database_url))
            print(f"  annotation DB: {database}", file=sys.stderr)
        except FileNotFoundError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

    common = dict(
        annotation_path=annotation_path, database=database, mz_tolerance_ppm=args.ppm,
        bbox=tuple(args.bbox) if args.bbox else None, top_n_peaks=args.top_n_peaks,
        stride=args.stride, progress=_make_progress(),
        use_headgroups=args.use_headgroups, fetch_reactions=args.fetch_reactions,
    )
    if args.stride > 1:
        print(f"  subsampling: every {args.stride}th pixel", file=sys.stderr)

    if stacked:
        from . import import_imzml_stack
        print(f"Stacking {len(paths)} imzML files as z-slices -> {out_dir}/", file=sys.stderr)
        data = import_imzml_stack(paths, z_spacing_um=args.z_spacing, **common)
    else:
        from . import import_imzml
        print(f"Processing {paths[0].name} -> {out_dir}/", file=sys.stderr)
        data = import_imzml(paths[0], **common)
    dataset = data.dataset

    dims = (f"3D volume: {dataset.z_slice_count} z-slices" if dataset.is_3d
            else "2D (single z-slice)")
    print(f"  {dims}", file=sys.stderr)

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
        extra = " + 3D" if dataset.is_3d else ""
        print(f"  rendering up to {args.max_images} ion/reaction diagrams ({kind}{extra}) ...",
              file=sys.stderr)
        n_img = _render_images(dataset, out_dir, args.max_images, render_html=args.render_html)

    resolved = sum(1 for lp in dataset.lipids if lp.lm_id or lp.generic_lm_id)
    print(
        f"Done [{dims}]: {n_ions} ions ({resolved} with LM IDs) across "
        f"{len(dataset.spatial_samples())} pixels; {n_rx} spatial reactions; "
        f"{n_img} images. Output in {out_dir}/",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
