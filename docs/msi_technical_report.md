# Technical Report — `msi` branch: Mass-Spectrometry-Imaging (MSI) support

## 1. Goal of this branch

`lipidmaps_py` on `main` ingests **tabular** lipidomics (CSV / MS-DIAL), standardizes
names (RefMet / LMSD), fetches reactions, and visualizes results in a Streamlit demo.
It has **no spatial capability**.

The `msi` branch adds **mass-spectrometry imaging** support end to end:

- Read **`.imzML` / `.ibd`** datasets, where every mass spectrum is tied to an
  `(x, y, z)` tissue **pixel**.
- Turn m/z peaks into **lipid/metabolite names** (imzML carries no names), from either
  a supplied annotation file (**MetaboLights MAF** or a simple `name,mz` CSV) or by
  **auto-annotating** observed peaks against a molecule database (**CoreMetabolome v3**).
- Reuse the existing **standardization + reaction** pipeline so spatial data flows into
  LM IDs, reactions and BioPAN unchanged.
- **Visualize** ion distributions across tissue (grid ion images + **Voronoi** maps),
  overlay **reactions over tissue** (reactant vs product vs log2-ratio), and stack
  multiple sections into **3D**.
- Do all of this from a **CLI** (`lipidmaps-msi`) for multi-GB files that can't go
  through a browser upload, and from the **Streamlit** app for small/demo data.

Design tenet: **reuse, don't fork.** MSI plugs into the existing `LipidDataset` /
`QuantifiedLipid` spine rather than a parallel model, so every downstream feature keeps
working.

## 2. What changed compared to `main`

28 files, ~+41k lines (most of which is the committed CoreMetabolome DB + the demo
imzML fixture). Grouped by area:

| Area | New | Modified |
|---|---|---|
| **Data models** | — | `data/models/sample.py` (additive fields + helpers) |
| **Pipeline** | — | `data/data_manager.py` (shared `_annotate_and_react`, candidate resolution) |
| **Ingestion (binary)** | `data/ingestion/imzml_reader.py` | `data/ingestion/__init__.py` |
| **Annotation** | `data/annotation/mz_annotator.py`, `data/annotation/db_provision.py`, `data/annotation/__init__.py` | — |
| **Importers / API** | — | `data_importer.py` (`import_imzml`, `import_imzml_stack`), `__init__.py` |
| **Spatial utils** | `data/utils/spatial.py`, `data/utils/spatial_html.py` (plotly), `data/utils/spatial_render.py` (matplotlib) | — |
| **CLI** | `msi_cli.py` (`lipidmaps-msi`) | `setup.py` (console script + `[msi]` extra) |
| **Streamlit** | `scripts/spatial_ui.py` | `scripts/streamlit_demo.py` (Spatial tab, imzML/MAF upload, stacking) |
| **Demo / data** | `scripts/make_synth_msi_demo.py`, `scripts/fetch_msi_demo.py`, `tests/data/inputs/demo/msi/*`, `tests/data/core_metabolome_v3.csv` | — |
| **Tests** | `tests/test_spatial_utils.py`, `tests/test_mz_annotator.py`, `tests/test_msi_cli.py`, `tests/data/test_imzml_ingestion.py` | — |

New optional dependencies (in the `[msi]` extra only — core stays lean):
`pyimzml`, `metaspace2020`, `matplotlib`, `plotly`.

New public API: `import_imzml`, `import_imzml_stack` (top-level), and the
`lipidmaps-msi` console script.

## 3. How we updated the base classes

All model changes are **additive and `Optional`** so the tabular workflow and the
downstream BioPAN importers are untouched (`None` for non-spatial data).

### `data/models/sample.py`

- **`PixelCoordinate(LipidmapsBaseModel)`** — `x: int`, `y: int`, `z: int = 0`. imzML
  coordinates are always 3-tuples (`z = 0`/`1` for a single layer), so 3D costs nothing
  to model.
- **`SampleMetadata.coordinates: Optional[PixelCoordinate] = None`** — the *only* field
  needed to make a "sample" spatial.
- **`QuantifiedLipid.annotation_candidates: Optional[List[Dict[str, str]]] = None`** —
  carries the candidate molecule names (+ ids) for a formula-level auto-annotated ion,
  so the DB names can drive standardization.
- **`LipidDataset`** gained:
  - `pixel_size_um: Optional[Tuple[float, float]]` — physical pixel size from imzML.
  - `z_spacing_um: Optional[float]` — section thickness for a stacked volume.
  - helpers: `is_spatial`, `spatial_samples()`, `z_layers()`, `z_slice_count`,
    **`is_3d`** (True only when `z_slice_count > 1` — a single layer at z=0/1 is 2D).

### The "pixel = sample" architecture

The key decision: **each pixel is a `sample`.** A pixel's `sample_name` encodes its
coordinate (`px_x0012_y0034_z01`), and its ion intensity lands in the existing
`QuantifiedLipid.values[sample_name]` dict. Consequences:

- Standardization, LM-ID lookup, reaction fetching, normalization and BioPAN work with
  **zero changes** — they already operate over `values` keyed by sample name.
- A "3D dataset" is just samples whose `coordinates.z` varies.

### `data/data_manager.py`

- Extracted the post-ingestion steps of `process_csv` into a shared
  **`DataManager._annotate_and_react(dataset)`**: RefMet → candidate resolution →
  headgroup generic-ID fill → LMSD name fill → reaction fetch → LMSD detail annotate.
  Both `process_csv` and `import_imzml*` call it, so tabular and spatial data are
  standardized identically.
- Added **`_resolve_candidate_lm_ids(dataset)`** (runs before the headgroup/LMSD fills):
  for auto-annotated ions labelled `formula [adduct]`, it batches the candidate molecule
  names through **RefMet** and adopts the best hit (prefers one yielding a specific
  `lm_id`, else one yielding a `standardized_name`), recording
  `lm_id_found_by = "candidate"`. This is how CoreMetabolome names become LM IDs and
  unlock reactions.

## 4. How we read the new file types

Three new inputs, each parsed into the same neutral shape (a list of `IonAnnotation`
= `{name, mz, adduct?, formula?, lm_id?, candidates?}` and, for the DB, a formula index).

### 4a. imzML / ibd (the binary MSI data) — §5 below

### 4b. MetaboLights **MAF** and simple annotation CSV — `parse_annotation_csv()`

MAF (`m_MTBLS…_maf.tsv`) is the curated per-study assignment file. The parser:

- **sniffs the delimiter** (tab vs comma — MAF is tab-delimited),
- maps columns case-insensitively:
  `metabolite_identification`→**name**, `mass_to_charge`→**m/z**,
  `chemical_formula`→formula, `database_identifier`→id (used as `lm_id` **only** when it
  starts with `LM`; HMDB/ChEBI ids are ignored, not mis-assigned),
- also accepts a plain `name,mz[,adduct,formula,lm_id]` CSV and METASPACE-style headers.

MAF m/z are **observed ion m/z**, so they match imzML peaks directly (no adduct math).
Validated on **MTBLS12782**: 367 annotations parsed, 89 resolved to LM IDs.

### 4c. CoreMetabolome v3 molecule database — `MetaboliteDatabase.load()`

A tab-delimited `id, name, formula, inchi` table (11,440 molecules) used for
auto-annotation when no MAF/CSV is supplied. Loaded once, **grouped by formula**
(isomers collapse into one entry with a list of candidate names), and each unique
formula's **monoisotopic mass** is computed from `formula_monoisotopic_mass()` (a
`([A-Z][a-z]?)(\d*)` tokenizer over a 15-element mass table — all DB formulas are
bracket-free).

### Annotation resolution order (CLI)

`--annotations <file>` → else a **sibling `m_*_maf.tsv`** auto-detected next to the
imzML → else **CoreMetabolome auto-annotation**.

## 5. How we read the binary file (`.imzML` / `.ibd`)

imzML is XML metadata (one `<spectrum>` per pixel) + a binary `.ibd` holding the m/z and
intensity arrays. We use **`pyimzml.ImzMLParser`**, which reads spectra **lazily** by
index (`getspectrum(i)`), so a multi-GB `.ibd` **streams** rather than loading whole.
All of this lives in `data/ingestion/imzml_reader.py` (`ImzMLIngestion.read`).

### Pipeline inside `read()`

1. **Parse metadata** — `ImzMLParser(path)`; extract polarity (for adduct defaults) and
   **pixel size** (`imzmldict["pixel size (x)/(y)"]`, `_pixel_size()`).
2. **Select pixels** — `_selected_indices()` applies an optional `bbox` crop and a
   `stride` (keep every Nth pixel) to bound work/memory on huge files.
3. **Annotate** — supplied annotations win; otherwise `_auto_annotate()` (see §6).
4. **Extract** — one lazy pass over the selected pixels, filling
   `ion_values[ion][pixel]`.
5. Return `ImzMLIngestionResult{pixels, ion_values, annotations, pixel_size}`.

### Continuous vs processed, and the fast path

- **Continuous** files share one m/z axis across all pixels (`_axis_if_continuous()`
  compares the first two spectra). We then:
  - match the DB against the **axis directly** — *no separate discovery pass*, and
  - precompute a **fixed channel index per ion** (`_channel_indices()`), so extraction
    is `intensities[channel]` instead of a per-pixel/per-ion `searchsorted`.
  This is the single biggest speed/robustness win — it removed the multi-minute silent
  hang on the Waters test file.
- **Processed** files have per-pixel m/z arrays; `_discover_processed()` pools peaks from
  a **sample** of pixels (≤ 2000) into ~1 mDa bins, then matches, then extracts with
  `_match_intensity()` (nearest peak within ppm). `_rank_matches()` keeps the top-N by
  intensity over a pixel sample.

### Memory & progress

- Memory cost is the per-ion × per-pixel `values` dicts ("pixel = sample"); `--bbox`,
  `--stride` and `--top-n-peaks` bound it. (See §10 for the full-resolution ceiling.)
- Every phase reports through a `progress(phase, done, total)` callback — the CLI prints
  a live `Parsing / Loading DB / Scanning peaks / Extracting ions  n/total (%)` line, so
  a long read is visibly working.

## 6. Auto-annotation (m/z → molecule) — `data/annotation/mz_annotator.py`

When there's no MAF/CSV, we reproduce a **METASPACE-style** formula+adduct match
(without isotope-pattern / FDR scoring — results are *putative*):

1. For each unique DB formula, compute monoisotopic mass; for each **adduct**
   (`Adduct{name, delta_mass, charge, polarity}`; positive `[M+H]+/[M+Na]+/[M+K]+`,
   negative `[M-H]-/[M+Cl]-`, chosen by imzML polarity), compute theoretical m/z.
2. `build_index()` sorts all `(theoretical_mz, formula, adduct)`; `annotate()` matches
   observed peaks with `searchsorted` within ppm.
3. Emit **one ion per matched `(formula, adduct)`**, labelled `formula [adduct]`, keeping
   **all isomeric candidate names** (m/z cannot distinguish isomers).
4. The pipeline then resolves LM IDs from those candidate names (§3).

Provisioning (`db_provision.resolve_metabolome_db`): `--database` → `$LIPIDMAPS_METABOLOME_DB`
→ user cache (`~/.cache/lipidmaps/`) → in-repo copy → optional download from a supplied
URL. No hardcoded URL (no stable public CoreMetabolome download exists).

## 7. Visualization & outputs

- **`data/utils/spatial.py`** (pure, testable): `values_to_grid` (dense z-slice array),
  `voronoi_regions` (`scipy.spatial.Voronoi` + finite-cell reconstruction + bbox clip,
  each region carries `x,y,z,value,polygon`), `lipid_spatial_series`, `z_slices`,
  `spatial_reactions` / `ratio_series` (reaction-over-tissue), `ion_display_name/full`
  (friendly molecule-name labels), `nice_length` (scale bars).
- **`data/utils/spatial_html.py`** (plotly): `grid_figure`, `voronoi_figure`,
  `scatter3d_figure` — rich hover (**pixel (x,y,z) + sample_name + value**), z in titles,
  a **µm scale bar** when pixel size is known. Shared by the app and the CLI's HTML.
- **`data/utils/spatial_render.py`** (matplotlib, headless): PNG equivalents incl.
  `save_scatter3d_png` for 3D.
- **Reaction over tissue**: pick a found reaction → reactant / product / **log2(product/
  reactant)** maps across the tissue.
- **3D**: `import_imzml_stack([...])` assigns each imzML file a z-layer; a stacked
  dataset renders a **3D volumetric point cloud** (interactive HTML + PNG + a Streamlit
  "3D volume" mode).
- **CLI outputs**: `processed_dataset.json`, `ions.csv`, `reactions.csv`, and
  `images/*.png` (+ `*.html` with `--html`, sharing one `plotly.min.js`). Filenames and
  titles lead with the **molecule name** when resolved.

## 8. CLI surface — `lipidmaps-msi`

```
lipidmaps-msi FILE.imzML [FILE2.imzML ...]        # >1 file → stacked as z-slices (3D)
  --annotations MAF/CSV        # else auto-detect sibling m_*_maf.tsv, else CoreMetabolome
  --database PATH  --database-url URL             # molecule DB provisioning
  --ppm 5  --bbox X0 Y0 X1 Y1  --stride N  --top-n-peaks N
  --z-spacing UM               # section thickness for a stack
  --no-fetch-reactions  --no-headgroups
  --html  --no-images  --max-images N  --no-json  -v
```

## 9. Testing

Fully offline (network calls mocked; committed synthetic fixture):
`tests/test_spatial_utils.py`, `tests/test_mz_annotator.py`, `tests/test_msi_cli.py`,
`tests/data/test_imzml_ingestion.py` — covering formula mass, adduct m/z, DB matching,
imzML reader (annotated + auto + MAF + stride + bbox), MAF parsing, `import_imzml`,
`import_imzml_stack` (3D), candidate resolution, 2D-vs-3D detection, display names, and
the CLI (data + PNG + HTML). **Full suite: 259 passing.** Also validated live on the
demo fixture and the real MTBLS12782 `7.imzML`.

## 10. Improvements we can still make

**Correctness / science**
- **Image registration** for 3D stacks — we currently assume sections are already
  x/y-aligned; add rigid/affine alignment across slices.
- **Intensity normalization** (TIC / RMS / median per pixel, or a reference ion) — values
  are currently raw; add optional normalization and a **shared colour scale** so ion maps
  are comparable.
- **Isotope-pattern / FDR scoring** to move auto-annotation from "putative match" toward
  METASPACE-grade confidence; expose per-ion score/FDR.
- Use MAF `database_identifier` (HMDB/ChEBI/KEGG) via a cross-reference to LM IDs instead
  of relying only on name→RefMet.

**Scale / performance**
- The "pixel = sample" model (pydantic object + dict per pixel) is demo-scale; for
  full-resolution datasets move to an **array-backed store** (ion × pixel matrices,
  memory-mapped) and lazy views, keeping the object API as a façade.
- Parallelize extraction; cache parsed peak lists; optional on-disk intermediate.
- Streamlit can't take multi-GB uploads — steer large data to the CLI (done), and add a
  server-side "point at a path" option.

**Features / UX**
- **Optical/histology overlay** and anatomical-region labels (needs a co-registered
  tissue image — not in imzML).
- Region/ROI selection → **reaction z-score comparison across regions** (scaffold exists
  in `aggregate_by_region`).
- Per-slice **montage** for 3D, and true volume rendering.
- Export ion-image **matrices** (`.npz`) for downstream analysis.
- Adduct set configurable from the CLI; per-polarity presets.

**Ops**
- A stable, licensed **CoreMetabolome download** (or bundle a slim subset) so
  auto-annotation works out of the box.
- Commit a small **real** MSI slice for CI (currently synthetic); wire an integration
  test behind a network/opt-in marker.
- Docs: add an MSI page to `mkdocs` and usage examples to `README`.
