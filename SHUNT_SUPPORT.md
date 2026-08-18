# Shunt / Alternate-Route Reaction Support

This document describes the changes made to support **shunt (alternate-route)
reactions** across the two projects — the Python package `lipidmaps_py` and the
LIPID MAPS web application at `/lipidmaps/lipidmaps/www` — and explains the
reasoning behind each change.

## Background

A **shunt** (or *alternate route*) is a pathway branch that produces the same or
related products via a different sequence of reactions. The canonical example is
the **"Cholesterol biosynthesis (Shunt)"** branch of sterol biosynthesis. Shunt
is treated as a **general concept**, not a cholesterol-only special case.

Two capabilities were required:

1. **Flag & display** — recognize shunt reactions as a *structured attribute* and
   label them distinctly in the BioPAN pathway graph / export.
2. **Match** — enable the matching pipeline to actually pair shunt /
   sterol-biosynthesis reactions to dataset species.

## Why changes were needed in *both* projects

Two facts, verified while investigating, shaped the design:

- **The shunt designation did not reach Python.** "Shunt" existed only as a text
  suffix in the legacy `biopan_reaction.class` column (`'Cholesterol
  biosynthesis (Shunt)'`), which the modern `/api/reactions` endpoint never
  reads. The API builds each reaction's `pathways[]` from the `reaction_pathways`
  table. So without an **API/DB change**, the Python client has no shunt signal
  to act on.

- **No existing matcher could pair sterol reactions.** The `SameStructureMatcher`
  pairs a reactant to a product only when their `(total_carbons,
  total_double_bonds)` are **equal**. Sterol biosynthesis *changes composition*
  (e.g. lanosterol C30 → cholesterol C27, with desaturations/reductions), and
  sterols parse as `StructureLevel.MOLECULAR` with no acyl chains. A new,
  identity-based matcher was therefore required in the **Python** project.

Hence the work spans both repos: the API/DB emits a structured shunt signal, and
the Python client consumes it, flags it through the pipeline, and matches the
reactions.

---

## Changes in `lipidmaps_py` (Python)

### 1. Shunt detection utility — `src/lipidmaps/data/utils/shunt.py` (new)
A single, general source of truth for "is this a shunt?". `is_shunt_pathway()`
and `reaction_is_shunt()` check, in precedence order:
1. a structured boolean flag (`is_shunt` / `shunt`) — what the API now emits;
2. a structured classification (`pathway_class` / `reaction_class`);
3. a text fallback (`SHUNT_MARKERS` in the pathway name or type).

*Why:* keeping detection in one place makes the concept general and future-proof.
An explicit structured `False` is authoritative, so a curator can override a
misleading pathway name.

### 2. `is_shunt` threaded through the models
- **`models/reaction.py`** — `ReactionData.is_shunt` is derived in
  `model_post_init` from `pathways` (honoring any structured value the API sends)
  and preserved through `filter_reaction()` (which the client always applies).
- **`models/species_reaction.py`** — `ClassReaction.is_shunt` carries the flag to
  the matching layer, orthogonal to the reaction mechanism.

*Why:* the flag must survive API parsing, the lm_main filtering step, and the
match-set JSON cache, so it lives on real Pydantic fields that serialize.

### 3. Named-compound matcher (enables sterol/shunt matching)
- **`models/species_reaction.py`** — new `ReactionType.NAMED_COMPOUND`, a
  `NAMED_COMPOUND_CLASSES` sterol allow-set, and routing in `reaction_type`
  (`_is_named_compound_reaction`) that diverts sterol/shunt reactions.
- **`matching/named_compound.py`** (new) — `NamedCompoundMatcher` pairs a reactant
  and product purely by **presence/identity**, ignoring composition. Registered
  in `matching/registry.py` and exported from `matching/__init__.py`.

*Why:* sterol reactions change composition, so composition-based matching cannot
work; they must be paired by compound identity. Routing is gated so that ordinary
chain-bearing reactions are completely unaffected (they still use
`SAME_STRUCTURE`, `FA_RELEASE`, etc.).

### 4. Export / display — `biopan_pathway_exporter.py`
- `_extract_class_reactions` `register()` sets `is_shunt` (OR-ing across reactions
  that collapse into one class edge) and `is_molspecies` for named-compound edges.
- `_get_result_pathways` adds a per-pathway `is_shunt`; `_build_reaction_edges`
  adds an edge-level `is_shunt`; `_format_pathway_labels` appends a `[shunt]`
  marker (strictly gated on the flag).

*Why:* the pathway graph / frontend can now style or filter shunt edges, and the
label makes shunt reactions visible without changing non-shunt output.

### 5. Tests
- `tests/test_shunt.py` — detection precedence, true/false cases.
- `tests/test_reaction.py` — `ReactionData.is_shunt` derivation + `filter_reaction`
  preservation.
- `tests/test_named_compound_matcher.py` — routing, presence-gated pairing, ordinary
  reactions unaffected, JSON round-trip.
- `tests/test_biopan_pathway_exporter.py` — shunt edge flagged (edge-level and
  per-pathway) with the `[shunt]` label; non-shunt edges not flagged.

---

## Changes in `/lipidmaps/lipidmaps/www` (API / DB)

### 1. Migration — `laravel/database/migrations/2026_08_17_000000_add_is_shunt_to_reaction_pathways.php` (new)
- Adds `is_shunt boolean not null default false` to `reaction_pathways`.
- Seeds it by matching shunt markers in the pathway `name` / `type`.

*Why:* "shunt" is a property of a pathway (an alternate route), so
`reaction_pathways` is its faithful home. Crucially, `ReactionsController::postReactions`
attaches **whole** `reaction_pathways` rows into each reaction's `pathways[]`, so
the new column reaches clients automatically as `pathways[].is_shunt` — **no
controller code change needed.**

### 2. Schema snapshot — `laravel/database/schema/pgsql-schema.sql`
Updated the `reaction_pathways` definition to include the new `is_shunt` column,
keeping the checked-in schema in sync with the migration.

---

## Data / curation follow-up (not code)

The migration seeds `is_shunt` from text markers, but the modern
wikipathways-sourced `reaction_pathways.name` may not carry the legacy "(Shunt)"
wording. Deciding which pathways/reactions are genuinely shunt needs curator input
(the legacy `biopan_reaction` row 162 is the known cholesterol case). The schema
and plumbing are ready regardless — curators simply set `is_shunt` per row, and
the flag flows end-to-end to the Python client and BioPAN output.

## Verification

- Python: `pytest` — full suite passes, including the new shunt/matcher/exporter
  tests. `python -c "import lipidmaps.data.utils.shunt, lipidmaps.data.matching"`
  confirms no import cycle.
- www: `php -l` clean on the migration. After `php artisan migrate`, hitting
  `/api/reactions` for a sterol LM ID returns `pathways[].is_shunt`.
