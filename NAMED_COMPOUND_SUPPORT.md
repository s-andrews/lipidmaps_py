# Named-compound (sterol / oxylipin) pathway support

## Summary

BioPAN was silently dropping every reaction between **specific compounds that have
no generic (class) form** — sterol biosynthesis/oxidation intermediates, oxylipins
/ eicosanoids, and cholesterol ⇄ cholesteryl-ester esterification. These species
were fetched, matched to reactions, and shown as "processed" on the data page, but
produced **zero edges** in the reaction/pathway graphs.

The root cause was an identity mismatch between the reaction side and the
measured-species side of the matcher. This document describes the failure and the
fix (two phases) in `lipidmaps_py`.

## Symptom

For an input containing e.g. `cholesterol`, `27-hydroxycholesterol`,
`Arachidonic acid`, `PGD2`, `PGJ2`, `PGG2`, `CE 20:4`:

- `summary.json` correctly lists them under `pathway.processed` (cholesterol alone
  had 13 matched reactions).
- But `lp_class_reaction.json` / `lp_species_reaction.json` came back with
  `"children": []`, and the reaction graph had no sterol/oxylipin edges.

## Root cause

The pipeline joins reactions to measured species in two places that keyed the same
compound **differently**:

1. **Reaction side** — `BioPANPathwayExporter._extract_class_reactions`
   (`src/lipidmaps/data/biopan_pathway_exporter.py`) keyed each reaction component
   by its generic `compound_headgroup`. Every sterol reports
   `compound_headgroup="ST"`, so *all* distinct sterol conversions
   (cholesterol → 27-hydroxycholesterol, cholesterol → other oxysterols,
   lanosterol steps, …) **collapsed into a single `ST → ST` self-loop**.
   FA-class oxylipins report `compound_headgroup="FA"` and were **excluded
   outright** by the FA / acyl-CoA cofactor filter
   (`headgroup not in {None, "FA", "acyl CoA"}`).

2. **Measured-species side** — `create_matcher_context`
   (`src/lipidmaps/data/matching/base.py`) builds its lookup by running the
   measured **name** through `ChainParser`, keyed by the parsed headgroup. The
   measured sterols parse as `headgroup="Cholesterol"`,
   `"27-Hydroxycholesterol"`, etc. — **never `"ST"`**.

So `NamedCompoundMatcher.get_species("ST")` returned `[]` → 0 sterol pairs → empty
trees. The `NAMED_COMPOUND` machinery
(`NamedCompoundMatcher`, `NAMED_COMPOUND_CLASSES`, `is_molspecies`) already existed
and worked in isolation (see `tests/test_named_compound_matcher.py`, which feeds it
specific-name classes like `ClassReaction(reactant_class="lanosterol",
product_class="cholesterol")`) — it was simply never reached, because extraction
handed it `ST → ST`.

### Why `lm_id` is the only reliable join key

The measured compound and the reaction component share a stable identifier only via
`lm_id`. The human-readable names do **not** match:

| endpoint | measured `standardized_name` | reaction `compound_name` | `lm_id` |
|---|---|---|---|
| cholesterol | `Cholesterol` | `Cholesterol` | `LMST01010001` |
| 27-OH-cholesterol | `27-Hydroxycholesterol` | `27-hydroxy-cholesterol` | `LMST01010057` |

`compound_abbrev` is also ambiguous — five different oxysterols all report
`ST 27:1;O2`. `lm_id` is the join.

## The fix

### Phase 1 — named-to-named reactions (sterols, oxylipins)

`_extract_class_reactions` now precomputes a map of **measured specific compounds
with no generic form**:

```python
measured_by_lm_id = { lipid.lm_id: display_name
                      for lipid in dataset.lipids
                      if lipid.lm_id and not lipid.generic_lm_id }
```

Before the generic headgroup logic, a **named-compound pass** fires when a reaction
links exactly one measured specific reactant to one measured specific product
(joined by `compound_lm_id ∈ measured_by_lm_id`). It registers an identity-based
edge keyed by the **measured display names**, marks it `named_compound=True`, and
stores both `lm_id`s.

Because the class node names now equal the measured display names, they parse to the
same headgroup the matcher context uses, so `get_species(...)` resolves and the
`NamedCompoundMatcher` emits the pair. This covers both `ST → ST` sterols and the
previously-excluded `FA → FA` oxylipins uniformly.

### Phase 2 — named-compound ⇄ class via a cofactor (cholesterol ⇄ CE)

For reactions where **one** side is a measured specific compound and the other is a
chain-bearing class plus an acyl-CoA / FA cofactor (cholesterol esterification /
hydrolysis), a new helper `_component_node(component, measured_by_lm_id)` resolves a
component to `(class, lm_id)`:

- a measured specific compound → keyed by its measured display name, carrying its
  `lm_id`;
- everything else → its generic headgroup (unchanged).

The generic 1:1 registration uses it, so the sterol side is keyed as `Cholesterol`
instead of the stranded `ST`, while `CE` stays a normal class. The reaction keeps
its inferred `compound_require` (FACoA add / FA release), so the existing
full-structure matchers pair the chainless cholesterol to specific `CE` species by
the acyl-chain difference (e.g. `Cholesterol → CE 20:4` is the arachidonoyl ester).

### Deterministic routing

`ClassReaction.reaction_type` (`src/lipidmaps/data/models/species_reaction.py`)
returns `NAMED_COMPOUND` first when the new stored `named_compound` flag is set,
so routing no longer depends on the endpoint names appearing in the hard-coded
`NAMED_COMPOUND_CLASSES` set (measured names such as `27-Hydroxycholesterol` are
not in it, and would otherwise fall through to composition matching and mis-pair at
the degenerate `(0,0)` level).

## Files changed

| File | Change |
|---|---|
| `src/lipidmaps/data/models/species_reaction.py` | `ClassReaction`: added stored `named_compound` flag (forces `NAMED_COMPOUND` routing) and `reactant_lm_id` / `product_lm_id` fields. |
| `src/lipidmaps/data/biopan_pathway_exporter.py` | `_extract_class_reactions`: build `measured_by_lm_id`; add the named-compound pass (phase 1); use `_component_node` in the generic 1:1 path (phase 2). New `_component_node` helper. `register` accepts `named_compound` / lm_ids. Class-level edge builder emits the specific lm_ids for `source_lm_id` / `target_lm_id`. |
| `tests/test_biopan_pathway_exporter.py` | Added `test_named_compounds_extract_by_lm_id_not_generic_headgroup`, `test_named_compound_reactions_match_and_produce_edges`, `test_named_compound_to_class_cofactor_reaction_links_by_identity`. |

## Verification

End-to-end run of the full pipeline on a demo input produced correct directed edges
with z-scores, `lm_id`s, and node shapes:

- **Cholesterol → 27-Hydroxycholesterol** (LMST, ellipse)
- **Arachidonic acid → PGG2**, **PGD2 → PGJ2** (LMFA, triangle)
- **Cholesterol → CE 20:4 / CE 16:0** and the reverse hydrolysis edges (phase 2)

Full test suite: **232 passed**. Existing sum-composition, shunt, and z-score
parity tests are unaffected — the generic class-reaction path is unchanged except
for substituting the node identity of measured specific compounds.

## Known limitations

- **LM-abbrev sterol input names.** The join relies on the measured display name
  being unparseable so it becomes its own `ChainParser` headgroup (true for RefMet
  names like `Cholesterol`, `27-Hydroxycholesterol`, `PGD2`). If a user inputs an
  LM-style abbrev such as `ST 27:1;O`, it parses to headgroup `ST` and would not
  join. This is a pre-existing ambiguity; address only if it appears in real data.
- **Multi-substrate named reactions** where more than one measured specific
  compound appears on a side fall through to the generic path (rare).

## Operational notes

- These changes live in the **`lipidmaps_py`** repository (separate from the
  `www` site repo); they must be committed and deployed there.
- Per-session reaction match-set caches (`config/reaction_match_set*.json`) and the
  comparison cache are rebuilt on reprocess, so existing BioPAN sessions must be
  re-run to pick up the change.
