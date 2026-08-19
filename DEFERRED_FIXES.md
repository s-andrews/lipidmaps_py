# Deferred fixes / items to revisit

Notes on suspected issues that were intentionally **not** changed (or were reverted),
pending a closer look. Each entry has enough detail to pick it back up.

---

## 1. `get_extended_facoa_names()` emits `CAR` instead of `CoA`

- **Status:** Reverted to the original `CAR` on 2026-07-01, at the user's request, to
  review in detail later. Currently the code is back to its pre-review state.
- **File:** `src/lipidmaps/data/utils/chain_parser.py`, function `get_extended_facoa_names()`.
- **Current (original) code:**
  ```python
  for carbons, double_bonds, _ in EXTENDED_FATTY_ACIDS:
      name = f"CAR {carbons}:{double_bonds}"
  ```

### Why it looks like a bug
- The function is named/documented as returning **fatty acyl-CoA** names, but it emits
  `CAR ...` (acylcarnitine), whereas its sibling `get_common_facoa_names()` emits
  `CoA {carbons}:{double_bonds}`.
- The matcher context that consumes these names only accepts acyl-CoA headgroups:
  `src/lipidmaps/data/matching/base.py` (`create_matcher_context`) keeps a parsed
  species only if `species.headgroup in ("acyl CoA", "CoA", "FACoA", "FaCoA", "FACOA")`.
- Verified with the parser (run inside the venv):
  ```
  'CAR 16:0'      -> headgroup 'CAR'        (NOT accepted -> silently dropped)
  'CoA 16:0'      -> headgroup 'acyl CoA'   (accepted)
  'FACoA 16:0'    -> headgroup 'acyl CoA'   (accepted)
  'acyl CoA 16:0' -> headgroup 'acyl CoA'   (accepted)
  ```
  So every name produced by `get_extended_facoa_names()` is dropped by the matcher and
  never enters the acyl-CoA (`facoa_species`) pool.

### The candidate fix (reverted, re-apply if confirmed correct)
```python
      name = f"CoA {carbons}:{double_bonds}"
```

### Why we're being cautious
- `CAR` is a legitimate LIPID MAPS abbreviation (acylcarnitine), so the intent may not
  be a simple typo — needs confirming whether the extended pool was meant to be
  acylcarnitines or acyl-CoAs.
- This changes the acyl-CoA pool that feeds **reaction matching**, which flows into the
  **BioPAN pathway exporter** (`biopan_pathway_exporter.py`) and its z-scores. It does
  **not** touch `biopan_exporter.py` (the display-file exporter).

### How to check the impact when we revisit
1. Grep for any downstream reliance on the literal `CAR ` prefix:
   `grep -rn "CAR " src/ scripts/` (note the space; test data legitimately contains `CAR` names).
2. Re-run the BioPAN parity/scoring tests and compare FACoA-addition reactions:
   `pytest tests/test_biopan_fa_pool_and_scoring.py tests/test_biopan_zscore_parity.py`.
3. Confirm whether `EXTENDED_FATTY_ACIDS` is meant to seed acyl-CoA or acylcarnitine names.
