# Pathway Scoring & Z-scores

`lipidmaps_py` computes two different "z-scores". They answer different questions
and are **not** interchangeable:

| Score | Where | Purpose |
| --- | --- | --- |
| **BioPAN reaction/pathway z-score** | `data/biopan_pathway_exporter.py` | The authoritative pathway-analysis score. A p-value from a t-test on the per-sample product/reactant ratio, mapped onto the standard normal. Drives BioPAN edge/pathway ranking and highlighting. |
| **`QuantitationAnalyzer.reaction_zscores`** | `data/quantitation.py` | A simple standardized difference (Cohen-like effect size). Used only by the Streamlit demo. Not a statistical test. |

There is also a **per-lipid** `QuantifiedLipid.zscore()` helper for row-wise
standardization of a single lipid's values across samples; it is unrelated to
pathway analysis and is described at the end.

---

## 1. BioPAN reaction / pathway z-score (authoritative)

This reproduces the legacy BioPAN R implementation
(`lib_pathway_analysis.r`: `get_lipid_pathway_zscore` / `extract_zscore`) so that
`lipidmaps-biopan` output matches the historical tool.

### 1.1 Inputs: per-sample product and reactant sums

For a reaction, each sample gets two numbers via `_sum_lipid_values()`:

- **product sum** = sum of the abundances of all matched product species in that sample
- **reactant sum** = sum of the abundances of all matched reactant species in that sample

Two entry points build these:

- `_score_result()` — **class level**: sums *all* matched reactant species and
  *all* matched product species for the reaction.
- `_score_pair()` — **species level**: a single reactant species and a single
  product species per pair.

Both then call `_compute_ratio_zscore(disease_products, disease_reactants,
control_products, control_reactants, alt, paired)`.

### 1.2 Building the ratio vectors — `_select_ratio_vectors()`

The per-sample score is the **ratio** `product_sum / reactant_sum` (a plain ratio,
not a log ratio). Before computing it, samples are filtered following the legacy
rules:

1. If **any** product or reactant sum is missing (`None`, `NaN`, or `Inf`), the
   reaction is unscorable and the z-score is forced to **0.0**.
2. Samples whose **reactant sum is 0** are dropped (division would be undefined):
   - *Unpaired, equal group sizes*: the **union** of zero-reactant positions is
     removed from both groups (keeps the groups aligned).
   - *Paired*: only positions where **both** groups are non-zero survive.
   - *Unpaired, unequal sizes*: each group is filtered independently.
3. If dropping leaves fewer than 2 usable samples in a group, the reaction is
   unscorable → z-score 0.0.

### 1.3 The t-test → z-score — `_compute_ratio_zscore()`

With the cleaned ratio vectors `disease` and `control`:

- If either group has ≤ 1 value → return **0.0**.
- **Paired** design (one-sample t-test on the per-sample differences):

  ```
  diff = disease - control
  t    = mean(diff) / sqrt( var(diff, ddof=1) / n )
  df   = n - 1
  ```

- **Unpaired** design (Welch's t-test, unequal variances):

  ```
  se1 = var(disease, ddof=1) / n1
  se2 = var(control, ddof=1) / n2
  t   = (mean(disease) - mean(control)) / sqrt(se1 + se2)
  df  = (se1 + se2)^2 / ( se1^2/(n1-1) + se2^2/(n2-1) )   # Welch–Satterthwaite
  ```

The t statistic is converted to a **one-sided p-value** using the Student-t CDF
(`scipy.special.stdtr`), with the direction chosen by `alt`:

```
alt == "greater":  p = 1 - stdtr(df, t)
alt == "less":     p = stdtr(df, t)
otherwise:         p = 2 * (1 - stdtr(df, |t|))   # two-sided
```

Finally the p-value is mapped onto the standard normal with the inverse normal
CDF (`scipy.special.ndtri`):

```
z = ndtri(1 - p)
```

and rounded to 3 decimals (`_round_zscore`). So a stronger, more significant
change in the expected direction yields a larger positive z.

> **Why `stdtr`/`ndtri` instead of `scipy.stats`?** The bare C ufuncs give
> identical results ~30× faster than `stats.ttest_*` + `stats.norm`, which matters
> because the score is evaluated ~12k times per view. When the variance is
> degenerate (e.g. all values within a group are equal, so `t`/`df` are non-finite),
> the code falls back to `_compute_ratio_zscore_scipy()` purely to stay bit-for-bit
> identical to the original on those rare edge cases.

### 1.4 Direction convention (mode → `alt`)

The direction is set when reaction edges are built
(`_build_reaction_edges`, `diagnose_reaction`):

```
alt = "greater"  if mode in {"active", "most_active"}  else "less"
```

- **active** pathway → test whether the disease product/reactant ratio is
  *greater* than control (positive z = reaction up-regulated in disease).
- **suppressed** → `alt = "less"`.

### 1.5 Thresholding and significance

The `threshold` parameter (default **0.05**, e.g. from `lipidmaps-biopan
--threshold 0.05`) is a **p-value** cutoff, converted to a z cutoff by
`_critical_zscore = norm.ppf(1 - threshold)` (0.05 → z ≈ 1.645). An edge is kept
as significant when its `z > critical_z` (`_select_significant_edges`). The same
cutoff filters multi-step pathway chains.

### 1.6 Multi-step pathway chains — `_chain_score()`

A pathway is a chain of reaction edges. Its score is the Stouffer combination of
the per-edge z-scores:

```
chain_z = sum(edge_z) / sqrt(number_of_edges)
```

### 1.7 Legacy substrate-consumption mode

`BioPANPathwayExporter(legacy_substrate_consumption=True)` (or `lipidmaps-biopan
--legacy-substrate-consumption`) reproduces the old tool's greedy substrate→product
pairing (`get_reaction_fa_coa` in the R code). Reactants are walked in ascending
`(total_carbons, total_double_bonds)` order; each claims still-available products,
and claimed products are removed from the pool. This changes which pairs feed
`_sum_lipid_values` and therefore the class-level z-score. It applies only to
sum-composition reactions and is cached separately (`..._legacy.json`). Leave it
off unless you specifically need parity with historical z-scores.

---

## 2. `QuantitationAnalyzer.reaction_zscores` (effect size, demo only)

Despite the name, this is **not** the BioPAN pathway z-score and is not used by
any export path — only the Streamlit demo calls it.

For each reaction it estimates a per-sample flux with
`get_reaction_flux_estimate(reaction, sample, method=...)`:

- `method="ratio"`   → `product_sum / reactant_sum`
- `method="difference"` → `product_sum - reactant_sum`

Then, for two groups, it computes a pooled-standard-deviation standardized
difference (Cohen's d style):

```
pooled_std = sqrt( ((n1-1)·s1^2 + (n2-1)·s2^2) / (n1 + n2 - 2) )   # ddof=1
z          = (mean1 - mean2) / pooled_std
```

It requires `n1 > 1` and `n2 > 1`; otherwise the score is `NaN`. Use it as a
quick effect-size indicator, not as a significance test — prefer the BioPAN
z-score (section 1) for pathway analysis.

```python
from lipidmaps import process_csv, create_analyzer

dataset = process_csv("data.csv")
analyzer = create_analyzer(dataset)
scores = analyzer.reaction_zscores("treated", "control", method="ratio")
# scores[reaction_id] -> {reaction_name, group1_mean, group2_mean, zscore, n1, n2}
```

---

## 3. Per-lipid `QuantifiedLipid.zscore()`

A row-wise helper that standardizes a single lipid's values across its samples:

```
z_s = (value_s - mean) / std        # mean/std over the non-missing values, ddof=1
```

Missing/`None` values map to `0.0`, and if fewer than two values are present or
the standard deviation is zero, all entries are `0.0`. This underlies the
`NormalizationMethod.ZSCORE` normalization and is independent of pathway analysis.
