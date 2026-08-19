"""Parity test for the optimized ratio z-score computation.

`_compute_ratio_zscore` was changed from `scipy.stats.ttest_ind` / `ttest_rel` +
`stats.norm.ppf` to a direct numpy + `scipy.special.stdtr` / `ndtri`
implementation (~30x faster, no per-call distribution-object construction). This
test pins that the new code returns the same z-scores as the original scipy.stats
formula across paired/unpaired designs, both one- and two-sided alternatives, and
equal/unequal group sizes.
"""

import math

import numpy as np
import pytest
from scipy import stats

from lipidmaps.data import BioPANPathwayExporter
from lipidmaps.data.models.sample import LipidDataset, QuantifiedLipid, SampleMetadata


def make_exporter() -> BioPANPathwayExporter:
    # _compute_ratio_zscore / _select_ratio_vectors / _round_zscore do not touch
    # the dataset, so a minimal one is enough to instantiate the exporter.
    dataset = LipidDataset(
        samples=[SampleMetadata(sample_name="s1", group="a")],
        lipids=[QuantifiedLipid(input_name="PC(34:1)", values={"s1": 1.0})],
    )
    return BioPANPathwayExporter(dataset=dataset)


def reference_zscore(disease_ratios, control_ratios, alt, paired):
    """The original scipy.stats implementation, used as the oracle."""
    if len(disease_ratios) <= 1 or len(control_ratios) <= 1:
        return 0.0
    try:
        if paired:
            result = stats.ttest_rel(disease_ratios, control_ratios, alternative=alt)
        else:
            result = stats.ttest_ind(disease_ratios, control_ratios, equal_var=False, alternative=alt)
    except Exception:
        return 0.0
    p_value = getattr(result, "pvalue", None)
    if p_value is None or math.isnan(p_value) or math.isinf(p_value):
        return 0.0
    z_score = float(stats.norm.ppf(1 - p_value))
    if math.isnan(z_score):
        return 0.0
    # mirror _round_zscore
    if z_score == 0:
        return 0.0
    return float(f"{round(z_score, 3):.3f}")


def new_zscore(exporter, disease_ratios, control_ratios, alt, paired):
    # Feed ratios directly by using unit reactants, so products == ratios and
    # _select_ratio_vectors performs no dropping.
    return exporter._compute_ratio_zscore(
        list(disease_ratios),
        [1.0] * len(disease_ratios),
        list(control_ratios),
        [1.0] * len(control_ratios),
        alt=alt,
        paired=paired,
    )


ALTERNATIVES = ["greater", "less", "two-sided"]


@pytest.mark.parametrize("paired", [False, True])
@pytest.mark.parametrize("alt", ALTERNATIVES)
def test_zscore_parity_random(paired, alt):
    exporter = make_exporter()
    rng = np.random.default_rng(1234)
    max_diff = 0.0
    for _ in range(500):
        if paired:
            n = int(rng.integers(2, 9))
            n1 = n2 = n
        else:
            n1 = int(rng.integers(2, 9))
            n2 = int(rng.integers(2, 9))
        disease = rng.normal(rng.uniform(-2, 2), rng.uniform(0.2, 2.0), n1).tolist()
        control = rng.normal(rng.uniform(-2, 2), rng.uniform(0.2, 2.0), n2).tolist()
        ref = reference_zscore(disease, control, alt, paired)
        got = new_zscore(exporter, disease, control, alt, paired)
        max_diff = max(max_diff, abs(ref - got))
    # Both round to 3 decimals from z-scores that agree to ~1e-13, so outputs are
    # effectively identical (a rounding-boundary case could differ by one ulp of
    # 0.001, which this tolerance still excludes in practice).
    assert max_diff <= 1e-9, f"max z-score difference {max_diff} for paired={paired} alt={alt}"


@pytest.mark.parametrize("alt", ALTERNATIVES)
def test_zscore_parity_unequal_groups(alt):
    exporter = make_exporter()
    disease = [1.2, 0.9, 1.5, 2.1, 0.7]
    control = [0.4, 0.5, 0.3]
    assert new_zscore(exporter, disease, control, alt, False) == reference_zscore(disease, control, alt, False)


def test_zscore_too_few_samples_returns_zero():
    exporter = make_exporter()
    assert new_zscore(exporter, [1.0], [2.0, 3.0], "greater", False) == 0.0


@pytest.mark.parametrize("alt", ALTERNATIVES)
@pytest.mark.parametrize("paired", [False, True])
def test_zscore_zero_variance_matches_scipy(alt, paired):
    exporter = make_exporter()
    # Identical values within each group -> degenerate variance. The fast path
    # defers to scipy here, so the result must match the original exactly
    # (including the ±inf / NaN edge behaviour).
    disease = [1.0, 1.0, 1.0]
    control = [2.0, 2.0, 2.0]
    got = new_zscore(exporter, disease, control, alt, paired)
    ref = reference_zscore(disease, control, alt, paired)
    assert got == ref or (math.isnan(got) and math.isnan(ref))
