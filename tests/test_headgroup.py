import sys
import types
import pytest

# Ensure package imports work during tests (use local src)
sys.path.insert(0, "src")

from lipidmaps.data.models.sample import QuantifiedLipid, LipidDataset, SampleMetadata
from lipidmaps.data.models.sample import Headgroup

# Ensure pydantic forward refs are resolved in this runtime (Pydantic v2)
for cls in (QuantifiedLipid, LipidDataset, Headgroup, SampleMetadata):
    try:
        cls.model_rebuild()
    except Exception:
        pass


def make_lipid(name, sample_vals=None, headgroup=None, generic_lm_id=None):
    return QuantifiedLipid(
        input_name=name,
        values=sample_vals or {},
        headgroup=headgroup,
        generic_lm_id=generic_lm_id,
    )


def test_generate_headgroups_creates_and_references_lipids():
    lip1 = make_lipid("PC(16:0/18:1)", {"S1": 10}, headgroup="PC", generic_lm_id="LMGP01010000")
    lip2 = make_lipid("PC(18:0/18:1)", {"S1": 5}, headgroup="PC", generic_lm_id="LMGP01010000")
    lip3 = make_lipid("TG(54:3)", {"S1": 3}, headgroup="TG", generic_lm_id="LMGL03010000")

    ds = LipidDataset(samples=[SampleMetadata(sample_name="S1", group="g")], lipids=[lip1, lip2, lip3])

    hgs = ds.generate_headgroups()

    names = {hg.name for hg in hgs}
    assert "PC" in names and "TG" in names

    pc = ds.get_headgroup("PC")
    assert pc is not None
    # ensure headgroup stores references to same QuantifiedLipid objects
    assert any(l is lip1 for l in pc.lipids)
    assert any(l is lip2 for l in pc.lipids)

    # ensure lm_ids aggregated
    assert "LMGP01010000" in pc.lm_ids


def test_headgroup_aggregation_sum_and_mean():
    lip1 = make_lipid("PC_A", {"S1": 2.0}, headgroup="PC")
    lip2 = make_lipid("PC_B", {"S1": 3.0}, headgroup="PC")
    sample = SampleMetadata(sample_name="S1", group="g")

    hg = Headgroup(name="PC", lm_ids=["L1"], lipids=[lip1, lip2])

    total = hg.aggregated_value_for_sample(sample, method="sum")
    mean = hg.aggregated_value_for_sample(sample, method="mean")

    assert pytest.approx(total, rel=1e-6) == 5.0
    assert pytest.approx(mean, rel=1e-6) == 2.5
