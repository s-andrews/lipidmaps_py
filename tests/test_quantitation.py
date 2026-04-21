"""
Tests for quantitation analysis functionality.

Tests cover:
1. Reaction component quantitation
2. Normalization methods
3. Unit tracking
4. Statistical analysis
5. get_values() method
6. Group-level quantitation
"""

import math
import pytest
import numpy as np
from unittest.mock import MagicMock

from lipidmaps.data.models.sample import (
    LipidDataset,
    QuantifiedLipid,
    SampleMetadata,
)
from lipidmaps.data.models.reaction import ReactionData, CompoundComponent
from lipidmaps.data.quantitation import (
    QuantitationAnalyzer,
    QuantitationConfig,
    NormalizationMethod,
    QuantitationUnit,
    create_analyzer,
)
from lipidmaps.data_importer import LipidData


@pytest.fixture
def sample_dataset():
    """Create a sample dataset for testing."""
    samples = [
        SampleMetadata(sample_name="Control_1", group="Control"),
        SampleMetadata(sample_name="Control_2", group="Control"),
        SampleMetadata(sample_name="Control_3", group="Control"),
        SampleMetadata(sample_name="Treatment_1", group="Treatment"),
        SampleMetadata(sample_name="Treatment_2", group="Treatment"),
        SampleMetadata(sample_name="Treatment_3", group="Treatment"),
    ]
    
    lipids = [
        QuantifiedLipid(
            input_name="PC(16:0/18:1)",
            standardized_name="PC 16:0/18:1",
            lm_id="LMGP01010001",
            values={
                "Control_1": 100.0,
                "Control_2": 110.0,
                "Control_3": 105.0,
                "Treatment_1": 200.0,
                "Treatment_2": 210.0,
                "Treatment_3": 195.0,
            },
        ),
        QuantifiedLipid(
            input_name="PE(18:0/20:4)",
            standardized_name="PE 18:0/20:4",
            lm_id="LMGP02010001",
            values={
                "Control_1": 50.0,
                "Control_2": 55.0,
                "Control_3": 52.0,
                "Treatment_1": 45.0,
                "Treatment_2": 48.0,
                "Treatment_3": 47.0,
            },
        ),
        QuantifiedLipid(
            input_name="TAG(16:0/18:1/18:2)",
            standardized_name="TAG 16:0/18:1/18:2",
            lm_id="LMGL03010001",
            generic_lm_id="LMGL03000000",
            values={
                "Control_1": 500.0,
                "Control_2": 520.0,
                "Control_3": 510.0,
                "Treatment_1": 550.0,
                "Treatment_2": 560.0,
                "Treatment_3": 545.0,
            },
        ),
        QuantifiedLipid(
            input_name="IS_PC(15:0/15:0)",  # Internal standard
            standardized_name="PC 15:0/15:0",
            values={
                "Control_1": 100.0,
                "Control_2": 100.0,
                "Control_3": 100.0,
                "Treatment_1": 100.0,
                "Treatment_2": 100.0,
                "Treatment_3": 100.0,
            },
        ),
    ]
    
    return LipidDataset(samples=samples, lipids=lipids)


@pytest.fixture
def analyzer(sample_dataset):
    """Create a QuantitationAnalyzer for testing."""
    return QuantitationAnalyzer(dataset=sample_dataset)


@pytest.fixture
def lipid_data(sample_dataset):
    """Create a LipidData object for testing."""
    return LipidData(dataset=sample_dataset)


class TestCreateAnalyzer:
    """Test factory function."""
    
    def test_create_analyzer(self, sample_dataset):
        """Test create_analyzer factory function."""
        analyzer = create_analyzer(sample_dataset)
        assert isinstance(analyzer, QuantitationAnalyzer)
        assert analyzer.dataset == sample_dataset


class TestQuantitationConfig:
    """Test quantitation configuration."""
    
    def test_default_config(self, analyzer):
        """Test default configuration values."""
        config = analyzer.config
        assert config.unit is None
        assert config.method is None
        assert config.normalization == NormalizationMethod.NONE
        assert config.internal_standard_lipid is None
    
    def test_set_config(self, analyzer):
        """Test setting configuration values."""
        analyzer.set_config(
            unit="pmol",
            method="LC-MS",
            notes="Test dataset",
        )
        assert analyzer.config.unit == "pmol"
        assert analyzer.config.method == "LC-MS"
        assert analyzer.config.notes == "Test dataset"


class TestNormalizationMethods:
    """Test normalization methods."""
    
    def test_normalize_total_lipid(self, analyzer):
        """Test total lipid normalization."""
        result = analyzer.normalize_total_lipid(scale_factor=100.0)
        
        # Check that normalization was applied
        assert "PC(16:0/18:1)" in result
        assert "Control_1" in result["PC(16:0/18:1)"]
        
        # Control_1 total = 100 + 50 + 500 + 100 = 750
        # PC value = 100 / 750 * 100 = 13.33...
        expected = (100.0 / 750.0) * 100.0
        assert abs(result["PC(16:0/18:1)"]["Control_1"] - expected) < 0.01
    
    def test_normalize_internal_standard(self, analyzer):
        """Test internal standard normalization."""
        result = analyzer.normalize_internal_standard("IS_PC(15:0/15:0)")
        
        # Internal standard is 100 in all samples, so values should be unchanged
        assert result["PC(16:0/18:1)"]["Control_1"] == 1.0  # 100/100
        assert result["TAG(16:0/18:1/18:2)"]["Control_1"] == 5.0  # 500/100
    
    def test_normalize_internal_standard_not_found(self, analyzer):
        """Test error when internal standard not found."""
        with pytest.raises(ValueError, match="not found"):
            analyzer.normalize_internal_standard("NonExistent")
    
    def test_normalize_log2(self, analyzer):
        """Test log2 transformation."""
        result = analyzer.normalize_log(base=2, offset=0.0)
        
        # log2(100) ≈ 6.644
        assert abs(result["PC(16:0/18:1)"]["Control_1"] - math.log2(100)) < 0.01
    
    def test_normalize_log10(self, analyzer):
        """Test log10 transformation."""
        result = analyzer.normalize_log(base=10, offset=0.0)
        
        # log10(100) = 2
        assert abs(result["PC(16:0/18:1)"]["Control_1"] - 2.0) < 0.01
    
    def test_normalize_log_with_offset(self, analyzer):
        """Test log transformation with offset."""
        result = analyzer.normalize_log(base=2, offset=1.0)
        
        # log2(100 + 1) = log2(101)
        assert abs(result["PC(16:0/18:1)"]["Control_1"] - math.log2(101)) < 0.01
    
    def test_normalize_median_center(self, analyzer):
        """Test median centering."""
        result = analyzer.normalize_median_center()
        
        # Check that values are centered (median subtracted)
        assert "PC(16:0/18:1)" in result
        # The median of Control_1 values (100, 50, 500, 100) is 100
        # So PC value should be 100 - 100 = 0
        assert abs(result["PC(16:0/18:1)"]["Control_1"] - 0.0) < 0.01
    
    def test_apply_normalization_in_place(self, sample_dataset):
        """Test in-place normalization."""
        analyzer = QuantitationAnalyzer(dataset=sample_dataset)
        original_value = sample_dataset.lipids[0].values["Control_1"]
        
        analyzer.apply_normalization(
            method=NormalizationMethod.LOG2,
            log_offset=0.0,
            in_place=True,
        )
        
        new_value = sample_dataset.lipids[0].values["Control_1"]
        expected = math.log2(original_value)
        assert abs(new_value - expected) < 0.01


class TestStatisticalAnalysis:
    """Test statistical analysis methods."""
    
    def test_calculate_fold_change(self, analyzer):
        """Test fold change calculation."""
        fc = analyzer.calculate_fold_change("Treatment", "Control", log2=True)
        
        # PC: Treatment mean ~202, Control mean ~105, log2 FC ~0.94
        assert "PC(16:0/18:1)" in fc
        assert fc["PC(16:0/18:1)"] > 0  # Up in treatment
        
        # PE: Treatment mean ~47, Control mean ~52, log2 FC ~ -0.15
        assert fc["PE(18:0/20:4)"] < 0  # Down in treatment
    
    def test_calculate_fold_change_non_log(self, analyzer):
        """Test non-log fold change."""
        fc = analyzer.calculate_fold_change("Treatment", "Control", log2=False)
        
        # PC: Treatment/Control ≈ 202/105 ≈ 1.92
        assert 1.8 < fc["PC(16:0/18:1)"] < 2.1
    
    def test_calculate_pvalue(self, analyzer):
        """Test p-value calculation."""
        pvals = analyzer.calculate_pvalue("Treatment", "Control", test="ttest")
        
        assert "PC(16:0/18:1)" in pvals
        # PC should be significantly different (large difference)
        assert pvals["PC(16:0/18:1)"] < 0.05
    
    def test_calculate_pvalue_mannwhitney(self, analyzer):
        """Test Mann-Whitney U test."""
        pvals = analyzer.calculate_pvalue(
            "Treatment", "Control", test="mannwhitney"
        )
        
        assert "PC(16:0/18:1)" in pvals
        assert not math.isnan(pvals["PC(16:0/18:1)"])
    
    def test_calculate_cv(self, analyzer):
        """Test CV calculation."""
        cvs = analyzer.calculate_cv(by_group=False)
        
        # Check all lipids have CV values
        assert len(cvs) == 4
        assert all(not math.isnan(v) for v in cvs.values())
        
        # Internal standard should have very low CV (all same values)
        assert cvs["IS_PC(15:0/15:0)"] == 0.0
    
    def test_calculate_cv_by_group(self, analyzer):
        """Test CV calculation by group."""
        cvs = analyzer.calculate_cv(by_group=True)
        
        assert "Control" in cvs
        assert "Treatment" in cvs
        assert "PC(16:0/18:1)" in cvs["Control"]
    
    def test_differential_analysis(self, analyzer):
        """Test differential analysis."""
        results = analyzer.differential_analysis(
            "Treatment",
            "Control",
            fc_threshold=0.5,
            pvalue_threshold=0.05,
        )
        
        # Should return list of dicts
        assert isinstance(results, list)
        assert len(results) == 4  # All lipids
        
        # Results should be sorted by p-value
        pvals = [r["p_value"] for r in results if not math.isnan(r["p_value"])]
        assert pvals == sorted(pvals)
        
        # Check result structure
        pc_result = next(r for r in results if r["lipid_name"] == "PC(16:0/18:1)")
        assert "log2_fold_change" in pc_result
        assert "p_value" in pc_result
        assert "significant" in pc_result
        assert "direction" in pc_result
        assert pc_result["direction"] == "up"


class TestGroupLevelQuantitation:
    """Test group-level quantitation methods."""
    
    def test_get_group_means(self, analyzer):
        """Test group mean calculation."""
        means = analyzer.get_group_means()
        
        assert "Control" in means
        assert "Treatment" in means
        
        # PC Control mean: (100 + 110 + 105) / 3 = 105
        assert abs(means["Control"]["PC(16:0/18:1)"] - 105.0) < 0.01
    
    def test_get_group_stds(self, analyzer):
        """Test group std calculation."""
        stds = analyzer.get_group_stds()
        
        assert "Control" in stds
        assert "Treatment" in stds
        
        # Check that internal standard has 0 std
        assert stds["Control"]["IS_PC(15:0/15:0)"] == 0.0
    
    def test_get_group_summary(self, analyzer):
        """Test comprehensive group summary."""
        summary = analyzer.get_group_summary()
        
        assert "Control" in summary
        pc_control = summary["Control"]["PC(16:0/18:1)"]
        
        assert "mean" in pc_control
        assert "std" in pc_control
        assert "min" in pc_control
        assert "max" in pc_control
        assert "median" in pc_control
        assert "n" in pc_control
        
        assert pc_control["n"] == 3
        assert pc_control["min"] == 100.0
        assert pc_control["max"] == 110.0
    
    def test_compare_groups(self, analyzer):
        """Test group comparison."""
        comparison = analyzer.compare_groups("Treatment", "Control")
        
        pc_comp = comparison["PC(16:0/18:1)"]
        assert "Treatment_mean" in pc_comp
        assert "Control_mean" in pc_comp
        assert "Treatment_std" in pc_comp
        assert "Control_std" in pc_comp
        assert "log2_fold_change" in pc_comp
        assert "p_value" in pc_comp
        assert "significant" in pc_comp


class TestReactionComponentQuantitation:
    """Test reaction component quantitation."""
    
    def test_get_value_for_reaction_component(self, analyzer, sample_dataset):
        """Test getting value for reaction component by LM ID."""
        # Use the PC lipid's LM ID
        value = analyzer.get_value_for_reaction_component(
            "LMGP01010001",
            "Control_1",
            method="sum",
        )
        
        assert value == 100.0  # PC value in Control_1
    
    def test_get_value_for_reaction_component_by_generic_id(self, analyzer):
        """Test getting value for reaction component by generic LM ID."""
        # TAG has generic_lm_id = LMGL03000000
        value = analyzer.get_value_for_reaction_component(
            "LMGL03000000",
            "Control_1",
            method="sum",
        )
        
        assert value == 500.0  # TAG value in Control_1
    
    def test_get_value_for_reaction_component_not_found(self, analyzer):
        """Test when component not found."""
        value = analyzer.get_value_for_reaction_component(
            "NONEXISTENT",
            "Control_1",
            method="sum",
        )
        
        assert value is None
    
    def test_get_values_for_reaction_component(self, analyzer):
        """Test getting all sample values for reaction component."""
        values = analyzer.get_values_for_reaction_component(
            "LMGP01010001",
            method="sum",
        )
        
        assert "Control_1" in values
        assert "Treatment_1" in values
        assert values["Control_1"] == 100.0
        assert values["Treatment_1"] == 200.0
    
    def test_get_value_with_compound_component(self, analyzer):
        """Test with CompoundComponent object."""
        component = CompoundComponent(compound_lm_id="LMGP01010001")
        value = analyzer.get_value_for_reaction_component(
            component,
            "Control_1",
            method="sum",
        )
        
        assert value == 100.0
    
    def test_get_reaction_flux_estimate(self, analyzer, sample_dataset):
        """Test reaction flux estimation."""
        # Create a simple reaction
        reaction = ReactionData(
            reaction_id=1,
            reactants=[CompoundComponent(compound_lm_id="LMGP01010001")],  # PC
            products=[CompoundComponent(compound_lm_id="LMGP02010001")],  # PE
        )
        
        flux = analyzer.get_reaction_flux_estimate(
            reaction,
            "Control_1",
            method="ratio",
        )
        
        # PE/PC = 50/100 = 0.5
        assert flux == 0.5


class TestUnitTracking:
    """Test unit tracking functionality."""
    
    def test_lipid_unit_fields(self):
        """Test that QuantifiedLipid has unit fields."""
        lipid = QuantifiedLipid(
            input_name="Test",
            values={"S1": 100.0},
            unit="pmol",
            measurement_method="LC-MS",
            quantitation_notes="Test notes",
        )
        
        assert lipid.unit == "pmol"
        assert lipid.measurement_method == "LC-MS"
        assert lipid.quantitation_notes == "Test notes"
    
    def test_dataset_quantitation_info(self, sample_dataset):
        """Test setting quantitation info on dataset."""
        sample_dataset.set_quantitation_info(
            unit="ng",
            method="GC-MS",
            notes="Sample dataset",
        )
        
        assert sample_dataset.quantitation_unit == "ng"
        assert sample_dataset.quantitation_method == "GC-MS"
        assert sample_dataset.quantitation_notes == "Sample dataset"
    
    def test_dataset_apply_to_lipids(self, sample_dataset):
        """Test applying quantitation info to all lipids."""
        sample_dataset.set_quantitation_info(
            unit="pmol",
            method="LC-MS",
            apply_to_lipids=True,
        )
        
        for lipid in sample_dataset.lipids:
            assert lipid.unit == "pmol"
            assert lipid.measurement_method == "LC-MS"


class TestGetValues:
    """Test get_values method."""
    
    def test_get_values_from_dataset(self, sample_dataset):
        """Test get_values on LipidDataset."""
        values = sample_dataset.get_values("PC(16:0/18:1)")
        
        assert values is not None
        assert "Control_1" in values
        assert values["Control_1"] == 100.0
    
    def test_get_values_case_insensitive(self, sample_dataset):
        """Test case-insensitive matching."""
        values = sample_dataset.get_values("pc(16:0/18:1)")
        
        assert values is not None
        assert "Control_1" in values
    
    def test_get_values_not_found(self, sample_dataset):
        """Test when lipid not found."""
        values = sample_dataset.get_values("NonExistent")
        
        assert values is None
    
    def test_get_values_from_lipid_data(self, lipid_data):
        """Test get_values on LipidData."""
        values = lipid_data.get_values("PC(16:0/18:1)")
        
        assert values is not None
        assert "Control_1" in values


class TestLipidDataIntegration:
    """Test integration with LipidData class."""
    
    def test_lipid_data_analyzer_property(self, lipid_data):
        """Test analyzer property on LipidData."""
        assert hasattr(lipid_data, "analyzer")
        assert isinstance(lipid_data.analyzer, QuantitationAnalyzer)

    def test_lipid_data_analyzer_is_cached(self, lipid_data):
        """Test analyzer property caches a single analyzer instance."""
        first = lipid_data.analyzer
        second = lipid_data.analyzer

        assert first is second
    
    def test_lipid_data_normalize(self, lipid_data):
        """Test normalize method on LipidData."""
        result = lipid_data.normalize("log2")
        
        assert "PC(16:0/18:1)" in result
    
    def test_lipid_data_fold_change(self, lipid_data):
        """Test fold change through LipidData."""
        fc = lipid_data.calculate_fold_change("Treatment", "Control")
        
        assert "PC(16:0/18:1)" in fc
        assert fc["PC(16:0/18:1)"] > 0
    
    def test_lipid_data_differential_analysis(self, lipid_data):
        """Test differential analysis through LipidData."""
        results = lipid_data.differential_analysis(
            "Treatment", "Control"
        )
        
        assert len(results) == 4
    
    def test_lipid_data_group_means(self, lipid_data):
        """Test group means through LipidData."""
        means = lipid_data.get_group_means()
        
        assert "Control" in means
        assert "Treatment" in means
    
    def test_lipid_data_set_quantitation_config(self, lipid_data):
        """Test setting quantitation config through LipidData."""
        lipid_data.set_quantitation_config(
            unit="pmol",
            method="LC-MS",
        )
        
        assert lipid_data.quantitation_config.unit == "pmol"
        assert lipid_data.quantitation_config.method == "LC-MS"



class TestQuantitationUnit:
    """Test QuantitationUnit enum."""
    
    def test_unit_values(self):
        """Test all unit enum values."""
        assert QuantitationUnit.PMOL.value == "pmol"
        assert QuantitationUnit.NMOL.value == "nmol"
        assert QuantitationUnit.NG.value == "ng"
        assert QuantitationUnit.AREA.value == "area"
        assert QuantitationUnit.RELATIVE.value == "relative"
        assert QuantitationUnit.PERCENT.value == "percent"


class TestNormalizationMethodEnum:
    """Test NormalizationMethod enum."""
    
    def test_method_values(self):
        """Test all normalization method values."""
        assert NormalizationMethod.NONE.value == "none"
        assert NormalizationMethod.TOTAL_LIPID.value == "total_lipid"
        assert NormalizationMethod.INTERNAL_STANDARD.value == "internal_standard"
        assert NormalizationMethod.LOG2.value == "log2"
        assert NormalizationMethod.LOG10.value == "log10"
        assert NormalizationMethod.MEDIAN_CENTER.value == "median_center"
        assert NormalizationMethod.ZSCORE.value == "zscore"


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_dataset(self):
        """Test with empty dataset."""
        dataset = LipidDataset(samples=[], lipids=[])
        analyzer = QuantitationAnalyzer(dataset=dataset)
        
        # Should not raise errors
        means = analyzer.get_group_means()
        assert means == {}
        
        fc = analyzer.calculate_fold_change("A", "B")
        assert fc == {}
    
    def test_single_sample_group(self):
        """Test with single sample per group."""
        samples = [
            SampleMetadata(sample_name="A", group="GroupA"),
            SampleMetadata(sample_name="B", group="GroupB"),
        ]
        lipids = [
            QuantifiedLipid(
                input_name="L1",
                values={"A": 100.0, "B": 200.0},
            ),
        ]
        dataset = LipidDataset(samples=samples, lipids=lipids)
        analyzer = QuantitationAnalyzer(dataset=dataset)
        
        # Should handle single-sample groups gracefully
        fc = analyzer.calculate_fold_change("GroupB", "GroupA")
        assert fc["L1"] == 1.0  # log2(200/100) = 1
    
    def test_missing_values(self):
        """Test handling of missing values."""
        samples = [
            SampleMetadata(sample_name="A", group="G"),
            SampleMetadata(sample_name="B", group="G"),
        ]
        lipids = [
            QuantifiedLipid(
                input_name="L1",
                values={"A": 100.0},  # Missing B
            ),
        ]
        dataset = LipidDataset(samples=samples, lipids=lipids)
        analyzer = QuantitationAnalyzer(dataset=dataset)
        
        # Should handle missing gracefully
        means = analyzer.get_group_means()
        assert means["G"]["L1"] == 100.0  # Only has A

    def test_zero_values(self):
        """Test handling of zero values in log transformation."""
        samples = [SampleMetadata(sample_name="A", group="G")]
        lipids = [
            QuantifiedLipid(
                input_name="L1",
                values={"A": 0.0},
            ),
        ]
        dataset = LipidDataset(samples=samples, lipids=lipids)
        analyzer = QuantitationAnalyzer(dataset)
        
        # Without offset, log(0) should give -inf
        result = analyzer.normalize_log(base=2, offset=0.0)
        assert result["L1"]["A"] == float('-inf')
        
        # With offset, should be log(0 + 1) = 0
        result = analyzer.normalize_log(base=2, offset=1.0)
        assert result["L1"]["A"] == 0.0
