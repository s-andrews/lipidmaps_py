"""Unit tests for the live-vs-dev reaction-notation normalizer."""

from scripts.compare_biopan_to_live import normalize_reaction


def test_sphingo_nacyl_reduction_aligns_live_and_dev():
    # Live labels sphingolipids by N-acyl only; dev carries the full backbone.
    # Both must reduce to headgroup + N-acyl.
    assert normalize_reaction("Cer(14:0)->SM(14:0)") == "CER14:0->SM14:0"
    assert normalize_reaction("Cer 18:1;O2/14:0 -> SM 18:1;O2/14:0") == "CER14:0->SM14:0"
    assert normalize_reaction("Cer(14:0)->SM(14:0)") == normalize_reaction("Cer 18:1;O2/14:0->SM 18:1;O2/14:0")


def test_cer1p_not_confused_with_cer():
    assert normalize_reaction("Cer(16:0)->Cer1P(16:0)") == "CER16:0->CER1P16:0"


def test_ether_prefix_matches_suffix():
    # Live's O-LPC prefix form must match dev's LPC O- suffix form.
    assert normalize_reaction("O-LPC(18:0)->O-PC(40:4)") == normalize_reaction("LPC O-18:0->PC O-40:4")


def test_plain_glycerophospholipid_unchanged():
    assert normalize_reaction("PC(36:0)->PS(36:0)") == "PC36:0->PS36:0"
    assert normalize_reaction("PC 36:0->PS 36:0") == "PC36:0->PS36:0"
