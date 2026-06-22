from lipidmaps.data.utils.chain_parser import infer_fa_from_lipids, parse_lipid


def test_parse_plasmalogen_pc_headgroup_from_inline_name():
    structure = parse_lipid("PC P-16:0_16:1")

    assert structure.headgroup == "PC P-"
    assert structure.linkage_type == "ether_vinyl"
    assert [chain.modifier for chain in structure.chains] == [None, None]
    assert [(chain.carbons, chain.double_bonds) for chain in structure.chains] == [(16, 0), (16, 1)]


def test_parse_plasmalogen_pe_headgroup_from_inline_name():
    structure = parse_lipid("PE P-18:0_20:4")

    assert structure.headgroup == "PE P-"
    assert structure.linkage_type == "ether_vinyl"
    assert [(chain.carbons, chain.double_bonds) for chain in structure.chains] == [(18, 0), (20, 4)]


def test_infer_fa_includes_mono_acyl_by_default():
    # Full-structure chains plus the single chain of a sum-level lyso species.
    result = infer_fa_from_lipids(["PC 16:0_18:1", "LPI 18:0"])
    assert "FA 16:0" in result and "FA 18:1" in result
    assert "FA 18:0" in result  # from the mono-acyl LPI


def test_infer_fa_excludes_mono_acyl_when_disabled():
    # Mono-acyl (lyso) species no longer contribute; full-structure chains do.
    result = infer_fa_from_lipids(["PC 16:0_18:1", "LPI 18:0"], include_mono_acyl=False)
    assert "FA 16:0" in result and "FA 18:1" in result
    assert "FA 18:0" not in result