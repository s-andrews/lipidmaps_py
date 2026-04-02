from lipidmaps.data.utils.chain_parser import parse_lipid


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