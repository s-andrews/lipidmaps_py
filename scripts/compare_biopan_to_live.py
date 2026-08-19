"""Compare dev BioPAN reaction z-scores against a live-server export TSV.

Live (legacy R) and dev use different node notations and an opposite
active/suppressed sign convention, so a raw string diff understates agreement.
This module normalizes both sides to a common key before comparing:

* ether linkage: live's prefix form ``O-LPC`` -> dev's suffix form ``LPC O-``;
* sphingolipids: live labels them by the N-acyl chain only (``Cer 14:0``) while
  dev carries the full backbone (``Cer 18:1;O2/14:0``) -> both reduce to
  headgroup + N-acyl (``CER14:0``);
* orientation: dev's sign is the data-correct one and live's export is inverted,
  so a match is ``dev == -live``.

Usage:
    python -m scripts.compare_biopan_to_live LIVE.tsv DATASET.json \
        --disease al_young --control aldr_old
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Prefer the in-repo source over any stale installed copy of the package.
_SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if _SRC_DIR.is_dir() and str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))

# Most-specific first so e.g. CER1P / DHCER are detected before CER.
SPHINGO_HEADGROUPS = ["DHCER1P", "CER1P", "DHCER", "DHSM", "HEX-CER", "GLCCER", "GALCER", "LACCER", "CER", "SM"]


def _normalize_node(raw: str) -> str:
    """Normalize a single lipid node to a comparison key."""
    s = raw.strip().replace(" ", "").replace("(", "").replace(")", "").upper()
    # Ether: prefix linkage (O-LPC) -> suffix linkage (LPCO-) to match dev.
    # The class is letters only; keep the chain digits after the moved linkage.
    s = re.sub(r"^([OP]-)([A-Z]+)", r"\2\1", s)
    # Sphingolipids: reduce to headgroup + N-acyl chain (the part after the
    # backbone), so live's N-acyl-only labels match dev's full-backbone labels.
    for headgroup in SPHINGO_HEADGROUPS:
        if s.startswith(headgroup):
            rest = s[len(headgroup):]
            if "/" in rest:
                rest = rest.split("/")[-1]
            return headgroup + rest
    return s


def normalize_reaction(reaction: str) -> str:
    """Normalize a full ``A->B`` reaction string to a comparison key."""
    return "->".join(_normalize_node(node) for node in reaction.split("->"))


def parse_live_tsv(path: Path) -> Dict[str, float]:
    """Parse single-edge reactions from a live export into ``key -> signed z``.

    Active sections contribute +z, suppressed sections -z (live's own
    convention; dev's correct convention is the negation of this).
    """
    live: Dict[str, float] = {}
    sign = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        low = line.lower()
        if low.startswith("lipid "):
            sign = 1.0 if "active" in low else (-1.0 if "suppress" in low else None)
            continue
        if line.startswith("Reactions chains") or not line.strip() or sign is None:
            continue
        parts = line.split("\t")
        if len(parts) < 2 or parts[0].count("->") != 1:
            continue
        try:
            z = float(parts[1])
        except ValueError:
            continue
        key = normalize_reaction(parts[0])
        live.setdefault(key, sign * z)
    return live


def build_dev_edges(dataset, disease: str, control: str) -> Dict[str, float]:
    """Build dev's signed reaction z-scores keyed by normalized reaction."""
    from lipidmaps.data import BioPANPathwayExporter

    exporter = BioPANPathwayExporter(dataset=dataset)
    result_set, _ = exporter.build_reaction_match_set(dataset)
    lookup = exporter._build_lipid_lookup(dataset)
    disease_samples = exporter._group_samples(dataset, disease)
    control_samples = exporter._group_samples(dataset, control)

    dev: Dict[str, float] = {}
    for result in result_set.results.values():
        if not result.has_pairs:
            continue
        rc = result.class_reaction
        dev[normalize_reaction(f"{rc.reactant_class}->{rc.product_class}")] = exporter._score_result(
            result, lookup, disease_samples, control_samples, alt="greater", paired=False
        )
        for pair in result.pairs:
            reactant = exporter._display_name_for_structure(pair.reactant, lookup)
            product = exporter._display_name_for_structure(pair.product, lookup)
            dev[normalize_reaction(f"{reactant}->{product}")] = exporter._score_pair(
                pair, lookup, disease_samples, control_samples, alt="greater", paired=False
            )
    return dev


def compare(live: Dict[str, float], dev: Dict[str, float]) -> Tuple[List, List]:
    """Return (matched, missing). Matched entries are (key, live_z, dev_z)."""
    matched = [(k, live[k], dev[k]) for k in live if k in dev]
    missing = [(k, live[k]) for k in live if k not in dev]
    return matched, missing


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("live_tsv", help="Live-server All_reaction export TSV")
    parser.add_argument("dataset_json", help="Dev processed_dataset.json (LipidDataset)")
    parser.add_argument("--disease", required=True, help="Condition-of-interest group")
    parser.add_argument("--control", required=True, help="Control group")
    args = parser.parse_args()

    from lipidmaps.data.models.sample import LipidDataset

    dataset = LipidDataset.model_validate_json(Path(args.dataset_json).read_text(encoding="utf-8"))
    live = parse_live_tsv(Path(args.live_tsv))
    dev = build_dev_edges(dataset, args.disease, args.control)
    matched, missing = compare(live, dev)

    exact = sum(1 for _, lz, dz in matched if abs((-lz) - dz) <= 0.0015)
    print(f"live single-edge reactions: {len(live)}")
    print(f"matched in dev: {len(matched)} | dev == -live (exact, |d|<=0.0015): {exact}/{len(matched)}")
    print(f"still missing: {len(missing)}")
    for key, lz in sorted(missing, key=lambda t: -abs(t[1])):
        print(f"  MISSING {key:34s} live={lz:+.3f}")


if __name__ == "__main__":
    main()
