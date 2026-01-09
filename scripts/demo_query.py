#!/usr/bin/env python3
"""Demo script for querying quantitative lipidomics CSV files.

Usage examples:
  # list samples and lipids
  python scripts/demo_query.py tests/data/inputs/small_demo.csv --list-samples --list-lipids

  # show detailed info for a lipid (input name or substring)
  python scripts/demo_query.py tests/data/inputs/small_demo.csv --lipid-info "PC(16:0/18:1)"

  # print a matrix of lipid x sample quantitative values
  python scripts/demo_query.py tests/data/inputs/small_demo.csv --table

This script uses the package's DataManager to process the CSV, so
it will run the same ingestion/RefMet annotation pipeline used
by the library.
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional

try:
    from lipidmaps.data.data_manager import DataManager
except Exception as e:
    print("Failed to import lipidmaps package. Make sure you installed it with `pip install -e .`", file=sys.stderr)
    raise


def list_samples(dataset) -> List[str]:
    return [s.sample_id for s in dataset.samples]


def list_lipids(dataset) -> List[str]:
    return [l.input_name for l in dataset.lipids]


def find_lipids(dataset, query: str):
    q = query.lower()
    return [l for l in dataset.lipids if q in (l.input_name or "").lower() or (l.standardized_name and q in l.standardized_name.lower())]


def print_lipid_info(lipid):
    recognized = bool(lipid.standardized_name or lipid.lm_id)
    print(f"Input name: {lipid.input_name}")
    print(f"Recognized: {recognized}")
    print(f"Standardized name: {lipid.standardized_name}")
    print(f"LM ID: {lipid.lm_id}")
    print(f"Values: {lipid.values}")


def print_table(dataset):
    samples = list_samples(dataset)
    header = ["Lipid"] + samples
    print("\t".join(header))
    for lipid in dataset.lipids:
        row = [lipid.input_name]
        for s in samples:
            v = lipid.values.get(s)
            row.append("" if v is None else str(v))
        print("\t".join(row))


def main(argv: Optional[List[str]] = None):
    parser = argparse.ArgumentParser(description="Query lipidomics CSV via lipidmaps DataManager")
    parser.add_argument("csv", nargs="?", default="tests/data/inputs/small_demo.csv", help="Path to CSV file")
    parser.add_argument("--list-samples", action="store_true", help="Print sample IDs found in the file")
    parser.add_argument("--list-lipids", action="store_true", help="Print input lipid names found in the file")
    parser.add_argument("--lipid-info", type=str, help="Show info for lipids matching this name or substring")
    parser.add_argument("--table", action="store_true", help="Print a tab-separated lipid x sample matrix of quantitative values")
    args = parser.parse_args(argv)

    mgr = DataManager()
    dataset = mgr.process_csv(args.csv)

    if args.list_samples:
        print("Samples:")
        for s in list_samples(dataset):
            print(" -", s)

    if args.list_lipids:
        print("Lipids:")
        for l in list_lipids(dataset):
            print(" -", l)

    if args.lipid_info:
        matches = find_lipids(dataset, args.lipid_info)
        if not matches:
            print(f"No lipids matching '{args.lipid_info}' found.")
        else:
            for m in matches:
                print("---")
                print_lipid_info(m)

    if args.table:
        print_table(dataset)


if __name__ == "__main__":
    main()
