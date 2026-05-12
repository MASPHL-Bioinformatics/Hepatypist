#!/usr/bin/env python
"""
Integration test: runs hepatypist, then compares results.tsv against a
pre-stored expected results.tsv in tests/expected/.

Compared columns: queryname, status, placement_gt, stage_failed.
Numeric/tree-derived columns (hamming pi, support values) are intentionally
excluded because they can vary slightly across runs.
"""

import sys
import argparse
import os
import pandas as pd

COMPARE_COLS = ["queryname", "status", "placement_gt", "stage_failed"]


def compare_results(expected_path, actual_path):
    exp = pd.read_csv(expected_path, sep="\t").fillna("")
    act = pd.read_csv(actual_path, sep="\t").fillna("")

    exp_queries = sorted(exp["queryname"].tolist())
    act_queries = sorted(act["queryname"].tolist())
    if exp_queries != act_queries:
        return False, (
            f"Query names differ.\n  Expected: {exp_queries}\n  Got:      {act_queries}"
        )

    cols = [c for c in COMPARE_COLS if c in exp.columns and c in act.columns]
    exp_sub = exp[cols].sort_values("queryname").reset_index(drop=True)
    act_sub = act[cols].sort_values("queryname").reset_index(drop=True)

    if exp_sub.equals(act_sub):
        return True, "OK"

    diff = exp_sub.compare(act_sub)
    return False, f"results.tsv differs from expected:\n{diff.to_string()}"


def main(sysargs=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description="Compare hepatypist results.tsv against expected output.",
        usage="python integration.py --expected <expected.tsv> --testdir <output_dir>",
    )
    parser.add_argument(
        "--expected",
        required=True,
        help="Path to the expected results.tsv (stored in tests/expected/)",
    )
    parser.add_argument(
        "--testdir",
        required=True,
        help="Path to the hepatypist output directory to evaluate",
    )

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    args = parser.parse_args(sysargs)

    actual_path = os.path.join(os.path.abspath(args.testdir), "results.tsv")
    expected_path = os.path.abspath(args.expected)

    if not os.path.isfile(actual_path):
        print(f"FAIL: results.tsv not found in {args.testdir}")
        sys.exit(1)

    if not os.path.isfile(expected_path):
        print(f"FAIL: expected file not found: {expected_path}")
        sys.exit(1)

    passed, msg = compare_results(expected_path, actual_path)
    if passed:
        print(f"PASS: results match {expected_path}")
    else:
        print(f"FAIL: {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
