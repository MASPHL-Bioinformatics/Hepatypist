#!/usr/bin/env python
# pipeline = "hepatypist"
"""
Collates output from:
    - blast genotyping top hit
    - per-query placement genotyping and statistics
    - reference statistics
    - tree visualization

Produces:
    - report.html  — interactive HTML report
    - results.tsv  — flat TSV with the same query outcome data

Created 24 Oct 2024
@author: Mary Godec
@contact: mary.godec@mass.gov
"""

from hepatypist.utils.report import *
from hepatypist.utils.tree import *

# PREPROCESSING output
queryfasta_qc = RESULTS_DIR + "qc/query_qc.txt"
preprocessing_failures = RESULTS_DIR + "qc/preprocessing_failures.tsv"

# BLAST output
blastparsed = os.path.join(RESULTS_DIR, "blast/blast_parsed.txt")
blast_failures = os.path.join(RESULTS_DIR, "blast/blast_failures.tsv")

# PLACEMENT output
aligned_queryref = RESULTS_DIR + "placement/ref_query_aligned.fasta"
tree_html_path = RESULTS_DIR + "placement/ref_query_tree.html"
placement_hits = RESULTS_DIR + "placement/placement_stats.tsv"

# REPORT output
report_html = RESULTS_DIR + "report.html"
report_tsv = RESULTS_DIR + "results.tsv"


rule make_report:
    input:
        aligned_queryref,
        placement_hits,
        ref_stats,
        tree_html_path,
        preprocessing_failures,
        blast_failures,
    output:
        report_html,
        report_tsv,
    run:
        make_report(
            aligned_queryref,
            placement_hits,
            ref_stats,
            tree_html_path,
            preprocessing_failures,
            blast_failures,
            report_html,
            report_tsv,
        )
