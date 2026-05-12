#!/usr/bin/env python
# pipeline = "hepatypist"
"""
Performs validation on reference and query fastas.
- Are they in valid fasta format?
- Do all headers have unique names?
- are all query seqs at least 250 bp long?
- Do all headers contain only approved alphanumeric/special chars?
- Do all sequences contain only approved chars?
- Replace any spaces in headers with "_".

If using user-supplied reference fasta:
    - performs validation check
    - calculates overall and per-clade reference statistics (i.e., # sequences, # unique sequences, nucleotide diversity)
    - checks for duplicate sequences, and if so, produces a deduplicated file for alignment/tree
    - aligns user reference fasta with MAFFT

Created 24 Oct 2024
@author: Mary Godec
@contact: mary.godec@mass.gov
"""
from hepatypist.utils.preprocessing import *

# INPUT definitions
queryfasta = config["query"]
reffasta = config["ref"]

# PREPROCESSING output
queryfasta_copy = RESULTS_DIR + "qc/query.fasta"
reffasta_copy = RESULTS_DIR + "qc/ref.fasta"
queryfasta_qc = RESULTS_DIR + "qc/query_qc.txt"
query_fmt = RESULTS_DIR + "qc/query_fmt.fasta"
reffasta_qc = RESULTS_DIR + "qc/ref_qc.txt"
preprocessing_failures = RESULTS_DIR + "qc/preprocessing_failures.tsv"

if config["ref_aligned"] != "":  # if user will use default ref
    reffasta_aligned = config["ref_aligned"]
    ref_stats = config["ref_stats"]
    ref_formatted = config["ref"]

else:  # if user specifies ref and therefore ref_aligned etc are not explicitly given in config...
    reffasta_aligned = RESULTS_DIR + "qc/ref_aligned.fasta"
    ref_stats = RESULTS_DIR + "qc/refstats.tsv"
    ref_formatted = RESULTS_DIR + "qc/ref_fmt.fasta"


rule copy_input:
    input:
        queryfasta,
        reffasta,
    output:
        queryfasta_copy,
        reffasta_copy,
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """


rule test_format_query:
    input:
        queryfasta_copy,
    output:
        queryfasta_qc,
        query_fmt,
        preprocessing_failures,
    run:
        fasta_check_pythonic(
            queryfasta,
            "query",
            QUERY_MINLEN,
            query_fmt,
            queryfasta_qc,
            dedup=False,
            failoutfile=preprocessing_failures,
        )


#### if user ref supplied
rule test_format_ref:
    input:
        reffasta_copy,
    output:
        reffasta_qc,
        ref_formatted,
    run:
        fasta_check_pythonic(
            reffasta_copy,
            "ref",
            REF_MINLEN,
            ref_formatted,
            reffasta_qc,
            dedup=REF_DEDUP,
        )


rule mafft_align_ref:
    input:
        ref_formatted,
    output:
        reffasta_aligned,
    shell:
        "mafft --auto {input} > {output}"


rule calculate_stats:
    input:
        reffasta_aligned,
    output:
        ref_stats,
    run:
        calc_ref_stats(reffasta_aligned, ref_stats)
