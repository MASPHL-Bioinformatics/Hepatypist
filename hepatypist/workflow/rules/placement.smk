#!/usr/bin/env python
"""
Per-query independent tree placement pipeline.

For each query that received a BLAST hit:
  1. Align the single query sequence against the reference alignment (MAFFT --add)
  2. Build a quick per-query tree with FastTree (-gtr -nt -fastest, no gamma)
  3. Midpoint-root the tree
  4. Extract placement statistics for that query

Individual query-ref trees are not saved; only per-query stats TSVs are kept,
then compiled into one placement_stats.tsv.

A separate combined tree (with all queries together, full GTR+gamma) is still
built for the report.html visualization.

Created 24 Oct 2024
@author: Mary Godec
@contact: mary.godec@mass.gov
"""
from hepatypist.utils.tree import *
from hepatypist.utils.blast import multifas_to_single

if config["ref_aligned"] != "":
    reffasta_aligned = config["ref_aligned"]
else:
    reffasta_aligned = RESULTS_DIR + "qc/ref_aligned.fasta"

blastparsed = os.path.join(RESULTS_DIR, "blast/blast_parsed.txt")
querymultifas = os.path.join(RESULTS_DIR, "blast/query_matches.fasta")

# Combined tree (all queries together) — used for visualization only
aligned_queryref = RESULTS_DIR + "placement/ref_query_aligned.fasta"
output_treepath = RESULTS_DIR + "placement/ref_query_tree.treefile"
rooted_treepath = RESULTS_DIR + "placement/ref_query_tree_rooted.treefile"
tree_image_path = RESULTS_DIR + "placement/ref_query_tree.png"
tree_html_path = RESULTS_DIR + "placement/ref_query_tree.html"

# Per-query placement intermediates
per_query_fas_dir = os.path.join(RESULTS_DIR, "placement/per_query_fastas/")
per_query_aln_dir = os.path.join(RESULTS_DIR, "placement/per_query_aligned/")
per_query_tree_dir = os.path.join(RESULTS_DIR, "placement/per_query_trees/")
per_query_stats_dir = os.path.join(RESULTS_DIR, "placement/per_query_stats/")

# Final compiled placement output
placement_hits = RESULTS_DIR + "placement/placement_stats.tsv"


# ── Per-query placement ───────────────────────────────────────────────────────


checkpoint split_blast_matches:
    """Split query_matches.fasta into one fasta per BLAST-passing query."""
    input:
        querymultifas,
    output:
        directory(per_query_fas_dir),
    run:
        os.makedirs(per_query_fas_dir, exist_ok=True)
        multifas_to_single(input[0], per_query_fas_dir)


def agg_per_query_stats(wildcards):
    FOUTD = checkpoints.split_blast_matches.get(**wildcards).output[0]
    (QN,) = glob_wildcards(FOUTD + "/{qn}.fasta")
    return expand(per_query_stats_dir + "/{qn}.tsv", qn=QN)


rule align_single_query:
    """Add one query to the reference alignment; preserve reference column positions."""
    input:
        qf=per_query_fas_dir + "/{qn}.fasta",
        rfa=reffasta_aligned,
    output:
        per_query_aln_dir + "/{qn}_aligned.fasta",
    shell:
        "mafft --add {input.qf} --keeplength --reorder {input.rfa} > {output}"


rule make_single_query_tree:
    """
    Build a quick per-query tree using -fastest (no gamma, no SPR rounds).
    Speed optimisation: we skip gamma rate variation and use faster heuristics
    since these trees are used only for placement stats, not for publication.
    """
    input:
        per_query_aln_dir + "/{qn}_aligned.fasta",
    output:
        temp(per_query_tree_dir + "/{qn}.treefile"),
    shell:
        "FastTree -gtr -nt -fastest {input} > {output}"


rule root_single_query_tree:
    input:
        per_query_tree_dir + "/{qn}.treefile",
    output:
        temp(per_query_tree_dir + "/{qn}_rooted.treefile"),
    run:
        root_tree(input[0], output[0])


rule place_single_query:
    input:
        rt=per_query_tree_dir + "/{qn}_rooted.treefile",
        bp=blastparsed,
    output:
        per_query_stats_dir + "/{qn}.tsv",
    params:
        support_threshold=config["high_support_threshold"],
    run:
        get_single_query_placement(
            input.rt, wildcards.qn, input.bp, output[0], params.support_threshold
        )


rule compile_placement_stats:
    """Aggregate all per-query placement TSVs into one file."""
    input:
        agg_per_query_stats,
    output:
        placement_hits,
    run:
        placement_cols = [
            "queryname",
            "querynode_support",
            "queryparent_support",
            "nn_name",
            "nn_gt",
            "nn_dist",
            "nn_support",
            "numgen_nngt",
            "monophyletic",
            "mrca_support",
            "blast_tophit",
            "blast_gt",
            "placement_gt",
            "warnings",
        ]
        if len(input) > 0:
            dfs = [pd.read_csv(f, sep="\t") for f in input]
            pd.concat(dfs, ignore_index=True).to_csv(output[0], sep="\t", index=False)
        else:
            pd.DataFrame(columns=placement_cols).to_csv(
                output[0], sep="\t", index=False
            )


# ── Combined tree for visualization ──────────────────────────────────────────


rule align_query_ref:
    """
    Align all BLAST-passing queries together with the reference (for visualization tree).
    If no queries passed BLAST (query_matches.fasta is empty), copies the reference
    alignment as-is so downstream tree rules still produce valid output.
    """
    input:
        rfa=reffasta_aligned,
        of=querymultifas,
    output:
        aligned_queryref,
    shell:
        """
        if [ -s {input.of} ]; then
            mafft --add {input.of} --keeplength --reorder {input.rfa} > {output}
        else
            cp {input.rfa} {output}
        fi
        """


rule make_tree:
    """Full GTR+gamma tree of all queries + reference, used only for report visualization."""
    input:
        aligned_queryref,
    output:
        output_treepath,
    shell:
        """
        FastTree -gtr -gamma -nt {input} > {output}
        """


rule root_tree:
    input:
        output_treepath,
    output:
        rooted_treepath,
    run:
        root_tree(output_treepath, rooted_treepath)


rule tree_viz:
    input:
        rt=rooted_treepath,
    output:
        ti=tree_image_path,
        th=tree_html_path,
    params:
        lowsupport=config["low_support_threshold"],
        highsupport=config["high_support_threshold"],
    run:
        make_treeviz(
            input.rt, output.ti, output.th, params.highsupport, params.lowsupport
        )
