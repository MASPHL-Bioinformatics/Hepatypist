#!/usr/bin/env python
import pandas as pd
from Bio import SeqIO


# ── Sequence utilities ────────────────────────────────────────────────────────


def strip_seq(s):
    return str(s).replace("\n", "").replace("\r", "")


def hamming_distance(s1, s2) -> int:
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def no_gap_hamming_distance(s1, s2) -> int:
    return sum(el1 != el2 for el1, el2 in zip(s1, s2) if "-" not in (el1, el2))


# ── Hamming pi per-query stats ────────────────────────────────────────────────


def per_query_medoid_stats(querygt, queryh, rsdfsubs, qra, is_single_ref=False):
    """
    Computes query hamming pi vs. genotype medoid.

    For single-ref genotypes the medoid ID is the one reference sequence itself
    (set in preprocessing.gt_metrics).  We still compute the distance but mark
    the genotype range as N/A and flag as SINGLE_REF so the user can see the
    raw query value without a spurious WITHIN/OUTSIDE verdict.
    """
    hammedoid_id = rsdfsubs.iloc[0]["gt_ham_medoid"]
    nghammedoid_id = rsdfsubs.iloc[0]["gt_ng_ham_medoid"]

    if hammedoid_id not in qra or nghammedoid_id not in qra or queryh not in qra:
        return {
            "queryname": queryh,
            "clade": querygt,
            "query_ham_pi": "",
            "gt_ham_pi_range": "",
            "ham_pi_flag": "",
            "query_ng_ham_pi": "",
            "gt_ng_ham_pi_range": "",
            "ng_ham_pi_flag": "",
        }

    queryseq = qra[queryh]
    hammedoid_seq = qra[hammedoid_id]
    nghammedoid_seq = qra[nghammedoid_id]
    avlen = len(queryseq)

    qcham = hamming_distance(queryseq, hammedoid_seq)
    qcngham = no_gap_hamming_distance(queryseq, nghammedoid_seq)

    query_ham_pi = "{:.3f}".format(float(qcham / avlen))
    query_ng_ham_pi = "{:.3f}".format(float(qcngham / avlen))

    if is_single_ref:
        gt_ham_pi_range = "N/A (single ref)"
        gt_ng_ham_pi_range = "N/A (single ref)"
        ham_pi_flag = "SINGLE_REF"
        ng_ham_pi_flag = "SINGLE_REF"
    else:
        gt_avg = str(rsdfsubs.iloc[0]["gt_avg_ham_pi"])
        gt_max = str(rsdfsubs.iloc[0]["max_ham_pi"])
        gt_ng_avg = str(rsdfsubs.iloc[0]["gt_avg_nogaps_ham_pi"])
        gt_ng_max = str(rsdfsubs.iloc[0]["max_nogaps_ham_pi"])
        gt_ham_pi_range = f"avg {gt_avg} / max {gt_max}"
        gt_ng_ham_pi_range = f"avg {gt_ng_avg} / max {gt_ng_max}"
        try:
            ham_pi_flag = "OUTSIDE" if float(query_ham_pi) > float(gt_max) else "WITHIN"
            ng_ham_pi_flag = (
                "OUTSIDE" if float(query_ng_ham_pi) > float(gt_ng_max) else "WITHIN"
            )
        except (ValueError, TypeError):
            ham_pi_flag = ""
            ng_ham_pi_flag = ""

    return {
        "queryname": queryh,
        "clade": querygt,
        "query_ham_pi": query_ham_pi,
        "gt_ham_pi_range": gt_ham_pi_range,
        "ham_pi_flag": ham_pi_flag,
        "query_ng_ham_pi": query_ng_ham_pi,
        "gt_ng_ham_pi_range": gt_ng_ham_pi_range,
        "ng_ham_pi_flag": ng_ham_pi_flag,
    }


# ── Build combined placement + hamming-pi dataframe ──────────────────────────

_PLACEMENT_RESULT_COLS = [
    "queryname",
    "placement_gt",
    "blast_gt",
    "query_ham_pi",
    "gt_ham_pi_range",
    "ham_pi_flag",
    "query_ng_ham_pi",
    "gt_ng_ham_pi_range",
    "ng_ham_pi_flag",
    "nn_gt",
    "nn_support",
    "querynode_support",
    "queryparent_support",
    "mrca_support",
    "monophyletic",
    "warnings",
]


def make_per_gt_df(queryrefaln, placement_stats, refdf):
    """
    Joins placement stats with per-query hamming pi metrics vs. the assigned
    genotype medoid.  Always calls per_query_medoid_stats; single-ref genotypes
    are handled inside that function.
    """
    if len(placement_stats) == 0:
        return pd.DataFrame(columns=_PLACEMENT_RESULT_COLS)

    qra = {r.description: strip_seq(r.seq) for r in SeqIO.parse(queryrefaln, "fasta")}
    queries = list(placement_stats["queryname"])
    ref_gts = list(set(refdf["clade"]))

    meddflist = []
    for q in queries:
        qgt = placement_stats[placement_stats["queryname"] == q].iloc[0]["placement_gt"]
        if qgt in ref_gts:
            refdfgt = refdf[refdf["clade"] == qgt]
            refdfsubs = refdfgt[
                [
                    "clade",
                    "ham_pi",
                    "max_ham_pi",
                    "nogaps_ham_pi",
                    "max_nogaps_ham_pi",
                    "num_seqs_gt",
                    "gt_ham_medoid",
                    "gt_ng_ham_medoid",
                ]
            ].rename(
                columns={
                    "ham_pi": "gt_avg_ham_pi",
                    "nogaps_ham_pi": "gt_avg_nogaps_ham_pi",
                }
            )
            is_single_ref = int(refdfsubs["num_seqs_gt"].iloc[0]) == 1
            meddict = per_query_medoid_stats(
                qgt, q, refdfsubs, qra, is_single_ref=is_single_ref
            )
        else:
            meddict = {
                "queryname": q,
                "clade": qgt,
                "query_ham_pi": "",
                "gt_ham_pi_range": "",
                "ham_pi_flag": "",
                "query_ng_ham_pi": "",
                "gt_ng_ham_pi_range": "",
                "ng_ham_pi_flag": "",
            }
        meddflist.append(meddict)

    medcomb = pd.DataFrame(meddflist)
    combined = pd.merge(placement_stats, medcomb, on="queryname", how="left")

    return combined[
        [
            "queryname",
            "placement_gt",
            "blast_gt",
            "query_ham_pi",
            "gt_ham_pi_range",
            "ham_pi_flag",
            "query_ng_ham_pi",
            "gt_ng_ham_pi_range",
            "ng_ham_pi_flag",
            "nn_gt",
            "nn_support",
            "querynode_support",
            "queryparent_support",
            "mrca_support",
            "monophyletic",
            "warnings",
        ]
    ].fillna("")


# ── Build full report dataframe (placements + failures) ───────────────────────

_ALL_COLS = [
    "queryname",
    "status",
    "stage_failed",
    "failure_reason",
    "placement_gt",
    "blast_gt",
    "query_ham_pi",
    "gt_ham_pi_range",
    "ham_pi_flag",
    "query_ng_ham_pi",
    "gt_ng_ham_pi_range",
    "ng_ham_pi_flag",
    "nn_gt",
    "nn_support",
    "querynode_support",
    "queryparent_support",
    "mrca_support",
    "monophyletic",
    "warnings",
]


def _build_full_report(placement_df, preproc_failures, blast_failures):
    rows = []

    for _, row in placement_df.iterrows():
        pgt = row.get("placement_gt", "")
        status = "UNASSIGNABLE" if pgt == "UNASSIGNABLE" else "SUCCESS"
        r = {"status": status, "stage_failed": "", "failure_reason": ""}
        r.update({c: row.get(c, "") for c in _ALL_COLS if c not in r})
        rows.append(r)

    for _, row in preproc_failures.iterrows():
        rows.append(
            {
                "queryname": row["queryname"],
                "status": "FAILED",
                "stage_failed": "preprocessing",
                "failure_reason": row["reason"],
                **{
                    c: ""
                    for c in _ALL_COLS
                    if c
                    not in ("queryname", "status", "stage_failed", "failure_reason")
                },
            }
        )

    for _, row in blast_failures.iterrows():
        rows.append(
            {
                "queryname": row["queryname"],
                "status": "FAILED",
                "stage_failed": "blast",
                "failure_reason": row["reason"],
                **{
                    c: ""
                    for c in _ALL_COLS
                    if c
                    not in ("queryname", "status", "stage_failed", "failure_reason")
                },
            }
        )

    return pd.DataFrame(rows, columns=_ALL_COLS).fillna("")


# ── HTML generation ──────────────────────────────────────────────────────────


def _style_cell(val):
    if val in ("SUCCESS", "WITHIN"):
        return "background-color:#d4edda;color:#155724;font-weight:bold"
    if val in ("FAILED",):
        return "background-color:#f8d7da;color:#721c24;font-weight:bold"
    if val in ("UNASSIGNABLE", "OUTSIDE"):
        return "background-color:#fff3cd;color:#856404;font-weight:bold"
    if val == "SINGLE_REF":
        return "background-color:#cce5ff;color:#004085;font-style:italic"
    return ""


def _df_to_styled_html(df):
    style_cols = [
        c for c in ("status", "ham_pi_flag", "ng_ham_pi_flag") if c in df.columns
    ]
    try:
        styled = df.style.map(_style_cell, subset=style_cols)
    except AttributeError:
        # pandas < 2.1 uses applymap
        styled = df.style.applymap(_style_cell, subset=style_cols)
    try:
        styled = styled.hide(axis="index")
    except AttributeError:
        # pandas < 1.4
        styled = styled.hide_index()
    return styled.to_html()


def generate_combined_interactive_html(full_df, ref_df, plot_html_path, output_path):
    with open(plot_html_path, "r") as f:
        plot_html = f.read()

    table1_html = _df_to_styled_html(full_df)
    table2_html = ref_df.to_html(index=False)

    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <style>
            body {{
                display: flex;
                flex-direction: row;
                font-family: Arial, sans-serif;
                padding: 0;
                margin: 0;
                overflow: hidden;
            }}
            .left {{
                width: 60%;
                overflow-y: scroll;
                padding: 10px;
                max-height: 100vh;
                margin-right: 30px;
                box-sizing: border-box;
                font-size: 0.75rem;
            }}
            .right {{
                width: 40%;
                overflow-y: auto;
                padding: 30px;
                margin-left: 50px;
                max-height: 100vh;
                box-sizing: border-box;
                border: 50px;
                border-style: solid;
                border-color: black;
                border-top: 0;
                border-bottom: 0;
                border-right: 0;
            }}
            .right img {{
                width: 100%;
                height: auto;
            }}
            table {{
                font-size: 0.75rem;
                border-collapse: collapse;
                border-spacing: 10px;
                margin-bottom: 20px;
            }}
            thead th {{
                background: #88CCF1;
                color: #FFF;
                font-family: 'Lato', sans-serif;
                font-size: 0.75rem;
                font-weight: 100;
                text-transform: uppercase;
            }}
            th, td {{
                border: 1px solid black;
                padding: 6px 10px;
                text-align: center;
            }}
        </style>
    </head>
    <body>
        <div class="left">
            <h2>Query genotype assignment and placement statistics</h2>
            <p style="font-size:0.7rem;color:#555">
                <b>ham_pi_flag</b>: WITHIN = query hamming pi ≤ genotype max;
                OUTSIDE = exceeds genotype max (potential outlier);
                SINGLE_REF = only one reference sequence for this genotype.
            </p>
            {table1_html}
            <h2>Per-genotype reference statistics</h2>
            {table2_html}
        </div>
        <div class="right">
            {plot_html}
        </div>
    </body>
    </html>
    """

    with open(output_path, "w") as f:
        f.write(html_template)
    return None


# ── Top-level entry point ────────────────────────────────────────────────────


def make_report(
    queryrefaln,
    placement_hits,
    refstats,
    tree_html_path,
    preprocessing_failures_path,
    blast_failures_path,
    report_html,
    report_tsv,
):
    placement_stats = pd.read_csv(placement_hits, sep="\t")
    refdf = pd.read_csv(refstats, sep="\t").fillna("")
    preproc_failures = pd.read_csv(preprocessing_failures_path, sep="\t")
    blast_failures = pd.read_csv(blast_failures_path, sep="\t")

    placement_stats_sub = placement_stats[
        [
            "queryname",
            "placement_gt",
            "warnings",
            "blast_gt",
            "nn_gt",
            "nn_support",
            "querynode_support",
            "queryparent_support",
            "monophyletic",
            "mrca_support",
        ]
    ].fillna("")

    placement_stats_improved = make_per_gt_df(queryrefaln, placement_stats_sub, refdf)

    full_df = _build_full_report(
        placement_stats_improved, preproc_failures, blast_failures
    )

    refdf_subset = refdf[
        [
            "clade",
            "num_seqs_gt",
            "ham_pi",
            "max_ham_pi",
            "nogaps_ham_pi",
            "max_nogaps_ham_pi",
        ]
    ].rename(
        columns={
            "ham_pi": "gt_avg_ham_pi",
            "max_ham_pi": "max_gt_ham_pi",
            "nogaps_ham_pi": "gt_avg_nogaps_ham_pi",
            "max_nogaps_ham_pi": "max_gt_nogaps_ham_pi",
            "num_seqs_gt": "n_seqs",
        }
    )

    display_df = full_df.copy()
    if (
        display_df["stage_failed"].eq("").all()
        and display_df["failure_reason"].eq("").all()
    ):
        display_df = display_df.drop(columns=["stage_failed", "failure_reason"])

    generate_combined_interactive_html(
        display_df, refdf_subset, tree_html_path, report_html
    )

    full_df.to_csv(report_tsv, sep="\t", index=False)

    return None
