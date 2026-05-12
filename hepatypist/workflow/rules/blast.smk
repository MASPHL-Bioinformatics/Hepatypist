#!/usr/bin/env python
# pipeline = "hepatypist"
"""
Creates blast database from reference fasta.
Runs BLAST for each query sequence (since initial QC was passed)
Outputs for each query:
    blast result text file
    match fasta

Produces combined overall tsv with "top hit" and genotypes

Created 24 Oct 2024
@author: Mary Godec
@contact: mary.godec@mass.gov
"""
from hepatypist.utils.blast import *
from Bio import SeqIO
import glob
import os
import pandas as pd


def fasta_to_dict(infasta):
    fasta = open(infasta)
    sequences = fasta.read()
    sequences = re.split(
        "^>", sequences, flags=re.MULTILINE
    )  # Only splits string at the start of a line.
    fasta.close()
    del sequences[0]  # the first in this list will be blank, so delete it
    fastadict = {}
    for fasta in sequences:
        header, sequence = fasta.split(
            "\n", 1
        )  # Split each fasta into header and sequence.
        header = (
            ">" + header + "\n"
        )  # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
        sequence = sequence.replace("\n", "")  # Take line breaks out of sequence.
        fastadict[header] = sequence
    return fastadict


def edit_fasta_headers(h):
    delim = "|" if "|" in h else " "
    newh = (
        h.split(delim)[0]
        .rstrip(" ")
        .replace(" ", "_")
        .replace(">", "")
        .replace("\n", "")
        .rstrip("_")
    )
    trunch = newh[
        :50
    ]  # complete filepath including folder layers should be <256, which is the limit of most systems - also cutting down query helps for tree display
    return trunch


SAMNAMES = [edit_fasta_headers(h) for h, s in fasta_to_dict(config["query"]).items()]

THCOLS = [
    "queryid",
    "matchid",
    "match_gt",
    "report_gt",
    "top_pid",
    "top_cov",
    "top_len",
    "seq",
]
BLAST_COLNAMES = [
    "qseqid",
    "qlen",
    "qseq",
    "sseqid",
    "slen",
    "sstart",
    "send",
    "sseq",
    "evalue",
    "length",
    "pident",
    "qcovs",
    "qcovshsp",
    "qcovus",
]


# PREPROCESSING output
queryfasta_copy = RESULTS_DIR + "qc/query.fasta"
queryfasta_qc = RESULTS_DIR + "qc/query_qc.txt"
query_fmt = RESULTS_DIR + "qc/query_fmt.fasta"
reffasta_copy = RESULTS_DIR + "qc/ref.fasta"  # conditional
reffasta_qc = RESULTS_DIR + "qc/ref_qc.txt"  # conditional
ref_stats = RESULTS_DIR + "qc/refstats.tsv"  # conditional

if config["ref_aligned"] != "":
    reffasta_aligned = config["ref_aligned"]
    ref_stats = config["ref_stats"]
    ref_formatted = config["ref"]
else:
    reffasta_aligned = RESULTS_DIR + "qc/ref_aligned.fasta"  # conditional
    reffasta_qc = RESULTS_DIR + "qc/ref_qc.txt"
    ref_formatted = RESULTS_DIR + "qc/ref_fmt.fasta"

# BLAST output
if config["ref_blastdb"] != "":
    blastdbpath = config["ref_blastdb"]
else:
    blastdbpath = RESULTS_DIR + "blast/blastdb/refdb"

blastdbndb = blastdbpath + ".ndb"
blastfasdir = os.path.join(RESULTS_DIR, "blast/blastinfiles/")
blastoutdir = os.path.join(RESULTS_DIR, "blast/blastoutfiles/")
blastparsed = os.path.join(RESULTS_DIR, "blast/blast_parsed.txt")
querymultifas = os.path.join(RESULTS_DIR, "blast/query_matches.fasta")
blast_failures = os.path.join(RESULTS_DIR, "blast/blast_failures.tsv")


wildcard_constraints:
    sn="|".join(SAMNAMES),


def agg_splitfas(wildcards):
    FOUTD = checkpoints.get_splitfas.get(**wildcards).output[0]
    (SN,) = glob_wildcards(FOUTD + "/{sn}.fasta")
    return expand(RESULTS_DIR + "blast/blastinfiles/{sn}.fasta", sn=SN)


def agg_blastout(wildcards):
    BOUTD = checkpoints.run_blast.get(**wildcards).output[0]
    (SN,) = glob_wildcards(BOUTD + "/{sn}.tsv")
    return expand(RESULTS_DIR + "blast/blastoutfiles/{sn}.tsv", sn=SN)


rule make_blastdb:
    input:
        ref_formatted,
    output:
        blastdbndb,
    params:
        bdir=blastdbpath,
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -parse_seqids -out {params.bdir}
        """


# https://stackoverflow.com/questions/56241962/snakemake-getting-checkpoint-and-aggregate-function-to-work
checkpoint get_splitfas:
    input:
        qf=query_fmt,
    output:
        directory(blastfasdir),
    run:
        os.makedirs(blastfasdir, exist_ok=True)
        multifas_to_single(input.qf, blastfasdir)


checkpoint run_blast:
    input:
        agg_splitfas,
    output:
        directory(blastoutdir),
    params:
        nthreads=config["nthreads"],
        bd=blastdbpath,
    run:
        os.makedirs(blastoutdir, exist_ok=True)
        for f in input:
            OF = blastoutdir + str(os.path.basename(f)).replace(".fasta", ".tsv")
            shell(
                "blastn -query {f} -out {OF} -db {params.bd} -num_threads {params.nthreads} -evalue 0.05 -qcov_hsp_perc 50 -outfmt '6 qseqid qlen qseq sseqid slen sstart send sseq evalue length pident qcovs qcovshsp qcovus'"
            )


rule parse_blast:
    input:
        agg_blastout,
    output:
        bp=blastparsed,
        qmf=querymultifas,
        bf=blast_failures,
    params:
        pid=config["pid_threshold"],
        pc=config["pc_threshold"],
        lenmin=config["blast_len_threshold"],
    run:
        rowlist = []
        fail_rows = []
        for f in input:
            sn = os.path.basename(f).replace(".tsv", "")
            if os.path.getsize(f) > 0:
                temp = pd.read_csv(f, sep="\t", names=BLAST_COLNAMES)
                if len(temp) > 0:
                    thr = get_top_hit(
                        temp, float(params.pid), float(params.pc), int(params.lenmin)
                    )
                    if len(thr) > 0:
                        rowlist.append(thr)
                    else:
                        fail_rows.append(
                            {
                                "queryname": sn,
                                "stage": "blast",
                                "reason": (
                                    f"no BLAST hits met thresholds "
                                    f"(pid≥{params.pid}, coverage≥{params.pc}, "
                                    f"len≥{params.lenmin})"
                                ),
                            }
                        )
                else:
                    fail_rows.append(
                        {
                            "queryname": sn,
                            "stage": "blast",
                            "reason": "no BLAST hits found in reference database",
                        }
                    )
            else:
                fail_rows.append(
                    {
                        "queryname": sn,
                        "stage": "blast",
                        "reason": "BLAST output file was empty",
                    }
                )

        pd.DataFrame(fail_rows, columns=["queryname", "stage", "reason"]).to_csv(
            output.bf, sep="\t", index=False
        )

        if len(rowlist) > 0:
            thdf = pd.DataFrame(columns=THCOLS, data=rowlist)
            thdf.to_csv(output.bp, sep="\t", index=False)

            reslist = [(h, s) for h, s in zip(thdf["queryid"], thdf["seq"])]
            with open(output.qmf, "w+") as qf:
                for h, s in reslist:
                    qf.write(">" + h + "\n")
                    qf.write(s + "\n")
        else:
            bpf = open(output.bp, "w+")
            bpf.close()
            outf = open(output.qmf, "w+")
            outf.close()
