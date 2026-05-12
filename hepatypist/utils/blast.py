#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO


# splits multifasta query into individual fasta files with one seq each
def multifas_to_single(qfmt, outdir):
    recs = SeqIO.parse(qfmt, "fasta")
    for r in recs:
        h = r.description  # s = r.seq
        # outfile = os.path.join(outdir, h, ".fasta")
        outfile = f"{outdir}/{h}.fasta"
        # with open(outfile, "w") as o:
        #    o.write(f">{h}\n{s}\n")
        o = open(outfile, "w")
        SeqIO.write(r, o, "fasta")
        o.close()
    return None


def get_top_hit(df, pid, pc, lenmin):
    # applies filters to find top hit
    tophitrow = []
    df[["pident", "qcovs"]] = df[["pident", "qcovs"]].apply(
        pd.to_numeric, errors="coerce"
    )
    # filter df on pid and pcov thresholds
    filt_df = df[(df["pident"] >= pid) & (df["qcovs"] >= pc)]

    # check if qlen is really long - that will cause super low qcovs output

    if len(filt_df) > 0:
        sortdf = filt_df.copy()
        sortdf = sortdf.sort_values(
            ["pident", "qcovs", "sseqid"], ascending=[False, False, True]
        )
        top_pid = sortdf["pident"].iloc[0]
        top_cov = sortdf["qcovs"].iloc[0]
        top_len = sortdf["length"].iloc[0]
        seq = sortdf["qseq"].iloc[0]
        qseqid = sortdf["qseqid"].iloc[0]
        sseqid = sortdf["sseqid"].iloc[0]

        # get top PID match(es)
        top_pid_df = sortdf.loc[sortdf["pident"] == top_pid]
        # check if other top-PID matches share the same percent coverage
        top_hits = top_pid_df[top_pid_df["qcovs"] == top_cov]

        # check to make sure there's only one genotype in top hits
        top_hit_genotypes = list(
            set(list([el.split("|")[-1] for el in top_hits["sseqid"]]))
        )

        # if len(top_hit_genotypes)==1:
        #    top_match_genotype = top_pid_df["sseqid"].iloc[0].split("|")[-1]
        # else: #somewhat unlikely situation here with multiple genotypes returning same PID and PCOV
        #    top_match_genotype=";".join(top_hit_genotypes)

        th_genotype = top_hit_genotypes[0]
        report_genotypes = "-".join(top_hit_genotypes)

        if (top_pid >= pid) & (top_cov >= pc) & (top_len >= lenmin):
            tophitrow = [
                qseqid,
                sseqid,
                th_genotype,
                report_genotypes,
                top_pid,
                top_cov,
                top_len,
                seq,
            ]

    return tophitrow


# def parse_blast_out(blastlist, blastparsed, querymultifas, pid, pc, lenmin):
#     # parse each blast output tsv and get "top hit". Combine all top hit results into one "top hit" tsv
#     # create an output multifasta containing top hit match for each query yielding one
#     colnames=["qseqid","qlen","qseq","sseqid","slen","sstart","send","sseq","evalue","length","pident","qcovs","qcovshsp","qcovus"]#qseq,top_match_genotype,report_genotypes,top_pid,top_cov,top_len,seq
#     thcols = ["queryid","match_gt","report_gt","top_pid","top_cov","top_len","seq"]
#     rowlist = []

#     #blastlist = glob.glob(f"{blastoutdir}/*.tsv")

#     if len(blastlist)==0:
#         print("error: no blast output files found")

#     for f in blastlist:
#         if (os.path.getsize(f) > 0):#check if output isn't empty (if it has a filesize > 0)
#             df = pd.read_csv(f, sep="\t", names = colnames)
#             thr = get_top_hit(df, float(pid), float(pc), int(lenmin))#pd.read_csv(f, sep="\t",names=colnames), pid, pc, lenmin)
#             if len(thr)>0: #if there is a top hit
#                 rowlist.append(thr)#tophitrow = [qseqid,th_genotype,report_genotypes,top_pid,top_cov,top_len,seq]

#     if len(rowlist)>0:
#         thdf = pd.DataFrame(columns=thcols, data=rowlist)
#         thdf.to_csv(blastparsed,sep="\t",index=False)#write top hits to parsed tsv file

#         #write to fasta
#         reslist = [f">{h}\n{s}" for h, s in zip(thdf['queryid'], thdf['seq'])]
#         with open(querymultifas, "w+") as qf:
#             for res in reslist:
#                 qf.write(res)
#         qf.close()

#     else:#create empty file
#         bpf = open(blastparsed, 'w+')
#         bpf.close()
#         outf = open(querymultifas, 'w+')
#         outf.close()

#     return None


# def parse_blast_out(blast_out_path, outdir, pid, pc):
#     """
#     Parses tabular BLAST output to report top hit/s, if any.
#     A top hit must be at least 100bp in length and have the highest percent identity (PID) value of all hits.
#     If multiple hits share the highest PID and PC, then all will be reported.
#     """
#     #where to define these?
#     colnames = ["qseqid","qlen","qseq","sseqid","slen","sstart","send","sseq","evalue","length","pident","qcovs","qcovshsp","qcovus"]

#     #check if output is empty (if it has a filesize of 0)
#     file_size = os.path.getsize(blast_out_path)
#     if (file_size == 0):
#         return ValueError("File is empty")

#     df = pd.read_csv(blast_out_path, sep="\t", names=colnames)
#     df[["pident", "qcovs"]] = df[["pident", "qcovs"]].apply(pd.to_numeric, errors='coerce')

#     #filter df on pid and pcov thresholds
#     filt_df = df[(df["pident"]>=pid) & (df["qcovs"]>=pc)]

#     # return error if there are no BLAST hits over specified thresholds
#     if len(filt_df)==0:
#         return ValueError(f"No matches above specified tresholds for pct ID ({pid}) and pct coverage ({pc})")

#     sortdf = filt_df.copy()
#     sortdf = sortdf.sort_values(['pident','qcovs', "sseqid"], ascending=[False,False,True])

#     top_pid = sortdf["pident"].iloc[0]
#     top_cov = sortdf["qcovs"].iloc[0]
#     print(f"top percent identity {top_pid}")
#     print(f"top percent coverage {top_cov}")

#     # get top PID match(es)
#     top_pid_df = sortdf.loc[sortdf["pident"] == top_pid]
#     # check if other top-PID matches share the same percent coverage
#     top_hits = top_pid_df[top_pid_df["qcovs"] == top_cov]

#     # check to make sure there's only one genotype in top hits
#     #top_hit_genotypes = list(set(list([el.split("|")[-1] for el in top_hits["sseqid"]])))
#     top_hit_genotypes = list(set(list([el.split("|")[-1] for el in top_hits["sseqid"]])))

#     if len(top_hit_genotypes)==1:
#         top_match_ID = top_pid_df["sseqid"].iloc[0]
#         top_match_genotype = top_match_ID.split("|")[-1]
#         print(f"Top hit ID: {top_match_ID}")
#         print(f"Top hit genotype: {top_match_genotype}")

#     else: #somewhat unlikely situation here with multiple genotypes returning same PID and PCOV
#         print("Multiple genotypes share same top percent identity and percent coverage")
#         print(f"{len(top_hit_genotypes)} genotypes share the same top percent identity and percent coverage:")
#         print(top_hit_genotypes)

#     th_genotype = top_hit_genotypes[0]
#     report_genotypes = "-".join(top_hit_genotypes)

#     print(f"blast genotype: {report_genotypes}")

#     #print top hit query sequence to output fasta
#     queryout_base = top_hits["qseqid"].iloc[0] + "_VP1-P2A_region.fasta"
#     query_header = ">" + top_hits["qseqid"].iloc[0]# + "|" + report_genotype
#     query_seq = top_hits["qseq"].iloc[0]

#     outfasta_path = os.path.join(outdir, queryout_base)

#     # #maybe should check if outfasta already exists - return error rather than append/overwrite?
#     # if os.path.exists(outfasta_path):
#     #     return ValueError("output fasta path already exists - delete to overwrite")

#     with open(outfasta_path, "w") as outf:
#         outf.seek(0) #to overwrite existing file
#         outf.write(query_header + "\n")
#         outf.write(query_seq + "\n")

#     #double-check
#     outf.close()

#     return outfasta_path, query_header, th_genotype
