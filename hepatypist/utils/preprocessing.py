#!/usr/bin/env python
import re
import itertools
import pandas as pd
import numpy as np
from Bio import SeqIO


# def read_fasta(fasta_filepath):
#     """
#     Reads fasta using BioPython's SeqIO module, which has built-in format checks
#     Splits headers on ", " to return just accession
#     Returns sequences as dictionary, in format dict[header] = seq
#     """
#     #fdict = SeqIO.to_dict(SeqIO.parse(fasta_filepath, "fasta"))
#     biod = SeqIO.to_dict(SeqIO.parse(fasta_filepath, "fasta"))
#     fdict = {h:strip_seq(s) for h,s in biod.items()}
#     #fixdict = {key.split(",")[0]:str(value.seq) for key, value in fdict.items()}
#     return fdict


def split_string_delim(s, delim, maxlen):
    # splits string at a given delimiter...
    # but if delimiter occurs multiple times... then splits at the occurrence of delimiter that gives the longest possible length without exceeding length threshold
    if delim not in s:
        return ValueError(f"delimiter {delim} not in string {s}")

    if len(s) <= maxlen:
        return s

    else:
        delim_ind = [i for i, letter in enumerate(s) if letter == delim and i <= maxlen]
        maxind = max(delim_ind)
        return s[:maxind]


def edit_fasta_headers(QorR, h, maxlen):
    if QorR == "query":
        delimchar = "|" if "|" in h else " "
        newh = (
            h.split(delimchar)[0]
            .rstrip(" ")
            .replace(" ", "_")
            .replace(">", "")
            .replace("\n", "")
            .rstrip("_")[:maxlen]
        )
        return newh
    else:
        fh = h.replace("\n", "").replace(">", "").replace(" ", "_")
        if len(fh) < maxlen:
            return fh
        else:
            hn, gt = fh.split("|")  # [-1]
            return hn[: (maxlen - len("|" + gt))].rstrip("_") + "|" + gt


def strip_seq(s):
    sedit = (
        str(s).replace("\n", "").replace("\r", "")
    )  # sedit = re.sub(r"\n", "", str(s))
    return sedit


def check_valid_headerchars_and_format(header):
    # checks to make sure headers contain only alphanumeric + allowed specials
    # raise error if any header contains disallowed characters; checks that all headers are unique.
    # ^ - Matches the start of the string
    # [a-zA-Z0-9._\-|] - Matches any of the allowed characters:
    #                   a-z (lowercase letters), A-Z (uppercase letters), 0-9 (digits),
    #                   . (dot), _ (underscore), - (dash), | (pipe)
    #                   + - Ensures one or more allowed characters are present
    #                   $ - Matches the end of the string
    pattern = r"^[a-zA-Z0-9._\-|]+$"
    return bool(re.match(pattern, header))


def check_valid_seqchars(seq):
    # Checks to see whether there are any invalid/not permitted chars in the sequence.
    # I arbitrarily threw in a limit of 20% ambiguity
    # I initially only allowed standard nucleotide (ATCG), gaps (-), and N characters, but expanded to anything FastTree tolerates.
    # Apparently FastTree doesn't like ambiguity, and issues warnings for R/N/etc and treats them as missing data, but will still run
    pattern = r"^[ACTGNRYKMBDHVSWactgnrykmbdhvsw-]+$"  # pattern = r'^[ACTGNactgn-]+$'
    return bool(re.match(pattern, str(seq)))


def check_ambig_pct(seq):
    ambigs = list("NRYKMBDHVWnrykmbdhvsw")
    seqambigs = [c for c in seq if c in ambigs]  # if len(seqambigs)>0:
    pctambig = float(len(seqambigs) / len(seq))
    return pctambig < 0.20


def check_seqlength(seq, lengthmin):
    if len(seq) >= lengthmin:
        return True
    else:
        return False


def _failure_reason(fmtdict, lengthmin):
    reasons = []
    if not fmtdict.get("replace_error", True):
        reasons.append("header could not be reformatted")
    if fmtdict.get("disallowed_header") == "Y":
        reasons.append("header contains disallowed characters")
    if fmtdict.get("disallowed_nt") == "Y":
        reasons.append("sequence contains disallowed nucleotide characters")
    if fmtdict.get("ambig_over_20pct") == "Y":
        reasons.append(">20% ambiguous nucleotides")
    if fmtdict.get("too_short") == "Y":
        reasons.append(f"sequence too short (minimum {lengthmin} bp)")
    return "; ".join(reasons)


def fasta_check_pythonic(
    infasta, QorR, lengthmin, fmtout, qcoutfile, dedup=None, failoutfile=None
):
    # new and improved. instead of for loop, uses more pythonic list comprhension to identify sequences that pass filter
    if QorR not in ["query", "ref"]:  # QorF must be either ref or query
        return ValueError("must specify query or ref")

    if dedup not in [
        None,
        True,
        False,
        "True",
        "true",
        "TRUE",
        "T",
        "False",
        "false",
        "FALSE",
        "F",
    ]:  # dedup must be unspecified, or string value of "True" or "False"
        return ValueError("dedup must be either unspecified or 'True'/'False'")

    try:
        fascheck = SeqIO.parse(infasta, "fasta")
    except Exception:  # FileExistsError, YouAreBeingMeanException) as e:
        return ValueError("Input fasta is not in valid format")

    fmtfasd = {}
    qcdict = {}
    maxheaderlen = 50

    for r in fascheck:
        fmtdict = {}  # header, header_fmt, disallowed_header, disallowed_nt, ambig_over_20_pct, too_short, replace_error
        h = r.description
        fmtdict["header"] = h

        strseq = strip_seq(r.seq)

        try:
            newh = edit_fasta_headers(QorR, h, maxheaderlen)
            fmtpass = True
        except Exception:  #:
            newh = h
            fmtpass = False

        fmtdict["header_fmt"] = newh
        fmtdict["replace_error"] = fmtpass

        headercheck = check_valid_headerchars_and_format(newh)
        fmtdict["disallowed_header"] = "N" if headercheck else "Y"

        seqcheck = check_valid_seqchars(strseq)
        fmtdict["disallowed_nt"] = "N" if seqcheck else "Y"

        ambigcheck = check_ambig_pct(strseq)
        fmtdict["ambig_over_20pct"] = "N" if ambigcheck else "Y"

        minlencheck = check_seqlength(strseq, lengthmin)
        fmtdict["too_short"] = "N" if minlencheck else "Y"

        qcdict[newh] = fmtdict

        # if "Y" not in fmtdict.values():#[fmtdict[fk] for fk in FAILLIST]:#if "Y" not in list(fmtdict.values())[2:5]:
        if all(
            [fmtpass, headercheck, seqcheck, ambigcheck, minlencheck]
        ):  # all(not x for x in [fmterr,headercheck,seqcheck,ambigcheck,minlencheck]):
            fmtfasd[newh] = strseq

    # deduplication of ref optional here can add back in as optional argument if desired
    if dedup in [True, "True", "TRUE", "true", "T"]:
        writed = dedup_ref(fmtfasd.copy())
    else:
        writed = fmtfasd.copy()

    ofi = open(fmtout, "w")
    for h, s in writed.items():
        ofi.write(f">{h}\n{s}\n")  # SeqIO.write(rf,of,'fasta')
    ofi.close()

    qdf = pd.DataFrame(qcdict).T
    qdf.to_csv(qcoutfile, sep="\t", index=False)

    if failoutfile is not None and QorR == "query":
        fail_rows = []
        for h, fmtd in qcdict.items():
            passed = all(
                [
                    fmtd.get("replace_error", True),
                    fmtd.get("disallowed_header") == "N",
                    fmtd.get("disallowed_nt") == "N",
                    fmtd.get("ambig_over_20pct") == "N",
                    fmtd.get("too_short") == "N",
                ]
            )
            if not passed:
                fail_rows.append(
                    {
                        "queryname": h,
                        "stage": "preprocessing",
                        "reason": _failure_reason(fmtd, lengthmin),
                    }
                )
        pd.DataFrame(fail_rows, columns=["queryname", "stage", "reason"]).to_csv(
            failoutfile, sep="\t", index=False
        )

    return None


def hamming_distance(s1, s2) -> int:
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def no_gap_hamming_distance(s1, s2) -> int:
    return sum(el1 != el2 for el1, el2 in zip(s1, s2) if "-" not in (el1, el2))


def diversity(fdict):
    # first check to make sure fasta is aligned (i.e., all seqs are same length)
    sumlen = sum([len(v) for v in fdict.values()])
    avlen = sumlen / len(fdict)
    # Boolean test if all strings are same length as average length
    res = all(len(x) == avlen for x in fdict.values())
    if not res:
        return ValueError("input fasta is not aligned")

    # iterate through dict and count all nonredundant, non-gap pairwise comparisons
    seq_pairs = list(itertools.product(fdict.keys(), fdict.keys()))
    nonred_pairs = [(a, b) for (a, b) in seq_pairs if a != b]

    nmat = len(fdict.keys())
    empmatr = np.empty((nmat, nmat))
    empmatr[:] = np.nan
    empmatr2 = np.empty((nmat, nmat))
    empmatr2[:] = np.nan

    keylist = sorted(list(fdict.keys()))
    keyindd = {}
    for i, k in enumerate(keylist):
        keyindd[k] = i

    named_matrix = {
        "row_names": keylist,
        "col_names": keylist,
        "ham": empmatr,
        "ngham": empmatr2,
    }

    hamlist = []
    ng_hamlist = []
    hammedoid = ""
    ngmedoid = ""

    for a, b in nonred_pairs:
        abham = hamming_distance(fdict[a], fdict[b])
        abngham = no_gap_hamming_distance(fdict[a], fdict[b])
        hamlist.append(abham)
        ng_hamlist.append(abngham)

        a_index = keyindd[a]
        b_index = keyindd[b]
        named_matrix["ham"][a_index, b_index] = abham
        named_matrix["ham"][b_index, a_index] = abham
        named_matrix["ngham"][a_index, b_index] = abngham
        named_matrix["ngham"][b_index, a_index] = abngham

    hamaveraged = {}
    nghamaveraged = {}
    for k in keylist:  # get all non-Nan vals in that key's row and find average
        k_index = named_matrix["row_names"].index(k)
        hamvals = named_matrix["ham"][k_index, :]
        non_nan_hams = hamvals[~np.isnan(hamvals)]
        nghamvals = named_matrix["ngham"][k_index, :]
        non_nan_nghams = nghamvals[~np.isnan(nghamvals)]
        kavgham = sum(non_nan_hams) / len(non_nan_hams)
        kavgngham = sum(non_nan_nghams) / len(non_nan_nghams)
        hamaveraged[k] = kavgham
        nghamaveraged[k] = kavgngham

    minavham = min(list(hamaveraged.values()))
    minngavham = min(list(nghamaveraged.values()))

    max_ham = max(list(hamaveraged.values()))
    max_ng_ham = max(list(nghamaveraged.values()))

    minhamkeys = sorted(list(set([k for k, v in hamaveraged.items() if v == minavham])))
    minngkeys = sorted(
        list(set([k for k, v in nghamaveraged.items() if v == minngavham]))
    )

    hammedoid = minhamkeys[0]
    ngmedoid = minngkeys[0]

    if len(hamlist) == 0:
        print(hamlist)
        return ValueError("Hamming dist list len 0")

    else:
        avg_ham = sum(hamlist) / len(hamlist)
        avg_ng_ham = sum(ng_hamlist) / len(ng_hamlist)
        ham_pi = "{:.3f}".format(avg_ham / avlen)
        ng_ham_pi = "{:.3f}".format(avg_ng_ham / avlen)
        max_ham_pi = "{:.3f}".format(max_ham / avlen)
        max_ng_ham_pi = "{:.3f}".format(max_ng_ham / avlen)
        avg_ham_100 = "{:.3f}".format(avg_ham / (avlen / 100))
        avg_ng_ham_100 = "{:.3f}".format(avg_ng_ham / (avlen / 100))

    return (
        avlen,
        ham_pi,
        ng_ham_pi,
        max_ham_pi,
        max_ng_ham_pi,
        avg_ham_100,
        avg_ng_ham_100,
        hammedoid,
        ngmedoid,
    )


def num_unique(fdict):
    # gets number of unique seqs in set of values
    return len(set(fdict.values()))


def gt_metrics(fd):
    gtd = {}
    genotypes = set([head.split("|")[-1] for head in fd.keys()])

    for g in genotypes:
        subdict = {k: v for k, v in fd.items() if k.endswith(f"|{g}")}
        if len(subdict) == 0:
            (
                avlen,
                ham_pi,
                ng_ham_pi,
                max_ham_pi,
                max_ng_ham_pi,
                avg_ham_100,
                avg_ng_ham_100,
                hammedoid,
                ngmedoid,
            ) = ["NA"] * 9
        elif len(subdict) == 1:
            single_key = list(subdict.keys())[0]
            avlen = len(list(subdict.values())[0])
            ham_pi = ng_ham_pi = max_ham_pi = max_ng_ham_pi = avg_ham_100 = (
                avg_ng_ham_100
            ) = "NA"
            hammedoid = single_key
            ngmedoid = single_key
        else:
            (
                avlen,
                ham_pi,
                ng_ham_pi,
                max_ham_pi,
                max_ng_ham_pi,
                avg_ham_100,
                avg_ng_ham_100,
                hammedoid,
                ngmedoid,
            ) = diversity(subdict)

        gtd[g] = [
            ham_pi,
            ng_ham_pi,
            max_ham_pi,
            max_ng_ham_pi,
            len(subdict),
            avlen,
            hammedoid,
            ngmedoid,
        ]

    return gtd


def non_zero_min(lt):
    nonzeroes = [x for x in lt if x != 0]
    return min(nonzeroes)


def calc_ref_stats(reffas, outfile):
    # calculates overall and per-clade seq and diversity stats, output in tsv format
    # align  = SeqIO.parse(reffas, "fasta")
    fd = {r.description: strip_seq(r.seq) for r in SeqIO.parse(reffas, "fasta")}
    gtd = gt_metrics(fd)

    gtddf = pd.DataFrame(
        gtd
    ).T.reset_index()  # avlen, ham_pi, ng_ham_pi, avg_ham_100, avg_ng_ham_100, avlen, hammedoid, ngmedoid
    # ham_pi, ng_ham_pi, max_ham_pi, max_ng_ham_pi, len(subdict), avlen, hammedoid, ngmedoid
    gtddf.columns = [
        "clade",
        "ham_pi",
        "nogaps_ham_pi",
        "max_ham_pi",
        "max_nogaps_ham_pi",
        "num_seqs_gt",
        "avlen",
        "gt_ham_medoid",
        "gt_ng_ham_medoid",
    ]

    sortgtddf = gtddf.sort_values(by="clade")

    sortgtddf.to_csv(outfile, sep="\t", index=False)

    return None


def dedup_ref(rd):  # , dedup_reffasta):
    # rd = read_fasta(reffasta)
    uniqueseqs = list(set(rd.values()))
    uniqueseqdict = {}
    for u in uniqueseqs:
        matches = {h: s for h, s in rd.items() if s == u}
        keylist = [h.split("|")[0] for h, s in matches.items()]
        gtset = list(set([h.split("|")[-1] for h, s in matches.items()]))

        if len(gtset) == 1:
            gt = gtset[0]
        elif len(gtset) > 1:
            gt = "_".join(gtset)
            print(f"{gt}: multiple gts for unique seq {u}")
        else:
            return ValueError(f"no gt found for unique seq {u}")

        newkey = "___".join(keylist) + "|" + gt
        uniqueseqdict[newkey] = u

    return uniqueseqdict

    # def dedup_ref(rd):#, dedup_reffasta):
    #     #rd = read_fasta(reffasta)
    #     uniqueseqs = list(set(rd.values()))
    #     uniqueseqdict = {}
    #     for u in uniqueseqs:
    #         keylist = []
    #         for h,s in rd.items():
    #             if s==u:
    #                 keylist.append(h.split("|")[0])
    #         #keylist = [h.split("|")[0] for h,s in rd.items() if s==u]
    #         gtlist = list(set([h.split("|")[1] for h,s in rd.items() if s==u]))
    #         if len(gtlist)==1:
    #             gt = gtlist
    #         else:
    #             gt = "_".join(gtlist)
    #         newkey = "___".join(keylist) + "|" + gt
    #         uniqueseqdict[newkey] = u

    return uniqueseqdict

    # with open(dedup_reffasta, "w") as f:
    #     for s,h in uniqueseqdict.items():
    #         f.write(f">{h}\n{s}\n")

    # return None
