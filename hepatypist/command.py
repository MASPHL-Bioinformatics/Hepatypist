#!/usr/bin/env python

from hepatypist import _program
from hepatypist import __version__
import hepatypist as _havtyper_pkg

import sys
import os
import argparse
import platform
import datetime as dt

import pandas as pd
from loguru import logger
from Bio import SeqIO
import snakemake


thisdir = os.path.abspath(os.path.dirname(os.path.realpath(sys.argv[0])))
cwd = os.path.abspath(os.getcwd())

# Locate package resources relative to the installed package directory.
# Using __file__ works for both regular and editable installs — unlike
# site.getsitepackages() + glob, which breaks for editable/path installs.
_pkg_dir = os.path.dirname(os.path.abspath(_havtyper_pkg.__file__))

defaultrefpath = os.path.join(_pkg_dir, "data", "reference_set.fasta")
defaultrefalnpath = os.path.join(_pkg_dir, "data", "reference_set_aligned.fasta")
defaultrefblastdb = os.path.join(_pkg_dir, "data", "blastdb", "refdb")
defaultrefstats = os.path.join(_pkg_dir, "data", "refstats.tsv")
snakefilepath = os.path.join(_pkg_dir, "workflow", "Snakefile")

libdir = os.path.dirname(snakefilepath)


# determines the min length threshold for BLAST top hit reporting
# user can supply a blast_len_threshold
# but if they didn't, calculates default threshold: 80% of the mean length of sequences in the ref fasta
def get_minlen(configval, reffasta):
    if configval > 0:
        return configval
    else:
        lenlist = []
        for r in SeqIO.parse(reffasta, "fasta"):
            lenlist.append(len(r.seq.replace("-", "")))

        meanlen = int(sum(lenlist) / len(lenlist))
        return int(0.80 * meanlen)


def get_maxqueries(argmax, reffasta):
    if argmax > 0:
        return argmax
    else:
        num_refseqs = len([r for r in SeqIO.parse(reffasta, "fasta")])
        return int(0.2 * num_refseqs)


def _tprint(msg=""):
    """Print a line to the real terminal regardless of sys.stdout redirects."""
    print(msg, flush=True)


def main(sysargs=sys.argv[1:]):
    logger.remove()  # remove default handler/logger configuration, may also require logger.remove(0)
    logname = "log-{date:%Y-%m-%d_%H-%M-%S}.log".format(
        date=dt.datetime.now()
    )  # capture start time in log file ame
    log_level = "DEBUG"
    log_format = "<green>{time:YYYY-MM-DD HH:mm:ss.SSS zz}</green> | <level>{level: <8}</level> | <yellow>Line {line: >4} ({file}):</yellow> <b>{message}</b>"

    parser = argparse.ArgumentParser(
        prog=_program,
        description="hepatypist: Subgenotype assignment of VP2-1A region",
        usage="""hepatypist <query> [options]""",
    )

    io_group = parser.add_argument_group("Input-Output options")
    io_group.add_argument(
        "--query",
        "-q",
        required=True,
        action="store",
        help="Query fasta file of sequences to analyze, maximum # can be specified with --max_query_num. Length of all query sequences must meet or exceed to the length of the shortest sequence in the reference set",
    )
    io_group.add_argument(
        "--ref",
        "-r",
        required=False,
        default=defaultrefpath,
        action="store",
        help="OPTIONAL: Reference fasta file with headers in format {> ID_NAME | CLADENUM}. If user doesn't supply own reference, Hepatypist will default its built-in reference fasta",
    )
    io_group.add_argument(
        "--refdedup",
        "-rd",
        required=False,
        default=False,
        action="store",
        help="OPTIONAL: If user specifies own reference, the --ref_dedup flag will deduplicate that reference before analysis",
    )
    io_group.add_argument(
        "--outdir",
        "-o",
        required=False,
        default=cwd,
        action="store",
        help="Output directory. Default: current working directory",
    )
    io_group.add_argument(
        "--max_query_num",
        "-maxq",
        required=False,
        default=0,
        type=int,
        help="Maximum number of query seqs allowed in query multifasta",
    )

    a_group = parser.add_argument_group("Analysis options")
    a_group.add_argument(
        "--nthreads",
        "-n",
        required=False,
        default=1,
        type=int,
        help="OPTIONAL: number of threads to run Snakemake. Default = 1",
    )
    a_group.add_argument(
        "--percent_id",
        "-pid",
        required=False,
        default=0.90,
        type=float,
        help="Optional alternate percent identity threshold for blast hit reporting. Default is 90 pct",
    )
    a_group.add_argument(
        "--percent_coverage",
        "-pc",
        required=False,
        default=0.50,
        type=float,
        help="Optional alternate percent coverage (qcovs) threshold for blast hit reporting. Default is 50 pct",
    )
    a_group.add_argument(
        "--blast_minlen",
        "-bl",
        required=False,
        default=0,
        type=int,
        help="Optional alternate minimum length for blast match top hit. Default value is 80pct of mean seq length in reference set.",
    )
    a_group.add_argument(
        "--physupport_thresh",
        "--pst",
        required=False,
        default=70,
        type=int,
        help='Optional alternate "acceptable" bootstrap/SH value for node support',
    )
    a_group.add_argument(
        "--lowsupport_thresh",
        "--lst",
        required=False,
        default=50,
        type=int,
        help='Minimum "acceptable" bootstrap/SH value, used for coloring support values',
    )

    # check for required arguments
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    logfile = os.path.join(str(args.outdir), logname)
    logger.add(
        logfile,
        level=log_level,
        format=log_format,
        colorize=False,
        backtrace=True,
        diagnose=True,
    )

    logger.info(f"Python {sys.version}")  # system and platform info
    platform_uname = list(platform.uname())
    logger.info(f"OS: {platform_uname[0]}")
    logger.info(f"node: {platform_uname[1]}")
    logger.info(f"release: {platform_uname[2]}")
    logger.info(f"version: {platform_uname[3]}")
    logger.info(f"machine: {platform_uname[4]}")
    for k, v in platform.freedesktop_os_release().items():
        logger.info(f"Distro {k}: {v}")

    if args.refdedup is True and args.ref == defaultrefpath:
        parser.error(
            "--refdedup/-rd flag should only be used if using non-default reference fasta"
        )

    ## do basic input file existence checks
    if not os.path.isfile(os.path.abspath(args.query)):
        print(f"couldn't find query {args.query}")
        sys.exit(0)

    if not os.path.isfile(os.path.abspath(args.ref)):
        print(f"couldn't find ref {args.ref}")
        sys.exit(0)

    query_count = sum(1 for _ in SeqIO.parse(os.path.abspath(args.query), "fasta"))
    if query_count > 25:
        print(
            f"Error: query fasta contains {query_count} sequences. "
            "Maximum allowed is 25. Please split into smaller batches."
        )
        sys.exit(1)

    if not os.path.isfile(snakefilepath):
        print(f"couldn't find snakefile {snakefilepath}")
        sys.exit(0)

    unformatted_outdir = str(os.path.abspath(args.outdir))
    formatted_outdir = (
        unformatted_outdir
        if unformatted_outdir.endswith("/")
        else unformatted_outdir + "/"
    )

    # Initialise config dict
    config = {}
    config["query"] = str(os.path.abspath(args.query))
    config["ref"] = str(os.path.abspath(args.ref))

    if args.ref == defaultrefpath:
        config["ref_aligned"] = defaultrefalnpath
        config["ref_blastdb"] = defaultrefblastdb
        config["ref_stats"] = defaultrefstats

    else:
        config["ref_aligned"] = ""
        config["ref_blastdb"] = ""
        config["ref_stats"] = ""

    config["ref_dedup"] = args.refdedup
    config["query_minlen"] = 100
    config["ref_minlen"] = 200
    config["pid_threshold"] = args.percent_id
    config["pc_threshold"] = args.percent_coverage
    config["blast_len_threshold"] = get_minlen(args.blast_minlen, args.ref)
    config["high_support_threshold"] = args.physupport_thresh
    config["low_support_threshold"] = args.lowsupport_thresh
    config["maxqueryseqs"] = get_maxqueries(args.max_query_num, args.ref)
    config["nthreads"] = str(args.nthreads)
    config["resultsdir"] = formatted_outdir
    config["configoutfile"] = f"{config['resultsdir']}config/config.yaml"
    config["logsdir"] = f"{config['resultsdir']}logs/"
    config["rulesdir"] = f"{libdir}/rules/"
    config["runinfo_version"] = __version__
    config["runinfo_rundir"] = cwd
    config["runinfo_platform"] = (
        f"{platform.platform()}"  # single user readable string with as much useful info as possible
    )
    # config["runinfo_osname"]=f"{os.name}"
    # config["runinfo_sysrelease"]=f"{platform.system()} {platform.release()}" # raw os name as supplied by OS itself
    # config["runinfo_dist"]=f"{" ".join(platform.dist())}"
    # config["runinfo_distversion"]=f"{platform.version()}"
    # config["runinfo_machine"]=f"{platform.machine()}"
    # config["runinfo_architecture"]=f"{platform.architecture()}"

    #### POTENTIAL TO-DO: SET CACHING HERE (or poss in CICD)

    # ── Terminal: start banner ────────────────────────────────────────────────
    _tprint(f"hepatypist {__version__}")
    _tprint(
        f"Query:    {os.path.abspath(args.query)}  ({query_count} sequence{'s' if query_count != 1 else ''})"
    )
    _tprint(f"Output:   {formatted_outdir}")
    _tprint(f"Threads:  {args.nthreads}")
    _tprint(f"Log:      {logfile}")
    _tprint("Running pipeline — full output in log file...")

    # ── Run Snakemake, routing ALL its output to the log file ─────────────────
    # We redirect at the file-descriptor level (os.dup2) so that even output
    # written directly to fd 1/2 by child processes goes to the log, not the
    # terminal.  Python-level sys.stdout/stderr are redirected too for safety.
    _saved_py_out = sys.stdout
    _saved_py_err = sys.stderr
    _saved_fd_out = os.dup(1)
    _saved_fd_err = os.dup(2)

    status = False
    try:
        with open(logfile, "a") as _lf:
            sys.stdout = _lf
            sys.stderr = _lf
            os.dup2(_lf.fileno(), 1)
            os.dup2(_lf.fileno(), 2)
            # https://snakemake.readthedocs.io/en/v7.0.0/api_reference/snakemake.html
            status = snakemake.snakemake(
                snakefilepath,
                printshellcmds=False,
                forceall=True,
                force_incomplete=True,
                config=config,
                cores=int(args.nthreads),
                lock=False,
            )
    finally:
        # Always restore terminal — even if Snakemake raised an exception
        sys.stdout = _saved_py_out
        sys.stderr = _saved_py_err
        os.dup2(_saved_fd_out, 1)
        os.dup2(_saved_fd_err, 2)
        os.close(_saved_fd_out)
        os.close(_saved_fd_err)

    # ── Terminal: summary ─────────────────────────────────────────────────────
    if status:
        logger.success("Finished script!")
        _tprint("")
        _tprint("Run complete.")
        results_tsv = f"{formatted_outdir}results.tsv"
        try:
            rdf = pd.read_csv(results_tsv, sep="\t")
            n_success = int((rdf["status"] == "SUCCESS").sum())
            n_unassign = int((rdf["status"] == "UNASSIGNABLE").sum())
            n_failed = int((rdf["status"] == "FAILED").sum())
            _tprint(f"  Sequences submitted:    {query_count}")
            _tprint(f"  Genotyped successfully: {n_success}")
            if n_unassign:
                _tprint(f"  Unassignable:           {n_unassign}  (warnings in report)")
            if n_failed:
                _tprint(f"  Failed (QC or BLAST):   {n_failed}  (details in report)")
            _tprint(f"  Report:  {formatted_outdir}report.html")
            _tprint(f"  Results: {results_tsv}")
        except Exception as _e:
            _tprint(f"  (Could not read results summary: {_e})")
        return 0
    else:
        logger.error("Snakemake workflow failure")
        _tprint("")
        _tprint("ERROR: Pipeline failed. Check log for details:")
        _tprint(f"  {logfile}")
        return 1


if __name__ == "__main__":
    main()
