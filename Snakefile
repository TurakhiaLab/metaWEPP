# Snakefile

import os
import csv, re, itertools
import subprocess
import pandas as pd
import shutil

from pathlib import Path
import sys, glob

from collections import defaultdict
import itertools, json

# 1) Load METAWEPP config
configfile: "config/config.yaml"

# 2) Constants from config
DIR = config.get("DIR")
KRAKEN_DB = config.get("KRAKEN_DB")

requested_rules = set()
for arg in sys.argv[1:]:
    if not arg.startswith("-"):
        requested_rules.add(arg)

rule all:
    input:
        f"results/{DIR}/classification_proportions.png",
        f"results/{DIR}/.split_read.done",
        f"results/{DIR}/run_wepp.txt",
        f"results/{DIR}/.wepp.done",
        f"results/{DIR}/.wepp_dashboard.done"

PATHOGENS = config["PATHOGENS"].split(",")
CLADE_LIST = config["CLADE_LIST"].split(",")
CLADE_IDX = [int(x) for x in config["CLADE_IDX"].split(",")]

# Check if the lengths of all WEPP variables for different pathogens are equal
if not (len(PATHOGENS) == len(CLADE_LIST) == len(CLADE_IDX)):
    print(
        f"ERROR: PATHOGENS, CLADE_LIST, CLADE_IDX must have equal length.\n"
        f"PATHOGENS={len(PATHOGENS)}, CLADE_LIST={len(CLADE_LIST)}, CLADE_IDX={len(CLADE_IDX)}"
    )
    sys.exit(1)

WEPP_ROOT        = "WEPP"
WEPP_WORKFLOW    = "WEPP/workflow/Snakefile"
WEPP_RESULTS_DIR = "WEPP/results"
WEPP_DATA_DIR    = "WEPP/data"
WEPP_CMD_LOG     = "WEPP/cmd_log"
CONFIG           = "config/config.yaml"
IS_SINGLE_END    = config.get("SEQUENCING_TYPE", "d").lower() in {"s", "n"}

RUNNING_TEST = "test" in requested_rules

FQ1 = ""
FQ2 = ""

def gzip_if_needed(file_path):
    """Gzips a file if it's not already gzipped, returns new path."""
    if not str(file_path).endswith(".gz"):
        print(f"Compressing {file_path}...", file=sys.stderr)
        try:
            subprocess.run(["gzip", str(file_path)], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error compressing {file_path}: {e}", file=sys.stderr)
            raise
        return Path(str(file_path) + ".gz")
    return file_path

def check_input_files():
    if not DIR:
        raise ValueError("No DIR folder specified.\n"
        "Call snakemake with, e.g.,  --config DIR=TEST_DIR"
        )
    if not KRAKEN_DB:
        raise ValueError( "Please provide the path to the Kraken2 database via\n"
        "  --config KRAKEN_DB=<folder path>"
        )

    fq_dir = Path("data") / config["DIR"]
    if not fq_dir.exists():
        raise FileNotFoundError(f"Input folder {DIR} does not exist")

    if IS_SINGLE_END:
        r1_files = sorted(fq_dir.glob("*.fastq*"))
        if len(r1_files) != 1:
            raise RuntimeError(
                f"{fq_dir} must contain exactly one .fastq file for single-end mode, "
                f"but found {len(r1_files)}."
            )
        fq1 = str(gzip_if_needed(r1_files[0]))
        fq2 = ""
    else: # Paired-end
        r1_files = sorted(fq_dir.glob("*_R1.fastq*"))
        if len(r1_files) != 1:
            raise RuntimeError(
                f"{fq_dir} must contain exactly one *_R1.fastq* file for paired-end mode, "
                f"but found {len(r1_files)}."
            )

        r2_files = sorted(fq_dir.glob("*_R2.fastq*"))
        if len(r2_files) != 1:
            raise RuntimeError(
                f"{fq_dir} must contain exactly one *_R2.fastq* file for paired-end mode, "
                f"but found {len(r2_files)}."
            )
        fq1 = str(gzip_if_needed(r1_files[0]))
        fq2 = str(gzip_if_needed(r2_files[0]))

    print(f"Input FASTQs chosen:\n  FQ1 = {fq1}\n  FQ2 = {fq2}", file=sys.stderr)
    return fq1, fq2

if not RUNNING_TEST:
    FQ1, FQ2 = check_input_files()
    
# ────────────────────────────────────────────────────────────────
# Adds new pathogens with interative script to data/pathogens_for_wepp
# ────────────────────────────────────────────────────────────────

PATHOGEN_ROOT = Path("data/pathogens_for_wepp")
haplotype_pathogens = []

rule add_pathogen:
    output:
        f"results/{DIR}/.add_pathogen.done"
    params:
        db=KRAKEN_DB
    threads: workflow.cores
    run:
        subprocess.run(
            ["python", "scripts/add_ref_mat.py", "--db", params.db, "--threads", str(threads)],
            check=True
        )

        if PATHOGEN_ROOT.exists() and PATHOGEN_ROOT.is_dir():
            for pathogen_dir in filter(Path.is_dir, PATHOGEN_ROOT.iterdir()):
                fasta_files = list(pathogen_dir.glob("*.[fF][aAn]*"))
                pb_files = list(pathogen_dir.glob("*.pb*"))

                num_fasta = len(fasta_files)
                num_pb = len(pb_files)

                if num_fasta == 1 and num_pb == 1:
                    haplotype_pathogens.append(pathogen_dir)
                    continue

                if num_fasta != 1:
                    print(f"\n[WARN]: Expected 1 FASTA file in {pathogen_dir}, but found {num_fasta}.\n", file=sys.stderr)
                
                if num_pb != 1:
                    print(f"\n[WARN]: Expected 1 .pb or .pb.gz file in {pathogen_dir}, but found {num_pb}.\n", file=sys.stderr)
        
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).touch()

rule kraken:
    input:
        r1 = FQ1,
        r2 = lambda wc: [] if IS_SINGLE_END else [FQ2],
        add_pathogen_done = rules.add_pathogen.output,
    output:
        kraken_report = f"results/{DIR}/kraken_report.txt",
        kraken_out    = f"results/{DIR}/kraken_output.txt",
    threads: workflow.cores
    params:
        db        = KRAKEN_DB,
        mode_flag = "" if IS_SINGLE_END else "--paired",
    shell:
        r"""
        kraken2 --db {params.db} --threads {threads} {params.mode_flag} \
                {input.r1} {input.r2} \
                --report {output.kraken_report} \
                --output {output.kraken_out}
        """

rule kraken_visualization:
    input:
        kraken_report = rules.kraken.output.kraken_report,
    output:
        png = f"results/{DIR}/classification_proportions.png"
    params:
        min_prop = float(config.get("MIN_PROP_FOR_WEPP", 0.01))
    shell:
        r"""
        python scripts/kraken_data_visualization.py \
            {input.kraken_report} {output.png} {params.min_prop}
        """

rule split_read:
    input:
        kraken_out = rules.kraken.output.kraken_out,
        kraken_report = rules.kraken.output.kraken_report,
        png = rules.kraken_visualization.output.png,
        r1 = FQ1,
        r2 = lambda wc: [] if IS_SINGLE_END else [FQ2],
    output:
        "results/{DIR}/pathogen_coverage.tsv",
        "results/{DIR}/.split_read.done"
    params:
        script = "scripts/split_read.py",
        sequencing_type = config.get("SEQUENCING_TYPE"),
    threads: workflow.cores
    run:
        output_dir = Path(output[0]).parent
        cmd = [
            "python",
            params.script,
            "-k",
            input.kraken_out,
            "-r",
            input.kraken_report,
            "--r1",
            input.r1,
            "-o",
            str(output_dir),
            "-t",
            str(threads),
            "-s",
            params.sequencing_type,
        ]

        if input.r2:
            cmd.extend(["--r2", input.r2[0]])

        subprocess.run(cmd, check=True)

        tsv_path = Path(output[0])
        if not tsv_path.exists():
            tsv_path.write_text("")

        Path(output[1]).touch()

rule prepare_for_wepp:
    input:
        "results/{DIR}/.split_read.done",
        coverage = "results/{DIR}/pathogen_coverage.tsv"
    output:
        "results/{DIR}/run_wepp.txt"
    run:
        MIN_DEPTH_FOR_WEPP = float(config.get("MIN_DEPTH_FOR_WEPP", 0))

        try:
            coverage_df = pd.read_csv(
                input.coverage,
                sep="\t",
                header=None,
                names=["pathogen", "depth"]
            )
        except pd.errors.EmptyDataError:
            coverage_df = pd.DataFrame(columns=["pathogen", "depth"])

        coverage_df = coverage_df.dropna()
        # If coverage_df is empty → create empty run_wepp.txt and exit
        if coverage_df.empty:
            Path(output[0]).write_text("")
            return

        # Keep Pathogens meeting depth requirement
        selected_pathogens = coverage_df[coverage_df["depth"] >= MIN_DEPTH_FOR_WEPP]["pathogen"].tolist()
        # No pathogens pass depth → empty run_wepp.txt
        if not selected_pathogens:
            Path(output[0]).write_text("")

        WEPP_METADATA = {
            PATHOGENS[i]: {
                "clade_list": CLADE_LIST[i],
                "clade_idx": CLADE_IDX[i]
            }
            for i in range(len(PATHOGENS))
        }

        DEFAULT_META = WEPP_METADATA.get("default", None)
        if DEFAULT_META is None:
            raise ValueError(
                "ERROR: PATHOGENS must include a 'default' entry in config."
            )
        
        # Write output
        with open(output[0], "w") as f:
            for pathogen in selected_pathogens:
                # Choose either specific metadata or default
                raw_meta = WEPP_METADATA.get(pathogen, DEFAULT_META)

                # Process clade list
                cl = (raw_meta["clade_list"] or "").replace(":", ",")
                clade_list_arg = f"CLADE_LIST={cl} " if cl else ""

                cl_idx = f"{raw_meta['clade_idx']}"

                # Source directories
                pathogen_dir = Path(f"data/pathogens_for_wepp/{pathogen}")
                result_dir = Path(f"results/{DIR}/{pathogen}")

                # Destination directory
                dest_dir = Path(WEPP_ROOT) / "data" / f"{DIR}_{pathogen}"
                dest_dir.mkdir(parents=True, exist_ok=True)

                # 1. Copy files from data/pathogens_for_wepp/{pathogen}
                for file in pathogen_dir.glob("*"):
                    if file.is_file():
                        shutil.copy2(file, dest_dir)

                # 2. Copy files from results/{pathogen} if present
                if result_dir.exists():
                    for file in result_dir.glob("*"):
                        if file.is_file():
                            shutil.copy2(file, dest_dir)

                tree_file = (
                    list(dest_dir.glob("*.pb")) +
                    list(dest_dir.glob("*.pb.gz"))
                )
                if not tree_file:
                    raise FileNotFoundError(
                        f"No *.pb or *.pb.gz TREE file found in {pathogen_dir}"
                    )
                tree_name = tree_file[0].name

                ref_file = (
                    list(dest_dir.glob("*.fa")) +
                    list(dest_dir.glob("*.fasta")) +
                    list(dest_dir.glob("*.fna"))
                )
                if not ref_file:
                    raise FileNotFoundError(
                        f"No reference FASTA (*.fa, *.fasta, *.fna) found in {pathogen_dir}"
                    )
                ref_name = ref_file[0].name

                taxonium_files = (
                    list(dest_dir.glob("*.jsonl")) +
                    list(dest_dir.glob("*.jsonl.gz"))
                )
                taxonium_arg = (
                    f"TAXONIUM_FILE={taxonium_files[0].name} "
                    if taxonium_files else ""
                )

                # Build command
                cmd = (
                    f"snakemake --snakefile {WEPP_WORKFLOW} --directory {WEPP_ROOT} "
                    f"--cores {workflow.cores} --use-conda "
                    f"--config DIR={DIR}_{pathogen} FILE_PREFIX=metaWEPP_run "
                    f"TREE={tree_name} REF={ref_name} "
                    f"SEQUENCING_TYPE={config.get('SEQUENCING_TYPE', 'd')} "
                    f"PRIMER_BED={config.get('PRIMER_BED')} "
                    f"MIN_AF={config.get('MIN_AF')} "
                    f"MIN_DEPTH={config.get('MIN_DEPTH')} "
                    f"MIN_Q={config.get('MIN_Q')} "
                    f"MIN_PROP={config.get('MIN_PROP')} "
                    f"MIN_LEN={config.get('MIN_LEN')} "
                    f"MAX_READS={config.get('MAX_READS')} "
                    f"{clade_list_arg}CLADE_IDX={cl_idx} {taxonium_arg}"
                    f"DASHBOARD_ENABLED={config.get('DASHBOARD_ENABLED')}\n"
                )

                f.write(cmd)

rule run_wepp:
    input:
        "results/{DIR}/run_wepp.txt"
    output:
        "results/{DIR}/.wepp.done"
    run:
        for line in open(input[0]):
            cmd = line.strip()
            if cmd:
                subprocess.check_call(cmd, shell=True)

        Path(output[0]).touch()

rule print_dashboard_instructions:
    input:
        f"results/{DIR}/.wepp.done"
    output:
        "results/{DIR}/.wepp_dashboard.done"
    run:
        if not config.get("DASHBOARD_ENABLED"):
            print("\n\n\n\nTo view the dashboard, re-run the following commands:\n")
            with open(f"results/{DIR}/run_wepp.txt", "r") as f:
                for line in f:
                    line = line.strip()
                    if "DASHBOARD_ENABLED" in line:
                        modified_line = line.replace("DASHBOARD_ENABLED=False", "DASHBOARD_ENABLED=True")
                    else:
                        modified_line = line + " DASHBOARD_ENABLED=True"
                    print(modified_line, "\n")
            print("\n")
        else:
            print("\nDashboard is already enabled. No action needed.")

        Path(output[0]).touch()


rule test:
    message: "Printing metaWEPP configuration help"
    shell:
        "python scripts/metawepp_help.py"
